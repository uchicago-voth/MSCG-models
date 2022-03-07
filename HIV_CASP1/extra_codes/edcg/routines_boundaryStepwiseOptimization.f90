!
!Copyright 2015 University of Chicago
!
!Licensed under the Apache License, Version 2.0 (the "License");
!you may not use this file except in compliance with the License.
!You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
!Unless required by applicable law or agreed to in writing, software
!distributed under the License is distributed on an "AS IS" BASIS,
!WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!See the License for the specific language governing permissions and
!limitations under the License.
!
! Provides routines for optimizing mappings contiguous in primary squence. 
! Used in linear_models .

module routines_boundaryStepwiseOptimization

        use env_kindtypes

        implicit none

        private

        public boundaryStepwiseDescent, boundaryStepwiseAnnealing

        contains

                subroutine boundaryStepwiseAnnealing(boundaries, siteCharges, logger, residualList, &
                                                     residual_offsets, residual_scalings, &
                                                     tdata, mapping_accumulator_type, residual_weights, &
                                                     startingAccept, tempRate)

                        use env_kindtypes,      only: si,dp
                        use obj_trajectory,     only: traj_training_data
                        use core_random,        only: gen_rand_bool, metropolisAccept, gen_rand_ordered_seq
                        use obj_accum,          only: accum, varAccum
                        use abs_obj_logger_real,only: logger_real
                        use obj_ll,             only: i_sp_dll
                        use core_filter,        only: boundariesToMapping
                        use fit_common,         only: fit_residuals

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),            intent(inout)  :: boundaries(:)
                        real    (dp),            intent(inout)  :: siteCharges(:)
                        real    (dp),            intent(inout)  :: residualList(:)
                        real    (dp),            intent(in   )  :: residual_offsets(:), residual_scalings(:)
                        class(logger_real),      intent(inout)  :: logger
                        type(traj_training_data),intent(in   )  :: tdata
                        character(*),            intent(in   )  :: mapping_accumulator_type
                        real    (dp),            intent(in   )  :: residual_weights(:)
                        real    (dp),            intent(in   )  :: tempRate
                        real    (dp),            intent(in   )  :: startingAccept
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! hardcoded constants controlling the annealing process. Eventually, should be made
                        ! adjustable via cmdline.
                        integer (si), parameter         :: max_SA_iterations   = 1000000
                        integer (si), parameter         :: tempModifyCount     = 1000
                        integer (si), parameter         :: maxChange_          = 2
                        integer (si), parameter         :: max_move_generation_attempts = 100
                        real    (dp), parameter         :: min_acceptance_rate = 0.20
                        !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!! local variables !!!
                        ! NOTE: "prop.* " variables are the proposed move variables.
                        integer (si)                    :: nCG, nFG !num of fine grained or coarse grained sites

                        !variables to hold mapping 
                        integer (si),allocatable        :: mapping(:), prevMappingCache(:)
                        integer (si),allocatable        :: boundaryChanges(:), prop_boundaries(:)

                        !variables to hold buffers so that we dont' have a huge amount of 
                        !memory reallocation
                        integer (si),allocatable        :: mappingChangesBuffer(:), writeBuffer(:)
                        integer (si)                    :: mappingChangesBufferEnd

                        !real    (dp),allocatable        :: siteCharges(:),  prop_siteCharges(:)
                        real    (dp),allocatable        :: prop_siteCharges(:)

                        real    (dp),allocatable        :: prop_residualList(:)

                        real    (dp)                    :: temperature, new_temperature, stdev
                        integer (si)                    :: SA_iteration

                        ! These accumulate the variance of accepted residuals (for calculating T change)
                        ! and the acceptance rates at each step.
                        type(varAccum)                  :: residual_accum
                        type(accum)                     :: moveAccept_accum

                        !array of linked lists to hold backmapping
                        type(i_sp_dLL),allocatable      :: backmapping(:)

                        !!! Derive bounds and allocate arrays.
                        nFG = tdata%trj%nAtoms
                        nCG = size(boundaries,1) + 1

                        allocate(         prop_residualList(size(residualList)))

                        !allocate(          siteCharges(nCG)   )
                        allocate(     prop_siteCharges(nCG)   )
                        allocate(     prevMappingCache(nFG)   )
                        allocate( mappingChangesBuffer(nFG*2) )
                        allocate( writeBuffer(maxChange_*2) )

                        !initialize stat accumulators for describing behavior at a temperature.
                        call residual_accum%reset()
                        call moveAccept_accum%reset()

                        !randomly initialize boundaries
                        boundaries = gen_rand_ordered_seq(size(boundaries),1,nFG-1,2)

                        !generate the mapping; all residual calculations are done using mappings, NOT boundaries.
                        mapping = boundariesToMapping(boundaries,nFG)

                        !completely calculate initial residual, populate backmapping.
                        call fit_residuals(            tdata = tdata,   &
                                                     mapping = mapping,        &
                                    mapping_accumulator_type = mapping_accumulator_type,&
                                            residual_weights = residual_weights,&
                                                   cgCharges = siteCharges,    &
                                                 backmapping = backmapping,    &
                                                residualList = residualList,   &
                                            residual_offsets = residual_offsets,&
                                           residual_scalings = residual_scalings,&
                                                        thin = .true.)

                        !copy over the estabilished state to the modulated state.
                        prop_residualList = residualList
                        prop_siteCharges  = siteCharges
                        prop_boundaries   = boundaries

                        !initialize runtime variables for the loop.
                        temperature = estimate_boundary_starting_temp(accept_ratio = startingAccept,&
                                                                  residual_weights = residual_weights,&
                                                                      nFG          = nFG,&
                                                                      nCG          = nCG,&
                                                                      tdata        = tdata,&
                                                   mapping_accumulator_type        = mapping_accumulator_type,&
                                                                      nIterations  = 100,&
                                                                      tolerance    = .01_dp,&
                                                                      control_param= 2.0_dp,&
                                                                      residual_offsets = residual_offsets,&
                                                                      residual_scalings= residual_scalings)

                        groupLoop: do SA_iteration=1,max_SA_iterations

                                !generates new possible boundaries, and keeps the mapping
                                !coherent with the proposed boundary move.
                                move_generate: block
                                        integer (si) :: attempt_count
                                        logical      :: attempt_status
                                        attempt_status = .false.
                                        attempt_count  = 0
                                        do while (.not. attempt_status)
                                                call genNBoundaryMoveMapping(boundaries     = prop_boundaries,&
                                                                             prevBoundaries = boundaries,&
                                                                             boundaryChanges= boundaryChanges,&
                                                                             mapping        = mapping,      &
                                                                             prevMappingCache        = prevMappingCache, &
                                                                             mappingChangesBuffer    = mappingChangesBuffer,&
                                                                             mappingChangesBufferEnd = mappingChangesBufferEnd,&
                                                                             writeBuffer    = writeBuffer,  &
                                                                             maxChange      = maxChange_,   &
                                                                             change_status  = attempt_status)

                                                attempt_count = attempt_count + 1
                                                if (attempt_count > max_move_generation_attempts) then
                                                        print*, "INFO: Annealing stopped due to move generation failure."
                                                        cycle groupLoop
                                                endif
                                        enddo
                                end block move_generate

                                !Fit charges and generate the residuals.
                                call fit_residuals(         tdata = tdata,&
                                                   mapping        = mapping,&
                                         mapping_accumulator_type = mapping_accumulator_type,&
                                                 residual_weights = residual_weights,&
                                                   cgCharges      = prop_siteCharges,&
                                                   backmapping    = backmapping,&
                                                   residualList   = prop_residualList,&
                                                 residual_offsets = residual_offsets,&
                                                residual_scalings = residual_scalings,&
                                                   update         = .true.,&
                                                   mappingCache   = prevMappingCache,&
                                                   mappingChanges = mappingChangesBuffer(1:mappingChangesBufferEnd),&
                                                   thin           = .true.)

                                !accept move based on criteria
                                if (metropolisAccept(prop_residualList(1),residualList(1),temperature)) then
                                        !log success
                                        call moveAccept_accum%add(.true.)

                                        !log residual
                                        call logger%add(prop_residualList(1))

                                        !store changes
                                        residualList = prop_residualList
                                        siteCharges  = prop_siteCharges
                                        boundaries(boundaryChanges)   = prop_boundaries(boundaryChanges)

                                        !store residual to calcualte variance
                                        call residual_accum%add(residualList(1))
                                else
                                        !log failure
                                        call moveAccept_accum%add(.false.)

                                        !undo the move in the mapping
                                        call updateMapping(&
                                                  mapping        = mapping,          &
                                                  mappingCache   = prevMappingCache, &
                                                  newBoundaries  = boundaries,       &
                                                  oldBoundaries  = prop_boundaries,  &
                                                  boundaryChanges= boundaryChanges,  &
                                                  mappingChangesBuffer    = mappingChangesBuffer,&
                                                  mappingChangesBufferEnd = mappingChangesBufferEnd,&
                                                  writeBuffer    = writeBuffer)

                                        !under move in backmapping
                                        call updateBackmapping(&
                                                  backmapping    = backmapping,      &
                                                  newMapping     = mapping,          &
                                                  oldMappingCache= prevMappingCache, &
                                                  changes        = mappingChangesBuffer(1:mappingChangesBufferEnd))

                                        !revert changes elsewhere
                                        prop_residualList = residualList 
                                        prop_siteCharges  = siteCharges
                                        prop_boundaries(boundaryChanges) = boundaries(boundaryChanges)
                                endif

                                if (mod(SA_iteration,tempModifyCount) == 0) then

                                        if (residual_accum%getSD() > tiny(stdev)) then
                                                stdev = residual_accum%getSD()
                                        else
                                                stdev = tiny(stdev)
                                        endif

                                        new_temperature = adaptiveTempStep(temperature,stdev,tempRate)

                                        if (new_temperature < temperature*.15_dp) then
                                                temperature = temperature*.15_dp
                                        else
                                                temperature = new_temperature
                                        endif

                                        print*, "INFO: SA loop acceptance/iteration: ", moveAccept_accum%getMean(), SA_iteration

                                        if ((temperature < 0) .or. &
                                            (moveAccept_accum%getMean() < min_acceptance_rate)) then
                                                print*, "INFO: Exiting SA loop. State: ", &
                                                        moveAccept_accum%getMean(), temperature
                                                exit groupLoop
                                        endif

                                        call residual_accum%reset()
                                        call moveAccept_accum%reset()
                                endif

                        enddo groupLoop

                endsubroutine boundaryStepwiseAnnealing

                subroutine boundaryStepwiseDescent(boundaries, siteCharges, logger, residualList, &
                                                   residual_offsets, residual_scalings, tdata, &
                                                   mapping_accumulator_type, residual_weights, max_SD_iterations)

                        use env_kindtypes,      only: si,dp
                        use obj_trajectory,     only: traj_training_data
                        use core_random,        only: gen_rand_bool, metropolisAccept, gen_rand_ordered_seq
                        use obj_accum,          only: accum, varAccum
                        use abs_obj_logger_real,only: logger_real
                        use obj_ll,             only: i_sp_dll
                        use core_filter,    only: boundariesToMapping
                        use fit_common,      only: fit_residuals

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),            intent(inout)  :: boundaries(:)
                        real    (dp),            intent(inout)  :: siteCharges(:)
                        real    (dp),            intent(inout)  :: residualList(:)
                        real    (dp),            intent(in   )  :: residual_offsets(:)
                        real    (dp),            intent(in   )  :: residual_scalings(:)
                        class(logger_real),      intent(inout)  :: logger
                        type(traj_training_data),intent(in   )  :: tdata
                        character(*),            intent(in   )  :: mapping_accumulator_type
                        real    (dp),            intent(in   )  :: residual_weights(:)
                        integer (si),            intent(in   )  :: max_SD_iterations
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!! local variables !!!
                        ! NOTE: "prop.* " variables are the proposed move variables.
                        integer (si)                    :: nCG, nFG !num of fine grained or coarse grained sites

                        !variables to hold mapping 
                        integer (si),allocatable        :: mapping(:), prevMappingCache(:)
                        integer (si),allocatable        :: prop_boundaries(:)
                        integer (si)                    :: boundaryChanges(1)

                        !variables to hold buffers so that we dont' have a huge amount of 
                        !memory reallocation
                        integer (si),allocatable        :: mappingChangesBuffer(:), writeBuffer(:)
                        integer (si)                    :: mappingChangesBufferEnd

                        !real    (dp),allocatable        :: siteCharges(:),  prop_siteCharges(:)
                        real    (dp),allocatable        :: prop_siteCharges(:)

                        real    (dp),allocatable        :: prop_residualList(:)

                        integer (si)                    :: SD_iteration, residue_index

                        ! Accumlates the number of accepted moves (to detect stationary points).
                        type(accum)                     :: moveAccept_accum

                        !array of linked lists to hold backmapping
                        type(i_sp_dLL),allocatable      :: backmapping(:)

                        !!! Derive bounds and allocate arrays.
                        nFG = tdata%trj%nAtoms
                        nCG = size(boundaries,1) + 1

                        allocate(         prop_residualList(size(residualList)))

                        allocate(     prop_siteCharges(nCG)   )
                        allocate(     prevMappingCache(nFG)   )
                        allocate( mappingChangesBuffer(nFG*2) )
                        allocate( writeBuffer(6)              )

                        !initialize stat accumulators for describing behavior at a temperature.
                        call moveAccept_accum%reset()

                        !generate the mapping; all residual calculations are done using mappings, NOT boundaries.
                        mapping = boundariesToMapping(boundaries,nFG)

                        !completely calculate initial residual, populate backmapping.
                        call fit_residuals(       tdata = tdata,   &
                                                mapping = mapping,        &
                               mapping_accumulator_type = mapping_accumulator_type,&
                                       residual_weights = residual_weights,&
                                              cgCharges = siteCharges,    &
                                            backmapping = backmapping,    &
                                           residualList = residualList,   &
                                       residual_offsets = residual_offsets,&
                                      residual_scalings = residual_scalings,&
                                                   thin = .true.)

                        !copy over the estabilished state to the modulated state.
                        prop_residualList = residualList
                        prop_siteCharges  = siteCharges
                        prop_boundaries   = boundaries

                        groupLoop: do SD_iteration=1,max_SD_iterations

                                residue_index = mod(SD_iteration,nCG - 1) + 1

                                boundaryChanges(1) = residue_index

                                !generates new possible boundaries, and keeps the mapping
                                !coherent with the proposed boundary move.
                                move_generate: block
                                        integer (si) :: attempt_count, max_move_generation_attempts
                                        logical      :: attempt_status, forward

                                        forward = gen_rand_bool()
                                        attempt_status = .false.
                                        attempt_count  = 0
                                        max_move_generation_attempts = 2

                                        do while (.not. attempt_status)
                                                call genSpecificBoundaryMoveMapping(&
                                                              boundaries     = prop_boundaries,&
                                                              prevBoundaries = boundaries,&
                                                              indexToMove    = residue_index,&
                                                              forward        = forward,      &
                                                              mapping        = mapping,      &
                                                              prevMappingCache        = prevMappingCache, &
                                                              mappingChangesBuffer    = mappingChangesBuffer,&
                                                              mappingChangesBufferEnd = mappingChangesBufferEnd,&
                                                              writeBuffer    = writeBuffer,  &
                                                              change_status  = attempt_status)

                                                forward = (.not. forward)

                                                attempt_count = attempt_count + 1
                                                if (attempt_count > max_move_generation_attempts) then
                                                        print*, "move generation failure."
                                                        exit groupLoop
                                                endif
                                        enddo
                                end block move_generate

                                !Fit charges and generate the residuals.
                                call fit_residuals(       tdata   = tdata,&
                                                   mapping        = mapping,&
                                         mapping_accumulator_type = mapping_accumulator_type,&
                                                 residual_weights = residual_weights,&
                                                   cgCharges      = prop_siteCharges,&
                                                   backmapping    = backmapping,&
                                                   residualList   = prop_residualList,&
                                                   residual_offsets=residual_offsets,&
                                                   residual_scalings=residual_scalings,&
                                                   update         = .true.,&
                                                   mappingCache   = prevMappingCache,&
                                                   mappingChanges = mappingChangesBuffer(1:mappingChangesBufferEnd),&
                                                   thin           = .true.)

                                !accept move based on criteria
                                if ((prop_residualList(1) - residualList(1)) < 0) then
                                        !log success
                                        call moveAccept_accum%reset()

                                        !log residual
                                        call logger%add(prop_residualList(1))

                                        !store changes
                                        residualList = prop_residualList
                                        siteCharges  = prop_siteCharges
                                        boundaries(boundaryChanges)   = prop_boundaries(boundaryChanges)
                                else
                                        !log failure
                                        call moveAccept_accum%add(1)

                                        !undo the move in the mapping
                                        call updateMapping(&
                                                  mapping        = mapping,          &
                                                  mappingCache   = prevMappingCache, &
                                                  newBoundaries  = boundaries,       &
                                                  oldBoundaries  = prop_boundaries,  &
                                                  boundaryChanges= boundaryChanges,  &
                                                  mappingChangesBuffer    = mappingChangesBuffer,&
                                                  mappingChangesBufferEnd = mappingChangesBufferEnd,&
                                                  writeBuffer    = writeBuffer)

                                        !under move in backmapping
                                        call updateBackmapping(&
                                                  backmapping    = backmapping,      &
                                                  newMapping     = mapping,          &
                                                  oldMappingCache= prevMappingCache, &
                                                  changes        = mappingChangesBuffer(1:mappingChangesBufferEnd))

                                        !revert changes elsewhere
                                        prop_residualList = residualList 
                                        prop_siteCharges  = siteCharges
                                        prop_boundaries(boundaryChanges) = boundaries(boundaryChanges)
                                endif

                                !exit if we've stationary
                                if (moveAccept_accum%getSum() > 5*size(boundaries) ) exit groupLoop

                        enddo groupLoop

                endsubroutine boundaryStepwiseDescent

                function genPossibleSingleBoundaryMove(previousBoundaries,indexToMove,forward)&
                                result(possibleMove)

                        use env_kindtypes,  only: si
                        use core_random,    only: gen_rand_int_omit

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)             :: previousBoundaries(:)
                        integer (si),intent(in)             :: indexToMove
                        logical,     intent(in),optional    :: forward
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value (STACK ALLOCATED)
                        integer (si)    :: possibleMove(size(previousBoundaries))

                        !default to the current mapping.
                        possibleMove = previousBoundaries

                        if (forward) then
                                possibleMove(indexToMove) = possibleMove(indexToMove) + 1
                        else
                                possibleMove(indexToMove) = possibleMove(indexToMove) - 1
                        endif
                endfunction genPossibleSingleBoundaryMove

                !generates a move from a distribution of 1:n subsequent moves.
                function genPossibleNBoundaryMove(previousBoundaries,maxChange,changes)&
                                result(possibleMove)

                        use env_kindtypes,      only: si
                        use core_random,        only: gen_rand_int, gen_rand_bool
                        use core_stat,          only: seq, expDistribution,&
                                                        randomSampleNoReplace
                        use core_sort,          only: qsort

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )                      :: previousBoundaries(:)
                        integer (si),intent(in   ),optional             :: maxChange
                        integer (si),intent(inout),allocatable,optional :: changes(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value (STACK ALLOCATED)
                        integer (si)                    :: possibleMove(size(previousBoundaries))

                        !optional argument counterparts
                        integer (si)                    :: maxChange_

                        !local variables
                        integer (si),allocatable        :: sample_indices(:)
                        integer (si)                    :: i, numMoves
                        logical                         :: forward

                        !!! parse options !!!!!!!!!!!!!!!!!!
                        if (present(maxChange)) then
                                if (maxChange >= 1) then
                                        maxChange_ = maxChange
                                else
                                        maxChange_ = 1
                                endif
                        else
                                maxChange_ = 1
                        endif

                        !default to the current mapping.
                        numMoves = expDistribution(maxChange_,size(previousBoundaries,1))

                        possibleMove = previousBoundaries

                        sample_indices = randomSampleNoReplace(numMoves,seq(1,size(previousBoundaries,1)))
                        call qsort(sample_indices)

                        do i=1,numMoves
                                forward = gen_rand_bool()
                                possibleMove = &
                                genPossibleSingleBoundaryMove(previousBoundaries = possibleMove,&
                                                                  indexToMove        = sample_indices(i),&
                                                                  forward            = forward)
                        enddo

                        if (present(changes)) then
                                !changes = sample_indices
                                call move_alloc(sample_indices,changes)
                        endif


                        return
                endfunction genPossibleNBoundaryMove

                pure function validateBoundaryMove(possibleMove,numTotalSites,movedIndex) result(status)
                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)             :: possibleMove(:)
                        integer (si),intent(in)             :: numTotalSites
                        integer (si),intent(in),optional    :: movedIndex
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        ! Return value
                        logical         :: status

                        !local variables
                        integer (si)    :: i

                        status = .true.

                        !If we know what was changed we can be faster
                        if (present(movedIndex)) then
                                !If we're the first boundary, make sure we don't go past the first
                                !residue
                                if (movedIndex == 1) then
                                        if (possibleMove(1) < 1) then
                                                status = .false.
                                                return
                                        elseif (possibleMove(movedIndex+1) == possibleMove(movedIndex)) then
                                                status = .false.
                                                return
                                        endif
                                !If we're the last boundary, make sure we don't go out of the protein
                                elseif (movedIndex == size(possibleMove)) then
                                        if (possibleMove(movedIndex) > numTotalSites) then
                                                status = .false.
                                                return
                                        elseif (possibleMove(movedIndex-1) == possibleMove(movedIndex)) then
                                                status = .false.
                                                return
                                        endif
                                !Test for overlapping boundaries.
                                elseif (possibleMove(movedIndex-1) == possibleMove(movedIndex) .or. &
                                                possibleMove(movedIndex+1) == possibleMove(movedIndex)) then
                                        status = .false.
                                        return
                                endif
                                return
                        !if we don't have movedIndex do global checks
                        else
                                !Test for duplicates
                                do i=1,size(possibleMove)-1
                                        if (possibleMove(i+1) - possibleMove(i) <= 1) then
                                                status = .false.
                                                return
                                        endif
                                enddo

                                !Check to seee if we went out of sane bounds for edges
                                if (possibleMove(1) < 1) then
                                        status = .false.
                                        return
                                elseif (possibleMove(size(possibleMove)) >= numTotalSites) then
                                        status = .false.
                                        return
                                endif
                        endif

                        return
                endfunction validateBoundaryMove

                subroutine genNBoundaryMoveMapping(boundaries, &
                                                   prevBoundaries,&
                                                   boundaryChanges,&
                                                   mapping, &
                                                   prevMappingCache,&
                                                   mappingChangesBuffer, &
                                                   mappingChangesBufferEnd,&
                                                   writeBuffer,            &
                                                   maxChange, &
                                                   change_status)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout)             :: boundaries(:)
                        integer (si),intent(inout)             :: prevBoundaries(:)
                        integer (si),intent(inout),allocatable :: boundaryChanges(:)

                        integer (si),intent(inout)    :: mapping(:),prevMappingCache(:)
                        integer (si),intent(inout)    :: mappingChangesBuffer(:), mappingChangesBufferEnd
                        integer (si),intent(inout)    :: writeBuffer(:)

                        integer (si),intent(in   )    :: maxChange
                        logical,     intent(  out)    :: change_status
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional argument counterparts
                        !automatic variable causes segfaults.
                        integer (si),allocatable             :: testBoundaries(:)

                        !local vars
                        integer (si)    :: nCG,nFG

                        nCG = size(boundaries)+1
                        nFG = size(mapping)

                        testBoundaries = genPossibleNBoundaryMove(previousBoundaries = boundaries,&
                                                                 maxChange = maxChange,           &
                                                                 changes   = boundaryChanges)


                        change_status = validateBoundaryMove(testBoundaries,nFG)

                        if (change_status) then
                                prevBoundaries = boundaries
                                boundaries = testBoundaries
                                call updateMapping(mapping,prevMappingCache,boundaries,prevBoundaries,&
                                              boundaryChanges,mappingChangesBuffer,&
                                              mappingChangesBufferEnd,writeBuffer)
                        endif
                endsubroutine genNBoundaryMoveMapping

                subroutine genSpecificBoundaryMoveMapping(boundaries, &
                                                          prevBoundaries,&
                                                          indexToMove,&
                                                          forward,&
                                                          mapping, &
                                                          prevMappingCache,&
                                                          mappingChangesBuffer, &
                                                          mappingChangesBufferEnd,&
                                                          writeBuffer,            &
                                                          change_status)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout)             :: boundaries(:)
                        integer (si),intent(inout)             :: prevBoundaries(:)
                        integer (si),intent(in   )             :: indexToMove
                        logical,     intent(in   )             :: forward

                        integer (si),intent(inout)             :: mapping(:),prevMappingCache(:)
                        integer (si),intent(inout)             :: mappingChangesBuffer(:), mappingChangesBufferEnd
                        integer (si),intent(inout)             :: writeBuffer(:)

                        logical,     intent(  out)             :: change_status
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional argument counterparts
                        integer (si)             :: testBoundaries(size(boundaries))

                        !local vars
                        integer (si)    :: nFG

                        nFG = size(mapping)

                        testBoundaries = genPossibleSingleBoundaryMove(previousBoundaries = boundaries,&
                                                                           indexToMove        = indexToMove,&
                                                                           forward            = forward)

                        change_status = validateBoundaryMove(testBoundaries,nFG)

                        if (change_status) then
                                prevBoundaries = boundaries
                                boundaries = testBoundaries
                                call updateMapping(mapping,prevMappingCache,boundaries,prevBoundaries,&
                                              [ indexToMove ],mappingChangesBuffer,&
                                              mappingChangesBufferEnd,writeBuffer)
                        endif
                endsubroutine genSpecificBoundaryMoveMapping

                !Modifies temperature according ref f Huang et al. (1986) as mentioned in
                !following cite. Adaptive method.
                !
                !Triki, Eric, Yann Collette, and Patrick Siarry. "A theoretical study on the
                !behavior of simulated annealing leading
                !to a new cooling schedule." European Journal of Operational Research 166.1 (2005): 77-92.
                pure function adaptiveTempStep(previousTemp, stdev, lambda) result(newTemp)
                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp), intent(in   )             :: previousTemp
                        real    (dp), intent(in   )             :: stdev
                        real    (dp), intent(in   ),optional    :: lambda
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)    :: newTemp

                        !local variables
                        real    (dp)    :: lambda_

                        if (present(lambda)) then
                                lambda_ = lambda
                        else
                                lambda_ = .6
                        endif

                        newTemp = previousTemp * (1 - previousTemp*lambda_/stdev)
                endfunction adaptiveTempStep

                subroutine sgenMappingChanges(buffer,endWriteMark,writeBuffer, newBoundaries,oldBoundaries)
                        !this is unfolded for performance
                        !it gets called in a huge number of places, and repeated memory allocation is quite
                        !slow.  Unfortunately, it suffers in its readibility; as such, extensive commenting is given.
                        
                        use env_kindtypes,      only: si
                        use core_stat,          only: buffer_seq

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        ! these must all be preallocated. No checks are done on size!

                        integer (si), intent(inout)             :: buffer(:)          !contains results in places (1:endWriteMark)
                        integer (si), intent(inout)             :: writeBuffer(:)     !temporary holding place for writes.
                        integer (si), intent(  out)             :: endWriteMark       !marks the last relevant value in buffer
                        integer (si), intent(in   )             :: newBoundaries(:), oldBoundaries(:)
                                                                                      !before and after boundaries given a move.

                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                :: subIter, lastChanged, newMark, oldMark
                        integer (si)                :: WB_EndWriteMark !marks end of meaningful sequence values

                        endWriteMark = 0

                        !we do first iteration separate as the sequences are slightly different
                        oldmark = oldboundaries(1)
                        newmark = newboundaries(1)
                        !if we moved backwards
                        if (oldmark > newmark) then
                                !generate changed indices given positions
                                call buffer_seq(writeBuffer,WB_EndWriteMark,newmark+1,oldmark)
                                !transfer positions; this could be optmized out to remove this copy if we need.
                                do subiter=1,WB_EndWriteMark
                                        endwritemark = endwritemark + 1
                                        buffer(endwritemark) = writebuffer(subiter)
                                enddo
                                !place marking so that future changes don't double count certain positions
                                lastchanged = oldmark - 1 
                        !if we moved forwards
                        elseif (oldmark < newmark) then
                                !generate changed indices given positions
                                call buffer_seq(writeBuffer,WB_EndWriteMark,oldmark+1,newmark)
                                !transfer positions
                                do subiter=1,WB_EndWriteMark
                                        endwritemark = endwritemark + 1
                                        buffer(endwritemark) = writebuffer(subiter)
                                enddo
                                !place marking so that future changes don't double count certain positions
                                lastchanged = newmark
                        else
                                lastchanged = 0
                        endif

                        !do the same for the rest of the iterations.
                        if (size(newBoundaries) > 1) then
                                block
                                        integer (si) :: iter
                                        do iter=2,size(newboundaries)
                                                oldmark = oldboundaries(iter)
                                                newmark = newboundaries(iter)
                                                if (oldmark > newmark) then
                                                        call buffer_seq(writeBuffer,WB_EndWriteMark,max(lastchanged+1,newmark+1),oldmark)
                                                        do subiter=1,WB_EndWriteMark
                                                                endwritemark = endwritemark + 1
                                                                buffer(endwritemark) = writebuffer(subiter)
                                                        enddo
                                                        lastchanged = oldmark 
                                                elseif (oldmark < newmark) then
                                                        call buffer_seq(writeBuffer,WB_EndWriteMark,max(lastchanged+1,oldmark+1),newMark)
                                                        do subiter=1,WB_EndWriteMark
                                                                endwritemark = endwritemark + 1
                                                                buffer(endwritemark) = writebuffer(subiter)
                                                        enddo
                                                        lastchanged = newmark
                                                endif
                                        enddo
                                end block
                        endif
                endsubroutine sgenMappingChanges

                subroutine updateMapping(mapping,mappingCache,newBoundaries,oldBoundaries,&
                                               boundaryChanges,mappingChangesBuffer,&
                                               mappingChangesBufferEnd,writeBuffer)

                        !memory is assumed preallocated. No checks are done.

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si), intent(inout) :: mapping(:), mappingCache(:)
                        integer (si), intent(in   ) :: newBoundaries(:), oldBoundaries(:), boundaryChanges(:)
                        integer (si), intent(inout) :: mappingChangesBuffer(:) !contains results in places (1:endWriteMark)
                        integer (si), intent(inout) :: mappingChangesBufferEnd
                        integer (si), intent(inout) :: writeBuffer(:) !temporary holding place for writes.
                        !!! End Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                :: iter
                        integer (si)                :: boundaryIndex, currentBoundary, changeIndex
                        integer (si)                :: CGid

                        call sgenMappingChanges(mappingChangesBuffer,mappingChangesBufferEnd,&
                                                writeBuffer,newBoundaries(boundaryChanges),&
                                                oldBoundaries(boundaryChanges))

                        do iter=1,mappingChangesBufferEnd
                                !save old mapping positions in the cache.
                                mappingCache(mappingChangesBuffer(iter)) = mapping(mappingChangesBuffer(iter))
                        enddo

                        boundaryIndex   = 1
                        currentBoundary = newBoundaries(boundaryIndex)
                        CGid            = 1

                        do iter=1,mappingChangesBufferEnd
                                changeIndex = mappingChangesBuffer(iter)
                                findLabel: do while (changeIndex > currentBoundary)
                                        if (boundaryIndex == size(newBoundaries)) then
                                                CGid = CGid + 1
                                                exit findLabel
                                        else
                                                CGid = CGid + 1
                                                boundaryIndex = boundaryIndex + 1
                                                currentBoundary = newBoundaries(boundaryIndex)
                                                cycle findLabel
                                        endif
                                enddo findLabel
                                mapping(changeIndex) = CGid
                        enddo
                endsubroutine updateMapping

                subroutine updateBackmapping(backmapping,newMapping,oldMappingCache,changes)
                        use env_kindtypes, only: si
                        use obj_ll,        only: i_sp_dll

                        implicit none

                        !!!!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        type(i_sp_dll),intent(inout)           :: backmapping(:)
                        integer (si),  intent(in   )           :: newMapping(:),oldMappingCache(:), changes(:)
                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)    :: FGidIndex, iter, changedID, old_Assignment, new_Assignment

                        if (size(changes) == 0) return

                        do FGidIndex=1,size(changes)
                                changedID = changes(FGidIndex)
                                old_Assignment = oldMappingCache(changedID)
                                new_Assignment = newMapping(changedID)

                                find: do iter=1,backmapping(old_Assignment)%getLength()
                                        if (backmapping(old_Assignment)%getCurrent() == changedID) then
                                                call backmapping(oldMappingCache(changes(FGidIndex)))%delete()
                                                exit find
                                        else
                                                call backmapping(old_Assignment)%next()
                                        endif
                                enddo find
                                call backmapping(new_assignment)%insert(changedID)
                        enddo
                endsubroutine updateBackmapping

                function estimate_boundary_starting_temp(accept_ratio,residual_weights,nFG,nCG,&
                                                         tdata,mapping_accumulator_type,nIterations,&
                                                         residual_offsets,residual_scalings,&
                                                         tolerance,control_param) & 
                                                         result(temp)

                        use env_kindtypes,      only: si, dp
                        use obj_ll,             only: i_sp_dll
                        use core_filter,        only: boundariesToMapping
                        use obj_trajectory,     only: traj_training_data
                        use core_random,        only: gen_rand_ordered_seq
                        use fit_common,         only: fit_residuals, num_residuals

                        implicit none

                        !!!!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),             intent(in   )            :: accept_ratio,residual_weights(:)
                        integer (si),             intent(in   )            :: nFG, nCG, nIterations
                        class(traj_training_data),intent(in   )            :: tdata
                        character(*),             intent(in   )            :: mapping_accumulator_type
                        real    (dp),             intent(in   )            :: residual_scalings(:)
                        real    (dp),             intent(in   )            :: residual_offsets(:)
                        real    (dp),             intent(in   ),optional   :: tolerance,control_param
                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)                    :: temp

                        !variables to hold mapping 
                        integer (si),allocatable        :: mapping(:), prevMappingCache(:)
                        integer (si),allocatable        :: boundaries(:), boundaryChanges(:), prop_boundaries(:)

                        !variables to hold buffers so that we dont' have a huge amount of 
                        !memory reallocation
                        integer (si),allocatable        :: mappingChangesBuffer(:), writeBuffer(:)
                        integer (si)                    :: mappingChangesBufferEnd

                        type(i_sp_dll),allocatable      :: backmapping(:)
                        real    (dp),  allocatable      :: siteCharges(:),residualList1(:),residualList2(:)

                        real    (dp),  allocatable      :: trans_start_record(:), trans_end_record(:)

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! hardcoded constants which finetune the process. 
                        integer (si), parameter         :: maxChange_          = 2
                        integer (si), parameter         :: max_move_generation_attempts = 130
                        integer (si), parameter         :: max_temp_optim_iter = 50
                        !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                    :: attempt_count, sample_iter
                        logical                         :: attempt_status


                        allocate(siteCharges(nCG))

                        allocate(residualList1(num_residuals + 1))
                        allocate(residualList2(num_residuals + 1))

                        allocate(trans_start_record(nIterations))
                        allocate(trans_end_record(nIterations))

                        allocate(boundaries(nCG-1))
                        allocate(prop_boundaries(nCG-1))

                        allocate(     prevMappingCache(nFG)   )
                        allocate( mappingChangesBuffer(nFG*2) )
                        allocate( writeBuffer(maxChange_*2) )


                        sample_iter = 1

                        sample_transitions: do while (.true.)
                                boundaries = gen_rand_ordered_seq(nCG - 1,1,nFG - 1,2)
                                prop_boundaries = boundaries

                                mapping = boundariesToMapping(boundaries,nFG)

                                call fit_residuals(tdata        = tdata,   &
                                                   mapping      = mapping,        &
                                       mapping_accumulator_type = mapping_accumulator_type,&
                                               residual_weights = residual_weights,&
                                                   cgCharges    = siteCharges,    &
                                                   backmapping  = backmapping,    &
                                                   residualList = residualList1,  &
                                                   residual_offsets=residual_offsets,&
                                                   residual_scalings=residual_scalings,&
                                                   thin         = .true.)

                                attempt_status = .false.
                                attempt_count = 0
                                do while (.not. attempt_status)
                                        call genNBoundaryMoveMapping(boundaries     = prop_boundaries,&
                                                                     prevBoundaries = boundaries,&
                                                                     boundaryChanges= boundaryChanges,&
                                                                     mapping        = mapping,      &
                                                                     prevMappingCache        = prevMappingCache, &
                                                                     mappingChangesBuffer    = mappingChangesBuffer,&
                                                                     mappingChangesBufferEnd = mappingChangesBufferEnd,&
                                                                     writeBuffer    = writeBuffer, &
                                                                     maxChange      = maxChange_,   &
                                                                     change_status  = attempt_status)

                                        attempt_count = attempt_count + 1
                                        if (attempt_count > max_move_generation_attempts) then
                                                !print*, "INFO: Move generation failure (temperature derivation)"
                                                cycle sample_transitions
                                        endif
                                enddo

                                call fit_residuals(tdata        = tdata,   &
                                                   mapping      = mapping,        &
                                       mapping_accumulator_type = mapping_accumulator_type,&
                                               residual_weights = residual_weights,&
                                                   cgCharges    = siteCharges,    &
                                                   backmapping  = backmapping,    &
                                                   residualList = residualList2,  &
                                                   residual_offsets=residual_offsets,&
                                                   residual_scalings=residual_scalings,&
                                                   thin         = .true.)

                                if (residualList1(1) < residualList2(1)) then
                                        trans_start_record(sample_iter) = residualList1(1)
                                        trans_end_record(sample_iter)   = residualList2(1)
                                else
                                        trans_start_record(sample_iter) = residualList2(1)
                                        trans_end_record(sample_iter)   = residualList1(1)
                                endif

                                sample_iter = sample_iter + 1
                                if (sample_iter > nIterations) exit sample_transitions

                        enddo sample_transitions

                        temp = (sum((trans_end_record - trans_start_record)) / max(size(trans_start_record,1),1))

                        temp_optimization: block

                                use core_stat, only: outlier_mask

                                real    (dp)                    :: tempRecord(3), control_param_
                                real    (dp)                    :: current_accept
                                logical,     allocatable        :: sample_mask(:)

                                tempRecord = temp

                                control_param_ = control_param

                                sample_mask =       outlier_mask(sample=trans_start_record,percentile=70,upper=.false.) &
                                              .and. outlier_mask(sample=trans_end_record,  percentile=70,upper=.false.)

                                temp_optim_loop: do sample_iter=1,max_temp_optim_iter
                                        current_accept = &
                                             estimate_acceptance_rate(temp,trans_start_record,trans_end_record,sample_mask)

                                        if (current_accept >= 1.0_dp) current_accept = .99_dp

                                        if (abs(current_accept - accept_ratio) < tolerance) exit temp_optim_loop

                                        call update_opt_temp(temp,accept_ratio,&
                                                             current_accept,control_param_)

                                        tempRecord(1) = tempRecord(2)
                                        tempRecord(2) = tempRecord(3)
                                        tempRecord(3) = temp

                                        if ((tempRecord(3) - tempRecord(2))*(tempRecord(1) - tempRecord(2)) > 0) then
                                                control_param_ = control_param_*2
                                        endif
                                        
                                enddo temp_optim_loop

                        endblock temp_optimization

                endfunction estimate_boundary_starting_temp

                !http://download.springer.com/static/pdf/34/art%253A10.1023%252FB%253ACOAP.0000044187.23143.bd.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1023%2FB%3ACOAP.0000044187.23143.bd&token2=exp=1477953154~acl=%2Fstatic%2Fpdf%2F34%2Fart%25253A10.1023%25252FB%25253ACOAP.0000044187.23143.bd.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Farticle%252F10.1023%252FB%253ACOAP.0000044187.23143.bd*~hmac=4c259e69f7bf8870b9e826be454a148dad149ffd5489ca684c4cfc67dc07b843
                function estimate_acceptance_rate(temp,trans_start,trans_end,mask) result(rate)

                        use env_kindtypes,      only: dp, qp

                        implicit none

                        !!!!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),  intent(in   )                      :: temp, trans_end(:), trans_start(:)
                        logical,       intent(in   ),optional             :: mask(:)
                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)                    :: rate

                        real    (qp)                    :: exp_sum_end,exp_sum_start
                        
                        if (present(mask)) then
                                exp_sum_end   = sum(exp(-(real(trans_end,qp)  / real(temp,qp))),mask=mask)
                                exp_sum_start = sum(exp(-(real(trans_start,qp)/ real(temp,qp))),mask=mask)
                        else
                                exp_sum_end   = sum(exp(-(real(trans_end,qp)  / real(temp,qp))))
                                exp_sum_start = sum(exp(-(real(trans_start,qp)/ real(temp,qp))))
                        endif

                        if ((abs(exp_sum_end) < tiny(exp_sum_end)) .or. (abs(exp_sum_start) < tiny(exp_sum_start))) then
                                rate = 100
                        else
                                rate = real(exp_sum_end / exp_sum_start,dp)
                        endif

                endfunction estimate_acceptance_rate

                pure subroutine update_opt_temp(temp,target_accept,accept,control_param)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!!!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),  intent(inout)    :: temp
                        real    (dp),  intent(in   )    :: target_accept, accept, control_param
                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        temp = temp * (log(accept) / log(target_accept))**(1/control_param)

                endsubroutine update_opt_temp

endmodule routines_boundaryStepwiseOptimization
