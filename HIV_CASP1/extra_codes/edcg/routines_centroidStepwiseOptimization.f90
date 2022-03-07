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
! Provides routines for optimizing the placement of centroids to optimize maps.
! Used in centroid_models.
module routines_centroidStepwiseOptimization

        implicit none

        private

        public centroidSimulatedAnnealing

        contains

                subroutine centroidSimulatedAnnealing(centroids, siteCharges, logger, residualList, &
                                                      residual_offsets, residual_scalings, &
                                                      tdata, mapping_accumulator_type, residual_weights, &
                                                      startingAccept, tempRate, cluster_tol, cluster_max_iter)

                        use env_kindtypes,      only: si,dp
                        use obj_trajectory,     only: traj_training_data
                        use core_random,        only: metropolisAccept, gen_rand_int
                        use core_kmeans,        only: random_centroid_init, kmeans, sget_labels,&
                                                        centroid_random_change
                        use obj_accum,          only: accum, varAccum
                        use abs_obj_logger_real,only: logger_real
                        use obj_ll,             only: i_sp_dll
                        use fit_common,         only: fit_residuals

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),            intent(inout)  :: centroids(:,:)
                        real    (dp),            intent(inout)  :: siteCharges(:)
                        real    (dp),            intent(inout)  :: residualList(:)
                        real    (dp),            intent(in   )  :: residual_scalings(:),residual_offsets(:)
                        class(logger_real),      intent(inout)  :: logger
                        type(traj_training_data),intent(in   )  :: tdata
                        character(*),            intent(in   )  :: mapping_accumulator_type
                        real    (dp),            intent(in   )  :: residual_weights(:)
                        real    (dp),            intent(in   )  :: tempRate
                        real    (dp),            intent(in   )  :: startingAccept
                        real    (dp),            intent(in   )  :: cluster_tol
                        integer (si),            intent(in   )  :: cluster_max_iter
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! hardcoded constants controlling the annealing process. Eventually, should be made
                        ! adjustable via cmdline.
                        integer (si), parameter         :: max_SA_iterations   = 100000
                        integer (si), parameter         :: tempModifyCount     = 500
                        integer (si), parameter         :: max_move_generation_attempts = 100
                        real    (dp), parameter         :: min_acceptance_rate = 0.003
                        !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!! local variables !!!
                        ! NOTE: "prop.* " variables are the proposed move variables.
                        integer (si)                    :: nCG, nFG !num of fine grained or coarse grained sites

                        !variables to hold mapping
                        integer (si),allocatable        :: mapping(:), prop_mapping(:), accept_mapping(:)

                        real    (dp),allocatable        :: prop_centroids(:,:), accept_centroids(:,:)

                        real    (dp),allocatable        :: prop_siteCharges(:), accept_siteCharges(:)

                        real    (dp),allocatable        :: prop_residualList(:),accept_residualList(:)

                        real    (dp),allocatable        :: distance_buffer(:,:)

                        real    (dp)                    :: temperature, new_temperature
                        integer (si)                    :: SA_iteration

                        ! These accumulate the variance of accepted residuals (for calculating T change)
                        ! and the acceptance rates at each step.
                        type(varAccum)                  :: residual_accum
                        type(accum)                     :: moveAccept_accum

                        !array of linked lists to hold backmapping
                        type(i_sp_dLL),allocatable      :: backmapping(:)

                        !!! Derive bounds and allocate arrays.
                        nFG = tdata%trj%nAtoms
                        nCG = size(centroids,1)

                        allocate(prop_residualList(size(residualList)))
                        allocate(prop_siteCharges(nCG))

                        !initialize stat accumulators for describing behavior at a temperature.
                        call residual_accum%reset()
                        call moveAccept_accum%reset()

                        !randomly initialize boundaries
                        centroids = random_centroid_init(tdata%trj%refAvg,nCG)

                        !generate the mapping; all residual calculations are done using mappings, NOT boundaries.
                        call sget_labels(mapping,tdata%trj%refAvg,centroids,distance_buffer)

                        !completely calculate initial residual, populate backmapping.
                        call fit_residuals(tdata        = tdata,   &
                                           mapping      = mapping,        &
                               mapping_accumulator_type = mapping_accumulator_type,&
                                       residual_weights = residual_weights,&
                                           cgCharges    = siteCharges,    &
                                           backmapping  = backmapping,    &
                                           residualList = residualList,   &
                                           residual_offsets=residual_offsets,&
                                           residual_scalings=residual_scalings,&
                                           thin         = .true.)

                        !copy over the estabilished state to the modulated state.
                        prop_residualList = residualList
                        prop_siteCharges  = siteCharges
                        prop_centroids    = centroids
                        prop_mapping      = mapping

                        !copy over the estabilished state to the best overall state
                        accept_residualList = residualList
                        accept_siteCharges  = siteCharges
                        accept_centroids    = centroids
                        accept_mapping      = mapping


                        !initialize runtime variables for the loop.
                        temperature = estimate_centroid_starting_temp(accept_ratio = startingAccept,&
                                                                  residual_weights = residual_weights,&
                                                                      nFG          = nFG,&
                                                                      nCG          = nCG,&
                                                                      tdata         = tdata,&
                                                          mapping_accumulator_type = mapping_accumulator_type,&
                                                                      nIterations  = 100,&
                                                                      tolerance    = .01_dp,&
                                                                      control_param= 2.0_dp,&
                                                                      residual_offsets = residual_offsets,&
                                                                      residual_scalings= residual_scalings)

                        groupLoop: do SA_iteration=1,max_SA_iterations

                                !generates new possible boundaries, and keeps the mapping
                                !coherent with the proposed boundary move.
                                move_generate: block
                                        integer (si) :: attempt_status
                                        integer (si) :: attempt_count, change_index

                                        attempt_status = 1
                                        attempt_count  = 0
                                        do while (attempt_status /= 0)
                                                change_index = gen_rand_int(lower=1,upper=nCG)
                                                call centroid_random_change(centroids       = prop_centroids,&
                                                                            index_to_change = change_index,&
                                                                            candidate_points= tdata%trj%refAvg,&
                                                                            max_tries       = max_move_generation_attempts,&
                                                                            change_status   = attempt_status)

                                                attempt_count = attempt_count + 1
                                                if (attempt_count > max_move_generation_attempts) then
                                                        print*, "Move generation failure.", attempt_count, max_move_generation_attempts
                                                        cycle
                                                endif
                                        enddo
                                endblock move_generate

                                !We run iterations of lloyds algorithm to make sure the centroids are accep
                                !close to the given centroids, and so that distant parts of the molecule aren't
                                !clustered together.
                                kmeans_mapping: block
                                        real    (dp),allocatable :: buffer(:,:)
                                        if (cluster_max_iter > 0) then
                                                call kmeans(         labels = prop_mapping,&
                                                              training_data = tdata%trj%refAvg,&
                                                                  centroids = prop_centroids,&
                                                             distortion_tol = cluster_tol,&
                                                                   max_iter = cluster_max_iter,&
                                                                 initialize = .false.)
                                        else
                                                call sget_labels(   labels = prop_mapping,&
                                                                    points = tdata%trj%refAvg,&
                                                                 centroids = prop_centroids,&
                                                                    buffer = buffer)
                                        endif
                                endblock kmeans_mapping

                                !Fit charges and generate the residuals.
                                call fit_residuals(tdata          = tdata,&
                                                   mapping        = prop_mapping,&
                                         mapping_accumulator_type = mapping_accumulator_type,&
                                                 residual_weights = residual_weights,&
                                                   cgCharges      = prop_siteCharges,&
                                                   backmapping    = backmapping,&
                                                   residualList   = prop_residualList,&
                                                 residual_offsets = residual_offsets,&
                                                residual_scalings = residual_scalings,&
                                                   update         = .false.,&
                                                   thin           = .true.)

                                !accept move based on criteria
                                if (metropolisAccept(prop_residualList(1),accept_residualList(1),temperature)) then
                                        !log success
                                        call moveAccept_accum%add(.true.)

                                        !log residual
                                        call logger%add(prop_residualList(1))

                                        !store changes
                                        accept_residualList = prop_residualList
                                        accept_siteCharges  = prop_siteCharges
                                        accept_centroids    = prop_centroids
                                        accept_mapping      = prop_mapping

                                        !store residual to calcualte variance
                                        call residual_accum%add(residualList(1))
                                        if (prop_residualList(1) < residualList(1)) then
                                                residualList = accept_residualList
                                                siteCharges  = accept_siteCharges
                                                centroids    = accept_centroids
                                                mapping      = accept_mapping
                                        endif
                                else
                                        !log failure
                                        call moveAccept_accum%add(.false.)

                                        !revert changes elsewhere
                                        prop_residualList = accept_residualList
                                        prop_siteCharges  = accept_siteCharges
                                        prop_centroids    = accept_centroids
                                endif

                                if (mod(SA_iteration,tempModifyCount) == 0) then

                                        new_temperature = adaptiveTempStep(temperature,residual_accum%getSD(),tempRate)

                                        if (new_temperature < temperature*.75_dp) then
                                                temperature = temperature*.75_dp
                                        else
                                                temperature = new_temperature
                                        endif

                                        !print*, "NEW RUNTIME TEMPERATURE", temperature
                                        print*, "LOG: centroid regression acceptance:", moveAccept_accum%getMean()
                                        if ((temperature < 0) .or. &
                                            (moveAccept_accum%getMean() < min_acceptance_rate)) then
                                                exit groupLoop
                                        endif

                                        call residual_accum%reset()
                                        call moveAccept_accum%reset()
                                endif

                        enddo groupLoop

                        residualList = accept_residualList
                        siteCharges  = accept_siteCharges
                        centroids    = accept_centroids
                        mapping      = accept_mapping

                endsubroutine centroidSimulatedAnnealing

                !subroutine spatialStepwiseDescent(boundaries, siteCharges, logger, residualList, &
                !                                filteredTraj, functionSpec, max_SD_iterations)

                !        use env_kindtypes,      only: si,dp
                !        use obj_trajectory,     only: extTraj
                !        use core_random,        only: gen_rand_bool, metropolisAccept, gen_rand_ordered_seq
                !        use obj_accum,          only: accum, varAccum
                !        use abs_obj_logger_real,only: logger_real
                !        use obj_ll,             only: i_sp_dll
                !        use core_filter,    only: boundariesToMapping

                !        implicit none

                !        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !        integer (si),       intent(inout)             :: boundaries(:)
                !        real    (dp),       intent(inout)             :: siteCharges(:)
                !        real    (dp),       intent(inout)             :: residualList(:)
                !        class(logger_real), intent(inout)             :: logger
                !        type(extTraj),      intent(in   )             :: filteredTraj
                !        real    (dp),       intent(in   )             :: functionSpec(:)
                !        integer (si),       intent(in   )             :: max_SD_iterations
                !        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !        !!! local variables !!!
                !        ! NOTE: "prop.* " variables are the proposed move variables.
                !        integer (si)                    :: nCG, nFG !num of fine grained or coarse grained sites

                !        !variables to hold mapping
                !        integer (si),allocatable        :: mapping(:), prevMappingCache(:)
                !        integer (si),allocatable        :: prop_boundaries(:)
                !        integer (si)                    :: boundaryChanges(1)

                !        !variables to hold buffers so that we dont' have a huge amount of
                !        !memory reallocation
                !        integer (si),allocatable        :: mappingChangesBuffer(:), writeBuffer(:)
                !        integer (si)                    :: mappingChangesBufferEnd

                !        !real    (dp),allocatable        :: siteCharges(:),  prop_siteCharges(:)
                !        real    (dp),allocatable        :: prop_siteCharges(:)

                !        real    (dp),allocatable        :: prop_residualList(:)

                !        integer (si)                    :: SD_iteration, residue_index

                !        ! Accumlates the number of accepted moves (to detect stationary points).
                !        type(accum)                     :: moveAccept_accum

                !        !array of linked lists to hold backmapping
                !        type(i_sp_dLL),allocatable      :: backmapping(:)

                !        !!! Derive bounds and allocate arrays.
                !        nFG = filteredTraj%nAtoms
                !        nCG = size(boundaries,1) + 1

                !        allocate(         prop_residualList(size(residualList)))

                !        allocate(     prop_siteCharges(nCG)   )
                !        allocate(     prevMappingCache(nFG)   )
                !        allocate( mappingChangesBuffer(nFG*2) )
                !        allocate( writeBuffer(6)              )

                !        !initialize stat accumulators for describing behavior at a temperature.
                !        call moveAccept_accum%reset()

                !        !generate the mapping; all residual calculations are done using mappings, NOT boundaries.
                !        mapping = boundariesToMapping(boundaries,nFG)

                !        !completely calculate initial residual, populate backmapping.
                !        call fit_residuals(filteredTraj = filteredTraj,   &
                !                           mapping      = mapping,        &
                !                           lambda       = functionSpec(1),&
                !                           cgCharges    = siteCharges,    &
                !                           backmapping  = backmapping,    &
                !                           residualList = residualList)

                !        !copy over the estabilished state to the modulated state.
                !        prop_residualList = residualList
                !        prop_siteCharges  = siteCharges
                !        prop_boundaries   = boundaries

                !        groupLoop: do SD_iteration=1,max_SD_iterations

                !                residue_index = mod(SD_iteration,nCG - 1) + 1

                !                boundaryChanges(1) = residue_index

                !                !generates new possible boundaries, and keeps the mapping
                !                !coherent with the proposed boundary move.
                !                move_generate: block
                !                        integer (si) :: attempt_count, max_move_generation_attempts
                !                        logical      :: attempt_status, forward

                !                        forward = gen_rand_bool()
                !                        attempt_status = .false.
                !                        attempt_count  = 0
                !                        max_move_generation_attempts = 2

                !                        do while (.not. attempt_status)
                !                                call genSpecificBoundaryMoveMapping(&
                !                                              boundaries     = prop_boundaries,&
                !                                              prevBoundaries = boundaries,&
                !                                              indexToMove    = residue_index,&
                !                                              forward        = forward,      &
                !                                              mapping        = mapping,      &
                !                                              prevMappingCache        = prevMappingCache, &
                !                                              mappingChangesBuffer    = mappingChangesBuffer,&
                !                                              mappingChangesBufferEnd = mappingChangesBufferEnd,&
                !                                              writeBuffer    = writeBuffer,  &
                !                                              change_status  = attempt_status)

                !                                forward = (.not. forward)

                !                                attempt_count = attempt_count + 1
                !                                if (attempt_count > max_move_generation_attempts) then
                !                                        print*, "move generation failure."
                !                                        exit groupLoop
                !                                endif
                !                        enddo
                !                end block move_generate

                !                !Fit charges and generate the residuals.
                !                call fit_residuals(filteredTraj   = filteredTraj,&
                !                                   mapping        = mapping,&
                !                                   lambda         = functionSpec(1),&
                !                                   cgCharges      = prop_siteCharges,&
                !                                   backmapping    = backmapping,&
                !                                   residualList   = prop_residualList,&
                !                                   update         = .true.,&
                !                                   mappingCache   = prevMappingCache,&
                !                                   mappingChanges = mappingChangesBuffer(1:mappingChangesBufferEnd))

                !                !accept move based on criteria
                !                if ((prop_residualList(1) - residualList(1)) < 0) then
                !                        !log success
                !                        call moveAccept_accum%reset()

                !                        !log residual
                !                        call logger%add(prop_residualList(1))

                !                        !store changes
                !                        residualList = prop_residualList
                !                        siteCharges  = prop_siteCharges
                !                        boundaries(boundaryChanges)   = prop_boundaries(boundaryChanges)
                !                else
                !                        !log failure
                !                        call moveAccept_accum%add(1)

                !                        !undo the move in the mapping
                !                        call updateMapping(&
                !                                  mapping        = mapping,          &
                !                                  mappingCache   = prevMappingCache, &
                !                                  newBoundaries  = boundaries,       &
                !                                  oldBoundaries  = prop_boundaries,  &
                !                                  boundaryChanges= boundaryChanges,  &
                !                                  mappingChangesBuffer    = mappingChangesBuffer,&
                !                                  mappingChangesBufferEnd = mappingChangesBufferEnd,&
                !                                  writeBuffer    = writeBuffer)

                !                        !under move in backmapping
                !                        call updateBackmapping(&
                !                                  backmapping    = backmapping,      &
                !                                  newMapping     = mapping,          &
                !                                  oldMappingCache= prevMappingCache, &
                !                                  changes        = mappingChangesBuffer(1:mappingChangesBufferEnd))

                !                        !revert changes elsewhere
                !                        prop_residualList = residualList
                !                        prop_siteCharges  = siteCharges
                !                        prop_boundaries(boundaryChanges) = boundaries(boundaryChanges)
                !                endif

                !                !exit if we've stationary
                !                if (moveAccept_accum%getSum() > 5*size(boundaries) ) exit groupLoop

                !        enddo groupLoop

                !        !End with absolute calculation of residual.
                !        call fit_residuals(filteredTraj   = filteredTraj,    &
                !                           mapping        = mapping,         &
                !                           lambda         = functionSpec(1), &
                !                           cgCharges      = siteCharges,     &
                !                           backmapping    = backmapping,     &
                !                           residualList   = residualList)
                !endsubroutine spatialStepwiseDescent

                !subroutine genNBoundaryMoveMapping(boundaries, &
                !                                   prevBoundaries,&
                !                                   boundaryChanges,&
                !                                   mapping, &
                !                                   prevMappingCache,&
                !                                   mappingChangesBuffer, &
                !                                   mappingChangesBufferEnd,&
                !                                   writeBuffer,            &
                !                                   maxChange, &
                !                                   change_status)

                !        use env_kindtypes,      only: si

                !        implicit none

                !        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !        integer (si),intent(inout)             :: boundaries(:)
                !        integer (si),intent(inout)             :: prevBoundaries(:)
                !        integer (si),intent(inout),allocatable :: boundaryChanges(:)

                !        integer (si),intent(inout)    :: mapping(:),prevMappingCache(:)
                !        integer (si),intent(inout)    :: mappingChangesBuffer(:), mappingChangesBufferEnd
                !        integer (si),intent(inout)    :: writeBuffer(:)

                !        integer (si),intent(in   )    :: maxChange
                !        logical,     intent(  out)    :: change_status
                !        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !        !optional argument counterparts
                !        integer (si)             :: testBoundaries(size(boundaries))

                !        !local vars
                !        integer (si)    :: nCG,nFG

                !        nCG = size(boundaries)+1
                !        nFG = size(mapping)

                !        testBoundaries = genPossibleNBoundaryMove(previousBoundaries = boundaries,&
                !                                                 maxChange = maxChange,           &
                !                                                 changes   = boundaryChanges)


                !        change_status = validateBoundaryMove(testBoundaries,nFG)

                !        if (change_status) then
                !                prevBoundaries = boundaries
                !                boundaries = testBoundaries
                !                call updateMapping(mapping,prevMappingCache,boundaries,prevBoundaries,&
                !                              boundaryChanges,mappingChangesBuffer,&
                !                              mappingChangesBufferEnd,writeBuffer)
                !        endif
                !endsubroutine genNBoundaryMoveMapping

                !subroutine genSpecificBoundaryMoveMapping(boundaries, &
                !                                          prevBoundaries,&
                !                                          indexToMove,&
                !                                          forward,&
                !                                          mapping, &
                !                                          prevMappingCache,&
                !                                          mappingChangesBuffer, &
                !                                          mappingChangesBufferEnd,&
                !                                          writeBuffer,            &
                !                                          change_status)

                !        use env_kindtypes,      only: si

                !        implicit none

                !        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !        integer (si),intent(inout)             :: boundaries(:)
                !        integer (si),intent(inout)             :: prevBoundaries(:)
                !        integer (si),intent(in   )             :: indexToMove
                !        logical,     intent(in   )             :: forward

                !        integer (si),intent(inout)             :: mapping(:),prevMappingCache(:)
                !        integer (si),intent(inout)             :: mappingChangesBuffer(:), mappingChangesBufferEnd
                !        integer (si),intent(inout)             :: writeBuffer(:)

                !        logical,     intent(  out)             :: change_status
                !        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !        !optional argument counterparts
                !        integer (si)             :: testBoundaries(size(boundaries))

                !        !local vars
                !        integer (si)    :: nFG

                !        nFG = size(mapping)

                !        testBoundaries = genPossibleSingleBoundaryMove(previousBoundaries = boundaries,&
                !                                                           indexToMove        = indexToMove,&
                !                                                           forward            = forward)

                !        change_status = validateBoundaryMove(testBoundaries,nFG)

                !        if (change_status) then
                !                prevBoundaries = boundaries
                !                boundaries = testBoundaries
                !                call updateMapping(mapping,prevMappingCache,boundaries,prevBoundaries,&
                !                              [ indexToMove ],mappingChangesBuffer,&
                !                              mappingChangesBufferEnd,writeBuffer)
                !        endif
                !endsubroutine genSpecificBoundaryMoveMapping

                !!Modifies temperature according ref f Huang et al. (1986) as mentioned in
                !!following cite. Adaptive method.
                !!
                !!Triki, Eric, Yann Collette, and Patrick Siarry. "A theoretical study on the
                !!behavior of simulated annealing leading
                !!to a new cooling schedule." European Journal of Operational Research 166.1 (2005): 77-92.
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

                function estimate_centroid_starting_temp(accept_ratio,residual_weights,nFG,nCG,&
                                                         tdata,mapping_accumulator_type,nIterations,&
                                                         residual_scalings, residual_offsets,&
                                                         tolerance,control_param)&
                                                         result(temp)

                        use env_kindtypes,      only: si, dp
                        use obj_trajectory,     only: traj_training_data
                        use core_kmeans,        only: sget_labels, random_centroid_init, centroid_random_change
                        use core_random,        only: gen_rand_int
                        use obj_ll,             only: i_sp_dLL
                        use fit_common,         only: fit_residuals, num_residuals

                        implicit none

                        !!!!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),             intent(in   )          :: accept_ratio,residual_weights(:)
                        integer (si),             intent(in   )          :: nFG, nCG, nIterations
                        class(traj_training_data),intent(in   )          :: tdata
                        character(*),             intent(in   )          :: mapping_accumulator_type
                        real    (dp),             intent(in   )          :: residual_scalings(:),residual_offsets(:)
                        real    (dp),             intent(in   ),optional :: tolerance,control_param
                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)                    :: temp

                        !variables to hold mapping
                        integer (si),allocatable        :: mapping(:)
                        real    (dp),allocatable        :: centroids(:,:), prop_centroids(:,:)

                        real    (dp),allocatable        :: distance_buffer(:,:)

                        real    (dp),allocatable        :: siteCharges(:), residualList1(:), residualList2(:)

                        real    (dp),allocatable        :: ref_struct(:,:)

                        real    (dp),allocatable        :: trans_start_record(:), trans_end_record(:)

                        type(i_sp_dLL),allocatable      :: backmapping(:)

                        integer (si)                    :: sample_iter

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! hardcoded constants which finetune the process.
                        integer (si), parameter         :: max_move_generation_attempts = 100
                        integer (si), parameter         :: max_temp_optim_iter = 50
                        !
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        allocate(ref_struct,source=tdata%trj%refAvg)

                        allocate(centroids(nCG,size(ref_struct,2)))
                        allocate(prop_centroids(nCG,size(ref_struct,2)))

                        allocate(siteCharges(nCG))

                        allocate(residualList1(num_residuals + 1))
                        allocate(residualList2(num_residuals + 1))

                        allocate(trans_start_record(nIterations))
                        allocate(trans_end_record(nIterations))

                        sample_transitions: do sample_iter=1,nIterations
                                centroids = random_centroid_init(ref_struct,nCG)
                                prop_centroids = centroids

                                call sget_labels(mapping,ref_struct,centroids,distance_buffer)

                                call fit_residuals(tdata        = tdata,   &
                                                   mapping      = mapping,        &
                                       mapping_accumulator_type = mapping_accumulator_type,&
                                               residual_weights = residual_weights,&
                                                   cgCharges    = siteCharges,    &
                                                   backmapping  = backmapping,    &
                                                   residualList = residualList1,  &
                                               residual_offsets = residual_offsets,&
                                              residual_scalings = residual_scalings,&
                                                   thin         = .true.)

                                move_attempt: block
                                        integer (si) :: attempt_status
                                        integer (si) :: attempt_count, change_index

                                        attempt_status = 1
                                        attempt_count  = 0
                                        do while (attempt_status /= 0)
                                                change_index = gen_rand_int(lower=1,upper=nCG)
                                                call centroid_random_change(centroids       = prop_centroids,&
                                                                            index_to_change = change_index,&
                                                                            candidate_points= ref_struct,&
                                                                            max_tries       = 100,&
                                                                            change_status   = attempt_status)

                                                attempt_count = attempt_count + 1
                                                if (attempt_count > max_move_generation_attempts) then
                                                        print*, "Move generation failure.", attempt_count, max_move_generation_attempts
                                                        cycle sample_transitions
                                                endif
                                        enddo
                                endblock move_attempt

                                call sget_labels(mapping,ref_struct,prop_centroids,distance_buffer)

                                call fit_residuals(tdata        = tdata,   &
                                                   mapping      = mapping,        &
                                       mapping_accumulator_type = mapping_accumulator_type,&
                                               residual_weights = residual_weights, &
                                                   cgCharges    = siteCharges,    &
                                                   backmapping  = backmapping,    &
                                                   residualList = residualList2,  &
                                               residual_offsets = residual_offsets,&
                                              residual_scalings = residual_scalings,&
                                                   thin         = .true.)

                                if (residualList1(1) < residualList2(1)) then
                                        trans_start_record(sample_iter) = residualList1(1)
                                        trans_end_record(sample_iter)   = residualList2(1)
                                else
                                        trans_start_record(sample_iter) = residualList2(1)
                                        trans_end_record(sample_iter)   = residualList1(1)
                                endif
                        enddo sample_transitions

                        optim_temp: block
                                real    (dp)                    :: temp_record(3), current_accept
                                real    (dp)                    :: control_param_

                                temp = 3*(sum((trans_start_record + trans_end_record)/2) / &
                                                                max(size(trans_start_record,1),1))

                                temp_record = temp

                                control_param_ = control_param

                                temp_optim_loop: do sample_iter=1,max_temp_optim_iter
                                        current_accept = &
                                             estimate_acceptance_rate(temp,trans_start_record,trans_end_record)

                                        if (current_accept >= 1.0_dp) current_accept = .99_dp

                                        if (abs(current_accept - accept_ratio) < tolerance) exit temp_optim_loop

                                        call update_opt_temp(temp,accept_ratio,&
                                                             current_accept,control_param_)

                                        temp_record(1) = temp_record(2)
                                        temp_record(2) = temp_record(3)
                                        temp_record(3) = temp

                                        if ((temp_record(3) - temp_record(2))*(temp_record(1) - temp_record(2)) > 0) then
                                                control_param_ = control_param_*2
                                        endif

                                enddo temp_optim_loop
                        endblock optim_temp

                endfunction estimate_centroid_starting_temp

                pure function estimate_acceptance_rate(temp,trans_start,trans_end) result(rate)

                        use env_kindtypes,      only: dp, qp

                        implicit none

                        !!!!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),  intent(in   )    :: temp, trans_end(:), trans_start(:)
                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)                    :: rate

                        real    (qp)                    :: exp_sum_end,exp_sum_start
                        
                        exp_sum_end   = sum(exp(-(real(trans_end,qp)  / real(temp,qp))))
                        exp_sum_start = sum(exp(-(real(trans_start,qp)/ real(temp,qp))))

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

endmodule routines_centroidStepwiseOptimization
