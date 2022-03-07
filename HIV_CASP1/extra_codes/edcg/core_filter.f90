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
! This module provides functions which often operate on mappings,
! such as those aimed at collapsing masses or charges per bead.
! Routines concercing representations of maps are also here.

module core_filter

        implicit none

        private

        public valueCollapse, mapProject, computeNumUniqueSorted, composeMapping, boundariesToMapping, &
               mappingToBoundaries, fragmentMapping, signDiscretize, canonicizeMapping, mapStack

        contains
                !Wrapped public variable collapse function.
                !This function _requires_ 2 arrays. Each iten in each array stands for the value
                !associated witha case, e.g. an atom.
                !       1: A partiioning of vaules given as an integer label for each index
                !               (e.g. element type or which cg bead this point is mapped to)
                !       2: Value associated with with each case (e.g. mass)
                !Optionally, one can provide additional information:
                !       3: Weight associated with each value (useful for e.g. mean). These are normalized
                !               to have a sum of 1 IN EACH GROUP.
                !       4: And if available, give the number of unique labels in the set. If not, it is
                !               derived.
                !Addtional flags to modify behavior:
                !       1: sort: If not given, the data is assumed to be sorted by label. If given and .true.,
                !               data is sorted using quicksort.
                !       2: validate: If given, checks are run along given data to make sure they are sane.
                !               e.g. that there are the same number of vaules as non-unique labels.
                !               NOT CURRENTLY IMPLEMENTED.
                !
                !This function then derives a mean to associate with each type.
                pure function valueCollapse(indexLabels,indexValues,indexWeights,nLabels,sort,perValueNorm,constGuard) &
                                                result(values)

                        use env_kindtypes,      only: dp, si
                        use core_sort,      only: indexQsort, permute

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)            :: indexLabels(:)  !partioning data.
                        real    (dp),intent(in)            :: indexValues(:)  !vaules to mean.

                        real    (dp),intent(in),optional   :: indexWeights(:) !prob. weights.

                        integer (si),intent(in),optional   :: nLabels         !# unique labels
                        logical,     intent(in),optional   :: sort            !whether to sort
                        logical,     intent(in),optional   :: perValueNorm
                        logical,     intent(in),optional   :: constGuard
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si),allocatable           :: indexLabels_(:) !not using auto arrays
                        real    (dp),allocatable           :: indexValues_(:) !  because they're
                        real    (dp),allocatable           :: indexWeights_(:)!  stack allocated.
                        integer (si)                       :: nLabels_
                        logical                            :: perValueNorm_
                        logical                            :: sort_, constGuard_
                        real    (dp),allocatable           :: values(:)

                        integer (si),allocatable           :: permutation(:) !not using auto arrays

                        !parse optional arguemnts.

                        if (present(sort)) then
                                sort_ = sort
                        else
                                sort_ = .true.
                        endif

                        if (present(constGuard)) then
                                constGuard_ = constGuard
                        else
                                constGuard_ = .false.
                        endif

                        !if weights arent given, assume equal weighting in each set.
                        if (present(indexWeights)) then
                                indexWeights_ = indexWeights !this will later get normalized
                        else
                                allocate(indexWeights_(size(indexValues)))
                                indexWeights_ = 1
                        endif

                        if (present(perValueNorm)) then
                                perValueNorm_ = perValueNorm
                        else
                                perValueNorm_ = .true.
                        endif

                        indexLabels_ = indexLabels
                        indexValues_ = indexValues

                        !If sorting is asked for, sort all the arrays we're going to pass.
                        if (sort_) then
                                allocate(permutation(size(indexLabels_)))
                                !derive sorting permutation based on the labels.
                                call indexQsort(indexLabels_,permutation)

                                !sort everything by this permutation.
                                call permute(indexLabels_,  permutation)
                                call permute(indexValues_,  permutation)
                                call permute(indexWeights_, permutation)
                        endif

                        !derive #labels if we need.
                        if (present(nLabels)) then
                                nLabels_ = nLabels
                        else
                                nLabels_ = computeNumUniqueSorted(indexLabels_)
                        endif

                        !call the actual collapse function.
                        values = valueCollapseSorted(indexLabels_, &
                                                     indexValues_, &
                                                     indexWeights_,&
                                                     nLabels_,     &
                                                     perValueNorm_,&
                                                     constGuard_)
                endfunction valueCollapse

                !Takes two maps, basis_map, and map to project, and returns a projected map. The projected map 
                !has an entry for each unique type in basis, as well as implied types from gaps; basis must be positive.
                !The entries in proj_map give the value of the map_to_project for that type (as dictated by positions in basis_map).
                !We only look at the first case for each type in basis mat-- do not rely on bevahior for map pairs which have
                !nonuniform values in map_to_project for types in basis_map.
                !
                !This is quite similar to valueCollapse, but assumes all values are the same, and acts on integers.
                !It is phrased differently because of difference use cases.
                !
                ! @basis_map:           map which controls the sites definitions. Must have positive values.
                ! @map_to_project:      map which has the values which are selected for proj_map
                ! @invalid_mark:        value to put down when a site is implied by sequence in basis_map, but not explicitly
                !                       present. defaults to -1.
                ! @sort:                wheter to sort the arguments. If presorted, it can be set to false. defailts to true.
                pure function mapProject(basis_map,map_to_project,invalid_mark_value,sort) result(proj_map)

                        use env_kindtypes,      only: si
                        use core_sort,          only: indexQsort, permute

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: basis_map(:)  !partioning data.
                        integer (si),intent(in   )          :: map_to_project(:)  !vaules to mean.
                        integer (si),intent(in   ),optional :: invalid_mark_value
                        logical,     intent(in   ),optional :: sort
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return variable
                        integer (si),allocatable           :: proj_map(:)

                        !local variables
                        integer (si),allocatable           :: basis_map_(:) !not using auto arrays
                        integer (si),allocatable           :: map_to_project_(:) !  because they're
                        integer (si)                       :: invalid_mark_value_
                        logical                            :: sort_

                        integer (si),allocatable           :: permutation(:) !not using auto arrays

                        if (present(invalid_mark_value)) then
                                invalid_mark_value_ = invalid_mark_value
                        else
                                invalid_mark_value_ = -1
                        endif

                        if (present(sort)) then
                                sort_ = sort
                        else
                                sort_ = .true.
                        endif

                        allocate(basis_map_,source=basis_map)
                        allocate(map_to_project_,source=map_to_project)

                        !If sorting is asked for, sort all the arrays we're going to pass.
                        if (sort_) then
                                allocate(permutation(size(basis_map_)))
                                !derive sorting permutation based on the labels.
                                call indexQsort(basis_map_,permutation)

                                !sort everything by this permutation.
                                call permute(basis_map_,  permutation)
                                call permute(map_to_project_,  permutation)
                        endif

                        !call the actual collapse function.
                        proj_map = mapProjectSorted(     basis_map = basis_map_, &
                                                    map_to_project = map_to_project_, &
                                                invalid_mark_value = invalid_mark_value_)

                endfunction mapProject

                !Returns the number of unique integers in a sorted array.
                pure function computeNumUniqueSorted(labels) result(numUnique)

                        use env_kindtypes, only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer  (si),intent(in)    :: labels(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer  (si)               :: numUnique !# of unique labels

                        !local variables
                        integer  (si)               :: lastLabel
                        integer  (si)               :: i

                        numUnique = 1
                        lastLabel = labels(1)

                        do i=2,size(labels)
                                 if (lastLabel /= labels(i)) then
                                         numUnique = numUnique + 1
                                         lastLabel = labels(i)
                                 endif
                        enddo

                        return

                endfunction computeNumUniqueSorted

                !This function returns the averages for each set of vales. It assumes that the
                !input is sortedd by label.
                pure function mapProjectSorted(basis_map, map_to_project, invalid_mark_value) result(proj_map)

                        use env_kindtypes, only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   ) :: basis_map(:)
                        integer (si),intent(in   ) :: map_to_project(:)
                        integer (si),intent(in   ) :: invalid_mark_value
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si),allocatable:: proj_map(:)

                        integer (si)            :: i               !iterator
                        integer (si)            :: record_index,prev_label !for knowing where to write.

                        !Since basis_map is canonical and sorted, the number of values in the final value
                        allocate(proj_map(maxval(basis_map)))
                        proj_map = invalid_mark_value

                        !initialize with the first case.
                        record_index = 1

                        prev_label   = basis_map(1)
                        proj_map(record_index) = map_to_project(1)

                        record_index = record_index + 1

                        do i=2,size(basis_map) !start at two because initialization counts the first.
                                if (prev_label /= basis_map(i)) then !we in a new site domain
                                        proj_map(record_index) = map_to_project(i) !record the projected value

                                        prev_label   = basis_map(i)
                                        record_index = record_index + 1
                                endif
                        enddo

                endfunction mapProjectSorted


                !This function returns the averages for each set of vales. It assumes that the
                !input is sortedd by label.
                pure function valueCollapseSorted(indexLabels, indexValues, indexWeights, nLabels, &
                                                  perValueNorm, constGuard) &
                                               result(values)

                        use env_kindtypes, only: dp, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in) :: indexLabels(:)
                        real    (dp),intent(in) :: indexValues(:)
                        integer (si),intent(in) :: nLabels
                        real    (dp),intent(in),optional&
                                                :: indexWeights(:)
                        logical,     intent(in),optional&
                                                :: perValueNorm
                        logical,     intent(in),optional&
                                                :: constGuard
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp),allocatable:: values(:)


                        !optional argument replacement variables
                        real    (dp),allocatable:: indexWeights_(:)
                        logical                 :: perValueNorm_, constGuard_

                        !local variables
                        real    (dp)            :: valueAccumulator   !accumulates values
                        real    (dp)            :: noweight_valueAccumulator   !accumulates values
                        real    (dp)            :: weightAccumulator  !accumulates weights (normalization)
                        real    (dp)            :: noweight_accumulator !accumulates weights (normalization)

                        integer (si)            :: i               !iterator
                        integer (si)            :: recordIndex,prevLabel !for knowing where to write.

                        if (present(indexWeights)) then
                                indexWeights_ = indexWeights
                        else
                                allocate(indexWeights_(size(indexLabels)))
                                indexWeights_ = 1
                        endif

                        if (present(perValueNorm)) then
                                perValueNorm_ = perValueNorm
                        else
                                perValueNorm_ = .true.
                        endif

                        if (present(constGuard)) then
                                constGuard_ = constGuard
                        else
                                constGuard_ = .false.
                        endif

                        allocate(values(nLabels))
                        values = 0

                        !initialize with the first case.
                        valueAccumulator  = indexValues(1)*indexWeights_(1)
                        noweight_valueAccumulator = indexValues(1)
                        noweight_accumulator = 1.0_dp
                        weightAccumulator = indexWeights_(1)
                        prevLabel         = indexLabels(1)
                        recordIndex       = 1

                        if (perValueNorm_) then
                                do i=2,(size(indexValues)+1) !start at two because initialization counts the first.
                                !                               got to +1 because use use it as a sentinel.
                                        if (i == (size(indexValues)+1)) then
                                                !sentinel case
                                                !record the previous set of values.
                                                if (weightAccumulator == 0 .and. constGuard_) then
                                                        values(recordIndex) = noweight_valueAccumulator/&
                                                                               noweight_accumulator
                                                else
                                                        values(recordIndex) = valueAccumulator/weightAccumulator
                                                endif
                                                recordIndex         = recordIndex+1

                                        elseif (prevLabel == indexLabels(i)) then
                                                !if we're in a not in a new label area, accumulate vaules.
                                                valueAccumulator  = valueAccumulator&
                                                                  + indexValues(i)*indexWeights_(i)

                                                noweight_valueAccumulator  = noweight_valueAccumulator&
                                                                           + indexValues(i)

                                                weightAccumulator = weightAccumulator + indexWeights_(i)
                                                noweight_accumulator = noweight_accumulator + 1.0_dp

                                        else
                                                !if we're in a new label area, store the last set.
                                                if (weightAccumulator == 0 .and. constGuard_) then
                                                        values(recordIndex) = noweight_valueAccumulator/&
                                                                               noweight_accumulator
                                                else
                                                        values(recordIndex) = valueAccumulator/weightAccumulator
                                                endif
                                                recordIndex         = recordIndex+1 !increment write spot

                                                !then continue with the case we're on.
                                                valueAccumulator  = indexValues(i)*indexWeights_(i)

                                                noweight_valueAccumulator  = indexValues(i)
                                                noweight_accumulator = 1.0_dp

                                                weightAccumulator = indexWeights_(i)
                                                !set the current label vaule.
                                                prevlabel = indexLabels(i)
                                        endif
                                enddo
                        else
                                do i=2,(size(indexValues)+1) !start at two because initialization counts the first.
                                !                               got to +1 because use use it as a sentinel.
                                        if (i == (size(indexValues)+1)) then
                                                !sentinel case
                                                !record the previous set of values.
                                                values(recordIndex) = valueAccumulator
                                                recordIndex         = recordIndex+1

                                        elseif (prevLabel == indexLabels(i)) then
                                                !if we're in a not in a new label area, accumulate vaules.
                                                valueAccumulator  = valueAccumulator  + indexValues(i)*indexWeights_(i)
                                        else
                                                !if we're in a new label area, store the last set.
                                                values(recordIndex) = valueAccumulator
                                                recordIndex         = recordIndex+1 !increment write spot

                                                !then continue with the case we're on.
                                                valueAccumulator  = indexValues(i)*indexWeights_(i)
                                                !set the current label vaule.
                                                prevlabel = indexLabels(i)
                                        endif
                                enddo

                                return
                        endif

                endfunction valueCollapseSorted

                !pure subroutine mappingToBoundaries(boundaries,nSites,mapping,nCg,boundaries_are_beginnings)
                pure subroutine mappingToBoundaries(boundaries,nSites,mapping,nCg)
                        !This subroutine takes a series of boundaries between different labels
                        !and returns an explicit mapping.

                        !assumes that the mapping is well formed _as_ boundaries (no e.g. 1111222221111)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si), intent(inout),allocatable  :: boundaries(:)
                        integer (si), intent(  out)              :: nSites
                        integer (si), intent(in   )              :: mapping(:)
                        integer (si), intent(in   )              :: nCg
                        !If the boundaries provided are start positions of each site
                        !logical,      intent(in),optional :: boundaries_are_beginnings
                        !!! End dummary arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional argument counterparts
                        !logical                         :: boundaries_are_beginnings_

                        !local variables
                        integer (si)                    :: boundaryIndex, i

                        !if (present(boundaries_are_beginnings)) then
                        !        boundaries_are_beginnings_ = boundaries_are_beginnings
                        !else
                        !        boundaries_are_beginnings_ = .false.
                        !endif

                        !allocate(mapping(numSites))

                        nSites = size(mapping)

                        !optional arguments guard.
                        if (allocated(boundaries )) then
                                if (.not. (size(boundaries ,1) .eq. nCg-1 )) then
                                        deallocate(boundaries )
                                        allocate(boundaries (nCg-1 ))
                                endif
                        else
                                allocate(boundaries (nCg-1 ))
                        endif

                        boundaryIndex = 1
                        do i=1,(size(mapping)-1)
                                if (mapping(i) /= mapping(i+1)) then
                                        boundaries(boundaryIndex) = i
                                        boundaryIndex = boundaryIndex + 1
                                endif
                        enddo

                        !if (.not. boundaries_are_beginnings_) then
                        !        siteAssignment = 1
                        !        boundaryIndex  = 1
                        !        do i=1,size(mapping)-1
                        !                if (mapping(i) /= mapping(i+1)) then
                        !
                        !                endif
                        !        enddo
                        !else
                        !        if (boundaries(1) /= 1) then
                        !                !This means we have a hidden cgSite.
                        !                !Honestly, this shouldn't happen. Give dead results.
                        !                mapping  = -1
                        !                return
                        !        endif

                        !        siteAssignment = 1
                        !        boundaryIndex  = 1
                        !        do i=1,numSites
                        !                if (boundaryIndex+1 > size(boundaries)) then
                        !                        !end case. finish the mapping.
                        !                        mapping(i) = siteAssignment
                        !                elseif (i < boundaries(boundaryIndex+1)) then
                        !                        mapping(i) = siteAssignment
                        !                else
                        !                        siteAssignment = siteAssignment + 1
                        !                        boundaryIndex = boundaryIndex + 1
                        !                        mapping(i) = siteAssignment
                        !                endif
                        !        enddo
                        !endif
                        !return

                endsubroutine mappingToBoundaries

                pure function boundariesToMapping(boundaries,numSites,boundaries_are_beginnings) result(mapping)
                        !This function takes a series of boundaries between different labels
                        !and returns an explicit mapping.

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si), intent(in)          :: boundaries(:)
                        integer (si), intent(in)          :: numSites

                        !If the boundaries provided are start positions of each site
                        logical,      intent(in),optional :: boundaries_are_beginnings
                        !!! End dummary arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si), allocatable       :: mapping(:)

                        !optional argument counterparts
                        logical                         :: boundaries_are_beginnings_
                        !integer (si)                    :: numSites_

                        !local variables
                        integer (si)                    :: boundaryIndex, siteAssignment, i

                        if (present(boundaries_are_beginnings)) then
                                boundaries_are_beginnings_ = boundaries_are_beginnings
                        else
                                boundaries_are_beginnings_ = .false.
                        endif

                        allocate(mapping(numSites))
                        if (.not. boundaries_are_beginnings_) then
                                siteAssignment = 1
                                boundaryIndex  = 1
                                do i=1,numSites
                                        if (boundaryIndex > size(boundaries,1)) then
                                                mapping(i) = size(boundaries,1) + 1
                                        elseif (i <= boundaries(boundaryIndex)) then
                                                mapping(i) = siteAssignment
                                        else
                                                boundaryIndex  = boundaryIndex  + 1
                                                siteAssignment = siteAssignment + 1
                                                mapping(i)     = siteAssignment
                                        endif
                                enddo
                        else
                                if (boundaries(1) /= 1) then
                                        !This means we have a hidden cgSite.
                                        !Honestly, this shouldn't happen. Give dead results.
                                        mapping  = -1
                                        return
                                endif

                                siteAssignment = 1
                                boundaryIndex  = 1
                                do i=1,numSites
                                        if (boundaryIndex+1 > size(boundaries)) then
                                                !end case. finish the mapping.
                                                mapping(i) = siteAssignment
                                        elseif (i < boundaries(boundaryIndex+1)) then
                                                mapping(i) = siteAssignment
                                        else
                                                siteAssignment = siteAssignment + 1
                                                boundaryIndex = boundaryIndex + 1
                                                mapping(i) = siteAssignment
                                        endif
                                enddo
                        endif
                        return

                endfunction boundariesToMapping

                pure function composeMapping(mapping1,mapping2) result(newMapping)
                        !compose two mappings a la mapping1(mapping2(x)) function notation.

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)          :: mapping1(:)
                        integer (si),intent(in)          :: mapping2(:)
                        !!! End dummary arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si),allocatable         :: newMapping(:)

                        !local variables
                        integer (si)                     :: i

                        allocate(newMapping(size(mapping2)))

                        do i=1,size(mapping2)
                                newMapping(i) = mapping1(mapping2(i))
                        enddo
                endfunction composemapping

                !this subroutine modifies a mapping such that it is canonical-- that all ids are
                !greater that 1, and that all IDs are occupied up to the max value.
                pure subroutine canonicizeMapping(mapping)

                        use env_kindtypes,      only: si
                        use core_sort,          only: indexQsort, permute, invert_permutation

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout)       :: mapping(:)
                        !!! End dummary arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si),allocatable         :: permutation(:)
                        integer (si)                     :: current_new_value, iter

                        allocate(permutation(size(mapping)))

                        call indexQsort(mapping,permutation)
                        call permute(mapping,permutation)

                        !do first element separately
                        current_new_value = 1
                        do iter=1,size(mapping)-1
                                mapping(iter) = current_new_value
                                if (mapping(iter) /= mapping(iter + 1)) then
                                        current_new_value = current_new_value + 1
                                endif
                        enddo

                        mapping(size(mapping)) = current_new_value

                        call permute(mapping,invert_permutation(permutation))

                endsubroutine canonicizeMapping

                !this function takes two mappings of the sample length, and creates a new mapping such that
                !all and only particles which differ in either mapping differ in the resulting mapping.
                !
                ! e.g. 1 1 2 2 3 3 4 4  | -> 1 2 3 4 5 6 7 8
                !      1 2 1 2 1 2 1 2  |
                !
                ! we assume mappings are canonical-- that they have at least 1 of every type in a given integer range
                ! starting at 1, i.e. [1-n]
                !
                ! This could be improved using a key-value store that accepted legnth 2 arrays.
                !
                pure function fragmentMapping(mapping1,mapping2,assume_canon) result(newMapping)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)          :: mapping1(:), mapping2(:)
                        logical,     intent(in),optional :: assume_canon
                        !!! End dummary arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si),allocatable         :: newMapping(:)

                        !local variables

                        integer (si),allocatable         :: mapping1_(:), mapping2_(:)

                        logical                          :: assume_canon_
                        integer (si)                     :: iter, id_1, id_2
                        integer (si)                     :: next_new_id
                        integer (si)                     :: num_ele_1, num_ele_2
                        integer (si),allocatable         :: lookup_table(:,:)

                        if (present(assume_canon)) then
                                assume_canon_ = assume_canon
                        else
                                assume_canon_ = .false.
                        endif

                        allocate(newMapping(size(mapping1)))

                        if (assume_canon_) then
                                num_ele_1 = maxval(mapping1)
                                num_ele_2 = maxval(mapping2)
                        else
                                mapping1_ = mapping1
                                mapping2_ = mapping2

                                call canonicizeMapping(mapping1_)
                                call canonicizeMapping(mapping2_)

                                num_ele_1 = maxval(mapping1_)
                                num_ele_2 = maxval(mapping2_)
                        endif

                        allocate(lookup_table(num_ele_1,num_ele_2))
                        !initialize to -1, signifying unknown new id
                        lookup_table = -1

                        next_new_id = 1

                        if (assume_canon_) then
                                do iter=1,size(mapping1)

                                        id_1 = mapping1(iter)
                                        id_2 = mapping2(iter)

                                        if ( lookup_table(id_1,id_2) == -1) then
                                                lookup_table(id_1,id_2) = next_new_id
                                                next_new_id = next_new_id + 1
                                        endif

                                        newMapping(iter) = lookup_table(id_1,id_2)
                                enddo
                        else
                                do iter=1,size(mapping1)

                                        id_1 = mapping1_(iter)
                                        id_2 = mapping2_(iter)

                                        if ( lookup_table(id_1,id_2) == -1) then
                                                lookup_table(id_1,id_2) = next_new_id
                                                next_new_id = next_new_id + 1
                                        endif

                                        newMapping(iter) = lookup_table(id_1,id_2)
                                enddo
                        endif
                endfunction fragmentMapping

                pure function signDiscretize(values)

                        use env_kindtypes, only: dp, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),         intent(in   ) :: values(:)
                        !logical,     optional,intent(in   ) :: zero_nearest_assign
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si),allocatable  :: signDiscretize(:)

                        integer (si)              :: label
                        logical                   :: positives_present
                        logical                   :: negatives_present
                        logical                   :: zeros_present

                        allocate(signDiscretize(size(values)))

                        positives_present = any(values >  0.0_dp)
                        negatives_present = any(values <  0.0_dp)
                        zeros_present     = any(values == 0.0_dp)

                        label = 1

                        if (positives_present) then
                                where (values > 0.0_dp)
                                        signDiscretize = label
                                endwhere

                                label = label + 1
                        endif

                        if (negatives_present) then
                                where (values < 0.0_dp)
                                        signDiscretize = label
                                endwhere

                                label = label + 1
                        endif

                        if (zeros_present) then
                                where (values == 0.0_dp)
                                        signDiscretize = label
                                endwhere
                        endif

                endfunction signDiscretize

                !returns a map which effectively concatenates the given maps,
                !ensuring that distinct ids are distinct.
                !Map1 comes first in output map, then map2, which is shifted as needed
                !to fulfill goals.
                !
                !Empty corner cases are handled.
                !
                pure function mapStack(map1, map2)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)         :: map1(:), map2(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !currently, we do the trivial 
                        integer (si),allocatable :: mapStack(:)

                        if ((size(map1) /= 0) .and. (size(map2) /= 0)) then
                                allocate(mapStack(size(map1) + size(map2)))

                                mapStack(1:size(map1)) = map1

                                if (minval(map2) <= 0) then
                                        mapStack(size(map1)+1:) = map2 + 1 - minval(map2) + maxval(map1)
                                else
                                        mapStack(size(map1)+1:) = map2 + maxval(map1)
                                endif
                        else if (size(map1) == 0) then
                                mapStack = map2
                        else if (size(map2) == 0) then
                                mapStack = map1
                        endif

                endfunction mapStack
endmodule core_filter

!unit test

!program main
!        use core_filter
!        use env_kindtypes
!
!        implicit none
!
!        real    (dp)       :: values(1:10) = (/ 1.2, 3.5, 1.1, 1.5 ,.2, .1, .6, .6, .9, .1/)
!        real    (dp)       :: weights(1:10)= (/ 1.2, 3.5, 1.1, 1.5 ,.2, .1, .6, .6, .9, .1/)
!        integer (si)       :: types(1:10)  = (/ 1, 1, 1, 1, 2, 2, 2 ,3 ,3 ,3 /)
!
!        real    (dp),allocatable :: ret(:)
!
!        ret = valueCollapse(types,values,weights)
!        print*, "return value",ret
!
!endprogram main

