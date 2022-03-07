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
! This module provides generic routines for sorting and permuting. The names of the functions dictate what 
! types they preferentially act on. Pivoting is done on the median value.
!
!   
! NOTE:
!   Currently, we use quick sort. Note that intent is in/out sometimes in places it may surprise you-- this is
!   often due to the recursive nature of the call considered: one such call will need (in the current style) 
!   such intent. This can be solved if needed by having wrappers before recursive calls are started.

module core_sort

        implicit none

        private

        public swap, permute, qsort, indexQsort, invert_permutation

        !generic quicksort. Sorts an array in place.
        interface qsort
                  module procedure qsort_int_si, qsort_real_dp
        end interface

        !index quicksort. modifies an index to report how it would sort an array.
        !This index can then be used with permute.
        interface indexQsort
                  module procedure indexQsort_int_si, indexQsort_real_dp, &
                                          indexQsort_partial_int_si, indexQsort_partial_real_dp
        end interface

        !permute an array based on the mapping given (as an integer array). Useful with
        !index sort.
        interface permute
                  module procedure permute_real_dp, permute_int_si, &
                                          permute_real_dp_target, permute_int_si_target
        end interface

        !swap two variables in an array.
        interface swap
                  module procedure swap_int_si, swap_real_dp
        end interface

        contains
                pure recursive subroutine indexQsort_real_dp(positions,record,worker)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !NOTE we need inout on positions even though we don't touch the initial--
                        !       we don't reallocate an ever recursion call.
                        real    (dp),intent(inout)          :: positions(:)
                        integer (si),intent(out)            :: record(:)
                        logical,     intent(in), optional   :: worker
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: divider
                        integer (si)               :: i
                        real    (dp),allocatable   :: positions_(:)

                        if (.not. present(worker)) then
                                !if this is our first go, then copy the array and initialize
                                !record
                                do i=1,size(record)
                                        record(i) = i
                                enddo

                                if (size(positions) <= 1) return

                                positions_ = positions

                                !note we use _positions
                                call index_fixed_partition_real_dp(positions_, record, divider)
                                call indexQsort_real_dp(positions_(:(divider-1)),record(:(divider-1)), .true.)
                                call indexQsort_real_dp(positions_(divider:),    record(divider:), .true.)
                        else
                                if (.not. worker) then
                                        !if this is our first go, then copy the array and initialize
                                        !record
                                        do i=1,size(positions)
                                                record(i) = i
                                        enddo

                                        if (size(positions) <= 1) return

                                        positions_ = positions

                                        !note we use _positions
                                        call index_fixed_partition_real_dp(positions_, record, divider)
                                        call indexQsort_real_dp(positions_(:(divider-1)),record(:(divider-1)), .true.)
                                        call indexQsort_real_dp(positions_(divider:),    record(divider:), .true.)
                                else
                                        if (size(positions) <= 1) return

                                        !note we don't use the _positions
                                        call index_fixed_partition_real_dp(positions, record, divider)
                                        call indexQsort_real_dp(positions(:(divider-1)),record(:(divider-1)), .true.)
                                        call indexQsort_real_dp(positions(divider:),    record(divider:), .true.)
                                endif
                        endif
                endsubroutine

                !returns the permutation array corresponding to an index sort.
                !
                !..
                ! @partial_size : indicates a partial sort returning the first partial_size indices
                pure subroutine indexQsort_partial_real_dp(positions,record,partial_num)

                        use env_kindtypes,      only: si,dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !NOTE we need inout on positions even though we don't touch the initial--
                        !       we don't reallocate an ever recursion call.
                        real    (dp),intent(inout)             :: positions(:)
                        integer (si),intent(inout),allocatable :: record(:)
                        integer (si),intent(in)                :: partial_num
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: divider, i
                        real    (dp),allocatable   :: positions_(:)
                        integer (si),allocatable   :: dummy_record(:)

                        if (size(positions) <= 1) then 
                                if (allocated(record)) then
                                        deallocate(record)
                                endif
                                allocate(record(1))
                                record(1) = 1
                                return
                        endif

                        !if this is our first go, then copy the array and initialize
                        !record
                        positions_ = positions
                        allocate(dummy_record(size(positions)))

                        do i=1,size(dummy_record)
                                dummy_record(i) = i
                        enddo

                        !note we use _positions
                        call index_fixed_partition_real_dp(positions_, dummy_record, divider)

                        if (partial_num < divider) then
                                call indexQsort_partial_real_dp_worker (positions_(:(divider-1)),&
                                                dummy_record(:(divider-1)), partial_num)
                        else
                                call indexQsort_partial_real_dp_worker(positions_(:(divider-1)),&
                                                dummy_record(:(divider-1)), partial_num)
                                call indexQsort_partial_real_dp_worker(positions_(divider:),&
                                                dummy_record(divider:), partial_num)
                        endif

                        record = dummy_record(1:partial_num)
                endsubroutine indexQsort_partial_real_dp

                !returns the permutation array corresponding to an index sort.
                !
                !..
                ! @partial_size : indicates a partial sort returning the first partial_size indices
                pure recursive subroutine indexQsort_partial_real_dp_worker(positions,record,partial_num)

                        use env_kindtypes,      only: si,dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !NOTE we need inout on positions even though we don't touch the initial--
                        !       we don't reallocate an ever recursion call.
                        real    (dp),intent(inout)             :: positions(:)
                        integer (si),intent(inout)             :: record(:)
                        integer (si),intent(in)                :: partial_num
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: divider

                        if (size(positions) <= 1) then 
                                return
                        endif

                        !note we don't use the _positions
                        call index_fixed_partition_real_dp(positions, record, divider)

                        if (partial_num < divider) then
                                call indexQsort_partial_real_dp_worker(positions(:(divider-1)),&
                                                record(:(divider-1)), partial_num)
                        else
                                call indexQsort_partial_real_dp_worker(positions(:(divider-1)),&
                                                record(:(divider-1)), partial_num)
                                call indexQsort_partial_real_dp_worker(positions(divider:),    &
                                                record(divider:), partial_num)
                        endif

                endsubroutine indexQsort_partial_real_dp_worker


                pure subroutine index_fixed_partition_real_dp(positions,record,marker)

                        use env_kindtypes,      only: si,dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout) :: positions(:)
                        integer (si),intent(inout) :: record(:)
                        integer (si),intent(out)   :: marker
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)       :: i,j
                        real    (dp)       :: pivot

                        pivot = positions(size(positions)/2) !fixed pivot
                        i = 0
                        j = size(positions) + 1

                        do
                                j = (j - 1)
                                do
                                        if (positions(j) <= pivot) exit
                                        j = (j - 1)
                                enddo
                                i = (i + 1)
                                do
                                        if (positions(i) >= pivot) exit
                                        i = (i + 1)
                                enddo
                                if (i < j) then
                                        call swap_real_dp(i,j,positions)
                                        call swap_int_si(i,j,record)
                                elseif (i == j) then
                                        marker = (i + 1)
                                        return
                                else
                                        marker = i
                                        return
                                endif
                        enddo

                endsubroutine

                !returns the permutation array corresponding to an index sort.
                !
                !..
                ! @partial_size : indicates a partial sort returning the first partial_size indices
                pure recursive subroutine indexQsort_int_si(positions,record,worker)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !NOTE we need inout on positions even though we don't touch the initial--
                        !       we don't reallocate an ever recursion call.
                        integer (si),intent(inout)          :: positions(:)
                        integer (si),intent(out)            :: record(:)
                        logical,     intent(in),   optional :: worker
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: divider, i
                        integer (si),allocatable   :: positions_(:)

                        if (.not. present(worker)) then
                                !if this is our first go, then copy the array and initialize
                                !record
                                do i=1,size(record)
                                        record(i) = i
                                enddo

                                if (size(positions) <= 1) return

                                positions_ = positions

                                !notej we use _positions
                                call index_fixed_partition_int_si(positions_, record, divider)
                                call indexQsort_int_si(positions_(:(divider-1)),record(:(divider-1)), .true.)
                                call indexQsort_int_si(positions_(divider:),    record(divider:), .true.)
                        else
                                if (.not. worker) then
                                        !if this is our first go, then copy the array and initialize
                                        !record
                                        do i=1,size(positions)
                                                record(i) = i
                                        enddo

                                        if (size(positions) <= 1) return

                                        positions_ = positions

                                        !note we use _positions
                                        call index_fixed_partition_int_si(positions_, record, divider)
                                        call indexQsort_int_si(positions_(:(divider-1)),record(:(divider-1)), .true.)
                                        call indexQsort_int_si(positions_(divider:),    record(divider:), .true.)
                                else
                                        if (size(positions) <= 1) return

                                        !note we don't use the _positions
                                        call index_fixed_partition_int_si(positions, record, divider)
                                        call indexQsort_int_si(positions(:(divider-1)),record(:(divider-1)), .true.)
                                        call indexQsort_int_si(positions(divider:),    record(divider:), .true.)
                                endif
                        endif
                endsubroutine

                !returns the permutation array corresponding to an index sort.
                !
                !..
                ! @partial_size : indicates a partial sort returning the first partial_size indices
                pure subroutine indexQsort_partial_int_si(positions,record,partial_num)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !NOTE we need inout on positions even though we don't touch the initial--
                        !       we don't reallocate an ever recursion call.
                        integer (si),intent(inout)             :: positions(:)
                        integer (si),intent(inout),allocatable :: record(:)
                        integer (si),intent(in)                :: partial_num
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: divider, i
                        integer (si),allocatable   :: positions_(:), dummy_record(:)

                        if (size(positions) <= 1) then 
                                if (allocated(record)) then
                                        deallocate(record)
                                endif
                                allocate(record(1))
                                record(1) = 1
                                return
                        endif

                        !if this is our first go, then copy the array and initialize
                        !record
                        positions_ = positions
                        allocate(dummy_record(size(positions)))

                        do i=1,size(dummy_record)
                                dummy_record(i) = i
                        enddo

                        !note we use _positions
                        call index_fixed_partition_int_si(positions_, dummy_record, divider)

                        if (partial_num < divider) then
                                call indexQsort_partial_int_si_worker (positions_(:(divider-1)),&
                                                dummy_record(:(divider-1)), partial_num)
                        else
                                call indexQsort_partial_int_si_worker(positions_(:(divider-1)),&
                                                dummy_record(:(divider-1)), partial_num)
                                call indexQsort_partial_int_si_worker(positions_(divider:),&
                                                dummy_record(divider:), partial_num)
                        endif

                        record = dummy_record(1:partial_num)
                endsubroutine indexQsort_partial_int_si

                !returns the permutation array corresponding to an index sort.
                !
                !..
                ! @partial_size : indicates a partial sort returning the first partial_size indices
                pure recursive subroutine indexQsort_partial_int_si_worker(positions,record,partial_num)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !NOTE we need inout on positions even though we don't touch the initial--
                        !       we don't reallocate an ever recursion call.
                        integer (si),intent(inout)             :: positions(:)
                        integer (si),intent(inout)             :: record(:)
                        integer (si),intent(in)                :: partial_num
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: divider

                        if (size(positions) <= 1) then 
                                return
                        endif

                        !note we don't use the _positions
                        call index_fixed_partition_int_si(positions, record, divider)

                        if (partial_num < divider) then
                                call indexQsort_partial_int_si_worker(positions(:(divider-1)),&
                                                record(:(divider-1)), partial_num)
                        else
                                call indexQsort_partial_int_si_worker(positions(:(divider-1)),&
                                                record(:(divider-1)), partial_num)
                                call indexQsort_partial_int_si_worker(positions(divider:),    &
                                                record(divider:), partial_num)
                        endif

                endsubroutine indexQsort_partial_int_si_worker

                pure subroutine index_fixed_partition_int_si(positions,record,marker)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout) :: positions(:)
                        integer (si),intent(inout) :: record(:)
                        integer (si),intent(out)   :: marker
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)       :: i,j, pivot
                        pivot = positions(size(positions)/2) !fixed pivot
                        i = 0
                        j = size(positions) + 1

                        do
                                j = (j - 1)
                                do
                                        if (positions(j) <= pivot) exit
                                        j = (j - 1)
                                enddo
                                i = (i + 1)
                                do
                                        if (positions(i) >= pivot) exit
                                        i = (i + 1)
                                enddo
                                if (i < j) then
                                        call swap_int_si(i,j,positions)
                                        call swap_int_si(i,j,record)
                                elseif (i == j) then
                                        marker = (i + 1)
                                        return
                                else
                                        marker = i
                                        return
                                endif
                        enddo

                endsubroutine

                pure recursive subroutine qsort_int_si(positions)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout) :: positions(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: pivot

                        if (size(positions) > 1) then !else tail case, return.
                                call fixed_partition_int_si(positions, pivot)
                                call qsort_int_si(positions(:(pivot-1)))
                                call qsort_int_si(positions(pivot:))
                        endif
                endsubroutine

                pure subroutine fixed_partition_int_si(positions,marker)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout) :: positions(:)
                        integer (si),intent(out)   :: marker
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)       :: i,j, pivot

                        pivot = positions(size(positions)/2) !fixed pivot
                        i = 0
                        j = size(positions) + 1

                        do
                                j = (j - 1)
                                do
                                        if (positions(j) <= pivot) exit
                                        j = (j - 1)
                                enddo
                                i = (i + 1)
                                do
                                        if (positions(i) >= pivot) exit
                                        i = (i + 1)
                                enddo
                                if (i < j) then
                                        call swap_int_si(i,j,positions)
                                elseif (i == j) then
                                        marker = (i + 1)
                                        return
                                else
                                        marker = i
                                        return
                                endif
                        enddo

                endsubroutine

                pure recursive subroutine qsort_real_dp(positions)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout) :: positions(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: pivot

                        if (size(positions) > 1) then !else tail case, return.
                                call fixed_partition_real_dp(positions, pivot)
                                call qsort_real_dp(positions(:(pivot-1)))
                                call qsort_real_dp(positions(pivot:))
                        endif
                endsubroutine

                pure subroutine fixed_partition_real_dp(positions,marker)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout) :: positions(:)
                        integer (si),intent(out)   :: marker
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        real    (dp)       :: pivot
                        integer (si)       :: i,j

                        pivot = positions(size(positions)/2) !fixed pivot
                        i = 0
                        j = size(positions) + 1

                        do
                                j = (j - 1)
                                do
                                        if (positions(j) <= pivot) exit
                                        j = (j - 1)
                                enddo
                                i = (i + 1)
                                do
                                        if (positions(i) >= pivot) exit
                                        i = (i + 1)
                                enddo
                                if (i < j) then
                                        call swap_real_dp(i,j,positions)
                                elseif (i == j) then
                                        marker = (i + 1)
                                        return
                                else
                                        marker = i
                                        return
                                endif
                        enddo

                endsubroutine

                !reorders positions given the indices in mapping.
                !Currently done in the original array. Mapping is real (dp).
                pure subroutine permute_int_si_target(positions,mapping,output)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)          :: positions(:)
                        integer (si),intent(in)          :: mapping(:)
                        integer (si),intent(out)         :: output(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)             :: i

                        do i=1,size(positions)
                                output(i) = positions(mapping(i))
                        enddo

                endsubroutine permute_int_si_target

                !reorders positions given the indices in mapping.
                !Currently done in the original array. Mapping is real (dp).
                pure subroutine permute_int_si(positions,mapping)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout)       :: positions(:)
                        integer (si),intent(in)          :: mapping(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)             :: i
                        integer (si),allocatable :: positions_(:)

                        allocate(positions_,source=positions)

                        do i=1,size(positions)
                                positions(i) = positions_(mapping(i))
                        enddo

                endsubroutine permute_int_si

                !reorders positions given the indices in mapping.
                !Currently not in the original array. We could be more clever about this if we need.
                pure subroutine permute_real_dp_target(positions,mapping,output)

                        use env_kindtypes,      only: dp, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in)          :: positions(:)
                        integer (si),intent(in)          :: mapping(:)
                        real    (dp),intent(out)         :: output(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)             :: i

                        do i=1,size(positions)
                                output(i) = positions(mapping(i))
                        enddo

                endsubroutine permute_real_dp_target

                !reorders positions given the indices in mapping.
                !Currently not in the original array. We could be more clever about this if we need.
                pure subroutine permute_real_dp(positions,mapping)

                        use env_kindtypes,      only: dp, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: positions(:)
                        integer (si),intent(in)             :: mapping(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)             :: i
                        real    (dp),allocatable :: positions_(:)

                        allocate(positions_,source=positions)

                        do i=1,size(positions)
                                positions(i) = positions_(mapping(i))
                        enddo

                endsubroutine permute_real_dp

                !reorders positions given the indices in permutation.
                !Currently done in the original array. Mapping is real (dp).
                pure function invert_permutation(permutation) result(inverse_perm)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )       :: permutation(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si),allocatable         :: inverse_perm(:)

                        integer (si)                     :: iter

                        allocate(inverse_perm(size(permutation,1)))

                        do iter=1,size(permutation)
                                inverse_perm(permutation(iter)) = iter
                        enddo

                endfunction invert_permutation


                pure subroutine swap_real_dp(i,j,array)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)    :: i,j
                        real    (dp),intent(inout) :: array(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        real    (dp)               :: temp

                        temp     = array(i)
                        array(i) = array(j)
                        array(j) = temp

                endsubroutine swap_real_dp

                pure subroutine swap_int_si(i,j,array)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)    :: i,j
                        integer (si),intent(inout) :: array(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)               :: temp

                        temp     = array(i)
                        array(i) = array(j)
                        array(j) = temp

                endsubroutine swap_int_si

endmodule core_sort

!unit test
!
!program main
!        use core_sort
!        use env_kindtypes
!
!        implicit none
!
!        integer (si)       :: array1(1:6)  = (/ 1, 3, 4, 2, 5, 6 /)
!        !real    (dp)       :: array2(1:6)  = (/ 6, 5, 4, 3, 2, 1 /)
!        integer (si), allocatable :: worker(:)
!
!
!
!        print*, "array1:",array1
!        !print*, "array2:",array2
!        !call qsort(array2, array1)
!        allocate(worker(size(array1)))
!        call indexQsort(array1,worker)
!        print*, "array1: return value",array1
!        print*, "worker: return value",worker
!        call indexQsort(array1,worker,2)
!        print*, "worker: return value",worker
!        call indexQsort(array1,worker,3)
!        print*, "worker: return value",worker
!        call indexQsort(array1,worker,4)
!        print*, "worker: return value",worker
!        call indexQsort(array1,worker,5)
!        print*, "worker: return value",worker
!        call indexQsort(array1,worker,6)
!        print*, "worker: return value",worker
!        !call permute(array1,worker)
!        !print*, "array1: return value",array1
!        !print*, "array2: return value",array2
!
!endprogram main
