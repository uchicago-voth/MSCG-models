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
! Provides interfaces and routines to random number generation and primitive stochastic decisions.
!

module core_random

        implicit none

        private

        public gen_rand_ordered_seq, gen_rand_int, gen_rand_bool, metropolisAccept, gen_rand_int_omit,&
                        sgenBSarray, gen_nonunif_int

        interface sgenBSarray
                  module procedure sgen_BS_3array_dp
        end interface

        contains
                !generates an ordered random sequence of integers numSeq long between integers minSeq and maxSeq
                !@:                                               ~~~~~~                       ~~~~~~     ~~~~~~
                function gen_rand_ordered_seq(numSeq,minSeq,maxSeq,sep)

                        use env_kindtypes, only: dp, si

                        implicit none

                        !!!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si),intent(in   )              :: numSeq
                        integer (si),intent(in   )              :: minSeq
                        integer (si),intent(in   )              :: maxSeq
                        integer (si),intent(in   ),optional     :: sep

                        !integer (si),intent(out)   :: seq(numSeq)

                        integer (si)               :: gen_rand_ordered_seq(numSeq)
                        integer (si)               :: sep_

                        !!!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp)       :: temp
                        integer (si)       :: tempInt
                        integer (si)       :: i, j
                        logical            :: diff

                        !take sep if valid.
                        sep_ = 1
                        if (present(sep)) then
                                if (sep >= 1) then 
                                        sep_ = sep
                                endif
                        endif

                        do i=1,numSeq
                                diff = .false.
                                ! genererate a random integer between minSeq and maxSeq and make sure we don't have it already
                                do while (diff.eqv..false.)
                                        call random_number(temp)
                                        tempInt = int(temp*(maxSeq-minSeq))+minSeq
                                        diff    = .true.
                                        if (i.gt.1) then
                                                do j=1,i-1
                                                        if (abs(tempInt - gen_rand_ordered_seq(j)) < sep_) then
                                                                diff = .false.
                                                                exit
                                                        endif
                                                enddo
                                        endif
                                enddo

                                gen_rand_ordered_seq(i) = tempInt
                                ! if we have more than one element we need to sort them in ascending order
                                if (i.gt.1) then
                                        do j=1,i-1
                                                !check to see if the current element is less than element j
                                                if (gen_rand_ordered_seq(i)<gen_rand_ordered_seq(j)) then
                                                        !shift the entire array from element j to i circularly 1 to the right
                                                        gen_rand_ordered_seq(j:i) = cshift(gen_rand_ordered_seq(j:i),-1)
                                                        exit
                                                endif
                                        enddo
                                endif
                        enddo

                endfunction gen_rand_ordered_seq

                !Generates a random boolean.
                function gen_rand_bool() result(randBool)

                        use env_kindtypes, only: dp

                        implicit none

                        !return vaule
                        logical                    :: randBool
                        real    (dp)               :: randReal

                        call random_number(randReal)

                        if (randReal < .5) then
                                randBool = .true.
                        else
                                randBool = .false.
                        endif

                endfunction gen_rand_bool

                !generates a random integer between upper and lower, inclusive.
                function gen_rand_int(upper,lower) result(randInt)

                        use env_kindtypes, only: si,dp

                        implicit none

                        !!!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)            :: upper
                        integer (si),intent(in),optional   :: lower
                        !!!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return vaule
                        integer (si)               :: randInt
                        real    (dp)               :: randReal

                        call random_number(randReal)

                        if (present(lower) ) then
                                randInt = nint(randReal*(upper-lower))+lower
                        else
                                randInt = nint(randReal*(upper))
                        endif

                endfunction gen_rand_int

                !returns a boolean based on the metropolis criterion:
                ! if newValue < oldValue -> return True
                ! else return true with probability  exp(-(newValue-oldValue)/temp)
                function metropolisAccept(newValue, oldValue, temp, maximization)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy Variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )           :: newValue
                        real    (dp),intent(in   )           :: oldValue
                        real    (dp),intent(in   )           :: temp
                        logical,     intent(in   ),optional  :: maximization
                        !!! End Dummy Variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        logical         :: metropolisAccept

                        !local variables
                        logical         :: maximization_
                        real    (dp)    :: newValue_, oldValue_
                        real    (dp)    :: randReal, scaledDiff

                        if (present(maximization)) then
                                maximization_ = maximization
                        else
                                maximization_ = .false.
                        endif

                        if (maximization_) then
                                newValue_ = -newValue
                                oldValue_ = -oldValue
                        else
                                newValue_ = newValue
                                oldValue_ = oldValue
                        endif

                        if (newValue < oldValue) then
                                metropolisAccept = .true.
                        else
                                call random_number(randReal)
                                scaledDiff = (newValue - oldValue)/temp

                                if (exp(- scaledDiff) > randReal) then
                                        metropolisAccept = .true.
                                        return
                                else
                                        metropolisAccept = .false.
                                        return
                                endif
                        endif

                endfunction metropolisAccept

                !Returns a random integer between upper and lower, inclusive, prohibiting returning a value in omit.
                function gen_rand_int_omit(upper,lower,omit) result(randInt)

                        use env_kindtypes, only: si

                        implicit none

                        !!!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in)            :: upper
                        integer (si),intent(in),optional   :: lower
                        integer (si),intent(in)            :: omit(:)
                        !!!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return vaule
                        integer (si)    :: randInt

                        !local vars
                        integer (si)    :: lower_
                        integer (si)    :: raw_randInt, shift

                        if (present(lower)) then
                                lower_ = lower
                        else
                                lower_ = 1
                        endif

                        !we generate an index looking at 
                        raw_randInt = gen_rand_int(upper-size(omit,1),lower_)

                        shift = boolSum(omit<=raw_randInt)

                        randInt = raw_randInt + shift

                        return

                endfunction gen_rand_int_omit

                pure function boolSum(boolArray)

                        use env_kindtypes, only: si

                        implicit none

                        !!!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        logical,intent(in)            :: boolArray(:)
                        !!!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si)    :: boolSum

                        !local vars
                        integer (si)    :: i

                        boolSum = 0
                        do i=1,size(boolArray)
                                if (boolArray(i)) boolSum = boolSum + 1
                        enddo

                        return
                endfunction boolSum

                subroutine sgen_BS_3array_dp(targetArray,sourceArray,margin)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout),allocatable :: targetArray(:,:,:)
                        integer (si),intent(in   )             :: sourceArray(:,:,:)
                        integer (si),intent(in   )             :: margin
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si),allocatable               :: samplingIndices(:)
                        integer (si)                           :: iter

                        if (allocated(targetArray)) then
                                if (any(.not.(shape(targetArray) == shape(sourceArray)))) then
                                        deallocate(targetArray)
                                        allocate(targetArray(size(sourceArray,1),&
                                                             size(sourceArray,2),&
                                                             size(sourceArray,3)))
                                endif
                        else
                                allocate(targetArray(size(sourceArray,1),&
                                                     size(sourceArray,2),&
                                                     size(sourceArray,3)))
                        endif

                        allocate(samplingIndices(size(sourceArray,margin)))

                        do iter=1,size(samplingIndices)
                                samplingIndices(iter) = gen_rand_int(size(samplingIndices),1)
                        enddo

                        if      (margin == 1) then
                                targetArray = sourceArray(samplingIndices,:,:)
                        else if (margin == 2) then
                                targetArray = sourceArray(:,samplingIndices,:)
                        else if (margin == 3) then
                                targetArray = sourceArray(:,:,samplingIndices)
                        else
                                deallocate(targetArray)
                                allocate(targetArray(0,0,0))
                        endif

                endsubroutine sgen_BS_3array_dp

                function gen_nonunif_int(pmf) result(chosen_index)

                        use env_kindtypes,      only: si,dp
                        use core_sort,          only: qsort

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   ) :: pmf(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)               :: chosen_index

                        real    (dp)               :: rand
                        real    (dp),allocatable   :: cmf(:)
                        integer (si)               :: iter

                        allocate(cmf,source=pmf)
                        call qsort(cmf)

                        !construct cmf
                        do iter=2,size(cmf)
                                cmf(iter) = cmf(iter - 1) + cmf(iter)
                        enddo

                        !normalize cmf
                        cmf = cmf/cmf(size(cmf))

                        call random_number(rand)

                        !if smaller than the smallest value, end.
                        if (rand < cmf(1)) then
                                chosen_index = 1
                                return
                        endif

                        !check to see if we're between the current index and the next largest.
                        do iter=size(cmf)-1,1
                                if (rand > cmf(iter)) then
                                        chosen_index = iter + 1
                                        return
                                endif
                        enddo

                        !should never get here.
                        chosen_index = -1

                endfunction gen_nonunif_int
endmodule core_random
