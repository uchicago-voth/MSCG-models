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
! Provides definitions for accumulator objects. These objects implement
! general online algorithms for use in e.g. loops, were a rolling mean/sd is desired.
!

module obj_accum

        use env_kindtypes, only: si, dp

        implicit none

        private

        public          :: accum, varAccum, accumV, accumM

        !basic accumulator for mean.
        type                           :: accum

                private

                integer (si)           :: nSamples
                real    (dp)           :: mean

                logical                :: tainted

                contains
                        private

                        procedure,public :: reset      => accum_reset
                        procedure        :: addReal_dp => accum_addReal_dp
                        procedure        :: addInt_si  => accum_addInt_si
                        procedure        :: addLogical => accum_addLogical
                        generic,  public :: add        => addReal_dp,&
                                                          addInt_si,&
                                                          addLogical
                        procedure,public :: getN       => accum_getN
                        procedure,public :: getSum     => accum_getSum
                        procedure,public :: getMean    => accum_getMean
        end type accum

        !extends the basic accumulator to calculate variances as it aggregates.
        type, extends(accum)           :: varAccum

                private

                real (dp)              :: var_accum

                contains
                        private

                        procedure,public :: reset      => varAccum_reset
                        procedure        :: addReal_dp => varAccum_addReal_dp
                        procedure        :: addInt_si  => varAccum_addInt_si
                        procedure        :: addLogical => varAccum_addLogical
                        procedure,public :: getVar     => varAccum_getVar
                        procedure,public :: getSD      => varAccum_getSD
        end type varAccum

        !basic accumulator for 1d arrays.
        type                           :: accumV

                private

                integer (si)             :: nSamples
                real    (dp),allocatable :: mean(:)
                logical                  :: tainted

                contains
                        private

                        procedure,public :: reset      => accumV_reset
                        procedure,public :: isTainted  => accumV_tainted
                        procedure        :: addReal_dp => accumV_addReal_dp
                        procedure        :: addInt_si  => accumV_addInt_si
                        procedure        :: addLogical => accumV_addLogical
                        generic,  public :: add        => addReal_dp,&
                                                          addInt_si,&
                                                          addLogical
                        procedure,public :: getDim     => accumV_getDim
                        procedure,public :: getN       => accumV_getN
                        procedure,public :: getSum     => accumV_getSum
                        procedure,public :: getMean    => accumV_getMean
        end type accumV

        !basic accumulator for 2d arrays.
        type                           :: accumM

                private

                integer (si)             :: nSamples
                real    (dp),allocatable :: mean(:,:)
                logical                  :: tainted

                contains
                        private

                        procedure,public :: reset      => accumM_reset
                        procedure,public :: isTainted  => accumM_tainted
                        procedure        :: addReal_dp => accumM_addReal_dp
                        procedure        :: addInt_si  => accumM_addInt_si
                        procedure        :: addLogical => accumM_addLogical
                        generic,  public :: add        => addReal_dp,&
                                                          addInt_si,&
                                                          addLogical
                        procedure,public :: getDim     => accumM_getDim
                        procedure,public :: getN       => accumM_getN
                        procedure,public :: getSum     => accumM_getSum
                        procedure,public :: getMean    => accumM_getMean
        end type accumM


        contains
                !!! Begin accum methods !!!

                elemental subroutine accum_reset(self)
                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accum),intent(inout)       :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        self%Nsamples = 0
                        self%mean     = 0
                        self%tainted  = .false.

                        return
                endsubroutine accum_reset

                elemental function accum_tainted(self)
                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accum),intent(in   )       :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        logical         :: accum_tainted

                        accum_tainted = self%tainted

                        return
                endfunction accum_tainted

                elemental subroutine accum_addReal_dp(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accum),intent(inout)      :: self
                        real    (dp),intent(in   )      :: newValue
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !local vars
                        real    (dp)    :: delta

                        self%nSamples = self%nSamples + 1
                        delta         = newValue  - self%mean
                        self%mean     = self%mean + delta/self%nSamples

                        return
                endsubroutine accum_addReal_dp

                elemental subroutine accum_addInt_si(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accum),intent(inout)      :: self
                        integer (si),intent(in   )      :: newValue
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        call self%add(real(newValue,dp))

                        return
                endsubroutine accum_addInt_si

                elemental subroutine accum_addlogical(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accum),intent(inout)      :: self
                        logical     ,intent(in   )      :: newValue
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        if (newValue) then
                                call self%add(1.0_dp)
                        else
                                call self%add(0.0_dp)
                        endif

                        return
                endsubroutine accum_addlogical

                !get number of samples.
                elemental function accum_getN(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accum),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp) :: accum_getN

                        accum_getN = self%nSamples

                        return
                endfunction accum_getN

                elemental function accum_getSum(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accum),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp) :: accum_getSum

                        accum_getSum = self%mean * real(self%nSamples,dp)

                        return
                endfunction accum_getSum

                elemental function accum_getMean(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accum),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp) :: accum_getMean

                        accum_getMean = self%mean

                        return
                endfunction accum_getMean

                !!! End accum methods !!!

                !!! Begin varAccum methods !!!

                elemental subroutine varAccum_reset(self)
                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(varAccum),intent(inout)       :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        self%nsamples  = 0
                        self%mean      = 0
                        self%var_accum = 0

                        return
                endsubroutine varAccum_reset

                !https://en.wikipedia.org/wiki/Algorithms_for_calculating_variancepure
                !Donald E. Knuth (1998). The Art of Computer Programming,
                !volume 2: Seminumerical Algorithms, 3rd edn., p. 232. Boston:
                !Addison-Wesley.
                elemental subroutine varAccum_addReal_dp(self,newValue)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(varAccum),intent(inout)      :: self
                        real       (dp),intent(in   )      :: newvalue
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !local vars
                        real    (dp)    :: delta

                        self%nSamples  = self%nSamples + 1
                        delta          = newValue  - self%mean
                        self%mean      = self%mean + delta/self%nSamples
                        self%var_accum = self%var_accum  + delta*(newValue - self%mean)

                        return
                endsubroutine varAccum_addReal_dp

                elemental subroutine varAccum_addInt_si(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(varAccum),intent(inout)      :: self
                        integer (si),   intent(in   )      :: newValue
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        call self%add(real(newValue,dp))

                        return
                endsubroutine varAccum_addInt_si

                elemental subroutine varAccum_addlogical(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(varAccum),intent(inout)      :: self
                        logical        ,intent(in   )      :: newValue
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        if (newValue) then
                                call self%add(1.0_dp)
                        else
                                call self%add(0.0_dp)
                        endif

                        return
                endsubroutine varAccum_addlogical

                elemental function varAccum_getVar(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(varAccum),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp) :: varAccum_getVar

                        if (self%nSamples > 1) then
                                varAccum_getVar = self%var_accum/(self%nSamples-1)
                        else
                                varAccum_getVar = 0
                        endif

                        return
                endfunction varAccum_getvar

                elemental function varAccum_getSD(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(varAccum),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp) :: varAccum_getSD

                        varAccum_getSD = sqrt(self%getVar())

                        return
                endfunction varAccum_getSD

                !!! End varAccum methods !!!

                !!! Begin accumV methods !!!

                pure subroutine accumV_reset(self)
                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(inout)       :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        self%Nsamples = 0

                        if (allocated(self%mean)) then
                                deallocate(self%mean)
                        endif

                        return
                endsubroutine accumV_reset

                pure function accumV_tainted(self)
                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(in   )       :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        logical         :: accumV_tainted

                        accumV_tainted = self%tainted

                        return
                endfunction accumV_tainted

                pure subroutine accumV_addReal_dp(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(inout)      :: self
                        real    (dp), intent(in   )      :: newValue(:)
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !local vars
                        real    (dp)    :: delta
                        integer (si)    :: i

                        if (.not. allocated(self%mean)) then
                                allocate(self%mean(size(newValue)))
                                self%mean = 0
                        else
                                if (size(self%mean) /= size(newValue)) then
                                        self%tainted = .true.
                                        self%mean = 0
                                        return
                                endif
                        endif

                        self%nSamples = self%nSamples + 1

                        do i=1,size(newValue)
                                delta         = newValue(i)  - self%mean(i)
                                self%mean(i)  = self%mean(i) + delta/self%nSamples
                        enddo

                        return
                endsubroutine accumV_addReal_dp

                pure subroutine accumV_addInt_si(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(inout)      :: self
                        integer (si), intent(in   )      :: newValue(:)
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        call self%add(real(newValue,dp))

                        return
                endsubroutine accumV_addInt_si

                pure subroutine accumV_addlogical(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(inout)      :: self
                        logical,      intent(in   )      :: newValue(:)
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !local vars
                        real    (dp),allocatable :: converted_value(:)
                        integer (si)             :: i

                        allocate(converted_value(size(newValue,1)))

                        do i=1,size(newValue)
                                if (newValue(i)) then
                                        converted_value(i) = 1.0_dp
                                else
                                        converted_value(i) = 0.0_dp
                                endif
                        enddo

                        call self%add(converted_value)

                        return
                endsubroutine accumV_addlogical

                !get number of samples.
                pure function accumV_getN(self)

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        integer (si) :: accumV_getN

                        accumV_getN = self%nSamples

                        return
                endfunction accumV_getN

                pure function accumV_getDim(self)
                        use env_kindtypes, only: si

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        integer (si) :: accumV_getDim

                        accumV_getDim = size(self%mean,1)

                        return
                endfunction accumV_getDim

                pure function accumV_getSum(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp),allocatable :: accumV_getSum(:)
                        
                        accumV_getSum = self%mean * real(self%nSamples,dp)

                        return
                endfunction accumV_getSum

                pure function accumV_getMean(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumV),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp),allocatable :: accumV_getMean(:)

                        accumV_getMean = self%mean

                        return
                endfunction accumV_getMean

                !!! M accumulator fuctions !!!

                pure subroutine accumM_reset(self)
                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(inout)       :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        self%Nsamples = 0

                        if (allocated(self%mean)) then
                                deallocate(self%mean)
                        endif

                        return
                endsubroutine accumM_reset

                pure function accumM_tainted(self)
                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(in   )       :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        logical         :: accumM_tainted

                        accumM_tainted = self%tainted

                        return
                endfunction accumM_tainted

                pure subroutine accumM_addReal_dp(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(inout)      :: self
                        real    (dp), intent(in   )      :: newValue(:,:)
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !local vars
                        real    (dp)    :: delta
                        integer (si)    :: i,j

                        if (.not. allocated(self%mean)) then
                                allocate(self%mean(size(newValue,1),size(newValue,2)))
                                self%mean = 0
                        else
                                if (any(shape(self%mean) /= shape(newValue))) then
                                        self%tainted = .true.
                                        self%mean = 0
                                        return
                                endif
                        endif

                        self%nSamples = self%nSamples + 1

                        do i=1,size(newValue,1)
                                do j=1,size(newValue,2)
                                        delta         = newValue(i,j)  - self%mean(i,j)
                                        self%mean(i,j)  = self%mean(i,j) + delta/self%nSamples
                                enddo
                        enddo

                        return
                endsubroutine accumM_addReal_dp

                pure subroutine accumM_addInt_si(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(inout)      :: self
                        integer (si), intent(in   )      :: newValue(:,:)
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        call self%add(real(newValue,dp))

                        return
                endsubroutine accumM_addInt_si

                pure subroutine accumM_addlogical(self,newValue)
                        use env_kindtypes, only: dp

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(inout)      :: self
                        logical,      intent(in   )      :: newValue(:,:)
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !local vars
                        real    (dp),allocatable :: converted_value(:,:)
                        integer (si)             :: i,j

                        allocate(converted_value(size(newValue,1),size(newValue,2)))

                        do i=1,size(newValue,1)
                                do j=1,size(newValue,2)
                                        if (newValue(i,j)) then
                                                converted_value(i,j) = 1.0_dp
                                        else
                                                converted_value(i,j) = 0.0_dp
                                        endif
                                enddo
                        enddo

                        call self%add(converted_value)

                        return
                endsubroutine accumM_addlogical

                !get number of samples.
                pure function accumM_getN(self)

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        integer (si) :: accumM_getN

                        accumM_getN = self%nSamples

                        return
                endfunction accumM_getN

                pure function accumM_getDim(self)
                        use env_kindtypes, only: si

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        integer (si),allocatable :: accumM_getDim(:)

                        accumM_getDim = shape(self%mean)

                        return
                endfunction accumM_getDim

                pure function accumM_getSum(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp),allocatable :: accumM_getSum(:,:)
                        
                        accumM_getSum = self%mean * real(self%nSamples,dp)

                        return
                endfunction accumM_getSum

                pure function accumM_getMean(self)
                        use env_kindtypes, only: dp

                        implicit none

                        !!! dummy arguments !!!!!!!!!!!!!!!!!!!!!!
                        class(accumM),intent(in   )      :: self
                        !!! end dummy arguments !!!!!!!!!!!!!!!!!!

                        !return val
                        real    (dp),allocatable :: accumM_getMean(:,:)

                        accumM_getMean = self%mean

                        return
                endfunction accumM_getMean

endmodule obj_accum

!program main
!        use env_kindtypes, only: si, dp
!        use obj_accum
!
!        implicit none
!
!        !real (dp)       :: a,b,c,d
!        real (dp)       :: array1(2,2) = reshape([1,2,3,4],[2,2])
!        real (dp)       :: array2(2,2) = reshape([3,2,5,2],[2,2])
!        real (dp)       :: array3(4) = [7,3,2,1]
!        integer (si)    :: array5(4) = [7,3,2,1]
!        integer (si)    :: array6(4) = [2,1,1,8]
!        type(AccumM)  :: accumulator
!
!        call accumulator%reset()
!        call accumulator%add(array1)
!        call accumulator%add(array2)
!
!        print*, accumulator%getN()
!        print*, accumulator%getSum()
!        print*, accumulator%getMean()
!!        print*, accumulator%getVar()
!!        print*, accumulator%getSD()
!endprogram
