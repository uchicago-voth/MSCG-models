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
! Provides definitions of real number loggers.
!

module obj_Rloggers

        use env_kindtypes,       only: si, sp, dp
        use abs_obj_logger_real, only: logger_real

        implicit none

        private

        integer (si),parameter :: MAXBUFFERSIZE = 2**10 !to stop uninit buffer runaway

        public          :: t_rLogger

        !memory efficient real logger, internel representation is sp.
        type, extends(logger_real) :: t_rLogger
                private

                logical                  :: tainted = .false.
                logical                  :: active  = .false.

                integer (si)             :: increment     = 20!length of log buffer
                integer (si)             :: length_buffer = 0 !length of log buffer
                integer (si)             :: buffer_place  = 0 !number of iterations
                real    (sp),allocatable :: buffer(:)         !buffer holding log

                contains
                        private 

                        procedure, public :: initialize   => t_rlogger_initialize
                        procedure, public :: configure    => t_rlogger_configure
                        procedure, public :: reset        => t_rlogger_reset
                        procedure, public :: trim         => t_rlogger_trim

                        procedure         :: addReal_dp   => t_rlogger_addReal_dp
                        procedure         :: addReal_sp   => t_rlogger_addReal_sp

                        procedure         :: sgetLog_dp   => t_rlogger_sgetLog_dp
                        procedure         :: sgetLog_sp   => t_rlogger_sgetLog_sp

                        procedure         :: setLog_dp    => t_rlogger_setLog_dp
                        procedure         :: setLog_sp    => t_rlogger_setLog_sp

                        procedure, public :: getLog_sp    => t_rlogger_getLog_sp

                        procedure, public :: isTainted    => t_rlogger_isTainted
                        procedure, public :: isActive     => t_rlogger_isActive
                        procedure, public :: getLength    => t_rlogger_getLength
        end type t_rLogger

        contains
                pure subroutine t_rlogger_initialize(self,buffer_chunk_size)

                        implicit none

                        class(t_rLogger),          intent(inout) :: self
                        integer (si),    optional, intent(in   ) :: buffer_chunk_size


                        if (present(buffer_chunk_size)) then
                                self%length_buffer = buffer_chunk_size
                        endif

                        if (self%length_buffer < 1 .or. self%length_buffer > MAXBUFFERSIZE) then
                                self%tainted = .true.
                                return
                        endif

                        allocate(self%buffer(self%length_buffer))

                        self%buffer_place = 0
                        self%tainted = .false.
                        self%active  = .true.

                endsubroutine t_rlogger_initialize

                pure subroutine t_rlogger_configure(self,buffer_chunk_size)

                        implicit none

                        class(t_rLogger),          intent(inout) :: self
                        integer (si),              intent(in   ) :: buffer_chunk_size

                        self%length_buffer = buffer_chunk_size
                        self%buffer_place = 0

                endsubroutine t_rlogger_configure

                pure subroutine t_rlogger_reset(self,shallow)

                        implicit none

                        class(t_rLogger),         intent(inout) :: self
                        logical,           optional,intent(in   ) :: shallow

                        logical   :: shallow_

                        !default arg guard
                        if (present(shallow)) then
                                shallow_ = shallow
                        else
                                shallow_ = .true.
                        endif

                        if (.not. shallow_) then
                                !careful deallocation
                                if (allocated(self%buffer)) then
                                        deallocate(self%buffer)
                                endif
                                self%length_buffer = 0
                        endif

                        self%buffer_place  = 0

                        self%tainted = .false.
                        self%active  = .false.

                endsubroutine t_rlogger_reset

                pure subroutine t_rlogger_trim(self)

                        implicit none

                        class(t_rLogger),        intent(inout) :: self

                        real    (sp), allocatable :: tmp_buffer(:)

                        if (.not. allocated(self%buffer)) then
                                self%tainted = .true.
                                return
                        endif

                        call move_alloc(self%buffer,tmp_buffer)

                        allocate(self%buffer(self%buffer_place))

                        self%buffer = tmp_buffer(1:(self%buffer_place))

                        self%length_buffer = self%buffer_place

                endsubroutine t_rlogger_trim

                pure subroutine t_rlogger_extend_buffer(self,amount)

                        implicit none

                        class(t_rLogger),        intent(inout) :: self
                        integer (si),            intent(in   ) :: amount

                        real    (si),allocatable :: tmp_buffer(:)

                        if (.not.allocated(self%buffer)) then
                                self%tainted = .true.
                                return
                        endif

                        call move_alloc(self%buffer,tmp_buffer)

                        self%length_buffer = self%length_buffer + amount

                        allocate(self%buffer(self%length_buffer))

                        self%buffer(1:(self%length_buffer-amount)) = tmp_buffer

                endsubroutine

                pure subroutine t_rlogger_addReal_sp(self,toAdd)

                        implicit none

                        class(t_rLogger),        intent(inout) :: self
                        real    (sp),              intent(in   ) :: toAdd

                        if (.not. allocated(self%buffer)) then
                                self%tainted = .true.
                                return
                        endif

                        if (self%buffer_place == self%length_buffer) then
                                call t_rlogger_extend_buffer(self, self%increment)
                        endif

                        self%buffer_place = self%buffer_place + 1

                        self%buffer(self%buffer_place) = toAdd

                endsubroutine t_rlogger_addReal_sp

                pure subroutine t_rlogger_addReal_dp(self,toAdd)

                        implicit none

                        class(t_rLogger),        intent(inout) :: self
                        real    (dp),              intent(in   ) :: toAdd

                        if (.not. allocated(self%buffer)) then
                                self%tainted = .true.
                                return
                        endif

                        if (self%buffer_place == self%length_buffer) then
                                call t_rlogger_extend_buffer(self, self%increment)
                        endif

                        self%buffer_place = self%buffer_place + 1

                        self%buffer(self%buffer_place) = real(toAdd,sp)

                endsubroutine t_rlogger_addReal_dp

                pure subroutine t_rlogger_sgetLog_sp(self,toSet)

                        implicit none

                        class(t_rLogger),         intent(in   ) :: self
                        real     (sp),allocatable,intent(inout) :: toSet(:)

                        if (self%buffer_place == 0) then
                                !careful deallocation
                                if (allocated(toSet)) then
                                        deallocate(toSet)
                                endif
                                allocate(toSet(0))
                                return
                        endif

                        toSet = self%buffer(1:self%buffer_place)

                endsubroutine t_rlogger_sgetLog_sp

                pure subroutine t_rlogger_sgetLog_dp(self,toSet)

                        implicit none

                        class(t_rLogger),         intent(in   ) :: self
                        real     (dp),allocatable,intent(inout) :: toSet(:)

                        if (self%buffer_place == 0) then
                                !careful deallocation
                                if (allocated(toSet)) then
                                        deallocate(toSet)
                                endif
                                allocate(toSet(0))
                                return
                        endif

                        toSet = real(self%buffer(1:self%buffer_place),dp)

                endsubroutine t_rlogger_sgetLog_dp

                pure function t_rlogger_getLog_sp(self)

                        implicit none

                        class(t_rLogger),       intent(in   ) :: self

                        real     (sp),allocatable             :: t_rlogger_getLog_sp(:)

                        if (self%buffer_place == 0) then
                                allocate(t_rlogger_getLog_sp(0))
                                return
                        endif

                        t_rlogger_getLog_sp = self%buffer(1:self%buffer_place)

                endfunction t_rlogger_getLog_sp

                pure function t_rlogger_isTainted(self)

                        implicit none

                        class(t_rLogger),       intent(in   ) :: self

                        logical         :: t_rlogger_isTainted

                        t_rlogger_isTainted = self%tainted

                endfunction t_rlogger_isTainted

                pure function t_rlogger_isActive(self)

                        implicit none

                        class(t_rLogger),       intent(in   ) :: self

                        logical         :: t_rlogger_isActive

                        t_rlogger_isActive = self%active

                endfunction t_rlogger_isActive

                pure function t_rlogger_getLength(self)

                        implicit none

                        class(t_rLogger),       intent(in   ) :: self

                        integer (si)         :: t_rlogger_getLength

                        t_rlogger_getLength = self%buffer_place

                endfunction t_rlogger_getLength

                pure subroutine t_rlogger_setLog_sp(self,toSet)

                        implicit none

                        class(t_rLogger),       intent(inout) :: self
                        real    (sp),           intent(in   ) :: toSet(:)

                        !if ( self%buffer_place /= 0 .and. self%length_buffer > 0) then
                        !        self%tainted = .true.
                        !        return
                        !endif

                        self%active  = .true.
                        self%tainted = .false.

                        self%buffer = toSet
                        self%buffer_place  = size(toSet)
                        self%length_buffer = size(toSet)

                endsubroutine t_rlogger_setLog_sp

                pure subroutine t_rlogger_setLog_dp(self,toSet)

                        implicit none

                        class(t_rLogger),       intent(inout) :: self
                        real    (dp),           intent(in   ) :: toSet(:)

                        !this happens more often than you would think.
                        !if ( self%buffer_place /= 0 .and. self%length_buffer > 0) then
                        !        self%tainted = .true.
                        !        return
                        !endif

                        self%active  = .true.
                        self%tainted = .false.

                        self%buffer = real(toSet,sp)
                        self%buffer_place = size(toSet)
                        self%length_buffer= size(toSet)

                endsubroutine t_rlogger_setLog_dp
endmodule obj_Rloggers
