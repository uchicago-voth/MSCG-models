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
! Provides _abstract_ definition of a logger, an object which stores a growing list.

module abs_obj_logger_real

        use env_kindtypes, only: si, sp, dp

        implicit none

        private

        public          :: logger_real

        type, abstract                 :: logger_real

                private

                contains
                        private

                        procedure(abs_init),         deferred, public :: initialize
                        !configures options of the logger
                        procedure(abs_configure),    deferred, public :: configure
                        procedure(abs_self_reset),   deferred, public :: reset
                        procedure(abs_self_inout_p), deferred, public :: trim

                        procedure(abs_add_real_dp),  deferred         :: addReal_dp
                        procedure(abs_add_real_sp),  deferred         :: addReal_sp
                        generic,                               public :: add    => &
                                        addReal_dp, addReal_sp

                        procedure(abs_sgetArray_dp), deferred         :: sgetLog_dp
                        procedure(abs_sgetArray_sp), deferred         :: sgetLog_sp
                        generic,                               public :: sgetLog => &
                                        sgetLog_dp , sgetLog_sp

                        procedure(abs_setArray_dp),  deferred         :: setLog_dp
                        procedure(abs_setArray_sp),  deferred         :: setLog_sp
                        generic,                               public :: setLog => &
                                        setLog_dp , setLog_sp

                        procedure(abs_getArray_sp),  deferred, public :: getLog_sp

                        !tainted logs are have untrustable content.
                        procedure(abs_tainted),      deferred, public :: isTainted
                        procedure(abs_tainted),      deferred, public :: isActive

                        procedure(abs_get_si),       deferred, public :: getLength
        end type logger_real

        abstract interface
            subroutine abs_self_inout_p(self)
                import :: logger_real
                class(logger_real),      intent(inout) :: self
            endsubroutine abs_self_inout_p
        end interface

        abstract interface
            subroutine abs_self_reset(self,shallow)
                import :: logger_real
                class(logger_real),         intent(inout) :: self
                logical,           optional,intent(in   ) :: shallow
            endsubroutine abs_self_reset
        end interface

        abstract interface
            subroutine abs_init(self,buffer_chunk_size)
                import :: logger_real,si
                class(logger_real),          intent(inout) :: self
                integer (si),      optional, intent(in   ) :: buffer_chunk_size
            endsubroutine abs_init
        end interface

        abstract interface
            subroutine abs_configure(self,buffer_chunk_size)
                import :: logger_real,si
                class(logger_real),          intent(inout) :: self
                integer (si),                intent(in   ) :: buffer_chunk_size
            endsubroutine abs_configure
        end interface


        abstract interface
            pure function abs_getArray_sp(self)
                import :: logger_real,sp
                class(logger_real),        intent(in   ) :: self

                real    (sp), allocatable                :: abs_getArray_sp(:)
            endfunction abs_getArray_sp
        end interface

        abstract interface
            pure subroutine abs_sgetArray_dp(self,toSet)
                import :: logger_real,dp
                class(logger_real),        intent(in   ) :: self
                real    (dp), allocatable, intent(inout):: toSet(:)
            endsubroutine abs_sgetArray_dp
        end interface

        abstract interface
            pure subroutine abs_sgetArray_sp(self,toSet)
                import :: logger_real,sp
                class(logger_real),        intent(in   ) :: self
                real     (sp), allocatable, intent(inout) :: toSet(:)
            endsubroutine abs_sgetArray_sp
        end interface

        abstract interface
            pure subroutine abs_add_real_dp(self,toAdd)
                import :: logger_real,dp
                class(logger_real), intent(inout) :: self
                real    (dp),  intent(in   ) :: toAdd
            endsubroutine abs_add_real_dp
        end interface

        abstract interface
            pure subroutine abs_add_real_sp(self,toAdd)
                import :: logger_real,sp
                class(logger_real), intent(inout) :: self
                real    (sp),  intent(in   ) :: toAdd
            endsubroutine abs_add_real_sp
        end interface

        abstract interface
            pure function abs_tainted(self)
                import :: logger_real,sp
                class(logger_real), intent(in   ) :: self
                logical abs_tainted
            endfunction abs_tainted
        end interface

        abstract interface
            pure function abs_get_si(self)
                import :: logger_real,si
                class(logger_real), intent(in   ) :: self
                integer (si) :: abs_get_si
            endfunction abs_get_si
        end interface

        abstract interface
            pure subroutine abs_setArray_sp(self,toSet)
                import :: logger_real,sp
                class(logger_real), intent(inout) :: self
                real    (sp),       intent(in   ) :: toSet(:)
            endsubroutine abs_setArray_sp
        end interface

        abstract interface
            pure subroutine abs_setArray_dp(self,toSet)
                import :: logger_real,dp
                class(logger_real), intent(inout) :: self
                real    (dp),       intent(in   ) :: toSet(:)
            endsubroutine abs_setArray_dp
        end interface

endmodule abs_obj_logger_real
