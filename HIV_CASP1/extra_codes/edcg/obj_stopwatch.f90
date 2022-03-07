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
! Provides an omp based performance stopwatch.

module obj_stopwatch

        use env_kindtypes, only: si, dp, dp_x

        implicit none

        private

        public          :: stopwatch

        !Basic trajectory object. Unlabeled.
        type                    :: stopwatch

                private

                real (dp)              :: start_time
                real (dp)              :: finish_time

                contains
                        procedure,public :: start      => start_stopwatch
                        procedure,public :: finish     => finish_stopwatch
                        procedure,public :: report     => report_stopwatch

        end type stopwatch

        contains
                subroutine start_stopwatch(self)

                        implicit none

                        interface
                                function omp_get_wtime()
                                        import dp_x
                                        real  (dp_x) :: omp_get_wtime
                                endfunction
                        end interface

                        class   (stopwatch),intent(inout)   :: self

                        self%start_time = omp_get_wtime()

                endsubroutine start_stopwatch

                subroutine finish_stopwatch(self)

                        implicit none

                        interface
                                function omp_get_wtime()
                                        import dp_x
                                        real (dp_x) :: omp_get_wtime
                                endfunction
                        end interface

                        class   (stopwatch),intent(inout)   :: self

                        self%finish_time = omp_get_wtime()

                endsubroutine finish_stopwatch

                subroutine report_stopwatch(self,file_pointer,prefix)

                        implicit none

                        class(stopwatch),intent(inout)          :: self
                        integer (si),    intent(in   ),optional :: file_pointer
                        character(*),    intent(in   ),optional :: prefix

                        real (dp)              :: total_time
                        character(:), allocatable :: style_string

                        total_time = self%finish_time - self%start_time

                        if (present(prefix)) then
                                style_string = '("'//prefix//'Total time elapsed:",f16.3)'
                        else
                                style_string = '("Total time elapsed:",f16.3)'
                        endif

                        if (present(file_pointer)) then
                                write(file_pointer,style_string) total_time
                        else
                                write(*,style_string) total_time
                        endif

                endsubroutine report_stopwatch
endmodule obj_stopwatch
