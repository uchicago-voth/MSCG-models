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
! This module contains general conversion utilities between types using virtual files.
! Only generic interfaces should be public.
!
module core_convert

        implicit none

        private

        character (len=1),dimension(*),parameter :: neg_digits_set &
                = ['-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' ]
        character (len=1),parameter               :: dot &
                = '.'

        public itoa, ftoa, atoi, atof, is_real, is_integer

        !integer to string conversion.
        interface itoa
                  module procedure integer_si_to_string
        end interface

        !float to integer conversion
        interface ftoa
                  module procedure real_dp_to_string
        end interface

        !string to integer conversion
        interface atoi
                  module procedure string_to_integer_si
        end interface

        !string to float conversion
        interface atof
                  module procedure string_to_real_dp
        end interface

        contains
                pure function integer_si_to_string(toConvert) result(toReturn)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si), intent(in   )     :: toConvert
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        character(:),allocatable :: toReturn

                        !temporary hold variable
                        character(range(toConvert)+2) :: tmp

                        write(tmp,'(i0)') toConvert

                        toReturn = trim(tmp)

                endfunction integer_si_to_string

                pure function real_dp_to_string(toConvert) result(toReturn)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )     :: toConvert
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        character(:),allocatable :: toReturn

                        !temporary hold variable
                        character(range(toConvert)+2) :: tmp

                        write(tmp,*) toConvert

                        toReturn = trim(tmp)

                endfunction real_dp_to_string

                elemental function string_to_integer_si(toConvert) result(toReturn)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character(*), intent(in   )     :: toConvert
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si)                    :: toReturn

                        read(toConvert,'(i10)')  toReturn

                endfunction string_to_integer_si

                elemental function string_to_real_dp(toConvert) result(toReturn)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character(*), intent(in   )     :: toConvert
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)                    :: toReturn

                        read(toConvert,*)  toReturn

                endfunction string_to_real_dp

                elemental function is_integer(toCheck)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character(*), intent(in   )     :: toCheck
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        logical         :: is_integer

                        integer (si)    :: iter

                        is_integer = .true.

                        do iter=1,len(trim(toCheck))
                                if (.not. string_in_set(toCheck(iter:iter),neg_digits_set)) then
                                        is_integer = .false.
                                        return
                                endif
                        enddo

                endfunction is_integer

                elemental function is_real(toCheck)

                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character(*), intent(in   )     :: toCheck
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical         :: is_real, dot_found

                        integer (si)    :: iter

                        is_real   = .true.
                        dot_found = .false.

                        outer: do iter=1,len(trim(toCheck))
                                if (.not. string_in_set(toCheck(iter:iter),neg_digits_set)) then
                                        if (toCheck(iter:iter) == dot) then
                                                if (dot_found) then
                                                        is_real = .false.
                                                        return
                                                else
                                                        dot_found = .true.
                                                        cycle outer
                                                endif
                                        endif
                                        is_real = .false.
                                        return
                                endif
                        enddo outer

                endfunction is_real

                pure function string_in_set(query,set)

                        use env_kindtypes,      only: si

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character(len=*),              intent(in   ) :: query
                        character(len=*),dimension(:), intent(in   ) :: set
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical                 :: string_in_set
                        integer (si)            :: iter

                        string_in_set = .false.

                        do iter=1,size(set)
                                !print*, "query" , query
                                !print*, "set(iter)" , set(iter)
                                !print*, "result" , (query == set(iter))
                                if (query == set(iter)) then
                                        string_in_set = .true.
                                        return
                                endif
                        enddo

                endfunction string_in_set


endmodule core_convert
