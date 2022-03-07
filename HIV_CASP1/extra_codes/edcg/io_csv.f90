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
! Provides IO routines for csv files, implemented as a generic interface.
!
! NOTE:
!   Currently, only write is implemented. This is done for 1-3 dim arrays.
!   Routines work with either file ids or filenames, real dp's or int si's.
!

module IO_csv

        implicit none

        private

        public                :: write_csv

        interface write_csv
                  module procedure openfile_write_csv_case_array_intsi,&
                                   openfile_write_csv_case_array_realdp,&
                                   openfile_write_csv_2array_intsi,&
                                   openfile_write_csv_2array_realdp,&
                                   openfile_write_csv_3array_intsi,&
                                   openfile_write_csv_3array_realdp,&
                                   filename_write_csv_case_array_intsi,&
                                   filename_write_csv_case_array_realdp,&
                                   filename_write_csv_2array_intsi,&
                                   filename_write_csv_2array_realdp,&
                                   filename_write_csv_3array_intsi,&
                                   filename_write_csv_3array_realdp
        end interface

        !tbd
        !interface write_csv
        !          module procedure
        !end interface

        contains
                !write 3d array as a csv with an initial index identifying a
                !sequence of 2d slices.
                subroutine openfile_write_csv_3array_realdp(file_handle,array)
                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: file_handle
                        real    (dp),intent(in   )          :: array(:,:,:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                :: i

                        do i=1,size(array,3)
                                call openfile_write_csv_2array_realdp(file_handle,array(:,:,i),i)
                        enddo

                        return
                endsubroutine openfile_write_csv_3array_realdp

                !write 3d array as a csv with an initial index identifying a
                !sequence of 2d slices.
                subroutine openfile_write_csv_3array_intsi(file_handle,array)
                        use env_kindtypes,      only: si

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: file_handle
                        integer (si),intent(in   )          :: array(:,:,:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                :: i

                        do i=1,(size(array,1))
                                call openfile_write_csv_2array_intsi(file_handle,array(i,:,:),i)
                        enddo

                        return
                endsubroutine openfile_write_csv_3array_intsi

                !write 2d array. with optional INT prefix.
                subroutine openfile_write_csv_2array_realdp(file_handle,array,prefix)
                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: file_handle
                        real    (dp),intent(in   )          :: array(:,:)
                        integer (si),intent(in   ),optional :: prefix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                :: i

                        do i=1,(size(array,1))
                                if (present(prefix)) then
                                        call openfile_write_csv_case_array_realdp(file_handle,array(i,:),prefix)
                                else
                                        call openfile_write_csv_case_array_realdp(file_handle,array(i,:))
                                endif
                        enddo

                        return
                endsubroutine openfile_write_csv_2array_realdp

                !write 2d array. with optional INT prefix.
                subroutine openfile_write_csv_2array_intsi(file_handle,array,prefix)
                        use env_kindtypes,      only: si

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: file_handle
                        integer (si),intent(in   )          :: array(:,:)
                        integer (si),intent(in   ),optional :: prefix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                :: i

                        do i=1,(size(array,1))
                                if (present(prefix)) then
                                        call openfile_write_csv_case_array_intsi(file_handle,array(i,:),prefix)
                                else
                                        call openfile_write_csv_case_array_intsi(file_handle,array(i,:))
                                endif
                        enddo

                        return
                endsubroutine openfile_write_csv_2array_intsi

                !write 1d array. with optional INT prefix.
                !There are two output styles depending on size.
                !(this is more important for integers to avoid whitespace.)
                subroutine openfile_write_csv_case_array_realdp(file_handle,case_array,prefix)
                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: file_handle
                        real    (dp),intent(in   )          :: case_array(:)
                        integer (si),intent(in   ),optional :: prefix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                :: i
                        integer (si),     parameter :: size_check   = 0
                        character (len=*),parameter :: real_short_format = "(f10.5)"  !note size occupied by signs
                        character (len=*),parameter :: real_long_format  = "(f35.5)"
                        character (len=*),parameter :: int_short_format = "(i7)"  !note size occupied by signs
                        character (len=*),parameter :: int_long_format  = "(i11)"

                        !apply the prefix.
                        if (present(prefix)) then
                                if (abs(prefix) <= size_check) then
                                        write(file_handle,int_short_format,advance="no") prefix
                                else
                                        write(file_handle,int_long_format, advance="no") prefix
                                endif
                                write(file_handle,'(a)',advance="no") ","
                        endif

                        !write all with separating commas except the last value.
                        do i=1,(size(case_array,1)-1)
                                if (abs(case_array(i)) <= size_check) then
                                        write(file_handle,real_short_format,advance="no") case_array(i)
                                else
                                        write(file_handle,real_long_format, advance="no") case_array(i)
                                endif
                                write(file_handle,'(a)',advance="no") ","
                        enddo

                        !write the last value. Done out of loop because of the lack of tailing comma.
                        if (abs(case_array(size(case_array,1))) <= size_check) then
                                write(file_handle,real_short_format,advance="yes") case_array(size(case_array,1))
                        else
                                write(file_handle,real_long_format, advance="yes") case_array(size(case_array,1))
                        endif

                        return
                endsubroutine openfile_write_csv_case_array_realdp

                !write 1d array. with optional INT prefix.
                !There are two output styles depending on size.
                !(this is more important for integers to avoid whitespace.)
                subroutine openfile_write_csv_case_array_intsi(file_handle,case_array,prefix)
                        use env_kindtypes,      only: si

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: file_handle
                        integer (si),intent(in   )          :: case_array(:)
                        integer (si),intent(in   ),optional :: prefix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                :: i
                        integer (si),     parameter :: size_check   = 999999
                        character (len=*),parameter :: short_format = "(i7)"  !note size occupied by signs
                        character (len=*),parameter :: long_format  = "(i11)"

                        !apply the prefix.
                        if (present(prefix)) then
                                if (abs(prefix) <= size_check) then
                                        write(file_handle,short_format,advance="no") prefix
                                else
                                        write(file_handle,long_format, advance="no") prefix
                                endif
                                write(file_handle,'(a)',advance="no") ","
                        endif

                        !write all with separating commas except the last value.
                        do i=1,(size(case_array,1)-1)
                                if (abs(case_array(i)) <= size_check) then
                                        write(file_handle,short_format,advance="no") case_array(i)
                                else
                                        write(file_handle,long_format, advance="no") case_array(i)
                                endif
                                write(file_handle,'(a)',advance="no") ","
                        enddo

                        !write the last value. Done out of loop because of the lack of tailing comma.
                        if (abs(case_array(i)) <= size_check) then
                                write(file_handle,short_format,advance="yes") case_array(size(case_array,1))
                        else
                                write(file_handle,long_format, advance="yes") case_array(size(case_array,1))
                        endif

                        return
                endsubroutine openfile_write_csv_case_array_intsi

                !!! One shot routines which write with a filename instead of an ID. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !They just point towards their file id based counterparts after opening the file, and then
                !close it.

                subroutine filename_write_csv_3array_realdp(filename,array)
                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character (len=*),intent(in)        :: filename
                        real    (dp),     intent(in   )     :: array(:,:,:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        integer (si)    :: file_handle, status

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "csv file couldn't be written."
                                return
                        endif

                        call openfile_write_csv_3array_realdp(file_handle,array)

                        close(file_handle)
                        return
                endsubroutine filename_write_csv_3array_realdp

                subroutine filename_write_csv_3array_intsi(filename,array)
                        use env_kindtypes,      only: si

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character (len=*),intent(in)        :: filename
                        integer (si),intent(in   )          :: array(:,:,:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        integer (si)    :: file_handle, status

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "csv file couldn't be written."
                                return
                        endif

                        call openfile_write_csv_3array_intsi(file_handle,array)

                        close(file_handle)
                        return
                endsubroutine filename_write_csv_3array_intsi

                subroutine filename_write_csv_2array_realdp(filename,array,prefix)
                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character (len=*),intent(in)        :: filename
                        real    (dp),intent(in   )          :: array(:,:)
                        integer (si),intent(in   ),optional :: prefix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        integer (si)    :: file_handle, status

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "csv file couldn't be written."
                                return
                        endif

                        if (present(prefix)) then
                                call openfile_write_csv_2array_realdp(file_handle,array,prefix)
                        else
                                call openfile_write_csv_2array_realdp(file_handle,array)
                        endif

                        close(file_handle)
                        return
                endsubroutine filename_write_csv_2array_realdp

                subroutine filename_write_csv_2array_intsi(filename,array,prefix)
                        use env_kindtypes,      only: si

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character (len=*),intent(in)        :: filename
                        integer (si),intent(in   )          :: array(:,:)
                        integer (si),intent(in   ),optional :: prefix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        integer (si)    :: file_handle, status

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "csv file couldn't be written."
                                return
                        endif

                        if (present(prefix)) then
                                call openfile_write_csv_2array_intsi(file_handle,array,prefix)
                        else
                                call openfile_write_csv_2array_intsi(file_handle,array)
                        endif

                        close(file_handle)
                        return
                endsubroutine filename_write_csv_2array_intsi

                subroutine filename_write_csv_case_array_realdp(filename,case_array,prefix)
                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character (len=*),intent(in)        :: filename
                        real    (dp),intent(in   )          :: case_array(:)
                        integer (si),intent(in   ),optional :: prefix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        integer (si)    :: file_handle, status

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "csv file couldn't be written."
                                return
                        endif

                        if (present(prefix)) then
                                call openfile_write_csv_case_array_realdp(file_handle,case_array,prefix)
                        else
                                call openfile_write_csv_case_array_realdp(file_handle,case_array)
                        endif

                        close(file_handle)
                        return
                endsubroutine filename_write_csv_case_array_realdp

                subroutine filename_write_csv_case_array_intsi(filename,case_array,prefix)
                        use env_kindtypes,      only: si

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character (len=*),intent(in)        :: filename
                        integer (si),intent(in   )          :: case_array(:)
                        integer (si),intent(in   ),optional :: prefix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        integer (si)    :: file_handle, status

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "csv file couldn't be written."
                                return
                        endif

                        if (present(prefix)) then
                                call openfile_write_csv_case_array_intsi(file_handle,case_array,prefix)
                        else
                                call openfile_write_csv_case_array_intsi(file_handle,case_array)
                        endif

                        close(file_handle)
                        return
                endsubroutine filename_write_csv_case_array_intsi

endmodule IO_csv

!program main
!        use env_kindtypes, only: si, dp
!        use IO_csv
!
!        implicit none
!
!        real (dp)       :: array3(2,2,2) = reshape([1,2,3,4,5,6,7,8],[2,2,2])
!        real (dp)       :: array2(2,4)   = reshape([1,2,3,4,5,6,7,8],[2,4])
!        real (dp)       :: array(8)      =         [1,2,3,4,5,6,7,8]
!        integer (si)        :: file_handle, status
!
!        print*, shape(array2)
!        print*, array3(:,1,1)
!        print*, array3(1,:,1)
!        print*, array3(1,1,:)
!        call write_csv("test.csv",array3)
!endprogram
!
