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
! Provides an interface to write xyz coordinate files (vmd style).
!
! NOTE:
!    These routines were developed in house, in contrast to the psf/dcd routines.

module IO_xyz

        implicit none

        private

        public                :: write_xyz_coord

        interface write_xyz_coord
                  module procedure openfile_write_xyz_coord, filename_write_xyz_coord,&
                                   openfile_write_xyz_frame, filename_write_xyz_frame
        end interface

        contains
                subroutine openfile_write_xyz_coord(file_handle,coord,comment) 
                        !write positions array to xyz open fileid.

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),     intent(in)          :: coord(:,:,:)
                        integer (si),     intent(in)          :: file_handle

                        character (len=*),intent(in),optional :: comment
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional variable counterparts
                        character (len=:),allocatable :: comment_

                        !local variables
                        integer (si)    :: numSites, spaceDim, numSteps
                        integer (si)    :: step, i, j

                        if (present(comment)) then
                                comment_ = comment
                        else
                                comment_ = '("CG coordinates from minimized charge/EDCG procedure")'
                        endif
                        numSites = size(coord,1)
                        spaceDim = size(coord,2)
                        numSteps = size(coord,3)

                        do step=1,numSteps
                                write(file_handle,'(i5)') numSites
                                write(file_handle, comment_)
                                do i=1,numSites
                                        write(file_handle,'("C  ", 3f8.3)') &
                                             (coord(i,j,step),j=1,spaceDim)
                                enddo
                        enddo

                        return

                endsubroutine openfile_write_xyz_coord

                subroutine openfile_write_xyz_frame(file_handle,frame,comment) 
                        !write positions array to xyz open fileid.

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),     intent(in)          :: frame(:,:)
                        integer (si),     intent(in)          :: file_handle

                        character (len=*),intent(in),optional :: comment
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional variable counterparts
                        character (len=:),allocatable :: comment_

                        !local variables
                        integer (si)    :: numSites, spaceDim
                        integer (si)    :: i, j

                        if (present(comment)) then
                                comment_ = comment
                        else
                                comment_ = '("CG coordinates from minimized charge/EDCG procedure")'
                        endif
                        numSites = size(frame,1)
                        spaceDim = size(frame,2)

                        write(file_handle,'(i5)') numSites
                        write(file_handle, comment_)
                        do i=1,numSites
                                write(file_handle,'("C  ", 3f8.3)') &
                                     (frame(i,j),j=1,spaceDim)
                        enddo

                        return

                endsubroutine openfile_write_xyz_frame

                subroutine filename_write_xyz_coord(filename,coord,comment) 
                        !write positions array to xyz file given a filename.

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),     intent(in)          :: coord(:,:,:)
                        character (len=*),intent(in),optional :: filename

                        character (len=*),intent(in),optional :: comment
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional variable counterparts
                        character (len=:),allocatable   :: comment_

                        !local variables
                        integer (si)                    :: numSites, spaceDim, numSteps
                        integer (si)                    :: step, i, j
                        integer (si)                    :: file_handle, status

                        if (len(filename) == 0) then
                                return
                        endif

                        if (present(comment)) then
                                comment_ = comment
                        else
                                comment_ = '("CG coordinates from minimized charge/EDCG procedure")'
                        endif

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "xyz file couldn't be written."
                                return
                        endif

                        numSites = size(coord,1)
                        spaceDim = size(coord,2)
                        numSteps = size(coord,3)

                        do step=1,numSteps
                                write(file_handle,'(i5)') numSites
                                write(file_handle, comment_)
                                do i=1,numSites
                                        write(file_handle,'("C  ", 3f8.3)') &
                                             (coord(i,j,step),j=1,spaceDim)
                                enddo
                        enddo

                        return

                endsubroutine filename_write_xyz_coord

                subroutine filename_write_xyz_frame(filename,frame,comment) 
                        !write positions array to xyz file given a filename.

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),     intent(in)          :: frame(:,:)
                        character (len=*),intent(in),optional :: filename

                        character (len=*),intent(in),optional :: comment
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional variable counterparts
                        character (len=:),allocatable   :: comment_

                        !local variables
                        integer (si)                    :: numSites, spaceDim
                        integer (si)                    :: i, j
                        integer (si)                    :: file_handle, status

                        if (len(filename) == 0) then
                                return
                        endif

                        if (present(comment)) then
                                comment_ = comment
                        else
                                comment_ = '("CG coordinates from minimized charge/EDCG procedure")'
                        endif

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "xyz file couldn't be written."
                                return
                        endif

                        numSites = size(frame,1)
                        spaceDim = size(frame,2)

                        write(file_handle,'(i5)') numSites
                        write(file_handle, comment_)
                        do i=1,numSites
                                write(file_handle,'("C  ", 3f8.3)') &
                                     (frame(i,j),j=1,spaceDim)
                        enddo

                        return

                endsubroutine filename_write_xyz_frame

endmodule IO_xyz
