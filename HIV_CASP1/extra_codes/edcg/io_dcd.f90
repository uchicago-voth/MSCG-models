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
! Provides routines for reading DCD files.
!
! NOTE:
!    It's not clear where these routines came from (at all). Some slight modernification 
!    was done, but these have mostly been left to their own devices.

module IO_dcd

        implicit none

        private

        public          :: read_dcd_header, read_dcd_step

        contains

                subroutine read_dcd_header(inputTraj, nAtomTotal, nSteps, filePointer)

                        use env_kindtypes,   only: sp_x, si
                        use env_endianess, only: nativeEndianess

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        character (*), intent(in)      :: inputTraj
                        integer (si),intent(out)            :: nAtomTotal
                        integer (si),intent(out)            :: nSteps
                        integer (si),intent(out)            :: filePointer !This is allocated 
                                                                                !here.

                        !!!!! End dummy arguments  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                        :: junk, i 
                        character (LEN=4)                        :: hdr
                        integer (si)                        :: startTimeStep
                        integer (si)                        :: delta_step
                        integer (si)                        :: lastTimeStep
                        integer (si)                        :: nFreeAtoms
                        real    (sp_x)                      :: delta_t
                        integer (si)                        :: nTitle

                        character (LEN=80) :: title

                        open(newunit=filePointer,file=inputTraj,access='stream')

                        !Read the initial integer (should be 84)
                        read(filePointer) junk
                        if (junk == 84) then
                                nativeEndianess=.true.
                                print*, 'INFO: (dcd file) file has native endianess'
                        else
                                nativeEndianess=.false.
                                print*, 'INFO: (dcd file) file has not native endianess'
                        endif
                        if (nativeEndianess.eqv..true.) then
                              !Read the header (should be 'CORD' char*4)
                                read(filePointer) hdr
                                !Read the number of steps (integer)
                                read(filePointer) nSteps
                                print*, "INFO: (dcd file) Number of steps in trajectory file:", nSteps
                                !Read the first time step (integer)
                                read(filePointer) startTimeStep
                                !Read the steps between trajectory writes
                                read(filePointer) delta_step
                                !Read the last time step (integer)
                                read(filePointer) lastTimeStep
                                !skip 4 blank integers
                                do i=1,4
                                   read(filePointer) junk
                                enddo
                                !read the number of free atoms (integer)
                                read(filePointer) nFreeAtoms
                                !read the time between steps (real single precision)
                                read(filePointer) delta_t
                                !skip some blank or uniteresting integers
                                do i=1,12
                                   read(filePointer) junk
                                enddo
                                !read the number of lines for the title (each line is 80 chars long)
                                read(filePointer) nTitle
                                !read the title lines
                                do i=1,nTitle
                                   read(filePointer) title
                                enddo
                                !read two more integers
                                read(filePointer) junk
                                read(filePointer) junk
                                !read the number of atoms
                                read(filePointer) nAtomTotal
                                print*, 'INFO: (dcd file) Number of atoms in trajectory file:', nAtomTotal
                                !read another integer
                                read(filePointer) junk
                        else
                                !Read the header (should be 'CORD' char*4)
                                read(filePointer) hdr
                                !Read the number of steps (integer)
                                read(filePointer) nSteps
                                nSteps = int_swap_endianess(nSteps)
                                print*, "INFO: (dcd file) Number of steps in trajectory file:", nSteps
                                !Read the first time step (integer)
                                read(filePointer) startTimeStep
                                startTimeStep = int_swap_endianess(StartTimeStep)
                                !Read the steps between trajectory writes
                                read(filePointer) delta_step
                                delta_step = int_swap_endianess(delta_step)
                                !Read the last time step (integer)
                                read(filePointer) lastTimeStep
                                lastTimeStep = int_swap_endianess(lastTimeStep)
                                !skip 4 blank integers
                                do i=1,4
                                   read(filePointer) junk
                                   junk = int_swap_endianess(junk)
                                enddo
                                !read the number of free atoms (integer)
                                read(filePointer) nFreeAtoms
                                nFreeAtoms = int_swap_endianess(nFreeAtoms)
                                !read the time between steps (real single precision)
                                read(filePointer) delta_t
                                delta_t = real_swap_endianess(delta_t)
                                !skip some blank or uniteresting integers
                                do i=1,12
                                   read(filePointer) junk
                                   junk = int_swap_endianess(junk)
                                enddo
                                !read the number of lines for the title (each line is 80 chars long)
                                read(filePointer) nTitle
                                nTitle = int_swap_endianess(nTitle)
                                !read the title lines
                                do i=1,nTitle
                                   read(filePointer) title
                                enddo
                                !read two uninteresting integers
                                read(filePointer) junk
                                junk=int_swap_endianess(junk)
                                read(filePointer) junk
                                junk=int_swap_endianess(junk)
                                !read the number of atoms
                                read(filePointer) nAtomTotal
                                nAtomTotal = int_swap_endianess(nAtomTotal)
                                print*, 'INFO: (dcd file) Number of atoms in trajectory file:', nAtomTotal
                                !read another integer
                                read(filePointer) junk
                                junk=int_swap_endianess(junk)
                        endif


                endsubroutine read_dcd_header

!                !subroutine to read the positions of one step in a NAMD dcd file
!                subroutine get_selected_dcd_coord(coord, nSelectedAtoms, selectedAtoms,&
!                                                        nAtomsTotal, axis, filePointer)
!
!                        use env_kindtypes, only: dp, si, sp_x, dp_x
!                        use endianess
!
!                        implicit none
!
!                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                        integer (si),  intent(in)  ::  nAtomsTotal
!                        integer (si),  intent(in)  ::  nSelectedAtoms
!                        integer (si),  intent(in)  ::  selectedAtoms(nSelectedAtoms)
!                        real    (dp),  intent(out) ::  coord(3,nSelectedAtoms)
!                        integer (si),  intent(in)  ::  filePointer
!                        real    (dp),  intent(out) ::  axis(3)
!
!                        !!!!! End dummy arguments  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                        integer (si)      :: i, j
!                        real    (sp_x)    :: tempPos
!                        integer (si)      :: junk
!                        integer (si)      :: counter
!                        real    (dp_x)    :: junkd
!
!                        if (nativeEndianess.eqv..true.) then
!                                read(filePointer) junk
!                                read(filePointer) axis(1)
!                                read(filePointer) junkd
!                                read(filePointer) axis(2)
!                                read(filePointer) junkd
!                                read(filePointer) junkd
!                                read(filePointer) axis(3)
!                                read(filePointer) junk
!                                do j=1,3
!                                      !nAtoms*4
!                                      counter = 1
!                                      read(filePointer) junk
!                                      do i=1, nAtomsTotal
!                                            read(filePointer) tempPos
!                                            if(i==selectedAtoms(counter)) then
!                                                   coord(j,counter) = real(tempPos,dp)
!                                                   counter          = counter + 1
!                                            endif
!                                      enddo
!                                      !nAtoms*4
!                                      read(filePointer) junk
!                                enddo
!                        else
!                                read(filePointer) junk
!                                junk = int_swap_endianess(junk)
!                                read(filePointer) axis(1)
!                                axis(1) = dble_swap_endianess(axis(1))
!                                read(filePointer) junkd
!                                junkd = dble_swap_endianess(junkd)
!                                read(filePointer) axis(2)
!                                axis(2) = dble_swap_endianess(axis(2))
!                                read(filePointer) junkd
!                                junkd = dble_swap_endianess(junkd)
!                                read(filePointer) junkd
!                                junkd = dble_swap_endianess(junkd)
!                                read(filePointer) axis(3)
!                                axis(3) = dble_swap_endianess(axis(3))
!                                read(filePointer) junk
!                                junk = int_swap_endianess(junk)
!                                do j=1,3
!                                   !nAtoms*4
!                                   counter = 1
!                                   read(filePointer) junk
!                                   junk = int_swap_endianess(junk)
!                                   do i=1, nAtomsTotal
!                                      read(filePointer) tempPos
!                                      if(i==selectedAtoms(counter)) then
!                                             tempPos           = real_swap_endianess(tempPos)
!                                             coord(j,counter)  = real(tempPos,dp)
!                                             counter           = counter + 1
!                                      endif
!                                   enddo
!                                   !nAtoms*4
!                                   read(filePointer) junk
!                                   junk = int_swap_endianess(junk)
!                                enddo
!                        endif
!
!                endsubroutine get_selected_dcd_coord


                !su     broutine to read the positions of one step in a NAMD dcd file
                subroutine read_dcd_step(coord,nAtoms, filePointer)

                        use env_kindtypes, only: dp, si, sp_x, dp_x
                        use env_endianess

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si),intent(in)    :: nAtoms
                        real    (dp),intent(out)   :: coord(nAtoms,3)
                        integer (si),intent(in)    :: filePointer

                        !!!!! End dummy arguments  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)      :: i, j, counter, junk
                        real    (sp_x)    :: tempPos
                        real    (dp_x)    :: axis(3), junkd

                        if (nativeEndianess.eqv..true.) then
                             read(filePointer) junk
                             read(filePointer) axis(1)
                             read(filePointer) junkd
                             read(filePointer) axis(2)
                             read(filePointer) junkd
                             read(filePointer) junkd
                             read(filePointer) axis(3)
                             read(filePointer) junk
                             do j=1,3
                                !nAtoms*4
                                counter = 1
                                read(filePointer) junk
                                do i=1, nAtoms
                                   read(filePointer) tempPos
                                   coord(i,j) = real(tempPos,dp)
                                enddo
                                !nAtoms*4
                                read(filePointer) junk
                             enddo
                        else
                             read(filePointer) junk
                             junk = int_swap_endianess(junk)
                             read(filePointer) axis(1)
                             axis(1) = dble_swap_endianess(axis(1))
                             read(filePointer) junkd
                             junkd = dble_swap_endianess(junkd)
                             read(filePointer) axis(2)
                             axis(2) = dble_swap_endianess(axis(2))
                             read(filePointer) junkd
                             junkd = dble_swap_endianess(junkd)
                             read(filePointer) junkd
                             junkd = dble_swap_endianess(junkd)
                             read(filePointer) axis(3)
                             axis(3) = dble_swap_endianess(axis(3))
                             read(filePointer) junk
                             junk = int_swap_endianess(junk)
                             do j=1,3
                                !nAtoms*4
                                counter = 1
                                read(filePointer) junk
                                junk = int_swap_endianess(junk)
                                do i=1, nAtoms
                                   read(filePointer) tempPos
                                   tempPos     = real_swap_endianess(tempPos)
                                   coord(i,j)  = real(tempPos,dp)
                                enddo
                                !nAtoms*4
                                read(filePointer) junk
                                junk = int_swap_endianess(junk)
                             enddo
                        endif

                endsubroutine read_dcd_step

                !!!!  Endianess Routines

                pure integer (si) function int_swap_endianess(intIn)

                        use env_kindtypes, only: si, singleb, c_si

                        implicit none

                        integer (si),intent(in)  ::   intIn

                        integer (si),    parameter             :: nBytes=c_si
                        integer (singleb), dimension(nBytes)   :: intArray
                        integer (singleb), dimension(nBytes)   :: intReversed

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! Transfer 64 bits of realIn to generic 64 bit INTEGER space:
                        intArray = transfer( intIn, intArray )
                        ! Transfer reversed order bytes to 64 bit REAL space
                        intReversed = intArray(nBytes:1:-1)
                        int_swap_endianess = transfer( intReversed, int_swap_endianess )

                endfunction int_swap_endianess

                pure real (sp) function real_swap_endianess(realIn)

                        use env_kindtypes, only: sp, si, singleb, c_sp

                        implicit none

                        real (sp), intent(in)            :: realIn

                        integer (si), parameter               :: nBytes=c_sp
                        integer (singleb),dimension(nBytes)   :: intIn
                        integer (singleb),dimension(nBytes)   :: intReversed

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! Transfer 64 bits of realIn to generic 64 bit INTEGER space:
                        intIn = transfer( realIn, intIn )
                        ! Transfer reversed order bytes to 64 bit REAL space
                        intReversed = intIn(nBytes:1:-1)
                        real_swap_endianess = transfer( intReversed, real_swap_endianess )

                endfunction real_swap_endianess

                pure real (dp) function dble_swap_endianess(dbleIn)

                        use env_kindtypes, only: dp, si, singleb, c_dp

                        implicit none

                        real (dp), intent(in)               ::  dbleIn

                        integer (si), parameter                  :: nBytes=c_dp
                        integer (singleb),dimension(nBytes)      :: intIn
                        integer (singleb),dimension(nBytes)      :: intReversed

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! Transfer 64 bits of realIn to generic 64 bit INTEGER space:
                        intIn = transfer( dbleIn, intIn )
                        ! Transfer reversed order bytes to 64 bit REAL space
                        intReversed = intIn(nBytes:1:-1)
                        dble_swap_endianess = transfer( intReversed, dble_swap_endianess )

                endfunction dble_swap_endianess

endmodule IO_dcd
