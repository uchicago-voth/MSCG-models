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
! Provides routines to read psf files. Note that the full generality of the 
! standard is not supported; only e.g. atoms types are read.
!
! NOTE:
!    The source of this routine is unknown.

module IO_psf

        implicit none

        private

        public                :: read_psf_file

        contains

                !Read atomic charges from some file
                subroutine read_psf_file(AtomPsfFile,nAtoms,nResidues,&
                                atomCharges,atomMasses,atomResidues,atomLabels,&
                                residueStartSites) 

                        use env_kindtypes,    only: si, dp

                        implicit none

                        !!! Dummy variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        
                        character (*), intent(in)               :: atomPsfFile

                        integer   (si),intent(out)              :: nAtoms           !Num of atoms.
                        integer   (si),intent(out)              :: nResidues

                        real      (dp),intent(out), allocatable :: atomCharges(:)
                        character (LEN=7),  intent(out), allocatable :: atomLabels(:)
                        integer   (si),intent(out), allocatable :: atomResidues(:)
                        real      (dp),intent(out), allocatable :: atomMasses(:)
                        integer   (si),intent(out), allocatable :: residueStartSites(:)

                        !!! End dummy variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer   (si)              :: atom
                        character (LEN=6)                :: check   !character to check if NATOM is in the line
                        character (LEN=8)                :: numChar !character to read number of atoms.  
                                                                    !must be converted to an integer
                        integer   (si)              :: resCount 
                        integer   (si)              :: fp_psf  !file pointer to psf file
                        integer   (si),allocatable  :: authentic_ResNumber(:)

                        resCount = 0

                        !open the psf file
                        open(newunit=fp_psf,file=atomPsfFile)

                        print*, "INFO: atom psf file being read: ",trim(atomPsfFile)

                        do

                                read(fp_psf,'(a8,2x,a6)') numChar, check

                                !if we read the number of atoms exit this do loop
                                if (check.eq.'NATOM ') then
                                        !this converts the character natoms_char to the integer natoms
                                        read(numChar,*) nAtoms
                                        write(*,*) "INFO: (psf file) Number of atoms present:", nAtoms

                                        !Now that we know the number of atoms we must allocate the arrays

                                        ! To hold the site charges.
                                        allocate( atomCharges(nAtoms) )

                                        ! To hold the residues of each site 
                                        ! Numerically different from the psf residue 
                                        ! number.
                                        allocate( atomResidues(nAtoms))

                                        ! To hold the mass of each site.
                                        allocate( atomMasses(nAtoms)   )

                                        ! To hold the type (e.g. atom type) of
                                        ! each site.
                                        allocate( atomLabels(nAtoms)   )

                                        ! To temporarily hold residue spec as held 
                                        ! in psf file.
                                        allocate( authentic_ResNumber(nAtoms)    )

                                        !Now we loop through the number of atoms and read the pertinent information
                                        resCount = 0

                                        !read in the data for each atom.
                                        do atom = 1,nAtoms

                                                !initialize the full string
                                                atomLabels(atom) = ''
                                                read(fp_psf,'(14x,i4, 6x,a4, 6x,f10.6, 4x,f10.4)') &
                                                                authentic_ResNumber(atom),    &
                                                                atomLabels(atom)(1:4),   &
                                                                atomCharges(atom), &
                                                                atomMasses(atom)             

                                                !count the number of unique residues.
                                                if (atom.eq.1) then
                                                        resCount = resCount + 1
                                                elseif (authentic_ResNumber(atom) .ne. &
                                                                authentic_ResNumber(atom-1)) then
                                                        resCount = resCount + 1
                                                endif

                                                !store residue index
                                                atomResidues(atom) = resCount

                                        enddo

                                        deallocate(authentic_ResNumber)

                                !if we get to the next section, exit the loop.
                                elseif (check.eq.'NBOND:') then
                                        exit
                                endif

                        enddo

                        close(fp_psf)

                        nResidues = resCount
                        print*, "INFO: (psf file) Number of residues:", nResidues

                        !Create array of residue start atom numbers
                        allocate(residueStartSites(nResidues))

                        resCount = 1
                        do atom = 1,nAtoms

                                if (atom==1) then
                                        residueStartSites(resCount) = atom
                                        resCount                                = resCount + 1
                                elseif (atomResidues(atom) .ne. atomResidues(atom-1)) then
                                        residueStartSites(resCount) = atom
                                        resCount                                = resCount + 1
                                endif

                        enddo

                        write(*,*) "INFO: (psf file) Total Charge of the system = ", sum(atomCharges)

                endsubroutine read_psf_file

endmodule IO_psf
