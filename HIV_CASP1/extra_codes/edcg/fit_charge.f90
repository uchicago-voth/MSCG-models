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
! This module contains the tools to fit charges based on an analytical approach.
!

module fit_charge

        implicit none

        private

        public compute_A_B_matrices, update_A_B_matrices, fit_charges, charge_residual

        contains

                pure subroutine compute_A_B_matrices(atomPos,cgPos,A,B)

                        use env_kindtypes, only: si, dp
                        use routines_math, only: computeDistanceMatrix

                        implicit none

                        !!!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in)    :: atomPos(:,:,:), cgPos(:,:,:)
                        real    (dp),intent(out),allocatable &
                                                        :: A(:,:), B(:,:)
                        !!!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable   :: ASummand(:,:), BSummand(:,:)
                        integer (si)               :: step
                        integer (si)               :: nAtoms, nSteps, nCg
                        integer (si),allocatable   :: shapeAtomPos(:)
                        integer (si),allocatable   :: shapeCgPos(:)

                        !Derive shape of the input matrices.
                        shapeAtomPos    = shape(atomPos)
                        nAtoms          = shapeAtomPos(1)
                        nSteps          = shapeAtomPos(3)

                        shapeCgPos      = shape(cgPos)
                        nCg             = shapeCgPos(1)

                        !allocate the matrices in which to accumulate distance sums.
                        if ( allocated(A) .eqv. .false. ) then
                                allocate(A(nCg,nCg))
                        endif
                        A = 0_dp

                        if ( allocated(B) .eqv. .false. ) then
                                allocate(B(nCg,nAtoms))
                        endif
                        B = 0_dp

                        !accumulate distances
                        do step=1,nSteps
                                ASummand = computeDistanceMatrix(cgPos(:,:,step))
                                A = A - ASummand

                                BSummand = computeDistanceMatrix(cgPos(:,:,step),atomPos(:,:,step))
                                B = B - BSummand
                        enddo

                endsubroutine compute_A_B_matrices

                !requires lapack
                !updates A and B matrices due to change of one boundary residue
                pure subroutine update_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B,nSteps,updatedBR)

                        use env_kindtypes, only: si, dp
                        !use env_openmp

                        implicit none

                        integer (si),intent(in)      :: nAtoms
                        integer (si),intent(in)      :: nCg
                        integer (si),intent(in)      :: nSteps
                        integer (si),intent(in)      :: updatedBR
                        real    (dp),intent(in)      :: atomPos(nAtoms,3,nSteps)
                        real    (dp),intent(in)      :: cgPos(nCg,3,nSteps)
                        real    (dp),intent(inout)   :: A(nCg,nCg)
                        real    (dp),intent(inout)   :: B(nCg,nAtoms)

                        real    (dp)   :: dist, temp
                        !loop indeces
                        integer (si)         :: cgSite1, cgSite2
                        integer (si)         :: j
                        !atom2 not currently used
                        !integer atom1, atom2
                        integer (si)         :: atom1
                        integer (si)         :: step

                        !zero the appropriate rows and columns
                        A(updatedBR,:)   = 0
                        A(updatedBR+1,:) = 0
                        A(:,updatedBR)   = 0
                        A(:,updatedBR+1) = 0
                        B(updatedBR,:)   = 0
                        B(updatedBR+1,:) = 0
                        !Compute the distance between the CG sites
                        do step=1,nSteps
                                do cgSite1 = updatedBR, updatedBR+1
                                        do cgSite2 = 1,nCg
                                                if (cgSite1 .ne. cgSite2) then
                                                        dist = 0
                                                        do j=1,3
                                                                temp = cgPos(cgSite1,j,step)-cgPos(cgSite2,j,step)
                                                                dist = dist + temp*temp
                                                        enddo
                                                        dist = sqrt(dist)
                                                        A(cgSite1,cgSite2) = A(cgSite1,cgSite2)-dist
                                                        !symmetrize the matrix
                                                        A(cgSite2,cgSite1) = A(cgSite1,cgSite2)
                                                endif
                                        enddo
                                enddo

                                do cgSite1 = updatedBR, updatedBR+1
                                        do atom1 = 1,nAtoms
                                                dist = 0
                                                do j=1,3
                                                        temp = cgPos(cgSite1,j,step)-atomPos(atom1,j,step)
                                                        dist = dist + temp*temp
                                                enddo
                                                B(cgSite1,atom1) = B(cgSite1,atom1)-sqrt(dist)
                                        enddo
                                enddo
                        enddo

                endsubroutine update_A_B_matrices


                !subroutine to compute CG charges from input matrices A and B.
                !The residual sum of squares is calculated with aid of all-atom matrix C
                subroutine fit_charges(A, B, C, atomCharges, cgCharges)

                !not pure becuase of lapack
                        use env_kindtypes, only: dp, si
                        use core_matrix,   only: lsolve

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),intent(in   )        :: A(:,:)
                        real    (dp),intent(in   )        :: B(:,:)
                        real    (dp),intent(in   )        :: C(:,:)
                        real    (dp),intent(in   )        :: atomCharges(:)
                        real    (dp),intent(  out)        :: cgCharges(:)
                        !!!!! End dummy arguments  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !we could overflow the stack if we aren't careful.
                        real    (dp)                 :: ATemp(size(A,1),&
                                                                   size(A,2))
                        real    (dp)                 :: BTemp(size(B,1),&
                                                                   size(B,2))
                        real    (dp)                 :: D(    size(A,1)-1,&
                                                                   size(A,1))

                        real    (dp)                 :: newB(size(B,1))

                        integer (si)                 :: nAtoms
                        integer (si)                 :: nCg
                        integer (si)                 :: i

                        integer (si),allocatable     :: shapeA(:),shapeB(:),shapeC(:)

                        shapeA  =  shape(A)
                        shapeB  =  shape(B)
                        shapeC  =  shape(C)

                        nCg     = shapeA(1)
                        nAtoms  = shapeC(1)

                        !First we need to modify A and B to have the correct matrix properties
                        !create D matrices
                        D=0
                        do i=1,nCg-1
                                D(i,i)=1
                                D(i,i+1)=-1
                        enddo
                        !multiply A by D0 and B by D1 giving new matrices A and B the correct behavior
                        ATemp(1:(nCg-1),:) = matmul(D,A)
                        BTemp(1:(nCg-1),:) = matmul(D,B)
                        !generate new matrices with last line having 1.0s forcing charge conservation
                        ATemp(nCg,:) = 1.0_dp
                        BTemp(nCg,:) = 1.0_dp

                        !Now use lapack routine to solve least squares problem A*cgCharges=B*atomCharges
                        newB = matmul(BTemp,atomCharges)

                        cgCharges = lsolve(Atemp,newB)

                endsubroutine fit_charges

                pure function charge_residual(charges_1,charges_2,dist_1, dist_2, dist_12, residual_offset, residual_scaling) result(residual)

                        use env_kindtypes, only: dp
                        use core_matrix,   only: quadForm

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),intent(in   )          :: charges_1(:),charges_2(:)
                        real    (dp),intent(in   )          :: dist_1(:,:), dist_2(:,:), dist_12(:,:)
                        real    (dp),intent(in   ),optional :: residual_offset, residual_scaling

                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp)                      :: residual_offset_, residual_scaling_
                        real    (dp)                      :: residual,cross_term

                        if (present(residual_offset)) then
                                residual_offset_ = residual_offset
                        else
                                residual_offset_ = 0
                        endif

                        if (present(residual_scaling)) then
                                residual_scaling_ = residual_scaling
                        else
                                residual_scaling_ = 1
                        endif

                        cross_term = -2 * dot_product( charges_1, matmul(dist_12,charges_2) )

                        residual =     quadForm(dist_1,charges_1) + quadForm(dist_2,charges_2) + cross_term

                        residual = (residual + residual_offset_) * residual_scaling_
                              
                endfunction 
endmodule fit_charge
