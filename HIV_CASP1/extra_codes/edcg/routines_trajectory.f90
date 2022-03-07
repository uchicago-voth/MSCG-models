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
! Provides more complex supporting routines for the trajectory objects.
!


module routines_trajectory

        use env_kindtypes, only: si,dp

        implicit none

        private

        public          :: compute_covar, align_to_avg, get_square_displacement_distmat

        contains
                subroutine compute_covar(coord,avgCoord,pcaMat,edcgNorm,edDof,isoPcaMat)

                        use env_kindtypes, only: si, dp
                        use io_csv,        only: write_csv
                        use core_matrix,   only: quadForm, spectrum, diagonal

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in)         :: coord(:,:,:)
                        real    (dp),intent(in)         :: avgcoord(:,:)
                        real    (dp),intent(out)        :: pcamat(:,:)
                        real    (dp),intent(out)        :: isopcamat(:,:)
                        real    (dp),intent(out)        :: edcgnorm
                        integer (si),intent(in)         :: eddof
                        !!!!! end dummy arguments  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                    :: natoms
                        integer (si)                    :: nsteps

                        real    (dp),allocatable        :: temppcamat(:,:) ! used in calculation of cov mat.
                        !indexing variables
                        integer (si)                    :: i, j, k, l, trace_index  ! do loops
                        integer (si)                    :: index1, index2           ! do loops
                        integer (si)                    :: step                     ! do loops (timestep index)

                        real    (dp),allocatable        :: eig_mat(:,:)         !! used in projection of cov matrix into
                                                                                 !  essential subspace
                        integer (si),allocatable        :: coord_shape(:)
                        real    (dp),allocatable        :: eigenvectors(:,:)
                        real    (dp),allocatable        :: eigenvalues(:)  !! to contain eigenvalues of covariance matrix

                        ! todo; why do we really need these?
                        !!! holds values for calculating normalization constants. !!!!!!!!!!!!!!!!!
                        real (dp)          :: avgedcgnorm
                        real (dp)          :: avgedcgnorm2
                        real (dp)          :: isotropic_variance

                        coord_shape = shape(coord)

                        natoms = coord_shape(1)
                        nsteps = coord_shape(3)

                        allocate(eigenvalues(3*natoms))
                        allocate(temppcamat(3*natoms,3*natoms))

                        ! Initialized pca matrices and summed averages
                        pcaMat       = 0
                        tempPcaMat   = 0
                        avgEDCGNorm  = 0
                        avgEDCGNorm2 = 0

                        ! Compute covariance matrix. Loop through timesteps and dimensions in trajectory.
                        ! Compute index{1,2} to locate entry in covariance matrix.
                        do step=1,nSteps
                                !loop thgouh atoms
                                do i=1,nAtoms
                                        ! Loop through the dimensions
                                        do j=1,3
                                                ! cov mat index
                                                index1 = (i-1)*3+j
                                                do k=1,nAtoms
                                                        do l=1,3
                                                                index2 = (k-1)*3+l
                                                                !Calculate covariance entry
                                                                tempPcaMat(index1,index2) = &
                                                                ( (coord(i,j,step) - avgCoord(i,j) )*&
                                                                (coord(k,l,step) - avgCoord(k,l)))

                                                                !Symmerterize matrix
                                                                tempPcaMat(index2,index1) = tempPcaMat(index1,index2)
                                                        enddo
                                                enddo
                                        enddo
                                enddo

                                ! Accumulate edcg nomes.
                                pcaMat       = pcaMat       + tempPcaMat
                        enddo

                        ! compute the normalization constant

                        avgEDCGNorm  = avgEDCGNorm  / real(nSteps,dp)
                        avgEDCGNorm2 = avgEDCGNorm2 / real(nSteps,dp)
                        edcgNorm     = avgEDCGNorm2 - avgEDCGNorm*avgEDCGNorm

                        !time average and finalize covariance
                        pcaMat = pcaMat/real(nSteps,dp)

                        !project onto the essential subspace if the given number of dimensions implies it.
                        if (3*edDof < size(pcaMat,1)) then

!                                call spectrum(         M = pcaMat,&
!                                                       d = eigenvalues,&
!                                                       Q = eigenvectors,&
!                                     num_top_eigenvalues = edDof)
                                call spectrum(         M = pcaMat,&
                                                       d = eigenvalues,&
                                                       Q = eigenvectors)
        
                                eigenvalues(1:(3*natoms-edDof)) = 0
                                eig_mat = diagonal(values=eigenvalues)
        
                                pcaMat = quadForm(eig_mat,eigenvectors,rightTranspose=.true.)
                        endif

                        isoPcaMat          = 0.0_dp
                        isotropic_variance = 0.0_dp
                        do i=1,nAtoms
                                do j=i,nAtoms
                                        do trace_index=1,3
                                                index1 = (i-1)*3+trace_index
                                                index2 = (j-1)*3+trace_index
                                                isotropic_variance = isotropic_variance &
                                                                + pcaMat(index1,index2)
                                        enddo
                                        isoPcaMat(i,j)     = isotropic_variance
                                        isoPcaMat(j,i)     = isotropic_variance
                                        isotropic_variance = 0.0_dp
                                enddo
                        enddo

                endsubroutine compute_covar

                function get_square_displacement_distmat(isoCovMat) result(distMat)

                        use env_kindtypes, only: si, dp

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )        :: isoCovMat(:,:)
                        !!!!! end dummy arguments  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp),allocatable          :: distMat(:,:)
                        integer (si)                      :: iter1, iter2, nSites

                        nSites = size(isoCovMat,1)

                        allocate(distMat(nSites,nSites))

                        do iter1=1,nSites-1
                        do iter2=iter1+1,nSites
                                distMat(iter1,iter2) = isoCovMat(iter1,iter1) + isoCovMat(iter2, iter2) - 2*isoCovMat(iter1,iter2)
                                distMat(iter2,iter1) = distMat(iter1,iter2)
                        enddo
                        enddo

                        do iter1=1,nSites
                                distMat(iter1,iter1) = 0
                        enddo

                endfunction get_square_displacement_distmat

                pure subroutine align_to_avg(coord,refAvg,center,status)

                        use env_kindtypes, only: si, dp

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp), intent(inout) :: coord(:,:,:) !! Mutable: this is what we are aligning.
                        real    (dp), intent(inout) :: refAvg(:,:)  !! Mutable: Fit is
                                                                    !   done iteratively, modifying
                                                                    !   the reference.
                        logical,      intent(in), optional&
                                                    :: center       !! Whether to center the protein before aligning.
                        integer (si), intent(out),optional&
                                                    :: status       !! Exit code
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional argument counterparts
                        integer (si)                :: out_status          !! Exit code

                        integer (si), parameter     :: maxIter = 100       !! Maximum number of iterations in
                                                                           !   structure alignment.
                        real    (dp), parameter     :: thresh  = 1E-5_dp   !! Threshold to stop alignment early.

                        integer (si)                :: nAtoms, nSteps, &   !! Derived parameters to control loops.
                                                            dimen
                        real    (dp)                :: residual            !! Fit residual for early exit.
                        integer (si)                :: iter                !! Current iteration.
                        integer (si)                :: step                !! timestep index.

                        integer (si), allocatable   :: shape_coord(:), &   !! Holds array shape of input.
                                                            shape_refAvg(:)
                        real    (dp), allocatable   :: newAvg(:,:)         !! Holds average structure during iteration.

                        logical                     :: center_             !! Whether to center the protein before aligning.

                        !parse optional arguments
                        if (present(center)) then
                                center_ = center
                        else
                                center_ = .true.
                        endif

                        !!!!! Start program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !Initialize status as failure.
                        out_status = 1

                        !get shape of arguments.
                        shape_coord  = shape(coord)
                        shape_refAvg = shape(refAvg)

                        if ( any( shape_coord(1:2).ne.shape_refAvg ) ) then
                                !cause failure.
                                coord      = -1
                                refAvg     = -1
                                out_status =  2
                                if( present(status) ) then
                                        status = out_status
                                endif
                                return
                        endif

                        nAtoms = shape_coord(1)
                        dimen  = shape_coord(2)
                        nSteps = shape_coord(3)

                        allocate(newAvg(nAtoms,dimen))

                        if (center_) then
                              call center_CoordSequence(coord)
                              call center_frame(refAvg)
                        endif

                        do iter = 1, maxIter

                                !zero the average.
                                newAvg = 0_dp

                                !orient each step and accum averge.
                                do step=1,nSteps

                                        !orient to reference using Kabsch algorithm
                                        call orient_to_reference(coord(:,:,step),refAvg)

                                        !average aligned coordinates, normalize later.
                                        newAvg = newAvg + coord(:,:,step)

                                enddo

                                !Normalize the positions of each atom
                                newAvg = newAvg/real(nSteps,dp)

                                !residual is just the RMSD of the average from this step and the last
                                residual =  compute_residual(newAvg,refAvg)

                                !update the reference structure.
                                refAvg = newAvg

                                !check convergence
                                if (residual < thresh) then
                                        out_status = 1
                                        exit
                                endif
                        enddo

                        deallocate(newAvg)

                        !success
                        out_status = 0

                        !if they asked for status, return it.
                        if( present(status) ) then
                                status = out_status
                        endif

                        return

                endsubroutine align_to_avg

                ! This subroutine calculates the optimal rotation matrix to rotate the new set
                ! of coordinates to the reference set using the Kabsch algorithm. This subroutine
                ! assumes that refCoord and alpha_pos have their COM at the origin.
                pure subroutine orient_to_reference(coord,refCoord)

                        use env_kindtypes, only: si, dp
                        use core_matrix,   only: svd

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: coord(:,:)    !! must be able to modify structure.
                        real    (dp),intent(in   )          :: refCoord(:,:) !! don't modify reference
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si)               :: i, j                !! loop indices

                        real    (dp),allocatable   :: singular_values(:)
                        real    (dp),allocatable   :: covar(:,:)
                        real    (dp),allocatable   :: right_eigv_T(:,:)
                        real    (dp),allocatable   :: left_eigv(:,:)

                        real    (dp)               :: rot(3,3)
                        real    (dp)               :: right_handed(3,3)   !! Matrix to enforce
                                                                          !  rotation, not reflection, from kabsch
                        real    (dp)               :: det

                        !The first step in the Kabsch algorithm is to compute a 3x3 covariance matrix
                        !between new set of coordinates and reference set of coordinates compute the
                        !covariance matrix between this time point and the reference

                        covar = matmul(transpose(coord),refCoord)

                        !Compute the determinant of covariance matrix
                        det   =  determinant_three_three(covar)

                        call svd(covar,singular_values,left_eigv,right_eigv_T)

                        !Create the rotation/reflection decision matrix (diagonal)
                        do i=1, 3
                             do j=1, 3
                                  if (i.eq.3 .and. j.eq.3 ) then
                                       right_handed(i,j) = abs(det)/det
                                  elseif (i.eq.j) then
                                       right_handed(i,j) = 1
                                  else
                                       right_handed(i,j) = 0
                                  endif
                             enddo
                        enddo

                        !Compute the rotation matrix
                        left_eigv = matmul(left_eigv,right_handed)
                        rot       = matmul(left_eigv,right_eigv_T)

                        !Rotate the coordinates by applying matrix.
                        coord     = matmul(coord,rot)

                endsubroutine orient_to_reference

                !calculate the determinant of a three by three matrix
                pure real (dp) function determinant_three_three(mat)

                        use env_kindtypes, only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real (dp),intent(in)     :: mat(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        determinant_three_three =   mat(1,1)*mat(2,2)*mat(3,3)&
                                                  - mat(1,1)*mat(2,3)*mat(3,2)&
                                                  - mat(1,2)*mat(2,1)*mat(3,3)&
                                                  + mat(1,2)*mat(2,3)*mat(3,1)&
                                                  + mat(1,3)*mat(2,1)*mat(3,2)&
                                                  - mat(1,3)*mat(2,2)*mat(3,1)

                endfunction determinant_three_three

                pure real (dp) function compute_residual(coord,refCoord)

                        use env_kindtypes, only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in)           :: coord(:,:)
                        real    (dp),intent(in)           :: refCoord(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)             :: nAtoms, i, j
                        real    (dp)             :: temp

                        nAtoms = size(coord,1)

                        compute_residual = 0
                        do i=1,nAtoms

                           do j=1,3
                              temp = coord(i,j)-refCoord(i,j)
                              compute_residual = compute_residual + temp*temp
                           enddo

                        enddo

                        compute_residual = real( sqrt( compute_residual/dble(nAtoms) ) )

                endfunction compute_residual

                !Center a coordinate frame to be zero.
                pure subroutine center_frame(coord)

                        use env_kindtypes, only: si, dp

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp), intent(inout)     :: coord(:,:) !! Mutable: this is what we centering
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        real    (dp)                    :: siteMeanPos(size(coord,2))
                        integer (si)                    :: i

                        !Center the trajectory.
                        siteMeanPos     = 0
                        do i=1,size(coord,1)
                               siteMeanPos = siteMeanPos + coord(i,:)
                        enddo
                        siteMeanPos = siteMeanPos/size(coord,1)

                        !Apply shift
                        coord = coord - spread(siteMeanPos,1,size(coord,1))

                endsubroutine center_frame

                !Center a coordinate frame to be zero.
                pure subroutine center_coordSequence(coordSequence)

                        use env_kindtypes, only: si, dp

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp), intent(inout)     :: coordSequence(:,:,:) !! Mutable: this is
                                                                                !   what we centering
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                    :: i

                        do i=1,size(coordSequence,3)
                                call center_frame(coordSequence(:,:,i))
                        enddo

                endsubroutine center_coordSequence

endmodule routines_trajectory
