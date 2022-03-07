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
! Provides more complex mathematical routines which aren't simple enough for a core module,
! but don't easily reside in other places. In general, this module should be kept as sparse as possible.
!
! NOTE:
!    If compilation order permits, consider decomposing this module and reloating routines.

module routines_math

        implicit none

        private

        public electroPot, computeDistanceMatrix, computeAsymmetricDistanceMatrix, diameter

        contains
                !Note that this takes in a nSites x ncoords arrray of positions.
                pure function computeDistanceMatrix(positions,secondaryPositions,ele_inverse) &
                                                result(distanceMatrix)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in)             :: positions(:,:)
                        real    (dp),intent(in),optional    :: secondaryPositions(:,:)
                        logical,     intent(in),optional    :: ele_inverse
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        ! Return value
                        real    (dp),allocatable   :: distanceMatrix(:,:)

                        !optional value counterparts
                        logical                             :: ele_inverse_

                        !parse optional arguments
                        if (present(ele_inverse)) then
                                ele_inverse_ = ele_inverse
                        else
                                ele_inverse_ = .false.
                        endif

                        !calculate the appropriate matrix.

                        if (.not. present(secondaryPositions) ) then
                                allocate(distanceMatrix(   &
                                         size(positions,1),&
                                         size(positions,1)))

                                call computeSymmetricDistanceMatrix(positions,distanceMatrix)

                                !invert elementwise if we need to.
                                if (ele_inverse_) distanceMatrix = 1.0_dp/distanceMatrix

                        else
                                allocate(distanceMatrix(   &
                                         size(positions,1),&
                                         size(secondaryPositions,1)))

                                call computeAsymmetricDistanceMatrix(positions,secondaryPositions,distanceMatrix)

                                !invert elementwise if we need to.
                                if (ele_inverse_) distanceMatrix = 1.0_dp/distanceMatrix

                        endif

                        return

                endfunction computeDistanceMatrix

                pure subroutine computeAsymmetricDistanceMatrix(positions,secondaryPositions,storage)
                        use env_kindtypes,      only: dp, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in)             :: positions(:,:)
                        real    (dp),intent(in)             :: secondaryPositions(:,:)
                        real    (dp),intent(out)          :: storage(:,:)           ! Array to do work in and
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        ! Internal variables
                        integer (si)               :: index1, index2
                        integer (si)               :: coordIndex, nDim
                        integer (si)               :: primaryNSites, secondaryNSites
                        real    (dp)               :: distance2, temp

                        !first dim is nSites, section is the dimension of the space.
                        primaryNSites   = size(positions,1)
                        secondaryNSites = size(secondaryPositions,1)

                        !if the dims don't match we'll crash.
                        nDim            = size(positions,2)

                        storage = 0

                        !loop though the sites.
                        do index1 = 1, primaryNSites
                                do index2 = 1, secondaryNSites
                                        distance2 = 0
                                        !iterate through the coordates of each point.
                                        do coordIndex=1,nDim
                                                temp = positions(index1,coordIndex) - secondaryPositions(index2,coordIndex)
                                                distance2 = distance2 + temp*temp
                                        enddo

                                        storage(index1,index2) = sqrt(distance2)
                                enddo
                        enddo

                endsubroutine computeAsymmetricDistanceMatrix

                pure subroutine computeSymmetricDistanceMatrix(positions,storage)

                        use env_kindtypes,      only: dp, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in)             :: positions(:,:)
                        real    (dp),intent(out)            :: storage(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        ! Internal variables
                        integer (si)               :: index1, index2
                        integer (si)               :: coordIndex, nDim
                        integer (si)               :: primaryNSites
                        real    (dp)               :: distance2, temp

                        !first dim is nSites, section is the dimension of the space.

                        primaryNSites = size(positions,1)
                        nDim          = size(positions,2)

                        storage = 0

                        !loop though the sites.
                        do index1 = 1, (primaryNSites - 1)
                                do index2 = (index1 + 1), primaryNSites
                                        distance2 = 0
                                        !iterate through the coordates of each point.
                                        do coordIndex=1,nDim
                                                temp = positions(index1,coordIndex) - positions(index2,coordIndex)
                                                distance2 = distance2 + temp*temp
                                        enddo

                                        storage(index1,index2) = sqrt(distance2)

                                        ! symmetrize the matrix
                                        storage(index2,index1) = storage(index1,index2)
                                enddo
                        enddo
                endsubroutine computeSymmetricDistanceMatrix

                pure function electroPot(distanceMatrix,charges) result(ESP)

                        use env_kindtypes,      only: si,dp

                        implicit none

                        !!! Dummy Variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),intent(in)  :: distanceMatrix(:,:)
                        real    (dp),intent(in)  :: charges(:)
                        real    (dp)             :: ESP

                        !!! End Dummy Variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable :: inverseDistanceMatrix(:,:)
                        real    (dp),allocatable :: tempMat(:,:)
                        real    (dp),allocatable :: atomChargesM(:,:)
                        integer (si)             :: iter1,iter2
                        integer (si),allocatable :: chargeShape(:)

                        chargeShape     = shape(charges)
                        allocate( atomChargesM(chargeShape(1),1) )

                        atomChargesM(:,1) = charges

                        allocate(inverseDistanceMatrix(size(distanceMatrix,1),size(distanceMatrix,1)))

                        do iter1=1,size(inverseDistanceMatrix,1)-1
                        do iter2=iter1+1,size(inverseDistanceMatrix,1)
                                inverseDistanceMatrix(iter1,iter2) = -1/distanceMatrix(iter1,iter2)
                                inverseDistanceMatrix(iter2,iter1) = inverseDistanceMatrix(iter1,iter2) 
                        enddo
                        enddo

                        do iter1 = 1, size(inverseDistanceMatrix,1)
                                inverseDistanceMatrix(iter1,iter1) = 0
                        enddo

                        !compute atom ESP for this step
                        tempMat = matmul(transpose(atomChargesM),&
                                        matmul(inverseDistanceMatrix,atomChargesM))

                        ESP = tempMat(1,1)

                endfunction electroPot

                pure function diameter(positions)

                        use env_kindtypes,      only: dp, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )      :: positions(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        ! Internal variables
                        integer (si)               :: index1, index2
                        integer (si)               :: primaryNSites

                        real    (dp)               :: diameter
                        real    (dp)               :: distance_to_compare

                        !first dim is nSites, section is the dimension of the space.

                        primaryNSites = size(positions,1)

                        diameter = 0

                        !loop though the sites.
                        do index1 = 1, (primaryNSites - 1)
                                do index2 = (index1 + 1), primaryNSites
                                        distance_to_compare = &
                                                distance(positions(index1,:),positions(index2,:))
                                        if (distance_to_compare > diameter) then
                                                diameter = distance_to_compare
                                        endif
                                enddo
                        enddo

                endfunction diameter

                pure function distance(x1,x2)

                        use env_kindtypes, only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )      :: x1(:),x2(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp)    :: distance

                        distance = sqrt(sum((x1 - x2)**2))

                endfunction distance

endmodule routines_math
