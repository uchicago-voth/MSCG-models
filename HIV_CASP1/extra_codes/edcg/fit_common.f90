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
! This module contains routines to calculate residuals, which are then used in 
! optimization procedures.

module fit_common

        use env_kindtypes

        implicit none

        private

        integer (si),public,parameter :: num_residuals = 3

        public fit_residuals, is_zero

        contains
                subroutine fit_residuals(tdata,mapping,mapping_accumulator_type,residual_weights,&
                                         cgCharges,backmapping,residualList,&
                                         residual_offsets, residual_scalings,&
                                         update, mappingCache, mappingChanges, thin)

                        use obj_trajectory,     only: labeledTraj, traj_training_data
                        use fit_edcg,           only: edcg_residual_spatial
                        use fit_charge,         only: compute_A_B_matrices, fit_charges, charge_residual
                        use routines_math,      only: computeDistanceMatrix
                        use obj_ll,             only: i_sp_dll
                        use core_filter,        only: valueCollapse

                        implicit none

                        !!!! Dummy variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        type (traj_training_data),intent(in   )             :: tdata
                        integer (si),             intent(inout)             :: mapping(:)
                        character(*),             intent(in   )             :: mapping_accumulator_type
                        real    (dp),             intent(in   )             :: residual_weights(:)
                        real    (dp),             intent(inout)             :: cgCharges(:)
                        type(i_sp_dLL),           intent(inout),allocatable :: backmapping(:)
                        logical,                  intent(in   ),optional    :: update
                        integer (si),             intent(inout),optional    :: mappingCache(:), mappingChanges(:)
                        real    (dp),             intent(inout)             :: residualList(:)
                        real    (dp),             intent(in   )             :: residual_offsets(:)
                        real    (dp),             intent(in   )             :: residual_scalings(:)
                        logical,                  intent(in   ),optional    :: thin
                        !!!! End Dummy variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        ! local variables
                        real    (dp), allocatable                :: A(:,:)
                        real    (dp), allocatable                :: B(:,:)
                        real    (dp), allocatable                :: mapped_cgCharges(:)
                        logical,      allocatable                :: residual_mask(:)
                        logical                                  :: update_, thin_

                        type (labeledTraj)                       :: cgTraj

                        if (present(update)) then
                                update_ = update
                        else
                                update_ = .false.
                        endif
                        
                        if (present(thin)) then
                                thin_ = thin
                        else
                                thin_ = .false.
                        endif

                        allocate(residual_mask(size(residual_weights)))
                        residual_mask(:) = is_zero(residual_weights)

                        select case (mapping_accumulator_type)
                        case ("cop")
                                cgTraj = tdata%trj%lmapTrajectory(map = mapping)
                        case ("com")
                                cgTraj = tdata%trj%lmapTrajectory(map = mapping,&
                                                     position_weights = tdata%trj%site_mass)
                        case ("coc")
                                cgTraj = tdata%trj%lmapTrajectory(map = mapping,&
                                                     position_weights = tdata%trj%atomCharges)
                        endselect

                        !the order of residuals goes fit charge, edcg, mapped charge (historical consistency)

                        if (update_) then

                                !Compute CG distance matrix and CG - AA distance matrix.
                                if (.not. (thin_ .and. residual_mask(1) .and. residual_mask(3))) then
                                        call compute_A_B_matrices(tdata%trj%coord,cgTraj%coord,A,B)
                                endif

                                if (.not. (thin_ .and. residual_mask(1))) then
                                        !Fit charges
                                        call fit_charges(A=A,&
                                                         B=B,&
                                                         C=tdata%stats%distanceMat,&
                                                         atomCharges=tdata%trj%atomCharges,&
                                                         cgCharges=cgCharges)

                                        !compute charge residual
                                        residualList(2) = charge_residual(charges_1 = cgCharges,&
                                                                          charges_2 = tdata%trj%atomCharges,&
                                                                          dist_1    = A,&
                                                                          dist_2    = tdata%stats%distanceMat,&
                                                                          dist_12   = B,&
                                                                          residual_offset  = residual_offsets(1),&
                                                                          residual_scaling = residual_scalings(1))
                                else
                                        residualList(2) = -1
                                endif

                                if (.not. (thin_ .and. residual_mask(3))) then
                                        mapped_cgCharges = valueCollapse(indexLabels = mapping,&
                                                                         indexValues = tdata%trj%atomCharges,&
                                                                             nLabels = size(cgCharges,1),&
                                                                        perValueNorm = .false.)

                                        !compute mapped charge residual
                                        residualList(4) = charge_residual(charges_1 = mapped_cgCharges,&
                                                                          charges_2 = tdata%trj%atomCharges,&
                                                                          dist_1    = A,&
                                                                          dist_2    = tdata%stats%distanceMat,&
                                                                          dist_12   = B,&
                                                                          residual_offset  = residual_offsets(3),&
                                                                          residual_scaling = residual_scalings(3))
                                else
                                        residualList(4) = -1
                                endif

                                !currently, we have to call edcg no matter what, as it updates the backmapping.
                                !we simply overwrite the residual value with -1 to be consistent.
                                call edcg_residual_spatial(         residual = residualList(3),       &
                                                             residual_offset = residual_offsets(2),   &
                                                            residual_scaling = residual_scalings(2),  &
                                                                     mapping = mapping,               &
                                                                         nCg = cgTraj%nAtoms,         &
                                                                      update = .true. ,               &
                                                                 backmapping = backmapping,           &
                                                                   isoPcaMat = tdata%stats%isoPcaMat,&
                                                                     changes = mappingChanges,        &
                                                                mappingCache = mappingCache)
                                if ((thin_ .and. residual_mask(2))) then
                                        residualList(3) = -1
                                endif
                        else
                                if (.not. (thin_ .and. residual_mask(1) .and. residual_mask(3))) then
                                        !Accumulate A and B
                                        call compute_A_B_matrices(tdata%trj%coord,cgTraj%coord,A,B)
                                endif

                                !Fit charges
                                if (.not. (thin_ .and. residual_mask(1))) then
                                        call fit_charges(A=A,&
                                                         B=B,&
                                                         C=tdata%stats%distanceMat,&
                                                         atomCharges=tdata%trj%atomCharges,&
                                                         cgCharges=cgCharges)

                                        !compute charge residual
                                        residualList(2) = charge_residual(charges_1 = cgCharges,&
                                                                          charges_2 = tdata%trj%atomCharges,&
                                                                          dist_1    = A,&
                                                                          dist_2    = tdata%stats%distanceMat,&
                                                                          dist_12   = B,&
                                                                          residual_offset  = residual_offsets(1),&
                                                                          residual_scaling = residual_scalings(1))
                                else
                                        residualList(2) = -1
                                endif

                                if (.not. (thin_ .and. residual_mask(3))) then
                                        mapped_cgCharges = valueCollapse(indexLabels = mapping,&
                                                                         indexValues = tdata%trj%atomCharges,&
                                                                             nLabels = size(cgCharges,1),&
                                                                        perValueNorm = .false.)

                                        !compute mapped charge residual
                                        residualList(4) = charge_residual(charges_1 = mapped_cgCharges,&
                                                                          charges_2 = tdata%trj%atomCharges,&
                                                                          dist_1    = A,&
                                                                          dist_2    = tdata%stats%distanceMat,&
                                                                          dist_12   = B,&
                                                                          residual_offset  = residual_offsets(3),&
                                                                          residual_scaling = residual_scalings(3))
                                else
                                        residualList(4) = -1
                                endif


                                !currently, we have to call edcg no matter what, as it updates the backmapping.
                                !we simply overwrite the residual value with -1 to be consistent.
                                call edcg_residual_spatial( residual   = residualList(3),      &
                                                       residual_offset = residual_offsets(2),  &
                                                      residual_scaling = residual_scalings(2), &
                                                               mapping = mapping,              &
                                                                   nCg = cgTraj%nAtoms,        &
                                                                update = .false.,              &
                                                           backmapping = backmapping,          &
                                                             isoPcaMat = tdata%stats%isoPcaMat)
                                if ((thin_ .and. residual_mask(2))) then
                                        residualList(3) = -1
                                endif
                        endif

                        !compute total residual
                        residualList(1) = residual_weights(1)*residualList(2) + &
                                          residual_weights(2)*residualList(3) + &
                                          residual_weights(3)*residualList(4)

                endsubroutine fit_residuals

                elemental function is_zero(num)

                        implicit none

                        real    (dp),intent(in   ) :: num

                        logical                    :: is_zero

                        real    (dp),parameter     :: tol = 10.0**(-14)

                        if (abs(num) < tol) then
                                is_zero = .true.
                        else
                                is_zero = .false.
                        endif

                endfunction is_zero

endmodule fit_common
