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
! This module contains routines to calculate the edcg type residuals.

module fit_edcg

        use env_kindtypes, only: si,dp

        implicit none

        private

        ! there is a possiblity of error here. We have a strong need for the ability to update and compute
        ! a residual anew.
        public          :: edcg_residual_spatial

        contains
                real (dp) function edcg_residual(boundaryRes,normalization,isoPcaMat)

                        use env_kindtypes, only: si, dp, qp

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),intent(in)   :: isoPcaMat(:,:)
                        integer (si),intent(in)   :: boundaryRes(:)
                        real    (dp),intent(in)   :: normalization

                        !!!!! End dummy arguments  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)    :: nRes
                        integer (si)    :: nCg
                        integer (si)    :: res1,res2
                        integer (si)    :: cg
                        integer (si)    :: startRes, stopRes
                        real    (qp)    :: edcg_residual_accum

                        integer (si), allocatable&
                                             :: boundaryres_shape(:)
                        integer (si), allocatable&
                                             :: isoPcaMat_shape(:)

                        !zero the residual
                        edcg_residual_accum = 0.0_qp

                        boundaryres_shape   = shape(boundaryRes)
                        isoPcaMat_shape     = shape(isoPcaMat)

                        !the number of cg sites is one more than 
                        !the number of boundary residues.
                        nCg  = boundaryres_shape(1) + 1
                        nRes = isoPcaMat_shape(1)

                        !Iterate across the boundary residues.
                        startRes = 1

                        do cg=1,nCg
                                if (cg<nCg) then
                                        stopRes = boundaryRes(cg)
                                else
                                        stopRes = nRes
                                endif
                                do res1 = startRes,stopRes
                                        do res2 = res1+1,stopRes
                                                        !print*,"comparing",res1,res2
                                                        edcg_residual_accum = edcg_residual_accum&
                                                                        + isoPcaMat(res1,res1)&
                                                                        - 2.0*isoPcaMat(res1,res2)&
                                                                        + isoPcaMat(res2,res2)
                                        enddo
                                enddo

                                startRes = stopRes+1
                        enddo

                        ! Until we settle on a normalization constant, avoid errors if we
                        ! choose ones that are very different.
                        edcg_residual = real(edcg_residual_accum/real(normalization,qp),dp)

                endfunction edcg_residual
                ! public function to return the edcg residual.
                subroutine edcg_residual_spatial(residual,residual_offset,residual_scaling,mapping,nCg,update,mappingCache,&
                                                changes,isoPcaMat,backMapping)

                        use env_kindtypes, only: si, dp, qp
                        use obj_ll,        only: i_sp_dll

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),  intent(inout)                 :: residual
                        real    (dp),  intent(in   )                 :: residual_offset
                        real    (dp),  intent(in   )                 :: residual_scaling
                        real    (dp),  intent(in   )                 :: isoPcaMat(:,:)
                        integer (si),  intent(in   )                 :: mapping(:)
                        integer (si),  intent(in   )                 :: nCg
                        logical,       intent(in   )                 :: update
                        integer (si),  intent(in   ),optional        :: mappingCache(:),changes(:)
                        type(i_sp_dll),intent(inout),allocatable     :: backmapping(:)
                        !!!!! End dummy arguments  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (qp)                   :: edcg_residual_accum
                        integer (si)                   :: iter
                        integer (si)                   :: changedIDindex, changedID, comparison_FGid
                        integer (si)                   :: new_assignment, old_assignment

                        !isoPcaMat_shape     = shape(isoPcaMat)

                        if (update) then
                                edcg_residual_accum = 0
                                !first, subtract terms which are no longer valid
                                do changedIDindex=1,size(changes)
                                        changedID = changes(changedIDindex)
                                        old_assignment = mappingCache(changedID)
                                        delete: do iter=1,backmapping(old_assignment)%getLength()
                                                comparison_FGid = backmapping(old_assignment)%getCurrent()
                                                if (comparison_FGid == changedID) then
                                                        call backmapping(old_Assignment)%delete()
                                                else
                                                        edcg_residual_accum = edcg_residual_accum&
                                                                               - isoPcaMat(changedID,changedID)&
                                                                               + 2.0*isoPcaMat(changedID,comparison_FGid)&
                                                                               - isoPcaMat(comparison_FGid,comparison_FGid)
                                                        call backmapping(old_Assignment)%next()
                                                endif
                                        enddo delete
                                enddo

                                do changedIDindex=1,size(changes)
                                        changedID = changes(changedIDindex)
                                        new_assignment = mapping(changedID)
                                        add: do iter=1,backmapping(new_assignment)%getLength()
                                                comparison_FGid = backmapping(new_assignment)%getCurrent()
                                                edcg_residual_accum = edcg_residual_accum&
                                                                       + isoPcaMat(changedID,changedID)&
                                                                       - 2.0*isoPcaMat(changedID,comparison_FGid)&
                                                                       + isoPcaMat(comparison_FGid,comparison_FGid)
                                                call backmapping(new_Assignment)%next()
                                        enddo add
                                        call backmapping(new_assignment)%prepend(changedID)
                                enddo

                                residual = real(edcg_residual_accum,dp)*residual_scaling + residual
                        else
                                !fresh residual calculation
                                edcg_residual_accum = 0

                                !NOTE THAT HERE WE ALSO BUILD OUR MAPPING!

                                if (allocated(backmapping)) then
                                        if (.not. (size(backmapping,1) .eq. nCg)) then
                                                deallocate(backmapping)
                                                allocate(backmapping(nCg))
                                        endif
                                else
                                        allocate(backmapping(nCg))
                                endif

                                do iter=1,nCg
                                        call backmapping(iter)%initialize()
                                enddo

                                do changedID=1,size(mapping)
                                        new_assignment = mapping(changedID)
                                        if (backmapping(new_assignment)%getLength() == 0) then
                                                call backmapping(new_assignment)%append(changedID)
                                        else
                                                do iter=1,(backmapping(new_assignment)%getLength())
                                                        comparison_FGid = backmapping(new_assignment)%getCurrent()
                                                        edcg_residual_accum = &
                                                             edcg_residual_accum&
                                                             + isoPcaMat(changedID,changedID)&
                                                             - 2.0*isoPcaMat(changedID,comparison_FGid)&
                                                             + isoPcaMat(comparison_FGid,comparison_FGid)
                                                        call backmapping(new_assignment)%next()
                                                enddo
                                                call backmapping(new_assignment)%append(changedID)
                                        endif
                                enddo

                                residual = (real(edcg_residual_accum,dp) + residual_offset) * residual_scaling
                        endif

                endsubroutine edcg_residual_spatial

                !generate a backmapping from an array. backmapping is a dLL.
                !pure subroutine generate_backmapping(backmapping,mapping,nCG)

                !        use env_kindtypes, only: si
                !        use obj_ll,        only: i_sp_dll

                !        implicit none

                !        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !        type(i_sp_dll),intent(inout),allocatable :: backmapping(:)
                !        integer (si),  intent(in   )             :: mapping(:)
                !        integer (si),  intent(in   )             :: nCG
                !        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !        integer (si)    :: FGid, iter

                !        m4_careful_allocate_1(backmapping,nCG)

                !        do iter=1,size(backmapping)
                !                call backmapping(iter)%initialize()
                !        enddo

                !        do FGid=1,size(mapping)
                !                !populate backmapping.
                !                call backmapping(mapping(FGid))%append(FGid)
                !        enddo

                !endsubroutine generate_backmapping
endmodule fit_edcg
