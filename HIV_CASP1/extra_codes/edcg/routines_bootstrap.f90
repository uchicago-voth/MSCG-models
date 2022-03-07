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
! This module provides routines for boostrapping.
!

module routines_bootstrap

        implicit none

        private

        public calc_bootstrap_CIs, solve_resampled_models

        contains

                subroutine calc_bootstrap_CIs(models,nModels,tdata,nBSReps,bca,nCAranks)

                        use env_kindtypes,      only: si
                        use core_random,        only: gen_rand_int
                        use obj_trajectory,     only: traj_training_data
                        use abs_obj_model,      only: model

                        implicit none

                        class(model),              intent(inout) :: models(:)
                        class(traj_training_data), intent(in   ) :: tdata
                        integer (si),              intent(in   ) :: nBSReps, nCAranks, nModels
                        logical,                   intent(in   ) :: bca

                        integer (si)                             :: nModels_safe

                        nModels_safe = min(size(models),nModels)

                        !if no reps are requested, return
                        if (nBSReps < 1) return

                        if (nCAranks == 1) then
                                if (bca) then
                                        call calc_jackknife_dists_omp(models(1:nModels_safe),tdata)
                                endif
                                call calc_bootstrap_dists_omp(models(1:nModels_safe),tdata,nBSReps)
                        else
                                !coarray inclusion goes here.
                                print*, 'Run with multiple mpi ranks, but MPI &
                                        &parallelism is not implemented. Stopping.'
                                stop
                        endif

                endsubroutine calc_bootstrap_CIs

                subroutine solve_resampled_models(result_models,optimism_dist,&
                                                this_model,tdata,nDraws,nCAranks)

                        use env_kindtypes,      only: si, dp
                        use core_random,        only: gen_rand_int
                        use obj_trajectory,     only: traj_training_data
                        use abs_obj_model,      only: model

                        implicit none

                        class(model),allocatable,  intent(  out) :: result_models(:)
                        real    (dp),allocatable,  intent(  out) :: optimism_dist(:)
                        class(model),              intent(in   ) :: this_model
                        class(traj_training_data), intent(in   ) :: tdata
                        integer (si),              intent(in   ) :: nDraws, nCAranks

                        if (nDraws == 0) then 
                                allocate(result_models(0),mold=this_model)
                                allocate(optimism_dist(0))
                                return
                        endif

                        if (nCAranks == 1) then
                                call solve_bootstrapped_models_smp(result_models,optimism_dist,&
                                                                 this_model,tdata,nDraws)
                        else
                                !coarray inclusion goes here.
                                print*, 'Run with multiple mpi ranks, but MPI &
                                        &parallelism is not implemented. Stopping.'
                                stop
                        endif

                endsubroutine solve_resampled_models

                subroutine calc_bootstrap_dists_omp(models,src_tdata,nReplicas)

                        use env_kindtypes,      only: si, si_x
                        use core_random,        only: gen_rand_int
                        use obj_trajectory,     only: traj_training_data
                        use abs_obj_model,      only: model

                        implicit none

                        interface
                                pure function omp_get_num_procs()
                                        import si_x
                                        integer (si_x) :: omp_get_num_procs
                                endfunction
                        end interface

                        class(model),              intent(inout) :: models(:)
                        class(traj_training_data), intent(in   ) :: src_tdata
                        integer (si),              intent(in   ) :: nreplicas

                        class(traj_training_data),allocatable    :: job_tdatas(:)

                        integer (si)                    :: default_chunk_size, chunk_size
                        integer (si)                    :: jobs_remaining

                        integer (si)                    :: prep_iter

                        !to avoid using extensive memory
                        default_chunk_size = omp_get_num_procs()*16

                        prep_loop: do prep_iter=1,size(models)
                                call models(prep_iter)%initializeCI(nReplicas)
                        enddo prep_loop


                        allocate(job_tdatas(default_chunk_size),mold=src_tdata)

                        jobs_remaining = nReplicas

                        job_loop: do while (.true.)

                                chunk_size = min(default_chunk_size,jobs_remaining)

                                call sgen_bootstrapped_trajs_omp(job_tdatas(1:chunk_size)%trj,src_tdata%trj)

                                !$omp parallel private (prep_iter)
                                !$omp do
                                do prep_iter=1,chunk_size
                                        call job_tdatas(prep_iter)%preprocess(src_tdata%preproc_design)
                                enddo
                                !$omp end do
                                !$omp end parallel

                                call calc_residuals_omp(models,job_tdatas(1:chunk_size))

                                jobs_remaining = jobs_remaining - chunk_size
                                if (jobs_remaining <= 0) then
                                        exit job_loop
                                endif

                        enddo job_loop

                endsubroutine calc_bootstrap_dists_omp

                subroutine sgen_bootstrapped_trajs_omp(result_trajs,src_traj)

                        use env_kindtypes,      only: si
                        use core_random,        only: gen_rand_int
                        use obj_trajectory,     only: traj

                        implicit none

                        class(traj),            intent(inout) :: result_trajs(:)
                        class(traj),            intent(in   ) :: src_traj

                        integer (si), allocatable             :: rindices(:,:)
                        integer (si)                          :: nReplicas

                        nReplicas = size(result_trajs)

                        allocate(rindices(src_traj%nSteps,nReplicas))

                        !This isn't pure, so likely won't be concurrent.
                        gen_random: block
                                integer (si) :: job_iter,step_iter

                                do job_iter =1,size(rindices,2)
                                do step_iter=1,size(rindices,1)
                                        rindices(step_iter,job_iter) = gen_rand_int(src_traj%nSteps,1)
                                enddo
                                enddo
                        end block gen_random

                        design_trajs: block
                                integer (si) :: job_iter
                                !$omp parallel private (job_iter)
                                !$omp do
                                do job_iter=1,size(result_trajs)
                                        call src_traj%bootstrapCopy(copy    = result_trajs(job_iter),&
                                                                    indices = rindices(:,job_iter))
                                enddo
                                !$omp end do
                                !$omp end parallel
                        end block design_trajs

                endsubroutine sgen_bootstrapped_trajs_omp

                subroutine sgen_jackknifed_trajs_omp(result_trajs,src_traj,indices,istatus)

                        use env_kindtypes,      only: si
                        use obj_trajectory,     only: traj
                        use core_stat,          only: seq_omit

                        implicit none

                        class(traj),               intent(inout) :: result_trajs(:)
                        class(traj),               intent(in   ) :: src_traj
                        integer (si),              intent(in   ) :: indices(:)
                        integer (si),optional,     intent(  out) :: istatus

                        integer (si),allocatable                 :: jk_indices(:)
                        integer (si)                             :: max_index

                        max_index = src_traj%nSteps

                        if (size(result_trajs) /= src_traj%nSteps) then
                                if (present(istatus)) istatus = 1
                                return
                        endif

                        design_trajs: block
                                integer (si) :: job_iter
                                !$omp parallel private (jk_indices)
                                !$omp do
                                do job_iter=1,size(result_trajs)
                                        jk_indices = seq_omit(upper = max_index, &
                                                              lower = 1,         &
                                                              omit  = [ indices(job_iter) ])
                                        call src_traj%bootstrapCopy(copy    = result_trajs(job_iter),&
                                                                    indices = jk_indices)
                                enddo
                                !$omp end do
                                !$omp end parallel
                        end block design_trajs

                        if (present(istatus)) istatus = 0

                endsubroutine sgen_jackknifed_trajs_omp

                subroutine calc_jackknife_dists_omp(models,src_tdata)

                        use env_kindtypes,      only: si, si_x
                        use core_random,        only: gen_rand_int
                        use obj_trajectory,     only: traj_training_data, labeledTraj
                        use abs_obj_model,      only: model
                        use core_stat,          only: seq

                        implicit none

                        interface
                                pure function omp_get_num_procs()
                                        import si_x
                                        integer (si_x) :: omp_get_num_procs
                                endfunction
                        end interface

                        class(model),              intent(inout) :: models(:)
                        class(traj_training_data), intent(in   ) :: src_tdata

                        class(traj_training_data),allocatable    :: job_tdatas(:)

                        integer (si)                    :: default_chunk_size, chunk_size
                        integer (si)                    :: projection_DOF
                        integer (si)                    :: jobs_remaining

                        integer (si)                    :: prep_iter, current_pos

                        integer (si),allocatable        :: chunk_indices(:)


                        prep_loop: do prep_iter=1,size(models)
                                call models(prep_iter)%initializeCI(src_tdata%trj%nSteps,jk=.true.)
                        enddo prep_loop

                        !to avoid using extensive memory
                        jobs_remaining     = src_tdata%trj%nSteps
                        default_chunk_size = omp_get_num_procs()*4

                        projection_DOF = src_tdata%stats%edDOF
                        current_pos = 1

                        allocate(job_tdatas(default_chunk_size),mold=src_tdata)

                        job_loop: do while (.true.)

                                chunk_size = min(default_chunk_size,jobs_remaining)

                                chunk_indices = seq(current_pos,current_pos+chunk_size - 1)
                                current_pos = current_pos + chunk_size

                                call sgen_jackknifed_trajs_omp(job_tdatas(1:chunk_size)%trj,src_tdata%trj,&
                                                               chunk_indices)

                                !$omp parallel private (prep_iter)
                                !$omp do
                                do prep_iter=1,chunk_size
                                        call job_tdatas(prep_iter)%preprocess(src_tdata%preproc_design)
                                enddo
                                !$omp end do
                                !$omp end parallel

                                call calc_residuals_omp(models,job_tdatas(1:chunk_size), jk=.true.)

                                jobs_remaining = jobs_remaining - chunk_size
                                if (jobs_remaining <= 0) then
                                        exit job_loop
                                endif

                        enddo job_loop

                endsubroutine calc_jackknife_dists_omp

                subroutine calc_residuals_omp(models,tdatas,jk)

                        use env_kindtypes,      only: si, dp
                        use obj_trajectory,     only: traj_training_data
                        use abs_obj_model,      only: model
                        use fit_common,         only: num_residuals

                        implicit none

                        class(model),             intent(inout) :: models(:)
                        class(traj_training_data),intent(in   ) :: tdatas(:)
                        logical,       optional,  intent(in   ) :: jk

                        real    (dp)                  :: residuals(num_residuals + 1)

                        logical                       :: jk_

                        integer (si)    :: tdata_iter, model_iter

                        if (present(jk)) then
                                jk_ = jk
                        else
                                jk_ = .false.
                        endif

                        do tdata_iter =1,size(tdatas)
                                !$omp parallel private(residuals)
                                !$omp do
                                do model_iter=1,size(models)
                                        residuals = models(model_iter)%recalcResidual(tdatas(tdata_iter))
                                        call models(model_iter)%addCIValues([ residuals(1) ], jk_)
                                enddo
                                !$omp end do
                                !$omp end parallel
                        enddo

                        !$omp parallel
                        !$omp do
                        do model_iter=1,size(models)
                                call models(model_iter)%calcCI(bca=.true.)
                        enddo
                        !$omp end do
                        !$omp end parallel
        
                endsubroutine calc_residuals_omp

                subroutine solve_bootstrapped_models_smp(result_models,optimism_dist,this_model,tdata,nBSReps)

                        use env_kindtypes,           only: si, si_x, dp
                        use obj_trajectory,          only: traj_training_data
                        use abs_obj_model,           only: model
                        use routines_modelSolve,     only: smp_solve_models

                        implicit none

                        interface
                                pure function omp_get_num_procs()
                                        import si_x
                                        integer (si_x) :: omp_get_num_procs
                                endfunction
                        end interface

                        class(model),  allocatable,intent(  out) :: result_models(:)
                        real    (dp),  allocatable,intent(  out) :: optimism_dist(:)

                        class(model),              intent(in   ) :: this_model
                        class(traj_training_data), intent(in   ) :: tdata
                        integer (si),              intent(in   ) :: nBSReps

                        real    (dp),             allocatable    :: model_residuals(:), native_residual

                        class(traj_training_data),allocatable    :: bs_tdatas(:)
                        class(model),             allocatable    :: bs_models(:), models(:)

                        integer (si)                    :: default_chunk_size, chunk_size
                        integer (si)                    :: jobs_remaining, current_pos

                        integer (si)                    :: model_iter, prep_iter

                        allocate(optimism_dist(nBSReps))

                        !to avoid using extensive memory
                        default_chunk_size = omp_get_num_procs()*16

                        jobs_remaining = nBSReps
                        current_pos = 1

                        allocate(models(1),source=this_model)

                        allocate(result_models(nBSreps),mold=this_model)

                        allocate(bs_tdatas(default_chunk_size),mold=tdata)

                        job_loop: do while (.true.)

                                chunk_size    = min(default_chunk_size,jobs_remaining)

                                call sgen_bootstrapped_trajs_omp(bs_tdatas(1:chunk_size)%trj,&
                                                                 tdata%trj)

                                !$omp parallel private (prep_iter)
                                !$omp do
                                do prep_iter=1,chunk_size
                                        call bs_tdatas(prep_iter)%preprocess(tdata%preproc_design)
                                enddo
                                !$omp end do
                                !$omp end parallel

                                call smp_solve_models(result_models = bs_models, &
                                                      models        = models, &
                                                     training_tdatas= bs_tdatas(1:chunk_size),  &
                                                      sort_by       = '', &
                                                      dedup         = .false.)

                                !$omp parallel do private(model_residuals,native_residual)
                                do model_iter=1,chunk_size
                                        native_residual = bs_models(model_iter)%getResidual()
                                        model_residuals = bs_models(model_iter)%recalcResidual(tdata)
                                        optimism_dist(current_pos + model_iter - 1) = native_residual - model_residuals(1)
                                enddo
                                !$omp end parallel do

                                !copy over trajs to output
                                !$omp parallel do
                                do model_iter=1,chunk_size
                                        result_models(current_pos + model_iter - 1) = bs_models(model_iter)
                                enddo
                                !$omp end parallel do

                                jobs_remaining = jobs_remaining - chunk_size
                                current_pos    = current_pos + chunk_size
                                if (jobs_remaining <= 0) then
                                        exit job_loop
                                endif

                        enddo job_loop

                endsubroutine solve_bootstrapped_models_smp

endmodule routines_bootstrap
