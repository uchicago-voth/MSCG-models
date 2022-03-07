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
! Describes routines needed for designing a set of models-- these are then 
! sent to be solved elsewhere.
!
! NOTE:
!    This could be extended to crate models with varying settings, though this isn't taken
!    advantage of now. It's harder to adjust it use a varity of types models, as then
!    we have to convert from an array of models to an array of containers holding models due to 
!    F limitations on polymorphic arrays.

module routines_modelDesign

        implicit none 
        private

        public createModelSet, populateModelSetNorms

        contains
                !Produces a set of configured models. Multiple invocation styles.
                !
                subroutine createModelSet(models,nModels,modelType,model_design)
                        use env_kindtypes,      only: si
                        use abs_obj_model,      only: model, model_config

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),intent(  out),allocatable  :: models(:)
                        character(*),intent(in   )              :: modelType
                        integer (si),intent(in   )              :: nModels
                        type(model_config),intent(in   )        :: model_design
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                        select case (modelType)
                        case ("linearDivision")
                                call gen_linear_models(models,nModels,model_design)
                        case ("centroid")
                                call gen_centroid_models(models,nModels,model_design)
                        case ("spectral")
                                call gen_spectral_models(models,nModels,model_design)
                        case default
                                print*, "Unknown model type in model creation &
                                &(this is an internal error; contact the developer)."
                                print*, "    attempted model type:"//modelType
                                stop
                        end select

                endsubroutine createModelSet

                subroutine gen_linear_models(models,nModels,model_design)

                        use env_kindtypes,              only: si
                        use abs_obj_model,              only: model, model_config
                        use obj_linearDivisionModel,    only: linearDivisionModel

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),intent(inout),allocatable  :: models(:)
                        integer (si),intent(in   )              :: nModels
                        type(model_config),intent(in   )        :: model_design
                        !!! End Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                            :: set

                        allocate(linearDivisionModel::models(nModels))

                        if (nModels == 0) return

                        select type (models)
                        type is (linearDivisionModel)
                        !$omp parallel private (set)
                        !$omp do
                        do set=1,nModels
                                call models(set)%configure(model_design)
                        enddo
                        !$omp end do
                        !$omp end parallel
                        end select

                endsubroutine gen_linear_models

                subroutine gen_centroid_models(models,nModels,model_design)

                        use env_kindtypes,        only: si
                        use abs_obj_model,        only: model, model_config
                        use obj_centroidModel,    only: centroidModel

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),intent(inout),allocatable  :: models(:)
                        integer (si),intent(in   )              :: nModels
                        type(model_config),intent(in   )        :: model_design
                        !!! End Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                            :: set

                        allocate(centroidModel::models(nModels))

                        if (nModels == 0) return

                        select type (models)
                        type is (centroidModel)
                        !$omp parallel private (set)
                        !$omp do
                        do set=1,nModels
                                call models(set)%configure(model_design)
                        enddo
                        !$omp end do
                        !$omp end parallel
                        end select

                endsubroutine gen_centroid_models

                subroutine gen_spectral_models(models,nModels,model_design)

                        use env_kindtypes,        only: si
                        use abs_obj_model,        only: model, model_config
                        use obj_spectralModel,    only: spectralModel

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),intent(inout),allocatable  :: models(:)
                        integer (si),intent(in   )              :: nModels
                        type(model_config),intent(in   )        :: model_design
                        !!! End Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                            :: set

                        allocate(spectralModel::models(nModels))

                        if (nModels == 0) return

                        select type (models)
                        type is (spectralModel)
                        !$omp parallel private (set)
                        !$omp do
                        do set=1,nModels
                                call models(set)%configure(model_design)
                        enddo
                        !$omp end do
                        !$omp end parallel
                        end select

                endsubroutine gen_spectral_models

                subroutine populateModelSetNorms(models,tdata,n_reps)

                        use env_kindtypes,        only: si
                        use abs_obj_model,        only: model
                        use obj_trajectory,       only: traj_training_data

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),             intent(inout)   :: models(:)
                        class(traj_training_data),intent(in   )   :: tdata
                        integer (si),             intent(in   )   :: n_reps
                        !!! End Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                              :: iter
                        logical, allocatable                      :: statuses(:)

                        if (size(models) < 1) return

                        allocate(statuses(size(models)))

                        if (n_reps == 0) then
                                call models(1)%populateNorms(tdata=tdata,deactivate_norms=.true.)
                        else
                                call models(1)%populateNorms(tdata=tdata,nSamples=n_reps)
                        endif

                        statuses(1) = .true.

                        if (size(models) > 1) then
                                do iter=2,size(models)
                                        call models(1)%copyNorms(target=models(iter),status=statuses(iter))
                                enddo
                        endif

                        if (.not. all(statuses)) then
                                print*, "Problem normalizing models. Stopping."
                                stop
                        endif

                endsubroutine populateModelSetNorms

endmodule routines_modelDesign
