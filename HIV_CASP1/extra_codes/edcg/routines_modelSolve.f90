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
!
! Describes routines needed for computing/optimizing a set of models. Future compatible for coarrays.
!

module routines_modelSolve

        implicit none

        private

        public solve_models, smp_solve_models

        contains
                !Produces a set of optimized models. Multiple invocation styles.
                !
                ! @models:         set of models to solve. 
                ! @result_count:   number top models to return (when parallel, we will want to filter)
                ! @compute_images: coarray indices to run with.
                ! @best_models:    array of top models.
                subroutine solve_models(models, result_count, result_models, &
                                        training_tdatas, dedup, sort_by, nCAranks)

                        use env_kindtypes,      only: si
                        use abs_obj_model,      only: model
                        use obj_trajectory,     only: traj_training_data

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),             intent(in   )             :: models(:)
                        integer (si),             intent(in   )             :: result_count
                        class(model),             intent(inout),allocatable :: result_models(:)
                        class(traj_training_data),intent(in   )             :: training_tdatas(:)
                        logical,                  intent(in   )             :: dedup
                        character(*),             intent(in   )             :: sort_by
                        integer (si),             intent(in   )             :: nCAranks
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        if (size(models) == 0) then
                                allocate(result_models(0),source=models)
                                return
                        endif

                        if (nCAranks == 1) then
                                 call smp_solve_models(result_models, models, result_count, &
                                                       training_tdatas, dedup, sort_by)
                        endif

                endsubroutine solve_models

                !Produces a array of optimized models. Every possible model and trajectory combindation
                ! is run, and if dedup is specified, all results sorted and deduplicated together to make 
                ! the results; this deduplication may or may not make sense in your application 
                ! (and may be expensive or blow the stack!)
                !
                ! @result_models:         Array of top models (if you filter, this holds the filtered results)
                ! @models:                Set of models to solve (these are solved in place)
                ! @result_count:          Number top models to return (when parallel, we will want to filter)
                ! @training_trajectories: Array of trajectories to train models on.
                ! @dedup:                 Deduplicate (and sort) result_models.

                !!! This can't be a function since you cannot assign to an allocatable polymorphic object !!!
                !!! in gfortran yet. !!!
                subroutine smp_solve_models(result_models, models, result_count, training_tdatas, dedup,&
                                            sort_by)

                        use env_kindtypes,      only: si, dp
                        use core_sort,          only: indexQsort
                        use obj_trajectory,     only: traj_training_data
                        use abs_obj_model,      only: model

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),             intent(  out),allocatable :: result_models(:)
                        class(model),             intent(in   )             :: models(:)
                        integer (si),             intent(in   ),optional    :: result_count
                        class(traj_training_data),intent(in   )             :: training_tdatas(:)
                        logical,                  intent(in   )             :: dedup
                        character(*),             intent(in   )             :: sort_by
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        class(model),allocatable :: sorted_models(:), expanded_models(:)

                        real    (dp),allocatable :: sort_values(:)
                        integer (si),allocatable :: bestModelIndices(:)
                        integer (si)             :: set, proxy_index, result_count_


                        !We need a model array for every model trajectory combination.
                        allocate(expanded_models(size(models)*size(training_tdatas)),source=models(1))

                        !Populate the extended set of models.

                        !$omp parallel private (proxy_index)
                        !$omp do
                        do set=1,(size(expanded_models))
                                proxy_index = floor(real(set-1)/size(training_tdatas))+1
                                expanded_models(set) = models(proxy_index)
                        enddo
                        !$omp end do
                        !$omp end parallel


                        !If we are going to deduplicate, record sort_values as we optimize the models.
                        if (dedup) allocate(sort_values(size(expanded_models)))

                        !$omp parallel private (proxy_index)
                        !$omp do
                        do set=1,size(expanded_models)

                                proxy_index = mod(set-1,size(training_tdatas))+1

                                call expanded_models(set)%initialize()
                                call expanded_models(set)%anneal(training_tdatas(proxy_index))
                                call expanded_models(set)%refine(training_tdatas(proxy_index))

                                if (dedup) then
                                        sort_values(set) = expanded_models(set)%getResidual()
                                endif
                        enddo
                        !$omp end do
                        !$omp end parallel

                        if (dedup) then

                                if (present(result_count)) then
                                        result_count_ = result_count
                                else
                                        result_count_ = size(sort_values)
                                endif

                                call indexQsort(sort_values,bestModelIndices,result_count_)

                                allocate(result_models(result_count_),source=expanded_models(1))

                                do set=1,result_count_
                                        result_models(set) = expanded_models(bestModelIndices(set))
                                enddo

                                call deduplicate_models(result_models,sorted_models)

                                deallocate(result_models)
                                allocate(result_models(size(sorted_models)),source=sorted_models(1))

                                select case (sort_by)
                                case ("residual","")
                                        ! we are already sorted by this. simply copy the models.
                                        do set=1,size(result_models)
                                                result_models(set) = sorted_models(set)
                                        enddo
                                case ("frequency")

                                        ! get corresponding frequencies, create a permutation, apply.

                                        ! gfortran bug. allocate on assignment fails here.
                                        deallocate(sort_values)
                                        allocate(sort_values(size(sorted_models)))
                                        sort_values = - sorted_models%getModelFreq()

                                        deallocate(bestModelIndices)
                                        allocate(bestModelIndices(size(sort_values)))

                                        call indexQsort(sort_values,bestModelIndices)

                                        do set=1,size(sorted_models)
                                                result_models(set) = sorted_models(bestModelIndices(set))
                                        enddo
                                case default
                                        print*, "WARNING: Unknown model sorting parameter. Defaulting &
                                                &to residual."
                                        do set=1,size(result_models)
                                                result_models(set) = sorted_models(set)
                                        enddo
                                end select
                        else 
                                allocate(result_models(size(expanded_models)),source=expanded_models(1))

                                do set=1,size(expanded_models)
                                        result_models(set) = expanded_models(set)
                                enddo
                        endif
                        
                endsubroutine smp_solve_models

                ! @models:         set of models to solve. 
                ! @result_count:   number top models to return (when parallel, we will want to filter)
                ! @compute_images: coarray indices to run with.
                ! @best_models:    array of top models.
                !subroutine coarray_solve_models&
                !                   (models, result_count, compute_images, best_models, training_tdatas)

                !        use env_kindtypes,      only: si, dp
                !        use obj_trajectory,     only: traj_training_data
                !        use abs_obj_model,      only: model

                !        implicit none

                !        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !        class(model),             intent(in   )             :: models(:)
                !        integer (si),             intent(in   )             :: result_count
                !        integer (si),             intent(in   ),optional    :: compute_images(:)
                !        class(model),             intent(inout),allocatable :: best_models(:)
                !        class(traj_training_data),intent(in   )             :: training_tdatas
                !        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !        !not implemented
                !        stop

                !endsubroutine coarray_solve_models

                pure subroutine deduplicate_models(models,target_models)

                        use env_kindtypes,      only: si
                        use abs_obj_model,      only: model

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),intent(in   )             :: models(:)
                        class(model),intent(inout),allocatable :: target_models(:)
                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)             :: iter

                        integer (si),allocatable :: buffer(:)
                        integer (si)             :: bufferEnd

                        allocate(buffer(size(models)))
                        bufferEnd = 1

                        buffer(1) = 1

                        do iter=2,size(models)
                                if (.not.( models(iter) == models(buffer(bufferEnd)) )) then
                                        bufferEnd = bufferEnd + 1
                                        buffer(bufferEnd) = iter
                                endif
                        enddo

                        allocate(target_models(bufferEnd),source=models(1))

                        do iter=1,bufferEnd
                                target_models(iter) = models(buffer(iter))
                        enddo

                        do iter=1,(bufferEnd-1)
                                call target_models(iter)%incrModelFreq(buffer(iter+1)-buffer(iter)-1)
                        enddo
                        call target_models(bufferEnd)%incrModelFreq(size(models)-buffer(bufferEnd))

                endsubroutine deduplicate_models

endmodule routines_modelSolve
