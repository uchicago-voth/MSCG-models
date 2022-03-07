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
! Provides _abstract_ definition of a model. Other modules define
! impementations, and provides the model_config object, which contains
! values used to parameterize a model. 
!
! NOTE:
!   Config objects reduce pain when rapidly prototyping what information is passed
!   during configuration.

module abs_obj_model

        use env_kindtypes, only: si, sp, dp
        use obj_trajectory,only: traj, traj_training_data

        implicit none

        private

        public          :: model, model_config, write_model, model_distance, write_model_stats

        interface write_model
                  module procedure write_model_array_filename
        end interface

        !This unifies the options used in configuring models. 
        !When possible, keep parameters outside this struct.
        type                           :: model_config

                !number of coarse-grained sites in model
                integer (si)    :: nSites = 1
                !ratio of residual types used when creating the total residual
                real    (dp),allocatable :: residual_weights(:)
                !max number of steps used when refining (e.g. steepest descent)
                integer (si)    :: refine_steps = 100
                !Initial acceptance probability of simulated annealing,
                !and control parameter for descent
                real    (dp)    :: anneal_accept = 1, anneal_param = 1
                !max number of simulated annealing steps
                integer (si)    :: anneal_steps = 10000
                !whether to log results for the model
                logical         :: record = .false.
                !Size of buffers in loggers for the model
                !these parameters aren't visible on the command line yet.
                integer (si)    :: size_buffer_refine = 50, size_buffer_anneal = 50
                !tolerances of how spatially distributed a cluster can be
                !(used in spatial centroid annealing)
                real    (dp)    :: cluster_density_tol = .05_dp

                integer (si)    :: kmeans_max_iter = 100
                integer (si)    :: kmeans_max_try  = 100

                !variables controling spectral models

                character(:),allocatable :: spectral_regularization_type

                !spectral edcg variables

                logical                  :: spectral_use_edcg_metric = .false.

                character(:),allocatable :: spectral_edcg_dist_conv_type

                character(:),allocatable :: spectral_edcg_kNN_type
                integer (si)             :: spectral_edcg_kNN = -1

                real    (dp)             :: spectral_edcg_dist_conv_param = 0

                !spectral spatial variables

                logical                  :: spectral_use_spatial_metric = .false.

                character(:),allocatable :: spectral_spatial_dist_conv_type

                character(:),allocatable :: spectral_spatial_kNN_type
                integer (si)             :: spectral_spatial_kNN = -1

                real    (dp)             :: spectral_spatial_dist_conv_param = 0

                !spectral pairVar variables

                logical                  :: spectral_use_pairVar_metric = .false.

                character(:),allocatable :: spectral_pairVar_dist_conv_type

                character(:),allocatable :: spectral_pairVar_kNN_type
                integer (si)             :: spectral_pairVar_kNN = -1

                real    (dp)             :: spectral_pairVar_dist_conv_param = 0

                !spectral charge variables

                logical                  :: spectral_use_charge_metric = .false.

                character(:),allocatable :: spectral_charge_dist_conv_type

                character(:),allocatable :: spectral_charge_kNN_type
                integer (si)             :: spectral_charge_kNN = -1

                real    (dp)             :: spectral_charge_dist_conv_param = 0

                !charge bead splitting
                logical                  :: mapping_charge_split = .false.

                !type of accumulation to use (e.g. center of mass) when using the map.
                character(:),allocatable :: map_accumulator_type

        end type model_config

        type, abstract                 :: model

                private

                real    (dp),allocatable :: normOffsets(:)
                real    (dp),allocatable :: normScalings(:)

                contains
                        procedure(abs_self),           deferred         :: reset

                        procedure(abs_anneal),         deferred         :: anneal
                        procedure(abs_refine),         deferred         :: refine

                        !recalculate residual given a new trajectory (no optimization)
                        procedure(abs_recalc),         deferred         :: recalcResidual

                        !prints the model description to a file
                        procedure(abs_self_write    ), deferred         :: summarize
                        procedure(abs_getInt_si ),     deferred         :: getNSites
                        procedure(abs_get_mapping),    deferred         :: getMapping

                        procedure(abs_is_map_charge_split),deferred     :: isMapChargeSplit
                        procedure(abs_has_parent_map), deferred         :: hasParentMap

                        procedure(abs_getArrayReal_dp),deferred         :: getSiteCharges
                        procedure(abs_getResidual),    deferred         :: getResidual

                        procedure(abs_tainted),        deferred         :: isTainted

                        procedure(abs_setTainted),     deferred         :: setTaint

                        procedure(abs_init),           deferred         :: initialize
                        procedure(abs_configure),      deferred         :: configure

                        !models which are identical are often coalesced;  these control 
                        !the number of occurances a model represents.
                        procedure(abs_e_getInt_si),    deferred         :: getModelFreq
                        procedure(abs_self_opt_si),    deferred         :: incrModelFreq
                        procedure(abs_self_opt_si),    deferred         :: decrModelFreq

                        !CI = confidence interval.
                        procedure(abs_init_ci),        deferred         :: initializeCI
                        procedure(abs_has_ci),         deferred         :: hasCI 

                        procedure(abs_calc_ci),        deferred         :: calcCI
                        procedure(abs_get_ci),         deferred         :: getCI

                        procedure(abs_addCIValues),    deferred         :: addCIValues

                        !generates a random mapping from the space of possible mappings for the model
                        procedure(abs_getRandomMapping),deferred,private :: getRandomMapping

                        procedure(abs_getMapAccumType), deferred        :: getMapAccumType

                        !translation applied to residuals.
                        procedure                                       :: getNormOffsets
                        !scaling applied to residuals (after translation)
                        procedure                                       :: getNormScalings
                        procedure                                       :: populateNorms
                        procedure                                       :: copyNorms

                        procedure,                              private :: eq    => generic_model_eq
                        procedure,                              private :: clone => generic_model_clone
                        generic                                         :: assignment(=) => clone
                        generic                                         :: operator(==)  => eq
        end type model

        abstract interface
            pure subroutine abs_self(self)

                import :: model
                implicit none

                class(model), intent(inout) :: self

            endsubroutine abs_self
        end interface

        abstract interface
            pure function abs_has_ci(self,bca)

                import :: model
                implicit none

                class(model),        intent(in   ) :: self
                logical,optional,    intent(in   ) :: bca

                logical                     :: abs_has_ci

            endfunction abs_has_ci
        end interface

        abstract interface
            subroutine abs_calc_ci(self,bca)

                import :: model
                implicit none

                class(model),         intent(inout) :: self
                logical,     optional,intent(in   ) :: bca

            endsubroutine abs_calc_ci
        end interface

        abstract interface
            pure subroutine abs_init_ci(self,dist_size,jk)

                import :: model,si
                implicit none

                class(model),          intent(inout) :: self
                integer (si),          intent(in   ) :: dist_size
                logical,     optional, intent(in   ) :: jk

            endsubroutine abs_init_ci
        end interface

        abstract interface
            function abs_get_ci(self,percentile,bca)

                import :: model,si,dp
                implicit none

                class(model),          intent(in   ) :: self
                real    (dp),          intent(in   ) :: percentile
                logical,     optional, intent(in   ) :: bca

                real    (dp)                         :: abs_get_ci

            endfunction abs_get_ci
        end interface

        !abstract interface
        !    pure function abs_get_lower_ci(self,percentile,bca)

        !        import :: model,si,dp
        !        implicit none

        !        class(model),          intent(in   ) :: self
        !        integer (si),          intent(in   ) :: percentile
        !        logical,     optional, intent(in   ) :: bca


        !        real    (dp)                         :: abs_get_lower_ci

        !    endfunction abs_get_lower_ci
        !end interface

        abstract interface
            pure subroutine abs_init(self)

                import :: model
                implicit none

                class(model),          intent(inout) :: self

            endsubroutine abs_init
        end interface

        abstract interface
            pure subroutine abs_configure(self,config)

                import :: model,model_config
                implicit none

                class(model),          intent(inout) :: self
                type(model_config),    intent(in   ) :: config

            endsubroutine abs_configure
        end interface

        abstract interface
            subroutine abs_anneal(self,tdata,acceptance,descentControl)

                import :: model,traj_training_data,dp
                implicit none

                class(model),                      intent(inout) :: self
                class(traj_training_data),         intent(in   ) :: tdata
                real (dp),                optional,intent(in   ) :: acceptance,descentControl

            endsubroutine abs_anneal
        end interface

        abstract interface
            subroutine abs_refine(self,tdata,maxSteps)

                import :: model,traj_training_data,si
                implicit none

                class(model),                      intent(inout) :: self
                class(traj_training_data),         intent(in   ) :: tdata
                integer (si),             optional,intent(in   ) :: maxSteps

            endsubroutine abs_refine
        end interface

        abstract interface
            function abs_recalc(self,tdata)

                import :: model,traj_training_data,dp
                implicit none

                class(model),             intent(inout) :: self
                class(traj_training_data),intent(in   ) :: tdata
                
                real    (dp), allocatable           :: abs_recalc(:)

            endfunction abs_recalc
        end interface

        abstract interface
            subroutine abs_self_traj(self,tdata)

                import :: model,traj_training_data
                implicit none

                class(model),             intent(inout) :: self
                class(traj_training_data),intent(in   ) :: tdata

            endsubroutine abs_self_traj
        end interface

        abstract interface
            subroutine abs_self_write(self,fileID,filename,verbose)

                import :: model,si
                implicit none

                class(model), intent(in   )          :: self
                integer (si), intent(in   ),optional :: fileID
                character (len=*),&
                              intent(in   ),optional :: filename
                logical,      intent(in   ),optional :: verbose

            endsubroutine abs_self_write
        end interface

        abstract interface
            pure function abs_getArrayReal_dp(self)

                import :: model,dp
                implicit none

                class(model), intent(in   ) :: self
                real    (dp), allocatable   :: abs_getArrayReal_dp(:)

            endfunction abs_getArrayReal_dp
        end interface

        abstract interface
            pure function abs_getArrayReal_sp(self)

                import :: model,sp
                implicit none

                class(model), intent(in   ) :: self
                real    (sp), allocatable   :: abs_getArrayReal_sp(:)

            endfunction abs_getArrayReal_sp
        end interface

        abstract interface
            pure subroutine abs_setArrayReal_sp(self,toSet)

                import :: model,sp
                implicit none

                class(model), intent(inout) :: self
                real    (sp), intent(in   ) :: toSet(:)

            endsubroutine abs_setArrayReal_sp
        end interface

        abstract interface
            pure subroutine abs_setArrayReal_dp(self,toSet)

                import :: model,dp
                implicit none

                class(model), intent(inout) :: self
                real    (dp), intent(in   ) :: toSet(:)

            endsubroutine abs_setArrayReal_dp
        end interface

        abstract interface
            pure function abs_getArrayInt_si(self)

                import :: model,si
                implicit none

                class(model), intent(in   ) :: self
                integer (si), allocatable   :: abs_getArrayInt_si(:)

            endfunction abs_getArrayInt_si
        end interface

        abstract interface
            pure function abs_get_mapping(self,parent_mapping,split_on_charge)

                import :: model,si
                implicit none

                class(model),          intent(in   ) :: self
                logical,      optional,intent(in   ) :: parent_mapping
                logical,      optional,intent(in   ) :: split_on_charge

                integer (si), allocatable   :: abs_get_mapping(:)

            endfunction abs_get_mapping
        end interface

        abstract interface
            pure function abs_has_parent_map(self)

                import :: model

                implicit none

                class(model),        intent(in   ) :: self

                logical                     :: abs_has_parent_map

            endfunction abs_has_parent_map
        end interface

        abstract interface
            pure function abs_is_map_charge_split(self)

                import :: model
                implicit none

                class(model),        intent(in   ) :: self

                logical                     :: abs_is_map_charge_split

            endfunction abs_is_map_charge_split
        end interface


        abstract interface
            pure subroutine abs_setArrayInt_si(self,toSet)

                import :: model,si
                implicit none

                class(model), intent(inout) :: self
                integer (si), intent(in   ) :: toSet(:)

            endsubroutine abs_setArrayInt_si
        end interface


        abstract interface
            pure function abs_getReal_dp(self)

                import :: model,dp
                implicit none

                class(model), intent(in   ) :: self
                real    (dp)                :: abs_getReal_dp

            endfunction abs_getReal_dp
        end interface

        abstract interface
            elemental function abs_getresidual(self)

                import :: model,dp
                implicit none

                class(model), intent(in   ) :: self
                real    (dp)                :: abs_getresidual

            endfunction abs_getResidual
        end interface

        abstract interface
            pure function abs_getInt_si(self)

                import :: model,si
                implicit none

                class(model), intent(in   ) :: self
                integer (si)                :: abs_getInt_si

            endfunction abs_getInt_si
        end interface

        abstract interface
            elemental function abs_e_getInt_si(self)

                import :: model,si
                implicit none

                class(model), intent(in   ) :: self
                integer (si)                :: abs_e_getInt_si

            endfunction abs_e_getInt_si
        end interface

        abstract interface
            pure subroutine abs_setInt_si(self,toSet)

                import :: model,si
                implicit none

                class(model), intent(inout) :: self
                integer (si), intent(in   ) :: toSet

            endsubroutine abs_setInt_si
        end interface

        abstract interface
            pure function abs_tainted(self)

                import :: model
                implicit none

                class(model), intent(in   ) :: self
                logical         :: abs_tainted

            endfunction abs_tainted
        end interface

        abstract interface
            pure subroutine abs_setTainted(self,toSet)

                import :: model
                implicit none

                class(model), intent(inout) :: self
                logical,      intent(in   ) :: toSet

            endsubroutine abs_setTainted
        end interface

        abstract interface
            pure function abs_getString(self)

                import :: model
                implicit none

                class(model), intent(in   ) :: self
                character(:),allocatable    :: abs_getString

            end function abs_getString
        end interface

        abstract interface
            pure subroutine abs_self_opt_si(self,amount)

                import :: model, si
                implicit none

                class(model), intent(inout)          :: self
                integer (si), intent(in   ),optional :: amount

            endsubroutine abs_self_opt_si
        end interface

        abstract interface
            pure subroutine abs_addArrayReal_dp(self,toAdd)

                import :: model,dp
                implicit none

                class(model), intent(inout)  :: self
                real    (dp), intent(in   )  :: toAdd(:)

            endsubroutine abs_addArrayReal_dp
        end interface

        abstract interface
            pure subroutine abs_addCIValues(self,toAdd,jk)

                import :: model,dp
                implicit none

                class(model),          intent(inout) :: self
                real    (dp),          intent(in   ) :: toAdd(:)
                logical,     optional, intent(in   ) :: jk

            endsubroutine abs_addCIValues
        end interface

        abstract interface
            function abs_getRandomMapping(self,tdata)

                import :: model,traj_training_data,si
                implicit none

                class(model),             intent(in   ) :: self
                class(traj_training_data),intent(in   ) :: tdata

                integer (si),allocatable :: abs_getRandomMapping(:)

            endfunction abs_getRandomMapping
        end interface

        abstract interface
            pure function abs_getMapAccumType(self)

                import :: model
                implicit none

                class(model),intent(in   )      :: self

                character(:),allocatable        :: abs_getMapAccumType

            endfunction abs_getMapAccumType
        end interface

        contains
                !If this isn't overridden, no copying occurs!
                pure subroutine generic_model_clone(lhs,rhs)

                        implicit none

                        class(model),intent(inout) :: lhs
                        class(model),intent(in   ) :: rhs

                        !This is actualyl to just stop warnings. This should never get called.
                        call lhs%setTaint(rhs%isTainted())
                        call lhs%setTaint(.true.)

                endsubroutine generic_model_clone

                !If this isn't overridden, equality is always false.
                pure function generic_model_eq(lhs,rhs)

                        implicit none

                        class(model),intent(in   ) :: lhs
                        class(model),intent(in   ) :: rhs

                        logical generic_model_eq

                        !dummy statement to stop compiler warnings about unused variable.
                        !forced by polymorphish and current structure.
                        if (allocated(lhs%normOffsets)) CONTINUE
                        if (allocated(rhs%normOffsets)) CONTINUE

                        generic_model_eq = .false.

                endfunction generic_model_eq

                subroutine write_model_array_fileid(models,fileid,verbose)

                        implicit none

                        class(model),         intent(in   ) :: models(:)
                        integer (si),         intent(in   ) :: fileid
                        logical,     optional,intent(in   ) :: verbose

                        logical                             :: verbose_
                        integer (si)                        :: set

                        if (present(verbose)) then
                                verbose_ = verbose
                        else
                                verbose_ = .false.
                        endif

                        do set=1,size(models)
                                call models(set)%summarize(fileID=fileid,verbose=verbose_)
                        enddo

                endsubroutine write_model_array_fileid

                subroutine write_model_array_filename(models,filename,verbose)

                        implicit none

                        class(model),     intent(in   ) :: models(:)
                        character (len=*),intent(in   ) :: filename
                        logical, optional,intent(in   ) :: verbose

                        logical                         :: verbose_
                        integer (si)                    :: file_handle, status

                        if (present(verbose)) then
                                verbose_ = verbose
                        else
                                verbose_ = .false.
                        endif

                        open(newunit = file_handle, file = filename, iostat = status)
                        if (status /= 0) then
                                print*, "Model couldn't be written."
                                return
                        endif

                        call write_model_array_fileid(models,file_handle,verbose_)

                        close(file_handle)

                endsubroutine write_model_array_filename

                !computes the distance between two models.
                pure function model_distance(model1,model2)

                        use core_stat,          only: assignment_difference

                        implicit none

                        class(model),     intent(in   ) :: model1, model2

                        real    (dp)                    :: model_distance

                        integer (si),allocatable        :: mapping1(:), mapping2(:)

                        allocate(mapping1,source=model1%getMapping())
                        allocate(mapping2,source=model2%getMapping())

                        model_distance = assignment_difference(mapping1,mapping2)

                endfunction model_distance

                !populates the normal values used in residual normalizations
                subroutine populateNorms(self,tdata,nSamples,deactivate_norms)

                        use obj_trajectory,     only: traj_training_data
                        use obj_accum,          only: varAccum
                        use obj_ll,             only: i_sp_dll
                        use fit_common,         only: fit_residuals, num_residuals

                        implicit none

                        class(model),             intent(inout) :: self
                        class(traj_training_data),intent(in   ) :: tdata
                        integer (si),  optional,  intent(in   ) :: nSamples
                        logical,       optional,  intent(in   ) :: deactivate_norms

                        !to hold accumulated residuals.
                        type(varAccum),allocatable      :: accums(:)

                        type(i_sp_dll),allocatable      :: backmapping(:)

                        integer (si)                    :: sample_iter,res_iter,nSamples_
                        integer (si),allocatable        :: mapping(:)
                        real    (dp),allocatable        :: residualList(:), cgCharges(:)
                        real    (dp),allocatable        :: zero_array(:), one_array(:)
                        logical                         :: deactivate_norms_

                        real    (dp),allocatable        :: equal_weights_array(:)

                        if (present(deactivate_norms)) then
                                deactivate_norms_ = deactivate_norms
                        else
                                deactivate_norms_ = .false.
                        endif

                        if (deactivate_norms_) then
                                allocate(self%normOffsets(num_residuals))
                                allocate(self%normScalings(num_residuals))
                                self%normOffsets  = 0
                                self%normScalings = 1
                        else
                                if (present(nSamples)) then
                                        nSamples_ = nSamples
                                else
                                        nSamples_ = 100
                                endif

                                allocate(residualList(num_residuals+1))
                                allocate(accums(num_residuals))
                                allocate(cgCharges(self%getNSites()))
                                allocate(backmapping(self%getNSites()))

                                allocate(zero_array(num_residuals))
                                allocate(one_array(num_residuals))
                                allocate(equal_weights_array(num_residuals))
                                zero_array = 0
                                one_array  = 1
                                equal_weights_array = (1.0_dp/3.0_dp)

                                call accums%reset()

                                do sample_iter=1,nSamples_
                                        mapping = self%getRandomMapping(tdata)
                                        call fit_residuals(      tdata = tdata,&
                                                               mapping = mapping,&
                                              mapping_accumulator_type = self%getMapAccumType(),&
                                                      residual_weights = equal_weights_array,&
                                                             cgCharges = cgCharges,&
                                                           backmapping = backmapping,&
                                                          residualList = residualList,&
                                                     residual_offsets  = zero_array,&
                                                     residual_scalings = one_array,&
                                                                update = .false.)
                                        do res_iter=1,num_residuals
                                                call accums(res_iter)%add(residualList(res_iter+1))
                                        enddo
                                enddo

                                self%normOffsets  = - accums%getMean()
                                self%normScalings = 1.0_dp / accums%getSD()
                        endif

                endsubroutine populateNorms

                pure subroutine copyNorms(self,target,status)

                        implicit none

                        class(model),         intent(in   ) :: self
                        class(model),         intent(inout) :: target
                        logical,     optional,intent(  out) :: status

                        !to hold accumulated residuals.

                        if (present(status)) status = .true.

                        if (allocated(self%normOffsets)) then
                                target%normOffsets = self%normOffsets
                        else
                                if (present(status)) status = .false.
                        endif

                        if (allocated(self%normScalings)) then
                                target%normScalings = self%normScalings
                        else
                                if (present(status)) status = .false.
                        endif
                        
                endsubroutine copyNorms

                pure function getNormOffsets(self)

                        implicit none

                        class(model),     intent(in   ) :: self

                        real    (dp), allocatable       :: getNormOffsets(:)

                        if (allocated(self%normOffsets)) then
                                getNormOffsets = self%normOffsets
                        else
                                allocate(getNormOffsets(0))
                        endif

                endfunction getNormOffsets

                pure function getNormScalings(self)

                        implicit none

                        class(model),     intent(in   ) :: self

                        real    (dp), allocatable       :: getNormScalings(:)

                        if (allocated(self%normScalings)) then
                                getNormScalings = self%normScalings
                        else
                                allocate(getNormScalings(0))
                        endif

                endfunction getNormScalings

                subroutine write_model_stats(models,optimism_dist,report_file)

                        use env_kindtypes,           only: si, dp
                        use core_convert,            only: itoa

                        implicit none

                        class(model),              intent(in   ) :: models(:)
                        real    (dp), optional,    intent(in   ) :: optimism_dist(:)
                        character(*),              intent(in   ) :: report_file

                        real    (dp), allocatable       :: model_similarities(:)
                        integer (si)                    :: iter1, iter2, index, file_handle
                        integer (si)                    :: n_models_1, n_models_2, status

                        !guard for zero model case
                        if (size(models) == 0) return

                        open(newunit = file_handle, file = report_file, iostat=status)
                        if (status /= 0) then
                                print*, "Model summary couldn't be written."
                                return
                        endif

                        if (present(optimism_dist)) then
                                if (size(optimism_dist) > 1) then
                                        write(file_handle,*) 'Resampled optimism distribution:', optimism_dist

                                        write(file_handle,*) 'Optimism estimate:', &
                                                sum(optimism_dist)/max(size(optimism_dist),1)
                                endif
                        endif

                        !guard for zero model case
                        n_models_1 = sum(models%getModelFreq())

                        if (n_models_1 > 1) then
                                allocate(model_similarities((n_models_1**2 - n_models_1)/2))

                                index = 1
                                do iter1 = 1 , size(models)
                                        n_models_1 = models(iter1)%getModelFreq()

                                        model_similarities(index:(  index &
                                                                  + number_pairwise_comp(n_models_1) &
                                                                  - 1                              )) = 0

                                        index = index + number_pairwise_comp(n_models_1)
                                enddo

                                if (size(models) > 1) then
                                        do iter1 = 1 , size(models) - 1
                                        do iter2 = iter1 + 1, size(models)
                                                n_models_1 = models(iter1)%getModelFreq()
                                                n_models_2 = models(iter2)%getModelFreq()

                                                model_similarities(index:(index + n_models_1*n_models_2 - 1)) = &
                                                        model_distance(models(iter1),models(iter2))

                                                index = index + (n_models_1*n_models_2)
                                        enddo 
                                        enddo
                                endif

                                write(file_handle,*) 'Clustering stability distribution:', model_similarities
                                write(file_handle,*) 'Clustering stability estimate:',     &
                                                  sum(model_similarities)/max(size(model_similarities),1)
                         else
                                write(file_handle,*) itoa(size(models))//' model(s) given; &
                                        &no model comparisons performed or recorded.'
                         endif

                         close(file_handle)

                endsubroutine 

                !returns the number of unique pairwise comparisons BETWEEN set1 and set2
                !or just in set1.
                ! @size_1 : size of set 1
                ! @size_2 : size of set 2 (optional)
                pure function number_pairwise_comp(size_1,size_2)

                        implicit none

                        integer (si),         intent(in   ) :: size_1
                        integer (si),optional,intent(in   ) :: size_2

                        integer (si)                        :: number_pairwise_comp

                        if (present(size_2)) then
                                number_pairwise_comp = size_1 * size_2
                                return
                        else
                                number_pairwise_comp = size_1 * (size_1 - 1) /2
                        endif

                endfunction number_pairwise_comp


endmodule abs_obj_model
