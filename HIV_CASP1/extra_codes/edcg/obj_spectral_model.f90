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
! Provides definition of a spectral model, which operates using spectral clustering methods
! to create noncontiguous mappings.
!

module obj_spectralModel

        use env_kindtypes, only: si, sp, dp
        use abs_obj_model, only: model, model_config
        use fit_common,    only: num_residuals
        use obj_Rloggers,  only: t_rLogger
        use routines_spectralOptimization,&
                           only: spectral_dist_config

        implicit none

        private

        public          :: spectralModel

        type, extends(model) :: spectralModel
                private

                !fund. state variables
                logical                  :: tainted = .false.
                logical                  :: record  = .false.

                integer (si)             :: frequency = 1

                !mapping state variables
                integer (si)             :: nSites = 0
                integer (si)             :: nAtoms = 0   !input number of sites

                real    (dp),allocatable :: minimum_projection(:,:)
                real    (dp),allocatable :: minimum_eigenvalues(:)
                integer (si),allocatable :: minimum_mapping(:)

                !when dealing with a mapped trajectory, this gives a map from
                !the parent trajectory to the mapped trajectory, so that 
                !we can get the full resolution mapping in a simple manner.
                integer (si),allocatable :: parent_mapping(:)
                real    (dp),allocatable :: parent_charges(:)

                real    (dp),allocatable :: domain_charges(:)

                real    (dp),allocatable :: minimum_cg_charges(:)
                real    (dp),allocatable :: mapped_cg_charges(:)

                logical                  :: postproc_split_on_charge = .false.

                !residual vars
                real    (dp)             :: residuals(num_residuals+1)      = -1
                real    (dp)             :: residual_weights(num_residuals) = -1

                !metric variables-- these control which how distance metrics are used
                !in the spectral clustering.
                type(spectral_dist_config) :: edcg_config, spatial_config, pairVar_config, &
                                              charge_config

                !regression control variables
                character(:),allocatable :: regularization_type

                integer (si)             :: kmeans_discretize_max_iter
                integer (si)             :: kmeans_discretize_max_try

                !bootstrapping variables

                real    (dp),allocatable :: bootstrap_distribution(:)
                integer (si)             :: bootstrap_dist_current_mark = 0
                logical                  :: bootstrap_distribution_issorted =.false.

                real    (dp)             :: bca_bias,bca_accel

                real    (dp),allocatable :: jackknife_distribution(:)
                integer (si)             :: jackknife_dist_current_mark = 0
                logical                  :: jackknife_distribution_issorted =.false.

                !holds type of accumulation used in mapping (e.g. center of mass).
                character(:),allocatable :: map_accumulator_type

                contains
                        procedure,public :: reset         => sm_reset
                        procedure,public :: anneal        => sm_sim_annealing
                        procedure,public :: refine        => sm_lucky_descent
                        procedure,public :: summarize     => sm_summarize

                        procedure,public :: recalcResidual=> sm_recalc_residual

                        procedure,public :: getAnnealLog  => sm_get_anneal_log
                        procedure,public :: getRefineLog  => sm_get_refine_log

                        procedure,public :: getMapping    => sm_get_mapping

                        procedure,public :: isMapChargeSplit => sm_is_mapping_charge_split
                        procedure,public :: hasParentMap  => sm_has_parent_mapping

                        procedure        :: setAnnealLog  => sm_set_anneal_log
                        procedure        :: setRefineLog  => sm_set_refine_log
                        procedure,public :: getSiteCharges =>sm_get_siteCharges
                        procedure,public :: getResidual    =>sm_get_residual
                        procedure,public :: getFullResidual=>sm_get_full_residual
                        procedure,public :: initialize    => sm_initialize
                        procedure,public :: configure     => sm_configure

                        procedure,public :: isTainted     => sm_isTainted

                        procedure        :: setTaint      => sm_setTainted

                        procedure,public :: getModelFreq  => sm_get_model_freq
                        procedure,public :: incrModelFreq => sm_incr_model_freq
                        procedure,public :: decrModelFreq => sm_decr_model_freq

                        procedure,public :: initializeCI  => sm_init_ci

                        procedure,public :: calcCI        => sm_calc_ci
                        procedure,public :: hasCI         => sm_has_ci

                        procedure,public :: getCI         => sm_get_ci

                        procedure,public :: getBias       => sm_get_bias

                        procedure,public :: addCIValues   => sm_add_ci_values

                        procedure,public :: getRandomMapping => sm_get_random_mapping
                        procedure,public :: getNSites     => get_nSites

                        procedure,public :: getMapAccumType => sm_get_map_accum_type

                        procedure        :: clone         => sm_clone
                        procedure        :: eq            => sm_eq
        end type spectralModel

        contains
                !In the spectral case, only a subset of these parameters are used.
                pure subroutine sm_configure(self,config)

                        implicit none

                        class(spectralModel),      intent(inout) :: self
                        type(model_config),        intent(in   ) :: config

                        call self%setTaint(.false.)

                        self%map_accumulator_type = config%map_accumulator_type

                        if (allocated(config%residual_weights)) then
                                if (size(config%residual_weights,1) == &
                                    size(self%residual_weights))    then
                                        self%residual_weights = config%residual_weights
                                else
                                        self%residual_weights = -1
                                        call self%setTaint(.true.)
                                endif
                        else
                                self%residual_weights = -1
                                call self%setTaint(.true.)
                        endif

                        self%nSites  = config%nSites
                        self%record  = config%record

                        self%kmeans_discretize_max_iter = config%kmeans_max_iter
                        self%kmeans_discretize_max_try = config%kmeans_max_try

                        self%regularization_type  = config%spectral_regularization_type

                        if (config%spectral_use_edcg_metric) then
                                self%edcg_config%active = .true.
                                self%edcg_config%dist_conv_type  = config%spectral_edcg_dist_conv_type
                                self%edcg_config%dist_conv_param = config%spectral_edcg_dist_conv_param
                                self%edcg_config%kNN             = config%spectral_edcg_kNN
                                self%edcg_config%kNN_type        = config%spectral_edcg_kNN_type
                        endif

                        if (config%spectral_use_spatial_metric) then
                                self%spatial_config%active = .true.
                                self%spatial_config%dist_conv_type  = config%spectral_spatial_dist_conv_type
                                self%spatial_config%dist_conv_param = config%spectral_spatial_dist_conv_param
                                self%spatial_config%kNN             = config%spectral_spatial_kNN
                                self%spatial_config%kNN_type        = config%spectral_spatial_kNN_type
                        endif

                        if (config%spectral_use_pairVar_metric) then
                                self%pairVar_config%active = .true.
                                self%pairVar_config%dist_conv_type  = config%spectral_pairVar_dist_conv_type
                                self%pairVar_config%dist_conv_param = config%spectral_pairVar_dist_conv_param
                                self%pairVar_config%kNN             = config%spectral_pairVar_kNN
                                self%pairVar_config%kNN_type        = config%spectral_pairVar_kNN_type
                        endif

                        if (config%spectral_use_charge_metric) then
                                self%charge_config%active = .true.
                                self%charge_config%dist_conv_type  = config%spectral_charge_dist_conv_type
                                self%charge_config%dist_conv_param = config%spectral_charge_dist_conv_param
                                self%charge_config%kNN             = config%spectral_charge_kNN
                                self%charge_config%kNN_type        = config%spectral_charge_kNN_type
                        endif

                        self%postproc_split_on_charge = config%mapping_charge_split

                endsubroutine sm_configure

                pure subroutine sm_initialize(self)

                        implicit none

                        class(spectralModel),intent(inout) :: self

                        allocate(self%minimum_cg_charges(self%nSites))

                        !if (self%record) then
                        !        call self%sa_logger%initialize()
                        !endif

                endsubroutine sm_initialize

                pure subroutine sm_reset(self)

                        implicit none

                        class(spectralModel), intent(inout) :: self

                        self%residuals = -1

                        !call self%sa_logger%reset(shallow=.true.)

                endsubroutine sm_reset

                !spectral routines ignore acceptance and descent control.
                subroutine sm_sim_annealing(self,tdata,acceptance,descentControl)

                        use obj_trajectory,     only: traj_training_data,labeledTraj
                        use routines_spectralOptimization,&
                                                only: spectralOptimization
                        use core_filter,        only: valueCollapse
                        use fit_common,         only: fit_residuals
                        use obj_ll,             only: i_sp_dll

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel),              intent(inout) :: self
                        class(traj_training_data),         intent(in   ) :: tdata
                        real (dp),                optional,intent(in   ) :: acceptance,descentControl
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        type(i_sp_dll),allocatable :: backmapping(:)

                        !dummy statement to stop compiler warnings about unused variable.
                        !forced by polymorphish and current structure.
                        if (present(acceptance)) CONTINUE
                        if (present(descentControl)) CONTINUE

                        if (self%nAtoms == 0) then
                                self%nAtoms = tdata%trj%nAtoms
                        endif

                        self%domain_charges = tdata%trj%atomCharges

                        if (tdata%trj%isMapped()) then
                                self%parent_mapping = tdata%trj%getParentMap()
                                self%parent_charges = tdata%trj%getParentCharges()
                        endif

                        call spectralOptimization(        mapping = self%minimum_mapping,&
                                                       projection = self%minimum_projection, &
                                                      eigenvalues = self%minimum_eigenvalues,&
                                                              nCG = self%nSites,&
                                                            tdata = tdata,&
                                                   regularization = self%regularization_type,&
                                                  kmeans_max_iter = self%kmeans_discretize_max_iter,&
                                                   kmeans_max_try = self%kmeans_discretize_max_try,&
                                                      edcg_config = self%edcg_config,&
                                                   spatial_config = self%spatial_config,&
                                                   pairVar_config = self%pairVar_config,&
                                                    charge_config = self%charge_config)

                        self%mapped_cg_charges = valueCollapse(indexLabels = self%minimum_mapping,&
                                                               indexValues = tdata%trj%atomCharges,&
                                                                   nLabels = self%nSites,&
                                                              perValueNorm = .false.)

                        call fit_residuals(tdata            = tdata,    &
                                           mapping          = self%minimum_mapping,&
                                   mapping_accumulator_type = self%getMapAccumType(),&
                                           residual_weights = self%residual_weights,&
                                           cgCharges        = self%minimum_cg_charges,&
                                           backmapping      = backmapping,     &
                                           residualList     = self%residuals,  &
                                           residual_offsets = self%getNormOffsets(),&
                                           residual_scalings= self%getNormScalings())

                endsubroutine sm_sim_annealing

                subroutine sm_lucky_descent(self,tdata,maxSteps)

                        use obj_trajectory, only: traj_training_data

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel),         intent(inout) :: self
                        class(traj_training_data),    intent(in   ) :: tdata
                        integer (si),        optional,intent(in   ) :: maxSteps
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !dummy statement to stop compiler warnings about unused variable.
                        !forced by polymorphish and current structure.
                        if (maxSteps == 0) CONTINUE

                        self%domain_charges = tdata%trj%atomCharges

                        if (tdata%trj%isMapped()) then
                                self%parent_mapping = tdata%trj%getParentMap()
                                self%parent_charges = tdata%trj%getParentCharges()
                        endif

                endsubroutine sm_lucky_descent

                pure subroutine sm_set_anneal_log(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(inout)          :: self
                        real    (sp),               intent(in   )          :: toSet(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !dummy statement to stop compiler warnings about unused variable.
                        !forced by polymorphish and current structure.
                        if (self%tainted) CONTINUE
                        if (size(toSet,1) == 0) CONTINUE

                endsubroutine sm_set_anneal_log

                pure function sm_get_anneal_log(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (sp),allocatable  :: sm_get_anneal_log(:)

                        !dummy statement to stop compiler warnings about unused variable.
                        !forced by polymorphish and current structure.
                        if (self%tainted) CONTINUE

                        allocate(sm_get_anneal_log(0))

                endfunction sm_get_anneal_log

                pure subroutine sm_set_refine_log(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(inout)          :: self
                        real    (sp),               intent(in   )          :: toSet(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !dummy statement to stop compiler warnings about unused variable.
                        !forced by polymorphish and current structure.
                        if (self%tainted) CONTINUE
                        if (size(toSet,1) == 0) CONTINUE

                endsubroutine sm_set_refine_log

                pure function sm_get_refine_log(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (sp),allocatable  :: sm_get_refine_log(:)

                        !dummy statement to stop compiler warnings about unused variable.
                        !forced by polymorphish and current structure.
                        if (self%tainted) CONTINUE

                        allocate(sm_get_refine_log(0))

                endfunction sm_get_refine_log

                elemental function sm_get_residual(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp)                 :: sm_get_residual

                        sm_get_residual = self%residuals(1)

                endfunction sm_get_residual

                pure function sm_get_full_residual(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp),allocatable  :: sm_get_full_residual(:)

                        sm_get_full_residual = self%residuals

                endfunction sm_get_full_residual

                pure function sm_get_mapping(self,parent_mapping,split_on_charge)

                        use core_filter,    only: boundariesToMapping, composeMapping, &
                                                  fragmentMapping, signDiscretize

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel),         intent(in   )          :: self
                        logical,             optional,intent(in   )          :: parent_mapping
                        logical,             optional,intent(in   )          :: split_on_charge
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si),allocatable :: sm_get_mapping(:)
                        logical                  :: parent_mapping_
                        logical                  :: do_split_on_charge
                        integer (si),allocatable :: charge_mapping(:)

                        if (present(parent_mapping)) then
                                parent_mapping_ = parent_mapping
                        else
                                parent_mapping_ = .false.
                        endif

                        if (present(split_on_charge)) then
                                do_split_on_charge = split_on_charge
                        else
                                do_split_on_charge = self%postproc_split_on_charge
                        endif


                        if (.not.(allocated(self%minimum_mapping))) then
                                allocate(sm_get_mapping(0))
                                return
                        endif

                        if (parent_mapping_ .and. allocated(self%parent_mapping)) then
                                sm_get_mapping = composeMapping(self%minimum_mapping,self%parent_mapping)

                                if (do_split_on_charge .and. allocated(self%parent_charges)) then
                                        charge_mapping = signDiscretize(self%parent_charges)
                                        sm_get_mapping = fragmentMapping(sm_get_mapping,charge_mapping)
                                endif
                        else
                                sm_get_mapping = self%minimum_mapping

                                if (do_split_on_charge .and. allocated(self%domain_charges)) then
                                        charge_mapping = signDiscretize(self%domain_charges)
                                        sm_get_mapping = fragmentMapping(self%minimum_mapping,charge_mapping)
                                endif
                        endif

                endfunction sm_get_mapping

                pure function sm_is_mapping_charge_split(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: sm_is_mapping_charge_split

                        sm_is_mapping_charge_split = self%postproc_split_on_charge

                endfunction sm_is_mapping_charge_split

                pure function sm_has_parent_mapping(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: sm_has_parent_mapping

                        if (allocated(self%parent_mapping)) then
                                sm_has_parent_mapping = .true.
                        else
                                sm_has_parent_mapping = .false.
                        endif

                endfunction sm_has_parent_mapping

                pure function sm_get_projection(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable :: sm_get_projection(:,:)

                        if (allocated(self%minimum_projection)) then
                                sm_get_projection = self%minimum_projection
                        else
                                allocate(sm_get_projection(0,0))
                        endif

                endfunction sm_get_projection

                pure function sm_get_eigenvalues(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable :: sm_get_eigenvalues(:)

                        if (allocated(self%minimum_eigenvalues)) then
                                sm_get_eigenvalues = self%minimum_eigenvalues
                        else
                                allocate(sm_get_eigenvalues(0))
                        endif

                endfunction sm_get_eigenvalues

                pure function sm_get_siteCharges(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp),allocatable :: sm_get_siteCharges(:)

                        sm_get_siteCharges = self%minimum_cg_charges

                endfunction sm_get_siteCharges

                pure function sm_isTainted(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: sm_isTainted

                        sm_isTainted = self%tainted

                endfunction sm_isTainted

                pure subroutine sm_setTainted(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(inout)          :: self
                        logical,                    intent(in   )          :: toSet
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        self%tainted = toSet

                endsubroutine sm_setTainted

                pure function sm_getNAtoms(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si)    :: sm_getNAtoms

                        sm_getNAtoms = self%nAtoms

                endfunction sm_getNAtoms

                pure function sm_getNSites(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si)    :: sm_getNSites

                        sm_getNSites = self%nSites

                endfunction sm_getNSites

                function sm_recalc_residual(self,tdata) result(residuals)

                        use fit_common,         only: fit_residuals
                        use obj_trajectory,     only: traj_training_data
                        use obj_ll,             only: i_sp_dll
                        use core_kmeans,        only: sget_labels

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel),           intent(inout) :: self
                        class(traj_training_data),      intent(in   ) :: tdata 
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp), allocatable       :: residuals(:)
                        type(i_sp_dll), allocatable     :: backmapping(:)

                        allocate(residuals(size(self%residuals)))

                        call fit_residuals(tdata        = tdata,             &
                                           mapping      = self%minimum_mapping,   &
                               mapping_accumulator_type = self%getMapAccumType(),&
                                           residual_weights = self%residual_weights,&
                                           cgCharges    = self%minimum_cg_charges,&
                                           backmapping  = backmapping,            &
                                           residualList = residuals,              &
                                           residual_offsets = self%getNormOffsets(),&
                                           residual_scalings = self%getNormScalings())

                endfunction sm_recalc_residual

                subroutine sm_summarize(self,fileID,filename,verbose)

                        use core_convert, only: itoa
                        use core_filter,  only: composeMapping

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(spectralModel), intent(in   )          :: self
                        integer (si),               intent(in   ),optional :: fileID
                        character (len=*),          intent(in   ),optional :: filename
                        logical,                    intent(in   ),optional :: verbose
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        logical         :: verbose_
                        integer (si)    :: file_handle, status

                        !optional arg guard
                        if (present(verbose)) then
                                verbose_ = verbose
                        else
                                verbose_ = .false.
                        endif

                        if (present(filename) .and. present(fileID)) then
                                print*, "Only one of filename or fileID may be specified in&
                                                &model print statements."
                                return
                        endif

                        if (present(filename)) then
                                open(newunit = file_handle, file = filename, iostat = status)
                                if (status /= 0) then
                                        print*, "Model summary file couldn't be written."
                                        return
                                endif
                        else
                                file_handle = fileID
                        endif

                        if (verbose_) then
                                write(fileID,*) "Model Summary"
                                write(fileID,*) "   Number of Models:       ", self%getModelFreq()
                                write(fileID,*) "   Combined Residual:       ", self%getResidual()
                                write(fileID,*) "   Individual Residuals:       ", self%getFullResidual()
                                write(fileID,*) "   Residual normalization offset/scaling:       ", &
                                                                self%getNormOffsets(), '/', self%getNormScalings()
                                if (self%hasCI()) then
                                        write(fileID,*) "   CI (90): [", self%getCI(5.0_dp), &
                                                                    ",", self%getCI(95.0_dp), &
                                                "] (", self%bootstrap_dist_current_mark, "samples )"
                                        if (self%hasCI(bca=.true.)) then
                                                write(fileID,*) "   CI (90): [", self%getCI(5.0_dp,bca=.true.), &
                                                                            ",", self%getCI(95.0_dp,bca=.true.), &
                                                        "] (", self%bootstrap_dist_current_mark, "samples )"
                                                write(fileID,*) "   Boostrap dist:", &
                                                                self%bootstrap_distribution(1:self%bootstrap_dist_current_mark)
                                        endif
                                        write(fileID,*) "   Bias est.:", self%getBias()
                                endif
                                write(fileID,*) "   Residual weights:   ", self%residual_weights

                                if (allocated(self%parent_mapping)) then
                                        write(fileID,*) "   Mapping:   ",       self%getMapping(parent_mapping=.true.)
                                        write(fileID,*) "   Proxy mapping:   ", self%minimum_mapping
                                        write(fileID,*) "   Proxy optimized site charges:   ", self%minimum_cg_charges
                                        write(fileID,*) "   Proxy mapped site charges:   ",    self%mapped_cg_charges
                                else
                                        write(fileID,*) "   Mapping:   ", self%getMapping()

                                        write(fileID,*) "   Optimized site charges:   ", &
                                                        self%minimum_cg_charges
                                        write(fileID,*) "   Mapped site charges:   ", &
                                                        self%mapped_cg_charges
                                endif

                                write(fileID,*) "   Graph eig. proj.:   ", self%minimum_projection
                                write(fileID,*) "   Graph eigval.:   ", self%minimum_eigenvalues
                                write(fileID,*) "   Optimized Site Charges:   ", &
                                                self%minimum_cg_charges
                                write(fileID,*) "   Mapped Site Charges:   ", &
                                                self%mapped_cg_charges

                                write(fileID,*) "End Model Summary"
                        else
                                write(fileID,*) "Model Summary"
                                write(fileID,*) "   Number of Models:       ", self%getModelFreq()
                                write(fileID,*) "   Residual:       ", self%getResidual()
                                if (self%hasCI()) then
                                        write(fileID,*) "   CI (90): [", self%getCI(5.0_dp), &
                                                                    ",", self%getCI(95.0_dp), &
                                                "] (", self%bootstrap_dist_current_mark, "samples )"
                                        if (self%hasCI(bca=.true.)) then
                                                write(fileID,*) "   CI (90): [", self%getCI(5.0_dp,bca=.true.), &
                                                                            ",", self%getCI(95.0_dp,bca=.true.), &
                                                        "] (", self%bootstrap_dist_current_mark, "samples )"
                                        endif
                                        write(fileID,*) "   Bias est.:", self%getBias()
                                endif
                                write(fileID,*) "   Residual weights:   ", self%residual_weights
                                write(fileID,*) "   Mapping:   ", self%getMapping()
                                write(fileID,*) "   Optimized Site Charges:   ", &
                                                self%minimum_cg_charges
                                write(fileID,*) "   Mapped Site Charges:   ", &
                                                self%mapped_cg_charges

                                write(fileID,*) "End Model Summary"
                        endif

                        if (present(filename)) then
                                close(fileID)
                        endif

                endsubroutine sm_summarize

                pure subroutine sm_clone(lhs,rhs)

                        implicit none

                        class(spectralModel),intent(inout) :: lhs
                        class(model),              intent(in   ) :: rhs

                        select type (REFrhs => rhs)
                        class is (spectralModel)

                        if (REFrhs%isTainted()) then
                                call lhs%setTaint(.true.)
                                return
                        else
                                call lhs%setTaint(.false.)
                        endif

                        lhs%nAtoms = REFrhs%nAtoms
                        lhs%nSites = REFrhs%nSites

                        lhs%frequency = REFrhs%frequency

                        if (allocated(REFrhs%minimum_mapping)) then
                                lhs%minimum_mapping = REFrhs%minimum_mapping
                        endif

                        if (allocated(REFrhs%minimum_eigenvalues)) then
                                lhs%minimum_eigenvalues = REFrhs%minimum_eigenvalues
                        endif

                        if (allocated(REFrhs%minimum_projection)) then
                                lhs%minimum_projection = REFrhs%minimum_projection
                        endif

                        if (allocated(REFrhs%minimum_cg_charges)) then
                                lhs%minimum_cg_charges        = REFrhs%minimum_cg_charges
                        endif

                        if (allocated(REFrhs%mapped_cg_charges)) then
                                lhs%mapped_cg_charges        = REFrhs%mapped_cg_charges
                        endif

                        lhs%residuals                 = REFrhs%residuals

                        lhs%residual_weights          = REFrhs%residual_weights

                        call REFrhs%copyNorms(lhs)

                        return
                        end select

                        !else
                        call lhs%setTaint(.true.)
                endsubroutine sm_clone

                pure subroutine sm_init_ci(self, dist_size, jk)

                        implicit none

                        class(spectralModel),          intent(inout) :: self
                        integer (si),                        intent(in   ) :: dist_size
                        logical,                   optional, intent(in   ) :: jk

                        logical         :: jk_

                        !optional arg guard
                        if (present(jk)) then
                                jk_ = jk
                        else
                                jk_ = .false.
                        endif

                        if (jk_) then
                                self%jackknife_distribution_isSorted = .false.

                                !careful allocate
                                if (allocated(self%jackknife_distribution)) then
                                        if (.not. (size(self%jackknife_distribution,1) .eq. dist_size)) then
                                                deallocate(self%jackknife_distribution)
                                                allocate(self%jackknife_distribution(dist_size))
                                        endif
                                else
                                        allocate(self%jackknife_distribution(dist_size))
                                endif

                                self%jackknife_dist_current_mark = 0
                        else
                                self%bootstrap_distribution_isSorted = .false.

                                !careful allocate
                                if (allocated(self%bootstrap_distribution)) then
                                        if (.not. (size(self%bootstrap_distribution,1) .eq. dist_size)) then
                                                deallocate(self%bootstrap_distribution)
                                                allocate(self%bootstrap_distribution(dist_size))
                                        endif
                                else
                                        allocate(self%bootstrap_distribution(dist_size))
                                endif

                                self%bootstrap_dist_current_mark = 0
                        endif

                endsubroutine sm_init_ci

                subroutine sm_calc_ci(self,bca)

                        use core_sort,          only: qsort
                        use core_stat,          only: get_bca_bias, get_bca_accel

                        implicit none

                        class(spectralModel),           intent(inout) :: self
                        logical,                    optional, intent(in   ) :: bca

                        logical                                  :: bca_

                        !optional arg guard
                        if (present(bca)) then
                                bca_ = bca
                        else
                                bca_ = .false.
                        endif

                        if (allocated(self%bootstrap_distribution)) then
                                call qsort(self%bootstrap_distribution(&
                                                1:self%bootstrap_dist_current_mark))
                                self%bootstrap_distribution_isSorted = .true.
                        endif

                        if (bca_) then
                                call qsort(self%jackknife_distribution(&
                                                1:self%jackknife_dist_current_mark))
                                self%bca_accel = get_bca_accel(&
                                                 self%jackknife_distribution(1:self%jackknife_dist_current_mark))
                                self%bca_bias  = get_bca_bias(&
                                                 self%bootstrap_distribution(1:self%bootstrap_dist_current_mark),&
                                                 self%residuals(1))
                        endif

                endsubroutine sm_calc_ci

                pure function sm_has_ci(self,bca)

                        implicit none

                        class(spectralModel),intent(in   ) :: self
                        logical,optional,    intent(in   ) :: bca

                        logical :: sm_has_ci

                        logical :: bca_

                        if (present(bca)) then
                                bca_ = bca
                        else
                                bca_ = .false.
                        endif

                        !this isn't ideal-- we should have more direct indications
                        !that these are valid.

                        if (bca_) then
                                sm_has_ci = (allocated(self%jackknife_distribution) .and. &
                                             self%bootstrap_distribution_isSorted)
                        else
                                sm_has_ci = self%bootstrap_distribution_isSorted
                        endif

                endfunction sm_has_ci

                function sm_get_ci(self, percentile, bca) result(bound)

                        use core_stat,    only: bca_transform_percentile
                        use core_convert, only: itoa

                        implicit none

                        class(spectralModel),intent(in   ) :: self
                        real    (dp),              intent(in   ) :: percentile
                        logical,         optional, intent(in   ) :: bca

                        real    (dp)                             :: bound

                        logical                                  :: bca_

                        integer (si)                             :: p_index
                        real    (dp)                             :: corr_perc

                        !optional arg guard
                        if (present(bca)) then
                                bca_ = bca
                        else
                                bca_ = .false.
                        endif

                        if (bca_) then
                                corr_perc = bca_transform_percentile(percentile,self%bca_accel,self%bca_bias)
                        else
                                corr_perc = percentile
                        endif

                        p_index = int((corr_perc/100.0)*(self%bootstrap_dist_current_mark))

                        if (self%hasCI()) then
                                if (p_index < 1) then
                                        print*, "warning: bootstrap index:",p_index, " was out of &
                                                &the range of expression of the bootstrapped sample. &
                                                &Substituting with 1."
                                        p_index = 1
                                else if (p_index > self%bootstrap_dist_current_mark) then
                                        print*, "warning: bootstrap index:",p_index, " was out of &
                                                &the range of expression of the bootstrapped sample. &
                                                &Substituting with "//itoa(self%bootstrap_dist_current_mark)//"."
                                        p_index = self%bootstrap_dist_current_mark
                                endif
                                bound = self%bootstrap_distribution(p_index)
                                return
                        else
                                bound = huge(bound)
                                return
                        endif

                endfunction sm_get_ci

                pure function sm_get_bias(self) result(estimate)

                        implicit none

                        class(spectralModel),intent(in   ) :: self

                        real    (dp)                             :: estimate

                        if (self%hasCI()) then
                                estimate = sum(self%bootstrap_distribution(1:self%bootstrap_dist_current_mark))/&
                                                self%bootstrap_dist_current_mark
                                estimate = estimate - self%residuals(1)
                                return
                        else
                                estimate = 0
                                return
                        endif

                endfunction sm_get_bias

                pure subroutine sm_add_ci_values(self, toAdd, jk) 
                                                               
                        implicit none                          

                        class(spectralModel),         intent(inout) :: self
                        real    (dp),                       intent(in   ) :: toAdd(:)
                        logical,                   optional,intent(in   ) :: jk

                        logical                                  :: jk_

                        integer (si)                             :: iter,copyIter

                        real    (dp),allocatable                 :: scratch(:)

                        !optional arg guard
                        if (present(jk)) then
                                jk_ = jk
                        else
                                jk_ = .false.
                        endif

                        !simple growing array approach.
                        if (jk_) then
                                do iter=1,size(toAdd)

                                        self%jackknife_dist_current_mark = self%jackknife_dist_current_mark +1 

                                        !if our write location is beyond the size of the array, expand the array.
                                        if (self%jackknife_dist_current_mark > size(self%jackknife_distribution)) then

                                                call move_alloc(self%jackknife_distribution, scratch)
                                                allocate(self%jackknife_distribution(&
                                                                2*(self%jackknife_dist_current_mark-1)))

                                                do copyIter=1,size(scratch)
                                                        self%jackknife_distribution(copyIter) = &
                                                                scratch(copyIter)
                                                enddo

                                        endif

                                        self%jackknife_distribution(self%jackknife_dist_current_mark) = toAdd(iter)
                                enddo

                                self%jackknife_distribution_isSorted = .false.
                        else
                                do iter=1,size(toAdd)

                                        self%bootstrap_dist_current_mark = self%bootstrap_dist_current_mark +1 

                                        !if our write location is beyond the size of the array, expand the array.
                                        if (self%bootstrap_dist_current_mark > size(self%bootstrap_distribution)) then

                                                call move_alloc(self%bootstrap_distribution, scratch)
                                                allocate(self%bootstrap_distribution(&
                                                                2*(self%bootstrap_dist_current_mark-1)))

                                                do copyIter=1,size(scratch)
                                                        self%bootstrap_distribution(copyIter) = &
                                                                scratch(copyIter)
                                                enddo

                                        endif

                                        self%bootstrap_distribution(self%bootstrap_dist_current_mark) = toAdd(iter)
                                enddo

                                self%bootstrap_distribution_isSorted = .false.
                        endif

                endsubroutine sm_add_ci_values

                elemental function sm_get_model_freq(self) result(count)

                        implicit none

                        class(spectralModel),intent(in   ) :: self

                        integer (si)            :: count

                        count = self%frequency

                endfunction sm_get_model_freq

                pure subroutine sm_incr_model_freq(self,amount)

                        implicit none

                        class(spectralModel),intent(inout)          :: self
                        integer (si)              ,intent(in   ),optional :: amount

                        if (present(amount)) then
                                self%frequency = self%frequency + amount
                        else 
                                self%frequency = self%frequency + 1
                        endif

                endsubroutine sm_incr_model_freq

                pure subroutine sm_decr_model_freq(self,amount)

                        implicit none

                        class(spectralModel),intent(inout)          :: self
                        integer (si)              ,intent(in   ),optional :: amount

                        if (present(amount)) then
                                self%frequency = self%frequency - amount
                        else 
                                self%frequency = self%frequency - 1
                        endif

                endsubroutine sm_decr_model_freq

                !For spectral models, it's essentially impossible to restrict our domain effectively without 
                !solving the entire model. We simply randomly sample models.
                function sm_get_random_mapping(self,tdata)

                        use core_random,    only: gen_rand_int
                        use obj_trajectory, only: traj_training_data

                        implicit none

                        class(spectralModel),           intent(in   ) :: self
                        class(traj_training_data),      intent(in   ) :: tdata

                        integer (si),allocatable                 :: sm_get_random_mapping(:)
                        integer (si)                             :: iter

                        if ((self%nSites == 0) .or. (tdata%trj%nAtoms == 0)) then
                                allocate(sm_get_random_mapping(0))
                        else
                                allocate(sm_get_random_mapping(tdata%trj%nAtoms))
                                do iter=1,size(sm_get_random_mapping)
                                        sm_get_random_mapping(iter) = gen_rand_int(lower=1,&
                                                                                   upper=self%nSites)
                                enddo
                        endif

                endfunction sm_get_random_mapping

                pure function sm_get_map_accum_type(self)

                        implicit none

                        class(spectralModel),   intent(in   ) :: self

                        character(:),allocatable              :: sm_get_map_accum_type

                        if (allocated(self%map_accumulator_type)) then
                                sm_get_map_accum_type = self%map_accumulator_type
                        else
                                allocate(character(len=0) :: sm_get_map_accum_type)
                        endif

                endfunction sm_get_map_accum_type

                pure function get_nSites(self)

                        implicit none

                        class(spectralModel),intent(in   ) :: self

                        integer (si)    :: get_nSites

                        get_nSites = self%nSites
                        
                endfunction get_nSites


                pure function sm_eq(lhs,rhs)

                        use core_stat,          only: assignment_difference

                        implicit none

                        class(spectralModel),intent(in   ) :: lhs
                        class(model),              intent(in   ) :: rhs

                        logical         :: sm_eq

                        !logical         :: lhs_charges_allocated, rhs_charges_allocated
                        logical         :: lhs_mapping_allocated, rhs_mapping_allocated
                        !logical         :: lhs_projection_allocated, rhs_projection_allocated
                        !logical         :: lhs_eigenvalues_allocated, rhs_eigenvalues_allocated

                        !this is lowered due to variation in the charge residual value based on labeling.
                        !numerical instability?
                        !real    (dp),parameter :: tolerance = real(10.0_dp**(-4),dp)

                        sm_eq = .true.


                        select type(rhs)
                        class is (spectralModel)


                        !if (any(abs(lhs%residuals - rhs%residuals) > tolerance)) then
                        !        sm_eq = .false.
                        !        return
                        !endif

                        !if (lhs%nAtoms /= rhs%nAtoms) then
                        !        sm_eq = .false.
                        !        return
                        !endif

                        !if (lhs%nSites /= rhs%nSites) then
                        !        sm_eq = .false.
                        !        return
                        !endif

                        !lhs_charges_allocated = allocated(lhs%minimum_cg_charges)
                        !rhs_charges_allocated = allocated(rhs%minimum_cg_charges)

                        !if (lhs_charges_allocated .and. rhs_charges_allocated) then
                        !        !if (any(abs(rhs%minimum_cg_charges - lhs%minimum_cg_charges ) > tolerance )) then
                        !        !        sm_eq = .false.
                        !        !        return
                        !        !endif
                        !else if (lhs_charges_allocated) then
                        !        sm_eq = .false.
                        !        return
                        !else if (rhs_charges_allocated) then
                        !        sm_eq = .false.
                        !        return
                        !endif

                        lhs_mapping_allocated = allocated(lhs%minimum_mapping)
                        rhs_mapping_allocated = allocated(rhs%minimum_mapping)

                        if (lhs_mapping_allocated .and. rhs_mapping_allocated) then
                                if (assignment_difference(rhs%minimum_mapping,lhs%minimum_mapping) > (10.0_dp**(-10))) then
                                        sm_eq = .false.
                                        return
                                endif
                        else if (lhs_mapping_allocated) then
                                sm_eq = .false.
                                return
                        else if (rhs_mapping_allocated) then
                                sm_eq = .false.
                                return
                        endif

                        !lhs_eigenvalues_allocated = allocated(lhs%minimum_eigenvalues)
                        !rhs_eigenvalues_allocated = allocated(rhs%minimum_eigenvalues)

                        !if (lhs_eigenvalues_allocated .and. rhs_eigenvalues_allocated) then
                        !        if (.not. all(shape(rhs%minimum_eigenvalues) == shape(lhs%minimum_eigenvalues))) then
                        !                sm_eq = .false.
                        !                return
                        !        endif
                        !        if (any(abs(rhs%minimum_eigenvalues - lhs%minimum_eigenvalues ) > tolerance )) then
                        !                sm_eq = .false.
                        !                return
                        !        endif
                        !else if (lhs_eigenvalues_allocated) then
                        !        sm_eq = .false.
                        !        return
                        !else if (rhs_eigenvalues_allocated) then
                        !        sm_eq = .false.
                        !        return
                        !endif

                        !lhs_projection_allocated = allocated(lhs%minimum_projection)
                        !rhs_projection_allocated = allocated(rhs%minimum_projection)

                        !if (lhs_projection_allocated .and. rhs_projection_allocated) then
                        !        if (.not. all(shape(rhs%minimum_projection) == shape(lhs%minimum_projection))) then
                        !                sm_eq = .false.
                        !                return
                        !        endif
                        !        if (any(abs(rhs%minimum_projection - lhs%minimum_projection ) > tolerance )) then
                        !                sm_eq = .false.
                        !                return
                        !        endif
                        !else if (lhs_projection_allocated) then
                        !        sm_eq = .false.
                        !        return
                        !else if (rhs_projection_allocated) then
                        !        sm_eq = .false.
                        !        return
                        !endif
                        end select

                        return

                endfunction sm_eq
endmodule obj_spectralModel
