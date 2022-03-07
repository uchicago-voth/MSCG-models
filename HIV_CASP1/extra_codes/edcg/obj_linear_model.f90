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
! Provides routines to generate mappings contigious in primary sequence.

module obj_linearDivisionModel

        use env_kindtypes, only: si, sp, dp
        use abs_obj_model, only: model, model_config
        use fit_common,    only: num_residuals
        use obj_Rloggers,  only: t_rLogger
        use obj_trajectory,only: traj


        implicit none

        private

        public          :: linearDivisionModel

        type, extends(model) :: linearDivisionModel
                private

                !fund. state variables
                logical                  :: tainted = .false.
                logical                  :: record  = .false.

                integer (si)             :: frequency = 1

                !mapping state variables
                integer (si)             :: nSites = 0
                integer (si)             :: nAtoms   = 0   !input number of sites
                integer (si),allocatable :: minimum_boundary_residues(:)
                integer (si),allocatable :: minimum_mapping(:)

                !when dealing with a mapped trajectory, this gives a map from
                !the parent trajectory to the mapped trajectory, so that 
                !we can get the full resolution mapping in a simple manner.
                integer (si),allocatable :: parent_mapping(:)
                real    (dp),allocatable :: parent_charges(:)

                real    (dp),allocatable :: domain_charges(:)
                real    (dp),allocatable :: minimum_cg_charges(:)
                real    (dp),allocatable :: mapped_cg_charges(:)

                logical                  :: postproc_split_on_charge

                !residual vars
                real    (dp)             :: residuals(num_residuals+1)      = -1
                real    (dp)             :: residual_weights(num_residuals) = -1

                !bootstrapping variables
                real    (dp),allocatable :: bootstrap_distribution(:)
                integer (si)             :: bootstrap_dist_current_mark = 0
                logical                  :: bootstrap_distribution_issorted =.false.

                real    (dp)             :: bca_bias,bca_accel

                real    (dp),allocatable :: jackknife_distribution(:)
                integer (si)             :: jackknife_dist_current_mark = 0
                logical                  :: jackknife_distribution_issorted =.false.

                !! simulated annelaing vars
                type    (t_rLogger)      :: sa_logger
                integer (si)             :: sa_logger_stride = 1
                real    (dp)             :: sa_accept,&    !initial temp
                                            sa_descent_control!control for SA process

                !! lucky descent vars
                type    (t_rLogger)      :: ld_logger
                integer (si)             :: ld_logger_stride = 1
                integer (si)             :: ld_max_steps     !max number of stepds

                !holds type of accumulation used in mapping (e.g. center of mass).
                character(:),allocatable :: map_accumulator_type

                contains
                        procedure,public :: reset         => ldm_reset
                        procedure,public :: anneal        => ldm_sim_annealing
                        procedure,public :: refine        => ldm_lucky_descent
                        procedure,public :: summarize     => ldm_summarize

                        procedure,public :: recalcResidual=> ldm_recalc_residual

                        procedure,public :: getAnnealLog  => ldm_get_anneal_log
                        procedure,public :: getRefineLog  => ldm_get_refine_log

                        procedure,public :: getMapping    => ldm_get_mapping
                        procedure,public :: getBoundaryResidues &
                                                          => ldm_get_boundaryResidues

                        procedure,public :: isMapChargeSplit => ldm_is_mapping_charge_split
                        procedure,public :: hasParentMap  => ldm_has_parent_mapping

                        procedure        :: setAnnealLog  => ldm_set_anneal_log
                        procedure        :: setRefineLog  => ldm_set_refine_log
                        procedure,public :: getSiteCharges =>ldm_get_siteCharges
                        procedure,public :: getResidual    =>ldm_get_residual
                        procedure,public :: getFullResidual=>ldm_get_full_residual
                        procedure,public :: initialize    => ldm_initialize
                        procedure,public :: configure     => ldm_configure

                        procedure,public :: isTainted     => ldm_isTainted

                        procedure        :: setTaint      => ldm_setTainted

                        procedure,public :: getModelFreq  => ldm_get_model_freq
                        procedure,public :: incrModelFreq => ldm_incr_model_freq
                        procedure,public :: decrModelFreq => ldm_decr_model_freq

                        procedure,public :: initializeCI  => ldm_init_ci

                        procedure,public :: calcCI        => ldm_calc_ci
                        procedure,public :: hasCI         => ldm_has_ci

                        procedure,public :: getCI         => ldm_get_ci
                        !procedure,public :: getLowerCI    => ldm_get_lower_ci

                        procedure,public :: getBias       => ldm_get_bias

                        procedure,public :: addCIValues   => ldm_add_ci_values

                        procedure,public :: getRandomMapping => ldm_get_random_mapping

                        procedure,public :: getMapAccumType => ldm_get_map_accum_type

                        procedure        :: getNSites     => get_nSites

                        procedure        :: clone         => ldm_clone
                        procedure        :: eq            => ldm_eq
        end type linearDivisionModel

        contains
                pure subroutine ldm_configure(self,config) 

                        implicit none

                        class(linearDivisionModel),intent(inout) :: self
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

                        self%ld_max_steps       = config%refine_steps
                        self%sa_accept          = config%anneal_accept
                        self%sa_descent_control = config%anneal_param

                        self%record = config%record

                        if (self%record) then
                                call self%sa_logger%configure(config%size_buffer_anneal)
                                call self%ld_logger%configure(config%size_buffer_refine)
                        endif

                        self%nSites = config%nSites

                        self%postproc_split_on_charge = config%mapping_charge_split

                endsubroutine ldm_configure

                pure subroutine ldm_initialize(self)

                        implicit none

                        class(linearDivisionModel),intent(inout) :: self

                        allocate(self%minimum_boundary_residues(self%nSites-1))
                        allocate(self%minimum_cg_charges(self%nSites))

                        if (self%record) then
                                call self%sa_logger%initialize()
                                call self%ld_logger%initialize()
                        endif

                endsubroutine ldm_initialize

                pure subroutine ldm_reset(self)

                        implicit none

                        class(linearDivisionModel), intent(inout) :: self

                        self%residuals = -1

                        call self%sa_logger%reset(shallow=.true.)
                        call self%ld_logger%reset(shallow=.true.)

                endsubroutine ldm_reset

                subroutine ldm_sim_annealing(self,tdata,acceptance,descentControl)

                        use obj_trajectory,     only: traj_training_data
                        use routines_boundaryStepwiseOptimization,&
                                                only: boundaryStepwiseAnnealing
                        use core_filter,        only: boundariesToMapping, valueCollapse
                        use fit_common,         only: fit_residuals
                        use obj_ll,             only: i_sp_dll

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(inout) :: self
                        class(traj_training_data),  intent(in   ) :: tdata
                        real (dp),     optional,    intent(in   ) :: acceptance,descentControl
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !temporary shim variable.
                        type(i_sp_dll),allocatable :: backmapping(:)

                        if (self%nAtoms == 0) then
                                self%nAtoms = tdata%trj%nAtoms
                        endif

                        if (present(acceptance)) then
                                self%sa_accept = acceptance
                        endif

                        if (present(descentControl)) then
                                self%sa_descent_control = descentControl
                        endif

                        self%domain_charges = tdata%trj%atomCharges

                        if (tdata%trj%isMapped()) then
                                self%parent_mapping = tdata%trj%getParentMap()
                                self%parent_charges = tdata%trj%getParentCharges()
                        endif

                        call boundaryStepwiseAnnealing( boundaries = self%minimum_boundary_residues,&
                                                       siteCharges = self%minimum_cg_charges, &
                                                            logger = self%sa_logger,&
                                                      residualList = self%residuals, &
                                                  residual_offsets = self%getNormOffsets(), &
                                                 residual_scalings = self%getNormScalings(), &
                                                             tdata = tdata, &
                                          mapping_accumulator_type = self%getMapAccumType(),&
                                                  residual_weights = self%residual_weights,&
                                                    startingAccept = self%sa_accept, &
                                                          tempRate = self%sa_descent_control)

                        self%minimum_mapping   = boundariesToMapping(self%minimum_boundary_residues,self%nAtoms)

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
                                           residualList     = self%residuals,    &
                                           residual_offsets = self%getNormOffsets(),&
                                           residual_scalings= self%getNormScalings())

                endsubroutine ldm_sim_annealing

                subroutine ldm_lucky_descent(self,tdata,maxSteps)

                        use obj_trajectory,     only: traj_training_data
                        use routines_boundaryStepwiseOptimization,&
                                                only: boundaryStepwiseDescent
                        use core_filter,        only: boundariesToMapping, valueCollapse
                        use fit_common,         only: fit_residuals
                        use obj_ll,             only: i_sp_dll

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(inout) :: self
                        class(traj_training_data),  intent(in   ) :: tdata
                        integer (si),  optional,    intent(in   ) :: maxSteps
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        type(i_sp_dll),allocatable :: backmapping(:)

                        if (self%nAtoms == 0) then
                                self%nAtoms = tdata%trj%nAtoms
                        endif

                        if (present(maxSteps)) then
                                self%ld_max_Steps = maxSteps
                        endif

                        self%domain_charges = tdata%trj%atomCharges

                        if (tdata%trj%isMapped()) then
                                self%parent_mapping = tdata%trj%getParentMap()
                                self%parent_charges = tdata%trj%getParentCharges()
                        endif

                        call boundaryStepwiseDescent( boundaries = self%minimum_boundary_residues, &
                                                     siteCharges = self%minimum_cg_charges, &
                                                          logger = self%ld_logger,&
                                                    residualList = self%residuals,&
                                                residual_offsets = self%getNormOffsets(),&
                                               residual_scalings = self%getNormScalings(),&
                                                           tdata = tdata, &
                                        mapping_accumulator_type = self%getMapAccumType(),&
                                                residual_weights = self%residual_weights, &
                                               max_SD_iterations = self%ld_max_Steps)

                        self%minimum_mapping = boundariesToMapping(self%minimum_boundary_residues,self%nAtoms)

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
                                           residualList     = self%residuals,    &
                                           residual_offsets = self%getNormOffsets(),&
                                           residual_scalings= self%getNormScalings())

                endsubroutine ldm_lucky_descent

                pure subroutine ldm_set_anneal_log(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(inout)          :: self
                        real    (sp),               intent(in   )          :: toSet(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        call self%sa_logger%setLog(toSet)

                endsubroutine ldm_set_anneal_log

                pure function ldm_get_anneal_log(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (sp),allocatable  :: ldm_get_anneal_log(:)

                        if (self%sa_logger%isActive()) then
                            call self%sa_logger%sgetLog(ldm_get_anneal_log)
                        else
                            allocate(ldm_get_anneal_log(0))
                        endif

                endfunction ldm_get_anneal_log

                pure subroutine ldm_set_refine_log(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(inout)          :: self
                        real    (sp),               intent(in   )          :: toSet(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        call self%ld_logger%setLog(toSet)

                endsubroutine ldm_set_refine_log

                pure function ldm_get_refine_log(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (sp),allocatable  :: ldm_get_refine_log(:)

                        if (self%ld_logger%isActive()) then
                            call self%ld_logger%sgetLog(ldm_get_refine_log)
                        else
                            allocate(ldm_get_refine_log(0))
                        endif

                endfunction ldm_get_refine_log

                elemental function ldm_get_residual(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp)                 :: ldm_get_residual

                        ldm_get_residual = self%residuals(1)

                endfunction ldm_get_residual

                pure function ldm_get_full_residual(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp),allocatable  :: ldm_get_full_residual(:)

                        ldm_get_full_residual = self%residuals

                endfunction ldm_get_full_residual

                pure function ldm_get_mapping(self,parent_mapping,split_on_charge)

                        use core_filter, only:boundariesToMapping, composeMapping, &
                                              fragmentMapping, signDiscretize

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel),          intent(in   )          :: self
                        logical,                   optional, intent(in   )          :: parent_mapping
                        logical,                   optional, intent(in   )          :: split_on_charge
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si),allocatable :: ldm_get_mapping(:)
                        logical                  :: parent_mapping_
                        logical                  :: do_split_on_charge
                        integer (si),allocatable :: working_mapping(:)

                        integer (si),allocatable :: charge_mapping(:)

                        if (present(parent_mapping)) then
                                parent_mapping_ = parent_mapping
                        else
                                parent_mapping_ = .false.
                        endif

                        if (.not.(allocated(self%minimum_boundary_residues))) then
                                allocate(ldm_get_mapping(0))
                                return
                        endif

                        if (present(split_on_charge)) then
                                do_split_on_charge = split_on_charge
                        else
                                do_split_on_charge = self%postproc_split_on_charge
                        endif

                        if (parent_mapping_ .and. allocated(self%parent_mapping)) then
                                working_mapping = boundariesToMapping(self%minimum_boundary_residues,&
                                                                          self%nAtoms)
                                ldm_get_mapping = composeMapping(working_mapping,self%parent_mapping)

                                if (do_split_on_charge .and. allocated(self%parent_charges)) then
                                        charge_mapping = signDiscretize(self%parent_charges)
                                        ldm_get_mapping = fragmentMapping(ldm_get_mapping,charge_mapping)
                                endif
                        else
                                ldm_get_mapping = boundariesToMapping(self%minimum_boundary_residues,&
                                                                  self%nAtoms)

                                if (do_split_on_charge .and. allocated(self%parent_charges)) then
                                        charge_mapping  = signDiscretize(self%domain_charges)
                                        ldm_get_mapping = fragmentMapping(ldm_get_mapping,charge_mapping)
                                endif
                        endif


                endfunction ldm_get_mapping

                pure function ldm_get_boundaryResidues(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si),allocatable :: ldm_get_boundaryResidues(:)

                        ldm_get_boundaryResidues = self%minimum_boundary_residues

                endfunction ldm_get_boundaryResidues

                pure function ldm_is_mapping_charge_split(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: ldm_is_mapping_charge_split

                        ldm_is_mapping_charge_split = self%postproc_split_on_charge

                endfunction ldm_is_mapping_charge_split

                pure function ldm_has_parent_mapping(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: ldm_has_parent_mapping

                        if (allocated(self%parent_mapping)) then
                                ldm_has_parent_mapping = .true.
                        else
                                ldm_has_parent_mapping = .false.
                        endif

                endfunction ldm_has_parent_mapping

                pure function ldm_get_siteCharges(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp),allocatable :: ldm_get_siteCharges(:)

                        ldm_get_siteCharges = self%minimum_cg_charges

                endfunction ldm_get_siteCharges

                pure function ldm_isTainted(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: ldm_isTainted

                        ldm_isTainted = self%tainted

                endfunction ldm_isTainted

                pure subroutine ldm_setTainted(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(inout)          :: self
                        logical,                    intent(in   )          :: toSet
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        self%tainted = toSet

                endsubroutine ldm_setTainted

                pure function ldm_getNAtoms(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si)    :: ldm_getNAtoms

                        ldm_getNAtoms = self%nAtoms

                endfunction ldm_getNAtoms

                pure function ldm_getNSites(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si)    :: ldm_getNSites

                        ldm_getNSites = self%nSites

                endfunction ldm_getNSites

                function ldm_recalc_residual(self,tdata) result(residuals)

                        use fit_common,         only: fit_residuals
                        use obj_trajectory,     only: traj_training_data
                        use obj_ll,             only: i_sp_dll
                        use core_filter,        only: boundariesToMapping

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(inout)       :: self
                        class(traj_training_data),  intent(in   )       :: tdata
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp), allocatable       :: residuals(:)

                        type(i_sp_dll), allocatable     :: backmapping(:)

                        allocate(residuals(size(self%residuals)))

                        self%minimum_mapping = boundariesToMapping(self%minimum_boundary_residues, &
                                                                   self%nAtoms)

                        call fit_residuals(tdata        = tdata,             &
                                           mapping      = self%minimum_mapping,   &
                               mapping_accumulator_type = self%getMapAccumType(),&
                                       residual_weights = self%residual_weights,  &
                                           cgCharges    = self%minimum_cg_charges,&
                                           backmapping  = backmapping,            &
                                           residualList = residuals,&
                                           residual_offsets= self%getNormOffsets(),&
                                           residual_scalings= self%getNormScalings())

                endfunction ldm_recalc_residual

                subroutine ldm_summarize(self,fileID,filename,verbose)

                        use core_convert, only: itoa
                        use core_filter,  only: composeMapping

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(linearDivisionModel), intent(in   )          :: self
                        integer (si),               intent(in   ),optional :: fileID
                        character (len=*),          intent(in   ),optional :: filename
                        logical,                    intent(in   ),optional :: verbose
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        logical         :: verbose_
                        integer (si)    :: file_handle, status

                        !default arg guard
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
                                        write(fileID,*) "   Mapping:   ", self%getMapping(parent_mapping=.true.)
                                        write(fileID,*) "   Proxy mapping:   ", self%minimum_mapping
                                        write(fileID,*) "   Proxy mapping boundaries:   ", &
                                                        self%minimum_boundary_residues
                                        write(fileID,*) "   Proxy optimized site charges:   ", &
                                                        self%minimum_cg_charges
                                        write(fileID,*) "   Proxy mapped site charges:   ", &
                                                        self%mapped_cg_charges
                                else
                                        write(fileID,*) "   Mapping:   ", self%getMapping()
                                        write(fileID,*) "   Mapping boundaries:   ", &
                                                        self%minimum_boundary_residues
                                        write(fileID,*) "   Optimized site charges:   ", &
                                                        self%minimum_cg_charges
                                        write(fileID,*) "   Mapped site charges:   ", &
                                                        self%mapped_cg_charges
                                endif

                                write(fileID,*) "   Simulated Annealing paramters:"
                                write(fileID,*) "      Initial Acceptance: ", self%sa_accept
                                write(fileID,*) "      Descent Parameter: ", self%sa_descent_control
                                write(fileID,*) "   Simulated Annealing iterations:  ",&
                                                       self%sa_logger%getLength()*self%sa_logger_stride
                                if (self%sa_logger%isActive()) then
                                        write(fileID,*) "   Simulated Annealing log (stride of "&
                                                            //itoa(self%sa_logger_stride)//") :   ",&
                                                               self%sa_logger%getLog_sp()
                                else
                                        write(fileID,*) "           (not logged)"
                                endif

                                write(fileID,*) "   Descent iterations:   ",&
                                                       self%ld_logger%getLength()*self%ld_logger_stride
                                write(fileID,*) "      Max iterations:", self%ld_max_steps

                                if (self%ld_logger%isActive()) then
                                        write(fileID,*) "   Descent log (stride of "&
                                                    //itoa(self%ld_logger_stride)//") :   ",&
                                                       self%ld_logger%getLog_sp()
                                else
                                        write(fileID,*) "           (not logged)"
                                endif
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
                                write(fileID,*) "   Boundary Residues:   ", &
                                                self%minimum_boundary_residues
                                write(fileID,*) "   Optimized Site Charges:   ", &
                                                self%minimum_cg_charges
                                write(fileID,*) "   Mapped Site Charges:   ", &
                                                self%mapped_cg_charges

                                write(fileID,*) "   Simulated Annealing paramters:"
                                write(fileID,*) "      Initial Acceptance: ", self%sa_accept
                                write(fileID,*) "      Descent Parameter: ", self%sa_descent_control
                                write(fileID,*) "   Simulated Annealing iterations:  ",&
                                                       self%sa_logger%getLength()*self%sa_logger_stride

                                write(fileID,*) "   Descent iterations:   ",&
                                                       self%ld_logger%getLength()*self%ld_logger_stride
                                write(fileID,*) "      Max iterations:", self%ld_max_steps
                                write(fileID,*) "End Model Summary"
                        endif

                        if (present(filename)) then
                                close(fileID)
                        endif

                endsubroutine ldm_summarize

                pure subroutine ldm_clone(lhs,rhs)

                        implicit none

                        class(linearDivisionModel),intent(inout) :: lhs
                        class(model),              intent(in   ) :: rhs

                        select type (REFrhs => rhs)
                        class is (linearDivisionModel)

                        if (REFrhs%isTainted()) then
                                call lhs%setTaint(.true.)
                                return
                        else
                                call lhs%setTaint(.false.)
                        endif

                        lhs%nAtoms = REFrhs%nAtoms
                        lhs%nSites = REFrhs%nSites

                        call lhs%setAnnealLog(REFrhs%getAnnealLog())
                        call lhs%setRefineLog(REFrhs%getRefineLog())

                        lhs%frequency = REFrhs%frequency

                        if (allocated(REFrhs%minimum_boundary_residues)) then
                                lhs%minimum_boundary_residues = REFrhs%minimum_boundary_residues
                        endif

                        if (allocated(REFrhs%minimum_cg_charges)) then
                                lhs%minimum_cg_charges        = REFrhs%minimum_cg_charges
                        endif

                        if (allocated(REFrhs%minimum_mapping)) then
                                lhs%minimum_mapping        = REFrhs%minimum_mapping
                        endif

                        if (allocated(REFrhs%mapped_cg_charges)) then
                                lhs%mapped_cg_charges        = REFrhs%mapped_cg_charges
                        endif

                        lhs%sa_accept = REFrhs%sa_accept
                        lhs%sa_descent_control = REFrhs%sa_descent_control
                        lhs%ld_max_steps = REFrhs%ld_max_steps

                        lhs%residuals                 = REFrhs%residuals

                        lhs%residual_weights          = REFrhs%residual_weights

                        call REFrhs%copyNorms(lhs)

                        return
                        end select

                        !else
                        call lhs%setTaint(.true.)
                endsubroutine ldm_clone

                pure subroutine ldm_init_ci(self, dist_size, jk)

                        implicit none

                        class(linearDivisionModel),          intent(inout) :: self
                        integer (si),                        intent(in   ) :: dist_size
                        logical,                   optional, intent(in   ) :: jk

                        logical         :: jk_

                        !default arg guard
                        if (present(jk)) then
                                jk_ = jk
                        else
                                jk_ = .false.
                        endif

                        if (jk_) then
                                self%jackknife_distribution_isSorted = .false.

                                !careful allocation
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

                                !careful allocation
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

                endsubroutine ldm_init_ci

                subroutine ldm_calc_ci(self,bca)

                        use core_sort,          only: qsort
                        use core_stat,          only: get_bca_bias, get_bca_accel

                        implicit none

                        class(linearDivisionModel),           intent(inout) :: self
                        logical,                    optional, intent(in   ) :: bca

                        logical                                  :: bca_

                        !default arg guards
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

                endsubroutine ldm_calc_ci

                pure function ldm_has_ci(self,bca)

                        implicit none

                        class(linearDivisionModel),intent(in   ) :: self
                        logical,optional,    intent(in   ) :: bca

                        logical :: ldm_has_ci

                        logical :: bca_

                        if (present(bca)) then
                                bca_ = bca
                        else
                                bca_ = .false.
                        endif

                        !this isn't ideal-- we should have more direct indications
                        !that these are valid.

                        if (bca_) then
                                ldm_has_ci = (allocated(self%jackknife_distribution) .and. &
                                             self%bootstrap_distribution_isSorted)
                        else
                                ldm_has_ci = self%bootstrap_distribution_isSorted
                        endif

                endfunction ldm_has_ci

                function ldm_get_ci(self, percentile, bca) result(bound)

                        use core_stat,    only: bca_transform_percentile
                        use core_convert, only: itoa

                        implicit none

                        class(linearDivisionModel),intent(in   ) :: self
                        real    (dp),              intent(in   ) :: percentile
                        logical,         optional, intent(in   ) :: bca

                        real    (dp)                             :: bound

                        logical                                  :: bca_

                        integer (si)                             :: p_index
                        real    (dp)                             :: corr_perc

                        !default arg guards
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

                endfunction ldm_get_ci

                pure function ldm_get_bias(self) result(estimate)

                        implicit none

                        class(linearDivisionModel),intent(in   ) :: self

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

                endfunction ldm_get_bias

                pure subroutine ldm_add_ci_values(self, toAdd, jk) 
                                                               
                        implicit none                          

                        class(linearDivisionModel),         intent(inout) :: self
                        real    (dp),                       intent(in   ) :: toAdd(:)
                        logical,                   optional,intent(in   ) :: jk

                        logical                                  :: jk_

                        integer (si)                             :: iter,copyIter

                        real    (dp),allocatable                 :: scratch(:)

                        !default arg guards
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

                endsubroutine ldm_add_ci_values

                elemental function ldm_get_model_freq(self) result(count)

                        implicit none

                        class(linearDivisionModel),intent(in   ) :: self

                        integer (si)            :: count

                        count = self%frequency

                endfunction ldm_get_model_freq

                pure subroutine ldm_incr_model_freq(self,amount)

                        implicit none

                        class(linearDivisionModel),intent(inout)          :: self
                        integer (si)              ,intent(in   ),optional :: amount

                        if (present(amount)) then
                                self%frequency = self%frequency + amount
                        else 
                                self%frequency = self%frequency + 1
                        endif

                endsubroutine ldm_incr_model_freq

                pure subroutine ldm_decr_model_freq(self,amount)

                        implicit none

                        class(linearDivisionModel),intent(inout)          :: self
                        integer (si)              ,intent(in   ),optional :: amount

                        if (present(amount)) then
                                self%frequency = self%frequency - amount
                        else 
                                self%frequency = self%frequency - 1
                        endif

                endsubroutine ldm_decr_model_freq

                function ldm_get_random_mapping(self,tdata)

                        use core_random,    only: gen_rand_ordered_seq
                        use core_filter,    only: boundariesToMapping
                        use obj_trajectory, only: traj_training_data

                        implicit none

                        class(linearDivisionModel),intent(in   ) :: self
                        class(traj_training_data), intent(in   ) :: tdata

                        integer (si),allocatable                 :: ldm_get_random_mapping(:)
                        integer (si),allocatable                 :: boundaries(:)

                        if ((self%nSites == 0) .or. (tdata%trj%nAtoms == 0)) then
                                allocate(ldm_get_random_mapping(0))
                        else
                                boundaries = gen_rand_ordered_seq(numSeq = self%nSites-1,&
                                                                  minSeq = 1,&
                                                                  maxSeq = tdata%trj%nAtoms)
                                ldm_get_random_mapping = boundariesToMapping(boundaries,tdata%trj%nAtoms)
                        endif

                endfunction ldm_get_random_mapping

                pure function ldm_get_map_accum_type(self)

                        implicit none

                        class(linearDivisionModel),   intent(in   ) :: self

                        character(:),allocatable              :: ldm_get_map_accum_type

                        if (allocated(self%map_accumulator_type)) then
                                ldm_get_map_accum_type = self%map_accumulator_type
                        else
                                allocate(character(len=0) :: ldm_get_map_accum_type)
                        endif

                endfunction ldm_get_map_accum_type

                pure function get_nSites(self)

                        implicit none

                        class(linearDivisionModel),intent(in   ) :: self

                        integer (si)    :: get_nSites

                        get_nSites = self%nSites
                        
                endfunction get_nSites

                pure function ldm_eq(lhs,rhs)

                        implicit none

                        class(linearDivisionModel),intent(in   ) :: lhs
                        class(model),              intent(in   ) :: rhs

                        logical         :: ldm_eq

                        logical         :: lhs_charges_allocated, rhs_charges_allocated        
                        logical         :: lhs_mapping_allocated, rhs_mapping_allocated        

                        real    (dp),parameter :: tolerance = real(10.0_dp**(-14),dp)

                        select type(rhs)
                        class is (linearDivisionModel)
                        if (any(abs(lhs%residuals - rhs%residuals) > tolerance)) then
                                ldm_eq = .false.
                                return
                        endif

                        if (lhs%nAtoms /= rhs%nAtoms) then
                                ldm_eq = .false.
                                return
                        endif

                        if (lhs%nSites /= rhs%nSites) then
                                ldm_eq = .false.
                                return
                        endif

                        lhs_charges_allocated = allocated(lhs%minimum_cg_charges)
                        rhs_charges_allocated = allocated(rhs%minimum_cg_charges)

                        if (lhs_charges_allocated .and. rhs_charges_allocated) then
                                if (any(abs(rhs%minimum_cg_charges - lhs%minimum_cg_charges ) > tolerance )) then
                                        ldm_eq = .false.
                                        return
                                endif
                        else if (lhs_charges_allocated) then
                                ldm_eq = .false.
                                return
                        else if (rhs_charges_allocated) then
                                ldm_eq = .false.
                                return
                        endif

                        lhs_mapping_allocated = allocated(lhs%minimum_mapping)
                        rhs_mapping_allocated = allocated(rhs%minimum_mapping)

                        if (lhs_mapping_allocated .and. rhs_mapping_allocated) then
                                if (any(.not. (rhs%minimum_mapping == lhs%minimum_mapping))) then
                                        ldm_eq = .false.
                                        return
                                endif
                        else if (lhs_mapping_allocated) then
                                ldm_eq = .false.
                                return
                        else if (rhs_mapping_allocated) then
                                ldm_eq = .false.
                                return
                        endif

                        ldm_eq = .true.
                        return
                        end select

                        ldm_eq = .false.
                        return

                endfunction ldm_eq
endmodule obj_linearDivisionModel
