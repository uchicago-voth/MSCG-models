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
! Provides the definition of a centroid model, which takes representative points in space and 
! assigns all nearby atoms to them as a means to generating noncontigous (in primary sequence) 
! maps.

module obj_centroidModel

        use env_kindtypes, only: si, sp, dp
        use abs_obj_model, only: model,model_config
        use fit_common,    only: num_residuals
        use obj_Rloggers,  only: t_rLogger
        use obj_trajectory,only: traj

        implicit none

        private

        public          :: centroidModel

        type, extends(model) :: centroidModel
                private

                !fund. state variables
                logical                  :: tainted = .false.
                logical                  :: record  = .false.

                integer (si)             :: frequency = 1

                !mapping state variables
                integer (si)             :: nSites = 0
                integer (si)             :: nAtoms = 0   !input number of sites
                real    (dp),allocatable :: minimum_centroids(:,:)
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

                !model specific paramters
                real    (dp)             :: cluster_density_tol
                integer (si)             :: kmeans_max_iter

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
                        procedure,public :: reset         => cm_reset
                        procedure,public :: anneal        => cm_sim_annealing
                        procedure,public :: refine        => cm_lucky_descent
                        procedure,public :: summarize     => cm_summarize

                        procedure,public :: recalcResidual=> cm_recalc_residual

                        procedure,public :: getAnnealLog  => cm_get_anneal_log
                        procedure,public :: getRefineLog  => cm_get_refine_log

                        procedure,public :: getMapping    => cm_get_mapping
                        procedure,public :: getCentroids  => cm_get_centroids

                        procedure,public :: isMapChargeSplit => cm_is_mapping_charge_split
                        procedure,public :: hasParentMap  => cm_has_parent_mapping

                        procedure        :: setAnnealLog  => cm_set_anneal_log
                        procedure        :: setRefineLog  => cm_set_refine_log
                        procedure,public :: getSiteCharges =>cm_get_siteCharges
                        procedure,public :: getResidual    =>cm_get_residual
                        procedure,public :: getFullResidual=>cm_get_full_residual
                        procedure,public :: initialize    => cm_initialize
                        procedure,public :: configure     => cm_configure

                        procedure,public :: isTainted     => cm_isTainted

                        procedure        :: setTaint      => cm_setTainted

                        procedure,public :: getModelFreq  => cm_get_model_freq
                        procedure,public :: incrModelFreq => cm_incr_model_freq
                        procedure,public :: decrModelFreq => cm_decr_model_freq

                        procedure,public :: initializeCI  => cm_init_ci

                        procedure,public :: calcCI        => cm_calc_ci
                        procedure,public :: hasCI         => cm_has_ci

                        procedure,public :: getCI         => cm_get_ci

                        procedure,public :: getBias       => cm_get_bias

                        procedure,public :: addCIValues   => cm_add_ci_values

                        procedure,public :: getRandomMapping => cm_get_random_mapping

                        procedure,public :: getMapAccumType => cm_get_map_accum_type

                        procedure,public :: getNSites     => get_nSites

                        procedure        :: clone         => cm_clone
                        procedure        :: eq            => cm_eq
        end type centroidModel

        contains
                pure subroutine cm_configure(self,config)

                        implicit none

                        class(centroidModel),intent(inout) :: self
                        type(model_config),  intent(in   ) :: config

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

                        self%cluster_density_tol= config%cluster_density_tol

                        self%kmeans_max_iter    = config%kmeans_max_iter

                        self%record = config%record

                        if (self%record) then
                                call self%sa_logger%configure(config%size_buffer_anneal)
                                call self%ld_logger%configure(config%size_buffer_refine)
                        endif

                        self%nSites = config%nSites

                        self%postproc_split_on_charge = config%mapping_charge_split

                endsubroutine cm_configure

                pure subroutine cm_initialize(self)

                        implicit none

                        class(centroidModel),intent(inout) :: self

                        allocate(self%minimum_centroids(self%nSites,3))
                        allocate(self%minimum_cg_charges(self%nSites))

                        if (self%record) then
                                call self%sa_logger%initialize()
                                call self%ld_logger%initialize()
                        endif

                endsubroutine cm_initialize

                pure subroutine cm_reset(self)

                        implicit none

                        class(centroidModel), intent(inout) :: self

                        self%residuals = -1

                        call self%sa_logger%reset(shallow=.true.)
                        call self%ld_logger%reset(shallow=.true.)

                endsubroutine cm_reset

                subroutine cm_sim_annealing(self,tdata,acceptance,descentControl)

                        use obj_trajectory,     only: traj_training_data,labeledTraj
                        use routines_centroidStepwiseOptimization,&
                                                only: centroidSimulatedAnnealing
                        use core_kmeans,        only: sget_labels
                        use core_filter,        only: valueCollapse
                        use fit_common,         only: fit_residuals
                        use obj_ll,             only: i_sp_dll

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel),               intent(inout) :: self
                        class(traj_training_data),          intent(in   ) :: tdata
                        real (dp),                optional, intent(in   ) :: acceptance,descentControl
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !buffer for final mapping calculation
                        real    (dp),  allocatable    :: dbuffer(:,:)
                        type(i_sp_dll),allocatable    :: backmapping(:)

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

                        call centroidSimulatedAnnealing(centroids = self%minimum_centroids,&
                                                     siteCharges  = self%minimum_cg_charges, &
                                                           logger = self%sa_logger,&
                                                     residualList = self%residuals, &
                                                 residual_offsets = self%getNormOffsets(), &
                                                residual_scalings = self%getNormScalings(), &
                                                            tdata = tdata, &
                                         mapping_accumulator_type = self%getMapAccumType(),&
                                                 residual_weights = self%residual_weights,&
                                                   startingAccept = self%sa_accept, &
                                                         tempRate = self%sa_descent_control,&
                                                      cluster_tol = self%cluster_density_tol,&
                                                 cluster_max_iter = self%kmeans_max_iter)

                        call sget_labels(self%minimum_mapping,tdata%trj%refAvg,&
                                         self%minimum_centroids,dbuffer)

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

                endsubroutine cm_sim_annealing

                subroutine cm_lucky_descent(self,tdata,maxSteps)

                        use obj_trajectory, only: traj_training_data
                        use routines_boundaryStepwiseOptimization,&
                                            only: boundaryStepwiseDescent
                        use core_filter,    only: boundariesToMapping

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel),              intent(inout) :: self
                        class(traj_training_data),         intent(in   ) :: tdata
                        integer (si),             optional,intent(in   ) :: maxSteps
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

                endsubroutine cm_lucky_descent

                pure subroutine cm_set_anneal_log(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(inout)          :: self
                        real    (sp),               intent(in   )          :: toSet(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        call self%sa_logger%setLog(toSet)

                endsubroutine cm_set_anneal_log

                pure function cm_get_anneal_log(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (sp),allocatable  :: cm_get_anneal_log(:)

                        if (self%sa_logger%isActive()) then
                            call self%sa_logger%sgetLog(cm_get_anneal_log)
                        else
                            allocate(cm_get_anneal_log(0))
                        endif

                endfunction cm_get_anneal_log

                pure subroutine cm_set_refine_log(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(inout)          :: self
                        real    (sp),               intent(in   )          :: toSet(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        call self%ld_logger%setLog(toSet)

                endsubroutine cm_set_refine_log

                pure function cm_get_refine_log(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (sp),allocatable  :: cm_get_refine_log(:)

                        if (self%ld_logger%isActive()) then
                            call self%ld_logger%sgetLog(cm_get_refine_log)
                        else
                            allocate(cm_get_refine_log(0))
                        endif

                endfunction cm_get_refine_log

                elemental function cm_get_residual(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp)                 :: cm_get_residual

                        cm_get_residual = self%residuals(1)

                endfunction cm_get_residual

                pure function cm_get_full_residual(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp),allocatable  :: cm_get_full_residual(:)

                        cm_get_full_residual = self%residuals

                endfunction cm_get_full_residual

                pure function cm_get_mapping(self,parent_mapping,split_on_charge)

                        use core_filter,    only: boundariesToMapping, composeMapping, &
                                                  fragmentMapping, signDiscretize

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel),         intent(in   )          :: self
                        logical,             optional,intent(in   )          :: parent_mapping
                        logical,             optional,intent(in   )          :: split_on_charge
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si),allocatable :: cm_get_mapping(:)
                        logical                  :: parent_mapping_
                        integer (si),allocatable :: charge_mapping(:)
                        logical                  :: do_split_on_charge

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
                                allocate(cm_get_mapping(0))
                                return
                        endif

                        if (parent_mapping_ .and. allocated(self%parent_mapping)) then
                                cm_get_mapping = composeMapping(self%minimum_mapping,self%parent_mapping)

                                if (do_split_on_charge .and. allocated(self%parent_charges)) then
                                        charge_mapping = signDiscretize(self%parent_charges)
                                        cm_get_mapping = fragmentMapping(cm_get_mapping,charge_mapping)
                                endif
                        else
                                cm_get_mapping = self%minimum_mapping

                                if (do_split_on_charge .and. allocated(self%parent_charges)) then
                                        charge_mapping  = signDiscretize(self%domain_charges)
                                        cm_get_mapping = fragmentMapping(self%minimum_mapping,charge_mapping)
                                endif
                        endif

                endfunction cm_get_mapping

                pure function cm_get_centroids(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable :: cm_get_centroids(:,:)

                        cm_get_centroids = self%minimum_centroids

                endfunction cm_get_centroids

                pure function cm_is_mapping_charge_split(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: cm_is_mapping_charge_split

                        cm_is_mapping_charge_split = self%postproc_split_on_charge

                endfunction cm_is_mapping_charge_split

                pure function cm_has_parent_mapping(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: cm_has_parent_mapping

                        if (allocated(self%parent_mapping)) then
                                cm_has_parent_mapping = .true.
                        else
                                cm_has_parent_mapping = .false.
                        endif

                endfunction cm_has_parent_mapping

                pure function cm_get_siteCharges(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real (dp),allocatable :: cm_get_siteCharges(:)

                        cm_get_siteCharges = self%minimum_cg_charges

                endfunction cm_get_siteCharges

                pure function cm_isTainted(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical :: cm_isTainted

                        cm_isTainted = self%tainted

                endfunction cm_isTainted

                pure subroutine cm_setTainted(self,toSet)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(inout)          :: self
                        logical,                    intent(in   )          :: toSet
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        self%tainted = toSet

                endsubroutine cm_setTainted

                pure function cm_getNAtoms(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si)    :: cm_getNAtoms

                        cm_getNAtoms = self%nAtoms

                endfunction cm_getNAtoms

                pure function cm_getNSites(self)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si)    :: cm_getNSites

                        cm_getNSites = self%nSites

                endfunction cm_getNSites

                function cm_recalc_residual(self,tdata) result(residuals)

                        use fit_common,         only: fit_residuals
                        use obj_trajectory,     only: traj_training_data
                        use obj_ll,             only: i_sp_dll
                        use core_kmeans,        only: sget_labels

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel),           intent(inout) :: self
                        class(traj_training_data),      intent(in   ) :: tdata
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp), allocatable       :: residuals(:)
                        real    (dp), allocatable       :: dbuffer(:,:)

                        type(i_sp_dll), allocatable     :: backmapping(:)

                        allocate(residuals(size(self%residuals)))

                        call sget_labels(self%minimum_mapping,tdata%trj%refAvg,&
                                         self%minimum_centroids,dbuffer)

                        call fit_residuals(tdata        = tdata,             &
                                           mapping      = self%minimum_mapping,   &
                               mapping_accumulator_type = self%getMapAccumType(),   &
                                       residual_weights = self%residual_weights,  &  
                                           cgCharges    = self%minimum_cg_charges,&
                                           backmapping  = backmapping,            &
                                           residualList = residuals,              &
                                           residual_offsets=self%getNormOffsets(),&
                                           residual_scalings=self%getNormScalings())

                endfunction cm_recalc_residual

                subroutine cm_summarize(self,fileID,filename,verbose)

                        use core_convert, only: itoa
                        use core_filter,  only: composeMapping

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(centroidModel), intent(in   )          :: self
                        integer (si),               intent(in   ),optional :: fileID
                        character (len=*),          intent(in   ),optional :: filename
                        logical,                    intent(in   ),optional :: verbose
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        logical         :: verbose_
                        integer (si)    :: file_handle, status

                        !optional arg guards
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
                                        write(fileID,*) "   Proxy mapping:   ", self%minimum_mapping
                                        write(fileID,*) "   Full resolution mapping:   ", &
                                                self%getMapping(parent_mapping=.true.)
        
                                        write(fileID,*) "   Proxy mapping centroids:   ", self%minimum_centroids

                                        write(fileID,*) "   Proxy optimized site charges:   ", &
                                                        self%minimum_cg_charges
                                        write(fileID,*) "   Proxy mapped site charges:   ", &
                                                        self%mapped_cg_charges
                                else
                                        write(fileID,*) "   Mapping:   ", self%minimum_mapping
                                        write(fileID,*) "   Mapping centroids:   ", self%minimum_centroids

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
                                write(fileID,*) "   Mapping:   ", self%minimum_mapping
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

                endsubroutine cm_summarize

                pure subroutine cm_clone(lhs,rhs)

                        implicit none

                        class(centroidModel),intent(inout) :: lhs
                        class(model),              intent(in   ) :: rhs

                        select type (REFrhs => rhs)
                        class is (centroidModel)

                        if (REFrhs%isTainted()) then
                                call lhs%setTaint(.true.)
                                return
                        else
                                call lhs%setTaint(.false.)
                        endif

                        lhs%nAtoms = REFrhs%nAtoms
                        lhs%nSites = REFrhs%nSites

                        lhs%frequency = REFrhs%frequency

                        call lhs%setAnnealLog(REFrhs%getAnnealLog())
                        call lhs%setRefineLog(REFrhs%getRefineLog())

                        if (allocated(REFrhs%minimum_centroids)) then
                                lhs%minimum_centroids = REFrhs%minimum_centroids
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
                endsubroutine cm_clone

                pure subroutine cm_init_ci(self, dist_size, jk)

                        implicit none

                        class(centroidModel),          intent(inout) :: self
                        integer (si),                        intent(in   ) :: dist_size
                        logical,                   optional, intent(in   ) :: jk

                        logical         :: jk_

                        !optional arg guards
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

                endsubroutine cm_init_ci

                subroutine cm_calc_ci(self,bca)

                        use core_sort,          only: qsort
                        use core_stat,          only: get_bca_bias, get_bca_accel

                        implicit none

                        class(centroidModel),           intent(inout) :: self
                        logical,                    optional, intent(in   ) :: bca

                        logical                                  :: bca_

                        !optional arg guards
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

                endsubroutine cm_calc_ci

                pure function cm_has_ci(self,bca)

                        implicit none

                        class(centroidModel),intent(in   ) :: self
                        logical,optional,    intent(in   ) :: bca

                        logical :: cm_has_ci

                        logical :: bca_

                        if (present(bca)) then
                                bca_ = bca
                        else
                                bca_ = .false.
                        endif

                        !this isn't ideal-- we should have more direct indications
                        !that these are valid.

                        if (bca_) then
                                cm_has_ci = (allocated(self%jackknife_distribution) .and. &
                                             self%bootstrap_distribution_isSorted)
                        else
                                cm_has_ci = self%bootstrap_distribution_isSorted
                        endif

                endfunction cm_has_ci

                function cm_get_ci(self, percentile, bca) result(bound)

                        use core_stat,    only: bca_transform_percentile
                        use core_convert, only: itoa

                        implicit none

                        class(centroidModel),intent(in   ) :: self
                        real    (dp),              intent(in   ) :: percentile
                        logical,         optional, intent(in   ) :: bca

                        real    (dp)                             :: bound

                        logical                                  :: bca_

                        integer (si)                             :: p_index
                        real    (dp)                             :: corr_perc

                        !optional arg guards
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

                endfunction cm_get_ci

                pure function cm_get_bias(self) result(estimate)

                        implicit none

                        class(centroidModel),intent(in   ) :: self

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

                endfunction cm_get_bias

                pure subroutine cm_add_ci_values(self, toAdd, jk) 
                                                               
                        implicit none                          

                        class(centroidModel),         intent(inout) :: self
                        real    (dp),                       intent(in   ) :: toAdd(:)
                        logical,                   optional,intent(in   ) :: jk

                        logical                                  :: jk_

                        integer (si)                             :: iter,copyIter

                        real    (dp),allocatable                 :: scratch(:)

                        !optional arg guards
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

                endsubroutine cm_add_ci_values

                elemental function cm_get_model_freq(self) result(count)

                        implicit none

                        class(centroidModel),intent(in   ) :: self

                        integer (si)            :: count

                        count = self%frequency

                endfunction cm_get_model_freq

                pure subroutine cm_incr_model_freq(self,amount)

                        implicit none

                        class(centroidModel),intent(inout)          :: self
                        integer (si)              ,intent(in   ),optional :: amount

                        if (present(amount)) then
                                self%frequency = self%frequency + amount
                        else 
                                self%frequency = self%frequency + 1
                        endif

                endsubroutine cm_incr_model_freq

                pure subroutine cm_decr_model_freq(self,amount)

                        implicit none

                        class(centroidModel),intent(inout)          :: self
                        integer (si)              ,intent(in   ),optional :: amount

                        if (present(amount)) then
                                self%frequency = self%frequency - amount
                        else 
                                self%frequency = self%frequency - 1
                        endif

                endsubroutine cm_decr_model_freq

                function cm_get_random_mapping(self,tdata)

                        use core_kmeans,    only: random_centroid_init,sget_labels
                        use obj_trajectory, only: traj_training_data

                        implicit none

                        class(centroidModel),     intent(in   ) :: self
                        class(traj_training_data),intent(in   ) :: tdata

                        integer (si),allocatable                 :: cm_get_random_mapping(:)
                        real    (dp),allocatable                 :: centroids(:,:)
                        real    (dp),allocatable                 :: buffer(:,:)

                        if ((self%nSites == 0) .or. (tdata%trj%nAtoms == 0)) then
                                allocate(cm_get_random_mapping(0))
                        else
                                centroids = random_centroid_init(tdata%trj%refAvg,self%nSites)
                                call sget_labels(cm_get_random_mapping,tdata%trj%refAvg,centroids,buffer)
                        endif

                endfunction cm_get_random_mapping

                pure function cm_get_map_accum_type(self)

                        implicit none

                        class(centroidModel),   intent(in   ) :: self

                        character(:),allocatable              :: cm_get_map_accum_type

                        if (allocated(self%map_accumulator_type)) then
                                cm_get_map_accum_type = self%map_accumulator_type
                        else
                                allocate(character(len=0) :: cm_get_map_accum_type)
                        endif

                endfunction cm_get_map_accum_type

                pure function get_nSites(self)

                        implicit none

                        class(centroidModel),intent(in   ) :: self

                        integer (si)    :: get_nSites

                        get_nSites = self%nSites
                        
                endfunction get_nSites

                pure function cm_eq(lhs,rhs)

                        use core_stat,          only: assignment_difference

                        implicit none

                        class(centroidModel),intent(in   ) :: lhs
                        class(model),              intent(in   ) :: rhs

                        logical         :: cm_eq

                        logical         :: lhs_charges_allocated, rhs_charges_allocated        
                        logical         :: lhs_mapping_allocated, rhs_mapping_allocated        

                        real    (dp),parameter :: tolerance = real(10.0_dp**(-14),dp)

                        select type(rhs)
                        class is (centroidModel)
                        if (any(abs(lhs%residuals - rhs%residuals) > tolerance)) then
                                cm_eq = .false.
                                return
                        endif

                        if (lhs%nAtoms /= rhs%nAtoms) then
                                cm_eq = .false.
                                return
                        endif

                        if (lhs%nSites /= rhs%nSites) then
                                cm_eq = .false.
                                return
                        endif

                        lhs_charges_allocated = allocated(lhs%minimum_cg_charges)
                        rhs_charges_allocated = allocated(rhs%minimum_cg_charges)

                        if (lhs_charges_allocated .and. rhs_charges_allocated) then
                                if (any(abs(rhs%minimum_cg_charges - lhs%minimum_cg_charges ) > tolerance )) then
                                        cm_eq = .false.
                                        return
                                endif
                        else if (lhs_charges_allocated) then
                                cm_eq = .false.
                                return
                        else if (rhs_charges_allocated) then
                                cm_eq = .false.
                                return
                        endif

                        lhs_mapping_allocated = allocated(lhs%minimum_mapping)
                        rhs_mapping_allocated = allocated(rhs%minimum_mapping)

                        if (lhs_mapping_allocated .and. rhs_mapping_allocated) then
                                if (0 /= assignment_difference(rhs%minimum_mapping,lhs%minimum_mapping)) then
                                        cm_eq = .false.
                                        return
                                endif
                        else if (lhs_mapping_allocated) then
                                cm_eq = .false.
                                return
                        else if (rhs_mapping_allocated) then
                                cm_eq = .false.
                                return
                        endif

                        cm_eq = .true.
                        return
                        end select

                        cm_eq = .false.
                        return

                endfunction cm_eq
endmodule obj_centroidModel
