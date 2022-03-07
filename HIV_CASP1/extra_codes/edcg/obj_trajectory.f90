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
! Defines trajectory objects, which hold CG and AA (atomistic) trajectories.
! Additionally, defintes traj_preproc_config objects, which hold the configurations
! options neccesary for trajectory preprocessing.

module obj_trajectory

        use env_kindtypes, only: si,dp

        implicit none

        private

        public          :: traj, labeledTraj, trajPreproc_config, traj_training_data

        !trajectory preprocessing config object. Succintly provides information for processing.
        type                             :: trajPreproc_config

                !controls for projecting the trajectory
                logical                   :: proj_flag
                integer (si)              :: proj_num_dof
                !character(:),allocatable  :: proj_type
                !character(:),allocatable  :: proj_dilation

                logical                   :: align_flag = .true.
                !character(:),allocatable  :: align_type

                !character(:),allocatable  :: gedcg_diss_mat_type
                character(len=256)        :: gedcg_diss_mat_type = 'squareDisplacement'

                logical                   :: pairwise_variance_distance_flag

        end type trajPreproc_config

        !Basic trajectory object. Unlabeled.
        type                             :: traj

                integer (si)              :: nAtoms      = -1 !Num of atoms.
                integer (si)              :: nDimensions = -1 !Num of dimensions
                integer (si)              :: nSteps      = -1 !Num of timesteps.
                integer (si)              :: stride      = -1 !Num skipped steps

                real    (dp),allocatable  :: refAvg(:,:)      !Reference struct.

                real    (dp),allocatable  :: coord(:,:,:)     !site coordinates.
                                                              !atom, dimension, step
                real    (dp),allocatable  :: atomCharges(:)

                logical,                 private  :: is_mapped_trajectory = .false.
                integer (si),allocatable,private  :: parent_mapping(:)
                real    (dp),allocatable,private  :: parent_charges(:)

                contains
                        procedure :: writeXYZfile         => write_xyz_trajectory
                        procedure :: writeCSVfile         => write_csv_trajectory
                        procedure :: writeAverageXYZfile  => write_xyz_average
                        procedure :: bootstrapCopy        => traj_bootstrapCopy

                        procedure :: iterativeAlign       => traj_align_to_avg

                        procedure :: strip                => traj_atom_strip

                        procedure :: genAvgStruct         => traj_gen_avg_coord

                        procedure :: genROGtraj           => traj_gen_rog_traj

                        procedure :: genDipoleTraj        => traj_omp_gen_dipole_traj
                        procedure :: genQuadrupoleTraj    => traj_omp_gen_quadrupole_traj

                        procedure :: mapTrajectory        => traj_map_trajectory

                        procedure :: stack                => traj_stack_traj

                        procedure :: isMapped             => isMapped
                        procedure :: getParentMap         => getParentMap
                        procedure :: getParentCharges     => getParentCharges

                        procedure :: reset                => traj_reset
                        final     :: traj_deallocate
        end type traj

        interface traj
                procedure traj_traj_list_construct
        end interface

        !Adds arrays of types and labels. Used for pulling data from
        !a trajectory file.
        type, extends(traj)             :: labeledTraj

                integer (si)              :: nRes = -1 

                !list different types of atoms present.
                character (LEN=7),allocatable  :: site_labels(:)
                integer (si),allocatable  :: site_type(:)
                integer (si),allocatable  :: site_res(:)
                real    (dp),allocatable  :: site_mass(:)

                integer (si),allocatable  :: residue_start_site(:)

                contains
                        procedure :: mapTrajectory      => ltraj_map_trajectory
                        procedure :: lmapTrajectory     => ltraj_ltraj_map_trajectory
                        procedure :: namedMapTrajectory => ltraj_ltraj_named_map_trajectory
                        procedure :: addPsfFile         => ltraj_add_psf_file
                        procedure :: addDcdFile         => ltraj_add_dcd_file
                        procedure :: bootstrapCopy      => ltraj_bootstrapCopy

                        procedure :: stack              => ltraj_stack_traj

                        procedure :: strip              => ltraj_atom_strip

                        procedure :: reset              => ltraj_reset
                        final     :: ltraj_deallocate
        end type labeledTraj

        interface labeledTraj
                procedure ltraj_file_construct, ltraj_ltraj_list_construct
        end interface

        !Adds arrays of types and labels. Used for pulling data from
        !a trajectory file.
        type                            :: traj_stats

                integer (si)              :: edDOF      = -1 !needed for pcamats.
                real    (dp)              :: edcgNorm 
                real    (dp)              :: atomEspVar = -1.0_dp

                real    (dp),allocatable  :: distanceMat(:,:)

                real    (dp),allocatable  :: pcaMat(:,:)        !covariance matrix of 3N particles.
                real    (dp),allocatable  :: isoPcaMat(:,:)     !trace collapsed cov mat.

                logical     ,private              :: is_calculated_refAvg_diameter = .false.
                real    (dp),private              :: refAvg_diameter
                real    (dp),private,allocatable  :: refAvg_distance_mat(:,:)
                real    (dp),private,allocatable  :: refAvg_distance_var_mat(:,:)

                real    (dp),private,allocatable  :: custom_distance_mat(:,:)
                character(:),private,allocatable  :: custom_distance_mat_type
                real    (dp),private,allocatable  :: custom_distance_mat_param(:)

                real    (dp),private,allocatable  :: sqChargeDiffMat(:,:)

                contains
                        procedure :: computeEDmat               => traj_stats_calculate_covar
                        procedure :: computeESPvar              => traj_stats_calculate_espvar

                        procedure :: computeAvgSDistanceMat     => traj_stats_calculate_avg_dmat
                        procedure :: computeAvgSDiameter        => traj_stats_calculate_avg_struct_diameter
                        procedure :: getAvgSDistanceMat         => traj_stats_get_avg_struct_dmat
                        procedure :: getAvgSDiameter            => traj_stats_get_avg_struct_diameter

                        procedure :: getAvgSDistanceVarMat      => traj_stats_get_avg_struct_var_dmat

                        procedure :: computeAccumulatedDistanceMatrix &
                                                                => traj_stats_calculate_accumDistanceMatrix
                        procedure :: computeCustomDistMat       => traj_stats_calculate_customDistMat
                        procedure :: getCustomDistMat           => traj_stats_get_customDistMat

                        procedure :: computePowChargeDiffMat     => traj_stats_calculate_pow_charge_diff_mat
                        procedure :: getPowChargeDiffMat         => traj_stats_get_pow_charge_diff_mat
        end type traj_stats

        !Struct to hold calculations related to a specific trajectory.
        !It should be passed as the fundamental object for trajectory
        !analysis.
        !
        !This divorces the large amount of statistics from 
        type                            :: traj_training_data

                type(traj_stats)                :: stats
                type(labeledTraj)               :: trj
                type(trajPreproc_config)        :: preproc_design

                contains
                        procedure :: preprocess => traj_training_data_preproc

        end type traj_training_data

        contains
                !trajectory deconstructor.
                pure subroutine traj_deallocate(self)

                        implicit none

                        type (traj), intent(inout)   :: self

                        if (allocated(self%refAvg))         deallocate(self%refAvg)
                        if (allocated(self%coord))          deallocate(self%coord)
                        if (allocated(self%atomCharges))    deallocate(self%atomCharges)
                        if (allocated(self%parent_mapping)) deallocate(self%parent_mapping)
                        if (allocated(self%parent_charges)) deallocate(self%parent_charges)

                endsubroutine traj_deallocate

                !trajectory deconstructor.
                pure subroutine ltraj_deallocate(self)

                        implicit none

                        type (labeledTraj), intent(inout)   :: self

                        if (allocated(self%refAvg))         deallocate(self%refAvg)
                        if (allocated(self%coord))          deallocate(self%coord)
                        if (allocated(self%atomCharges))    deallocate(self%atomCharges)
                        if (allocated(self%parent_mapping)) deallocate(self%parent_mapping)
                        if (allocated(self%parent_charges)) deallocate(self%parent_charges)

                        if (allocated(self%site_labels))        deallocate(self%site_labels)
                        if (allocated(self%site_type))          deallocate(self%site_type)
                        if (allocated(self%site_res))           deallocate(self%site_res)
                        if (allocated(self%site_mass))          deallocate(self%site_mass)
                        if (allocated(self%residue_start_site)) deallocate(self%residue_start_site)

                endsubroutine ltraj_deallocate

                !trajectory deconstructor.
                pure subroutine traj_reset(self)

                        implicit none

                        class (traj), intent(inout)   :: self

                        !Destructors aren't polymorphic, so we have to select
                        !on typ here.
                        select type (self)
                        type is (traj)
                                call traj_deallocate(self)

                                self%nAtoms      = -1
                                self%nDimensions = -1
                                self%nSteps      = -1
                                self%stride      = -1
                        end select

                endsubroutine traj_reset

                !trajectory deconstructor.
                pure subroutine ltraj_reset(self)

                        implicit none

                        class (labeledTraj), intent(inout)   :: self

                        !Destructors aren't polymorphic, so we have to select
                        !on typ here.
                        select type (self)
                        type is (labeledTraj)
                                call ltraj_deallocate(self)

                                self%nAtoms      = -1
                                self%nDimensions = -1
                                self%nSteps      = -1
                                self%stride      = -1
                                self%nRes        = -1
                        end select

                endsubroutine ltraj_reset

                !trajectory deconstructor.
                function traj_traj_list_construct(trj_list) result(trj)

                        use env_kindtypes, only: si

                        implicit none

                        type (traj),intent(in   )   :: trj_list(:)
                        type (traj)                 :: trj

                        integer (si) :: iter, stack_status

                        if (size(trj_list) <= 0) then
                                return
                        endif

                        trj = trj_list(1)

                        if (size(trj_list) > 1) then
                                do iter=2,size(trj_list)
                                        call trj%stack(trj_list(iter),stack_status)
                                        if (stack_status > 0) then
                                                call trj%reset()
                                                return
                                        endif
                                enddo
                        endif

                endfunction traj_traj_list_construct

                !trajectory deconstructor.
                function ltraj_file_construct(psf_file,dcd_file) result(ltraj)

                        implicit none

                        character(*),intent(in   )   :: dcd_file,psf_file
                        type (labeledTraj)           :: ltraj

                        call ltraj%addPsfFile(psf_file)
                        call ltraj%addDcdFile(dcd_file)

                endfunction ltraj_file_construct

                !trajectory deconstructor.
                function ltraj_ltraj_list_construct(ltraj_list) result(ltraj)

                        use env_kindtypes, only: si

                        implicit none

                        type (labeledTraj),intent(in   ) :: ltraj_list(:)
                        type (labeledTraj)               :: ltraj

                        integer (si) :: iter, stack_status

                        if (size(ltraj_list) <= 0) then
                                return
                        endif

                        ltraj = ltraj_list(1)

                        if (size(ltraj_list) > 1) then
                                do iter=2,size(ltraj_list)
                                        call ltraj%stack(ltraj_list(iter),stack_status)
                                        if (stack_status > 0) then
                                                call ltraj%reset()
                                                return
                                        endif
                                enddo
                        endif

                endfunction ltraj_ltraj_list_construct

                !check if trajectory is derived by mapping a different directory
                pure function isMapped(self)

                        implicit none

                        class   (traj), intent(in   )   :: self

                        logical         :: isMapped

                        isMapped = self%is_mapped_trajectory

                endfunction isMapped

                !get map that sends parent sites to this trajectory's sites.
                pure function getParentMap(self)

                        implicit none

                        class   (traj), intent(in   )   :: self

                        integer (si),allocatable :: getParentMap(:)

                        if (self%is_mapped_Trajectory .and. allocated(self%parent_mapping)) then
                                getParentMap = self%parent_mapping
                        else
                                allocate(getParentMap(0))
                        endif

                endfunction getParentMap

                !get map that sends parent sites to this trajectory's sites.
                pure function getParentCharges(self)

                        implicit none

                        class   (traj), intent(in   )   :: self

                        real    (dp),allocatable :: getParentCharges(:)

                        if (self%is_mapped_Trajectory .and. allocated(self%parent_charges)) then
                                getParentCharges = self%parent_charges
                        else
                                allocate(getParentCharges(0))
                        endif

                endfunction getParentCharges

                !compute the covariance matrix in a trajectory object.
                !If we're given filterDof, use that to filter the pca matrix.
                !if not, default to no filtering.
                subroutine traj_stats_calculate_covar(self,trj,filterDof)

                        use routines_trajectory, only: compute_covar

                        implicit none

                        class   (traj_stats),intent(inout)          :: self
                        class   (traj),      intent(in   )          :: trj
                        integer (si),        intent(in   ),optional :: filterDof

                        if( .not. present(filterDof) ) then
                                self%edDof = trj%nAtoms
                        else if (filterDof < 1) then
                                self%edDof = trj%nAtoms
                        else
                                self%edDof = filterDof
                        endif

                        !allocate or clear pca matrices
                        if (allocated(self%isoPcaMat)) then
                                self%isoPcaMat = 0
                        else
                                allocate(self%isoPcaMat(trj%nAtoms,&
                                                trj%nAtoms))
                                self%isoPcaMat = 0
                        endif

                        if (allocated(self%pcaMat)) then
                                self%pcaMat = 0
                        else
                                allocate(self%pcaMat(3*trj%nAtoms,&
                                                3*trj%nAtoms))
                                self%pcaMat = 0
                        endif

                        !calculate pca matrix
                        call compute_covar(trj%coord,   &
                                           trj%refAvg,  &
                                           self%pcaMat,  &
                                           self%edcgNorm,&
                                           self%edDof,   &
                                           self%isoPcaMat)

                endsubroutine traj_stats_calculate_covar

                pure function traj_gen_avg_coord(self) result(avg_coord)

                        use env_kindtypes,      only: si, dp
                        use obj_accum,          only: accumM

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj),       intent(in   )          :: self
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp), allocatable :: avg_coord(:,:)

                        !local variables
                        integer (si)    :: step
                        type(accumM)    :: meanStructAccum

                        if (.not. allocated(self%coord)) then
                                allocate(avg_coord(0,0))
                                return
                        endif

                        call meanStructAccum%reset()
                        do step=1,size(self%coord,3)
                                call meanStructAccum%add(self%coord(:,:,step))
                        enddo

                        avg_coord = meanStructAccum%getMean()

                endfunction traj_gen_avg_coord

                pure subroutine traj_align_to_avg(self,center,use_cached_avg)

                        use routines_trajectory, only: align_to_avg
                        use obj_accum,           only: accumM

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj),       intent(inout)          :: self
                        logical,              intent(in   ),optional :: center
                        logical,              intent(in   ),optional :: use_cached_avg
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional var couterparts
                        logical         :: center_
                        logical         :: use_cached_avg_

                        !local variables
                        integer (si)    :: step
                        type(accumM)    :: meanStructAccum

                        if (present(center)) then
                                center_ = center
                        else
                                center_ = .true.
                        endif

                        if (present(use_cached_avg)) then
                                use_cached_avg_ = use_cached_avg
                        else
                                use_cached_avg_ = .false.
                        endif

                        if (.not. use_cached_avg_) then
                                call meanStructAccum%reset()
                                do step=1,size(self%coord,3)
                                        call meanStructAccum%add(self%coord(:,:,step))
                                enddo
                                self%refAvg = meanStructAccum%getMean()
                        endif

                        call align_to_avg(coord  = self%coord, &
                                          refAvg = self%refAvg,&
                                          center = center_)

                endsubroutine traj_align_to_avg

                !Modifies a trajectory in place to remove atoms which
                !have a .false. entry in @mask.
                pure subroutine traj_atom_strip(self,mask)

                        use core_stat, only: which

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj),       intent(inout)          :: self
                        logical,              intent(in   )          :: mask(:)
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si),allocatable        :: indices(:)

                        if (.not. any(mask)) return

                        indices = which(mask)

                        self%nAtoms = size(indices)
                        if (allocated(self%refAvg)) deallocate(self%refAvg)

                        if (allocated(self%coord)) self%coord = self%coord(indices,:,:)

                        if (allocated(self%atomCharges)) self%atomCharges &
                                = self%atomCharges(indices)

                        self%is_mapped_trajectory = .false.

                        if (allocated(self%parent_mapping)) deallocate(self%parent_mapping)
                        if (allocated(self%parent_charges)) deallocate(self%parent_charges)

                endsubroutine traj_atom_strip

                !Adds a trajectory to the current trajectory, e.g. adds a molecule.
                !
                !If add_trj is malformed, self will likely become malformed.
                !The routing should not segfault, though.
                !
                !Stride of added trj is ignored.
                !
                ! @self:        trajectory to be modified
                ! @add_trj:     trajectory to add to @self
                ! @status:      return code explaining results
                pure subroutine traj_stack_traj(self,add_trj,status)

                        use routines_trajectory, only: align_to_avg
                        use obj_accum,           only: accumM

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj),         intent(inout) :: self
                        class   (traj),         intent(in   ) :: add_trj
                        integer (si),  optional,intent(inout) :: status
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !temp variables
                        real    (dp),allocatable    :: scratch_real_3array(:,:,:)
                        real    (dp),allocatable    :: scratch_real_2array(:,:)
                        real    (dp),allocatable    :: scratch_real_1array(:)

                        integer (si)                :: comb_nAtoms

                        if (present(status)) status = 0

                        !Check to see if new traj is incompatible
                        if (self%nSteps /= add_trj%nSteps) then
                                if (present(status)) status = 1
                                return
                        endif

                        if (self%nDimensions /= add_trj%nDimensions) then
                                if (present(status)) status = 2
                                return
                        endif

                        comb_nAtoms = self%nAtoms + add_trj%nAtoms

                        !Move over coordinates if valid.
                        if (allocated(self%coord) .and. allocated(add_trj%coord)) then
                                !temporarily store host traj's coords
                                call move_alloc(self%coord,scratch_real_3array)

                                allocate(self%coord(comb_nAtoms,self%nDimensions,self%nSteps))

                                !transfer old coordinates
                                self%coord(1:self%nAtoms,:,:) = scratch_real_3array

                                deallocate(scratch_real_3array)

                                !transfer add trj coordinates
                                self%coord(self%nAtoms+1:,:,:) = add_trj%coord

                        else
                                if (present(status)) status = 3
                                return
                        endif

                        !Move over charges if valid. Ommited sets charges are set at 0
                        !if other charges are present.
                        if (allocated(self%atomCharges) .and. allocated(add_trj%atomCharges)) then
                                !temporarily store host traj's atomCharges
                                call move_alloc(self%atomCharges,scratch_real_1array)

                                allocate(self%atomCharges(comb_nAtoms))

                                !transfer old charges
                                self%atomCharges(1:self%nAtoms) = scratch_real_1array

                                deallocate(scratch_real_1array)

                                !transfer add trj coordinates
                                self%atomCharges(self%nAtoms+1:) = add_trj%atomCharges

                        else if (allocated(self%atomCharges)) then
                                !when we have host changes, but not add charges

                                !temporarily store host traj's atomCharges
                                call move_alloc(self%atomCharges,scratch_real_1array)

                                allocate(self%atomCharges(comb_nAtoms))

                                !transfer old charges
                                self%atomCharges(1:self%nAtoms) = scratch_real_1array

                                deallocate(scratch_real_1array)

                                !transfer add trj coordinates
                                self%atomCharges(self%nAtoms+1:) = 0

                        else if (allocated(add_trj%atomCharges)) then
                                !when we have add charges, but no host charges

                                allocate(self%atomCharges(comb_nAtoms))

                                !transfer old charges
                                self%atomCharges(1:self%nAtoms) = 0

                                !transfer add trj coordinates
                                self%atomCharges(self%nAtoms+1:) = add_trj%atomCharges

                        endif

                        !Move over refAvg if BOTH trajectories have it.
                        if (allocated(self%refAvg) .and. allocated(add_trj%refAvg)) then
                                !temporarily store host traj's refAvg
                                call move_alloc(self%refAvg,scratch_real_2array)

                                allocate(self%refAvg(comb_nAtoms,self%nDimensions))

                                !transfer old refAvg
                                self%refAvg(1:self%nAtoms,:) = scratch_real_2array

                                deallocate(scratch_real_2array)

                                !transfer add trj refAvg
                                self%refAvg(self%nAtoms+1:,:) = add_trj%refAvg

                        else
                                !We have to chuck the average we have.
                                if (allocated(self%refAvg)) deallocate(self%refAvg)
                        endif

                        !the usage of these variables is questionable, so it's set to false.
                        !how could be backreference to two distinct trajectories?
                        self%is_mapped_trajectory = .false.

                        !Update scalar values
                        self%nAtoms = comb_nAtoms

                endsubroutine traj_stack_traj

                !Adds a trajectory to the current trajectory, e.g. adds a molecule, for ltrajs.
                !
                !If add_trj is malformed, self will likely become malformed.
                !The routing should not segfault, though.
                !
                !Stride of added trj is ignored.
                !
                ! @self:        trajectory to be modified
                ! @add_trj:     trajectory to add to @self
                ! @status:      return code explaining results
                pure subroutine ltraj_stack_traj(self,add_trj,status)

                        use routines_trajectory, only: align_to_avg
                        use obj_accum,           only: accumM

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (labeledtraj),  intent(inout) :: self
                        class   (traj),         intent(in   ) :: add_trj
                        integer (si),  optional,intent(inout) :: status
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !temp variables
                        real    (dp),     allocatable  :: scratch_real_1array(:)
                        character (LEN=4),allocatable  :: scratch_string_1array(:)
                        integer (si),     allocatable  :: scratch_integer_1array(:)

                        !old nAtoms
                        integer (si)                :: old_self_nAtoms, substatus

                        if (present(status)) status = 0

                        select type (add_trj)
                        type is (traj)
                                !there's no reasonable way to combine a traj into an 
                                !ltraj. It is better to fail here.
                                if (present(status)) status = 100
                                return
                        type is (labeledTraj)

                                old_self_nAtoms = self%nAtoms

                                call traj_stack_traj(self,add_trj,substatus)

                                if (substatus /= 0) then
                                        if (present(status)) status = substatus
                                        return
                                endif

                                if (allocated(self%site_labels) .and. &
                                    allocated(add_trj%site_labels)) then
                                        call move_alloc(self%site_labels,scratch_string_1array)

                                        allocate(self%site_labels(self%nAtoms))

                                        !transfer old charges
                                        self%site_labels(1:old_self_nAtoms) = scratch_string_1array

                                        deallocate(scratch_string_1array)

                                        !transfer add trj coordinates
                                        self%site_labels(old_self_nAtoms+1:) = add_trj%site_labels
                                else if (allocated(self%site_labels)) then
                                        call move_alloc(self%site_labels,scratch_string_1array)

                                        allocate(self%site_labels(self%nAtoms))

                                        !transfer old charges
                                        self%site_labels(1:old_self_nAtoms) = scratch_string_1array

                                        deallocate(scratch_string_1array)

                                        !transfer add trj coordinates
                                        self%site_labels(old_self_nAtoms+1:) = "UNKN"

                                else if (allocated(add_trj%site_labels)) then

                                        allocate(self%site_labels(self%nAtoms))

                                        !transfer old charges
                                        self%site_labels(1:old_self_nAtoms) = "UNKN"

                                        !transfer add trj coordinates
                                        self%site_labels(old_self_nAtoms+1:) = add_trj%site_labels
                                endif

                                if (allocated(self%site_type) .and. &
                                    allocated(add_trj%site_type)) then
                                        call move_alloc(self%site_type,scratch_integer_1array)

                                        allocate(self%site_type(self%nAtoms))

                                        !transfer old charges
                                        self%site_type(1:old_self_nAtoms) = scratch_integer_1array

                                        deallocate(scratch_integer_1array)

                                        !transfer add trj coordinates
                                        self%site_type(old_self_nAtoms+1:) = add_trj%site_type
                                else if (allocated(self%site_type)) then
                                        call move_alloc(self%site_type,scratch_integer_1array)

                                        allocate(self%site_type(self%nAtoms))

                                        !transfer old charges
                                        self%site_type(1:old_self_nAtoms) = scratch_integer_1array

                                        !transfer add trj coordinates
                                        self%site_type(old_self_nAtoms+1:) = maxval(scratch_integer_1array) + 1

                                        deallocate(scratch_integer_1array)

                                else if (allocated(add_trj%site_type)) then

                                        allocate(self%site_type(self%nAtoms))

                                        !transfer old charges
                                        self%site_type(1:old_self_nAtoms) = maxval(add_trj%site_type) + 1

                                        !transfer add trj coordinates
                                        self%site_type(old_self_nAtoms+1:) = add_trj%site_type
                                endif

                                !handle nRes, site_res and residue_start_site at the same time.

                                if (allocated(self%site_res) .and. &
                                    allocated(add_trj%site_res)) then

                                        !set site_res
                                        call move_alloc(self%site_res,scratch_integer_1array)

                                        allocate(self%site_res(self%nAtoms))

                                        self%site_res(1:old_self_nAtoms) = scratch_integer_1array

                                        deallocate(scratch_integer_1array)

                                        self%site_res(old_self_nAtoms+1:) = add_trj%site_res

                                        !set residue_start_site

                                        call move_alloc(self%residue_start_site,scratch_integer_1array)

                                        allocate(self%residue_start_site(self%nRes + add_trj%nRes))

                                        self%residue_start_site(1:self%nRes) = scratch_integer_1array

                                        deallocate(scratch_integer_1array)

                                        self%residue_start_site(self%nRes+1:) = &
                                                add_trj%residue_start_site + old_self_nAtoms

                                        !set nRes
                                        self%nRes = self%nRes + add_trj%nRes

                                else if (allocated(self%site_res)) then

                                        !set site_res
                                        call move_alloc(self%site_res,scratch_integer_1array)

                                        allocate(self%site_res(self%nAtoms))

                                        self%site_res(1:old_self_nAtoms) = scratch_integer_1array

                                        self%site_res(old_self_nAtoms+1:) = maxval(scratch_integer_1array) + 1

                                        deallocate(scratch_integer_1array)

                                        !set residue_start_site
                                        call move_alloc(self%residue_start_site,scratch_integer_1array)

                                        allocate(self%residue_start_site(self%nRes + 1))

                                        self%residue_start_site(1:self%nRes) = scratch_integer_1array

                                        deallocate(scratch_integer_1array)

                                        self%residue_start_site(self%nRes+1:) = old_self_nAtoms + 1

                                        !set nRes
                                        self%nRes = self%nRes + 1

                                else if (allocated(add_trj%site_res)) then

                                        !set site_res
                                        allocate(self%site_res(self%nAtoms))

                                        self%site_res(1:old_self_nAtoms) = 1

                                        self%site_res(old_self_nAtoms+1:) = add_trj%site_res + 1

                                        !set residue_start_site
                                        allocate(self%residue_start_site(add_trj%nRes + 1))

                                        self%residue_start_site(1) = 1

                                        self%residue_start_site(2:) = old_self_nAtoms + add_trj%residue_start_site

                                        !set nRes
                                        self%nRes = add_trj%nRes + 1
                                endif

                                !combine masses if we can. We do not make up masses, as they 
                                !are likely to violate assumptions of mass routines, or are hard to 
                                !detect.
                                if (allocated(self%site_mass) .and. &
                                    allocated(add_trj%site_mass)) then
                                        call move_alloc(self%site_mass,scratch_real_1array)

                                        allocate(self%site_mass(self%nAtoms))

                                        !transfer old charges
                                        self%site_mass(1:old_self_nAtoms) = scratch_real_1array

                                        deallocate(scratch_real_1array)

                                        !transfer add trj coordinates
                                        self%site_mass(old_self_nAtoms+1:) = add_trj%site_mass
                                endif

                        class default
                                if (present(status)) status = 101
                                return
                        end select

                endsubroutine ltraj_stack_traj

                subroutine ltraj_add_psf_file(self,psfFile)

                        use IO_psf, only: read_psf_file

                        implicit none

                        class   (labeledTraj),intent(inout) :: self
                        character (*),   intent(in)    :: psfFile

                        !An ugly call, but psfIO can't be aware of traj
                        !and be referenced here.
                        call read_psf_file(AtomPsfFile = psfFile,      &
                                        nAtoms      = self%nAtoms,     &
                                        nResidues   = self%nRes,       &
                                        atomCharges = self%atomCharges,&
                                        atomMasses  = self%site_mass,  &
                                        atomResidues= self%site_res,   &
                                        atomLabels  = self%site_labels,&
                                        residueStartSites              &
                                            = self%residue_start_site  &
                                        )

                        self%nDimensions = 3

                endsubroutine ltraj_add_psf_file

                subroutine ltraj_add_dcd_file(self,atomDcdFile,compute_meanStruct)

                        use env_kindtypes,      only: si
                        use io_dcd,             only: read_dcd_header, read_dcd_step
                        use obj_accum,          only: accumM

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        class(labeledTraj),intent(inout)          :: self
                        character (*),     intent(in   )          :: atomDcdFile
                        logical,           intent(in   ),optional :: compute_meanStruct

                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !optional arg couterparts
                        logical        :: compute_meanStruct_

                        !local vars
                        integer (si)   :: step!, i
                        integer (si)   :: fp_dcd
                        type(accumM)   :: meanStructAccum

                        if (present(compute_meanStruct)) then
                                compute_meanStruct_ = compute_meanStruct
                        else
                                compute_meanStruct_ = .false.
                        endif

                        ! Read atom dcd header and grab nAtoms, nSteps
                        ! This opens the dcd file too.
                        call read_dcd_header(atomDcdFile,&
                                             self%nAtoms,&
                                             self%nSteps,&
                                             fp_dcd)

                        print*,"INFO: (trajectory) number of Steps:", self%nSteps

                        !allocate space for reading the full trajectory.
                        allocate( self%coord(self%nAtoms, 3, self%nSteps) )

                        !Read in the entire trajectory into our trajectory object.
                        do step = 1, (self%nSteps)
                                call read_dcd_step(self%coord(:,:,step), self%nAtoms, fp_dcd)
                        enddo

                        if (compute_meanStruct_) then
                                if (allocated(self%refAvg)) then
                                        deallocate(self%refAvg)
                                endif

                                allocate(self%refAvg(size(self%coord,1),size(self%coord,2)))

                                call meanStructAccum%reset()
                                do step=1,size(self%coord,3)
                                        call meanStructAccum%add(self%coord(:,:,step))
                                enddo

                                self%refAvg = meanStructAccum%getMean()
                        endif

                        close(fp_dcd)

                endsubroutine ltraj_add_dcd_file

                ! function to coarse a trajectory through mapping or filtering.
                function ltraj_ltraj_named_map_trajectory(self,stride,map_type,pass_res) result(mappedTrajectory)

                        use env_kindtypes,      only: si
                        use core_filter,        only: fragmentMapping, signDiscretize

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (labeledTraj),intent(in   )          :: self
                        character(*),       intent(in   )          :: map_type    
                        integer (si),       intent(in   ),optional :: stride      !number of steps to skip between
                        logical,            intent(in   ),optional :: pass_res
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        type    (labeledTraj)                   :: mappedTrajectory

                        !optional argument mirrors.
                        integer (si)   :: stride_
                        logical        :: pass_res_

                        if (.not. present(stride)) then
                                stride_ = 1
                        elseif (stride <= 0) then
                                stride_ = 1
                        else
                                stride_ = stride
                        endif

                        if (present(pass_res)) then
                                pass_res_ = pass_res
                        else
                                pass_res_ = .true.
                        endif

                        select case (map_type)
                        case ("carbon_alpha")
                                mappedTrajectory = ltraj_ltraj_map_trajectory_CA(self=self,&
                                                                               stride=stride_,&
                                                                             pass_res=pass_res_)
                        case ("residue_cop")
                                mappedTrajectory = ltraj_ltraj_map_trajectory(self=self,&
                                                                               map=self%site_res,&
                                                                            stride=stride_,&
                                                                          pass_res=pass_res_)
                        case ("residue_com")
                                mappedTrajectory = ltraj_ltraj_map_trajectory(self=self,&
                                                                               map=self%site_res,&
                                                                            stride=stride_,&
                                                                  position_weights=self%site_mass,&
                                                                          pass_res=pass_res_)
                        case ("residue_coc")
                                mappedTrajectory = ltraj_ltraj_map_trajectory(self=self,&
                                                                               map=self%site_res,&
                                                                            stride=stride_,&
                                                                  position_weights=self%atomCharges,&
                                                                          pass_res=pass_res_)
                        case ("subresidue_coc")
                                subresidue_coc: block 
                                        integer (si),allocatable :: mapping(:)

                                        allocate(mapping,&
                                                 source=fragmentMapping(signDiscretize(self%atomCharges),&
                                                                        self%site_res))
                                        mappedTrajectory = ltraj_ltraj_map_trajectory(self = self,&
                                                                                       map = mapping,&
                                                                                    stride = stride_,&
                                                                          position_weights = self%atomCharges,&
                                                                                  pass_res = .false.)
                                endblock subresidue_coc
                        case ("subresidue_com")
                                subresidue_com: block 
                                        integer (si),allocatable :: mapping(:)

                                        allocate(mapping,&
                                                 source=fragmentMapping(signDiscretize(self%atomCharges),&
                                                                        self%site_res))
                                        mappedTrajectory = ltraj_ltraj_map_trajectory(self = self,&
                                                                                       map = mapping,&
                                                                                    stride = stride_,&
                                                                          position_weights = self%site_mass,&
                                                                                  pass_res = .false.)
                                endblock subresidue_com
                        case ("identity")
                                call ltraj_stride_copy(to = mappedTrajectory,&
                                                     from = self,&
                                                   stride = stride_)
                        case default
                                print*, "Unknown map type requested: "//map_type//"."
                                stop
                        end select

                endfunction ltraj_ltraj_named_map_trajectory

                !Modifies a trajectory in place to remove atoms which
                !have a .false. entry in @mask.
                !The start places for residues are destroyed in this process.
                pure subroutine ltraj_atom_strip(self,mask)

                        use core_stat, only: which

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (labeledTraj), intent(inout)         :: self
                        logical,               intent(in   )         :: mask(:)
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si),allocatable        :: indices(:)

                        if (.not. any(mask)) return

                        call traj_atom_strip(self,mask)

                        indices = which(mask)

                        if (allocated(self%site_labels)) &
                                self%site_labels = self%site_labels(indices)

                        if (allocated(self%site_type)) &
                                self%site_type = self%site_type(indices)

                        if (allocated(self%site_mass)) &
                                self%site_mass = self%site_mass(indices)

                        !note that we don't canonicize this. There can be gaps
                        !in the numbers which correspond to present residues.
                        if (allocated(self%site_res)) &
                                self%site_res = self%site_res(indices)

                        !these are nontrivial to reconstruct. For now, we simply remove them.
                        !implementation requires us to rescan our residue topology.
                        self%nRes = -1

                        if (allocated(self%residue_start_site)) deallocate(self%residue_start_site)

                endsubroutine ltraj_atom_strip

                pure subroutine ltraj_stride_copy(to,from,stride)

                        use env_kindtypes,      only: si
                        use core_stat,          only: seq

                        implicit none

                        class (labeledTraj),intent(inout)       :: to
                        class (labeledTraj),intent(in   )       :: from
                        integer (si),       intent(in   )       :: stride

                        integer (si),allocatable :: step_seq(:)

                        !normal traj values
                        to%nAtoms = from%nAtoms
                        to%nSteps = from%nSteps
                        to%stride = stride

                        to%is_mapped_trajectory = from%is_mapped_trajectory

                        if (to%is_mapped_trajectory) then
                                to%parent_mapping = from%parent_mapping
                                to%parent_charges = from%parent_charges
                        endif

                        step_seq = seq(1,from%nSteps,stride)

                        to%coord = from%coord(:,:,step_seq)

                        to%nRes  = from%nRes

                        if (allocated(from%site_labels))        to%site_labels = from%site_labels
                        if (allocated(from%site_type))          to%site_type   = from%site_type
                        if (allocated(from%site_res))           to%site_res    = from%site_res
                        if (allocated(from%site_mass))          to%site_mass   = from%site_mass
                        if (allocated(from%residue_start_site)) to%residue_start_site = from%residue_start_site

                endsubroutine ltraj_stride_copy


                subroutine traj_stats_calculate_accumDistanceMatrix(self,trj)

                        use env_kindtypes,      only: si, dp
                        use routines_math,      only: computeDistanceMatrix

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (traj_stats), intent(inout)     :: self
                        class (traj),       intent(in   )     :: trj
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                    :: step
                        real    (dp),allocatable        :: distanceMatrix(:,:) !will be nAtoms:nAtoms
                        real    (dp),allocatable        :: accumDistanceMatrix(:,:) !will be nAtoms:nAtoms

                        allocate(accumDistanceMatrix(trj%nAtoms,trj%nAtoms))
                        accumDistanceMatrix = 0

                        do step=1,trj%nSteps
                                distanceMatrix      = computeDistanceMatrix(trj%coord(:,:,step))
                                accumDistanceMatrix = accumDistanceMatrix - distanceMatrix
                        enddo

                        self%distanceMat = accumDistanceMatrix

                endsubroutine traj_stats_calculate_accumDistanceMatrix

                subroutine traj_stats_calculate_customDistMat(self,distance_type,mix_parameters)

                        use env_kindtypes,      only: dp
                        use routines_trajectory,only: get_square_displacement_distmat

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(traj_stats),     intent(inout)          :: self
                        character (len=*),    intent(in   )          :: distance_type
                        real    (dp),         intent(in   ),optional :: mix_parameters(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        select case (distance_type)
                        case ("squareDisplacement")
                                if (present(mix_parameters)) then
                                        print*, "Custom distance matrix calculation: ignoring mix parameters, &
                                                &as square displacement matrices do not take parameters."
                                endif
                                if (allocated(self%isoPcaMat)) then
                                        self%custom_distance_mat = get_square_displacement_distmat(self%isoPcaMat)
                                        self%custom_distance_mat_type = distance_type
                                else
                                        print*, "FATAL ERROR: Custom distance matrix calculation: &
                                                &no isotropic variance matrix present at time of calcution."
                                        stop
                                endif
                        case default
                                print*, "FATAL ERROR: Custom distance matrix calculation: &
                                        &Unknown distance type in custom distance matrix calculation."
                                print*, "    attempted matrix type:"//distance_type
                                stop
                        end select

                endsubroutine traj_stats_calculate_customDistMat

                subroutine traj_stats_calculate_pow_charge_diff_mat(self,trj,power)

                        use routines_trajectory,only: get_square_displacement_distmat

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(traj_stats),     intent(inout)          :: self
                        class(traj),           intent(in   )          :: trj
                        integer (si),          intent(in   )          :: power
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)            :: iter1, iter2
                        logical                 :: do_alloc

                        do_alloc = .false.

                        if (allocated(self%sqchargediffmat)) then
                                if (.not. ((size(self%sqchargediffmat,1) == size(trj%atomCharges)) .and. &
                                           (size(self%sqchargediffmat,2) == size(trj%atomCharges)))) then
                                        deallocate(self%sqChargeDiffMat)
                                        do_alloc = .true.
                                endif
                        else
                                do_alloc=.true.
                        endif

                        if (do_alloc) then
                                allocate(self%sqChargeDiffMat(size(trj%atomCharges),&
                                                              size(trj%atomCharges)))
                        endif

                        do iter1=1,size(self%sqChargeDiffMat,1)
                        do iter2=iter1+1,size(self%sqChargeDiffMat,1)

                                self%sqChargeDiffMat(iter1,iter2) = &
                                        abs(trj%atomCharges(iter1) - trj%atomCharges(iter1))**power

                                self%sqChargeDiffMat(iter2,iter1) = self%sqChargeDiffMat(iter1,iter2)

                        enddo
                        enddo

                        do iter1=1,size(self%sqChargeDiffMat,1)
                                self%sqChargeDiffMat(iter1,iter1) = 0
                        enddo

                endsubroutine traj_stats_calculate_pow_charge_diff_mat

                pure function traj_stats_get_pow_charge_diff_mat(self) result(mat)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(traj_stats),     intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable        :: mat(:,:)

                        if (allocated(self%sqChargeDiffMat)) then
                                mat = self%sqChargeDiffMat
                        else
                                allocate(mat(0,0))
                        endif

                endfunction traj_stats_get_pow_charge_diff_mat

                pure function traj_stats_get_customDistMat(self) result(dist_mat)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj_stats), intent(in   )          :: self
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp), allocatable                 :: dist_mat(:,:)

                        if (allocated(self%custom_distance_mat)) then
                                dist_mat = self%custom_distance_mat
                        else
                                allocate(dist_mat(0,0))
                        endif

                endfunction traj_stats_get_customDistMat

                pure subroutine traj_stats_calculate_espvar(self,trj,save_DM)

                        use env_kindtypes,      only: si, dp
                        use routines_math,      only: electroPot, computeDistanceMatrix

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj_stats),intent(inout)             :: self
                        class   (traj),      intent(in   )             :: trj 
                        logical,             intent(in   ),optional    :: save_DM !whether to save distance 
                                                                                  !matrix
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variables
                        integer (si)                    :: step
                        real    (dp),allocatable        :: distanceMatrix(:,:) !will be nAtoms:nAtoms
                        real    (dp),allocatable        :: accumDistanceMatrix(:,:) !will be nAtoms:nAtoms
                        real    (dp)                    :: atomEsp
                        real    (dp)                    :: atomEspAvg
                        real    (dp)                    :: atomEspAvg2

                        allocate(accumDistanceMatrix(trj%nAtoms,trj%nAtoms))
                        accumDistanceMatrix = 0

                        atomEspAvg  = 0
                        atomEspAvg2 = 0
                        do step=1,trj%nSteps
                                distanceMatrix      = computeDistanceMatrix(trj%coord(:,:,step))

                                atomEsp             = electroPot(distanceMatrix,trj%atomCharges)

                                atomEspAvg          = atomEspAvg  + atomEsp
                                atomEspAvg2         = atomEspAvg2 + atomEsp**2

                                accumDistanceMatrix = accumDistanceMatrix - distanceMatrix
                        enddo

                        self%atomEspVar =    atomEspAvg2 /real(trj%nSteps,dp)   &
                                          - (atomEspAvg  /real(trj%nSteps,dp))**2

                        if (present(save_DM)) then
                                if (save_DM) then
                                        self%distanceMat = accumDistanceMatrix
                                endif
                        endif

                endsubroutine traj_stats_calculate_espvar

                function traj_map_trajectory(self,map,stride,map_charges,position_weights) &
                                result(mappedTrajectory)

                        use env_kindtypes,      only: si, dp
                        use core_filter,    only: valueCollapse,computeNumUniqueSorted,&
                                                        composeMapping, boundariesToMapping
                        use core_sort,          only: permute, indexQsort

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (traj),       intent(in)          :: self

                        integer (si),       intent(in)          :: map(:)      !!coarse-grained map
                        integer (si),       intent(in),optional :: stride      !! number of steps to skip between
                                                                               !  transferring frames +1.
                        logical,            intent(in),optional :: map_charges !! whether to map charges.
                        real    (dp),       intent(in),optional :: position_weights(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        type    (traj)             :: mappedTrajectory

                        integer (si),allocatable   :: workingMapping(:)
                        logical                    :: map_charges_

                        !local variables
                        integer (si)   :: step

                        !optional argument mirrors.
                        integer (si)   :: stride_

                        integer (si),allocatable   :: permutation_record(:)
                        integer (si),allocatable   :: work_array_int(:)
                        real    (dp),allocatable   :: work_array_real(:)
                        integer (si)               :: numCGsites, spaceDim

                        real    (dp),allocatable   :: reordered_weights(:)

                        if (present(stride)) then
                                if (stride >= 1) then
                                        stride_ = stride
                                else
                                        stride_ = 1
                                endif
                        else
                                stride_ = 1
                        endif

                        if (present(map_charges)) then
                                map_charges_ = map_charges
                        else
                                map_charges_ = .true.
                        endif

                        if (.not.(size(map) == self%nAtoms)) then
                                print*, "ERROR: Attempting to map a trajectory using a mapping of &
                                        &non-understandable length: ", size(map)
                                print*, "Stopping."
                                stop
                        endif

                        workingMapping = map

                        allocate(permutation_record(size(workingMapping)))

                        !We first need the number of sites for memory allocation.
                        ! This requires sorting, and we save the permutation.
                        call indexQsort(workingMapping,permutation_record)
                        call    permute(workingMapping,permutation_record)

                        numCGsites  =  computeNumUniqueSorted(workingMapping)

                        allocate(mappedTrajectory%coord( numCGsites,size(self%coord,2),self%nSteps/stride_))
                        allocate(mappedTrajectory%refAvg(numCGsites,size(self%coord,2)))

                        mappedTrajectory%refAvg = 0

                        allocate(work_array_real(self%nAtoms))
                        allocate(work_array_int(self%nAtoms))

                        if (present(position_weights)) then
                                allocate(reordered_weights,source=position_weights)
                                call permute(reordered_weights,permutation_record)
                        else
                                allocate(reordered_weights(self%nAtoms))
                                reordered_weights =1
                        endif

                        do step = 1, (self%nSteps), stride_ !for every step
                                do spaceDim = 1, (size(self%coord,2)) !for each dimension
                                        !permute the atomistic order to the same as the labels
                                        call permute(positions = self%coord(:,spaceDim,step),&
                                                       mapping = permutation_record,&
                                                        output = work_array_real)

                                        !collapse based on the labels and value.
                                        !we know mapping is already sorted.
                                        mappedTrajectory%coord(:,spaceDim,step) =&
                                                valueCollapse(indexLabels = workingMapping,&
                                                              indexValues = work_array_real,&
                                                             indexWeights = reordered_weights,&
                                                                     sort = .false.,&
                                                             perValueNorm = .true.,&
                                                               constGuard = .true.)
                                enddo
                                mappedTrajectory%refAvg =  mappedTrajectory%refAvg(:,:) &
                                                         + mappedTrajectory%coord(:,:,step)
                        enddo
                        !complete average structure
                        mappedTrajectory%nSteps      = self%nSteps/stride_
                        mappedTrajectory%refAvg      = mappedTrajectory%refAvg/real(mappedTrajectory%nSteps,dp)
                        mappedTrajectory%nDimensions = self%nDimensions

                        mappedTrajectory%nAtoms      = numCGsites

                        call permute(positions=self%atomCharges,&
                                       mapping=permutation_record,&
                                        output=work_array_real)

                        if (map_charges_) then
                                charge_map: block
                                        real    (dp),allocatable :: ones(:)

                                        allocate(ones(size(work_array_real,1)))
                                        ones = 1.0_dp

                                        !Collapse charges per residue.
                                        mappedTrajectory%atomCharges = &
                                                valueCollapse(indexLabels = workingMapping,&
                                                              indexValues = work_array_real,&
                                                                     sort = .false.,&
                                                             indexWeights = ones,&
                                                             perValueNorm = .false.)
                                endblock charge_map
                        endif

                        mappedTrajectory%stride      = stride_

                        !encode parent mapping information
                        mappedTrajectory%is_mapped_trajectory = .true.
                        mappedTrajectory%parent_mapping = map
                        mappedTrajectory%parent_charges = self%atomCharges

                endfunction traj_map_trajectory

                function ltraj_map_trajectory(self,map,stride,map_charges,position_weights) &
                                result(mappedTrajectory)

                        use env_kindtypes,      only: si, dp
                        use core_filter,    only: valueCollapse,computeNumUniqueSorted,&
                                                        composeMapping, boundariesToMapping
                        use core_sort,          only: permute, indexQsort

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (labeledTraj),intent(in)          :: self

                        integer (si),       intent(in)          :: map(:)      !!coarse-grained map
                        integer (si),       intent(in),optional :: stride      !! number of steps to skip between
                                                                               !  transferring frames +1.
                        logical,            intent(in),optional :: map_charges !! whether to map charges.
                        real    (dp),       intent(in),optional :: position_weights(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        type    (traj)             :: mappedTrajectory

                        integer (si),allocatable   :: workingMapping(:)
                        logical                    :: map_charges_

                        !local variables
                        integer (si)   :: step

                        !optional argument mirrors.
                        integer (si)   :: stride_

                        integer (si),allocatable   :: permutation_record(:)
                        integer (si),allocatable   :: work_array_int(:)
                        real    (dp),allocatable   :: work_array_real(:)
                        integer (si)               :: numCGsites, spaceDim

                        real    (dp),allocatable   :: reordered_weights(:)

                        if (present(stride)) then
                                if (stride >= 1) then
                                        stride_ = stride
                                else
                                        stride_ = 1
                                endif
                        else
                                stride_ = 1
                        endif

                        if (present(map_charges)) then
                                map_charges_ = map_charges
                        else
                                map_charges_ = .true.
                        endif

                        !We need to detect if this mapping is defined on the atomistic scale, or the 
                        !residue scale. This is done by size inspection. If on the residue scale, we convert it 
                        !to the atomistic scale.

                        if (size(map) == self%nRes) then
                                !convert from residue to atomistic scale.
                                allocate(workingMapping,source=composeMapping(map,self%site_res))
                        else if (size(map) == self%nAtoms) then
                                allocate(workingMapping,source=map)
                        else
                                print*, "ERROR: Attempting to map a trajectory using a mapping of &
                                        &non-understandable length: ", size(map)
                                print*, "Stopping."
                                stop
                        endif

                        allocate(permutation_record(size(workingMapping)))

                        !We first need the number of sites for memory allocation.
                        ! This requires sorting, and we save the permutation.
                        call indexQsort(workingMapping,permutation_record)
                        call    permute(workingMapping,permutation_record)
                        numCGsites  =  computeNumUniqueSorted(workingMapping)

                        allocate(mappedTrajectory%coord( numCGsites,size(self%coord,2),self%nSteps/stride_))
                        allocate(mappedTrajectory%refAvg(numCGsites,size(self%coord,2)))

                        mappedTrajectory%refAvg = 0

                        allocate(work_array_real(self%nAtoms))
                        allocate(work_array_int(self%nAtoms))

                        if (present(position_weights)) then
                                allocate(reordered_weights,source=position_weights)
                                call permute(reordered_weights,permutation_record)
                        else
                                allocate(reordered_weights(self%nAtoms))
                                reordered_weights =1
                        endif

                        do step = 1, (self%nSteps), stride_ !for every step
                                do spaceDim = 1, (size(self%coord,2)) !for each dimension
                                        !permute the atomistic order to the same as the labels
                                        call permute(positions=self%coord(:,spaceDim,step),&
                                                       mapping=permutation_record,&
                                                        output=work_array_real)

                                        !collapse based on the labels and value.
                                        !we know the mapping is sorted.
                                        mappedTrajectory%coord(:,spaceDim,step) =&
                                                valueCollapse(indexLabels = workingMapping,&
                                                              indexValues = work_array_real,&
                                                             indexWeights = reordered_weights,&
                                                                     sort = .false.,&
                                                             perValueNorm = .true.,&
                                                               constGuard = .true.)
                                enddo
                                mappedTrajectory%refAvg =  mappedTrajectory%refAvg(:,:) &
                                                         + mappedTrajectory%coord(:,:,step)
                        enddo
                        !complete average structure
                        mappedTrajectory%nSteps      = self%nSteps/stride_
                        mappedTrajectory%refAvg      = mappedTrajectory%refAvg/real(mappedTrajectory%nSteps,dp)
                        mappedTrajectory%nDimensions = self%nDimensions

                        mappedTrajectory%nAtoms      = numCGsites

                        call permute(positions=self%atomCharges,&
                                       mapping=permutation_record,&
                                        output=work_array_real)

                        if (map_charges_) then
                                charge_map: block
                                        real    (dp),allocatable :: ones(:)

                                        allocate(ones(size(work_array_real,1)))
                                        ones = 1.0_dp

                                        !Collapse charges per residue.
                                        mappedTrajectory%atomCharges = &
                                                valueCollapse(indexLabels = workingMapping,&
                                                              indexValues = work_array_real,&
                                                                     sort = .false.,&
                                                             indexWeights = ones,&
                                                             perValueNorm = .false.)
                                endblock charge_map
                        endif

                        mappedTrajectory%stride      = stride_

                        !encode parent mapping information
                        mappedTrajectory%is_mapped_trajectory = .true.
                        mappedTrajectory%parent_mapping = map
                        mappedTrajectory%parent_charges = self%atomCharges

                endfunction ltraj_map_trajectory

                function ltraj_ltraj_map_trajectory(self,map,stride,map_charges,position_weights,pass_res)&
                                result(mappedTrajectory)

                        use env_kindtypes,      only: si, dp
                        use core_filter,    only: valueCollapse,computeNumUniqueSorted,&
                                                        composeMapping, boundariesToMapping
                        use core_sort,          only: permute, indexQsort
                        use core_convert,       only: itoa

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (labeledTraj),intent(in)          :: self

                        integer (si),       intent(in)          :: map(:)      !!coarse-grained map
                        integer (si),       intent(in),optional :: stride      !! number of steps to skip between
                                                                               !  transferring frames +1.
                        logical,            intent(in),optional :: map_charges !! whether to map charges.
                        real    (dp),       intent(in),optional :: position_weights(:)
                        logical,            intent(in),optional :: pass_res    !! whether to map charges.
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        type    (labeledTraj)      :: mappedTrajectory

                        integer (si),allocatable   :: workingMapping(:)
                        logical                    :: map_charges_

                        !local variables
                        integer (si)   :: step

                        !optional argument mirrors.
                        integer (si)   :: stride_

                        integer (si),allocatable   :: permutation_record(:)
                        integer (si),allocatable   :: work_array_int(:)
                        real    (dp),allocatable   :: work_array_real(:)
                        integer (si)               :: numCGsites, spaceDim, iter

                        real    (dp),allocatable   :: reordered_weights(:)
                        logical                    :: pass_res_

                        if (present(pass_res)) then
                                pass_res_ = pass_res
                        else
                                pass_res_ = .false.
                        endif

                        if (present(stride)) then
                                if (stride >= 1) then
                                        stride_ = stride
                                else
                                        stride_ = 1
                                endif
                        else
                                stride_ = 1
                        endif

                        if (present(map_charges)) then
                                map_charges_ = map_charges
                        else
                                map_charges_ = .true.
                        endif

                        !We need to detect if this mapping is defined on the atomistic scale, or the 
                        !residue scale. This is done by size inspection. If on the residue scale, we convert it 
                        !to the atomistic scale.

                        if (size(map) == self%nRes) then
                                !convert from residue to atomistic scale.
                                allocate(workingMapping,source=composeMapping(map,self%site_res))
                        else if (size(map) == self%nAtoms) then
                                allocate(workingMapping,source=map)
                        else
                                print*, "ERROR: Attempting to map a trajectory using a mapping of &
                                        &non-understandable length: ", size(map)
                                print*, "Stopping."
                                stop
                        endif

                        allocate(permutation_record(size(workingMapping)))

                        !We first need the number of sites for memory allocation.
                        ! This requires sorting, and we save the permutation.
                        call indexQsort(workingMapping,permutation_record)
                        call    permute(workingMapping,permutation_record)
                        numCGsites  =  computeNumUniqueSorted(workingMapping)

                        allocate(mappedTrajectory%coord( numCGsites,size(self%coord,2),self%nSteps/stride_))
                        allocate(mappedTrajectory%refAvg(numCGsites,size(self%coord,2)))

                        mappedTrajectory%refAvg = 0

                        allocate(work_array_real(self%nAtoms))
                        allocate(work_array_int(self%nAtoms))

                        if (present(position_weights)) then
                                allocate(reordered_weights,source=position_weights)
                                call permute(reordered_weights,permutation_record)
                        else
                                allocate(reordered_weights(self%nAtoms))
                                reordered_weights =1
                        endif


                        do step = 1, (self%nSteps), stride_ !for every step
                                do spaceDim = 1, (size(self%coord,2)) !for each dimension
                                        !permute the atomistic order to the same as the labels
                                        call permute(positions=self%coord(:,spaceDim,step),&
                                                       mapping=permutation_record,&
                                                        output=work_array_real)

                                        !collapse based on the labels and value.
                                        !we know mapping is already sorted.
                                        mappedTrajectory%coord(:,spaceDim,step) =&
                                                valueCollapse(indexLabels = workingMapping,&
                                                              indexValues = work_array_real,&
                                                             indexWeights = reordered_weights,&
                                                                     sort = .false.,&
                                                             perValueNorm = .true.,&
                                                               constGuard = .true.)
                                enddo
                                mappedTrajectory%refAvg =  mappedTrajectory%refAvg(:,:) &
                                                         + mappedTrajectory%coord(:,:,step)
                        enddo
                        !complete average structure
                        mappedTrajectory%nSteps      = self%nSteps/stride_
                        mappedTrajectory%refAvg      = mappedTrajectory%refAvg/real(mappedTrajectory%nSteps,dp)
                        mappedTrajectory%nDimensions = self%nDimensions

                        mappedTrajectory%nAtoms      = numCGsites

                        !this is quite primitive in terms of assigning residue labels--
                        !the problem is inherently not valid generally. We warn, but we purposely 
                        !aren't doing inteligent matching here.
                        if (pass_res_) then
                                mappedTrajectory%nRes        = self%nRes
                                !Now all residues are 1 bead each.
                                mappedTrajectory%residue_start_site = (/ (iter, iter = 1, self%nRes) /)
                                mappedTrajectory%site_res           = (/ (iter, iter = 1, self%nRes) /)
                        else
                                mappedTrajectory%nRes        = -1
                        endif

                        call permute(positions=self%atomCharges,&
                                       mapping=permutation_record,&
                                        output=work_array_real)

                        if (map_charges_) then
                                charge_map: block
                                        real    (dp),allocatable :: ones(:)

                                        allocate(ones(size(work_array_real,1)))
                                        ones = 1.0_dp

                                        !Collapse charges per residue.
                                        mappedTrajectory%atomCharges = &
                                                valueCollapse(indexLabels = workingMapping,&
                                                              indexValues = work_array_real,&
                                                                     sort = .false.,&
                                                             indexWeights = ones,&
                                                             perValueNorm = .false.)
                                endblock charge_map
                        endif

                        !map masses.

                        mass_map: block
                                real    (dp),allocatable :: ones(:)

                                allocate(ones(size(work_array_real,1)))
                                ones = 1.0_dp

                                call permute(positions=self%site_mass,&
                                               mapping=permutation_record,&
                                                output=work_array_real)

                                !Collapse charges per residue.
                                mappedTrajectory%site_mass = &
                                        valueCollapse(indexLabels = workingMapping,&
                                                      indexValues = work_array_real,&
                                                             sort = .false.,&
                                                     indexWeights = ones,&
                                                     perValueNorm = .false.)
                        endblock mass_map

                        site_names: block
                                integer                 :: iter
                                character(1),parameter  :: res_prefix = 'C'

                                allocate(mappedTrajectory%site_labels(numCGsites))

                                do iter=1,size(mappedTrajectory%site_labels)
                                        mappedTrajectory%site_labels(iter) = res_prefix//trim(itoa(iter))
                                enddo
                        endblock site_names

                        mappedTrajectory%stride      = stride_

                        !encode parent mapping information
                        mappedTrajectory%is_mapped_trajectory = .true.
                        mappedTrajectory%parent_mapping = map
                        mappedTrajectory%parent_charges = self%atomCharges

                endfunction ltraj_ltraj_map_trajectory

                pure function ltraj_ltraj_map_trajectory_CA(self,stride,map_charges,pass_res) &
                                result(mappedTrajectory)

                        use env_kindtypes,      only: si, dp
                        use core_filter,    only: valueCollapse
                        use core_convert,    only: itoa

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (labeledTraj),intent(in)          :: self

                        integer (si),       intent(in),optional :: stride !number of steps to skip between
                                                                          !transferring frames +1.
                        logical,            intent(in),optional :: map_charges !Whether to collapse charges.
                        logical,            intent(in),optional :: pass_res !Whether to collapse charges.
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        type    (labeledTraj)                   :: mappedTrajectory

                        !local variables
                        integer (si)   :: step
                        integer (si)   :: atom1, i

                        !real    (dp)   :: caAvg(3)
                        integer (si)   :: caIndex

                        !optional argument mirrors.
                        integer (si)   :: stride_
                        logical        :: map_charges_, pass_res_

                        if (present(pass_res)) then
                                pass_res_ = pass_res
                        else
                                pass_res_ = .true.
                        endif

                        if (present(stride)) then
                                if (stride >= 1) then
                                        stride_ = stride
                                else
                                        stride_ = 1
                                endif
                        else
                                stride_ = 1
                        endif

                        if (present(map_charges)) then
                                map_charges_  = map_charges
                        else
                                map_charges_  = .true.
                        endif

                        !If we're using carbon alpha based filtering
                        allocate( mappedTrajectory%coord( self%nRes,3, self%nSteps/stride_) )
                        allocate( mappedTrajectory%refAvg(self%nRes,3) )

                        mappedTrajectory%refAvg = 0

                        do step = 1, (self%nSteps), stride_
                                !pull out C alphas from the trajectory and put them into
                                !our analysis trajectory
                                caIndex = 1
                                do atom1 = 1, self%nAtoms
                                        if (self%site_labels(atom1) .eq. "CA") then
                                                mappedTrajectory%coord(caIndex,:,step) = &
                                                                self%coord(atom1,:,step)
                                                caIndex = caIndex + 1
                                        endif
                                enddo

                                !Provide an average configuration. Note normalization after loop.
                                mappedTrajectory%refAvg =  mappedTrajectory%refAvg(:,:) &
                                                         + mappedTrajectory%coord(:,:,step)
                        enddo
                        mappedTrajectory%nSteps      = self%nSteps / stride_
                        !complete average structure
                        mappedTrajectory%refAvg      = mappedTrajectory%refAvg/real(mappedTrajectory%nSteps,dp)

                        mappedTrajectory%nDimensions = self%nDimensions

                        !This is a one site per residue mapping.
                        mappedTrajectory%nAtoms      = self%nRes

                        if (pass_res_) then
                                mappedTrajectory%nRes        = self%nRes
                                !Now all residues are 1 bead each.
                                mappedTrajectory%residue_start_site = (/ (i, i = 1, self%nRes) /)
                                mappedTrajectory%site_res           = (/ (i, i = 1, self%nRes) /)
                        else
                                mappedTrajectory%nRes = -1
                        endif

                        !Copy over stride
                        mappedTrajectory%stride      = stride_

                        if (map_charges_) then
                                charge_map: block
                                        real    (dp),allocatable :: ones(:)

                                        allocate(ones(size(self%site_res,1)))
                                        ones = 1.0_dp

                                        !Collapse charges per residue.
                                        mappedTrajectory%atomCharges = &
                                                valueCollapse(indexLabels = self%site_res,&
                                                              indexValues = self%atomCharges,&
                                                                     sort = .true.,&
                                                             indexWeights = ones,&
                                                             perValueNorm = .false.)
                                endblock charge_map
                        endif

                        !map masses.
                        mass_map: block
                                real    (dp),allocatable :: ones(:)

                                allocate(ones(size(self%site_res,1)))
                                ones = 1.0_dp

                                !Collapse charges per residue.
                                mappedTrajectory%site_mass = &
                                        valueCollapse(indexLabels = self%site_res,&
                                                      indexValues = self%site_mass,&
                                                             sort = .true.,&
                                                     indexWeights = ones,&
                                                     perValueNorm = .false.)
                        endblock mass_map

                        site_names: block
                                integer                 :: iter
                                character(1),parameter  :: res_prefix = 'C'

                                allocate(mappedTrajectory%site_labels(self%nRes))

                                do iter=1,size(mappedTrajectory%site_labels)
                                        mappedTrajectory%site_labels(iter) = res_prefix//trim(itoa(iter))
                                enddo
                        endblock site_names

                        !encode parent mapping information
                        mappedTrajectory%is_mapped_trajectory = .true.
                        mappedTrajectory%parent_mapping = self%site_res
                        mappedTrajectory%parent_charges = self%atomCharges

                endfunction ltraj_ltraj_map_trajectory_CA

                subroutine write_xyz_trajectory(self,filename,comment)
                        !Writes an xyz trajectory (VMD format).

                        use IO_xyz,             only: write_xyz_coord

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (traj),     intent(in)          :: self
                        character (len=*),intent(in)          :: filename
                        character (len=*),intent(in),optional :: comment
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variable

                        if (present(comment)) then
                                call write_xyz_coord(filename,self%coord)
                        else
                                call write_xyz_coord(filename,self%coord,comment)
                        endif

                        return

                endsubroutine write_xyz_trajectory

                subroutine write_xyz_average(self,filename,comment)
                        !Writes an xyz trajectory (VMD format).

                        use IO_xyz,             only: write_xyz_coord

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (traj),     intent(in)          :: self
                        character (len=*),intent(in)          :: filename
                        character (len=*),intent(in),optional :: comment
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        if (allocated(self%refAvg)) then
                                if (present(comment)) then
                                        call write_xyz_coord(filename,self%refAvg)
                                else
                                        call write_xyz_coord(filename,self%refAvg,comment)
                                endif
                        else
                                print*, "WARNING: Could not write average structure as array is not allocated."
                        endif

                endsubroutine write_xyz_average

                subroutine write_csv_trajectory(self,filename)

                        use IO_csv,             only: write_csv

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (traj),     intent(in)          :: self
                        character (len=*),intent(in)          :: filename
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local variable

                        call write_csv(filename,self%coord)

                endsubroutine write_csv_trajectory

                pure subroutine traj_bootstrapCopy(self,copy,indices)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (traj),     intent(in   )          :: self
                        class (traj),     intent(inout)          :: copy
                        integer (si),     intent(in   )          :: indices(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        copy%nAtoms = self%nAtoms
                        copy%nSteps = size(indices)
                        copy%stride = 1

                        !if the indices are valid and self has coords, copy.
                        if (allocated(self%coord) .and. &
                            maxval(indices) <= size(self%coord,3)) then
                                copy%coord = self%coord(:,:,indices)
                        !else, destroy coords in copy to indicate failure.
                        elseif (allocated(copy%coord)) then
                                deallocate(copy%coord)
                        endif

                        !average is no loner valid.
                        if (allocated(copy%refAvg)) then
                                deallocate(copy%refAvg)
                        endif

                        if (allocated(self%atomCharges)) then
                                copy%atomCharges = self%atomCharges
                        endif

                endsubroutine traj_bootstrapCopy

                pure subroutine ltraj_bootstrapCopy(self,copy,indices)

                        implicit none

                        !!!!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class (labeledtraj),     intent(in   )  :: self
                        class (traj),            intent(inout)  :: copy
                        integer (si),            intent(in   )  :: indices(:)
                        !!!!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        select type (copy)
                        class is (traj)

                                copy%nAtoms = self%nAtoms
                                copy%nSteps = size(indices)
                                copy%stride = 1

                                if (allocated(self%coord) .and. &
                                    maxval(indices) <= size(self%coord,3)) then
                                        copy%coord = self%coord(:,:,indices)
                                elseif (allocated(copy%coord)) then
                                        deallocate(copy%coord)
                                endif

                                !average is no loner valid.
                                if (allocated(copy%refAvg)) then
                                        deallocate(copy%refAvg)
                                endif

                                if (allocated(self%atomCharges)) then
                                        copy%atomCharges = self%atomCharges
                                endif

                        class is (labeledTraj)

                                !careful deallocate
                                if (allocated(copy%refAvg)) then
                                        deallocate(copy%refAvg)
                                endif

                                copy%nAtoms = self%nAtoms
                                copy%nSteps = size(indices)
                                copy%stride = 1

                                if (allocated(self%coord) .and. &
                                    maxval(indices) <= size(self%coord,3)) then
                                        copy%coord = self%coord(:,:,indices)
                                elseif (allocated(copy%coord)) then
                                        deallocate(copy%coord)
                                endif

                                if (allocated(self%atomCharges)) then
                                        copy%atomCharges = self%atomCharges
                                endif

                                if (allocated(self%site_labels)) then
                                        copy%site_labels = self%site_labels
                                endif 

                                if (allocated(self%site_type)) then
                                        copy%site_type = self%site_type
                                endif 

                                if (allocated(self%site_res)) then
                                        copy%site_res = self%site_res
                                endif 

                                if (allocated(self%site_mass)) then
                                        copy%site_mass = self%site_mass
                                endif 

                                if (allocated(self%residue_start_site)) then
                                        copy%residue_start_site = self%residue_start_site
                                endif 

                                return
                        end select

                endsubroutine ltraj_bootstrapCopy

                subroutine traj_stats_calculate_avg_dmat(self,trj,variance)

                        use routines_math,      only: computeDistanceMatrix
                        use obj_accum,          only: accum, varAccum

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj_stats),intent(inout)          :: self
                        class   (traj),      intent(in   )          :: trj
                        logical,             intent(in   ),optional :: variance
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical                    :: variance_
                        class(accum),allocatable   :: dist_accum(:,:)
                        real    (dp),allocatable   :: dist_mat(:,:)
                        integer (si)               :: frame_iter, var_iter_1, var_iter_2

                        if (present(variance)) then
                                variance_ = variance
                        else
                                variance_ = .false.
                        endif

                        if (variance_) then
                                allocate(varAccum :: dist_accum(trj%nAtoms,trj%nAtoms))
                        else
                                allocate(   accum :: dist_accum(trj%nAtoms,trj%nAtoms))
                        endif

                        call dist_accum%reset()

                        do frame_iter=1,size(trj%coord,3)
                                dist_mat = computeDistanceMatrix(trj%coord(:,:,frame_iter))
                                !$omp parallel do private(var_iter_1, var_iter_2)
                                do var_iter_1=1,trj%nAtoms
                                do var_iter_2=1,trj%nAtoms
                                        call dist_accum(var_iter_1,var_iter_2)%&
                                                add(dist_mat(var_iter_1,var_iter_2))
                                enddo
                                enddo
                                !$omp end parallel do
                        enddo

                        self%refAvg_distance_mat = dist_accum%getMean()

                        !we need type select because the compiler isn't smart enough to know the subtype.
                        if (variance_) then
                                select type (dist_accum)
                                type is (varAccum)
                                        self%refAvg_distance_var_mat = dist_accum%getVar()
                                end select
                        endif

                endsubroutine traj_stats_calculate_avg_dmat

                pure subroutine traj_stats_calculate_avg_struct_diameter(self,trj)

                        use routines_math,      only: diameter

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj_stats),intent(inout)          :: self
                        class   (traj),      intent(in   )          :: trj
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        if (allocated(trj%refAvg)) then
                                self%refAvg_diameter = diameter(trj%refAvg)
                                self%is_calculated_refAvg_diameter = .true.
                        else
                                self%is_calculated_refAvg_diameter = .false.
                        endif

                endsubroutine traj_stats_calculate_avg_struct_diameter

                pure function traj_stats_get_avg_struct_diameter(self)

                        use env_kindtypes,      only: dp 

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj_stats),intent(in   )          :: self
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp)    :: traj_stats_get_avg_struct_diameter

                        if (self%is_calculated_refAvg_diameter) then
                                traj_stats_get_avg_struct_diameter = self%refAvg_diameter
                        else
                                traj_stats_get_avg_struct_diameter = 0
                        endif

                endfunction traj_stats_get_avg_struct_diameter

                pure function traj_stats_get_avg_struct_dmat(self)

                        use env_kindtypes,      only: dp 

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj_stats),intent(in   )          :: self
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable    :: traj_stats_get_avg_struct_dmat(:,:)

                        if (allocated(self%refAvg_distance_mat)) then
                                traj_stats_get_avg_struct_dmat = self%refAvg_distance_mat
                        else
                                allocate(traj_stats_get_avg_struct_dmat(0,0))
                        endif

                endfunction traj_stats_get_avg_struct_dmat

                pure function traj_stats_get_avg_struct_var_dmat(self)

                        use env_kindtypes,      only: dp 

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj_stats),intent(in   )          :: self
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable    :: traj_stats_get_avg_struct_var_dmat(:,:)

                        if (allocated(self%refAvg_distance_var_mat)) then
                                traj_stats_get_avg_struct_var_dmat = self%refAvg_distance_var_mat
                        else
                                allocate(traj_stats_get_avg_struct_var_dmat(0,0))
                        endif

                endfunction traj_stats_get_avg_struct_var_dmat

                subroutine traj_stats_gen_statistics(self,trj,design)

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(traj_stats),        intent(inout)   :: self
                        class(traj),              intent(in   )   :: trj
                        type (trajPreproc_config),intent(in   )   :: design
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        ! compute ESP var
                        call self%computeESPvar(trj=trj,save_DM=.true.)

                        call self%computeavgsdistancemat(trj=trj,&
                                                         variance=design%pairwise_variance_distance_flag)

                        call self%computeavgsdiameter(trj=trj)

                        ! compute the covariance matrix
                        call self%computeedmat(trj=trj,filterdof=design%proj_num_dof)

                        ! compute kernel dist mat
                        call self%computecustomdistmat(distance_type=design%gedcg_diss_mat_type)

                        ! compute chargediff mat
                        call self%computePowChargeDiffMat(trj=trj,power=2)

                        ! Print reference structure
                        call trj%writeAverageXYZfile(filename='average_struct.xyz')

                endsubroutine traj_stats_gen_statistics

                subroutine traj_training_data_preproc(self,design)

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj_training_data),intent(inout)          :: self
                        class   (trajPreproc_config),intent(in   )          :: design
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        self%preproc_design = design

                        if (self%preproc_design%align_flag) then
                                call self%trj%iterativealign(center=.true.)
                        endif

                        call traj_stats_gen_statistics(self%stats,self%trj,self%preproc_design)

                endsubroutine traj_training_data_preproc

                function traj_omp_gen_dipole_traj(self,section_mapping,origin,site_weights)&
                               result(dipole_traj)

                        use core_multipole,     only: dipole
                        use core_stat,          only: which

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj), intent(in   )           :: self
                        integer (si),   intent(in   )           :: section_mapping(:)
                        real    (dp),   intent(in   ),optional  :: origin(:)
                        real    (dp),   intent(in   ),optional  :: site_weights(:)
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        type    (traj)                  :: dipole_traj
                        integer (si), allocatable       :: working_mapping(:)
                        integer (si), allocatable       :: selection(:)
                        integer (si)                    :: orig_min_map_value, num_sites
                        integer (si)                    :: site_index, step_index

                        orig_min_map_value = minval(section_mapping)

                        if (orig_min_map_value /= 1) then
                                working_mapping = section_mapping - orig_min_map_value + 1
                        else
                                working_mapping = section_mapping
                        endif

                        num_sites = maxval(working_mapping)

                        dipole_traj%nAtoms     = num_sites
                        dipole_traj%nDimensions= 3          !there are 3 values in each dipole.
                        dipole_traj%nSteps     = self%nSteps
                        dipole_traj%stride     = 1          !we force full stride right now.

                        allocate(dipole_traj%coord(dipole_traj%nAtoms,&
                                            dipole_traj%nDimensions,&
                                            dipole_traj%nSteps))

                        if (present(origin)) then
                                do step_index=1,dipole_traj%nSteps
                                do site_index=1,num_sites
                                        selection = which(section_mapping == site_index)
                                        dipole_traj%coord(site_index,:,step_index) = &
                                                dipole(points = self%coord(selection,:,step_index),&
                                                       values = self%atomCharges(selection),&
                                                       origin = origin)
                                enddo
                                enddo
                        elseif (present(site_weights)) then
                                do step_index=1,dipole_traj%nSteps
                                do site_index=1,num_sites
                                        selection = which(section_mapping == site_index)
                                        dipole_traj%coord(site_index,:,step_index) = &
                                                dipole(points = self%coord(selection,:,step_index),&
                                                       values = self%atomCharges(selection),&
                                               origin_weights = site_weights(selection))
                                enddo
                                enddo
                        else
                                do step_index=1,dipole_traj%nSteps
                                do site_index=1,num_sites
                                        selection = which(section_mapping == site_index)
                                        dipole_traj%coord(site_index,:,step_index) = &
                                                dipole(points = self%coord(selection,:,step_index),&
                                                       values = self%atomCharges(selection))
                                enddo
                                enddo
                        endif

                endfunction traj_omp_gen_dipole_traj

                function traj_omp_gen_quadrupole_traj(self,section_mapping,origin, site_weights) result(quadrupole_traj)

                        use core_multipole,     only: quadrupole
                        use core_stat,          only: which

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj), intent(in   )          :: self
                        integer (si),   intent(in   )          :: section_mapping(:)
                        real    (dp),   intent(in   ),optional :: origin(:)
                        real    (dp),   intent(in   ),optional :: site_weights(:)
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        type    (traj)                  :: quadrupole_traj
                        integer (si), allocatable       :: working_mapping(:)
                        real    (dp), allocatable       :: quad(:,:)
                        integer (si), allocatable       :: selection(:)
                        integer (si)                    :: orig_min_map_value, num_sites
                        integer (si)                    :: site_index, step_index

                        orig_min_map_value = minval(section_mapping)

                        if (orig_min_map_value /= 1) then
                                working_mapping = section_mapping - orig_min_map_value + 1
                        else
                                working_mapping = section_mapping
                        endif

                        num_sites = maxval(working_mapping)

                        quadrupole_traj%nAtoms     = num_sites
                        quadrupole_traj%nDimensions= 9          !there are 3 values in each quadrupole.
                        quadrupole_traj%nSteps     = self%nSteps
                        quadrupole_traj%stride     = 1          !we force full stride right now.

                        allocate(quadrupole_traj%coord(quadrupole_traj%nAtoms,&
                                            quadrupole_traj%nDimensions,&
                                            quadrupole_traj%nSteps))

                        if (present(origin)) then
                                do step_index=1,quadrupole_traj%nSteps
                                do site_index=1,num_sites
                                        selection = which(section_mapping == site_index)
                                        quad      = quadrupole(points = self%coord(selection,:,step_index),&
                                                               values = self%atomCharges(selection),&
                                                               origin = origin)
                                        quadrupole_traj%coord(site_index,:,step_index) = &
                                                                        reshape(quad,[size(quad)])
                                enddo
                                enddo
                        elseif (present(site_weights)) then
                                do step_index=1,quadrupole_traj%nSteps
                                do site_index=1,num_sites
                                        selection = which(section_mapping == site_index)
                                        quad      = quadrupole(points = self%coord(selection,:,step_index),&
                                                               values = self%atomCharges(selection),&
                                                       origin_weights = site_weights(selection))

                                        quadrupole_traj%coord(site_index,:,step_index) = &
                                                                        reshape(quad,[size(quad)])
                                enddo
                                enddo
                        else
                                do step_index=1,quadrupole_traj%nSteps
                                do site_index=1,num_sites
                                        selection = which(section_mapping == site_index)
                                        quad      = quadrupole(points = self%coord(selection,:,step_index),&
                                                               values = self%atomCharges(selection))
                                        quadrupole_traj%coord(site_index,:,step_index) = &
                                                                        reshape(quad,[size(quad)])
                                enddo
                                enddo
                        endif 

                endfunction traj_omp_gen_quadrupole_traj

                !generates an array of radius of gyration values for subsets of atoms in a trajectory.
                !assums the section_mapping is contiguous.
                function traj_gen_rog_traj(self,section_mapping) result(rog_traj)

                        use core_stat,     only: rog, which

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class   (traj), intent(in   )   :: self
                        integer (si),   intent(in   )   :: section_mapping(:)
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        type    (traj)                  :: rog_traj

                        integer (si), allocatable       :: working_mapping(:)
                        integer (si), allocatable       :: selection(:)
                        integer (si)                    :: orig_min_map_value, num_sites
                        integer (si)                    :: step_index, site_index

                        orig_min_map_value = minval(section_mapping)

                        if (orig_min_map_value /= 1) then
                                working_mapping = section_mapping - orig_min_map_value + 1
                        else
                                working_mapping = section_mapping
                        endif

                        num_sites = maxval(working_mapping)

                        rog_traj%nAtoms     = num_sites
                        rog_traj%nDimensions= 1          !there are 3 values in each rog.
                        rog_traj%nSteps     = self%nSteps
                        rog_traj%stride     = 1          !we force full stride right now.

                        allocate(rog_traj%coord(num_sites,1,rog_traj%nSteps))

                        do step_index=1,rog_traj%nSteps
                        do site_index=1,num_sites
                                selection = which(section_mapping == site_index)
                                rog_traj%coord(site_index,:,step_index) = &
                                        rog(self%coord(selection,:,step_index))
                        enddo
                        enddo

                endfunction traj_gen_rog_traj

endmodule obj_trajectory
