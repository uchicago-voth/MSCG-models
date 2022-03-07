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
! Provides routines for reading command line options.
! Basic argument checking is done here, as static typing makes multitype argument
! handeling later difficult.

module obj_commandline

        use env_kindtypes,    only: si, dp
        use abs_obj_model,    only: model_config
        use obj_trajectory,   only: trajPreproc_config
        use routines_multipoleAnalysis,&
                              only: multipoleAnalysis_config

        implicit none

        private

        public runParameters, get_option_value_string, sget_option_value_static_string_array

        type                    :: runParameters

                private

                !Flags to hold whether options were changed from the default.

                logical :: model_log_flag    = .false.
                logical :: bs_log_flag       = .false.
                logical :: resampling_summary_log_flag = .false.
                logical :: model_summary_log_flag = .false.

                logical :: delta_step_flag   = .false.
                logical :: atom_dcd_flag     = .false.
                logical :: n_cg_flag         = .false.

                logical :: traj_preproc_map_type_flag = .false.

                logical :: atom_psf_flag     = .false.
                logical :: ed_dof_flag       = .false.
                logical :: residual_weights_flag &
                                             = .false.
                logical :: nsets_flag        = .false.
                logical :: modelType_flag    = .false.
                logical :: sa_accept_flag    = .false.
                logical :: sa_descent_control_flag&
                                             = .false.
                logical :: max_sa_steps_flag = .false.
                logical :: max_sd_steps_flag = .false.

                logical :: linearModel_flag  = .false.
                logical :: spatialModel_flag =.false.
                logical :: spectralModel_flag=.false.

                logical :: n_res_bs_reps_Flag    = .false.
                logical :: n_res_bs_models_Flag  = .false.

                logical :: n_bs_reps_Flag = .false.
                logical :: n_mc_norm_reps_Flag = .false.
                logical :: use_bca_ci_flag = .false.

                logical :: model_sort_by_flag = .false.

                logical :: runConfig      = .false.

                logical :: cluster_density_tol_flag = .false.
                logical :: kmeans_max_iter_flag = .false.
                logical :: kmeans_max_try_flag  = .false.
                logical :: spectral_regularization_type_flag = .false.

                logical :: spectral_use_edcg_metric_flag = .false.
                logical :: spectral_edcg_dist_conv_type_flag = .false.
                logical :: spectral_edcg_dist_conv_param_flag= .false.
                logical :: spectral_edcg_kNN_flag      = .false.
                logical :: spectral_edcg_kNN_type_flag = .false.

                logical :: spectral_use_spatial_metric_flag = .false.
                logical :: spectral_spatial_dist_conv_type_flag = .false.
                logical :: spectral_spatial_dist_conv_param_flag= .false.
                logical :: spectral_spatial_kNN_flag      = .false.
                logical :: spectral_spatial_kNN_type_flag = .false.

                logical :: spectral_use_pairVar_metric_flag = .false.
                logical :: spectral_pairVar_dist_conv_type_flag = .false.
                logical :: spectral_pairVar_dist_conv_param_flag= .false.
                logical :: spectral_pairVar_kNN_flag      = .false.
                logical :: spectral_pairVar_kNN_type_flag = .false.

                logical :: spectral_use_charge_metric_flag = .false.
                logical :: spectral_charge_dist_conv_type_flag = .false.
                logical :: spectral_charge_dist_conv_param_flag= .false.
                logical :: spectral_charge_kNN_flag      = .false.
                logical :: spectral_charge_kNN_type_flag = .false.

                logical :: mapping_charge_split_flag  = .false.

                logical :: multipole_analysis_flag = .false.

                logical :: map_accumulator_type_flag = .false.

                !Command line options.
                character (:), allocatable,public :: atomPsfFile, atomDcdFile, modelLog, bsLog, modelType
                character (:), allocatable,public :: resampling_summary_log, model_summary_log
                integer   (si),            public :: deltaStep, n_res_bs_reps, n_res_bs_models, n_bs_reps
                integer   (si),            public :: n_mc_norm_reps, nSets
                logical,                   public :: use_bca_ci
                character (:), allocatable,public :: model_sort_by
                character (:), allocatable,public :: traj_preproc_map_type

                type(model_config),             public :: model_design
                type(trajPreproc_config),       public :: trajPreproc_design
                type(multipoleAnalysis_config), public :: multipoleAnalysis_design

                contains
                        procedure       :: parse_command_line   => parse_command_line
                        procedure       :: check_values         => check_set_values
        end type runParameters


        contains

                !trivial function to determine if a given string is a option.
                elemental function is_option(to_check)

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character(*), intent(in   )     :: to_check
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical         :: is_option

                        is_option = .false.

                        if (len_trim(to_check) == 0) then
                                return
                        endif

                        if (to_check(1:1) == '-') then
                                is_option = .true.
                        endif

                endfunction is_option

                !extracts the maximum length of nonzero characters in array of strings
                pure function max_string_length(string_array,subset_max_index)

                        use env_kindtypes,    only: si

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character(*),   intent(in   )          :: string_array(:)
                        integer,        intent(in   ),optional :: subset_max_index
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si) :: max_string_length

                        integer (si) :: local_length
                        integer (si) :: iter

                        max_string_length = 0

                        if (present(subset_max_index)) then
                                do iter=1,min(size(string_array),subset_max_index)
                                        local_length = len_trim(string_array(iter))
                                        if (local_length > max_string_length) then
                                                max_string_length = local_length
                                        endif
                                enddo
                        else
                                do iter=1,size(string_array)
                                        local_length = len_trim(string_array(iter))
                                        if (local_length > max_string_length) then
                                                max_string_length = local_length
                                        endif
                                enddo
                        endif

                endfunction max_string_length

                function get_option_value_string(arg_index,increment_index) result(string_option)

                        use env_kindtypes,    only: si

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),   intent(inout)          :: arg_index
                        logical,        intent(in   ),optional :: increment_index
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        character(:),allocatable               :: string_option  ! Holds read argument.

                        !return value
                        character(255)                         :: buffer! Holds read argument.

                        logical         :: increment_index_

                        !optional arg guards
                        if (present(increment_index)) then
                                increment_index_ = increment_index
                        else
                                increment_index_ = .false.
                        endif

                        call get_command_argument(arg_index,buffer)

                        string_option = trim(buffer)

                        !the argument is already a string, so no more processing required.

                        if (increment_index_) then
                                arg_index = arg_index + 1
                        endif

                endfunction get_option_value_string

                function get_option_value_integer_si(arg_index,increment_index) result(integer_option)

                        use env_kindtypes,    only: si
                        use core_convert,     only: atoi, is_integer

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),   intent(inout)          :: arg_index
                        logical,        intent(in   ),optional :: increment_index
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        character (LEN=128)      :: arg_value  ! Holds read argument.

                        character (LEN=128)     :: parent_arg_value  ! Holds read argument.

                        integer (si)            :: integer_option

                        logical                 :: increment_index_

                        !optional arg guards
                        if (present(increment_index)) then
                                increment_index_ = increment_index
                        else
                                increment_index_ = .false.
                        endif

                        call get_command_argument(arg_index,arg_value)

                        if (is_integer(arg_value)) then
                                integer_option = atoi(arg_value)
                        else
                                call get_command_argument(arg_index - 1 ,parent_arg_value)
                                print*, 'ERROR: Bad argument ', trim(parent_arg_value)
                                print*, 'ERROR: Expected integer type argument. Got ', arg_value
                                print*, 'ERROR: Erroroneous argument type. Stopping.'
                                stop
                        endif

                        if (increment_index_) then
                                arg_index = arg_index + 1
                        endif

                endfunction get_option_value_integer_si

                function get_option_value_real_dp(arg_index,increment_index) result(real_option)

                        use env_kindtypes,    only: si
                        use core_convert,     only: atof, is_real

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),   intent(inout)          :: arg_index
                        logical,        intent(in   ),optional :: increment_index
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        character (LEN=128)     :: arg_value  ! Holds read argument.

                        character (LEN=128)     :: parent_arg_value  ! Holds read argument.

                        real    (dp)            :: real_option

                        logical                 :: increment_index_

                        !optional arg guards
                        if (present(increment_index)) then
                                increment_index_ = increment_index
                        else
                                increment_index_ = .false.
                        endif

                        call get_command_argument(arg_index,arg_value)

                        if (is_real(arg_value)) then
                                real_option = atof(arg_value)
                        else
                                call get_command_argument(arg_index - 1 ,parent_arg_value)
                                print*, 'ERROR: Bad argument ', trim(parent_arg_value)
                                print*, 'ERROR: Expected real number type argument. Got ', arg_value
                                print*, 'ERROR: Erroroneous argument type. Stopping.'
                                stop
                        endif

                        if (increment_index_) then
                                arg_index = arg_index + 1
                        endif

                endfunction get_option_value_real_dp

                !this has to be a subroutine, not a function, to stop ICEs in gfortran.
                !Specific error unknown.
                subroutine sget_option_value_static_string_array(string_array_option,arg_index,increment_index)

                        use env_kindtypes,    only: si

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character(*),allocatable,intent(  out)          :: string_array_option(:)
                        integer (si),            intent(inout)          :: arg_index
                        logical,                 intent(in   ),optional :: increment_index
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value

                        !gives the lengths of buffers for filenames.
                        integer (si),parameter   :: string_option_len = 4096 !this seems to be the only
                                                                                  !safe length.

                        !gives the buffer parameter for the number of accepted strings.
                        integer (si),parameter   :: string_array_option_buffer_size = 64
                        character(string_option_len),&
                                     allocatable :: string_array_option_buffer(:)
                        integer (si)             :: buffer_write_mark
                        integer (si)             :: internal_arg_index
                        character (LEN=4096)     :: arg_value         ! Holds individual argument.
                        character (LEN=4096)     :: parent_arg_value  ! Holds individual argument.

                        logical                  :: increment_index_

                        !optional arg guards
                        if (present(increment_index)) then
                                increment_index_ = increment_index
                        else
                                increment_index_ = .false.
                        endif

                        allocate(string_array_option_buffer(string_array_option_buffer_size))

                        !this is move around reversibly inside the routine.
                        internal_arg_index = arg_index

                        buffer_write_mark = 1
                        outer: do while (.true.)
                                call get_command_argument(internal_arg_index,arg_value)

                                if (.not. is_option(arg_value) .and. &
                                    (command_argument_count() >= internal_arg_index)) then
                                        if (buffer_write_mark <= string_array_option_buffer_size) then
                                                string_array_option_buffer(buffer_write_mark) = arg_value
                                                buffer_write_mark = buffer_write_mark + 1

                                                internal_arg_index = internal_arg_index + 1
                                        else
                                                print*, 'ERROR: Reading a string array command line &
                                                        &parameter overan the given buffer of size ', &
                                                        string_array_option_buffer_size
                                                print*, 'ERROR: Stopping.'
                                                stop
                                        endif
                                else
                                        exit outer
                                endif
                        enddo outer

                        if (buffer_write_mark == 1) then
                                call get_command_argument(arg_index - 1 ,parent_arg_value)
                                print*, 'ERROR: Absent or bad arguments for ', trim(parent_arg_value)
                                print*, 'ERROR: Expected to populate an string array from command line arguments, &
                                        &but found no valid non-option strings.'
                                print*, 'ERROR: Stopping.'
                                stop
                        else
                                !all is valid, copy our trimmed results into the return value.
                                return_transfer: block
                                        integer (si) :: required_len_buffer_size
                                        integer (si) :: iter

                                        required_len_buffer_size = max_string_length(string_array_option_buffer,buffer_write_mark-1)

                                        if (allocated(string_array_option)) deallocate(string_array_option)
                                        allocate(string_array_option(buffer_write_mark-1))

                                        do iter=1,size(string_array_option)
                                                string_array_option(iter) = string_array_option_buffer(iter)(1:required_len_buffer_size)
                                        enddo
                                endblock return_transfer
                        endif

                        if (increment_index_) then
                                arg_index = internal_arg_index
                        endif

                endsubroutine sget_option_value_static_string_array

                function get_option_value_real_array_dp(arg_index,increment_index) result(real_option)

                        use env_kindtypes,    only: si, dp
                        use core_convert,     only: atof, is_real

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),   intent(inout)          :: arg_index
                        logical,        intent(in   ),optional :: increment_index
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp),allocatable :: real_option(:)

                        integer (si),parameter   :: real_option_buffer_size = 30
                        real    (dp)             :: real_option_buffer(real_option_buffer_size)
                        integer (si)             :: buffer_write_mark

                        integer (si)             :: internal_arg_index
                        character (LEN=128)      :: arg_value  ! Holds individual argument.
                        character (LEN=128)      :: parent_arg_value  ! Holds individual argument.

                        logical                  :: increment_index_

                        !optional arg guards
                        if (present(increment_index)) then
                                increment_index_ = increment_index
                        else
                                increment_index_ = .false.
                        endif

                        internal_arg_index = arg_index

                        buffer_write_mark = 1

                        outer: do while (.true.)
                                call get_command_argument(internal_arg_index,arg_value)

                                if (is_real(arg_value)) then
                                        if (buffer_write_mark <= real_option_buffer_size) then
                                                real_option_buffer(buffer_write_mark) = atof(arg_value)
                                                buffer_write_mark = buffer_write_mark + 1

                                                internal_arg_index = internal_arg_index + 1
                                        else
                                                print*, 'ERROR: Reading a real array command line &
                                                &parameter overan the given buffer of size', real_option_buffer_size
                                                print*, 'ERROR: Stopping.'
                                                stop
                                        endif
                                else
                                        exit outer
                                endif
                        enddo outer

                        if (buffer_write_mark == 1) then
                                call get_command_argument(arg_index - 1 ,parent_arg_value)
                                print*, 'ERROR: Absent or bad arguments for ', trim(parent_arg_value)
                                print*, 'ERROR: expected to populate an array from command line arguments, &
                                        &but found no valid numbers.'
                                stop
                        else
                                real_option = real_option_buffer(1:buffer_write_mark-1)
                        endif

                        if (increment_index_) then
                                arg_index = internal_arg_index
                        endif

                endfunction get_option_value_real_array_dp

                subroutine parse_command_line(self)

                        use env_kindtypes,    only: si
                        use core_convert,      only: itoa

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(runParameters),intent(out) :: self
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !to hold command line arguments

                        integer   (si)             :: i    ! Iteration variable
                        character (LEN=256)        :: arg  ! Holds read argument.

                        i = 1
                        do

                                arg = get_option_value_string(i,increment_index=.true.)

                                select case (arg)

                                case ('--linear-model')
                                        self%linearModel_flag =.true.

                                case ('--spatial-model')
                                        self%spatialModel_flag =.true.

                                case ('--spectral-model')
                                        self%spectralModel_flag =.true.

                                case ('--num-CG-sites')
                                        self%model_design%nSites = get_option_value_integer_si(i,&
                                                                                increment_index=.true.)
                                        self%n_cg_flag=.true.

                                case ('--residual-weights')
                                        self%model_design%residual_weights = get_option_value_real_array_dp(i,&
                                                                                increment_index=.true.)
                                        self%residual_weights_flag=.true.

                                case ('--input-psf')
                                        self%atomPsfFile = get_option_value_string(i,increment_index=.true.)
                                        self%atom_psf_flag=.true.

                                case ('--traj-preproc-map-type')
                                        self%traj_preproc_map_type = &
                                                get_option_value_string(i,increment_index=.true.)
                                        self%traj_preproc_map_type_flag =.true.

                                case ('--bs-log')
                                        self%bsLog    = get_option_value_string(i,increment_index=.true.)
                                        self%bs_log_flag=.true.

                                case ('--model-summary-log')
                                        self%model_summary_log = get_option_value_string(i,increment_index=.true.)
                                        self%model_summary_log_flag = .true.

                                case ('--resampling-summary-log')
                                        self%resampling_summary_log = get_option_value_string(i,increment_index=.true.)
                                        self%resampling_summary_log_flag =.true.

                                case ('--model-log')
                                        self%modelLog = get_option_value_string(i,increment_index=.true.)
                                        self%model_log_flag=.true.

                                case ('--input-dcd')
                                        self%atomDcdFile = get_option_value_string(i,increment_index=.true.)
                                        self%atom_dcd_flag=.true.

                                case ('--stride')
                                        self%deltaStep = get_option_value_integer_si(i,increment_index=.true.)
                                        self%delta_step_flag = .true.

                                case ('--sim-ann-steps')
                                        self%model_design%anneal_steps = get_option_value_integer_si(i,&
                                                                                        increment_index=.true.)
                                        self%max_sa_steps_flag=.true.

                                case ('--sim-ann-acceptance')
                                        self%model_design%anneal_accept = get_option_value_real_dp(i,&
                                                                                        increment_index=.true.)
                                        self%sa_accept_flag =.true.

                                case ('--sim-ann-control')
                                        self%model_design%anneal_param = get_option_value_real_dp(i,&
                                                                                        increment_index=.true.)
                                        self%sa_descent_control_flag =.true.

                                case ('--max-descent-steps')
                                        self%model_design%refine_steps = get_option_value_integer_si(i,increment_index=.true.)
                                        self%max_sd_steps_flag=.true.

                                case ('--num-models')
                                        self%nSets = get_option_value_integer_si(i,increment_index=.true.)
                                        if (self%nSets < 1) then
                                                print*, "WARNING: You specified an illegal number of sets:"//itoa(self%nSets)//'.'
                                                self%nSets = 1
                                                print*, "WARNING: Changing this to "//itoa(self%nSets)//'.'
                                        endif
                                        self%nsets_flag=.true.

                                case ('--num-cov-dof')
                                        self%trajPreproc_design%proj_num_dof = &
                                                get_option_value_integer_si(i,increment_index=.true.)

                                        if (self%trajPreproc_design%proj_num_dof > 0) then
                                                self%trajPreproc_design%proj_flag = .true.
                                        else
                                                self%trajPreproc_design%proj_flag = .false.
                                        endif

                                        self%ed_dof_flag=.true.

                                case ('--num-res-bs-reps')
                                        self%n_res_bs_reps    = get_option_value_integer_si(i,increment_index=.true.)
                                        self%n_res_bs_reps_Flag=.true.

                                case ('--num-res-bs-models')
                                        self%n_res_bs_models    = get_option_value_integer_si(i,increment_index=.true.)
                                        self%n_res_bs_models_Flag=.true.

                                case ('--num-bs-reps')
                                        self%n_bs_reps     = get_option_value_integer_si(i,increment_index=.true.)
                                        self%n_bs_reps_Flag=.true.

                                case ('--num-mc-norm-reps')
                                        self%n_mc_norm_reps     = get_option_value_integer_si(i,increment_index=.true.)
                                        self%n_mc_norm_reps_Flag=.true.

                                case ('--use-bca-ci')
                                        block
                                                integer (si) :: bca_ci_int
                                                bca_ci_int = get_option_value_integer_si(i,increment_index=.true.)
                                                if (bca_ci_int == 0) then
                                                        self%use_bca_ci = .false.
                                                else
                                                        self%use_bca_ci = .true.
                                                endif
                                        endblock

                                        self%use_bca_ci_flag =.true.

                                case ('--cluster-density-tol')
                                        self%model_design%cluster_density_tol = get_option_value_real_dp(i,increment_index=.true.)
                                        self%cluster_density_tol_Flag=.true.

                                case ('--kmeans-max-iter')
                                        self%model_design%kmeans_max_iter = get_option_value_integer_si(i,increment_index=.true.)
                                        self%kmeans_max_iter_flag =.true.

                                case ('--kmeans-max-try')
                                        self%model_design%kmeans_max_try = get_option_value_integer_si(i,increment_index=.true.)
                                        self%kmeans_max_try_flag =.true.

                                case ('--spectral-regularization-type')
                                        self%model_design%spectral_regularization_type = get_option_value_string(i,&
                                                                                                increment_index=.true.)
                                        self%spectral_regularization_type_flag =.true.

                                case ('--spectral-edcg-dist-conv-type')
                                        self%model_design%spectral_edcg_dist_conv_type  = get_option_value_string(i,&
                                                                                                increment_index=.true.)
                                        self%spectral_edcg_dist_conv_type_flag =.true.

                                        self%model_design%spectral_use_edcg_metric = .true.

                                case ('--spectral-spatial-dist-conv-type')
                                        self%model_design%spectral_spatial_dist_conv_type = get_option_value_string(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_spatial_dist_conv_type_flag =.true.

                                        self%model_design%spectral_use_spatial_metric = .true.

                                case ('--spectral-pairVar-dist-conv-type')
                                        self%model_design%spectral_pairVar_dist_conv_type = get_option_value_string(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_pairVar_dist_conv_type_flag =.true.

                                        self%model_design%spectral_use_pairVar_metric = .true.

                                case ('--spectral-charge-dist-conv-type')
                                        self%model_design%spectral_charge_dist_conv_type = get_option_value_string(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_charge_dist_conv_type_flag =.true.

                                        self%model_design%spectral_use_charge_metric = .true.

                                case ('--spectral-edcg-dist-conv-param')
                                        self%model_design%spectral_edcg_dist_conv_param = get_option_value_real_dp(i,&
                                                                                                increment_index=.true.)
                                        self%spectral_edcg_dist_conv_param_flag =.true.

                                        self%model_design%spectral_use_edcg_metric = .true.

                                case ('--spectral-spatial-dist-conv-param')
                                        self%model_design%spectral_spatial_dist_conv_param = get_option_value_real_dp(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_spatial_dist_conv_param_flag =.true.

                                        self%model_design%spectral_use_spatial_metric = .true.

                                case ('--spectral-pairVar-dist-conv-param')
                                        self%model_design%spectral_pairVar_dist_conv_param = get_option_value_real_dp(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_pairVar_dist_conv_param_flag =.true.

                                        self%model_design%spectral_use_pairVar_metric = .true.

                                case ('--spectral-charge-dist-conv-param')
                                        self%model_design%spectral_charge_dist_conv_param = get_option_value_real_dp(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_charge_dist_conv_param_flag =.true.

                                        self%model_design%spectral_use_charge_metric = .true.

                                case ('--spectral-edcg-kNN')
                                        self%model_design%spectral_edcg_kNN = get_option_value_integer_si(i,&
                                                                                                increment_index=.true.)
                                        self%spectral_edcg_kNN_flag =.true.

                                        self%model_design%spectral_use_edcg_metric = .true.

                                case ('--spectral-spatial-kNN')
                                        self%model_design%spectral_spatial_kNN = get_option_value_integer_si(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_spatial_kNN_flag =.true.

                                        self%model_design%spectral_use_spatial_metric = .true.

                                case ('--spectral-pairVar-kNN')
                                        self%model_design%spectral_pairVar_kNN = get_option_value_integer_si(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_pairVar_kNN_flag =.true.

                                        self%model_design%spectral_use_pairVar_metric = .true.

                                case ('--spectral-charge-kNN')
                                        self%model_design%spectral_charge_kNN = get_option_value_integer_si(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_charge_kNN_flag =.true.

                                        self%model_design%spectral_use_charge_metric = .true.

                                case ('--spectral-edcg-kNN-type')
                                        self%model_design%spectral_edcg_kNN_type = get_option_value_string(i,&
                                                                                                increment_index=.true.)
                                        self%spectral_edcg_kNN_type_flag =.true.

                                        self%model_design%spectral_use_edcg_metric = .true.

                                case ('--spectral-spatial-kNN-type')
                                        self%model_design%spectral_spatial_kNN_type = get_option_value_string(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_spatial_kNN_type_flag =.true.

                                        self%model_design%spectral_use_spatial_metric = .true.

                                case ('--spectral-pairVar-kNN-type')
                                        self%model_design%spectral_pairVar_kNN_type = get_option_value_string(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_pairVar_kNN_type_flag =.true.

                                        self%model_design%spectral_use_pairVar_metric = .true.

                                case ('--spectral-charge-kNN-type')
                                        self%model_design%spectral_charge_kNN_type = get_option_value_string(i,&
                                                                                                increment_index=.true.)

                                        self%spectral_charge_kNN_type_flag =.true.

                                        self%model_design%spectral_use_charge_metric = .true.

                                case ('--mapping-charge-split')
                                        self%mapping_charge_split_flag = .true.

                                        self%model_design%mapping_charge_split =.true.

                                case ('--model-sort-by')
                                        self%model_sort_by = get_option_value_string(i,increment_index=.true.)

                                        self%model_sort_by_flag =.true.

                                case ('--multipole-analysis')
                                        self%multipoleAnalysis_design%do_multipole_analysis = .true.

                                        self%multipole_analysis_flag =.true.

                                case ('--map-accumulator-type')
                                        self%model_design%map_accumulator_type = get_option_value_string(i,increment_index=.true.)

                                        self%map_accumulator_type_flag = .true.

                                case ('-f','--run-config')
                                        self%runConfig=.true.
                                        !do nothing else for now.

                                case default
                                        print '(a,a,/)', 'FATAL ERROR: Unrecognized command-line option: ', arg
                                        print*, 'INFO: Accepted options:  &
                                                &--linear-model &
                                                &--spatial-model &
                                                &--spectral-model &
                                                &--num-CG-sites &
                                                &--residual-weights &
                                                &--input-psf &
                                                &--bs-log &
                                                &--model-summary-log &
                                                &--resampling-summary-log &
                                                &--model-log &
                                                &--input-dcd &
                                                &--stride &
                                                &--sim-ann-steps &
                                                &--sim-ann-acceptance &
                                                &--sim-ann-control &
                                                &--max-descent-steps &
                                                &--num-models &
                                                &--num-cov-dof &
                                                &--num-res-bs-reps &
                                                &--num-res-bs-models &
                                                &--num-bs-reps &
                                                &--num-mc-norm-reps &
                                                &--use-bca-ci &
                                                &--cluster-density-tol &
                                                &--kmeans-max-iter &
                                                &--spectral-regularization-type &
                                                &--spectral-distance-conv-type &
                                                &--spectral-distance-conv-param.'
                                        stop
                                end select

                                if (i > command_argument_count()) exit
                        enddo

                endsubroutine parse_command_line

                subroutine check_set_values(self,apply_defaults)

                        use core_convert,     only: itoa, ftoa

                        implicit none

                        !!!!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(runParameters),intent(inout)          :: self
                        logical,             intent(in   ),optional :: apply_defaults
                        !!!!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        logical         :: apply_defaults_

                        !optional arg guards
                        if (present(apply_defaults)) then
                                apply_defaults_ = apply_defaults
                        else
                                apply_defaults_ = .false.
                        endif

                        if (.not.(self%spatialModel_flag &
                                  .or. self%linearModel_flag &
                                  .or. self%spectralModel_flag)) then
                                write(*,*) "WARNING: No model type was specified on the command line. Set this&
                                                &with --spatial-model, --linear-model or --spectral-model."
                                if (apply_defaults_) then
                                        self%modelType = 'linearDivision'
                                        write(*,*) 'INFO: Using a linear division type model for analysis.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        elseif (count([self%spatialModel_flag,self%linearModel_flag,self%spectralModel_flag])&
                                      /= 1) then
                                write(*,*) "ERROR: Conflicting model types specified. &
                                                &Check your use of --spatial-model, --linear-model &
                                                &and --spectral-model."
                                stop
                        elseif (self%linearModel_flag) then
                                write(*,*) "INFO: Using a linear model."
                                self%modelType = 'linearDivision'
                        elseif (self%spatialModel_flag)     then
                                write(*,*) "INFO: Using a centroid based spatial model."
                                self%modelType = 'centroid'
                        elseif (self%spectralModel_flag)     then
                                write(*,*) "INFO: Using a spectral model."
                                self%modelType = 'spectral'

                                if (.not. (self%model_design%spectral_use_edcg_metric    .or. &
                                           self%model_design%spectral_use_spatial_metric .or. &
                                           self%model_design%spectral_use_pairVar_metric .or. &
                                           self%model_design%spectral_use_charge_metric ))     then
                                        write(*,*) "INFO: No metric for spectral clustering specified. Enabling &
                                                   &edcg metric by default."
                                        self%model_design%spectral_use_edcg_metric = .true.
                                        self%spectral_use_edcg_metric_flag = .true.
                                endif
                        endif

                        if (self%n_cg_flag.eqv..false.) then
                                write(*,*) 'ERROR: Must provide number of CG sites using command line argument -ncg [number of CG sites]'
                                stop
                        endif

                        if (self%atom_psf_flag.eqv..false.) then
                                write(*,*) 'ERROR: Must provide a psf file using command line argument -psf [psf file name]'
                                stop
                        endif

                        if (self%atom_dcd_flag.eqv..false.) then
                                write(*,*) 'ERROR: Must provide a atom dcd file using command line argument -dcd [atom dcd file name]'
                                stop
                        endif

                        if (self%traj_preproc_map_type_flag .eqv..false.) then
                                write(*,*) "WARNING: The preprocesing option controlling how the trajectory &
                                           &is mapped before the model is made was not specified. This can be &
                                           &set with --traj-preproc-map-type [type], with [type] as &
                                           &carbon_alpha, identity, or residue_com."
                                if (apply_defaults_) then
                                        self%traj_preproc_map_type = 'carbon_alpha'
                                        write(*,*) 'INFO: Using trajectory preprocessing type'&
                                                  //self%traj_preproc_map_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%model_log_flag.eqv..false.) then
                                write(*,*) "WARNING: A file to record the native models was not set via the command line.&
                                           &Manually set this with the command line argument &
                                           &--model-log [logfile]"
                                if (apply_defaults_) then
                                        self%modelLog = 'models.log'
                                        write(*,*) 'INFO: Using model log file '//self%modelLog//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%bs_log_flag.eqv..false.) then
                                write(*,*) "WARNING: A log file for analysis of the algorithm itself via &
                                           &boostrapping and jacknifing was not set via the command line. &
                                           &Manually set this with the command line argument &
                                           &--bs-log [logfile]"
                                if (apply_defaults_) then
                                        self%bsLog = 'bs.log'
                                        write(*,*) 'INFO: Using log file '//self%bsLog//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%model_summary_log_flag.eqv..false.) then
                                write(*,*) "WARNING: A log file for summarizing the model statistics was not set. &
                                           &If no model was done, this is fine. Else, set this log &
                                           &with --model-summary-log [logfile]."
                                if (apply_defaults_) then
                                        self%model_summary_log = 'model_summary.log'
                                        write(*,*) 'INFO: Using model log file '//self%model_summary_log//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%resampling_summary_log_flag.eqv..false.) then
                                write(*,*) "WARNING: A log file for summarizing the resampling was not set. &
                                           &If no resampling was done, this is fine. Else, set this log &
                                           &with --resampling-summary-log [logfile]."
                                if (apply_defaults_) then
                                        self%resampling_summary_log = 'reasmpling_summary.log'
                                        write(*,*) 'INFO: Using resampling log file '//self%resampling_summary_log//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%delta_step_flag.eqv..false.) then
                                write(*,*) "WARNING: A sampling step size (stride) was not set via the command line. &
                                           &Manually set this with the command line argument &
                                           &--stride [step size]"
                                if (apply_defaults_) then
                                        self%deltaStep = 1
                                        write(*,*) 'INFO: Using a default stride (step size) of '//itoa(self%deltaStep)//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%residual_weights_flag.eqv..false.) then
                                write(*,*) "WARNING: Residual weights were not set via the command line. &
                                           &Provide weights be specifying the corresponding numbers in order. &
                                           &First is fit charge, second is EDCG, third is mapped charged. &
                                           &--residual-weights [coefficients]"
                                if (apply_defaults_) then
                                        self%model_design%residual_weights = [ 0.0_dp, 1.0_dp, 0.0_dp ]
                                        write(*,*) 'INFO: Using a default charge residual weight of '&
                                                   //ftoa(self%model_design%residual_weights(1))//' '&
                                                   //ftoa(self%model_design%residual_weights(2))//' '&
                                                   //ftoa(self%model_design%residual_weights(3))//' '&
                                                   //'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        else
                                if (size(self%model_design%residual_weights,1) /= 3) then
                                        write(*,*) "ERROR: Residual weights were set, but not the right number of them. &
                                                   &Expected 3, got "//itoa(size(self%model_design%residual_weights,1))//'.'
                                        write(*,*) "ERROR: Stopping."
                                        stop
                                endif
                                if (abs(sum(self%model_design%residual_weights) - 1.0_dp) > (10.0_dp ** (-12))) then
                                        if (sum(self%model_design%residual_weights) < (10.0_dp**(-12))) then
                                                write(*,*) "WARNING: Residual weights were set, but did not sum to 1. "
                                        endif
                                        write(*,*) "WARNING: Residual weights were set, but did not sum to 1. "
                                        write(*,*) "WARNING: Residual weights (sum): ", self%model_design%residual_weights, &
                                                   '(', ftoa(sum(self%model_design%residual_weights)), ')'

                                        self%model_design%residual_weights = self%model_design%residual_weights/sum(self%model_design%residual_weights)
                                        write(*,*) "WARNING: Normalizing given residual_weights to", self%model_design%residual_weights, '.'
                                endif
                        endif

                        if (self%max_sa_steps_flag.eqv..false.) then
                                write(*,*) "WARNING: Maxium number of simulated annealing steps to do per set was not set. &
                                           &Manually set this with the command line argument &
                                           &--sim-ann-steps [number of steps]"
                                if (apply_defaults_) then
                                        self%model_design%anneal_steps = 10000
                                        write(*,*) 'INFO: Using default '//itoa(self%model_design%anneal_steps)//' simulated annealing steps.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%sa_accept_flag.eqv..false.) then
                                write(*,*) "WARNING: Initial uphill acceptance rate for simulated annealing wasn't set. &
                                           &This can be set via --sim-ann-acceptance [acceptance decimal] "
                                self%model_design%anneal_accept = .95_dp
                                if (apply_defaults_) then
                                        write(*,*) 'INFO: Using a default '//ftoa(self%model_design%anneal_accept)//' for the acceptance ratio.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%sa_descent_control_flag.eqv..false.) then
                                write(*,*) "WARNING: The parameter which controls the rate of simulated annealing temperature relaxation &
                                           &was not set via --sim-ann-control [control value]."
                                if (apply_defaults_) then
                                        self%model_design%anneal_param = .0007
                                        write(*,*) 'INFO: Using default '//ftoa(self%model_design%anneal_param)//' for relaxation value.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%max_sd_steps_flag.eqv..false.) then
                                write(*,*) "WARNING: Maxium number of refine/descent steps to do per set was not set. &
                                           &Manually set this with the command line argument &
                                           &--max-descent-steps [number of steps]"
                                if (apply_defaults_) then
                                        self%model_design%refine_steps = 2000
                                        write(*,*) 'INFO: Using default '//itoa(self%model_design%refine_steps)//' maximum steepest descent steps.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%nsets_flag.eqv..false.) then
                                write(*,*) "WARNING: Number of starting conditions to run native optimization from was not set. &
                                           &Manually set this with the command line argument &
                                           &--num-models [number of models]"
                                if (apply_defaults_) then
                                        self%nSets    = 15
                                        write(*,*) 'INFO: Using default '//itoa(self%nSets)//' repitions.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%ed_dof_flag.eqv..false.) then
                                write(*,*) 'WARNING: Number of modes to use in the essential subspace projection was not set. &
                                           &Manually set this with the command line argument &
                                           &--num-cov-dog [degrees of freedom]. Set this to -1 to use all degrees.'
                                if (apply_defaults_) then

                                        self%trajPreproc_design%proj_num_dof = 3*self%model_design%nSites - 6
                                        self%trajPreproc_design%proj_flag    = .true.

                                        write(*,*) 'INFO: Using a default 3*nSites-6 degrees of freedom.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        else
                                self%trajPreproc_design%proj_flag    = .false.
                        endif

                        if (self%n_res_bs_reps_Flag .eqv..false.) then
                                write(*,*) 'WARNING: Number of bootstrapping repitions to use in residual confidence interval&
                                           & (NOT GENERAL BOOTSTRAPPING)&
                                           & calculation was not set.&
                                           & Manually set this with the command line argument &
                                           & --num-residual-bs-reps [number of reps].'
                                if (apply_defaults_) then
                                        self%n_res_bs_reps = 0
                                        write(*,*) 'INFO: Using a default of '//itoa(self%n_res_bs_reps)//'reps.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%n_res_bs_models_Flag .eqv. .false.) then
                                write(*,*) 'WARNING: Number of top unique models to calculate bootstrapped&
                                           & confidence intervals for the residual (NOT GENERAL BOOTSTRAPPING)&
                                           & was not specified. This is a costly calculation.&
                                           & Manually set this with the command line argument&
                                           & --num_res_bs_models [number of reps].'
                                if (apply_defaults_) then
                                        self%n_res_bs_models = 0
                                        write(*,*) 'INFO: Calculating for a default '//itoa(self%n_res_bs_models)//'model.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%n_bs_reps_Flag .eqv. .false.) then
                                write(*,*) 'WARNING: Number of draws to use for general bootstrap calculations&
                                           & was not specified. This is a very costly calculation!&
                                           & Manually set this with the command line argument&
                                           & --num_bs_reps [number of reps].'
                                if (apply_defaults_) then
                                        self%n_res_bs_models = 0 
                                        write(*,*) 'INFO: Calculating for a default '//itoa(self%n_res_bs_models)//'model.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%n_mc_norm_reps_flag .eqv. .false.) then
                                write(*,*) 'WARNING: Number of Monte Carlo repitions used to calculate &
                                           &normalization constants not specified. Specify this with &
                                           &--num-mc-norm-reps [number of reps].'
                                if (apply_defaults_) then
                                        self%n_mc_norm_reps = 0
                                        write(*,*) 'INFO: Using a default '//itoa(self%n_mc_norm_reps)//' repititions.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%use_bca_ci_flag .eqv. .false.) then
                                write(*,*) 'WARNING: BCA type confidence intervals were not enabled nor disabled.'
                                if (apply_defaults_) then
                                        self%use_bca_ci = .false.
                                        write(*,*) 'INFO: Disabling BCA type confidence intervals as per default.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%cluster_density_tol_flag .eqv. .false.) then
                                write(*,*) 'WARNING: A tolerance for the spatial density of clusters accepted in moves &
                                           &was not set. This only matters if you are using a spatial centroid model. &
                                           &Set it with --cluster-density-tol [tolerance].'
                                if (apply_defaults_) then
                                        self%model_design%cluster_density_tol = .05_dp
                                        write(*,*) 'INFO: Setting the cluster spatial density tolerance &
                                        &to '//ftoa(self%model_design%cluster_density_tol)//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%kmeans_max_iter_flag .eqv. .false.) then
                                write(*,*) 'WARNING: A maximum for allowed kmeans iterations was not set. &
                                           &This matters for centroid spatial and spectral models. Set it &
                                           &with --max-kmeans-iter [number of iterations].'
                                if (apply_defaults_) then
                                        self%model_design%kmeans_max_iter = 100
                                        write(*,*) 'INFO: Setting the max number of kmeans iterations &
                                        &to '//itoa(self%model_design%kmeans_max_iter)//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%kmeans_max_try_flag .eqv. .false.) then
                                write(*,*) "WARNING: A maximum for allowed kmeans tries was not set. &
                                           &This isn't the number of kmeans iterations. &
                                           &This matters for spectral models. Set it &
                                           &with --max-kmeans-try [number of tryations]."
                                if (apply_defaults_) then
                                        self%model_design%kmeans_max_try = 150
                                        write(*,*) 'INFO: Setting the max number of kmeans tries &
                                        &to '//itoa(self%model_design%kmeans_max_try)//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_regularization_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no spectral regularization was set. This matters if using &
                                           &a spectral model. Set it via --spectral-regularization-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_regularization_type = 'unregularized'
                                        write(*,*) 'INFO: Setting spectral regularization type &
                                        &to '//self%model_design%spectral_regularization_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_edcg_dist_conv_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no conversion type transforming a distance matrix &
                                           &to an affinity matrix was set. This matters if using spectral clustering. &
                                           &Set this via --spectral-edcg-dist-conv-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_edcg_dist_conv_type = 'trans'
                                        write(*,*) 'INFO: Setting distance conversion type &
                                        &to '//self%model_design%spectral_edcg_dist_conv_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_spatial_dist_conv_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no conversion type transforming a distance matrix &
                                           &to an affinity matrix FOR THE REAL SPACE DISTANCE was set. &
                                           &This matters if using spectral clustering and using REAL SPACE DISTANCES. &
                                           &Set this via --spectral-spatial-dist-conv-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_spatial_dist_conv_type = 'trans'
                                        write(*,*) 'INFO: Setting spatial distance conversion type &
                                        &to '//self%model_design%spectral_spatial_dist_conv_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_pairVar_dist_conv_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no conversion type transforming a distance matrix &
                                           &to an affinity matrix for pairwise variance similarity was set. &
                                           &Set this via --spectral-pairVar-dist-conv-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_pairVar_dist_conv_type = 'trans'
                                        write(*,*) 'INFO: Setting pairVar distance conversion type &
                                        &to '//self%model_design%spectral_pairVar_dist_conv_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_charge_dist_conv_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no conversion type transforming a distance matrix &
                                           &to an affinity matrix for charge similarity was set. &
                                           &Set this via --spectral-charge-dist-conv-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_charge_dist_conv_type = 'trans'
                                        write(*,*) 'INFO: Setting charge distance conversion type &
                                        &to '//self%model_design%spectral_charge_dist_conv_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_edcg_dist_conv_param_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no parameter controlling the distance matrix to affinity matrix &
                                        &transformation was given. This matters if using a spectral model, and the effect &
                                        &of this parameter depends on the conversion method used. &
                                        &Set this via --spectral-edcg-dist-conv-param [parameter].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_edcg_dist_conv_param = 0
                                        write(*,*) 'INFO: Setting distance conversion param &
                                        &to '//ftoa(self%model_design%spectral_edcg_dist_conv_param)//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_spatial_dist_conv_param_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no parameter controlling the distance matrix to affinity matrix &
                                        &transformation FOR THE REAL SPACE DISTANCE was given. This matters if using &
                                        &a spectral model and using real space values, and the effect &
                                        &of this parameter depends on the conversion &
                                        &method used. Set this via --spectral-spatial-dist-conv-param [parameter].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_spatial_dist_conv_param = 0
                                        write(*,*) 'INFO: Setting spatial distance conversion param &
                                        &to '//ftoa(self%model_design%spectral_spatial_dist_conv_param)//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_pairVar_dist_conv_param_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no parameter controlling the distance matrix to affinity matrix &
                                        &transformation pairwise variance was given. &
                                        &Set this via --spectral-pairVar-dist-conv-param [parameter].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_pairVar_dist_conv_param = 0
                                        write(*,*) 'INFO: Setting pairVar distance conversion param &
                                        &to '//ftoa(self%model_design%spectral_pairVar_dist_conv_param)//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_charge_dist_conv_param_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no parameter controlling the distance matrix to affinity matrix &
                                        &transformation for charge difference was given. &
                                        &Set this via --spectral-charge-dist-conv-param [parameter].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_charge_dist_conv_param = 0
                                        write(*,*) 'INFO: Setting charge distance conversion param &
                                        &to '//ftoa(self%model_design%spectral_charge_dist_conv_param)//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif


                        if (self%spectral_edcg_kNN_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no type for kNN based filtering set (used in spectral clustering). &
                                           &Set this via --spectral-edcg-kNN-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_edcg_kNN_type = 'mutual'
                                        write(*,*) 'INFO: Setting k Nearest Neighbor type &
                                        &to '//self%model_design%spectral_edcg_kNN_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_spatial_kNN_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no type for spatial kNN based filtering set (used in spectral &
                                           &clustering). Set this via --spectral-spatial-kNN-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_spatial_kNN_type = 'mutual'
                                        write(*,*) 'INFO: Setting spatial k Nearest Neighbor type &
                                        &to '//self%model_design%spectral_spatial_kNN_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_pairVar_kNN_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no type for pairVar kNN based filtering set (used in spectral &
                                           &clustering). Set this via --spectral-pairVar-kNN-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_pairVar_kNN_type = 'mutual'
                                        write(*,*) 'INFO: Setting pairVar k Nearest Neighbor type &
                                        &to '//self%model_design%spectral_pairVar_kNN_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_charge_kNN_type_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no type for charge kNN based filtering set (used in spectral &
                                           &clustering). Set this via --spectral-charge-kNN-type [type].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_charge_kNN_type = 'mutual'
                                        write(*,*) 'INFO: Setting charge k Nearest Neighbor type &
                                        &to '//self%model_design%spectral_charge_kNN_type//'.'
                                else
                                        stop
                                endif
                        endif

                        if (self%spectral_edcg_kNN_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no value for setting k nearest neighbors (used in &
                                           &spectral clustering). Set this via --spectral-edcg-kNN [#].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_edcg_kNN = -1
                                        write(*,*) 'INFO: Setting k nearest-neighbors &
                                        &to '//itoa(self%model_design%spectral_edcg_kNN)//'. &
                                        &This deactivates kNN filtering.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_spatial_kNN_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no value for setting k nearest neighbors for spatial evaluation &
                                           &(used in spectral clustering). Set this via --spectral-spatial-kNN [#].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_spatial_kNN = -1
                                        write(*,*) 'INFO: Setting spatial k nearest-neighbors &
                                        &to '//itoa(self%model_design%spectral_spatial_kNN)//'. &
                                        &This deactivates kNN filtering.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_pairVar_kNN_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no value for setting k nearest neighbors for pairVar evaluation &
                                           &(used in spectral clustering). Set this via --spectral-pairVar-kNN [#].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_pairVar_kNN = -1
                                        write(*,*) 'INFO: Setting pairVar k nearest-neighbors &
                                        &to '//itoa(self%model_design%spectral_pairVar_kNN)//'. &
                                        &This deactivates kNN filtering.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%spectral_charge_kNN_flag .eqv. .false.) then
                                write(*,*) 'WARNING: no value for setting k nearest neighbors for charge evaluation &
                                           &(used in spectral clustering). Set this via --spectral-charge-kNN [#].'
                                if (apply_defaults_) then
                                        self%model_design%spectral_charge_kNN = -1
                                        write(*,*) 'INFO: Setting charge k nearest-neighbors &
                                        &to '//itoa(self%model_design%spectral_charge_kNN)//'. &
                                        &This deactivates kNN filtering.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%mapping_charge_split_flag .eqv. .false.) then
                                write(*,*) "WARNING: Whether to split each bead determined by the first round of &
                                           &of mapping generation wasn\'t specified. Enable splitting with &
                                           &--model-charge-split."
                                if (apply_defaults_) then
                                        self%model_design%mapping_charge_split = .false.
                                        write(*,*) 'INFO: Deactivating charge splitting.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%model_sort_by_flag .eqv. .false.) then
                                if (apply_defaults_) then
                                        self%model_sort_by = 'residual'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                        if (self%multipoleAnalysis_design%do_multipole_analysis) then
                                write(*,*) "INFO: Performing multipole analysis."
                                self%multipoleAnalysis_design%site_dipole_traj_filename = &
                                        "site_dipole_traj.csv"
                                self%multipoleAnalysis_design%site_quadrupole_traj_filename = &
                                        "site_quadrupole_traj.csv"

                                self%multipoleAnalysis_design%cg_site_dipole_traj_filename = &
                                        "cg_site_dipole_traj.csv"
                                self%multipoleAnalysis_design%cg_site_quadrupole_traj_filename = &
                                        "cg_site_quadrupole_traj.csv"

                                self%multipoleAnalysis_design%system_dipole_traj_filename = &
                                        "system_dipole_traj.csv"
                                self%multipoleAnalysis_design%system_quadrupole_traj_filename = &
                                        "system_quadrupole_traj.csv"

                                self%multipoleAnalysis_design%cg_system_dipole_traj_filename = &
                                        "cg_system_dipole_traj.csv"
                                self%multipoleAnalysis_design%cg_system_quadrupole_traj_filename = &
                                        "cg_system_quadrupole_traj.csv"

                                self%multipoleAnalysis_design%nosplit_cg_system_dipole_traj_filename = &
                                        "nosplit_cg_system_dipole_traj.csv"
                                self%multipoleAnalysis_design%nosplit_cg_system_quadrupole_traj_filename = &
                                        "nosplit_cg_system_quadrupole_traj.csv"

                        endif

                        if (self%map_accumulator_type_flag .eqv. .false.) then
                                if (apply_defaults_) then
                                        self%model_design%map_accumulator_type = 'com'
                                        write(*,*) 'WARNING: No acculator type for applying mappings was specified &
                                                   &e.g. com (center of mass). Set this with --map-accumulator-type [type].'
                                        write(*,*) 'INFO: Map accumulator type set to '//self%model_design%map_accumulator_type//'.'
                                else
                                        write(*,*) 'ERROR: Lacking vital option and no default. Stopping.'
                                        stop
                                endif
                        endif

                endsubroutine check_set_values

endmodule obj_commandline

