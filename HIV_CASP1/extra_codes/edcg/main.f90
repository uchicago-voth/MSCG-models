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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                       !!
!! Coarse Grained Mapping Generation Code                                                !!
!!                                                                                       !!
!! This code is property of the Voth Group at the University of Chicago (c) 05-01-2016   !!
!! and is released under the Apache license (see LICENSE file).                          !!
!!                                                                                       !!
!! Origial Author: Martin McCullagh (martin.mccullagh@colostate.edu)                     !!
!! Current Maintainer: Aleksander Durumeric (alekepd@gmail.com)                          !!
!!                                                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program map_optimize

        use env_kindtypes,              only: dp
        use obj_commandline,            only: runParameters
        use obj_stopwatch,              only: stopwatch
        use obj_trajectory,             only: traj, labeledTraj, traj_training_data
        use abs_obj_model,              only: model, write_model, write_model_stats
        use routines_modelDesign,       only: createModelSet, populateModelSetNorms
        use routines_modelSolve,        only: solve_models
        use routines_bootstrap,         only: calc_bootstrap_CIs, solve_resampled_models
        use routines_multipoleAnalysis, only: multipole_analysis

        implicit none

        !!! full objects.
        type(labeledTraj)          :: full_traj
        type(traj)                 :: cg_traj
        type(traj_training_data)   :: training_data
        type(runParameters)        :: cmd_line
        type(stopwatch)            :: timer
        class(model),allocatable   :: models(:), best_models(:), bs_models(:)
        real (dp),   allocatable   :: optimism_record(:)

        call timer%start()

        call cmd_line%parse_command_line()
        call cmd_line%check_values(apply_defaults=.true.)

        !add files 
        call full_traj%addPsfFile(cmd_line%atomPsfFile)
        call full_traj%addDcdFile(cmd_line%atomDcdFile)

        !map trajectory as a preprocessing step before applying models.
        training_data%trj = full_traj%namedMapTrajectory(map_type=cmd_line%traj_preproc_map_type)

        !preprocess trajectory (generated cov mats, etc).
        call training_data%preprocess(cmd_line%trajPreproc_design)

        !create total set of models we will run.
        call createModelSet(models        = models,&
                            nModels       = cmd_line%nSets,&
                            modelType     = cmd_line%modelType,&
                            model_design  = cmd_line%model_design)

        !calculate normalization constants we need later.
        call populateModelSetNorms(models = models,&
                                   tdata  = training_data,&
                                   n_reps = cmd_line%n_mc_norm_reps)

        !calculate models on bootstrapped trajectories for CI intervals.
        call solve_resampled_models(result_models = bs_models,&
                                    optimism_dist = optimism_record,&
                                    this_model    = models(1),&
                                    tdata         = training_data,&
                                    nDraws        = cmd_line%n_bs_reps,&
                                    nCAranks      = 1)

        call write_model_stats(models           = bs_models,&
                               optimism_dist    = optimism_record,&
                               report_file      = cmd_line%resampling_summary_log)

        call write_model(bs_models,cmd_line%bsLog,verbose=.true.)

        !solve legitimate model repititions.
        call solve_models(models          = models,&
                          result_count    = size(models),&
                          result_models   = best_models, &
                          training_tdatas = [ training_data ],&
                          dedup           = .true.,&
                          sort_by         = cmd_line%model_sort_by,&
                          nCAranks        = 1)


        !Calculate CIs for _residual values_.
        call calc_bootstrap_CIs(models   = best_models,&
                                nModels  = cmd_line%n_res_bs_models,&
                                tdata    = training_data,&
                                nBSreps  = cmd_line%n_res_bs_reps,&
                                bca      = cmd_line%use_bca_ci,&
                                nCAranks = 1)

        call write_model(best_models,cmd_line%modelLog,verbose=.true.)

        call write_model_stats(models      = best_models,&
                               report_file = cmd_line%model_summary_log)

        !apply the map our optimal model gave us to out original trajectory.
        select case (best_models(1)%getMapAccumType())
        case ("cop")
                cg_traj = full_traj%mapTrajectory(map = best_models(1)%getMapping(parent_mapping=.true.))
        case ("com")
                cg_traj = full_traj%mapTrajectory(map = best_models(1)%getMapping(parent_mapping=.true.),&
                                         position_weights = full_traj%site_mass)
        case ("coc")
                cg_traj = full_traj%mapTrajectory(map = best_models(1)%getMapping(parent_mapping=.true.),&
                                         position_weights = full_traj%atomCharges)
        endselect

        call cg_traj%writeXYZfile("minCgtraj.xyz")
        call multipole_analysis(best_models(1:1),full_traj,cmd_line%multipoleAnalysis_design,1)

        call timer%finish()
        call timer%report(prefix=' INFO: ')

endprogram map_optimize
