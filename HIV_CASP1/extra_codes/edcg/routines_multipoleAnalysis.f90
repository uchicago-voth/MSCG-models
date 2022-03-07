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
! Provides routines for multipole analysis of systems and subsystems
! of trajectories.
!

module routines_multipoleAnalysis

        implicit none

        private

        public multipole_analysis, multipoleAnalysis_config

        !gives values that control filenames and (later) parameters for
        !multipole analysis. Keeps function signatures sane.
        type         :: multipoleAnalysis_config
                logical                  :: do_multipole_analysis = .false.

                !hold trajectories of multipoles of subsets of AA trajectories.
                !not active if no charge splitting (they could be made active).
                character(:),allocatable :: site_dipole_traj_filename
                character(:),allocatable :: site_quadrupole_traj_filename

                !hold trajectories of multipoles of subsets of CG trajectories--
                !determined by how the beads were charge split. If not charge split,
                !these aren't used.
                character(:),allocatable :: cg_site_dipole_traj_filename
                character(:),allocatable :: cg_site_quadrupole_traj_filename

                !Overall system multipoles of AA system.
                character(:),allocatable :: system_dipole_traj_filename
                character(:),allocatable :: system_quadrupole_traj_filename

                !Overall system multipoles of production CG system.
                character(:),allocatable :: cg_system_dipole_traj_filename
                character(:),allocatable :: cg_system_quadrupole_traj_filename

                !Overall system multipoles of non-charge-split CG system.
                !Not active if no charge splitting.
                character(:),allocatable :: nosplit_cg_system_dipole_traj_filename
                character(:),allocatable :: nosplit_cg_system_quadrupole_traj_filename
        end type        multipoleAnalysis_config

        contains
                !Wrapper function for multipole analysis. Currently, we only take one 
                !model to analyze, but we can expand later.
                !
                ! @models:      array of models to analyze (only size 1 arrays allowed currently).
                ! @trajectory:  trajectory to use for analysis
                ! @design:      options to control analysis
                ! @nCAranks:    coarray ranks (futher compatibility)
                subroutine multipole_analysis(models, trajectory, design, nCAranks)

                        use env_kindtypes,      only: si
                        use abs_obj_model,      only: model
                        use obj_trajectory,     only: labeledTraj

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),                  intent(in   ) :: models(:)
                        class(labeledTraj),            intent(in   ) :: trajectory
                        type(multipoleAnalysis_config),intent(in   ) :: design
                        integer (si),                  intent(in   ) :: nCAranks
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        if (size(models) == 0) then
                                return
                        endif

                        if (nCAranks <= 1) then
                                 call multipole_analysis_smp(models,trajectory,design)
                        endif

                endsubroutine multipole_analysis

                !smp function for multipole analysis. Currently, we only take one 
                !model to analyze, but we can expand later.
                !
                ! @models:      array of models to analyze (only size 1 arrays allowed currently).
                ! @trajectory:  trajectory to use for analysis
                ! @design:      options to control analysis
                subroutine multipole_analysis_smp(models,trajectory,design)

                        use env_kindtypes,      only: si
                        use abs_obj_model,      only: model
                        use obj_trajectory,     only: traj, labeledTraj
                        use core_filter,        only: canonicizeMapping, mapProject

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        class(model),                  intent(in   ) :: models(:)
                        class(labeledTraj),            intent(in   ) :: trajectory
                        type(multipoleAnalysis_config),intent(in   ) :: design
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        type(labeledtraj)        :: cg_traj
                        type(traj)               :: stat_traj
                        integer (si),allocatable :: fine_map(:), coarse_map(:), link_map(:)

                        !check if we should do anything.
                        if (.not. design%do_multipole_analysis) then
                                return
                        endif

                        if (size(models) > 1) then
                                print*, "ERROR: Multipole analysis attempted on more than one model. Exiting."
                                stop
                        endif

                        !Much of this is just a list of unconnected commands to produce information to write.

                        !if we are splitting on charge, then we have additional actions.
                        if (models(1)%hasParentMap() .and. models(1)%isMapChargeSplit()) then
                                !do per bead analysis

                                !first, get non charge map of non charge split-- we calculate the dipole for each of 
                                !_these beads_.
                                coarse_map = models(1)%getMapping(parent_mapping=.true.,split_on_charge=.false.)
                                call canonicizeMapping(coarse_map)

                                !get the charge split map too.
                                fine_map   = models(1)%getMapping(parent_mapping=.true.,split_on_charge=.true.)

                                !per bead AA dipole.
                                stat_traj     = trajectory%genDipoleTraj(section_mapping = coarse_map,&
                                                                            site_weights = trajectory%site_mass)
                                call stat_traj%writeCSVfile(design%site_dipole_traj_filename)

                                !per bead AA quadrupole.
                                stat_traj     = trajectory%genQuadrupoleTraj(section_mapping = coarse_map,&
                                                                                site_weights = trajectory%site_mass)
                                call stat_traj%writeCSVfile(design%site_quadrupole_traj_filename)

                                !Now, get statistics for the split coarse CG beads from charge splitting.

                                !get the charge split CG trajectory.
                                select case (models(1)%getMapAccumType())
                                case ("cop")
                                        cg_traj  = trajectory%lmapTrajectory(fine_map)
                                case ("com")
                                        cg_traj  = trajectory%lmapTrajectory(fine_map,position_weights=trajectory%site_mass)
                                case ("coc")
                                        cg_traj  = trajectory%lmapTrajectory(fine_map,position_weights=trajectory%atomCharges)
                                endselect

                                !get map determining which charge split beads belong to coarse CG beads.
                                link_map = mapProject(basis_map      = fine_map,&
                                                      map_to_project = coarse_map)

                                !per coarse CG bead dipoles.
                                stat_traj     = cg_traj%genDipoleTraj(section_mapping = link_map,&
                                                                         site_weights = cg_traj%site_mass)
                                call stat_traj%writeCSVfile(design%cg_site_dipole_traj_filename)

                                !per coarse CG bead quadrupoles.
                                stat_traj     = cg_traj%genQuadrupoleTraj(section_mapping = link_map,&
                                                                             site_weights = cg_traj%site_mass)
                                call stat_traj%writeCSVfile(design%cg_site_quadrupole_traj_filename)

                                !get non charge split cg traj to do SYSTEM (not bead) analysis as comparison to later
                                !derived values.
                                !get the charge split CG trajectory.
                                select case (models(1)%getMapAccumType())
                                case ("cop")
                                        cg_traj  = trajectory%lmapTrajectory(coarse_map)
                                case ("com")
                                        cg_traj  = trajectory%lmapTrajectory(coarse_map,position_weights=trajectory%site_mass)
                                case ("coc")
                                        cg_traj  = trajectory%lmapTrajectory(coarse_map,position_weights=trajectory%atomCharges)
                                endselect

                                !Derive map for full system analysis.
                                deallocate(coarse_map)
                                allocate(coarse_map(cg_traj%nAtoms))
                                coarse_map = 1

                                !get coarse CG system multiple trajectories.
                                stat_traj     = cg_traj%genDipoleTraj(section_mapping = coarse_map,&
                                                                         site_weights = cg_traj%site_mass)
                                call stat_traj%writeCSVfile(design%nosplit_cg_system_dipole_traj_filename)

                                stat_traj     = cg_traj%genQuadrupoleTraj(section_mapping = coarse_map,&
                                                                             site_weights = cg_traj%site_mass)
                                call stat_traj%writeCSVfile(design%nosplit_cg_system_quadrupole_traj_filename)

                                deallocate(coarse_map)
                        endif

                        !Here, we get total system multiples for the fine CG and AA representations.

                        allocate(coarse_map(trajectory%nAtoms))
                        coarse_map = 1

                        !get AA system multipoles
                        stat_traj     = trajectory%genDipoleTraj(section_mapping = coarse_map,&
                                                                    site_weights = trajectory%site_mass)
                        call stat_traj%writeCSVfile(design%system_dipole_traj_filename)

                        stat_traj     = trajectory%genQuadrupoleTraj(section_mapping = coarse_map,&
                                                                        site_weights = trajectory%site_mass)
                        call stat_traj%writeCSVfile(design%system_quadrupole_traj_filename)

                        !get CG traj.
                        fine_map   = models(1)%getMapping(parent_mapping=.true.)

                        select case (models(1)%getMapAccumType())
                        case ("cop")
                                cg_traj  = trajectory%lmapTrajectory(fine_map)
                        case ("com")
                                cg_traj  = trajectory%lmapTrajectory(fine_map,position_weights=trajectory%site_mass)
                        case ("coc")
                                cg_traj  = trajectory%lmapTrajectory(fine_map,position_weights=trajectory%atomCharges)
                        endselect

                        deallocate(coarse_map)
                        allocate(coarse_map(cg_traj%nAtoms))
                        coarse_map = 1

                        !get CG system multipoles
                        stat_traj     = cg_traj%genDipoleTraj(section_mapping = coarse_map,&
                                                                 site_weights = cg_traj%site_mass)
                        call stat_traj%writeCSVfile(design%cg_system_dipole_traj_filename)

                        stat_traj     = cg_traj%genQuadrupoleTraj(section_mapping = coarse_map,&
                                                                     site_weights = cg_traj%site_mass)
                        call stat_traj%writeCSVfile(design%cg_system_quadrupole_traj_filename)

                endsubroutine multipole_analysis_smp
endmodule routines_multipoleAnalysis
