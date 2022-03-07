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
! Provides routines for manipulating mappings in a nontrivial manner.
!

module obj_mapping

        use env_kindtypes, only: si, dp
        use obj_trajectory,only: labeledTraj

        implicit none

        private

        public          :: mapping_postProc_config, postProcess_mapping

        type            :: mapping_postProc_config
                !subdivide each bead 
                logical           :: charge_subdivide
                logical,parameter :: translate_aminoacid_to_atom = .true.
        endtype

        type                                    :: map

                logical                  :: coherent = .true.

                !subdivide each bead 
                integer (si)             :: num_sites = 0

                integer (si),allocatable :: raw_mapping(:)

                contains
                        procedure :: getNSites          => get_num_sites
                        procedure :: getRawMapping      => get_raw_mapping
                        procedure :: setRawMapping      => map_set_raw_mapping

                        procedure :: isCoherent         => is_coherent
                        procedure :: setCoherent        => set_coherent
        endtype

        !extends map to have charge qualities
        type, extends(map)                  :: charge_map

                real    (dp),allocatable :: modified_charges(:)

                contains
                        procedure :: getMappedCharges   => get_mapped_charges

                        procedure :: setModifiedCharges => set_modified_charges
                        procedure :: getModifiedCharges => get_modified_charges
        endtype

        !this serves to represent a fully featured map.
        type, extends(charge_map)           :: attr_map
        endtype

        contains

                pure function get_num_sites(self)
                        implicit none

                        class(map),              intent(in   ) :: self
                        integer (si)            :: get_num_sites

                        get_num_sites = self%num_sites
                endfunction get_num_sites

                pure function get_raw_mapping(self,canonical)
                        implicit none

                        class(map),              intent(in   ) :: self
                        logical,   optional,     intent(in   ) :: canonical

                        integer (si),allocatable :: get_raw_mapping(:)

                        if (allocated(self%mapping)) then
                                get_raw_mapping = self%raw_mapping
                        else
                                allocate(get_raw_mapping(0))
                        endif

                endfunction get_raw_mapping

                pure subroutine map_set_raw_mapping(self,raw_mapping,nSites,&
                                                    assume_canonical,stay_coherent)
                        implicit none

                        class(map),           intent(in   ) :: self
                        integer (si),optional,intent(in   ) :: raw_mapping(:)
                        logical,     optional,intent(in   ) :: assume_canonical
                        logical,     optional,intent(in   ) :: stay_coherent

                        if (present(stay_coherent)) then
                                stay_coherent_ = stay_coherent
                        else
                                stay_coherent_ = .false.
                        endif

                        if (present(assume_canonical)) then
                                assume_canonical_ = assume_canonical
                        else
                                assume_canonical_ = .false.
                        endif

                endsubroutine map_set_raw_mapping

                pure function is_coherent(self)
                        implicit none

                        class(map),              intent(in   ) :: mapping

                        logical         :: is_coherent

                        is_coherent = 

                endfunction is_coherent

                !If this isn't overridden, no copying occurs!
                pure function postprocess_mapping(mapping,design,traj) result(processed_mapping)

                        implicit none

                        class(attr_map),              intent(in   ) :: mapping
                        type(mapping_postProc_config),intent(in   ) :: design
                        class(labeledTraj),           intent(in   ) :: traj

                        type(attr_map)                  :: processed_mapping
                        integer (si), allocatable       :: working_mapping(:), domain_charges(:)

                        if (design%translate_aminoacid_to_atom) then
                                working_mapping  = composeMapping(mapping,traj%site_res)
                                domain_charges   = traj%atomCharges
                        else
                                f
                        endif

                endsubroutine postProcess_mapping
endmodule obj_mappingProcess
