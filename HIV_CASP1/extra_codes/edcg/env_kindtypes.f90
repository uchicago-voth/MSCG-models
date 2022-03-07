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
! this module stores the kind definitions used.
!

module env_kindtypes

        !Single and double precision floats definined with
        !decimal precision statements.

        !integer, parameter:: dp=kind(0.d0)
        integer, parameter  :: sp = selected_real_kind(6)
        integer, parameter  :: dp = selected_real_kind(15)
        integer, parameter  :: qp = selected_real_kind(32)

        !Single and double precision integers defined via
        !ranges of fitable values.
        integer, parameter  :: si = selected_int_kind(7)
        integer, parameter  :: di = selected_int_kind(12)

        !Number of bytes each float occupies (needed for
        !endian conversion
        integer, parameter   :: c_sp = storage_size(1.0_sp)
        integer, parameter   :: c_dp = storage_size(1.0_dp)

        !Number of bytes each float occupies (needed for
        !endian conversion
        integer, parameter   :: c_si = (storage_size(1_si)/8)
        integer, parameter   :: c_di = (storage_size(1_di)/8)

        !Single byte integer
        integer, parameter   :: singleb = selected_int_kind(1)


        !!! NOTE !!!

        !these values are hypothetically used when the choice of precision is fixed
        !by external constraint. Their continued usage is strong discouraged, as
        !this code is no longer designed to allow for reparameterization of e.g. dp to sp.

        !Some things just can't change when switchign from double to single prec.
        !these are marked with _x.
        integer, parameter  :: sp_x = selected_real_kind(6)
        integer, parameter  :: dp_x = selected_real_kind(15)

        !Single and double precision integers defined via
        !ranges of fitable values.
        integer, parameter  :: si_x = selected_int_kind(7)
        integer, parameter  :: di_x = selected_int_kind(12)


endmodule env_kindtypes
