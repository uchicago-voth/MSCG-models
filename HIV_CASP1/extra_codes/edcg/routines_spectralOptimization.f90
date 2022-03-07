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
! Provides routines for optimizing a spectral based mapping.
! Used in spectral_models.
module routines_spectralOptimization

        use env_kindtypes, only: si, dp

        implicit none

        private

        public spectralOptimization, spectral_dist_config

        !for each metric considered, we have a few options which repeatedly come up.
        type                           :: spectral_dist_config

                !whether this metric is active in the calculation.
                logical                       :: active = .false.
                !how to modify the distance to create and affinity value, if needed.
                character(len=:),allocatable  :: dist_conv_type
                !how to parameterize distance to aff conversion.
                real    (dp)                  :: dist_conv_param
                !how to do preprocess k Nearest Neighbor (kNN) analysis
                character(len=:),allocatable  :: kNN_type
                !how many kNN to keep.
                integer (si)                  :: kNN

        end type spectral_dist_config

        contains

                !performs spectral clustering. 
                subroutine spectralOptimization(mapping,projection,eigenvalues,nCG,tdata,regularization,&
                                                kmeans_max_iter,kmeans_max_try,edcg_config,spatial_config,&
                                                pairVar_config,charge_config)

                        use env_kindtypes,              only: si,dp
                        use abs_obj_logger_real,        only: logger_real
                        use obj_ll,                     only: i_sp_dll
                        use core_graph,                 only: dist_to_lapl, sget_spectral_projection
                        use core_kmeans,                only: kmeans
                        use obj_trajectory,             only: traj_training_data
                        use core_sort,                  only: qsort
                        use core_convert,               only: itoa

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),            intent(inout),allocatable :: mapping(:)
                        real    (dp),            intent(inout),allocatable :: projection(:,:)
                        real    (dp),            intent(inout),allocatable :: eigenvalues(:)
                        integer (si),            intent(in   )             :: nCG
                        type(traj_training_data),intent(in   )             :: tdata
                        character(*),            intent(in   )             :: regularization
                        integer (si),            intent(in   )             :: kmeans_max_iter
                        integer (si),            intent(in   )             :: kmeans_max_try
                        !these give configuration on whether to include specific metrics, and
                        !how to process them.
                        type(spectral_dist_config),intent(in   )       :: edcg_config
                        type(spectral_dist_config),intent(in   )       :: spatial_config
                        type(spectral_dist_config),intent(in   )       :: pairVar_config
                        type(spectral_dist_config),intent(in   )       :: charge_config
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !buffer length to hold types (needs length since it is an array of strings).
                        integer (si),parameter                        :: string_buffer_len = 128
                        !index for transferring active metrics to clustering routines.
                        integer (si)                                  :: fill_index
                        !number of clustering distance metrics in use
                        integer (si)                                  :: num_active_metrics
                        real    (dp),          allocatable            :: dists(:,:,:)
                        !distance conversion types (e.g. exp)
                        character(len=string_buffer_len),allocatable  :: dist_convs(:)
                        !distance conversion parameters
                        real    (dp),          allocatable            :: dist_conv_params(:)
                        !k Nearest Neighbor types
                        character(len=string_buffer_len),allocatable  :: kNN_types(:)
                        !number of K nearest neighbors
                        integer (si),          allocatable            :: kNNs(:)

                        !get the count of the number of metrics we need to handle.
                        num_active_metrics = count( [ edcg_config%active,    spatial_config%active , &
                                                      pairVar_config%active, charge_config%active   ] )

                        if (num_active_metrics == 0) then
                                print*, "ERROR: No metrics were active in spectral clustering. Exiting."
                                stop
                        endif

                        !allocate the variables to send to underlying graph routines-- those routines do 
                        !_not_ know about e.g. 'edcg' vs 'spatial'.

                        allocate(dists(tdata%trj%nAtoms,tdata%trj%nAtoms,num_active_metrics))

                        allocate(dist_convs(num_active_metrics))
                        allocate(dist_conv_params(num_active_metrics))

                        allocate(kNN_types(num_active_metrics) )
                        allocate(kNNs(num_active_metrics))

                        !fill these variables
                        fill_index = 1

                        if (edcg_config%active) then
                                dists(:,:,fill_index)        = tdata%stats%getCustomDistMat()

                                dist_convs(fill_index)       = edcg_config%dist_conv_type
                                dist_conv_params(fill_index) = edcg_config%dist_conv_param
                                kNN_types(fill_index)        = edcg_config%kNN_type
                                kNNs(fill_index)             = edcg_config%kNN

                                fill_index = fill_index + 1
                        endif

                        !NOTE THAT VALUES ARE SOMETIMES ADJUSTED AS THEY ARE PASSED (e.g. diameter)
                        if (spatial_config%active) then
                                dists(:,:,fill_index)        = tdata%stats%getAvgSDistanceMat()

                                dist_convs(fill_index)  = spatial_config%dist_conv_type
                                dist_conv_params(fill_index) = spatial_config%dist_conv_param*&
                                                                        tdata%stats%getAvgSDiameter()
                                kNN_types(fill_index)        = spatial_config%kNN_type
                                kNNs(fill_index)             = spatial_config%kNN

                                fill_index = fill_index + 1
                        endif

                        if (pairVar_config%active) then
                                dists(:,:,fill_index)        = tdata%stats%getAvgSDistanceVarMat()

                                dist_convs(fill_index)       = pairVar_config%dist_conv_type
                                dist_conv_params(fill_index) = pairVar_config%dist_conv_param
                                kNN_types(fill_index)        = pairVar_config%kNN_type
                                kNNs(fill_index)             = pairVar_config%kNN

                                fill_index = fill_index + 1
                        endif

                        if (charge_config%active) then
                                dists(:,:,fill_index)        = tdata%stats%getPowChargeDiffMat()

                                dist_convs(fill_index)       = charge_config%dist_conv_type
                                dist_conv_params(fill_index) = charge_config%dist_conv_param
                                kNN_types(fill_index)        = charge_config%kNN_type
                                kNNs(fill_index)             = charge_config%kNN

                                !not needed, as it's the last config given.
                                !fill_index = fill_index + 1
                        endif

                        !we divide along the regularization cases.
                        select case (regularization)
                        case ("unregularized")
                                projection_compute: block
                                        real    (dp),allocatable        :: laplacian(:,:)

                                        laplacian = dist_to_lapl(distance_mats = dists,&
                                                                    dist_convs = dist_convs,&
                                                                   conv_params = dist_conv_params,&
                                                                          kNNs = kNNs,&
                                                                     kNN_types = kNN_types,&
                                                                     normalize = '')

                                        call sget_spectral_projection(projection,eigenvalues,laplacian,nCG)
                                endblock projection_compute
                        case ("sym")
                                sym_projection_compute: block
                                        real    (dp),allocatable        :: sym_laplacian(:,:)
                                        real    (dp)                    :: lapl_surplus

                                        sym_laplacian = dist_to_lapl(distance_mats = dists,&
                                                                        dist_convs = dist_convs,&
                                                                       conv_params = dist_conv_params,&
                                                                              kNNs = kNNs,&
                                                                         kNN_types = kNN_types,&
                                                                         normalize = 'sym')

                                        lapl_surplus = abs(sum(sym_laplacian) - real(size(sym_laplacian,1),dp))
                                        if (lapl_surplus < 1) then
                                                print*,"WARNING: The computed laplacian appears to have very small off &
                                                       &diagonal elements. This can crash lapack routines; check your conversion of &
                                                       &distances to affinities. Continuing."
                                                print*, "INFO: laplacian: ", sym_laplacian
                                        endif

                                        call sget_spectral_projection(projection,eigenvalues,sym_laplacian,nCG,.true.)
                                endblock sym_projection_compute
                        case ("rw")
                                rw_projection_compute: block
                                        real    (dp),allocatable        :: rw_laplacian(:,:)
                                        real    (dp)                    :: lapl_surplus

                                        rw_laplacian = dist_to_lapl(distance_mats = dists,&
                                                                       dist_convs = dist_convs,&
                                                                      conv_params = dist_conv_params,&
                                                                             kNNs = kNNs,&
                                                                        kNN_types = kNN_types,&
                                                                        normalize = 'rw')

                                        lapl_surplus = abs(sum(rw_laplacian) - real(size(rw_laplacian,1),dp))
                                        if (lapl_surplus < 1) then
                                                print*,"WARNING: The computed laplacian appears to have very small off &
                                                       &diagonal elements. This can crash lapack routines; check your conversion of &
                                                       &distances to affinities. Continuing."
                                                print*, "INFO: laplacian: ", rw_laplacian
                                        endif

                                        call sget_spectral_projection(  projection = projection,&
                                                                       eigenvalues = eigenvalues,&
                                                                         laplacian = rw_laplacian,&
                                                                      n_dimensions = nCG, &
                                                                 degree_regularize = .true.)
                                endblock rw_projection_compute
                        case default
                                print*, "ERROR: The passed spectral regularization parameter &
                                        &("//regularization//") wasn't understood. &
                                        &Check your command line arguments. &
                                        &'unregularized' and 'sym' are the currently accepted options."
                                STOP
                        end select

                        discretize_compute: block
                                real    (dp)    :: kmeans_sites(nCG,size(projection,2))

                                call kmeans(       labels = mapping,&
                                            training_data = projection,&
                                                centroids = kmeans_sites,&
                                                 max_iter = kmeans_max_iter,&
                                                 n_trials = kmeans_max_try,&
                                                initialize = .true.)
                        endblock discretize_compute

                        !make sure we have the expected number of sites. Else, future computations fail.
                        check_n_sites: block
                                integer (si),allocatable   :: mapping_(:)
                                integer (si)               :: iter, nSites

                                mapping_ = mapping
                                call qsort(mapping_)

                                nSites = 1
                                do iter=1,(size(mapping_)-1) 
                                        if (mapping_(iter) /= mapping_(iter+1)) then
                                                nSites = nSites + 1
                                        endif
                                enddo

                                if (nSites /= nCG) then
                                        print*, "ERROR: The number of sites in the derived model is "//&
                                                &itoa(nSites)//" but the model has "//itoa(nCG)//" sites."
                                        print*, "ERROR:" ,Mapping
                                        STOP
                                endif
                                        
                        endblock check_n_sites

                endsubroutine spectralOptimization

endmodule routines_spectralOptimization
