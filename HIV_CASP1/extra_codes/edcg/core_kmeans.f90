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
! This module provides geneal routines for dealing with kmeans problems.

module core_kmeans

        implicit none

        private

        public kmeans, centroid_random_change, random_centroid_init, sget_labels

        contains
                ! naively associates points with nearest centroid. O(n2) cost.
                ! @points:    points to be associated with cluster
                ! @centroids: points to be associated with cluster
                !
                ! returns classes: array of integer assignments of each cluster.
                pure subroutine sget_labels(labels, points, centroids, buffer)

                        use env_kindtypes,      only: si, dp
                        use routines_math,      only: computeAsymmetricDistanceMatrix
                        use core_stat,          only: which_min

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )             :: points(:,:)
                        real    (dp),intent(in   )             :: centroids(:,:)
                        real    (dp),intent(inout),allocatable :: buffer(:,:)
                        integer (si),intent(inout),allocatable :: labels(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)            :: nCentroids, nPoints, iter

                        nPoints    = size(points,1)
                        nCentroids = size(centroids,1)

                        !verify/correct the labels dimensions
                        if (allocated(labels)) then
                                if (size(labels)  /= nPoints ) then
                                        deallocate(labels)
                                        allocate(labels(nPoints))
                                endif
                        else
                                allocate(labels(nPoints))
                        endif

                        !verify/correct the buffer dimensions
                        if (allocated(buffer)) then

                                if (.not. (nPoints    == size(buffer,1) .and. &
                                           nCentroids == size(buffer,2))) then
                                        deallocate(buffer)
                                        allocate(buffer(nPoints,nCentroids))
                                endif
                        else
                                allocate(buffer(nPoints,nCentroids))
                        endif

                        call computeAsymmetricDistanceMatrix(points,centroids,buffer)

                        !This isn't ideal; we're going to thrash the cache a bit, but 
                        !either this or the previous loop has to do so.
                        do iter=1,size(buffer,1)
                                labels(iter) = which_min(buffer(iter,:))
                        enddo

                endsubroutine sget_labels

                ! derives centroids of a labeled data set.
                ! @points:  point positions
                ! @labels:  labels of points
                ! @nLabels: max number of unique labels
                !
                ! returns centroids: real 2 array of centroid positions
                pure function get_centroids(points, labels, nLabels) result(centroids)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )             :: points(:,:)
                        integer (si),intent(in   )             :: labels(:)
                        integer (si),intent(in   )             :: nLabels
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return
                        real    (dp)            :: centroids(nLabels,size(points,2))

                        integer (si)            :: label_count(nLabels), curr_label
                        integer (si)            :: id_iter, dimension_iter

                        centroids   = 0
                        label_count = 0

                        do dimension_iter=1,size(points,2)
                        do id_iter=1,size(points,1)
                                curr_label = labels(id_iter)

                                label_count(curr_label) = label_count(curr_label) + 1
                                centroids(curr_label,dimension_iter) =   centroids(curr_label,dimension_iter)&
                                                                       + points(id_iter,dimension_iter)
                        enddo
                        enddo

                        do dimension_iter=1,size(points,2)
                        do curr_label=1,nLabels
                                centroids(curr_label,dimension_iter) =   centroids(curr_label,dimension_iter)&
                                                                       / max(label_count(curr_label)/size(points,2),1)
                        enddo
                        enddo

                endfunction get_centroids

                ! optimizes clusters according to lloyds algorithm, with special termination conditions
                ! as specified by zhang et al for spatial EDCG. 
                ! !!! centroids are input; the labels are for output only (inout retains allocations).
                !
                ! @labels:              resulting labels from clustering.
                ! @training_data:       2 array of data to optimize on.
                ! @centroids:           2 array of inital centroid positions
                ! @max_iter (optional): bound on number of iterations to complete
                ! @intialize(optional): whether to intialize the centroid locations.
                !
                ! modifies centroids in place.
                subroutine kmeans(labels,training_data,centroids,distortion_tol,improvement_tol,&
                                  n_trials,max_iter,initialize)

                        use env_kindtypes,      only: si, dp
                        use core_stat,          only: which_min

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout),allocatable :: labels(:)
                        real    (dp),intent(in   )             :: training_data(:,:)
                        real    (dp),intent(inout)             :: centroids(:,:)
                        real    (dp),intent(in   ),optional    :: distortion_tol
                        real    (dp),intent(in   ),optional    :: improvement_tol
                        integer (si),intent(in   ),optional    :: n_trials
                        integer (si),intent(in   ),optional    :: max_iter
                        logical,     intent(in   ),optional    :: initialize
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)             :: max_iter_, n_trials_
                        integer (si)             :: iter, trial_iter, nLabels

                        integer (si),allocatable :: trial_labels(:)
                        real    (dp),allocatable :: dbuffer(:,:), trial_centroids(:,:)
                        real    (dp)             :: drift, old_drift, distortion_tol_
                        real    (dp)             :: residual, old_residual, trial_residual
                        real    (dp)             :: improvement_tol_
                        logical                  :: initialize_

                        if (present(initialize)) then
                                initialize_ = initialize
                        else
                                initialize_ = .true.
                        endif

                        if (initialize_) then
                                if (present(n_trials)) then
                                        if (n_trials > 0) then
                                                n_trials_ = n_trials
                                        else
                                                n_trials_ = 40
                                        endif
                                else
                                        n_trials_ = 40
                                endif
                        else
                                n_trials_ = 1
                        endif

                        if (present(max_iter)) then
                                max_iter_ = max_iter
                        else
                                max_iter_ = huge(max_iter)
                        endif

                        if (present(distortion_tol)) then
                                distortion_tol_ = distortion_tol
                        else
                                distortion_tol_ = real(10000,dp)*tiny(distortion_tol)
                        endif

                        if (present(improvement_tol)) then
                                improvement_tol_ = improvement_tol
                        else
                                improvement_tol_ = real(1000,dp)*tiny(improvement_tol)
                        endif

                        nLabels = size(centroids,1)
                        allocate(trial_labels(size(training_data,1)))

                        !!!! prep for kmeans !!!!!!!!!!!!!!!!
                        if (initialize_) then
                                trial_centroids = random_centroid_init(training_data,size(centroids,1))
                        else
                                trial_centroids = centroids
                        endif

                        call sget_labels(trial_labels,&
                                         training_data, &
                                         trial_centroids, &
                                         dbuffer)

                        drift = get_cluster_drift(training_data,&
                                                  trial_labels,&
                                                  nLabels,&
                                                  trial_centroids,&
                                                  .true.) 

                        trial_residual = get_kmeans_residual(training_data,&
                                                       trial_labels,&
                                                       nLabels,&
                                                       trial_centroids)

                        residual  = trial_residual
                        centroids = trial_centroids
                        labels    = trial_labels

                        old_drift = -1
                        convergence_trials: do trial_iter=1, n_trials_

                                if (initialize_) then
                                        trial_centroids = random_centroid_init(training_data,size(centroids,1))
                                else
                                        trial_centroids = centroids
                                endif

                                iter = 1

                                optimize: do

                                        !get assignments.
                                        call sget_labels(trial_labels,&
                                                         training_data, &
                                                         trial_centroids, &
                                                         dbuffer)

                                        drift  = get_cluster_drift(training_data,&
                                                                   trial_labels,&
                                                                   nLabels,&
                                                                   trial_centroids,&
                                                                   .true.) 


                                        !refresh the centroid positions.
                                        trial_centroids = get_centroids(training_data, trial_labels, nLabels)


                                        !get stats
                                        old_residual = trial_residual
                                        trial_residual = get_kmeans_residual(training_data,&
                                                                       trial_labels,&
                                                                       nLabels,&
                                                                       trial_centroids)

                                        iter = iter + 1

                                        !exit conditions
                                        if ((                        iter > max_iter_        ) .or.&
                                            (                       drift < distortion_tol_  ) .or.&
                                            (abs(trial_residual - old_residual) < improvement_tol_)) then
                                                !print*, "iter", iter, max_iter_       
                                                !print*, "drift", drift, distortion_tol_ 
                                                !print*, "residual", abs(trial_residual - old_residual), improvement_tol_
                                                !print*, "Exiting kmeans loop on iter ", iter
                                                exit optimize
                                        endif
                                enddo optimize
                                if (residual > trial_residual) then
                                        centroids = trial_centroids
                                        labels    = trial_labels
                                        residual  = trial_residual
                                        !print*, "Accepting new kmeans state", trial_iter
                                endif
                        enddo convergence_trials

                endsubroutine kmeans

                ! intializes cluster locations by assigning them to random data points.
                !  impure due to stochastic choice.
                ! @points : data points from which to choose locations (npoints, dimension)
                ! @numCluster: number of clusters to choose.
                !
                ! returns the list of centroids.
                function random_centroid_init(points, numClust) result(centroids)

                        use env_kindtypes,      only: si, dp
                        use core_random,        only: gen_rand_ordered_seq

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )          :: points(:,:)
                        integer (si),intent(in   )          :: numClust
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp),allocatable    :: centroids(:,:)
                        integer (si),allocatable    :: centroid_indices(:)
                        integer (si)                :: iter

                        allocate(centroids(numClust,size(points,2)))

                        centroid_indices = gen_rand_ordered_seq(numClust,1,size(points,1))

                        do iter=1,size(centroid_indices)
                                centroids(iter,:) = points(centroid_indices(iter),:)
                        enddo

                endfunction random_centroid_init

                ! randomly updates in place one of the centroid locations to a new data point.
                !  checks to make sure we aren't overlapping with another cluster.
                ! @centroids: array to choose centroid and update.
                ! @candidate_points: possible new positions for the centroid
                ! @max_tries: maximum number of tries to make the move
                ! @tol      : the taxicab distance used to tell if two points are too close
                ! @status   : exit status. 0 on success, >0 on failure to make a move.
                subroutine centroid_random_change(centroids, index_to_change, candidate_points, max_tries, tol, change_status)
                        use env_kindtypes,      only: si, dp
                        use core_random,        only: gen_rand_int_omit, gen_rand_int
                        use core_stat,          only: seq, seq_omit

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: centroids(:,:)
                        integer (si),intent(in   )          :: index_to_change
                        real    (dp),intent(in   )          :: candidate_points(:,:)
                        integer (si),intent(in   )          :: max_tries
                        real    (dp),intent(in   ),optional :: tol
                        integer (si),intent(  out),optional :: change_status
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                :: iter
                        integer (si)                :: try
                        logical                     :: success
                        real    (dp)                :: tol_, residual
                        real    (dp),allocatable    :: candidate_centroid(:)
                        integer (si)                :: candidate_index, centroid_index
                        integer (si),allocatable    :: unchanged_centroid_indices(:)

                        if (present(tol)) then
                                tol_ = tol
                        else
                                tol_ = .01
                        endif

                        !choose index of centroid we want to change
                        !index_to_change = gen_rand_int(lower=1,upper=size(centroids,1))

                        !generate an array of the unchanged centroids for later collision comparison
                        allocate(unchanged_centroid_indices,&
                                 source=seq_omit(lower=1,&
                                                 upper=size(centroids,1),&
                                                 omit=[ index_to_change ]))

                        try = 0
                        success = .false.
                        !loop to try moves
                        do while (.not. success)

                                try = try + 1
                                if (try > max_tries) then
                                        print*, "Couldn't randomly select new k-means centroids."
                                        if (present(change_status)) change_status = 1
                                        stop
                                endif

                                candidate_index = gen_rand_int(lower=1,upper=size(candidate_points,1))

                                candidate_centroid = candidate_points(candidate_index,:)

                                success = .true.
                                do iter=1,size(unchanged_centroid_indices)
                                        centroid_index = unchanged_centroid_indices(iter)
                                        residual = sum(abs(centroids(centroid_index,:) - candidate_centroid))

                                        if (residual < tol_) then
                                                success = .false.
                                        endif
                                enddo
                        enddo

                        if (present(change_status)) change_status = 0

                        !success
                        centroids(index_to_change,:) = candidate_centroid

                endsubroutine centroid_random_change

                ! calculates how far the center of the clusters are from the recorded centroids (needed in 
                !  the traditional spatial EDCG method.
                ! @points   : training data use in analysis.
                ! @labels   : labels of training data
                ! @nlabels  : number of unique labels (i.e. clusters)
                ! @centroids: centroids to do the comparison of drift with.
                ! @taxicab  : whether to use taxicab distance or euclidean distance.
                !
                ! returns the drift value.
                pure function get_cluster_drift(points, labels, nLabels, centroids, taxicab) result(drift)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )          :: points(:,:), centroids(:,:)
                        integer (si),intent(in   )          :: labels(:)
                        integer (si),intent(in   )          :: nLabels
                        logical,     intent(in   ),optional :: taxicab
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)                :: drift

                        real    (dp)                :: normalization

                        real    (dp),allocatable    :: cluster_means(:,:)
                        integer (si)                :: iter

                        logical                     :: taxicab_

                        if (present(taxicab)) then
                                taxicab_ = taxicab
                        else
                                taxicab_ = .false.
                        endif

                        !calculate normalization constant
                        normalization = 3 * nLabels

                        !generate the means of the current partitions
                        allocate(cluster_means,source=get_centroids(points,labels,nLabels))

                        if (taxicab_) then
                                !taxicab distance case
                                drift = sum(abs(centroids - cluster_means)) /normalization
                        else
                                !euclidean distance case
                                cluster_means = (centroids - cluster_means)**2

                                drift = 0
                                do iter=1,size(centroids,1)
                                        drift = drift + sqrt(sum(cluster_means(iter,:)))
                                enddo

                                drift = drift / normalization
                        endif

                endfunction get_cluster_drift

                ! calculates in-cluster dispersion: the sum of squared distances of each point to its cluster's 
                ! centroid.
                ! @points   : training data use in analysis.
                ! @labels   : labels of training data
                ! @nlabels  : number of unique labels (i.e. clusters)
                ! @centroids: centroids
                !
                ! returns the dispersion
                pure function get_kmeans_residual(points, labels, nLabels, centroids) result(residual)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )          :: points(:,:), centroids(:,:)
                        integer (si),intent(in   )          :: labels(:)
                        integer (si),intent(in   )          :: nLabels
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)                :: residual, sq_dist

                        integer (si)                :: iter

                        residual = 0

                        do iter=1,size(labels)
                                sq_dist = sum((points(iter,:) - centroids(labels(iter),:))**2)
                                residual = residual + sq_dist
                        enddo

                endfunction get_kmeans_residual

                !! throwaway test function for kmeans testing.
                !function kmeans_test(training_data, numClust) result(centroids)

                !        use env_kindtypes,      only: si, dp

                !        implicit none

                !        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !        real    (dp),intent(in   )          :: training_data(:,:)
                !        integer (si),intent(in   )          :: numClust
                !        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !        real    (dp),              allocatable :: centroids(:,:)

                !        allocate(centroids(numClust,size(training_data,2)))

                !        centroids = 0
                !endfunction  kmeans_test

endmodule core_kmeans
