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
! This module provides geneal routines for dealing with spectal graph problems.

module core_graph

        implicit none

        private

        public dist_to_lapl, sget_spectral_projection

        contains
                !performs a spectral projection of the data. Returned eigenvectors (which correspond to the project
                !coordinates) and eigenvalues are in ascending order, skipping the first constant eigenvalue.
                subroutine sget_spectral_projection(projection,eigenvalues,laplacian,n_dimensions,row_normalize,degree_regularize)

                        use env_kindtypes,      only: si, dp
                        use core_matrix,        only: spectrum

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(  out),allocatable                 :: projection(:,:)
                        real    (dp),intent(  out),allocatable,optional        :: eigenvalues(:)
                        real    (dp),intent(in   )                             :: laplacian(:,:)
                        integer (si),intent(in   )                             :: n_dimensions
                        logical,     intent(in   )            ,optional        :: row_normalize
                        logical,     intent(in   )            ,optional        :: degree_regularize
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable                               :: eigenvalues_(:)
                        real    (dp),allocatable                               :: eigenvectors(:,:)
                        integer (si)                                           :: n_dimensions_
                        integer (si)                                           :: iter
                        logical                                                :: row_normalize_, degree_regularize_
                        real    (dp)                                           :: norm

                        if (present(row_normalize)) then
                                row_normalize_ = row_normalize
                        else
                                row_normalize_ = .false.
                        endif

                        !use guarded var: n_dimensions_
                        if (n_dimensions == -1 .or. n_dimensions > size(laplacian,1) - 1) then
                                n_dimensions_ = size(laplacian,1) -1
                        else
                                n_dimensions_ = n_dimensions
                        endif

                        if (present(degree_regularize)) then
                                degree_regularize_ = degree_regularize
                        else
                                degree_regularize_ = .false.
                        endif

                        !do the spectral decomposition
                        if (degree_regularize_) then
                                spectrum_work: block
                                        integer (si)             :: iter
                                        real    (dp),allocatable :: degree_matrix(:,:)

                                        allocate(degree_matrix(size(laplacian,1),size(laplacian,2)))
                                        degree_matrix = 0
                                        do iter=1,size(laplacian,1)
                                                degree_matrix(iter,iter) = laplacian(iter,iter)
                                        enddo

                                        call spectrum(       M = laplacian,&
                                                             B = degree_matrix,&
                                                             d = eigenvalues_,&
                                                             Q = eigenvectors,&
                                        num_bottom_eigenvalues = n_dimensions + 1)
                                endblock spectrum_work
                        else
                                call spectrum(       M = laplacian,&
                                                     d = eigenvalues_,&
                                                     Q = eigenvectors,&
                                num_bottom_eigenvalues = n_dimensions + 1)
                        endif

                        !filter and reverse eigenvalues
                        allocate(projection(size(laplacian,1),n_dimensions_))
                        do iter=1,n_dimensions_
                                projection(:,iter) = eigenvectors(:,iter+1)
                        enddo

                        if (row_normalize_) then
                                do iter=1,size(projection,1)
                                        norm = sqrt(sum(projection(iter,:)**2))
                                        if (norm > 10*tiny(norm)) then
                                                projection(iter,:) = projection(iter,:) / norm
                                        else
                                                projection(iter,:) = 1.0_dp/size(projection(iter,:))
                                        endif
                                enddo
                        endif

                        !transfer allocation if eigenvalues arg present
                        if (present(eigenvalues)) then
                                allocate(eigenvalues(n_dimensions_ + 1))
                                do iter=1,n_dimensions_+1
                                        eigenvalues(iter) = eigenvalues_(iter)
                                enddo
                        endif

                endsubroutine sget_spectral_projection


                !Converts an array of distance matrices (as a 3 array) into a laplacian.
                !
                !Arguments are _arrays_ of conversion parameters (except normalize).
                !
                !@distance_mats(:,:,:)  :: Primary distance matrix between points. Required.
                !@normalize             :: Whether to normalize the laplacian: accepts 'sym',
                !                       ::     'lw', '' (unnormalized). Defaults to unnormalized.
                !@dist_convs(:)         :: How to convert between distance and affinity. Accepts
                !                       ::     'trans', 'exp' 
                !@conv_params(:)        :: When converting to affinity, what real number parameterizes the 
                !                       ::     conversion.
                !@kNNs(:)               :: Number of nearest neighbors for kmeans. 0 -> no filtering.
                !@kNN_types(:)          :: Method for kNN: accepts 'mutual', 'max'
                function dist_to_lapl(distance_mats,normalize,dist_convs,conv_params,&
                                      kNNs,kNN_types) result(lapl)

                        use env_kindtypes,      only: si,dp
                        use core_matrix,        only: quadForm

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )                       :: distance_mats(:,:,:)
                        real    (dp),intent(in   )                       :: conv_params(:)
                        character(*),intent(in   )                       :: dist_convs(:)
                        character(*),intent(in   )                       :: normalize
                        integer (si),intent(in   )                       :: kNNs(:)
                        character(*),intent(in   ),dimension(:)          :: kNN_types
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable        :: lapl(:,:)
                        real    (dp),allocatable        :: affs(:,:,:)
                        real    (dp),allocatable        :: mod_degree(:,:)
                        integer (si)                    :: iter, dist_to_aff_iter

                        !check to see that the sizes of parameters were valid. Exit if not.
                        if ( (size(kNNs)      /= size(distance_mats,3))  .or. &
                             (size(kNN_types) /= size(distance_mats,3))) then
                                allocate(lapl(0,0))
                                return
                        endif

                        allocate(affs,mold=distance_mats)

                        do dist_to_aff_iter=1,size(affs,3)

                                !negative sign makes results positive
                                affs(:,:,dist_to_aff_iter) = &
                                        - dist_to_unnorm_lapl(distance_mat = distance_mats(:,:,dist_to_aff_iter),&
                                                                 dist_conv = dist_convs(dist_to_aff_iter),&
                                                                conv_param = conv_params(dist_to_aff_iter),&
                                                                   no_diag = .true.)

                                !if present and nonblank, do kNN filtering.
                                if ((len_trim(kNN_types(dist_to_aff_iter)) > 0) .and. &
                                    (              kNNs(dist_to_aff_iter)  > 0)) then
                                        select case (kNN_types(dist_to_aff_iter))
                                                case ("mutual")
                                                        call aff_kNN_filter_mutual(affs(:,:,dist_to_aff_iter),&
                                                                                   kNNs(dist_to_aff_iter))
                                                case ("max")
                                                        call aff_kNN_filter_max(affs(:,:,dist_to_aff_iter),&
                                                                                kNNs(dist_to_aff_iter))
                                                case default
                                                        allocate(lapl(0,0))
                                                        return
                                        endselect
                                endif

                        enddo

                        !make the resulting laplacian.
                        lapl = - product(affs,dim=3)

                        call add_degree_entries(lapl,-1.0_dp)

                        !if no normalization is requested, we are done.
                        if (len_trim(normalize) == 0) then
                                return
                        endif

                        !If normalization is requested, we do it. Note that normalized 
                        !clustering doesn't always use these explicitly, so they aren't always used in practice.
                        select case (normalize)
                                case ("sym")
                                        allocate(mod_degree(size(lapl,1),size(lapl,1)))
                                        do iter=1,size(lapl,1)
                                                mod_degree(iter,:)    = 0
                                                mod_degree(iter,iter) = 1.0_dp/sqrt(lapl(iter,iter))
                                        enddo
                                        lapl = quadForm(lapl,mod_degree)
                                        return
                                case ("rw")
                                        allocate(mod_degree(size(lapl,1),size(lapl,1)))
                                        do iter=1,size(lapl,1)
                                                mod_degree(iter,:)    = 0
                                                mod_degree(iter,iter) = 1.0_dp/lapl(iter,iter)
                                        enddo
                                        lapl = matmul(mod_degree,lapl)
                                        return
                                case default
                                        if (allocated(lapl)) then
                                                deallocate(lapl)
                                        else
                                                allocate(lapl(0,0))
                                        endif
                                        return
                        end select

                endfunction dist_to_lapl

                !
                ! @distance_mat:
                !
                ! returns graph laplacian.
                pure function dist_to_unnorm_lapl(distance_mat,dist_conv,conv_param,&
                                                  no_diag) result(lapl)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )          :: distance_mat(:,:)
                        character(*),intent(in   )          :: dist_conv
                        real    (dp),intent(in   )          :: conv_param
                        logical,     intent(in   ),optional :: no_diag
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable        :: lapl(:,:)
                        real    (dp),allocatable        :: aff(:,:)
                        logical                         :: no_diag_

                        if (present(no_diag)) then
                                no_diag_ = no_diag
                        else
                                no_diag_ = .false.
                        endif

                        select case (dist_conv)
                                case ("trans")
                                        aff = dist_to_aff_translation(distance_mat,conv_param)
                                case ("exp")
                                        aff = dist_to_aff_expDecay(distance_mat,conv_param)
                                case default
                                        allocate(lapl(0,0))
                                        return
                        end select


                        lapl = - aff

                        if (no_diag_) then
                                call add_degree_entries(lapl,0.0_dp)
                        else
                                call add_degree_entries(lapl,-1.0_dp)
                        endif

                endfunction dist_to_unnorm_lapl

                pure function dist_to_aff_translation(distance_mat,offset) result(aff)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )      :: distance_mat(:,:)
                        real    (dp),intent(in   )      :: offset
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable        :: aff(:,:)

                        real    (dp)                    :: true_offset

                        true_offset = maxval(distance_mat) + offset

                        aff = true_offset - distance_mat

                endfunction dist_to_aff_translation

                pure function dist_to_aff_expDecay(distance_mat,normalization) result(aff)

                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )      :: distance_mat(:,:)
                        real    (dp),intent(in   )      :: normalization
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable        :: aff(:,:)

                        !explicity alloc is neccesary because of gfortran bug (generic function?)
                        allocate(aff,mold=distance_mat)

                        aff = exp( - distance_mat / normalization)

                endfunction dist_to_aff_expDecay

                pure subroutine add_degree_entries(graph_matrix,factor)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: graph_matrix(:,:)
                        real    (dp),intent(in   ),optional :: factor
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)    :: iter
                        real    (dp)    :: factor_

                        if (present(factor)) then
                                factor_ = factor
                        else
                                factor_ = 1
                        endif

                        do iter=1,size(graph_matrix,1)
                                graph_matrix(iter,iter) = 0
                                graph_matrix(iter,iter) = factor_ * real(sum(graph_matrix(iter,:)),dp)
                        enddo

                endsubroutine add_degree_entries

                pure subroutine aff_kNN_filter_max(aff,kNN)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: aff(:,:)
                        integer (si),intent(in   )          :: kNN
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                        :: iter

                        do iter=1,size(aff,1)
                                call top_n_filter(aff(:,1),kNN)
                        enddo

                        call symmeterize_matrix_max(aff)

                endsubroutine aff_kNN_filter_max

                !does kNN on a affinity matrix WITH 0 ON THE DIAGONAL.
                pure subroutine aff_kNN_filter_mutual(aff,kNN)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: aff(:,:)
                        integer (si),intent(in   )          :: kNN
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                        :: iter

                        do iter=1,size(aff,1)
                                call top_n_filter(aff(:,1),kNN)
                        enddo

                        call filter_matrix_and(aff)

                endsubroutine aff_kNN_filter_mutual

                !does kNN on a affinity matrix WITH 0 ON THE DIAGONAL.
                pure subroutine top_n_filter(array,N)

                        use env_kindtypes,      only: si, dp
                        use core_sort,          only: indexQsort

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: array(:)
                        integer (si),intent(in   )          :: N
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),parameter              :: replacement = 0.0_dp
                        integer (si),allocatable            :: sort_indices(:)
                        integer (si)                        :: iter

                        allocate(sort_indices(size(array,1)))

                        call indexQsort(array,sort_indices)

                        do iter=1,(size(array,1) - N)
                                array(sort_indices(iter)) = replacement
                        enddo

                endsubroutine top_n_filter

                pure subroutine filter_matrix_and(M)

                        use env_kindtypes,      only: si,dp

                        implicit none

                        real    (dp),intent(inout)    :: M(:,:)

                        integer (si)                  :: iter1, iter2

                        do iter1=1,size(M,1)-1
                        do iter2=iter1+1,size(M,1)
                                call real_pair_positive_and(M(iter1,iter2),M(iter2,iter1))
                        enddo
                        enddo

                endsubroutine filter_matrix_and

                pure subroutine symmeterize_matrix_max(M)

                        use env_kindtypes,      only: si,dp

                        implicit none

                        real    (dp),intent(inout)    :: M(:,:)

                        integer (si)                  :: iter1, iter2

                        do iter1=1,size(M,1)-1
                        do iter2=iter1+1,size(M,1)
                                call real_pair_positive_max(M(iter1,iter2),M(iter2,iter1))
                        enddo
                        enddo

                endsubroutine symmeterize_matrix_max

                pure subroutine real_pair_positive_and(r_1,r_2)

                        use env_kindtypes,      only: dp
                        use core_sort,          only: indexQsort

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: r_1,r_2
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        if (.not.((r_1 > 0.0_dp) .and. (r_2 > 0.0_dp))) then
                                r_1 = 0
                                r_2 = 0
                        endif

                endsubroutine real_pair_positive_and

                pure subroutine real_pair_positive_max(r_1,r_2)

                        use env_kindtypes,      only: dp
                        use core_sort,          only: indexQsort

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)          :: r_1,r_2
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        r_1 = max( r_1, r_2 )
                        r_2 = r_1

                endsubroutine real_pair_positive_max

endmodule core_graph
