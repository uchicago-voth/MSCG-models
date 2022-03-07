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
! This module provides general routines for calculating multipole moments of
! three dimensional points.

module core_multipole

        use env_kindtypes, only: si, dp

        implicit none

        private

        public dipole, quadrupole, get_center_of_scalar_value, inertial_matrix

        contains

                !returns dipole of points and values, usually charges, with respect to a origin.
                ! @points: array of locations: (id, dimension)
                ! @values: array of values (usually charges).
                ! @origin: origin of coordinate system from which to calculate dipole.
                !               Takes precidence over origin_weights.
                ! @origin_weights: weights for points when calculating origin.
                ! @buffer: a buffer, which if not given, is allocated and destroyed internally
                !               (for performance). Dims must match @points
                !
                ! If .not. ( present(origin) .or. present(origin_weights) ) then the points are weighted
                ! equally to find the origin.
                !
                function dipole(points, values, origin, origin_weights, buffer) !result(dipole)

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )                      :: points(:,:)
                        real    (dp),intent(in   )                      :: values(:)
                        real    (dp),intent(in   ),            optional :: origin(:)
                        real    (dp),intent(in   ),            optional :: origin_weights(:)
                        real    (dp),intent(inout),allocatable,optional :: buffer(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable                :: dipole(:)

                        real    (dp),            parameter      :: tol = 1e-14 !tolerance for checking if
                                                                               !floats equal zero.
                        real    (dp),allocatable                :: buffer_(:,:), origin_(:)
                        logical                                 :: use_ext_buffer, zero_sum
                        real    (dp),allocatable                :: weights(:)

                        !we check to see if an external bufer is available, and use it if so.
                        use_ext_buffer = .false.
                        if (present(buffer)) then
                                if (allocated(buffer)) then
                                        use_ext_buffer = .true.
                                endif
                        endif
                        if (use_ext_buffer) then
                                call move_alloc(buffer,buffer_)
                        else
                                allocate(buffer_,mold=points)
                        endif

                        !if the monompole is 0, then we don't need to translate points.
                        zero_sum = .false.
                        if (abs(sum(values)) < tol) zero_sum = .true.

                        if (.not. zero_sum) then
                                if (present(origin)) then
                                        origin_ = origin
                                else if (present(origin_weights)) then
                                        weights = origin_weights/sum(origin_weights)
                                        origin_ = get_center_of_scalar_value(points,weights,buffer_)
                                else
                                        !get center of points if we need it.
                                        allocate(weights(size(values)))
                                        weights = 1.0_dp/size(values)
                                        origin_ = get_center_of_scalar_value(points,weights,buffer_)
                                endif
                        else
                                origin_ = [ 0.0_dp ] !so that we can do the simple conditional below.
                        endif

                        !check to see if origin is essentially zero
                        ! (this means we can avoid translating points.
                        if (zero_sum .or. sum(abs(origin_)) < tol) then
                                !in this case, the dipole is invariant to translation.
                                dipole = get_center_of_scalar_value(points,values)
                        else
                                !here, we must translate. We use a given buffer if possible.
                                dipole = get_center_of_scalar_value(points - spread(origin_,1,size(points,1)),&
                                                                    values,&
                                                                    buffer_)
                        endif

                endfunction dipole

                !returns quadrupole of points and values, usually charges, with respect to a origin.
                ! @points: array of locations: (id, dimension)
                ! @values: array of values (usually charges).
                ! @origin: origin of coordinate system from which to calculate quadrupole.
                !               Overides origin_weights.
                ! @origin_weights: weights used to calculate origin. Are normalized internally.
                ! @buffer: a buffer, which if not give, is allocated and destroyed internally
                !               (for performance). Dims must match @points
                !
                ! If .not. ( present(origin) .or. present(origin_weights) ) then the points are weighted
                ! equally to find the origin.
                !
                function quadrupole(points, values, origin, origin_weights, buffer)

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )                      :: points(:,:)
                        real    (dp),intent(in   )                      :: values(:)
                        real    (dp),intent(in   ),            optional :: origin(:)
                        real    (dp),intent(in   ),            optional :: origin_weights(:)
                        real    (dp),intent(inout),allocatable,optional :: buffer(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable                :: quadrupole(:,:)

                        real    (dp),            parameter      :: tol = 1e-14 !for checking if floats
                                                                               !are essentially zero
                        real    (dp),allocatable                :: buffer_(:,:), origin_(:)
                        logical                                 :: use_ext_buffer
                        real    (dp),allocatable                :: weights(:), prod(:,:)

                        integer (si)                            :: iter1

                        !check if an external buffer is given, and if so use it.
                        use_ext_buffer = .false.
                        if (present(buffer)) then
                                if (allocated(buffer)) then
                                        use_ext_buffer = .true.
                                endif
                        endif
                        if (use_ext_buffer) then
                                call move_alloc(buffer,buffer_)
                        else
                                allocate(buffer_,mold=points)
                        endif

                        !calculate the origin if needed
                        if (present(origin)) then
                                origin_ = origin
                        else if (present(origin_weights)) then
                                weights = origin_weights/sum(origin_weights)
                                origin_ = get_center_of_scalar_value(points,weights,buffer_)
                        else
                                !use center of points as origin if not given
                                allocate(weights(size(values)))
                                weights = 1.0_dp/size(values)
                                origin_ = get_center_of_scalar_value(points,weights,buffer_)
                        endif

                        allocate(quadrupole(size(points,2),size(points,2)))
                        quadrupole = 0

                        allocate(prod(size(points,2),size(points,2)))

                        !formulas for dipoles can be found:
                        !http://www.pa.msu.edu/~duxbury/courses/phy481/Fall2009/Lecture14.pdf

                        if (sum(abs(origin_)) < tol) then
                                do iter1=1,size(points,1)
                                        call sgen_outer_product(prod,points(iter1,:))
                                        quadrupole = quadrupole + values(iter1)*3*prod
                                        call add_diag(quadrupole,values(iter1)*sum(points(iter1,:)**2))
                                enddo
                        else
                                buffer_ = points - spread(origin_,1,size(points,1))
                                do iter1=1,size(buffer_,1)
                                        call sgen_outer_product(prod,buffer_(iter1,:))
                                        quadrupole = quadrupole + values(iter1)*3*prod
                                        call add_diag(quadrupole,-values(iter1)*sum(buffer_(iter1,:)**2))
                                enddo
                        endif

                        quadrupole = quadrupole * 0.5_dp

                        if (use_ext_buffer) call move_alloc(buffer,buffer_)

                endfunction quadrupole

                !returns moment of intertia of points and masses with respect to a origin.
                ! @points: array of locations: (id, dimension)
                ! @values: array of values (usually charges).
                ! @origin: origin of coordinate system from which to calculate quadrupole.
                ! @buffer: a buffer, which if not give, is allocated and destroyed internally
                !               (for performance). Dims must match @points
                function inertial_matrix(points, masses, origin, buffer) result(im)

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )                      :: points(:,:)
                        real    (dp),intent(in   )                      :: masses(:)
                        real    (dp),intent(in   ),            optional :: origin(:)
                        real    (dp),intent(inout),allocatable,optional :: buffer(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable                :: im(:,:)

                        real    (dp),            parameter      :: tol = 1e-14 !for checking if floats
                                                                               !are essentially zero
                        real    (dp),allocatable                :: buffer_(:,:)
                        real    (dp),allocatable                :: origin_(:)
                        logical                                 :: use_ext_buffer
                        real    (dp)                            :: xyz_sq(3)

                        !check if an external buffer is given, and if so use it.
                        use_ext_buffer = .false.
                        if (present(buffer)) then
                                if (allocated(buffer)) then
                                        use_ext_buffer = .true.
                                endif
                        endif
                        if (use_ext_buffer) then
                                call move_alloc(buffer,buffer_)
                        else
                                allocate(buffer_,mold=points)
                        endif

                        !calculate the origin if needed
                        if (present(origin)) then
                                origin_ = origin
                        else
                                !use center of points as origin if not given
                                origin_ = get_center_of_scalar_value(points,masses,buffer_)
                        endif

                        allocate(im(size(points,2),size(points,2)))

                        !formulas for dipoles can be found:
                        if (sum(abs(origin_)) < tol) then
                                xyz_sq = sum(points**2*&
                                                 spread(masses,2,size(points,2)),&
                                                                                   dim=2)
                                im(1,1) = xyz_sq(2) + xyz_sq(3)
                                im(2,2) = xyz_sq(1) + xyz_sq(3)
                                im(3,3) = xyz_sq(1) + xyz_sq(2)

                                im(1,2) = - sum(points(:,1)*points(:,2)*masses)
                                im(2,1) = im(1,2)

                                im(1,3) = - sum(points(:,1)*points(:,3)*masses)
                                im(3,1) = im(1,3)

                                im(2,3) = - sum(points(:,2)*points(:,3)*masses)
                                im(3,2) = im(2,3)
                        else
                                buffer_ = points - spread(origin_,1,size(points,1))
                                xyz_sq = sum(buffer_**2*&
                                                spread(masses,2,size(buffer_,2)),&
                                                                                   dim=2)
                                im(1,1) = xyz_sq(2) + xyz_sq(3)
                                im(2,2) = xyz_sq(1) + xyz_sq(3)
                                im(3,3) = xyz_sq(1) + xyz_sq(2)

                                im(1,2) = - sum(buffer_(:,1)*buffer_(:,2))
                                im(2,1) = im(1,2)

                                im(1,3) = - sum(buffer_(:,1)*buffer_(:,3))
                                im(3,1) = im(1,3)

                                im(2,3) = - sum(buffer_(:,2)*buffer_(:,3))
                                im(3,2) = im(2,3)
                        endif

                endfunction inertial_matrix

                !computes a weighted average per dimension, WITHOUT normalizing weights.
                ! @points: points to perform average on. (id,dimension)
                ! @values: weights to use in weighted average.
                ! @buffer: buffer to accelerate execution (prior allocation). 
                !               dim must match (points)
                ! @collapse_dim: which dimension to perform the average along (previously 
                !               in these argument desriptions it was assumed 1 for id.)
                function get_center_of_scalar_value(points,values,buffer,collapse_dim) result(center)

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )                      :: points(:,:)
                        real    (dp),intent(in   )                      :: values(:)
                        real    (dp),intent(inout),allocatable,optional :: buffer(:,:)
                        integer (si),intent(in   ),            optional :: collapse_dim
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable        :: center(:)

                        real    (dp),allocatable        :: buffer_(:,:)
                        integer (si)                    :: collapse_dim_
                        logical                         :: use_ext_collapse_dim, use_ext_buffer
                        integer (si)                    :: iter1

                        !optional argument parsing. 
                        use_ext_collapse_dim = .false.
                        if (present(collapse_dim)) then
                                if (collapse_dim < 3 .and. collapse_dim > 0) then
                                        use_ext_collapse_dim = .true.
                                endif
                        endif

                        if (use_ext_collapse_dim) then
                                collapse_dim_ = collapse_dim
                        else
                                collapse_dim_ = 1
                        endif

                        !if we can use the external buffer, use it.

                        use_ext_buffer = .false.
                        if (present(buffer)) then
                                if (allocated(buffer)) then
                                        use_ext_buffer = .true.
                                endif
                        endif

                        if (use_ext_buffer) then
                                call move_alloc(buffer,buffer_)
                        else
                                allocate(buffer_,mold=points)
                        endif

                        allocate(center(size(buffer_,2)))

                        if (collapse_dim_ == 1) then
                                !do weighted mean, depending on which dimensions we should collapse on.
                                do iter1=1,size(buffer_,1)
                                        buffer_(iter1,:) = points(iter1,:) * values(iter1)
                                enddo

                                center = sum(buffer_,1)
                        else
                                !do weighted mean, depending on which dimensions we should collapse on.
                                do iter1=1,size(buffer_,2)
                                        buffer_(:,iter1) = points(:,iter1) * values(iter1)
                                enddo

                                center = sum(buffer_,2)
                        endif

                        !return buffer to the original handle if needed

                        if (use_ext_buffer) call move_alloc(buffer_,buffer)

               endfunction get_center_of_scalar_value

               !add a value to the diagonal of a matrix
               pure subroutine add_diag(matrix,toadd)

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)      :: matrix(:,:)
                        real    (dp),intent(in   )      :: toadd
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                    :: iter1

                        do iter1=1,min(size(matrix,1),size(matrix,2))
                                matrix(iter1,iter1) = matrix(iter1,iter1) + toadd
                        enddo

                endsubroutine add_diag

                !give the outer product of a vector with itself.
               pure subroutine sgen_outer_product(prod,vec)

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(inout)      :: prod(:,:)
                        real    (dp),intent(in   )      :: vec(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)                    :: iter1,iter2

                        prod = 0
                        do iter1=1,size(vec)
                        do iter2=1,size(vec)
                                prod(iter1,iter2) = vec(iter1) * vec(iter2)
                        enddo
                        enddo

                endsubroutine sgen_outer_product

endmodule core_multipole
