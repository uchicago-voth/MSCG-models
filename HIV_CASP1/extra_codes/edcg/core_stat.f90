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
! This module contains general purpose statistics routines. They are pure when possible.
!

module core_stat

        implicit none

        private

        public var, sd, seq, seq_omit, buffer_seq, expDistribution, randomSampleNoReplace, &
               get_bca_accel, get_bca_bias, bca_transform_percentile, assignment_difference,&
               which_min, outlier_mask, rog, which

        contains
                !compute the naive sample variance of an array. stack allocated.
                !this should be changed to a stable variance formula.
                pure function var(arg)
                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp), intent(in   )     :: arg(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)    :: var

                        !local variables
                        real    (dp)    :: mean
                        real    (dp)    :: arg_(size(arg))

                        arg_ = arg

                        mean = sum(arg_)/size(arg_)

                        var = sum((arg_ - mean)**2) * (1/(size(arg_)-1))

                        return
                endfunction var

                !compute the sample standard deviation of an array. Pure when possible.
                pure function sd(arg)
                        use env_kindtypes,      only: dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp), intent(in   )     :: arg(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)    :: sd

                        !local sdiables
                        real    (dp)    :: mean
                        real    (dp)    :: arg_(size(arg))

                        arg_ = arg

                        mean = sum(arg_)/size(arg_)

                        sd = sqrt(sum((arg_ - mean)**2) * (1/(size(arg_)-1)))

                        return
                endfunction sd

                pure function which(array)
                        use env_kindtypes,      only: si

                        implicit none

                        logical,intent(in   )   :: array(:)

                        integer (si), allocatable :: which(:)

                        integer (si)              :: iter1, write_mark
                        integer (si)              :: true_count

                        true_count = 0

                        do iter1=1,size(array)
                                if (array(iter1)) true_count = true_count + 1
                        enddo

                        allocate(which(true_count))

                        write_mark = 1
                        do iter1=1,size(array)
                                if (array(iter1)) then
                                        which(write_mark) = iter1
                                        write_mark = write_mark + 1
                                endif
                        enddo

                endfunction which

                pure function which_min(x) result(which)
                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp), intent(in   )     :: x(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si)    :: which

                        integer (si)    :: iter
                        real    (dp)    :: best

                        which = 1
                        best  = x(1)

                        do iter=2,size(x)
                                if (x(iter) < best) then
                                        which = iter
                                        best  = x(iter)
                                endif
                        enddo

                endfunction which_min

                !return a array which contains a sequence of integers.
                pure function seq(lower,upper,stride)

                        use env_kindtypes,      only: si
                        use core_random,        only: gen_rand_int, gen_rand_int_omit

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: lower, upper
                        integer (si),intent(in   ),optional :: stride
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si),allocatable:: seq(:)

                        !local vars
                        integer (si)            :: stride_
                        integer (si)            :: upper_
                        integer (si)            :: index, i, remainder

                        if (present(stride)) then
                                if (stride >= 1) then
                                        stride_ = stride
                                else
                                        stride_ = 1
                                endif
                        else
                                stride_ = 1
                        endif

                        upper_ = upper

                        if (stride_ /= 1) then
                                remainder = mod(upper_-lower,stride_)

                                if (remainder > 0) then
                                        upper_ = upper_ - remainder
                                endif
                        endif

                        allocate(seq(((upper_-lower)/stride_)+1))

                        index = 1
                        do i=lower,upper_,stride_
                                seq(index) = i
                                index = index + 1
                        enddo

                        return
                endfunction seq

                !write a integer sequence inclusive from lower to upper in buffer.
                !marks last written location with endWrite. No checks are done.
                pure subroutine buffer_seq(buffer,endWrite,lower,upper)
                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(inout)          :: buffer(:)
                        integer (si),intent(  out)          :: endWrite
                        integer (si),intent(in   )          :: lower, upper
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !local vars
                        integer (si)            :: i, index

                        index = 1
                        do i=lower,upper
                                buffer(index) = i
                                index = index + 1
                        enddo

                        endWrite = upper - lower + 1
                endsubroutine buffer_seq

                !samples from a POSITIVE BOUNDED exponential distribution 
                function expDistribution(maxChoice,base) result(choice)
                        use env_kindtypes,      only: si
                        use core_random,        only: gen_rand_int

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si), intent(in   )     :: maxChoice,base
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si)    :: choice

                        !local vars
                        integer (si)    :: denom, choiceIndex, i

                        denom = 0
                        do i=1,maxChoice
                                denom = denom + base**i
                        enddo

                        choiceIndex = gen_rand_int(1,denom)

                        choice = 1
                        do i=maxChoice,1,-1
                                if (i <denom) then
                                        choice = i
                                        return
                                endif
                        enddo
                endfunction expDistribution

                !Samples nSamples from a given domain without replacement.
                function randomSampleNoReplace(nSamples,domain) result(sample)
                        use env_kindtypes,      only: si
                        use core_random,        only: gen_rand_int, gen_rand_int_omit

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )      :: domain(:)
                        integer (si),intent(in   )      :: nSamples
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value (STACK ALLOCATED)
                        integer (si)    :: sample(nSamples)

                        !local vars
                        integer (si)    :: omit_list(nSamples-1) !STACK ALLOCATED
                        integer (si)    :: index, i

                        if (nSamples == 1) then
                                index = gen_rand_int(size(domain,1),1)
                                sample(1) = domain(index)
                                omit_list(1) = index
                        else
                                do i=1,nSamples-1
                                        index = gen_rand_int_omit(size(domain,1),1,omit_list(:i-1))
                                        sample(i) = domain(index)
                                        omit_list(i) = index
                                enddo
                                index = gen_rand_int_omit(size(domain,1),1,omit_list(:i-1))
                                sample(i) = domain(index)
                        endif

                        return
                endfunction randomSampleNoReplace

                !return a sequence from lower to upper, inclusive, omitting numbers in omit.
                pure function seq_omit(lower,upper,omit) result(seq)
                        use env_kindtypes,      only: si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        integer (si),intent(in   )          :: lower, upper
                        integer (si),intent(in   )          :: omit(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        integer (si),allocatable:: seq(:)

                        !local vars
                        integer (si)            :: index, i

                        allocate(seq(upper-lower+1-size(omit)))

                        index=1
                        do i=lower,upper
                                if (.not. any(i == omit)) then
                                        seq(index) = i
                                        index = index + 1
                                endif
                        enddo

                endfunction seq_omit

                pure function get_bca_accel(jk_dist) result(accel)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        real    (dp),intent(in   ) :: jk_dist(:)

                        real    (dp)               :: accel

                        real    (dp)               :: jk_mean, denom
                        integer (si)               :: iter

                        jk_mean = sum(jk_dist)/max(size(jk_dist),1)

                        accel = 0

                        do iter=1,size(jk_dist)
                                accel = accel + (jk_mean - jk_dist(iter))**3
                        enddo 

                        denom = 0

                        do iter=1,size(jk_dist)
                                denom = denom + (jk_mean - jk_dist(iter))**2
                        enddo 

                        accel = accel / ( 6* denom **(3.0_dp/2.0_dp) )

                endfunction get_bca_accel

                !assumes that bs_dist is sorted.
                function get_bca_bias(bs_dist,estimate) result(bias)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )   :: bs_dist(:)
                        real    (dp),intent(in   )   :: estimate
                        !!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp)                    :: bias
                        real    (dp)                    :: population_prop

                        integer (si)                    :: iter

                        do iter=1,size(bs_dist)
                                if (bs_dist(iter) > estimate) then
                                        exit
                                endif
                        enddo 

                        population_prop = real(iter,dp) / max(size(bs_dist),1)

                        bias = normal_01_cdf_inv(population_prop)

                        !print*, "BIAS VALUES", population_prop, bias

                endfunction get_bca_bias

                !r8poly_value_horner evaluates a polynomial using Horner's method.
                !  Original Author: John Burkardt (2014)
                !  Modified for syntax, minor passing conventions, purity.
                !  @c(:):, @x:
                !    The polynomial 
                !
                !      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
                !
                !    is to be evaluated at the value @x.
                pure function r8poly_value_horner (c, x)

                        use env_kindtypes,      only: si, dp

                        implicit none
                        
                        !!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   ) :: c(0:)
                        real    (dp),intent(in   ) :: x
                        !!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)               :: r8poly_value_horner
                        
                        integer (si)               :: i
                        real    (dp)               :: scratch
                        
                        scratch = c(size(c)-1)
                        do i = (size(c)-2), 0, -1
                                scratch = scratch * x + c(i)
                        end do
                        
                        r8poly_value_horner = scratch
                endfunction

                !! normal_01_cdf computes the standard normal CDF.
                !  Original Author: John Burkardt (1999)
                !  Modified for syntax, minor passing conventions, purity.
                !  @x gives the domain value to calculate the probability for.
                !  Reference:
                !
                !    AG Adams,
                !    Algorithm 39,
                !    Areas Under the Normal Curve,
                !    Computer Journal,
                !    Volume 12, pages 197-198, 1969.
                !
                pure function normal_01_cdf(x) result(cdf)

                        use env_kindtypes,      only: dp

                        implicit none

                        real (dp), intent(in   ) :: x

                        !return value
                        real (dp)            :: cdf
                        
                        real (dp), parameter :: a1 = 0.398942280444D+00
                        real (dp), parameter :: a2 = 0.399903438504D+00
                        real (dp), parameter :: a3 = 5.75885480458D+00
                        real (dp), parameter :: a4 = 29.8213557808D+00
                        real (dp), parameter :: a5 = 2.62433121679D+00
                        real (dp), parameter :: a6 = 48.6959930692D+00
                        real (dp), parameter :: a7 = 5.92885724438D+00
                        real (dp), parameter :: b0 = 0.398942280385D+00
                        real (dp), parameter :: b1 = 3.8052D-08
                        real (dp), parameter :: b2 = 1.00000615302D+00
                        real (dp), parameter :: b3 = 3.98064794D-04
                        real (dp), parameter :: b4 = 1.98615381364D+00
                        real (dp), parameter :: b5 = 0.151679116635D+00
                        real (dp), parameter :: b6 = 5.29330324926D+00
                        real (dp), parameter :: b7 = 4.8385912808D+00
                        real (dp), parameter :: b8 = 15.1508972451D+00
                        real (dp), parameter :: b9 = 0.742380924027D+00
                        real (dp), parameter :: b10= 30.789933034D+00
                        real (dp), parameter :: b11= 3.99019417011D+00

                        real (dp)            :: q,y

                        !
                        !  |X| <= 1.28.
                        !
                        if ( abs ( x ) <= 1.28_dp ) then
                          y = 0.5_dp * x * x
                          q = 0.5_dp - abs(x) * (a1 - a2*y/( y + a3 - a4/( y + a5 &
                            + a6/(y+a7) ) ) )
                        !
                        !  1.28 < |X| <= 12.7
                        !
                        else if ( abs ( x ) <= 12.7_dp ) then
                          y = 0.5_dp * x * x
                          q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
                            + b2 / ( abs ( x ) + b3 &
                            + b4 / ( abs ( x ) - b5 &
                            + b6 / ( abs ( x ) + b7 &
                            - b8 / ( abs ( x ) + b9 &
                            + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
                        !
                        !  12.7 < |X|
                        !
                        else
                          q = 0.0_dp
                        end if
                        !
                        !  Take account of negative X.
                        !
                        if ( x < 0.0_dp ) then
                          cdf = q
                        else
                          cdf = 1.0_dp - q
                        end if
                endfunction

                !! normal_01_cdf_inv inverts the standard normal CDF.
                !  Original Author: John Burkardt (2014)
                !  Modified for syntax, minor passing conventions, purity.
                !  @p gives the probability to invert.
                !  Reference:
                !
                !    Michael Wichura,
                !    Algorithm AS241:
                !    The Percentage Points of the Normal Distribution,
                !    Applied Statistics,
                !    Volume 37, Number 3, pages 477-484, 1988.
                !
                pure function normal_01_cdf_inv ( p ) result(x)
                
                        use env_kindtypes,      only: dp

                        implicit none

                        !!!! Dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp), intent(in   ) :: p
                        !!!! End dummy Arguments !!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp) x
                        
                        !constants used in estimation
                        real    (dp), parameter :: a(8) = (/ &
                                3.3871328727963666080D+00, &
                                1.3314166789178437745D+02, &
                                1.9715909503065514427D+03, &
                                1.3731693765509461125D+04, &
                                4.5921953931549871457D+04, &
                                6.7265770927008700853D+04, &
                                3.3430575583588128105D+04, &
                                2.5090809287301226727D+03 /)
                        real    (dp), parameter :: b(8) = (/ &
                                1.0D+00, &
                                4.2313330701600911252D+01, &
                                6.8718700749205790830D+02, &
                                5.3941960214247511077D+03, &
                                2.1213794301586595867D+04, &
                                3.9307895800092710610D+04, &
                                2.8729085735721942674D+04, &
                                5.2264952788528545610D+03 /)
                        real    (dp), parameter :: c(8) = (/ &
                                1.42343711074968357734D+00, &
                                4.63033784615654529590D+00, &
                                5.76949722146069140550D+00, &
                                3.64784832476320460504D+00, &
                                1.27045825245236838258D+00, &
                                2.41780725177450611770D-01, &
                                2.27238449892691845833D-02, &
                                7.74545014278341407640D-04 /)
                        real    (dp), parameter :: const1 = 0.180625_dp
                        real    (dp), parameter :: const2 = 1.6_dp
                        real    (dp), parameter :: d(8) = (/ &
                                1.0D+00, &
                                2.05319162663775882187D+00, &
                                1.67638483018380384940D+00, &
                                6.89767334985100004550D-01, &
                                1.48103976427480074590D-01, &
                                1.51986665636164571966D-02, &
                                5.47593808499534494600D-04, &
                                1.05075007164441684324D-09 /)
                        real    (dp), parameter :: e(8) = (/ &
                                6.65790464350110377720D+00, &
                                5.46378491116411436990D+00, &
                                1.78482653991729133580D+00, &
                                2.96560571828504891230D-01, &
                                2.65321895265761230930D-02, &
                                1.24266094738807843860D-03, &
                                2.71155556874348757815D-05, &
                                2.01033439929228813265D-07 /)
                        real    (dp), parameter :: f(8) = (/ &
                                1.0D+00, &
                                5.99832206555887937690D-01, &
                                1.36929880922735805310D-01, &
                                1.48753612908506148525D-02, &
                                7.86869131145613259100D-04, &
                                1.84631831751005468180D-05, &
                                1.42151175831644588870D-07, &
                                2.04426310338993978564D-15 /)
                        real    (dp), parameter :: split1 = 0.425_dp
                        real    (dp), parameter :: split2 = 5.0_dp

                        !work variables
                        real    (dp)            :: q, r
                        
                        if ( p <= 0.0_dp ) then
                                x = - huge ( x )
                                return
                        end if
                        
                        if ( 1.0_dp <= p ) then
                                x = huge ( x )
                                return
                        end if
                        
                        q = p - 0.5_dp
                        
                        if ( abs (q) <= split1 ) then
                                r = const1 - q * q
                                x = q * r8poly_value_horner ( a, r ) / r8poly_value_horner ( b, r )
                        else
                                if ( q < 0.0_dp ) then
                                        r = p
                                else
                                        r = 1.0_dp - p
                                end if
                              
                                if ( r <= 0.0_dp ) then
                                        x = huge ( x )
                                else
                                        r = sqrt ( - log(r) )
                              
                                        if ( r <= split2 ) then
                                                r = r - const2
                                                x = r8poly_value_horner ( c, r ) / r8poly_value_horner ( d, r )
                                        else
                                                r = r - split2
                                                x = r8poly_value_horner ( e, r ) / r8poly_value_horner ( f, r )
                                        end if
                                end if
                              
                                if ( q < 0.0_dp ) then
                                        x = -x
                                end if
                        end if
                endfunction

                function bca_transform_percentile(percentile,accel,bias) result(perc)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        real    (dp), intent(in   ) :: percentile
                        real    (dp), intent(in   ) :: accel, bias

                        real    (dp)                :: perc

                        real    (dp)                :: arg, z_alpha

                        z_alpha = normal_01_cdf_inv(percentile/100.0_dp)

                        !print*, "percentile", percentile, "z_alpha", z_alpha, "bias", bias, "accel", accel

                        arg  = bias + ( (bias + z_alpha) / (1 - accel*(bias + z_alpha)))

                        perc = int(100.0_dp * normal_01_cdf(arg),si)

                        !print*, "arg", arg, "perc", perc

                endfunction

                !fully generic comparison of two mappings with no constaints on numbers of 
                !assignments. as
                pure function assignment_difference(assignment1,assignment2) result(tally)

                        use env_kindtypes,      only: si

                        implicit none

                        integer (si), intent(in   ) :: assignment1(:), assignment2(:)

                        !return value
                        integer (si)                :: tally

                        integer (si)                :: iter1, iter2
                        logical                     :: status1, status2

                        tally = 0

                        do iter1= 1,        size(assignment1)-1
                        do iter2= iter1 + 1,size(assignment1)
                                status1 = (assignment1(iter1) == assignment1(iter2))
                                status2 = (assignment2(iter1) == assignment2(iter2))
                                if ((       status1  .and. (.not. status2)) .or. &
                                    ((.not. status1) .and.        status2 )) then
                                        tally = tally + 1
                                endif
                        enddo
                        enddo

                endfunction assignment_difference

                !returns a mask excluding outsiders. O(n log n).
                !
                ! @sample : data to derive mask from
                ! @percentile: the percentile, outside of which is filtered:
                !       e.g. 70 with upper=.true. -> top 30 discarded.
                !       e.g. 70 with upper=.true. and lower=.true -> top 30 and bottom
                !            30 discarded (so ~40% remaining).
                ! @upper : whether to filter out the high values.
                ! @lower : whether to filter out the log values.
                function outlier_mask(sample,percentile,upper,lower) result(mask)

                        use env_kindtypes,      only: si, dp
                        use core_sort,          only: indexQsort

                        implicit none

                        real    (dp),          intent(in   ) :: sample(:)
                        integer (si),          intent(in   ) :: percentile
                        logical,     optional, intent(in   ) :: upper, lower

                        !return value
                        logical,     allocatable             :: mask(:)

                        logical                              :: upper_, lower_
                        real    (dp)                         :: n_max_val, n_min_val
                        real    (dp),allocatable             :: sample_copy(:)
                        integer (si)                         :: n_val
                        integer (si),allocatable             :: sort_indices(:)

                        if (present(upper)) then
                                upper_ = upper
                        else
                                upper_ = .true.
                        endif

                        if (present(lower)) then
                                lower_ = lower
                        else
                                lower_ = .true.
                        endif

                        allocate(mask(size(sample)))

                        !exit if no filtering is requested.
                        if (.not. (upper_ .or. lower_)) then
                                mask = .true.
                                return
                        endif

                        !convert percentile value into the n-th top or bottom value to use as an
                        !inclusive threshold. Conservative on including more of the original data.
                        n_val = floor(real(percentile,dp)/100.0_dp*size(sample))

                        !exit if we have nothing to trim.
                        if (n_val == 0) then
                                mask = .true.
                                return
                        endif

                        sample_copy = sample

                        !derive the thresholds to filter with using partial sorting.
                        if (upper_) then
                                !indexQsort has intent in/out on its first argument.
                                call indexQsort(sample_copy,sort_indices,n_val)
                                n_max_val = sample(sort_indices(n_val))
                        else
                                n_max_val = 0 !avoid spurious possible uninitialized variable warnings
                        endif

                        if (lower_) then
                                !indexQsort cannot have a constant/temporary variable as its array.
                                sample_copy = -1.0_dp* sample_copy 
                                call indexQsort( sample_copy,sort_indices,n_val)
                                n_min_val = sample(sort_indices(n_val))
                        else
                                n_min_val = 0 !avoid spurious possible uninitialized variable warnings
                        endif

                        !Fill mask values.
                        if (upper_ .and. lower_) then
                                where (sample > n_max_val .or. sample < n_min_val)
                                        mask = .false.
                                elsewhere
                                        mask = .true.
                                endwhere
                        elseif (upper_) then
                                where (sample > n_max_val)
                                        mask = .false.
                                elsewhere
                                        mask = .true.
                                endwhere
                        elseif (lower_) then
                                where (sample < n_min_val)
                                        mask = .false.
                                elsewhere
                                        mask = .true.
                                endwhere
                        endif


                endfunction outlier_mask

                pure function rog(points,do_center)

                        use env_kindtypes,      only: si,dp

                        implicit none

                        !given (ncases,ndim)
                        real    (dp),          intent(in   ) :: points(:,:)
                        logical,      optional,intent(in   ) :: do_center

                        real    (dp),allocatable             :: rog

                        integer (si)                         :: iter1

                        logical                              :: do_center_

                        real    (dp),allocatable             :: centered_points(:,:)
                        real    (dp),allocatable             :: center_point(:)

                        if (present(do_center)) then
                                do_center_ = do_center
                        else
                                do_center_ = .true.
                        endif

                        if (do_center_) then

                                allocate(centered_points,mold=points)

                                center_point = sum(points,1)/size(points,1)

                                do iter1=1,size(points,1)
                                        centered_points(iter1,:) = points(iter1,:) - center_point
                                enddo

                                rog = sqrt(sum(centered_points**2)/size(points,1))
                        else
                                rog = sqrt(sum(points**2)/size(points,1))
                        endif

                endfunction rog

endmodule core_stat
