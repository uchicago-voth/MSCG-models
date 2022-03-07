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
! This module provides geneal routines for dealing with standard matrix problems. Eventually,
! conditional compilation here should provide compatibility for e.g. MKL compilation, removing complex
! conditionals from the code.

module core_matrix

        implicit none

        private

        public lsolve, spectrum, quadForm, svd, &
                quadratic_min_symmetric_equality_LDL_dp, diagonal, diagonal_worker


        !interface for solving exact linear systems
        interface lsolve
                module procedure lsolve_v,lsolve_m
        end interface

        !solves eigenvalue problems.
        interface spectrum
                module procedure spectrum_symmetric_wrapper, gen_spectrum_symmetric_wrapper
        end interface

        !solves quadratic expressions, e.g. t(A) %*% B %*% A
        interface quadForm
                module procedure quadForm_mat, quadForm_vector
        end interface

        !solves singular value decomposition problems
        interface svd
                module procedure singular_value_decomp_dp
        end interface

        contains
                !linear system A x = b with x as a 1D array (wraps matrix version, uses a copy).
                ! @A         : design matrix
                ! @x         : matrix to solve
                ! @b         : system parameters
                ! @symmetric : whether A is symmetric
                pure function lsolve_v(A,b,symmetric) result(x)

                        use env_kindtypes,      only: dp_x

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )          :: A(:,:)
                        real    (dp_x),intent(in   )          :: b(:)
                        logical,       intent(in   ),optional :: symmetric
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !lapack call variables
                        real    (dp_x),allocatable            :: x_mat(:,:), b_mat(:,:), x(:)
                        logical                               :: symmetric_

                        if (present(symmetric)) then
                                symmetric_ = symmetric
                        else
                                symmetric_ = .false.
                        endif

                        allocate(b_mat(size(b),1))

                        b_mat(:,1) = b

                        x_mat = lsolve_m(A,b_mat,symmetric_)

                        x = x_mat(:,1)

                endfunction lsolve_v

                !linear system A x = b with x as a 2D array.
                ! @A         : design matrix
                ! @x         : matrix to solve
                ! @b         : system parameters
                ! @symmetric : whether A is symmetric
                pure function lsolve_m(A,B,symmetric) result(X)

                        use env_kindtypes,      only: dp_x

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )          :: A(:,:)
                        real    (dp_x),intent(in   )          :: B(:,:)
                        logical,       intent(in   ),optional :: symmetric
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !lapack call variables
                        real    (dp_x),allocatable            :: X(:,:)
                        logical                               :: symmetric_

                        if (present(symmetric)) then
                                symmetric_ = symmetric
                        else
                                symmetric_ = .false.
                        endif

                        if (symmetric_) then
                                X = solve_symmetric_ldl_dp(A,B)
                        else
                                X = solve_general_lu_p_dp(A,B)
                        endif

                endfunction lsolve_m

                !provides a direct solution to an equality constrained

                !quadratic minimization problem
                !of the following form:
                !
                ! @:
                ! arg min_q [ (1/2)(q*)Bq - bq ] s.t. Aq = constr
                !
                ! This is done via langrangian multipliers and LDL factorization.
                ! This routine is for real(dp_x).
                !
                ! notes for this algorithm may be found at:
                ! http://www.math.uh.edu/~rohop/fall_06/Chapter3.pdf
                !
                ! Note that we have constraints on B

                pure function quadratic_min_symmetric_equality_LDL_dp(A,q_B,l_b,constr) result(minVector)

                        use env_kindtypes,      only: dp_x

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )          :: A(:,:), q_B(:,:), l_b(:)
                        real    (dp_x),intent(in   ),optional :: constr(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp_x),allocatable      :: minVector(:)

                        !optional arg counterparts
                        real    (dp_x),allocatable      :: constr_(:)

                        !local variables
                        real    (dp_x),allocatable      :: KKT(:,:)
                        real    (dp_x),allocatable      :: toSolve(:,:), compositeTarget(:,:)

                        if (present(constr)) then
                                constr_ = constr
                        else
                                allocate(constr_(size(A,1)))
                                constr_ = [ 0 ]
                        endif

                        !allocate memory for the extended problem matrices.
                        allocate(KKT((size(q_B,1)+size(A,1)),&
                                      size(q_B,1)+size(A,1)))

                        allocate(compositeTarget(size(l_b,1)+size(constr),1))

                        !setup extended KKT matrix
                        KKT = 0

                        KKT(1:size(q_B,1),1:size(q_B,2)) = q_B

                        KKT((size(q_B,1)+1):(size(q_B,1)+1+size(A,1)),&
                             1:size(A,2))                = A

                        KKT(1:size(A,2),&
                             (size(q_B,2)+1):(size(q_B,2)+1+size(A,1))) &
                                                         = transpose(A)

                        compositeTarget(1:size(l_b),1)   = l_b
                        compositeTarget((size(l_b)+1):(size(l_b) + 1 + size(constr_)),1) = constr_

                        toSolve = solve_symmetric_ldl_dp(KKT,compositeTarget)

                        minVector = toSolve(1:size(l_b),1)

                endfunction quadratic_min_symmetric_equality_ldl_dp

                !Interfaces with lapack to solve a symmetric determined linear

                !system via LDL.
                ! @:
                ! Ax = B
                pure function solve_symmetric_ldl_dp(A,B) result(B_)

                        use env_kindtypes,      only: si_x, dp_x

                        implicit none

                        interface
                                pure subroutine dsytrf( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
                                        import si_x, dp_x
                                        character,           intent(in   ) :: UPLO
                                        integer (KIND=si_x), intent(in   ) :: N,LDA,LWORK
                                        integer (KIND=si_x), intent(  out) :: IPIV(*)
                                        real    (KIND=dp_x), intent(inout) :: A(LDA,*)
                                        real    (KIND=dp_x), intent(  out) :: WORK(*)
                                        integer (KIND=si_x), intent(  out) :: INFO
                                endsubroutine dsytrf
                        end interface

                        interface
                                pure subroutine dsytrs( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
                                        import si_x, dp_x
                                        character,           intent(in   ) :: UPLO
                                        integer (KIND=si_x), intent(in   ) :: N, NRHS, LDA, LDB
                                        integer (KIND=si_x), intent(  out) :: IPIV(*)
                                        real    (KIND=dp_x), intent(inout) :: A(LDA,*), B(LDB,*)
                                        integer (KIND=si_x), intent(  out) :: INFO
                                endsubroutine dsytrs
                        end interface

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )          :: A(:,:), B(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value

                        !local vars
                        real    (dp_x),allocatable      :: B_(:,:)
                        real    (dp_x),allocatable      :: A_(:,:)

                        real    (dp_x)                  :: scratch(1)

                        !lapack call variables
                        character                       :: uplo
                        integer (si_x)                  :: info, lwork
                        integer (si_x),allocatable      :: ipiv(:)
                        real    (dp_x),allocatable      :: work(:)

                        !set upper option
                        uplo ='U'

                        !lapack calls are destructive, so we make copies.

                        !copy A, B
                        B_ = B
                        A_ = A

                        !first we factor A -> LDL* via dsytrs

                        !allocate relevant dsytrs variables

                        allocate(ipiv(size(A_))) !to contain info about pivots

                        !Now ask what the optimal work size is.
                        !DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
                        call dsytrf(uplo,    &
                                 size(A_,1), &
                                 A_,         &
                                 size(A_,1), &
                                 ipiv,       &
                                 scratch,    & !scratch(1) gets lwork.
                                 -1,         & !This indicates that this is a work size query.
                                 info)

                        lwork = nint(scratch(1),si_x)

                        !allocate work array accordingly.
                        allocate(work(lwork))

                        !!!! do factoring call. !!!!!!!!!!!!

                        !DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
                        call dsytrf(uplo,    &
                                 size(A_,1), &
                                 A_,         &
                                 size(A_,1), &
                                 ipiv,       &
                                 work,       &
                                 lwork,      &
                                 info)

                        !!! end factoring call !!!!!!!!!!!!!

                        if (info /= 0) then
                                B_ = 0 !if not successful bail
                                return
                        endif

                        !allocate relevant dsytrs variables

                        !DSYTRS (UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
                        call dsytrs(uplo,     &
                                size(A_,1),   &
                                size(B,2),    &
                                A_,           &
                                size(A_,1),   &
                                ipiv,         &
                                B_,        & !gets overridden with answer
                                size(B_,1),&
                                info)

                        if (info /= 0) then
                                B_ = 0 !if not successful bail
                                return
                        endif

                        !all good

                        return

                endfunction

                !Solves a general EXACT pivot system by LU with pivoting (lu,p)
                ! @:
                !A*X=B.
                pure function solve_general_lu_p_dp(A,B) result(B_)

                        use env_kindtypes,      only: si_x, dp_x

                        implicit none

                        interface
                                pure subroutine dgesv(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
                                        import si_x, dp_x
                                        integer (KIND=si_x), intent(in   ) :: N,NRHS,LDA,LDB
                                        integer (KIND=si_x), intent(  out) :: IPIV(*)
                                        real    (KIND=dp_x), intent(inout) :: A(LDA,N)
                                        real    (KIND=dp_x), intent(inout) :: B(LDB,NRHS)
                                        integer (KIND=si_x), intent(  out) :: INFO
                                endsubroutine dgesv
                        end interface

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )          :: A(:,:)
                        real    (dp_x),intent(in   )          :: B(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp_x),allocatable      :: A_(:,:),B_(:,:)

                        !lapack call variables
                        integer (si_x)                  :: info, lda, ldb, n, nrhs
                        integer (si_x),allocatable      :: ipiv(:)

                        !lapack calls are destructive, so we make copies.

                        !copy A,B
                        allocate(A_,source=A)
                        allocate(B_,source=B)

                        !prepare parameters for dgesv call.
                        n    = size(A,1)
                        lda  = size(A,2)
                        ldb  = size(B,1)
                        nrhs = size(B,2)
                        allocate(ipiv(n))

                        call dgesv(n,   &
                                   nrhs,&
                                   A_,  &
                                   lda, &
                                   ipiv,&
                                   B_,  &
                                   ldb, &
                                   info)

                        return

                endfunction

                ! spectrum provides the eigenvalues and optionally eigenvectors of matrix M.
                ! @M : Matrix to be decomposed.
                ! @D : Returned array of eigenvalues
                ! @Q : Returned 2 array of eigenvectors (optional)
                ! @num_top_eigenvalues : Return partial results corresponding to @ top eigenvalues (optional)
                ! @num_bottom_eigenvalues : Return partial results corresponding to @ bottom eigenvalues (optional)
                !
                ! NOTE: num_top_eigenvalues takes precidence over num_bottom_eigenvalues if both are present.
                !
                ! This wrapper allocates variables if neccesary, and can be expanded to check arguments.
                pure subroutine spectrum_symmetric_wrapper(M,d,Q,num_top_eigenvalues,num_bottom_eigenvalues)

                        use env_kindtypes,      only: dp_x, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )             :: M(:,:)
                        real    (dp_x),intent(  out),allocatable :: d(:)
                        real    (dp_x),intent(  out),allocatable,optional,target&
                                                                 :: Q(:,:)
                        integer (si),  intent(in   ),optional    :: num_top_eigenvalues
                        integer (si),  intent(in   ),optional    :: num_bottom_eigenvalues
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp_x),target                    :: dummy_2array(0,0)
                        integer (si)                             :: num_eigenvalues
                        integer (si),allocatable                 :: eig_range(:)
                        real    (dp_x),pointer                   :: p_eigenvectors(:,:)

                        !We analyze arguments and allocate as neccesary.

                        if (present(num_top_eigenvalues)) then
                                num_eigenvalues = num_top_eigenvalues
                        else if (present(num_bottom_eigenvalues)) then
                                num_eigenvalues = num_bottom_eigenvalues
                        else
                                num_eigenvalues = size(M,1)
                        endif

                        !allocate output variables if needed and present
                        if (.not.(allocated(D))) then
                                allocate(D(num_eigenvalues))
                        else
                                if (size(D,1) /= num_eigenvalues) then
                                        deallocate(D)
                                        allocate(D(num_eigenvalues))
                                endif
                        endif

                        if (present(Q)) then
                                if (.not.(allocated(Q))) then
                                        allocate(Q(size(M,1),num_eigenvalues))
                                        p_eigenvectors => Q
                                else
                                        if ((size(Q,1) /= size(M,1)) .or.&
                                            (size(Q,2) /= num_eigenvalues))      then
                                                deallocate(Q)
                                                allocate(Q(size(M,2),num_eigenvalues))
                                        endif
                                        p_eigenvectors => Q
                                endif
                        else
                                p_eigenvectors => dummy_2array
                        endif

                        !if partial, create the rang eof eigenvalues to be passed.
                        if (present(num_top_eigenvalues)) then
                                allocate(eig_range(2))
                                eig_range(1) = size(M,1) - num_top_eigenvalues + 1
                                eig_range(2) = size(M,1)
                                !print *, 'eig range', eig_range
                        else if (present(num_bottom_eigenvalues)) then
                                allocate(eig_range(2))
                                eig_range(1) = 1
                                eig_range(2) = num_bottom_eigenvalues
                                !print *, 'eig range', eig_range
                        else
                                allocate(eig_range(0))
                        endif

                        !call work routine.
                        call spectrum_symmetric(M,D,p_eigenvectors,eig_range)

                endsubroutine spectrum_symmetric_wrapper

                ! spectrum provides the eigenvalues and optionally eigenvectors of matrix M.
                ! @M : Matrix to be decomposed.
                ! @eigenvalues : Returned array of eigenvalues (ascending)
                ! @eigenvectors : Returned 2 array of eigenvectors (if size == 0, ignored) (optional)
                ! @eig_range : array of upper and lower bounds, which delimit the indices of returned eigenvalue/eigenvectors.
                !                    (if size /= 2 ignored) (optional)
                !
                ! This assumeds D and Q to be appropriately allocated.
                pure subroutine spectrum_symmetric(M,eigenvalues,eigenvectors,eig_range)


                        use env_kindtypes,      only: si_x, dp_x

                        implicit none

                        interface
                                pure subroutine dsyevr(jobz, range, uplo, n, A, lda, vl, vu,  il, iu, &
                                                       abstol, m, w, Z, ldz, isuppz, work, lwork, iwork, liwork, info)
                                        import si_x, dp_x
                                        character,           intent(in   ) :: jobz, range, uplo
                                        integer (KIND=si_x), intent(in   ) :: n,lda,lwork
                                        real    (KIND=dp_x), intent(inout) :: A(lda,*)
                                        real    (KIND=dp_x), intent(in   ) :: vl, vu
                                        integer (KIND=si_x), intent(in   ) :: il, iu
                                        real    (KIND=dp_x), intent(in   ) :: abstol
                                        integer (KIND=si_x), intent(in   ) :: ldz
                                        integer (KIND=si_x), intent(  out) :: isuppz(*)
                                        integer (KIND=si_x), intent(  out) :: m
                                        real    (KIND=dp_x), intent(  out) :: Z(ldz,*)
                                        integer (KIND=si_x), intent(  out) :: info
                                        real    (KIND=dp_x), intent(  out) :: work(*), w(*)
                                        integer (KIND=si_x), intent(  out) :: iwork(*)
                                        integer (KIND=si_x), intent(in   ) :: liwork
                                endsubroutine dsyevr
                        end interface

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )                 :: M(:,:)
                        real    (dp_x),intent(  out),         target :: eigenvalues(:)
                        real    (dp_x),intent(  out),optional,target :: eigenvectors(:,:)
                        integer (si_x),intent(in   ),optional        :: eig_range(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!! local variables
                        real      (dp_x),allocatable    :: M_(:,:)

                        !!! values for diagonalization routine. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character (len=1)               :: jobz,uplo,range !! paramters for lapack call.
                        integer   (si_x)                :: lwork,liwork    !! size of double work array for

                                                                           !   lapack work to be computed below
                                                                           !   (depends on natoms)
                        real      (dp_x), allocatable   :: work(:)         !! work array for dsyevd subroutine
                        integer   (si_x), allocatable   :: iwork(:)        !! work array for dsyevd subroutine
                        integer   (si_x)                :: info            !! output integer (si_x) from dsyevd

                                                                           !   which describes if
                                                                           !   routine completed successfully
                        integer   (si_x)                :: n,lda,ldz       !! output integer (si_x) from dsyevd

                        integer   (si_x), allocatable   :: isuppz(:)
                        real      (dp_x)                :: vl, vu
                        integer   (si_x)                :: il, iu, num_eigenvalues
                        real      (dp_x),pointer        :: p_eigenvectors(:,:)
                        real      (dp_x),pointer        :: p_eigenvalues(:)

                        ! lapack is destructive, make copies.
                        allocate(M_,source=M)

                        ! set lapack options

                        !detect if we're doing a partial dcomposition. 'A' is full, 'I' is partial.
                        !set routine indices (il iu) if needed.
                        !if doing a partial decomp, we are passed an eigenvalue array of size < rank.
                        !however, lapack always needs a full range return array (according to standard).
                        !use pointers to make such tricks readable.
                        range= 'A'
                        if (present(eig_range)) then
                                if (size(eig_range) == 2) then
                                        range = 'I'
                                        il = eig_range(1)
                                        iu = eig_range(2)
                                        allocate(p_eigenvalues(size(M,1)))
                                else
                                        p_eigenvalues => eigenvalues
                                endif
                        endif

                        ! If we are given space to put the eigenvectors, mark them to be computed.
                        if (present(eigenvectors)) then
                                if (size(eigenvectors) /= 0) then
                                        p_eigenvectors => eigenvectors
                                        jobz='V'

                                else
                                        allocate(p_eigenvectors(1,1))
                                        jobz='N'
                                endif
                        else
                                allocate(p_eigenvectors(1,1))
                                jobz='N'
                        endif

                        uplo = 'U'  !! Parameters for lapack routine. Use the upper
                                    !  half of the matrix (U)


                        n   = size(M,1)
                        lda = size(M,1)

                        ldz = size(p_eigenvectors,1)
                        !print*, "LDZ", ldz, shape(p_eigenvectors)

                        allocate(isuppz(2*size(M,1)))

                        ! workspace query
                        allocate(work(1))
                        allocate(iwork(1))

                        call dsyevr(  jobz = jobz, &
                                     range = range, &
                                      uplo = uplo, &
                                         n = n, &
                                         A = M_, &
                                       lda = lda, &
                                        vl = vl,&
                                        vu = vu,&
                                        il = il,&
                                        iu = iu,&
                                    abstol = 0.0_dp_x,&
                                         m = num_eigenvalues,&
                                         W = p_eigenvalues,&
                                         Z = p_eigenvectors,&
                                       ldz = ldz,&
                                    isuppz = isuppz,&
                                      work = work, &
                                     lwork = -1, &
                                     iwork = iwork, &
                                    liwork = -1, &
                                      info = info)

                        lwork  = nint(work(1))
                        liwork =      iwork(1)

                        !!allocate arrays for lapack subroutine
                        deallocate(work)
                        allocate(work(lwork))
                        deallocate(iwork)
                        allocate(iwork(liwork))

                        call dsyevr(  jobz = jobz, &
                                     range = range, &
                                      uplo = uplo, &
                                         n = n, &
                                         A = M_, &
                                       lda = lda, &
                                        vl = vl,&         !ignored currently
                                        vu = vu,&         !ignored currently
                                        il = il,&         !ignored currently
                                        iu = iu,&         !ignored currently
                                    abstol = 0.0_dp_x,&
                                         m = num_eigenvalues,&
                                         W = p_eigenvalues,&
                                         Z = p_eigenvectors,&
                                       ldz = ldz,&
                                    isuppz = isuppz,&
                                      work = work, &
                                     lwork = lwork, &
                                     iwork = iwork, &
                                    liwork = liwork, &
                                      info = info)

                        if (range == 'I') then
                                eigenvalues = p_eigenvalues(1:num_eigenvalues)
                        endif

                endsubroutine spectrum_symmetric

                ! spectrum provides the eigenvalues and optionally eigenvectors of matrix M.
                !
                ! Form: M x = \lambda B x
                !
                ! @M : Matrix to be decomposed.
                ! @B : Modifier matrix
                ! @d : Returned array of eigenvalues
                ! @Q : Returned 2 array of eigenvectors (optional)
                ! @num_top_eigenvalues : Return partial results corresponding to @ top eigenvalues (optional)
                ! @num_bottom_eigenvalues : Return partial results corresponding to @ bottom eigenvalues (optional)
                !
                ! NOTE: num_top_eigenvalues takes precidence over num_bottom_eigenvalues if both are present.
                !
                ! This wrapper allocates variables if neccesary, and can be expanded to check arguments.
                pure subroutine gen_spectrum_symmetric_wrapper(M,B,d,Q,num_top_eigenvalues,num_bottom_eigenvalues)

                        use env_kindtypes,      only: dp_x, si

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )             :: M(:,:)
                        real    (dp_x),intent(in   )             :: B(:,:)
                        real    (dp_x),intent(  out),allocatable :: d(:)
                        real    (dp_x),intent(  out),allocatable,optional,target&
                                                                 :: Q(:,:)
                        integer (si),  intent(in   ),optional    :: num_top_eigenvalues
                        integer (si),  intent(in   ),optional    :: num_bottom_eigenvalues
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp_x),target                    :: dummy_2array(0,0)
                        integer (si)                             :: num_eigenvalues
                        integer (si),allocatable                 :: eig_range(:)
                        real    (dp_x),pointer                   :: p_eigenvectors(:,:)

                        !We analyze arguments and allocate as neccesary.

                        if (present(num_top_eigenvalues)) then
                                num_eigenvalues = num_top_eigenvalues
                        else if (present(num_bottom_eigenvalues)) then
                                num_eigenvalues = num_bottom_eigenvalues
                        else
                                num_eigenvalues = size(M,1)
                        endif

                        !allocate output variables if needed and present
                        if (.not.(allocated(D))) then
                                allocate(D(num_eigenvalues))
                        else
                                if (size(D,1) /= num_eigenvalues) then
                                        deallocate(D)
                                        allocate(D(num_eigenvalues))
                                endif
                        endif

                        if (present(Q)) then
                                if (.not.(allocated(Q))) then
                                        allocate(Q(size(M,1),num_eigenvalues))
                                        p_eigenvectors => Q
                                else
                                        if ((size(Q,1) /= size(M,1)) .or.&
                                            (size(Q,2) /= num_eigenvalues))      then
                                                deallocate(Q)
                                                allocate(Q(size(M,2),num_eigenvalues))
                                        endif
                                        p_eigenvectors => Q
                                endif
                        else
                                p_eigenvectors => dummy_2array
                        endif

                        !if partial, create the rang eof eigenvalues to be passed.
                        if (present(num_top_eigenvalues)) then
                                allocate(eig_range(2))
                                eig_range(1) = size(M,1) - num_top_eigenvalues + 1
                                eig_range(2) = size(M,1)
                                !print *, 'eig range', eig_range
                        else if (present(num_bottom_eigenvalues)) then
                                allocate(eig_range(2))
                                eig_range(1) = 1
                                eig_range(2) = num_bottom_eigenvalues
                                !print *, 'eig range', eig_range
                        else
                                allocate(eig_range(0))
                        endif

                        !call work routine.
                        call gen_spectrum_symmetric(M,B,D,p_eigenvectors,eig_range)

                endsubroutine gen_spectrum_symmetric_wrapper

                ! spectrum provides the eigenvalues and optionally eigenvectors of matrix M.
                !
                ! Form: M x = \lambda B x
                !
                ! @M : Matrix to be decomposed.
                ! @B : Modifier matrix
                ! @eigenvalues : Returned array of eigenvalues (ascending)
                ! @eigenvectors : Returned 2 array of eigenvectors (if size == 0, ignored) (optional)
                ! @eig_range : array of upper and lower bounds, which delimit the indices of

                !              returned eigenvalue/eigenvectors.
                !                    (if size /= 2 ignored) (optional)
                !
                ! This assumeds D and Q to be appropriately allocated.
                pure subroutine gen_spectrum_symmetric(M,B,eigenvalues,eigenvectors,eig_range)


                        use env_kindtypes,      only: si_x, dp_x

                        implicit none

                        interface
                                pure subroutine dsygvx(itype, jobz, range, uplo, n, A, lda, B, ldb, vl, vu,&
                                                       il, iu, abstol, m, w, Z, ldz, work, lwork, iwork,&
                                                       ifail, info)
                                        import si_x, dp_x
                                        integer (si_x), intent(in   ) :: itype
                                        character,      intent(in   ) :: jobz, range, uplo
                                        integer (si_x), intent(in   ) :: n,lda,ldb,lwork
                                        real    (dp_x), intent(inout) :: A(lda,*)
                                        real    (dp_x), intent(inout) :: B(ldb,*)
                                        real    (dp_x), intent(in   ) :: vl, vu
                                        integer (si_x), intent(in   ) :: il, iu
                                        real    (dp_x), intent(in   ) :: abstol
                                        integer (si_x), intent(in   ) :: ldz
                                        integer (si_x), intent(  out) :: m
                                        real    (dp_x), intent(  out) :: Z(ldz,*)
                                        real    (dp_x), intent(  out) :: work(*), w(*)
                                        integer (si_x), intent(  out) :: iwork(*)
                                        integer (si_x), intent(  out) :: ifail(*)
                                        integer (si_x), intent(  out) :: info
                                endsubroutine dsygvx
                        end interface

                        interface
                                pure function dlamch(cmach)
                                        import dp_x
                                        character,      intent(in   ) :: cmach
                                        real    (dp_x)                :: dlamch
                                endfunction dlamch
                        end interface

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp_x),intent(in   )                 :: M(:,:), B(:,:)
                        real    (dp_x),intent(  out),         target :: eigenvalues(:)
                        real    (dp_x),intent(  out),optional,target :: eigenvectors(:,:)
                        integer (si_x),intent(in   ),optional        :: eig_range(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !!! local variables
                        real      (dp_x),allocatable   :: M_(:,:), B_(:,:)

                        !!! values for diagonalization routine. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        character (len=1)               :: jobz,uplo,range !! paramters for lapack call.
                        integer   (si_x)                :: lwork           !! size of double work array for

                                                                           !   lapack work to be computed below
                                                                           !   (depends on natoms)
                        real      (dp_x),allocatable    :: work(:)         !! work array for dsyevd subroutine
                        integer   (si_x),allocatable    :: iwork(:)        !! work array for dsyevd subroutine
                        integer   (si_x)                :: info            !! output integer (si_x) from dsyevd

                                                                           !   which describes if
                                                                           !   routine completed successfully
                        integer   (si_x),allocatable    :: ifail(:)        !! Marks eigenvectors with failed convergence

                        integer   (si_x)                :: n,lda,ldb,ldz   !! output integer (si_x) from dsyevd

                        real      (dp_x)                :: vl, vu
                        integer   (si_x)                :: il, iu, num_eigenvalues
                        real      (dp_x),pointer        :: p_eigenvectors(:,:)
                        real      (dp_x),pointer        :: p_eigenvalues(:)

                        integer   (si_x),parameter      :: itype = 1       !! Specifies the type of eigenvalue.
                                                                           !  Set to M x = \lambda B x

                        ! lapack is destructive, make copies.
                        allocate(M_,source=M)
                        allocate(B_,source=B)

                        !!!!!!!!!!!! set lapack options !!!!!!!!!!!!

                        !detect if we're doing a partial dcomposition. 'A' is full, 'I' is partial.
                        !set routine indices (il iu) if needed.
                        !if doing a partial decomp, we are passed an eigenvalue array of size < rank.
                        !however, lapack always needs a full range return array (according to standard).
                        !use pointers to make such tricks readable.
                        range= 'A'
                        if (present(eig_range)) then
                                if (size(eig_range) == 2) then
                                        range = 'I'
                                        il = eig_range(1)
                                        iu = eig_range(2)
                                        allocate(p_eigenvalues(size(M,1)))
                                else
                                        p_eigenvalues => eigenvalues
                                endif
                        endif

                        ! If we are given space to put the eigenvectors, mark them to be computed.
                        if (present(eigenvectors)) then
                                if (size(eigenvectors) /= 0) then
                                        p_eigenvectors => eigenvectors
                                        jobz='V'

                                else
                                        allocate(p_eigenvectors(1,1))
                                        jobz='N'
                                endif
                        else
                                allocate(p_eigenvectors(1,1))
                                jobz='N'
                        endif

                        uplo = 'U'  !! Parameters for lapack routine. Use the upper
                                    !  half of the matrix (U)


                        n   = size(M,1)
                        lda = size(M,1)
                        ldb = size(B,1)

                        ldz = size(p_eigenvectors,1)
                        !print*, "LDZ", ldz, shape(p_eigenvectors)

                        allocate(iwork(5*n))
                        allocate(ifail(n))

                        ! workspace query
                        allocate(work(1))

                        call dsygvx( itype = itype,&
                                      jobz = jobz,&
                                     range = range,&
                                      uplo = uplo,&
                                         n = n,&
                                         A = M_,&
                                       lda = lda,&
                                         B = B_,&
                                       ldb = ldb,&
                                        vl = vl,&
                                        vu = vu,&
                                        il = il,&
                                        iu = iu,&
                                    abstol = dlamch('C'),&
                                         m = num_eigenvalues,&
                                         W = p_eigenvalues,&
                                         Z = p_eigenvectors,&
                                       ldz = ldz,&
                                      work = work,&
                                     lwork = -1,&
                                     iwork = iwork,&
                                     ifail = ifail,&
                                      info = info)

                        lwork  = nint(work(1))

                        !!allocate arrays for lapack subroutine
                        deallocate(work)
                        allocate(work(lwork))

                        call dsygvx( itype = itype,&
                                      jobz = jobz, &
                                     range = range, &
                                      uplo = uplo, &
                                         n = n, &
                                         A = M_, &
                                       lda = lda, &
                                         B = B_,&
                                       ldb = ldb,&
                                        vl = vl,&         !ignored currently
                                        vu = vu,&         !ignored currently
                                        il = il,&
                                        iu = iu,&
                                    abstol = dlamch('C'),&
                                         m = num_eigenvalues,&
                                         W = p_eigenvalues,&
                                         Z = p_eigenvectors,&
                                       ldz = ldz,&
                                      work = work, &
                                     lwork = lwork, &
                                     iwork = iwork, &
                                     ifail = ifail, &
                                      info = info)

                        if (range == 'I') then
                                eigenvalues = p_eigenvalues(1:num_eigenvalues)
                        endif

                endsubroutine gen_spectrum_symmetric


                ! values takes precidence over n and length
                pure function diagonal(n,length,values) result(returnMatrix)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   ),optional :: n
                        integer (si),intent(in   ),optional :: length
                        real    (dp),intent(in   ),optional :: values(:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),allocatable            :: returnMatrix(:,:),content(:)

                        if (present(values)) then
                                allocate(returnMatrix(size(values),size(values)))
                                call diagonal_worker(values,returnMatrix)
                                return
                        elseif (present(n) .and. present(length)) then
                                allocate(returnMatrix(length,length))
                                allocate(content(length))
                                content = n
                                call diagonal_worker(content,returnMatrix)
                                return
                        else
                                allocate(returnMatrix(1,1))
                                returnMatrix = 0
                                return
                        endif

                 endfunction diagonal

                 pure subroutine diagonal_worker(content,matrix)

                        use env_kindtypes,      only: si, dp

                        implicit none

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   ) :: content(:)
                        real    (dp),intent(  out) :: matrix(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        integer (si)               :: i,j, limit

                        limit = size(content,1)

                        !Fortran trick for diagonal matrices. Not the most efficient.
                        ForAll(i = 1:limit, j = 1:limit) matrix(i,j) = real((i/j)*(j/i),dp)*content(i)

                 endsubroutine diagonal_worker

                 !computes t(A)*X*A type produces with minimal copying.
                 pure function quadForm_mat(center,side,symmetric,rightTranspose) result(toReturn)

                        use env_kindtypes,      only: dp, si_x, dp_x

                        implicit none

                        !http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3.html
                        interface
                                pure subroutine dgemm(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
                                        import si_x, dp_x
                                        character,     intent(in   ) :: transA, transB
                                        integer (si_x),intent(in   ) :: m, n, k, lda, ldb, ldc
                                        real    (dp_x),intent(in   ) :: alpha, beta, A(lda, *), B(ldb,*)
                                        real    (dp_x),intent(inout) :: C(ldc, *)
                                endsubroutine dgemm
                        end interface

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )          :: center(:,:)
                        real    (dp),intent(in   )          :: side(:,:)
                        logical     ,intent(in   ),optional :: symmetric, rightTranspose
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp_x),allocatable  :: toReturn(:,:)

                        !optional argument counterparts
                        logical                     :: symmetric_
                        logical                     :: rightTranspose_

                        !local variables
                        real    (dp_x),allocatable  :: scratchMatrix(:,:)

                        !lapack call variables
                        character                   :: transA,transB
                        integer (si_x)              :: m, n, k, lda, ldb, ldc
                        real    (dp_x)              :: alpha, beta

                        if (present(symmetric)) then
                                symmetric_ = symmetric
                        else
                                symmetric_ = .false.
                        endif

                        if (present(rightTranspose)) then
                                rightTranspose_ = rightTranspose
                        else
                                rightTranspose_ = .false.
                        endif

                        !currently, we treat the symmetric case identically.
                        if (symmetric_) then
                                return
                        else
                                if (rightTranspose_) then
                                        allocate(scratchMatrix(size(side,1),size(center,2)))
                                        allocate(toReturn(size(side,1),size(center,2)))

                                        beta         =  0  !we don't have C.
                                        alpha        =  1  !we don't want scaling.

                                        !do [side * center] first
                                        !   [  A  *   B   ]
                                        transA = 'n' !don't transpose side
                                        transB = 'n' !don't transpose center
                                        m      = size(side,1)   !num rows of side
                                        n      = size(center,2) !num columns center
                                        k      = size(side,2)   !num columns side
                                        lda    = size(side,1)   !LD side
                                        ldb    = size(center,1) !LD center
                                        ldc    = size(scratchMatrix,1) !LD center
                                        call dgemm(transA,transB,m,n,k,alpha,&
                                                        side,  lda,&
                                                        center,ldb,&
                                                        beta,&
                                                        scratchMatrix,ldc)

                                        !then do [side * Center] * t(size)
                                        !        [     A       ] * t(B)
                                        transA = 'n' !don't transpose toReturn
                                        transB = 't' !transpose side

                                        m      = size(scratchMatrix,1)   !num rows of A
                                        n      = size(side,1)            !num columns t(B)
                                        k      = size(scratchMatrix,2)   !num columns A
                                        lda    = size(scratchMatrix,1)   !LD A

                                        ldb    = size(side,1)            !LD B
                                        ldc    = size(toReturn,1)        !LD C
                                        call dgemm(transA,transB,m,n,k,alpha,&
                                                        scratchMatrix,  lda,&
                                                        side,ldb,&
                                                        beta,&
                                                        toReturn,ldc)
                                        return
                                else
                                        allocate(scratchMatrix(size(side,2),size(center,2)))
                                        allocate(toReturn(size(side,2),size(center,2)))

                                        beta         =  0  !we don't have C.
                                        alpha        =  1  !we don't want scaling.

                                        !do [t(side) * center] first
                                        transA = 't' !transpose side
                                        transB = 'n' !don't transpose center
                                        m      = size(side,2)   !num rows of t(side)
                                        n      = size(center,2) !num columns center
                                        k      = size(side,1)   !num columns t(side)
                                        lda    = size(side,1)   !LD side
                                        ldb    = size(center,1) !LD center
                                        ldc    = size(scratchMatrix,1) !LD center
                                        call dgemm(transA,transB,m,n,k,alpha,&
                                                        side,  lda,&
                                                        center,ldb,&
                                                        beta,&
                                                        scratchMatrix,ldc)

                                        !then do [t(A) * B] * A
                                        transA = 'n' !don't transpose toReturn
                                        !transB = 'n' !don't transpose center (already set)
                                        m      = size(scratchMatrix,1)   !num rows of t(side)
                                        n      = size(side,2) !num columns center
                                        k      = size(scratchMatrix,2)   !num columns t(side)
                                        lda    = size(scratchMatrix,1)   !LD side
                                        ldb    = size(side,1) !LD center
                                        ldc    = size(toReturn,1) !LD center
                                        call dgemm(transA,transB,m,n,k,alpha,&
                                                        scratchMatrix,  lda,&
                                                        side,ldb,&
                                                        beta,&
                                                        toReturn,ldc)
                                        return
                                endif
                        endif

                 endfunction quadForm_mat

                 !computes t(A)*X*A type produces with minimal copying, where A is a dim=1 array.

                 !we use pointer to avoid copies, but to present two dimensional versions of 
                 !one dimensional arguments.
                 pure function quadForm_vector(center,side,symmetric) result(toReturn)

                        use env_kindtypes,      only: dp

                        implicit none
                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )          :: center(:,:)
                        real    (dp),intent(in   ),target   :: side(:)
                        logical     ,intent(in   ),optional :: symmetric
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        !return value
                        real    (dp)                :: toReturn

                        real    (dp)                :: toReturn_mat(1,1)

                        !optional argument counterparts
                        logical                     :: symmetric_

                        logical,parameter           :: right_transpose = .false.

                        !pointers
                        real    (dp),pointer        :: s2_side(:,:)

                        !hopefully the compiler will optimize out this copy.
                        !we can't use pointers because of purity.
                        allocate(s2_side(size(side,1),1))
                        s2_side(1:size(side,1),1) = side(1:size(side,1))

                        if (present(symmetric)) then
                                symmetric_ = symmetric
                        else
                                symmetric_ = .false.
                        endif

                        toReturn_mat = quadForm_mat(center,s2_side,symmetric,right_transpose)

                        toReturn = toReturn_mat(1,1)

                 endfunction quadForm_vector


                 !computes SVD of the following style (real arithmetic):
                 ! M = U E V*

                 ! Currently, no partial calculation is accepted.
                 pure subroutine singular_value_decomp_dp(inputMatrix,singular,left,right_trans)

                        use env_kindtypes,      only: si, dp, si_x, dp_x

                        implicit none

                        !http://www.netlib.org/lapack/explore-html/d1/d7e/
                        !group__double_g_esing.html#ga84fdf22a62b12ff364621e4713ce02f2
                        interface
                                pure subroutine dgesvd(jobu, jobvt, m, n, A, lda, S, U, ldu, VT, &
                                                                ldvt, work, lwork, info)
                                        import si_x, dp_x
                                        character,      intent(in   )   :: jobu, jobvt
                                        integer (si_x), intent(in   )   :: m,n,lda,ldu,lwork,ldvt
                                        integer (si_x), intent(  out)   :: info
                                        real    (dp_x), intent(inout)   :: A(lda,*)
                                        real    (dp_x), intent(  out)   :: work(*), S(*), U(ldu,*), &
                                                                                VT(ldvt,*)
                                endsubroutine dgesvd
                        end interface

                        !!! Dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        real    (dp),intent(in   )                              :: inputMatrix(:,:)
                        real    (dp),intent(inout),         allocatable         :: singular(:)
                        real    (dp),intent(inout),optional,allocatable,target  :: left(:,:)
                        real    (dp),intent(inout),optional,allocatable,target  :: right_trans(:,:)
                        !!! End dummy arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        real    (dp),target        :: dummyMatrix(1,1)

                        real    (dp),allocatable        :: inputMatrix_(:,:)
                        real    (dp),pointer,contiguous :: left_pointer(:,:)
                        real    (dp),pointer,contiguous :: right_pointer(:,:)

                        !lapack variables
                        real    (dp),allocatable   :: work(:)         !! for lapack routine
                        integer (si)               :: lwork

                        character                  :: jobu, jobvt
                        integer (si)               :: m, n

                        integer (si)               :: LDA, LDU, LDVT
                        integer (si)               :: info

                        allocate(inputMatrix_,source=inputMatrix)

                        m    = size(inputMatrix,1)
                        n    = size(inputMatrix,2)
                        LDA  = size(inputMatrix,1)
                        LDU  = size(inputMatrix,1)
                        LDVT = size(inputMatrix,2)

                        if (present(left)) then
                                jobu = 'A'

                                if (allocated(left)) then
                                        if (.not. ((size(left,1) .eq. m) .and. &
                                                  (size(left,2) .eq. m))) then
                                                deallocate(left)
                                                allocate(left(m,m))
                                        endif
                                else
                                        allocate(left(m,m))
                                endif

                                left_pointer => left
                        else
                                jobu = 'N'
                                left_pointer => dummyMatrix
                        endif

                        if (present(right_trans)) then
                                jobvt = 'A'

                                if (allocated(right_trans)) then
                                        if (.not. ((size(right_trans,1) .eq. n) .and. &
                                                  (size(right_trans,2) .eq. n))) then
                                                deallocate(right_trans)
                                                allocate(right_trans(n,n))
                                        endif
                                else
                                        allocate(right_trans(n,n))
                                endif

                                right_pointer => right_trans
                        else
                                jobvt = 'N'
                                right_pointer => dummyMatrix
                        endif

                        if (allocated(singular)) then
                                if (.not. (size(singular,1) .eq. min(m,n))) then
                                        deallocate(singular)
                                        allocate(singular(min(m,n)))
                                endif
                        else
                                allocate(singular(min(m,n)))
                        endif

                        !workspace query
                        allocate(work(1))
                        call dgesvd(jobu,jobvt,m,n,&
                                        inputMatrix_,LDA,&
                                        singular,&
                                        left_pointer,m,&
                                        right_pointer,n,&
                                        work,-1,info)

                        lwork = nint(work(1))
                        deallocate(work)
                        allocate(work(lwork))

                        !production call
                        call dgesvd(jobu,jobvt,m,n,&
                                        inputMatrix_,LDA,&
                                        singular,&
                                        left_pointer,m,&
                                        right_pointer,n,&
                                        work,lwork,info)

                 endsubroutine singular_value_decomp_dp
endmodule core_matrix

!unit test

!program main
!
!        use env_kindtypes
!        use core_matrix
!
!        implicit none
!
!        !yes, this is ugly. but it works as a simple unit test. Results were derived in R.
!
!        real    (dp)       :: array1(10,10)  = reshape( &
!                        [0.231376278912649_dp, 0.693346086423844_dp, 0.160741497762501_dp, &
!        0.746368581661955_dp, 0.68386468430981_dp, 0.741752249421552_dp, 0.778853614814579_dp, &
!        0.877396091585979_dp, 0.972667377907783_dp, 0.493809894658625_dp, 0.0548648966941983_dp, &
!        0.898453368572518_dp, 0.270436847116798_dp, 0.587363301077858_dp, 0.747960405657068_dp, &
!        0.875459946459159_dp, 0.212156609399244_dp, 0.987908352166414_dp, 0.0394798037596047_dp, &
!        0.428808676544577_dp, 0.81942753167823_dp, 0.16654721647501_dp, 0.509482031222433_dp, &
!        0.445587576832622_dp, 0.108940435806289_dp, 0.320269274292514_dp, 0.00429537380114198_dp, &
!        0.127538143889979_dp, 0.622021368471906_dp, 0.0621072174981236_dp, 0.0540236292872578_dp, &
!        0.0169521861243993_dp, 0.923818471143022_dp, 0.788286265917122_dp, 0.864530677441508_dp, &
!        0.48182024131529_dp, 0.141209222143516_dp, 0.363187932409346_dp, 0.786711658118293_dp, &
!        0.884133638814092_dp, 0.549317098455504_dp, 0.146840059198439_dp, 0.0287633682601154_dp, &
!        0.0379053701180965_dp, 0.832561332266778_dp, 0.964146455517039_dp, 0.544074802426621_dp, &
!        0.359897785121575_dp, 0.547797162784263_dp, 0.175326750148088_dp, 0.357881496893242_dp, &
!        0.959622764727101_dp, 0.638440294889733_dp, 0.942666088230908_dp, 0.352900360012427_dp, &
!        0.399023131700233_dp, 0.727050063433126_dp, 0.456970115657896_dp, 0.846911886008456_dp, &
!        0.9209067102056_dp, 0.965000302530825_dp, 0.210625957464799_dp, 0.39397593983449_dp, &
!        0.201783254742622_dp, 0.42156618530862_dp, 0.368297443957999_dp, 0.536969932727516_dp, &
!        0.0602388503029943_dp, 0.188702642219141_dp, 0.69248143886216_dp, 0.661275706021115_dp, &
!        0.802070792997256_dp, 0.741438026772812_dp, 0.81903372448869_dp, 0.343384377192706_dp, &
!        0.6662566002924_dp, 0.750419952441007_dp, 0.304543149191886_dp, 0.138954193564132_dp, &
!        0.510839792899787_dp, 0.406359357060865_dp, 0.123135861009359_dp, 0.938722611404955_dp, &
!        0.950078680878505_dp, 0.744901102036238_dp, 0.456361019518226_dp, 0.597216539783403_dp, &
!        0.797275994671509_dp, 0.541320266202092_dp, 0.44853696343489_dp, 0.758698848774657_dp, &
!        0.140811997931451_dp, 0.887924589216709_dp, 0.838451216463_dp, 0.543776330072433_dp, &
!        0.233618144877255_dp, 0.653859071433544_dp, 0.542026714654639_dp, 0.731736903078854_dp, &
!        0.0315509843640029_dp], [10,10])
!        real    (dp)       :: symm_array1(10,10)  = reshape( &
!                        [2.72330747358501_dp, 0.778620538767427_dp, 1.25721962377429_dp, &
!        0.166886705905199_dp, 1.38195694913156_dp, 1.12927551520988_dp, 0.711528190411627_dp, &
!        0.98657527891919_dp, 1.10746893333271_dp, 0.875894732307643_dp, 0.778620538767427_dp, &
!        2.92687488859519_dp, 1.70472226897255_dp, 0.977211354067549_dp, 0.149245712673292_dp, &
!        1.25045907567255_dp, 0.840626067714766_dp, 0.972975023789331_dp, 0.852923363447189_dp, &
!        0.739653498400003_dp, 1.25721962377429_dp, 1.70472226897255_dp, 2.85886506922543_dp, &
!        1.55976403784007_dp, 0.191622453741729_dp, 0.731647792970762_dp, 0.662440445739776_dp, &
!        1.23254952859133_dp, 0.273822978371754_dp, 0.896642932668328_dp, 0.166886705905199_dp, &
!        0.977211354067549_dp, 1.55976403784007_dp, 3.24274884117767_dp, 0.750555692240596_dp, &
!        0.124686117051169_dp, 1.00826012645848_dp, 0.725064822006971_dp, 1.11448277416639_dp, &
!        1.21676548290998_dp, 1.38195694913156_dp, 0.149245712673292_dp, 0.191622453741729_dp, &
!        0.750555692240596_dp, 2.12951065739617_dp, 0.482213373063132_dp, 0.621653700247407_dp, &
!        1.75129548530094_dp, 1.39098125230521_dp, 0.706736625405028_dp, 1.12927551520988_dp, &
!        1.25045907567255_dp, 0.731647792970762_dp, 0.124686117051169_dp, 0.482213373063132_dp, &
!        3.3759339354001_dp, 0.782613976858556_dp, 1.59810881828889_dp, 1.30673734867014_dp, &
!        0.621303273132071_dp, 0.711528190411627_dp, 0.840626067714766_dp, 0.662440445739776_dp, &
!        1.00826012645848_dp, 0.621653700247407_dp, 0.782613976858556_dp, 3.41408031992614_dp, &
!        0.647184573113918_dp, 1.52004841901362_dp, 0.254456369904801_dp, 0.98657527891919_dp, &
!        0.972975023789331_dp, 1.23254952859133_dp, 0.725064822006971_dp, 1.75129548530094_dp, &
!        1.59810881828889_dp, 0.647184573113918_dp, 2.79017546027899_dp, 1.42070761485957_dp, &
!        0.796630741329864_dp, 1.10746893333271_dp, 0.852923363447189_dp, 0.273822978371754_dp, &
!        1.11448277416639_dp, 1.39098125230521_dp, 1.30673734867014_dp, 1.52004841901362_dp, &
!        1.42070761485957_dp, 2.44817775348201_dp, 1.07244482659735_dp, 0.875894732307643_dp, &
!        0.739653498400003_dp, 0.896642932668328_dp, 1.21676548290998_dp, 0.706736625405028_dp, &
!        0.621303273132071_dp, 0.254456369904801_dp, 0.796630741329864_dp, 1.07244482659735_dp, &
!        3.82914637541398_dp], [10,10])
!        real    (dp)       :: symm_array2(10,10)  = reshape( &
!        [4.61698451275729_dp, 0.0292803538594762_dp,&
!        1.05197120932541_dp, 0.684688056121411_dp,&
!        -1.28562942551886_dp, -0.721686073481445_dp,&
!        3.1229627538451_dp, 1.55058048087287_dp,&
!        -0.0473611148721462_dp, 2.79081630236308_dp,&
!        0.0292803538594762_dp, 9.08864873310607_dp,&
!        0.315996228756272_dp, -0.3358660843886_dp,&
!        -3.86608803678688_dp, -0.699104984523349_dp,&
!        0.262858902703926_dp, 0.364713663743865_dp,&
!        -0.00874620991940618_dp, -4.21964397508402_dp,&
!        1.05197120932541_dp, 0.315996228756272_dp,&
!        5.27266256831365_dp, 0.0609474408552645_dp,&
!        0.313958036314005_dp, 0.549243308305777_dp,&
!        1.72968776339159_dp, 0.683551557356109_dp,&
!        1.26408725872342_dp, -1.11952292913463_dp,&
!        0.684688056121411_dp, -0.3358660843886_dp,&
!        0.0609474408552645_dp, 6.81536084843332_dp,&
!        -2.33530070104547_dp, 2.32078921670246_dp,&
!        1.13348185036571_dp, 2.10717762953615_dp,&
!        -0.632557918910687_dp, 0.401051718202871_dp,&
!        -1.28562942551886_dp, -3.86608803678688_dp,&
!        0.313958036314005_dp, -2.33530070104547_dp,&
!        5.57058903440522_dp, 0.330311560139098_dp,&
!        -3.42036426186318_dp, -2.79756895819969_dp,&
!        -0.481146687362744_dp, 0.171635577925558_dp,&
!        -0.721686073481445_dp, -0.699104984523349_dp,&
!        0.549243308305777_dp, 2.32078921670246_dp,&
!        0.330311560139098_dp, 4.55033995092687_dp,&
!        -0.944702525086404_dp, 4.68786013394113_dp,&
!        -1.67123659268596_dp, 0.256951070239753_dp,&
!        3.1229627538451_dp, 0.262858902703926_dp,&
!        1.72968776339159_dp, 1.13348185036571_dp,&
!        -3.42036426186318_dp, -0.944702525086404_dp,&
!        6.4671709852884_dp, 2.53330000354259_dp,&
!        1.22915904840649_dp, 0.940181348209596_dp,&
!        1.55058048087287_dp, 0.364713663743865_dp,&
!        0.683551557356109_dp, 2.10717762953615_dp,&
!        -2.79756895819969_dp, 4.68786013394113_dp,&
!        2.53330000354259_dp, 12.7740384490773_dp,&
!        -2.77617942518654_dp, 1.90268001626435_dp,&
!        -0.0473611148721462_dp, -0.00874620991940618_dp,&
!        1.26408725872342_dp, -0.632557918910687_dp,&
!        -0.481146687362744_dp, -1.67123659268596_dp,&
!        1.22915904840649_dp, -2.77617942518654_dp,&
!        2.2474803967473_dp, -0.823167365477737_dp,&
!        2.79081630236308_dp, -4.21964397508402_dp,&
!        -1.11952292913463_dp, 0.401051718202871_dp,&
!        0.171635577925558_dp, 0.256951070239753_dp,&
!        0.940181348209596_dp, 1.90268001626435_dp,&
!        -0.823167365477737_dp, 5.72623257232743_dp], [10,10])
!        real    (dp)       :: array2(10)    = &
!                        [0.554870725376531_dp, 0.999600998358801_dp, 0.925282047828659_dp, &
!        0.0222018531057984_dp, 0.709053714759648_dp, 0.157962263328955_dp, 0.961079649860039_dp, &
!        0.490833985153586_dp, 0.90480772475712_dp, 0.448791508330032_dp]
!        real    (dp)       :: linearsystem_answer(10)    = &
!                        [ -7.08315230128404_dp, 1.09714160983019_dp, -0.75419920275508_dp, &
!        -2.96988094568714_dp, &
!        4.96835477574137_dp, 7.78513857623552_dp, -1.77575087869568_dp, -3.98011759268542_dp, &
!        3.00859758803731_dp, 0.68489503313388_dp ]
!        real    (dp)       :: eigenvalues_answer(10) = &
!                        [ -0.0613118396280499_dp, 0.439347456217176_dp, 1.17272160139944_dp, &
!        1.27202809968562_dp, 2.13957292629257_dp, 2.48763498121484_dp, 3.38370377046507_dp, &
!        3.49800818887043_dp, 4.08477086611459_dp, 11.322344723849_dp ]
!        real    (dp),allocatable       :: array3(:),array4(:,:),array5(:,:)
!
!        real    (dp)                   :: tol = 1e-12
!        real    (dp)                   :: diff, diff2
!
!        !lapack routine variables
!        real    (dp),allocatable     :: work(:)
!        real    (dp)                 :: workQuery(1)
!        integer (si)                 :: info
!        integer (si)                 :: lwork
!        integer (si)                 :: i
!
!        integer (si)                 :: nAtoms
!        integer (si)                 :: nCg
!
!        print*, "Linear system vector test:"
!
!        array3=lsolve(array1,array2)
!        diff = sum(abs(array3-linearsystem_answer))/size(array3)
!
!        if (diff .lt. tol) then
!                print*, ':: Test passed. Difference: ', diff
!        else
!                print*, ':: Test failed. Difference: ', diff
!        endif
!
!        print*, "Eigendecomposition test:"
!
!        call spectrum(symm_array1,array3,array4)
!        diff  = sum(abs(array3-eigenvalues_answer))/size(array3)
!        diff2 = sum(abs(matmul(array4,transpose(array4))-diagonal(1.0_dp,size(array4,1))))/size(array4,1)
!
!        if ((diff .lt. tol) .and. (diff2 .lt. tol)) then
!                print*, ':: Test passed. Eigenvalue difference: ', diff
!                print*, ':: Test passed. Eigenvector difference: ', diff2
!        else
!                print*, ':: Test failed. Eigenvalue difference: ', diff
!                print*, ':: Test failed. Eigenvector difference: ', diff2
!                print*, eigenvalues_answer
!                print*, array3
!        endif
!
!        call spectrum(M=symm_array1,&
!                      B=symm_array2,&
!                      d=array3,&
!                      Q=array4,&
!                      num_bottom_eigenvalues=3)
!
!        print*, "d", array3
!        print*, "Q", array4
!
!        !if ((diff .lt. tol) .and. (diff2 .lt. tol)) then
!        !        print*, ':: Test passed. Eigenvalue difference: ', diff
!        !        print*, ':: Test passed. Eigenvector difference: ', diff2
!        !else
!        !        print*, ':: Test failed. Eigenvalue difference: ', diff
!        !        print*, ':: Test failed. Eigenvector difference: ', diff2
!        !        print*, eigenvalues_answer
!        !        print*, array3
!        !endif
!
!        !array4 = quadForm(symm_array1,array1,rightTranspose=.true.)
!        !print*,array4
!
!        !array4 = quadForm(symm_array1,array1,rightTranspose=.false.)
!        !print*,array4
!!
!!        print*, array1
!!        print*, array3
!!        print*, array5
!
!endprogram main
