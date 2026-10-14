module gwr_linear_algebra
    ! BLAS/LAPACK counterparts of the NumPy operations used in model.py:
    ! np.linalg.solve, np.linalg.svd, matrix @ matrix and matrix @ vector.
    use gwr_precision, only: dp
    implicit none
    private

    integer, parameter, public :: adjoint_weighted = 1
    integer, parameter, public :: adjoint_euclidean = 2

    public :: solve_linear_system
    public :: condition_number_2
    public :: build_tikhonov_system

    interface
        ! Explicit Fortran-77 library interfaces provide compile-time argument
        ! checking while retaining compatibility with a Fortran-95 caller.
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
            double precision, intent(inout) :: a(lda, *)
            double precision, intent(inout) :: b(ldb, *)
        end subroutine dgesv

        subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, &
                          ldvt, work, lwork, info)
            character(len=1), intent(in) :: jobu, jobvt
            integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
            integer, intent(out) :: info
            double precision, intent(inout) :: a(lda, *)
            double precision, intent(out) :: s(*)
            double precision, intent(out) :: u(ldu, *)
            double precision, intent(out) :: vt(ldvt, *)
            double precision, intent(inout) :: work(*)
        end subroutine dgesvd

        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, &
                         ldb, beta, c, ldc)
            character(len=1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            double precision, intent(in) :: alpha, beta
            double precision, intent(in) :: a(lda, *), b(ldb, *)
            double precision, intent(inout) :: c(ldc, *)
        end subroutine dgemm

        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, &
                         y, incy)
            character(len=1), intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            double precision, intent(in) :: alpha, beta
            double precision, intent(in) :: a(lda, *), x(*)
            double precision, intent(inout) :: y(*)
        end subroutine dgemv
    end interface

contains

    subroutine solve_linear_system(matrix, rhs, solution, info)
        ! Python counterpart: np.linalg.solve(matrix, rhs).
        ! DGESV overwrites both matrix and right-hand side, hence the copies.
        real(dp), intent(in) :: matrix(:, :)
        real(dp), intent(in) :: rhs(:)
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: info
        real(dp), allocatable :: matrix_copy(:, :), rhs_copy(:, :)
        integer, allocatable :: pivots(:)
        integer :: size_n, stat

        info = 0
        size_n = size(rhs)
        if (size(matrix, 1) /= size_n .or. size(matrix, 2) /= size_n) then
            info = -100
            return
        end if

        if (size(solution) /= size_n) then
            info = -101
            return
        end if
        allocate(matrix_copy(size_n, size_n), rhs_copy(size_n, 1), &
                 pivots(size_n), stat=stat)
        if (stat /= 0) then
            info = stat
            return
        end if
        matrix_copy = matrix
        rhs_copy(:, 1) = rhs

        call dgesv(size_n, 1, matrix_copy, size_n, pivots, rhs_copy, &
                   size_n, info)
        if (info == 0) solution = rhs_copy(:, 1)
    end subroutine solve_linear_system

    subroutine condition_number_2(matrix, condition, info)
        ! Python counterpart:
        ! singular_values = np.linalg.svd(matrix, compute_uv=False)
        ! condition = singular_values[0] / singular_values[-1]
        real(dp), intent(in) :: matrix(:, :)
        real(dp), intent(out) :: condition
        integer, intent(out) :: info
        real(dp), allocatable :: matrix_copy(:, :), singular_values(:)
        real(dp), allocatable :: work(:)
        real(dp) :: workspace_query(1), unused_u(1, 1), unused_vt(1, 1)
        integer :: rows, columns, lwork, stat

        info = 0
        condition = 0.0_dp
        rows = size(matrix, 1)
        columns = size(matrix, 2)
        allocate(matrix_copy(rows, columns), &
                 singular_values(min(rows, columns)), stat=stat)
        if (stat /= 0) then
            info = stat
            return
        end if
        matrix_copy = matrix

        ! LAPACK workspace query.  The first DGESVD call performs no SVD.
        lwork = -1
        call dgesvd('N', 'N', rows, columns, matrix_copy, rows, &
                    singular_values, unused_u, 1, unused_vt, 1, &
                    workspace_query, lwork, info)
        if (info /= 0) return

        lwork = max(1, int(workspace_query(1)))
        allocate(work(lwork), stat=stat)
        if (stat /= 0) then
            info = stat
            return
        end if
        matrix_copy = matrix
        call dgesvd('N', 'N', rows, columns, matrix_copy, rows, &
                    singular_values, unused_u, 1, unused_vt, 1, work, &
                    lwork, info)
        if (info /= 0) return

        if (singular_values(size(singular_values)) <= 0.0_dp) then
            condition = huge(1.0_dp)
        else
            condition = singular_values(1) / &
                singular_values(size(singular_values))
        end if
    end subroutine condition_number_2

    subroutine build_tikhonov_system(matrix, data, beta, adjoint_mode, &
                                      lhs, rhs, info)
        ! Python counterpart: the beta > 0 branch of analyze().  Keeping the
        ! two adjoints in separate BLAS calls prevents an accidental change of
        ! the mathematical model during a later migration.
        real(dp), intent(in) :: matrix(:, :), data(:), beta
        integer, intent(in) :: adjoint_mode
        real(dp), intent(out) :: lhs(:, :), rhs(:)
        integer, intent(out) :: info
        integer :: i, size_n

        info = 0
        size_n = size(data)
        if (size(matrix, 1) /= size_n .or. size(matrix, 2) /= size_n) then
            info = -100
            return
        end if
        if (beta <= 0.0_dp) then
            info = -101
            return
        end if
        if (adjoint_mode /= adjoint_weighted .and. &
            adjoint_mode /= adjoint_euclidean) then
            info = -102
            return
        end if

        if (size(lhs, 1) /= size_n .or. size(lhs, 2) /= size_n .or. &
            size(rhs) /= size_n) then
            info = -103
            return
        end if
        lhs = 0.0_dp
        rhs = 0.0_dp

        if (adjoint_mode == adjoint_weighted) then
            ! Consistent weighted adjoint: (beta I + A*A) alpha = A*y.
            ! Python: lhs = beta*I + A @ A; rhs = A @ data.
            call dgemm('N', 'N', size_n, size_n, size_n, 1.0_dp, &
                       matrix, size_n, matrix, size_n, 0.0_dp, lhs, size_n)
            call dgemv('N', size_n, size_n, 1.0_dp, matrix, size_n, &
                       data, 1, 0.0_dp, rhs, 1)
        else
            ! Euclidean comparison: (beta I + A^T*A) alpha = A^T*y.
            ! Python: lhs = beta*I + A.T @ A; rhs = A.T @ data.
            call dgemm('T', 'N', size_n, size_n, size_n, 1.0_dp, &
                       matrix, size_n, matrix, size_n, 0.0_dp, lhs, size_n)
            call dgemv('T', size_n, size_n, 1.0_dp, matrix, size_n, &
                       data, 1, 0.0_dp, rhs, 1)
        end if

        do i = 1, size_n
            lhs(i, i) = lhs(i, i) + beta
        end do
    end subroutine build_tikhonov_system

end module gwr_linear_algebra

