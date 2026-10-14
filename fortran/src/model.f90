module gwr_model
    ! Discrete model corresponding to Problem and build_problem() in
    ! python/src/gwr_parameteridentifikation/model.py.
    use gwr_precision, only: dp
    implicit none
    private

    type, public :: problem_data
        integer :: n = 0
        ! Pointer components keep this data type compatible with Fortran 95.
        ! destroy_problem() owns the corresponding deallocation.
        real(dp), pointer :: x(:) => null()
        real(dp), pointer :: weights(:) => null()
        real(dp), pointer :: matrix(:, :) => null()
        real(dp), pointer :: data(:) => null()
        real(dp), pointer :: exact_lambda(:) => null()
    end type problem_data

    public :: build_problem
    public :: destroy_problem

contains

    subroutine build_problem(n, problem, info)
        ! Build exactly the same nodes, trapezoidal weights, matrix and
        ! reference vectors as Python build_problem(n).  Bounds 0:n mirror
        ! the mathematical indices and Python's entries 0,...,n.
        integer, intent(in) :: n
        type(problem_data), intent(out) :: problem
        integer, intent(out) :: info
        integer :: i, j, stat

        info = 0
        nullify(problem%x, problem%weights, problem%matrix, problem%data, &
                problem%exact_lambda)
        if (n < 1) then
            info = -1
            return
        end if

        problem%n = n
        allocate(problem%x(0:n), problem%weights(0:n), &
                 problem%matrix(0:n, 0:n), problem%data(0:n), &
                 problem%exact_lambda(0:n), stat=stat)
        if (stat /= 0) then
            info = stat
            call destroy_problem(problem)
            return
        end if

        do i = 0, n
            problem%x(i) = real(i, dp) / real(n, dp)
            problem%weights(i) = 1.0_dp / real(n, dp)
        end do
        problem%weights(0) = 1.0_dp / (2.0_dp * real(n, dp))
        problem%weights(n) = problem%weights(0)

        do j = 0, n
            do i = 0, n
                ! Python: np.exp(np.outer(x, x)) * weights[np.newaxis, :]
                problem%matrix(i, j) = problem%weights(j) * &
                    exp(problem%x(i) * problem%x(j))
            end do
        end do

        do i = 0, n
            problem%data(i) = (exp(problem%x(i) + 1.0_dp) - 1.0_dp) / &
                (problem%x(i) + 1.0_dp)
            problem%exact_lambda(i) = exp(problem%x(i))
        end do
    end subroutine build_problem

    subroutine destroy_problem(problem)
        ! Allocatable components are a Fortran 2003 feature.  The explicit
        ! cleanup is therefore part of the strict Fortran-95 implementation.
        type(problem_data), intent(inout) :: problem

        if (associated(problem%x)) deallocate(problem%x)
        if (associated(problem%weights)) deallocate(problem%weights)
        if (associated(problem%matrix)) deallocate(problem%matrix)
        if (associated(problem%data)) deallocate(problem%data)
        if (associated(problem%exact_lambda)) deallocate(problem%exact_lambda)
        nullify(problem%x, problem%weights, problem%matrix, problem%data, &
                problem%exact_lambda)
        problem%n = 0
    end subroutine destroy_problem

end module gwr_model

