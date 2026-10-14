module gwr_performance
    ! Phase-separated timings equivalent to measure() in
    ! python/src/gwr_parameteridentifikation/benchmark.py and to
    ! gwrparameter.measure_performance() in MATLAB.
    use gwr_precision, only: dp
    use gwr_model, only: problem_data, build_problem, destroy_problem
    use gwr_linear_algebra, only: solve_linear_system, &
        build_tikhonov_system, adjoint_weighted
    implicit none
    private

    ! Non-default SYSTEM_CLOCK kinds are a Fortran-2003 extension supported by
    ! gfortran.  They are required here for wall-clock resolution below 1 ms.
    ! The mathematical implementation itself remains strict Fortran 95.
    integer, parameter :: clock_kind = selected_int_kind(18)

    type, public :: benchmark_result
        integer :: n = 0
        real(dp) :: beta = 0.0_dp
        real(dp) :: assembly_seconds = 0.0_dp
        real(dp) :: direct_solve_seconds = 0.0_dp
        real(dp) :: tikhonov_setup_seconds = 0.0_dp
        real(dp) :: tikhonov_solve_seconds = 0.0_dp
    end type benchmark_result

    public :: measure_performance

contains

    subroutine measure_performance(n, beta, repeats, result, info, message)
        integer, intent(in) :: n, repeats
        real(dp), intent(in) :: beta
        type(benchmark_result), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(problem_data) :: problem, temporary_problem
        real(dp), allocatable :: samples(:), solution(:), lhs(:, :), rhs(:)
        integer :: repeat_index, stat
        integer(kind=clock_kind) :: start_count, end_count, count_rate

        info = 0
        message = ''
        result%n = n
        result%beta = beta
        if (n < 1 .or. repeats < 1 .or. beta <= 0.0_dp) then
            info = -1
            message = 'n, beta and repeats must be positive'
            return
        end if

        call system_clock(count_rate=count_rate)
        if (count_rate <= 0_clock_kind) then
            info = -2
            message = 'system_clock has no usable rate'
            return
        end if

        allocate(samples(repeats), solution(n + 1), lhs(n + 1, n + 1), &
                 rhs(n + 1), stat=stat)
        if (stat /= 0) then
            info = stat
            message = 'benchmark workspace allocation failed'
            return
        end if

        call build_problem(n, problem, info)
        if (info /= 0) then
            message = 'reference problem construction failed'
            return
        end if
        call build_tikhonov_system(problem%matrix, problem%data, beta, &
                                   adjoint_weighted, lhs, rhs, info)
        if (info /= 0) then
            message = 'initial Tikhonov setup failed'
            call destroy_problem(problem)
            return
        end if

        ! Warm-up and timings for problem assembly.  Destruction is included,
        ! matching the lifetime of the temporary arrays returned by Python.
        call build_problem(n, temporary_problem, info)
        if (info == 0) call destroy_problem(temporary_problem)
        if (info /= 0) then
            message = 'assembly warm-up failed'
            call destroy_problem(problem)
            return
        end if
        do repeat_index = 1, repeats
            call system_clock(start_count)
            call build_problem(n, temporary_problem, info)
            if (info == 0) call destroy_problem(temporary_problem)
            call system_clock(end_count)
            if (info /= 0) then
                message = 'timed problem assembly failed'
                call destroy_problem(problem)
                return
            end if
            samples(repeat_index) = elapsed_seconds( &
                start_count, end_count, count_rate)
        end do
        result%assembly_seconds = median(samples)

        ! DGESV copies its inputs inside solve_linear_system, just as
        ! numpy.linalg.solve preserves the arrays passed by the caller.
        call solve_linear_system(problem%matrix, problem%data, solution, info)
        if (info /= 0) then
            message = 'direct-solve warm-up failed'
            call destroy_problem(problem)
            return
        end if
        do repeat_index = 1, repeats
            call system_clock(start_count)
            call solve_linear_system( &
                problem%matrix, problem%data, solution, info)
            call system_clock(end_count)
            if (info /= 0) then
                message = 'timed direct solve failed'
                call destroy_problem(problem)
                return
            end if
            samples(repeat_index) = elapsed_seconds( &
                start_count, end_count, count_rate)
        end do
        result%direct_solve_seconds = median(samples)

        ! Allocate temporary outputs for every setup, corresponding to the
        ! newly returned lhs/rhs arrays in MATLAB and Python.
        call build_tikhonov_allocating(problem, beta, info)
        if (info /= 0) then
            message = 'Tikhonov-setup warm-up failed'
            call destroy_problem(problem)
            return
        end if
        do repeat_index = 1, repeats
            call system_clock(start_count)
            call build_tikhonov_allocating(problem, beta, info)
            call system_clock(end_count)
            if (info /= 0) then
                message = 'timed Tikhonov setup failed'
                call destroy_problem(problem)
                return
            end if
            samples(repeat_index) = elapsed_seconds( &
                start_count, end_count, count_rate)
        end do
        result%tikhonov_setup_seconds = median(samples)

        call solve_linear_system(lhs, rhs, solution, info)
        if (info /= 0) then
            message = 'Tikhonov-solve warm-up failed'
            call destroy_problem(problem)
            return
        end if
        do repeat_index = 1, repeats
            call system_clock(start_count)
            call solve_linear_system(lhs, rhs, solution, info)
            call system_clock(end_count)
            if (info /= 0) then
                message = 'timed Tikhonov solve failed'
                call destroy_problem(problem)
                return
            end if
            samples(repeat_index) = elapsed_seconds( &
                start_count, end_count, count_rate)
        end do
        result%tikhonov_solve_seconds = median(samples)

        call destroy_problem(problem)
    end subroutine measure_performance

    subroutine build_tikhonov_allocating(problem, beta, info)
        type(problem_data), intent(in) :: problem
        real(dp), intent(in) :: beta
        integer, intent(out) :: info
        real(dp), allocatable :: temporary_lhs(:, :), temporary_rhs(:)
        integer :: size_n, stat

        size_n = problem%n + 1
        allocate(temporary_lhs(size_n, size_n), temporary_rhs(size_n), &
                 stat=stat)
        if (stat /= 0) then
            info = stat
            return
        end if
        call build_tikhonov_system(problem%matrix, problem%data, beta, &
                                   adjoint_weighted, temporary_lhs, &
                                   temporary_rhs, info)
    end subroutine build_tikhonov_allocating

    real(dp) function elapsed_seconds(start_count, end_count, count_rate)
        integer(kind=clock_kind), intent(in) :: start_count, end_count
        integer(kind=clock_kind), intent(in) :: count_rate

        elapsed_seconds = real(end_count - start_count, dp) / &
            real(count_rate, dp)
    end function elapsed_seconds

    real(dp) function median(values)
        real(dp), intent(in) :: values(:)
        real(dp), allocatable :: sorted(:)
        real(dp) :: current
        integer :: i, j, number

        number = size(values)
        allocate(sorted(number))
        sorted = values
        do i = 2, number
            current = sorted(i)
            j = i - 1
            do while (j >= 1)
                if (sorted(j) <= current) exit
                sorted(j + 1) = sorted(j)
                j = j - 1
            end do
            sorted(j + 1) = current
        end do
        if (mod(number, 2) == 0) then
            median = 0.5_dp * (sorted(number / 2) + &
                               sorted(number / 2 + 1))
        else
            median = sorted((number + 1) / 2)
        end if
    end function median

end module gwr_performance

