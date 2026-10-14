module gwr_evaluation
    ! High-level analysis corresponding to AnalysisResult and analyze() in
    ! python/src/gwr_parameteridentifikation/model.py.
    use gwr_precision, only: dp
    use gwr_model, only: problem_data, build_problem, destroy_problem
    use gwr_linear_algebra, only: solve_linear_system, condition_number_2, &
        build_tikhonov_system, adjoint_weighted, adjoint_euclidean
    implicit none
    private

    type, public :: analysis_result
        integer :: n = 0
        real(dp) :: beta = 0.0_dp
        integer :: adjoint_mode = adjoint_weighted
        real(dp) :: condition = 0.0_dp
        real(dp) :: direct_error = 0.0_dp
        real(dp) :: direct_relative_error = 0.0_dp
        real(dp) :: regularized_error = 0.0_dp
        real(dp) :: regularized_relative_error = 0.0_dp
    end type analysis_result

    public :: analyze

contains

    subroutine analyze(n, beta, adjoint_mode, result, info, message)
        ! Execute the common reference workflow: build the problem, determine
        ! cond_2(A), solve A*alpha=y and optionally solve one Tikhonov system.
        integer, intent(in) :: n, adjoint_mode
        real(dp), intent(in) :: beta
        type(analysis_result), intent(out) :: result
        integer, intent(out) :: info
        character(len=*), intent(out) :: message
        type(problem_data) :: problem
        real(dp) :: lhs(n + 1, n + 1), rhs(n + 1)
        real(dp) :: alpha_direct(n + 1), alpha_regularized(n + 1)
        real(dp) :: exact_norm

        info = 0
        message = ''
        result%n = n
        result%beta = beta
        result%adjoint_mode = adjoint_mode

        if (beta < 0.0_dp) then
            info = -1
            message = 'beta must be non-negative'
            return
        end if
        if (adjoint_mode /= adjoint_weighted .and. &
            adjoint_mode /= adjoint_euclidean) then
            info = -2
            message = 'invalid adjoint mode'
            return
        end if

        call build_problem(n, problem, info)
        if (info /= 0) then
            message = 'building the discrete problem failed'
            return
        end if

        ! Both language versions compute the condition before DGESV can modify
        ! a matrix copy.  This also keeps diagnostics independent of the solve.
        call condition_number_2(problem%matrix, result%condition, info)
        if (info /= 0) then
            message = 'DGESVD failed while computing the condition number'
            call destroy_problem(problem)
            return
        end if

        call solve_linear_system(problem%matrix, problem%data, &
                                 alpha_direct, info)
        if (info /= 0) then
            message = 'DGESV failed for the direct system'
            call destroy_problem(problem)
            return
        end if

        exact_norm = maxval(abs(problem%exact_lambda))
        result%direct_error = maxval(abs(alpha_direct - &
                                          problem%exact_lambda))
        result%direct_relative_error = result%direct_error / exact_norm

        if (beta > 0.0_dp) then
            ! build_tikhonov_system contains the explicit A versus transpose(A)
            ! choice; analyze itself treats both modes identically afterwards.
            call build_tikhonov_system(problem%matrix, problem%data, beta, &
                                       adjoint_mode, lhs, rhs, info)
            if (info /= 0) then
                message = 'building the Tikhonov system failed'
                call destroy_problem(problem)
                return
            end if
            call solve_linear_system(lhs, rhs, alpha_regularized, info)
            if (info /= 0) then
                message = 'DGESV failed for the Tikhonov system'
                call destroy_problem(problem)
                return
            end if
            result%regularized_error = maxval(abs( &
                alpha_regularized - problem%exact_lambda))
            result%regularized_relative_error = &
                result%regularized_error / exact_norm
        end if
        call destroy_problem(problem)
    end subroutine analyze

end module gwr_evaluation

