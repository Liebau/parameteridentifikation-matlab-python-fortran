program test_reference
    ! Fortran equivalents of the dense-reference tests in
    ! python/tests/test_model.py.  The low-rank tests remain Python/MATLAB-only.
    use gwr_precision, only: dp
    use gwr_model, only: problem_data, build_problem, destroy_problem
    use gwr_linear_algebra, only: adjoint_weighted, adjoint_euclidean
    use gwr_evaluation, only: analysis_result, analyze
    implicit none

    integer :: failures

    failures = 0
    call test_trapezoidal_problem(failures)
    call test_stable_direct_case(failures)
    call test_weighted_regularization(failures)
    call test_adjoint_conventions(failures)
    call test_problem_lifecycle(failures)

    if (failures /= 0) then
        write(*, '(I0,A)') failures, ' Fortran reference test(s) failed.'
        stop 1
    end if
    write(*, '(A)') 'All Fortran reference tests passed.'

contains

    subroutine test_trapezoidal_problem(failure_count)
        ! Python: test_builds_trapezoidal_matrix.
        integer, intent(inout) :: failure_count
        type(problem_data) :: problem
        real(dp) :: expected_x(0:2), expected_weights(0:2)
        real(dp) :: expected_matrix(0:2, 0:2)
        integer :: i, j, info

        call build_problem(2, problem, info)
        call assert_true(info == 0, 'n=2 problem construction', failure_count)
        if (info /= 0) return

        expected_x = (/ 0.0_dp, 0.5_dp, 1.0_dp /)
        expected_weights = (/ 0.25_dp, 0.5_dp, 0.25_dp /)
        do j = 0, 2
            do i = 0, 2
                expected_matrix(i, j) = expected_weights(j) * &
                    exp(expected_x(i) * expected_x(j))
            end do
        end do
        call assert_close(maxval(abs(problem%x - expected_x)), 0.0_dp, &
                          1.0e-15_dp, 'n=2 nodes', failure_count)
        call assert_close(maxval(abs(problem%weights - expected_weights)), &
                          0.0_dp, 1.0e-15_dp, 'n=2 weights', failure_count)
        call assert_close(maxval(abs(problem%matrix - expected_matrix)), &
                          0.0_dp, 1.0e-15_dp, 'n=2 matrix', failure_count)
        call destroy_problem(problem)
    end subroutine test_trapezoidal_problem

    subroutine test_stable_direct_case(failure_count)
        ! Python: test_reproduces_stable_direct_case.
        integer, intent(inout) :: failure_count
        type(analysis_result) :: result
        integer :: info
        character(len=256) :: message

        call analyze(4, 0.0_dp, adjoint_weighted, result, info, message)
        call assert_true(info == 0, 'n=4 direct solve: '//trim(message), &
                         failure_count)
        if (info /= 0) return
        call assert_close(result%condition, 9.100210765820118e5_dp, &
                          2.0e-7_dp, 'n=4 condition', failure_count, .true.)
        call assert_close(result%direct_error, 1.089640687693824_dp, &
                          2.0e-9_dp, 'n=4 direct error', failure_count, .true.)
    end subroutine test_stable_direct_case

    subroutine test_weighted_regularization(failure_count)
        ! Python: test_weighted_regularization_improves_fine_grid.
        integer, intent(inout) :: failure_count
        type(analysis_result) :: result
        integer :: info
        character(len=256) :: message

        call analyze(32, 1.0e-6_dp, adjoint_weighted, result, info, message)
        call assert_true(info == 0, 'weighted regularization: '//trim(message), &
                         failure_count)
        if (info /= 0) return
        call assert_close(result%regularized_error, &
                          1.641158222633532e-2_dp, 2.0e-7_dp, &
                          'n=32 weighted regularized error', failure_count, &
                          .true.)
        call assert_true(result%regularized_error < 2.0e-2_dp, &
                         'weighted regularization accuracy', failure_count)
    end subroutine test_weighted_regularization

    subroutine test_adjoint_conventions(failure_count)
        ! Python: test_adjoint_conventions_remain_explicit.
        integer, intent(inout) :: failure_count
        type(analysis_result) :: weighted, euclidean
        integer :: info
        character(len=256) :: message

        call analyze(32, 1.0e-3_dp, adjoint_weighted, weighted, info, message)
        call assert_true(info == 0, 'weighted adjoint: '//trim(message), &
                         failure_count)
        call analyze(32, 1.0e-3_dp, adjoint_euclidean, euclidean, info, message)
        call assert_true(info == 0, 'euclidean adjoint: '//trim(message), &
                         failure_count)
        if (info /= 0) return
        call assert_true(weighted%regularized_error < 6.0e-2_dp, &
                         'weighted adjoint expected error', failure_count)
        call assert_true(euclidean%regularized_error > 1.0_dp, &
                         'euclidean adjoint differs clearly', failure_count)
        call assert_true(euclidean%regularized_error > &
                         10.0_dp * weighted%regularized_error, &
                         'adjoint conventions remain distinct', failure_count)
    end subroutine test_adjoint_conventions

    subroutine test_problem_lifecycle(failure_count)
        ! Exercise pointer ownership repeatedly and verify the error path.
        integer, intent(inout) :: failure_count
        type(problem_data) :: problem
        integer :: iteration, info

        call build_problem(0, problem, info)
        call assert_true(info == -1, 'reject n=0', failure_count)
        call destroy_problem(problem)

        do iteration = 1, 20
            call build_problem(8, problem, info)
            call assert_true(info == 0, 'repeated problem construction', &
                             failure_count)
            if (info /= 0) return
            call assert_true(associated(problem%matrix), &
                             'matrix associated after construction', &
                             failure_count)
            call destroy_problem(problem)
            call assert_true(.not. associated(problem%matrix), &
                             'matrix disassociated after destruction', &
                             failure_count)
        end do
    end subroutine test_problem_lifecycle

    subroutine assert_true(condition, label, failure_count)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failure_count

        if (.not. condition) then
            failure_count = failure_count + 1
            write(*, '(A,A)') 'FAIL: ', trim(label)
        end if
    end subroutine assert_true

    subroutine assert_close(actual, expected, tolerance, label, &
                            failure_count, relative)
        real(dp), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failure_count
        logical, intent(in), optional :: relative
        real(dp) :: difference, scale

        difference = abs(actual - expected)
        scale = 1.0_dp
        if (present(relative)) then
            if (relative) scale = max(abs(expected), 1.0_dp)
        end if
        if (difference > tolerance * scale) then
            failure_count = failure_count + 1
            write(*, '(A,A,2(A,ES24.16))') 'FAIL: ', trim(label), &
                ', actual=', actual, ', expected=', expected
        end if
    end subroutine assert_close

end program test_reference

