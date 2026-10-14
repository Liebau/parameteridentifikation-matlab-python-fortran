program run_reference
    ! Console/CSV counterpart of
    ! python/src/gwr_parameteridentifikation/run_reference.py.
    use gwr_precision, only: dp
    use gwr_linear_algebra, only: adjoint_weighted, adjoint_euclidean
    use gwr_evaluation, only: analysis_result, analyze
    implicit none

    ! Direct cross-language equality is restricted to stable cases.  The
    ! regularized grids match the Python and MATLAB reference runners.
    integer, parameter :: direct_sizes(3) = (/ 1, 2, 4 /)
    integer, parameter :: beta_1e3_sizes(6) = (/ 1, 2, 4, 8, 16, 32 /)
    integer, parameter :: beta_1e6_sizes(7) = (/ 1, 2, 4, 8, 16, 32, 128 /)
    character(len=64), parameter :: csv_path = 'reference_fortran.csv'
    integer :: csv_unit, io_status

    csv_unit = 20
    open(unit=csv_unit, file=trim(csv_path), status='replace', &
         action='write', iostat=io_status)
    if (io_status /= 0) then
        write(*, '(A,1X,A)') 'Cannot open CSV output:', trim(csv_path)
        stop 1
    end if
    write(csv_unit, '(A)') &
        'implementation,method,n,beta,condition,error,relative_error'

    call run_direct_cases(csv_unit)
    call run_regularized_cases(csv_unit, beta_1e3_sizes, 1.0e-3_dp, &
                               adjoint_weighted, 'weighted')
    call run_regularized_cases(csv_unit, beta_1e3_sizes, 1.0e-3_dp, &
                               adjoint_euclidean, 'euclidean')
    call run_regularized_cases(csv_unit, beta_1e6_sizes, 1.0e-6_dp, &
                               adjoint_weighted, 'weighted')
    call run_regularized_cases(csv_unit, beta_1e6_sizes, 1.0e-6_dp, &
                               adjoint_euclidean, 'euclidean')

    close(csv_unit)
    write(*, '(/,A,1X,A)') 'CSV written to', trim(csv_path)

contains

    subroutine run_direct_cases(unit_number)
        integer, intent(in) :: unit_number
        type(analysis_result) :: result
        integer :: i, info
        character(len=256) :: message

        write(*, '(A)') 'Direct solution without regularization'
        write(*, '(A5,3A16)') 'n', 'cond(A)', 'error', 'relative error'
        do i = 1, size(direct_sizes)
            call analyze(direct_sizes(i), 0.0_dp, adjoint_weighted, &
                         result, info, message)
            call require_success(info, message)
            write(*, '(I5,3ES16.6)') result%n, result%condition, &
                result%direct_error, result%direct_relative_error
            write(unit_number, '(A,I0,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3)') &
                'fortran,direct,', result%n, ',0.0,', result%condition, ',', &
                result%direct_error, ',', result%direct_relative_error
        end do
    end subroutine run_direct_cases

    subroutine run_regularized_cases(unit_number, sizes, beta, mode, mode_name)
        ! Equivalent to Python _print_regularized(), with an additional CSV
        ! row containing the unrounded values used by compare_with_python.py.
        integer, intent(in) :: unit_number, sizes(:), mode
        real(dp), intent(in) :: beta
        character(len=*), intent(in) :: mode_name
        type(analysis_result) :: result
        integer :: i, info
        character(len=256) :: message

        write(*, '(/,A,ES10.2,A,A)') 'Tikhonov: beta=', beta, &
            ', adjoint=', trim(mode_name)
        write(*, '(A5,3A16)') 'n', 'cond(A)', 'error', 'relative error'
        do i = 1, size(sizes)
            call analyze(sizes(i), beta, mode, result, info, message)
            call require_success(info, message)
            write(*, '(I5,3ES16.6)') result%n, result%condition, &
                result%regularized_error, result%regularized_relative_error
            write(unit_number, '(A,A,A,I0,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3,A,ES24.16E3)') &
                'fortran,', trim(mode_name), ',', result%n, ',', beta, ',', &
                result%condition, ',', result%regularized_error, ',', &
                result%regularized_relative_error
        end do
    end subroutine run_regularized_cases

    subroutine require_success(info, message)
        integer, intent(in) :: info
        character(len=*), intent(in) :: message

        if (info /= 0) then
            write(*, '(A,I0,A,A)') 'Numerical operation failed (info=', &
                info, '): ', trim(message)
            stop 1
        end if
    end subroutine require_success

end program run_reference

