program run_benchmark
    ! Single-thread benchmark matching Python benchmark.py and MATLAB
    ! run_benchmark.m.  Thread count and CPU affinity are controlled by the
    ! calling shell, not by the numerical source code.
    use gwr_precision, only: dp
    use gwr_performance, only: benchmark_result, measure_performance
    implicit none

    integer, parameter :: sizes(5) = (/ 64, 128, 256, 512, 1024 /)
    integer, parameter :: repeats = 9
    real(dp), parameter :: beta = 1.0e-6_dp
    character(len=*), parameter :: output_path = 'benchmark_fortran.csv'
    character(len=*), parameter :: thread_mode = 'single-cpu0'
    integer, parameter :: thread_count = 1
    type(benchmark_result) :: result
    integer :: i, info, output_unit, io_status
    character(len=256) :: message

    output_unit = 30
    open(unit=output_unit, file=output_path, status='replace', &
         action='write', iostat=io_status)
    if (io_status /= 0) stop 1
    call write_header(output_unit)

    write(*, '(A)') 'Fortran dense benchmark (median after warm-up)'
    write(*, '(A6,4A16)') 'n', 'assembly', 'direct solve', &
        'Tikh. setup', 'Tikh. solve'
    do i = 1, size(sizes)
        call measure_performance( &
            sizes(i), beta, repeats, result, info, message)
        call require_success(info, message)
        call write_result(output_unit, result, thread_mode, thread_count)
        call print_result(result)
    end do
    close(output_unit)
    write(*, '(A,1X,A)') 'Benchmark written to', output_path

contains

    subroutine write_header(unit_number)
        integer, intent(in) :: unit_number
        write(unit_number, '(A)') 'implementation,thread_mode,'// &
            'thread_count,n,beta,assembly_seconds,direct_solve_seconds,'// &
            'tikhonov_setup_seconds,tikhonov_solve_seconds'
    end subroutine write_header

    subroutine write_result(unit_number, measured, mode, threads)
        integer, intent(in) :: unit_number, threads
        type(benchmark_result), intent(in) :: measured
        character(len=*), intent(in) :: mode

        write(unit_number, '(A,A,A,I0,A,I0,5(A,ES24.16E3))') &
            'Fortran/gfortran,', trim(mode), ',', threads, ',', measured%n, &
            ',', measured%beta, ',', measured%assembly_seconds, ',', &
            measured%direct_solve_seconds, ',', &
            measured%tikhonov_setup_seconds, ',', &
            measured%tikhonov_solve_seconds
    end subroutine write_result

    subroutine print_result(measured)
        type(benchmark_result), intent(in) :: measured
        write(*, '(I6,4ES16.6)') measured%n, measured%assembly_seconds, &
            measured%direct_solve_seconds, measured%tikhonov_setup_seconds, &
            measured%tikhonov_solve_seconds
    end subroutine print_result

    subroutine require_success(status, text)
        integer, intent(in) :: status
        character(len=*), intent(in) :: text
        if (status /= 0) then
            write(*, '(A,I0,A,A)') 'Benchmark failed (info=', status, &
                '): ', trim(text)
            stop 1
        end if
    end subroutine require_success

end program run_benchmark

