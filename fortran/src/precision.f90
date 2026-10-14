module gwr_precision
    implicit none

    ! Match the double precision used by numpy.float64 in the Python
    ! reference (python/src/gwr_parameteridentifikation/model.py).
    ! The external D... BLAS/LAPACK routines use Fortran DOUBLE PRECISION.
    ! Deriving dp from a double literal guarantees the same kind at their
    ! explicitly declared interfaces.
    integer, parameter :: dp = kind(1.0d0)

end module gwr_precision

