C     =================================================================
C     File: toyprob.f
C     =================================================================

C     =================================================================
C     Module: Spectral Projected Gradient Method. Problem definition.
C     =================================================================

C     Last update of any of the component of this module:

C     March 14, 2008.

C     Users are encouraged to download periodically updated versions of
C     this code at the TANGO Project web page:
C
C     www.ime.usp.br/~egbirgin/tango/

C     *****************************************************************
C     *****************************************************************

      subroutine inip(n, x)

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n

        ! ARRAY ARGUMENTS
        double precision, intent(inout) :: x(n)

        ! PARAMETERS
        integer, parameter :: nmax = 100000

        ! COMMON ARRAYS
        double precision l(nmax), u(nmax)

        ! LOCAL SCALARS
        integer :: i

        ! COMMON BLOCKS
        common /bounds/ l, u

        ! Number of variables
        n = 10

        ! Initial point
        do i = 1, n
            x(i) = 60.0d0
        end do

        ! Bound constraints (or any other constraints that define a
        ! convex set)
        do i = 1, n
            l(i) = -100.0d0
            u(i) = 50.0d0
        end do

      end subroutine inip

      subroutine evalf(n, x, f, flag)

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: f
        integer, intent(out) :: flag

        ! LOCAL SCALARS
        integer :: i

        flag = 0

        f = 0.0d0
        do i = 1, n
            f = f + x(i) ** 2
        end do

      end subroutine evalf

      subroutine evalg(n, x, g, flag)

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        double precision, intent(in) :: x(n)
        double precision, intent(out) :: g(n)
        integer, intent(out) :: flag

        ! LOCAL SCALARS
        integer :: i

        flag = 0

        do i = 1, n
            g(i) = 2.0d0 * x(i)
        end do

      end subroutine evalg

      subroutine proj(n, x, flag)

        ! SCALAR ARGUMENTS
        integer, intent(in) :: n
        double precision, intent(inout) :: x(n)
        integer, intent(out) :: flag

        ! PARAMETERS
        integer, parameter :: nmax = 100000

        ! COMMON ARRAYS
        double precision l(nmax), u(nmax)

        ! LOCAL SCALARS
        integer :: i

        ! COMMON BLOCKS
        common /bounds/ l, u

        flag = 0

        do i = 1, n
            x(i) = max(l(i), min(x(i), u(i)))
        end do

      end subroutine proj

      program spgma

        use adolc
        implicit none

        integer :: n
        double precision :: x(100000), f
        integer :: flag

        call inip(n, x)
        call evalf(n, x, f, flag)
        write(*,*) "f =", f

      end program spgma
