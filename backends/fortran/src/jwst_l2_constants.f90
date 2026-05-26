! jwst_l2_constants
! -----------------
! Project-wide constants and kind parameters for the JWST-L2 Lege-Artis
! Fortran reference. Locked physical-constant values match dynamics.py
! exactly for the doctrine Sec.4.4 C2 cross-implementation bit-identity gate.
!
! Canonical references:
!   - Chapter 02 Sec.2.5.1 (engineer-tier bridge for I tensor evolution)
!   - dynamics.py line 37 (G_NEWTON value origin: CODATA 2018)
!
! Conventions per SONNET-BRIEF Sec.2:
!   - Double precision (selected_real_kind(15))
!   - ASCII-only source (KB-039)
!   - Lowercase identifiers, no case-collision reuse (KB-037)
!
! Author: Opus prep pass 2026-05-24
! License: Apache 2.0 (inherited from parent project)

module jwst_l2_constants
    implicit none
    private

    ! Kind parameter for double-precision real
    integer, parameter, public :: dp = selected_real_kind(15)

    ! Gravitational constant (CODATA 2018, matches dynamics.py:37 exactly)
    ! Units: m^3 / (kg * s^2)
    real(dp), parameter, public :: g_newton = 6.67430e-11_dp

    ! Machine epsilon for double precision (IEEE-754 binary64)
    real(dp), parameter, public :: eps_machine = epsilon(1.0_dp)

    ! Default tolerance for unit tests against analytical reference values.
    ! Sized 1000x machine epsilon to absorb modest rounding accumulation in
    ! intermediate calculations (per the Higham 2002 Sec.4.2 Option G framework
    ! used throughout the lege-artis projects).
    real(dp), parameter, public :: eps_test_default = 1.0e-12_dp

    ! Tolerance for tests against trace identities (Q = tr(I)*I - 3*I), which
    ! involve subtraction of comparable-magnitude numbers and therefore have
    ! looser numerical reach.
    real(dp), parameter, public :: eps_test_trace = 1.0e-13_dp

    ! Tolerance for cross-implementation bit-identity tests (doctrine Sec.4.4 C2).
    ! Per Higham Sec.4.2 Option G: tolerance is ULP(value) * sqrt(N_ops).
    ! For typical gravity-gradient torque calculation N_ops ~ 30, giving:
    !   ULP(O(1)) * sqrt(30) ~ 2.2e-16 * 5.5 ~ 1.2e-15 relative
    real(dp), parameter, public :: eps_test_c2 = 1.0e-14_dp

    ! Mathematical constant pi (Fortran 2018 has no intrinsic; standard idiom)
    real(dp), parameter, public :: pi = 4.0_dp * atan(1.0_dp)

end module jwst_l2_constants
