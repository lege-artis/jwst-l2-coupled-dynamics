! test_trace_identity  (D2)
! --------------------------
! Tests the trace-free quadrupole identity Q = tr(I) * I_3 - 3 * I.
!
! Tests:
!   T1. Thin-rod sanity case: I = diag(1/3, 1/3, 0).
!       Expect Q = diag(-1/3, -1/3, 2/3) with tr(Q) = 0.
!   T2..T101. 100 random diagonal SPD matrices (eigenvalues from a
!       deterministic LCG in [0.1, 10]): verify |tr(Q)| < eps_test_trace.
!
! Acceptance: all 101 cases PASS.
! Exit 0 on all PASS, exit 1 on any FAIL.
!
! Canonical reference: Chapter 03 Sec.3.5 (trace identity).
!
! Author: Sonnet Phase A 2026-05-24

program test_trace_identity
    use jwst_l2_constants, only: dp, eps_test_trace, eps_test_default
    use jwst_l2_potential,  only: quadrupole_from_inertia
    implicit none

    integer  :: n_pass, n_fail
    real(dp) :: I_body(3, 3), Q_body(3, 3)
    real(dp) :: tr_Q, maxdiff, eig1, eig2, eig3
    integer  :: icase
    integer(kind=8) :: seed

    n_pass = 0
    n_fail = 0

    ! --- T1: thin-rod case ---
    ! I = diag(1/3, 1/3, 0), representing a unit-mass unit-length thin rod along z.
    ! Analytical Q: tr(I) = 2/3, Q_ii = tr(I) - 3*I_ii.
    !   Q_11 = 2/3 - 3*(1/3) = 2/3 - 1 = -1/3
    !   Q_22 = 2/3 - 3*(1/3) = -1/3
    !   Q_33 = 2/3 - 3*0 = 2/3
    I_body = 0.0_dp
    I_body(1, 1) = 1.0_dp / 3.0_dp
    I_body(2, 2) = 1.0_dp / 3.0_dp
    I_body(3, 3) = 0.0_dp

    Q_body = quadrupole_from_inertia(I_body)
    tr_Q = Q_body(1,1) + Q_body(2,2) + Q_body(3,3)
    call report("T1 thin-rod |tr(Q)| < eps_trace", abs(tr_Q) < eps_test_trace, abs(tr_Q))

    ! Check diagonal elements against analytical values
    maxdiff = max(abs(Q_body(1,1) - (-1.0_dp/3.0_dp)), &
                  abs(Q_body(2,2) - (-1.0_dp/3.0_dp)), &
                  abs(Q_body(3,3) - ( 2.0_dp/3.0_dp)))
    call report("T1 thin-rod Q diagonal elements vs analytical", &
                maxdiff < eps_test_default, maxdiff)

    ! --- T2..T101: 100 random diagonal SPD matrices ---
    ! Use a simple linear congruential generator (seed fixed for reproducibility).
    ! LCG: x_{n+1} = (1664525 * x_n + 1013904223) mod 2^32
    seed = 42_8
    do icase = 1, 100
        ! Draw 3 eigenvalues in [0.1, 10]
        eig1 = 0.1_dp + 9.9_dp * lcg_next(seed)
        eig2 = 0.1_dp + 9.9_dp * lcg_next(seed)
        eig3 = 0.1_dp + 9.9_dp * lcg_next(seed)

        I_body = 0.0_dp
        I_body(1, 1) = eig1
        I_body(2, 2) = eig2
        I_body(3, 3) = eig3

        Q_body = quadrupole_from_inertia(I_body)
        tr_Q = Q_body(1,1) + Q_body(2,2) + Q_body(3,3)

        if (abs(tr_Q) >= eps_test_trace) then
            ! Only print failures to keep output clean
            write (*, '(a,i3,a,es12.4)') "  T", icase+1, " FAIL |tr(Q)|=", abs(tr_Q)
            n_fail = n_fail + 1
        else
            n_pass = n_pass + 1
        end if
    end do

    ! Summary line for the 100 random cases
    if (n_fail == 0) then
        write (*, '(a,i3,a)') "  T2..T101 all 100 random SPD cases: ", n_pass - 2, " PASS"
    end if

    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_trace_identity: ", n_pass, " PASS / ", n_fail, " FAIL"
    if (n_fail > 0) stop 1

contains

    !-----------------------------------------------------------------------
    ! report: print test result and accumulate counters.
    !-----------------------------------------------------------------------
    subroutine report(label, passed, residual)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: passed
        real(dp),         intent(in) :: residual
        if (passed) then
            n_pass = n_pass + 1
            write (*, '(a,a,a,es12.4,a)') "  ", trim(label), "  (residual=", residual, ")  PASS"
        else
            n_fail = n_fail + 1
            write (*, '(a,a,a,es12.4,a)') "  ", trim(label), "  (residual=", residual, ")  FAIL"
        end if
    end subroutine report

    !-----------------------------------------------------------------------
    ! lcg_next: linear congruential generator, returns uniform in [0, 1).
    ! Numerical Recipes 3rd ed Sec.7.1 constants.
    !-----------------------------------------------------------------------
    real(dp) function lcg_next(iseed)
        integer(kind=8), intent(inout) :: iseed
        iseed = mod(1664525_8 * iseed + 1013904223_8, 4294967296_8)
        lcg_next = real(iseed, dp) / 4294967296.0_dp
    end function lcg_next

end program test_trace_identity
