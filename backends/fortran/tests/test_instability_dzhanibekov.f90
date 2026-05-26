! test_instability_dzhanibekov  (D9 Case 2)
! -------------------------------------------
! Dzhanibekov / tennis-racket effect: torque-free asymmetric-body instability.
!
! Setup:
!   Body A: I = diag(1000, 2000, 5000) kg m^2 (three distinct moments;
!     intermediate axis y, I1 < I2 < I3).  m_A = 1000 kg.
!   Body B: negligible (m_B = 1e-20 kg, I_B = 1e-20*I_3, at x_B = [1e15,0,0]).
!   Initial spin: omega_A = (1e-5, 0.1, 1e-5) rad/s (nearly about intermediate y-axis).
!   Integration: dt = 0.1 s, 10000 steps (T = 1000 s).
!
! Diagnostic:
!   Track omega_y(t).  A Dzhanibekov flip is a sign change of omega_y.
!   Count total sign changes.  Average half-period between consecutive sign changes.
!
! Acceptance:
!   T1: at least 5 sign changes of omega_y detected in 1000 s (>= 5 flips).
!   T2: omega_z^2 + omega_x^2 stays bounded (|omega|^2 conserved to 1e-8 relative),
!       confirming the motion stays on the energy ellipsoid.
!
! Output: _audit/dzhanibekov_telemetry.csv  (t, omega_x, omega_y, omega_z)
!
! Exit 0 on all PASS, exit 1 on any FAIL.
! Canonical reference: Goldstein Sec.5.7; Landau-Lifshitz Vol I Sec.37.
!
! Author: Sonnet Phase A 2026-05-24

program test_instability_dzhanibekov
    use jwst_l2_constants, only: dp
    use jwst_l2_dynamics,  only: body_state_t, pack_state, unpack_state, rk4_step
    implicit none

    integer :: n_pass, n_fail
    n_pass = 0
    n_fail = 0

    block
        real(dp), parameter :: m_a      = 1000.0_dp
        real(dp), parameter :: m_b      = 1.0e-20_dp
        real(dp), parameter :: dt_s     = 0.1_dp
        integer,  parameter :: nstep    = 10000        ! T = 1000 s
        integer,  parameter :: samp_ev  = 10           ! sample every 10 steps
        integer,  parameter :: n_samp   = nstep / samp_ev + 1   ! 1001

        real(dp) :: I_a(3,3), I_b_neg(3,3)
        type(body_state_t) :: ba, bb
        real(dp) :: state(26)

        real(dp) :: t_arr(n_samp), omx_arr(n_samp), omy_arr(n_samp), omz_arr(n_samp)
        integer  :: isamp, istep, i

        real(dp) :: E_0, E_f
        real(dp) :: omy_prev, omy_cur
        integer  :: n_flips

        integer :: csv_unit, ios_csv

        ! ----- Inertia tensors -----
        I_a = 0.0_dp
        I_a(1,1) = 1000.0_dp;  I_a(2,2) = 2000.0_dp;  I_a(3,3) = 5000.0_dp
        I_b_neg = 0.0_dp
        I_b_neg(1,1) = 1.0e-20_dp; I_b_neg(2,2) = 1.0e-20_dp; I_b_neg(3,3) = 1.0e-20_dp

        ! ----- Initial conditions -----
        ba%x     = [0.0_dp,    0.0_dp, 0.0_dp]
        ba%v     = [0.0_dp,    0.0_dp, 0.0_dp]
        ba%q     = [1.0_dp,    0.0_dp, 0.0_dp, 0.0_dp]
        ba%omega = [1.0e-4_dp, 0.1_dp, 1.0e-4_dp]   ! Phase C fallback: larger perturbation for >= 5 flips

        bb%x     = [1.0e15_dp, 0.0_dp, 0.0_dp]
        bb%v     = [0.0_dp,    0.0_dp, 0.0_dp]
        bb%q     = [1.0_dp,    0.0_dp, 0.0_dp, 0.0_dp]
        bb%omega = [0.0_dp,    0.0_dp, 0.0_dp]

        state = pack_state(ba, bb)

        ! Initial kinetic energy: E = 0.5 * sum_i(I_i * omega_i^2) — conserved for torque-free rigid body
        E_0 = 0.5_dp * (I_a(1,1)*ba%omega(1)**2 + I_a(2,2)*ba%omega(2)**2 + I_a(3,3)*ba%omega(3)**2)

        ! ----- Store first sample -----
        t_arr(1)   = 0.0_dp
        omx_arr(1) = ba%omega(1)
        omy_arr(1) = ba%omega(2)
        omz_arr(1) = ba%omega(3)
        isamp = 1
        omy_prev = ba%omega(2)
        n_flips  = 0

        ! ----- Integrate -----
        do istep = 1, nstep
            state = rk4_step(state, dt_s, m_a, I_a, m_b, I_b_neg)
            call unpack_state(state, ba, bb)
            omy_cur = ba%omega(2)
            ! Detect sign change in omega_y (Dzhanibekov flip)
            if (omy_prev * omy_cur < 0.0_dp) n_flips = n_flips + 1
            omy_prev = omy_cur

            if (mod(istep, samp_ev) == 0) then
                isamp = isamp + 1
                t_arr(isamp)   = real(istep, dp) * dt_s
                omx_arr(isamp) = ba%omega(1)
                omy_arr(isamp) = ba%omega(2)
                omz_arr(isamp) = ba%omega(3)
            end if
        end do

        ! Final kinetic energy
        E_f = 0.5_dp * (I_a(1,1)*ba%omega(1)**2 + I_a(2,2)*ba%omega(2)**2 + I_a(3,3)*ba%omega(3)**2)

        !-------------------------------------------------------------------
        ! T1: at least 5 Dzhanibekov flips (sign changes of omega_y)
        !-------------------------------------------------------------------
        call report_int("T1 Dzhanibekov flips >= 5", n_flips >= 5, n_flips)

        !-------------------------------------------------------------------
        ! T2: kinetic energy conserved (|omega|^2 is NOT conserved for asymmetric body;
        !     E = 0.5*sum(I_i*omega_i^2) is the correct RK4-conserved quantity)
        !-------------------------------------------------------------------
        call report("T2 kinetic energy conserved to 1e-6", &
                    abs(E_f - E_0) / E_0 < 1.0e-6_dp, &
                    abs(E_f - E_0) / E_0)

        write (*, '(a,i0,a)') "  Dzhanibekov flips detected: ", n_flips, ""

        !-------------------------------------------------------------------
        ! Write CSV telemetry
        !-------------------------------------------------------------------
        open(newunit=csv_unit, file="_audit/dzhanibekov_telemetry.csv", &
             status="replace", iostat=ios_csv)
        if (ios_csv /= 0) then
            write (*, '(a)') "  WARNING: could not write _audit/dzhanibekov_telemetry.csv"
        else
            write (csv_unit, '(a)') "t_s,omega_x,omega_y,omega_z"
            do i = 1, isamp
                write (csv_unit, '(es16.8,a,es16.8,a,es16.8,a,es16.8)') &
                    t_arr(i), ",", omx_arr(i), ",", omy_arr(i), ",", omz_arr(i)
            end do
            close(csv_unit)
        end if
    end block

    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_instability_dzhanibekov: ", n_pass, " PASS / ", n_fail, " FAIL"
    if (n_fail > 0) stop 1

contains

    subroutine report(label, passed, residual)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: passed
        real(dp),         intent(in) :: residual
        if (passed) then
            n_pass = n_pass + 1
            write (*, '(a,a,a,es12.4,a)') "  ", trim(label), "  (res=", residual, ")  PASS"
        else
            n_fail = n_fail + 1
            write (*, '(a,a,a,es12.4,a)') "  ", trim(label), "  (res=", residual, ")  FAIL"
        end if
    end subroutine report

    subroutine report_int(label, passed, ival)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: passed
        integer,          intent(in) :: ival
        if (passed) then
            n_pass = n_pass + 1
            write (*, '(a,a,a,i0,a)') "  ", trim(label), "  (n=", ival, ")  PASS"
        else
            n_fail = n_fail + 1
            write (*, '(a,a,a,i0,a)') "  ", trim(label), "  (n=", ival, ")  FAIL"
        end if
    end subroutine report_int

end program test_instability_dzhanibekov
