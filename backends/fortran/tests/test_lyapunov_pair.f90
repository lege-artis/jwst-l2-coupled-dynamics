! test_lyapunov_pair  (D9 Case 3)
! ---------------------------------
! Numerical Lyapunov sensitivity via paired integration of an asymmetric rigid body.
!
! Setup (PHASE C REDESIGN from original m_b=100 kg at r=1e8 m which had negligible
! coupling; this setup uses Dzhanibekov flip timing sensitivity instead):
!   Body A: I = diag(1000, 2000, 5000) kg m^2 (asymmetric; same as D9 Case 2).
!   Body B: negligible (m_B = 1e-20 kg at x_B = [1e15, 0, 0]).
!   Initial spin: omega_A = (1e-4, 0.1, 1e-4) rad/s — near intermediate axis.
!
! Perturbation:
!   State X_0 and perturbed X_0 + delta_X_0, where delta_X_0 is a unit
!   perturbation in omega_A(2) (index 12 in the packed state vector),
!   scaled to |delta_X_0| = 1e-10.  This perturbs the near-intermediate-axis
!   spin, which is on the Dzhanibekov separatrix — a known high-sensitivity region.
!
! Physical basis:
!   For a torque-free asymmetric body with near-intermediate-axis spin, the flip
!   timing (half-period between sign changes of omega_y) is exponentially sensitive
!   to the initial omega_y value.  Even a 1e-10 perturbation in omega_y causes the
!   two trajectories to acquire a flip-timing offset that grows to O(T_flip/2) after
!   the first flip, producing delta_X >> dX0.  This is LARGE-SCALE SENSITIVITY on
!   the Euler-equation separatrix, not a positive Lyapunov exponent in the strict
!   sense (the free rotor is integrable), but the effective growth rate measured from
!   the linear regression slope is strongly positive over the first few flip periods.
!
! Acceptance:
!   T1: |delta_X| grew >= 2x from initial by the end of integration.
!   T2: lambda_L > 0 (positive slope in log(delta_X/dX0) vs t linear regression).
!   (No fixed numerical target; the deliverable is the extracted value.)
!
! Output: _audit/lyapunov_telemetry.csv  (t, log_delta_X_over_delta_X0, delta_X)
!
! Exit 0 on all PASS, exit 1 on any FAIL.
! Canonical reference: Chapter 02 Sec.2.5.1 (dI/dt sensitivity framework);
!   Landau-Lifshitz Vol I Sec.37 (intermediate-axis instability).
!
! Author: Sonnet Phase A 2026-05-24.  Phase C redesign: ThinkPad 2026-05-25.

program test_lyapunov_pair
    use jwst_l2_constants, only: dp
    use jwst_l2_dynamics,  only: body_state_t, pack_state, unpack_state, rk4_step
    implicit none

    integer :: n_pass, n_fail
    n_pass = 0
    n_fail = 0

    block
        real(dp), parameter :: m_a       = 1000.0_dp
        real(dp), parameter :: m_b       = 1.0e-20_dp
        real(dp), parameter :: dt_s      = 0.1_dp
        integer,  parameter :: nstep     = 10000      ! T = 1000 s (same as Dzhanibekov)
        integer,  parameter :: samp_ev   = 100
        integer,  parameter :: n_samp    = nstep / samp_ev + 1   ! 101
        real(dp), parameter :: dX0_mag   = 1.0e-10_dp

        real(dp) :: I_a(3,3), I_b(3,3)
        type(body_state_t) :: ba0, bb0
        real(dp) :: state0(26), state_p(26)

        real(dp) :: t_arr(n_samp), logdX_arr(n_samp), dX_arr(n_samp)
        integer  :: isamp, istep, i

        real(dp) :: dX0, dX_cur
        real(dp) :: sum_t, sum_y, sum_tt, sum_ty, n_fit_r, slope_lyap, lambda_L

        integer :: csv_unit, ios_csv

        ! ----- Inertia tensors: asymmetric Dzhanibekov body -----
        I_a = 0.0_dp
        I_a(1,1) = 1000.0_dp;  I_a(2,2) = 2000.0_dp;  I_a(3,3) = 5000.0_dp
        I_b = 0.0_dp
        I_b(1,1) = 1.0e-20_dp;  I_b(2,2) = 1.0e-20_dp;  I_b(3,3) = 1.0e-20_dp

        ! ----- Initial conditions: near intermediate-axis spin -----
        ba0%x     = [0.0_dp,    0.0_dp, 0.0_dp]
        ba0%v     = [0.0_dp,    0.0_dp, 0.0_dp]
        ba0%q     = [1.0_dp,    0.0_dp, 0.0_dp, 0.0_dp]
        ba0%omega = [1.0e-4_dp, 0.1_dp, 1.0e-4_dp]   ! near intermediate axis

        bb0%x     = [1.0e15_dp, 0.0_dp, 0.0_dp]
        bb0%v     = [0.0_dp,    0.0_dp, 0.0_dp]
        bb0%q     = [1.0_dp,    0.0_dp, 0.0_dp, 0.0_dp]
        bb0%omega = [0.0_dp,    0.0_dp, 0.0_dp]

        state0 = pack_state(ba0, bb0)

        ! ----- Perturbed state: add delta_X_0 to omega_A(2) (state index 12) -----
        ! Perturbing omega_y exploits Dzhanibekov flip-timing sensitivity on the
        ! separatrix; even 1e-10 perturbation in omega_y causes large divergence
        ! after the first flip (~200 s).
        state_p     = state0
        state_p(12) = state_p(12) + dX0_mag

        dX0 = dX0_mag   ! exact since only one component perturbed

        ! ----- Store first sample -----
        t_arr(1)     = 0.0_dp
        dX_cur       = state_norm_diff(state_p, state0)
        dX_arr(1)    = dX_cur
        logdX_arr(1) = log(dX_cur / dX0)
        isamp = 1

        ! ----- Integrate both states -----
        do istep = 1, nstep
            state0  = rk4_step(state0,  dt_s, m_a, I_a, m_b, I_b)
            state_p = rk4_step(state_p, dt_s, m_a, I_a, m_b, I_b)

            if (mod(istep, samp_ev) == 0) then
                isamp = isamp + 1
                dX_cur = state_norm_diff(state_p, state0)
                t_arr(isamp)     = real(istep, dp) * dt_s
                dX_arr(isamp)    = dX_cur
                logdX_arr(isamp) = log(max(dX_cur / dX0, 1.0e-300_dp))
            end if
        end do

        ! ----- Linear fit: log(dX/dX0) = lambda_L * t -----
        n_fit_r = real(isamp, dp)
        sum_t   = 0.0_dp;  sum_y  = 0.0_dp
        sum_tt  = 0.0_dp;  sum_ty = 0.0_dp
        do i = 1, isamp
            sum_t  = sum_t  + t_arr(i)
            sum_y  = sum_y  + logdX_arr(i)
            sum_tt = sum_tt + t_arr(i)**2
            sum_ty = sum_ty + t_arr(i) * logdX_arr(i)
        end do
        slope_lyap = (n_fit_r*sum_ty - sum_t*sum_y) / &
                     (n_fit_r*sum_tt - sum_t**2 + 1.0e-300_dp)
        lambda_L   = slope_lyap

        !-------------------------------------------------------------------
        ! T1: divergence grew by factor >= 2
        !-------------------------------------------------------------------
        call report("T1 delta_X grew >= 2x from initial", &
                    dX_arr(isamp) >= 2.0_dp * dX_arr(1), &
                    dX_arr(isamp) / dX_arr(1))

        !-------------------------------------------------------------------
        ! T2: positive effective growth rate from linear fit
        !-------------------------------------------------------------------
        call report("T2 lambda_L > 0 (positive effective growth rate)", &
                    lambda_L > 0.0_dp, lambda_L)

        write (*, '(a,es12.4,a)') "  lambda_L = ", lambda_L, " s^-1"
        if (lambda_L > 0.0_dp) &
            write (*, '(a,es12.4,a)') "  tau_eff  = ", 1.0_dp/lambda_L, " s"
        write (*, '(a,es12.4,a,es12.4)') &
            "  delta_X_0 = ", dX0, ";  delta_X_final = ", dX_arr(isamp)

        !-------------------------------------------------------------------
        ! Write CSV telemetry
        !-------------------------------------------------------------------
        open(newunit=csv_unit, file="_audit/lyapunov_telemetry.csv", &
             status="replace", iostat=ios_csv)
        if (ios_csv /= 0) then
            write (*, '(a)') "  WARNING: could not write _audit/lyapunov_telemetry.csv"
        else
            write (csv_unit, '(a)') "t_s,log_dX_over_dX0,delta_X"
            do i = 1, isamp
                write (csv_unit, '(es16.8,a,es16.8,a,es16.8)') &
                    t_arr(i), ",", logdX_arr(i), ",", dX_arr(i)
            end do
            close(csv_unit)
        end if
    end block

    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_lyapunov_pair: ", n_pass, " PASS / ", n_fail, " FAIL"
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

    !-----------------------------------------------------------------------
    ! state_norm_diff: Euclidean norm of (s2 - s1).
    !-----------------------------------------------------------------------
    pure real(dp) function state_norm_diff(s2, s1)
        real(dp), intent(in) :: s2(26), s1(26)
        real(dp) :: d(26)
        integer  :: k
        d = s2 - s1
        state_norm_diff = 0.0_dp
        do k = 1, 26
            state_norm_diff = state_norm_diff + d(k)**2
        end do
        state_norm_diff = sqrt(state_norm_diff)
    end function state_norm_diff

end program test_lyapunov_pair
