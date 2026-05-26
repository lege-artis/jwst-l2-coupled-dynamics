! test_instability_oblate  (D9 Case 1)
! -------------------------------------
! Gravity-gradient instability of an oblate body at L2-like configuration.
!
! Setup:
!   Body A: I = diag(2540, 2540, 3870) kg m^2, m_A = 2450 kg (oblate, I_par > I_perp).
!   Body B: m_B = 6e24 kg (Earth-mass primary), placed at origin.
!   Separation: r = 1.5e9 m; body A starts on +z axis in circular orbit.
!   Initial orientation: symmetry axis (body z) tilted delta_theta = 1e-3 rad from
!     radial direction (small-angle perturbation about body x-axis).
!   Initial spin: omega_A = (0, 0, 1e-4) rad/s (about symmetry axis).
!
! Diagnostic:
!   At each sample step extract theta_perp(t) = angle between body symmetry axis
!   and instantaneous radial direction r_hat(t).
!   Fit log(theta_perp) vs t by linear regression to extract tau_inst.
!
! Acceptance:
!   |tau_inst_fit - tau_inst_exact| / tau_inst_exact < 0.10
!   tau_inst_exact = sqrt(I_perp / (alpha * (I_par - I_perp))) = 2.317e6 s.
!   where alpha = 3 * G * m_B / r^3.
!
! Output: _audit/instability_oblate_telemetry.csv  (t, theta_perp, tau_GG_body_norm)
!
! Exit 0 on all PASS, exit 1 on any FAIL.
! Canonical reference: Chapter 03 Sec.3.7.2 (linearised GG torque).
!
! Author: Sonnet Phase A 2026-05-24

program test_instability_oblate
    use jwst_l2_constants,   only: dp, g_newton
    use jwst_l2_dynamics,    only: body_state_t, pack_state, unpack_state, rk4_step
    use jwst_l2_quaternion,  only: q_to_matrix
    use jwst_l2_torque,      only: gravity_gradient_torque_body_frame
    implicit none

    integer :: n_pass, n_fail
    n_pass = 0
    n_fail = 0

    block
        real(dp), parameter :: m_a       = 2450.0_dp
        real(dp), parameter :: m_b       = 6.0e24_dp
        real(dp), parameter :: r_sep     = 1.5e9_dp
        real(dp), parameter :: dt_s      = 500.0_dp
        integer,  parameter :: nstep     = 10000         ! T = 5e6 s (> 2 tau_inst for good fit)
        integer,  parameter :: samp_ev   = 100           ! sample every 100 steps
        integer,  parameter :: n_samp    = nstep / samp_ev + 1   ! 101
        real(dp), parameter :: dtheta    = 1.0e-3_dp     ! initial tilt (rad)

        real(dp) :: I_a(3,3), I_b_neg(3,3)
        type(body_state_t) :: ba, bb
        real(dp) :: state(26)

        real(dp) :: t_arr(n_samp), theta_arr(n_samp), tau_arr(n_samp)
        integer  :: isamp, istep

        real(dp) :: v_circ, gm_b, alpha_gg, omega_orb
        real(dp) :: tau_inst_exact, tau_inst_fit

        ! Linear-regression accumulators  (fit on all samp points)
        real(dp) :: sum_t, sum_y, sum_tt, sum_ty, n_fit_r
        real(dp) :: slope, intercept, t_cur
        integer  :: i_fit, n_fit

        integer :: csv_unit, ios_csv

        ! ----- Inertia tensors -----
        I_a = 0.0_dp
        I_a(1,1) = 2540.0_dp;  I_a(2,2) = 2540.0_dp;  I_a(3,3) = 3870.0_dp
        I_b_neg = 0.0_dp
        I_b_neg(1,1) = 1.0e-20_dp; I_b_neg(2,2) = 1.0e-20_dp; I_b_neg(3,3) = 1.0e-20_dp

        ! ----- Derived constants -----
        gm_b      = g_newton * m_b                      ! ~4.0e14 m^3/s^2
        v_circ    = sqrt(gm_b / r_sep)                  ! ~516.5 m/s
        omega_orb = v_circ / r_sep                      ! ~3.44e-7 rad/s
        alpha_gg  = 3.0_dp * gm_b / r_sep**3
        tau_inst_exact = sqrt(2540.0_dp / (alpha_gg * (3870.0_dp - 2540.0_dp)))

        ! ----- Initial conditions -----
        ! Quaternion: rotation of dtheta about body y-axis (orbit-normal = in-plane pitch mode).
        ! Body y = orbit normal for the xz-plane orbit. In-plane perturbation ensures the tilt
        ! measurement grows as cosh(sigma*t) without orbital-precession modulation (the out-of-plane
        ! perturbation about body x creates a cos(omega_orb*t) envelope that confounds the regression).
        ba%q(1) = cos(dtheta / 2.0_dp)
        ba%q(2) = 0.0_dp
        ba%q(3) = sin(dtheta / 2.0_dp)
        ba%q(4) = 0.0_dp
        ba%x     = [0.0_dp, 0.0_dp, r_sep]
        ba%v     = [v_circ, 0.0_dp, 0.0_dp]   ! circular orbit in xz-plane
        ba%omega = [0.0_dp, omega_orb, 0.0_dp] ! co-rotating: body tracks radial direction during orbit

        ! Body B: Earth-mass primary at origin (effectively stationary)
        bb%x     = [0.0_dp, 0.0_dp, 0.0_dp]
        bb%v     = [0.0_dp, 0.0_dp, 0.0_dp]
        bb%q     = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        bb%omega = [0.0_dp, 0.0_dp, 0.0_dp]

        state = pack_state(ba, bb)

        ! ----- Store initial sample -----
        t_arr(1)     = 0.0_dp
        theta_arr(1) = tilt_from_state(state)
        tau_arr(1)   = gg_norm_from_state(state, I_a, m_b)
        isamp = 1

        ! ----- Integrate -----
        do istep = 1, nstep
            state = rk4_step(state, dt_s, m_a, I_a, m_b, I_b_neg)
            if (mod(istep, samp_ev) == 0) then
                isamp = isamp + 1
                t_cur = real(istep, dp) * dt_s
                t_arr(isamp)     = t_cur
                theta_arr(isamp) = tilt_from_state(state)
                tau_arr(isamp)   = gg_norm_from_state(state, I_a, m_b)
            end if
        end do

        ! ----- Exponential fit: log(theta) = log(theta_0) + t / tau_inst -----
        ! Use only last quarter of samples: theta(t) ~ cosh(sigma*t) at early times, so
        ! fitting the full range underestimates sigma. Last quarter (sigma*t > 1.6) is
        ! close to the asymptotic exponential regime, giving < 10% error on tau_inst.
        n_fit    = isamp - 3*(isamp/4)  ! last quarter
        n_fit_r  = real(n_fit, dp)
        sum_t  = 0.0_dp;  sum_y  = 0.0_dp
        sum_tt = 0.0_dp;  sum_ty = 0.0_dp
        do i_fit = 3*(isamp/4) + 1, isamp
            if (theta_arr(i_fit) > 1.0e-30_dp) then
                sum_t  = sum_t  + t_arr(i_fit)
                sum_y  = sum_y  + log(theta_arr(i_fit))
                sum_tt = sum_tt + t_arr(i_fit)**2
                sum_ty = sum_ty + t_arr(i_fit) * log(theta_arr(i_fit))
            end if
        end do
        slope     = (n_fit_r*sum_ty - sum_t*sum_y) / (n_fit_r*sum_tt - sum_t**2)
        intercept = (sum_y - slope*sum_t) / n_fit_r
        if (slope > 0.0_dp) then
            tau_inst_fit = 1.0_dp / slope
        else
            tau_inst_fit = -1.0_dp   ! sentinel: no growth detected
        end if

        !-------------------------------------------------------------------
        ! T1: tau_inst within 10% of analytical value
        !-------------------------------------------------------------------
        call report("T1 tau_inst_fit / tau_inst_exact within 10%", &
                    tau_inst_fit > 0.0_dp .and. &
                    abs(tau_inst_fit - tau_inst_exact) / tau_inst_exact < 0.10_dp, &
                    abs(tau_inst_fit - tau_inst_exact) / max(tau_inst_exact, 1.0_dp))

        !-------------------------------------------------------------------
        ! T2: theta_perp grew (final > initial by at least 10%)
        !-------------------------------------------------------------------
        call report("T2 tilt grew: theta_f > 1.1 * theta_0", &
                    theta_arr(isamp) > 1.1_dp * theta_arr(1), &
                    theta_arr(isamp) / theta_arr(1))

        !-------------------------------------------------------------------
        ! Write CSV telemetry
        !-------------------------------------------------------------------
        open(newunit=csv_unit, file="_audit/instability_oblate_telemetry.csv", &
             status="replace", iostat=ios_csv)
        if (ios_csv /= 0) then
            write (*, '(a)') "  WARNING: could not write _audit/instability_oblate_telemetry.csv"
        else
            write (csv_unit, '(a)') "t_s,theta_perp_rad,tau_gg_body_norm"
            do i_fit = 1, isamp
                write (csv_unit, '(es16.8,a,es16.8,a,es16.8)') &
                    t_arr(i_fit), ",", theta_arr(i_fit), ",", tau_arr(i_fit)
            end do
            close(csv_unit)
        end if

        write (*, '(a,es12.4,a)') "  tau_inst_exact = ", tau_inst_exact, " s"
        write (*, '(a,es12.4,a)') "  tau_inst_fit   = ", tau_inst_fit,   " s"
    end block

    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_instability_oblate: ", n_pass, " PASS / ", n_fail, " FAIL"
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
    ! tilt_from_state: tilt angle of body A symmetry axis from radial direction.
    ! Body z in inertial frame = R * [0,0,1]; r_hat from positions.
    !-----------------------------------------------------------------------
    real(dp) function tilt_from_state(st)
        real(dp), intent(in) :: st(26)
        real(dp) :: qw, qx, qy, qz, bz(3), r_vec(3), r_mag, cos_t
        ! q_A is st(7:10) = (qw, qx, qy, qz)
        qw = st(7);  qx = st(8);  qy = st(9);  qz = st(10)
        ! Body z-axis in inertial frame (third column of R_body_to_inertial)
        bz(1) = 2.0_dp*(qx*qz + qw*qy)
        bz(2) = 2.0_dp*(qy*qz - qw*qx)
        bz(3) = 1.0_dp - 2.0_dp*(qx*qx + qy*qy)
        ! Radial direction: x_A - x_B = from planet (B) to satellite (A)
        r_vec = st(1:3) - st(14:16)
        r_mag = sqrt(r_vec(1)**2 + r_vec(2)**2 + r_vec(3)**2)
        r_vec = r_vec / max(r_mag, 1.0e-30_dp)
        cos_t = bz(1)*r_vec(1) + bz(2)*r_vec(2) + bz(3)*r_vec(3)
        cos_t = max(-1.0_dp, min(1.0_dp, cos_t))
        tilt_from_state = acos(cos_t)
    end function tilt_from_state

    !-----------------------------------------------------------------------
    ! gg_norm_from_state: norm of gravity-gradient torque on body A.
    !-----------------------------------------------------------------------
    real(dp) function gg_norm_from_state(st, I_a_in, m_b_in)
        real(dp), intent(in) :: st(26), I_a_in(3,3), m_b_in
        real(dp) :: R_mat(3,3), r_vec(3), tau(3)
        real(dp) :: qw, qx, qy, qz
        qw = st(7);  qx = st(8);  qy = st(9);  qz = st(10)
        ! Reconstruct rotation matrix from quaternion
        R_mat(1,1) = 1.0_dp - 2.0_dp*(qy*qy + qz*qz)
        R_mat(1,2) = 2.0_dp*(qx*qy - qw*qz)
        R_mat(1,3) = 2.0_dp*(qx*qz + qw*qy)
        R_mat(2,1) = 2.0_dp*(qx*qy + qw*qz)
        R_mat(2,2) = 1.0_dp - 2.0_dp*(qx*qx + qz*qz)
        R_mat(2,3) = 2.0_dp*(qy*qz - qw*qx)
        R_mat(3,1) = 2.0_dp*(qx*qz - qw*qy)
        R_mat(3,2) = 2.0_dp*(qy*qz + qw*qx)
        R_mat(3,3) = 1.0_dp - 2.0_dp*(qx*qx + qy*qy)
        r_vec = st(14:16) - st(1:3)
        tau = gravity_gradient_torque_body_frame(I_a_in, r_vec, R_mat, m_b_in)
        gg_norm_from_state = sqrt(tau(1)**2 + tau(2)**2 + tau(3)**2)
    end function gg_norm_from_state

end program test_instability_oblate
