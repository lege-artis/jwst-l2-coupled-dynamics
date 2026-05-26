! test_symmetric_top  (D8 -- SECOND C2 GATE)
! -------------------------------------------
! Torque-free Euler precession of an axisymmetric rigid body.
!
! Setup (matches tests/fixtures/symmetric_top_trajectory.json):
!   I_A = diag(2540, 2540, 3870) kg m^2, m_A = 2450 kg  (JWST-like oblate body)
!   I_B = 1e-20 * I_3, m_B = 1e-20 kg                   (negligible secondary)
!   omega_initial = (1e-3, 0, 1e-2) rad/s, q_initial = (1,0,0,0)
!   x_B_initial = (1e15, 0, 0) m  (far enough to suppress gravity-gradient torque)
!   dt = 0.05 s,  12000 RK4 steps  (window 600 s),  sample every 60 steps -> 201 points.
!
! Tests:
!   T1: Euler precession frequency.
!       omega_x(t) = omega_perp * cos(lambda * t).  At t=600 s the phase is
!       lambda*600 ~ pi.  atan2(omega_y_f, omega_x_f) gives ~-pi; unwrap by
!       +2*pi to recover total phase.  lambda_est = theta_unwrapped / 600.
!       |lambda_est - lambda_exact| / lambda_exact < 1e-6.
!       lambda_exact = (I_par - I_perp) / I_perp * omega_z = 5.236220...e-3 rad/s.
!   T2: Total mechanical energy conserved to 1e-11 relative.
!   T3: Total angular momentum magnitude conserved to 1e-9 relative.
!   T4 (C2 gate): q and omega_body of body A at steps 0, 60, 600, 6000, 12000
!       match the Python fixture values at eps_test_c2 * sqrt(step + 1) tolerance.
!
! Output: _audit/fortran_symmetric_top_output.json  (all 201 sampled states).
!
! Exit 0 on all PASS, exit 1 on any FAIL.
! Canonical reference: Landau-Lifshitz Vol I Sec.36 (torque-free symmetric top).
! Python oracle: tests/fixtures/symmetric_top_trajectory.json.
!
! Author: Sonnet Phase A 2026-05-24

program test_symmetric_top
    use jwst_l2_constants,   only: dp, eps_test_c2
    use jwst_l2_dynamics,    only: body_state_t, pack_state, unpack_state, rk4_step
    use jwst_l2_diagnostics, only: total_kinetic_energy, total_potential_energy, &
                                   total_angular_momentum
    implicit none

    integer :: n_pass, n_fail

    ! Compile-time constants
    integer,  parameter :: nstep_k    = 12000
    integer,  parameter :: sample_ev  = 60
    integer,  parameter :: n_samp     = nstep_k / sample_ev + 1   ! 201
    real(dp), parameter :: dt_s       = 0.05_dp
    real(dp), parameter :: t_tot_s    = real(nstep_k, dp) * dt_s  ! 600 s
    real(dp), parameter :: pi_c       = 4.0_dp * atan(1.0_dp)

    ! Working variables
    real(dp) :: m_a, m_b
    real(dp) :: I_a(3,3), I_b(3,3)
    type(body_state_t) :: ba, bb
    real(dp) :: state(26)
    real(dp) :: traj_q(4, n_samp), traj_om(3, n_samp)
    integer  :: isamp, istep
    real(dp) :: E_ke, E_tr, E_rot, V_pe, E_0, E_f
    real(dp) :: L_vec_0(3), L_vec_f(3), L0_mag, Lf_mag
    real(dp) :: lambda_exact, lambda_est, theta_f
    real(dp) :: q_ref(4), om_ref(3), resid_c2, tol_c2
    integer  :: chk_step, chk_idx

    n_pass = 0
    n_fail = 0

    ! ----- Inertia tensors -----
    I_a = 0.0_dp
    I_a(1,1) = 2540.0_dp;  I_a(2,2) = 2540.0_dp;  I_a(3,3) = 3870.0_dp
    I_b = 0.0_dp
    I_b(1,1) = 1.0e-20_dp; I_b(2,2) = 1.0e-20_dp; I_b(3,3) = 1.0e-20_dp

    m_a = 2450.0_dp
    m_b = 1.0e-20_dp

    ! ----- Initial conditions (from fixture setup) -----
    ba%x     = [0.0_dp,    0.0_dp, 0.0_dp]
    ba%v     = [0.0_dp,    0.0_dp, 0.0_dp]
    ba%q     = [1.0_dp,    0.0_dp, 0.0_dp, 0.0_dp]
    ba%omega = [1.0e-3_dp, 0.0_dp, 1.0e-2_dp]

    bb%x     = [1.0e15_dp, 0.0_dp, 0.0_dp]
    bb%v     = [0.0_dp,    0.0_dp, 0.0_dp]
    bb%q     = [1.0_dp,    0.0_dp, 0.0_dp, 0.0_dp]
    bb%omega = [0.0_dp,    0.0_dp, 0.0_dp]

    state = pack_state(ba, bb)

    ! ----- Initial diagnostics -----
    call total_kinetic_energy(state, m_a, I_a, m_b, I_b, E_tr, E_rot, E_ke)
    V_pe    = total_potential_energy(state, m_a, m_b)
    E_0     = E_ke + V_pe
    L_vec_0 = total_angular_momentum(state, m_a, I_a, m_b, I_b)
    L0_mag  = sqrt(L_vec_0(1)**2 + L_vec_0(2)**2 + L_vec_0(3)**2)

    ! ----- Store first sample -----
    traj_q(:, 1)  = ba%q
    traj_om(:, 1) = ba%omega
    isamp = 1

    ! ----- Integrate 12000 RK4 steps; sample every 60 -----
    do istep = 1, nstep_k
        state = rk4_step(state, dt_s, m_a, I_a, m_b, I_b)
        if (mod(istep, sample_ev) == 0) then
            isamp = isamp + 1
            call unpack_state(state, ba, bb)
            traj_q(:, isamp)  = ba%q
            traj_om(:, isamp) = ba%omega
        end if
    end do

    ! Final state (ba, bb already set by last unpack_state inside loop)
    call total_kinetic_energy(state, m_a, I_a, m_b, I_b, E_tr, E_rot, E_ke)
    V_pe    = total_potential_energy(state, m_a, m_b)
    E_f     = E_ke + V_pe
    L_vec_f = total_angular_momentum(state, m_a, I_a, m_b, I_b)
    Lf_mag  = sqrt(L_vec_f(1)**2 + L_vec_f(2)**2 + L_vec_f(3)**2)

    !-----------------------------------------------------------------------
    ! T1: Euler precession frequency
    ! omega_x(t) = omega_perp*cos(lambda*t), omega_y(t) = omega_perp*sin(lambda*t).
    ! At t=600 s, lambda*600 ~ pi.  atan2 returns ~-pi (just past the branch cut);
    ! add 2*pi to recover the accumulated phase > 0, then divide by t_tot_s.
    !-----------------------------------------------------------------------
    lambda_exact = (3870.0_dp - 2540.0_dp) / 2540.0_dp * 1.0e-2_dp
    theta_f = atan2(ba%omega(2), ba%omega(1))
    if (theta_f < 0.0_dp) theta_f = theta_f + 2.0_dp * pi_c
    lambda_est = theta_f / t_tot_s
    call report("T1 precession |dlambda/lambda| < 1e-6", &
                abs(lambda_est - lambda_exact) / lambda_exact < 1.0e-6_dp, &
                abs(lambda_est - lambda_exact) / lambda_exact)

    !-----------------------------------------------------------------------
    ! T2: Total mechanical energy conserved to 1e-11 relative
    !-----------------------------------------------------------------------
    call report("T2 energy |dE/E0| < 1e-11", &
                abs(E_f - E_0) / abs(E_0) < 1.0e-11_dp, &
                abs(E_f - E_0) / abs(E_0))

    !-----------------------------------------------------------------------
    ! T3: Total angular momentum magnitude conserved to 1e-9 relative
    !-----------------------------------------------------------------------
    call report("T3 angular momentum |dL/L0| < 1e-9", &
                abs(Lf_mag - L0_mag) / L0_mag < 1.0e-9_dp, &
                abs(Lf_mag - L0_mag) / L0_mag)

    !-----------------------------------------------------------------------
    ! T4 (C2 gate): element-wise comparison to Python fixture at 5 checkpoints.
    ! Tolerance per checkpoint: eps_test_c2 * sqrt(step + 1).
    ! Fixture source: tests/fixtures/symmetric_top_trajectory.json (hardcoded).
    !-----------------------------------------------------------------------

    ! --- step 0 (isamp index 1) ---
    chk_step = 0;  chk_idx = 1
    q_ref  = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    om_ref = [1.0e-3_dp, 0.0_dp, 1.0e-2_dp]
    tol_c2 = eps_test_c2 * sqrt(real(chk_step + 1, dp))
    resid_c2 = max(maxval(abs(traj_q(:, chk_idx) - q_ref)), &
                   maxval(abs(traj_om(:, chk_idx) - om_ref)))
    call report("T4a C2 step=0 max|delta| < eps*sqrt(1)", &
                resid_c2 < tol_c2, resid_c2)

    ! --- step 60 (isamp index 2) ---
    chk_step = 60;  chk_idx = 2
    q_ref  = [9.998863771748902e-1_dp,  1.499822597231296e-3_dp, &
              1.1780344929140052e-5_dp,  1.4999437771927375e-2_dp]
    om_ref = [9.998766215153628e-4_dp,  1.570801537506176e-5_dp,  1.0e-2_dp]
    tol_c2 = eps_test_c2 * sqrt(real(chk_step + 1, dp))
    resid_c2 = max(maxval(abs(traj_q(:, chk_idx) - q_ref)), &
                   maxval(abs(traj_om(:, chk_idx) - om_ref)))
    call report("T4b C2 step=60 max|delta| < eps*sqrt(61)", &
                resid_c2 < tol_c2, resid_c2)

    ! --- step 600 (isamp index 11) ---
    chk_step = 600;  chk_idx = 11
    q_ref  = [9.886592317033457e-1_dp,  1.4823360073739807e-2_dp, &
              1.1666758037951975e-3_dp,  1.4943838339410465e-1_dp]
    om_ref = [9.876872484248281e-4_dp,  1.5644136057638386e-4_dp,  1.0e-2_dp]
    tol_c2 = eps_test_c2 * sqrt(real(chk_step + 1, dp))
    resid_c2 = max(maxval(abs(traj_q(:, chk_idx) - q_ref)), &
                   maxval(abs(traj_om(:, chk_idx) - om_ref)))
    call report("T4c C2 step=600 max|delta| < eps*sqrt(601)", &
                resid_c2 < tol_c2, resid_c2)

    ! --- step 6000 (isamp index 101) ---
    chk_step = 6000;  chk_idx = 101
    q_ref  = [6.468969106237606e-2_dp,   3.4828491135168775e-2_dp, &
              3.483092276898071e-2_dp,    9.966890321933733e-1_dp]
    om_ref = [-6.98149373310077e-8_dp,   9.999999975629443e-4_dp,   1.0e-2_dp]
    tol_c2 = eps_test_c2 * sqrt(real(chk_step + 1, dp))
    resid_c2 = max(maxval(abs(traj_q(:, chk_idx) - q_ref)), &
                   maxval(abs(traj_om(:, chk_idx) - om_ref)))
    call report("T4d C2 step=6000 max|delta| < eps*sqrt(6001)", &
                resid_c2 < tol_c2, resid_c2)

    ! --- step 12000 (isamp index 201) ---
    chk_step = 12000;  chk_idx = 201
    q_ref  = [-9.892041013786311e-1_dp,  4.532742399079596e-6_dp, &
              -6.492510858840163e-2_dp,   1.313772281254755e-1_dp]
    om_ref = [-9.99999990251751e-4_dp,  -1.3962987432107976e-7_dp,  1.0e-2_dp]
    tol_c2 = eps_test_c2 * sqrt(real(chk_step + 1, dp))
    resid_c2 = max(maxval(abs(traj_q(:, chk_idx) - q_ref)), &
                   maxval(abs(traj_om(:, chk_idx) - om_ref)))
    call report("T4e C2 step=12000 max|delta| < eps*sqrt(12001)", &
                resid_c2 < tol_c2, resid_c2)

    !-----------------------------------------------------------------------
    ! JSON output for cross-test harness
    !-----------------------------------------------------------------------
    call write_json(traj_q, traj_om)

    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_symmetric_top: ", n_pass, " PASS / ", n_fail, " FAIL"
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
    ! write_json: emit 201 trajectory samples to _audit/fortran_symmetric_top_output.json.
    ! Format mirrors tests/fixtures/symmetric_top_trajectory.json for diff comparison.
    !-----------------------------------------------------------------------
    subroutine write_json(q_arr, om_arr)
        real(dp), intent(in) :: q_arr(4, n_samp), om_arr(3, n_samp)
        integer  :: junit, ios, i, step_i
        real(dp) :: t_i

        open(newunit=junit, file="_audit/fortran_symmetric_top_output.json", &
             status="replace", iostat=ios)
        if (ios /= 0) then
            write (*, '(a)') "  WARNING: could not open _audit/fortran_symmetric_top_output.json"
            return
        end if

        write (junit, '(a)') '{'
        write (junit, '(a)') '  "spec_id": "fortran_symmetric_top_output.v0.1",'
        write (junit, '(a)') '  "trajectory": ['
        do i = 1, n_samp
            step_i = (i - 1) * sample_ev
            t_i    = real(step_i, dp) * dt_s
            write (junit, '(a,i0,a,es22.15,a)') &
                '    {"step": ', step_i, ', "t": ', t_i, ','
            write (junit, '(a,es22.15,a,es22.15,a,es22.15,a,es22.15,a)') &
                '     "q": [', q_arr(1,i), ', ', q_arr(2,i), ', ', &
                               q_arr(3,i), ', ', q_arr(4,i), '],'
            write (junit, '(a,es22.15,a,es22.15,a,es22.15,a)') &
                '     "omega_body": [', om_arr(1,i), ', ', om_arr(2,i), &
                                 ', ', om_arr(3,i), ']'
            if (i < n_samp) then
                write (junit, '(a)') '    },'
            else
                write (junit, '(a)') '    }'
            end if
        end do
        write (junit, '(a)') '  ]'
        write (junit, '(a)') '}'
        close(junit)
    end subroutine write_json

end program test_symmetric_top
