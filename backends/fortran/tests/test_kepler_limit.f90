! test_kepler_limit  (D7)
! -----------------------
! Validates jwst_l2_dynamics in the Kepler limit:
!   I_A = I_B = eps*I_3 (isotropic tiny value -> zero gravity-gradient torque,
!   non-singular solve3x3).
!
! Configuration:
!   m_A = 5.972e24 kg (Earth-mass primary)
!   m_B = 7.350e22 kg (Moon-mass secondary)
!   a = 384400 km, eccentricity e = 0.1 (non-zero LRL for LRL-conservation test)
!   N = 1000 RK4 steps over one analytical orbital period T_orbital.
!
! Tests:
!   T1: Orbital period / position closure.
!       After N steps = 1 analytical period, |r_rel_final - r_rel_init| / r_p < 1e-6.
!   T2: Total mechanical energy conserved to 1e-10 relative.
!   T3: LRL vector magnitude conserved to 1e-8 relative.
!       (LRL = v_rel x (r_rel x v_rel) - mu * r_hat; |LRL| = mu*e for Kepler)
!
! Exit 0 on all PASS, exit 1 on any FAIL.
! Canonical reference: Chapter 03 Sec.3.8.1 (Kepler reduction).
!
! Author: Sonnet Phase A 2026-05-24

program test_kepler_limit
    use jwst_l2_constants,   only: dp, g_newton
    use jwst_l2_dynamics,    only: body_state_t, pack_state, unpack_state, rk4_step
    use jwst_l2_diagnostics, only: total_kinetic_energy, total_potential_energy
    implicit none

    integer :: n_pass, n_fail
    n_pass = 0
    n_fail = 0

    block
        real(dp) :: m_a, m_b, a_orb, ecc, mu
        real(dp) :: r_p, v_p, T_orbit, dt
        real(dp) :: I_eps(3,3), state(26)
        type(body_state_t) :: ba_0, bb_0, ba_f, bb_f
        real(dp) :: E_tr, E_rot, E_ke, V_pe, E_mech_0, E_mech_f
        real(dp) :: r_rel_0(3), r_rel_f(3), r_p_0, dr_close
        real(dp) :: lrl_0, lrl_f, lrl_rdiff
        integer  :: nstep, istep

        ! ----- Physical constants -----
        m_a   = 5.972e24_dp
        m_b   = 7.350e22_dp
        a_orb = 3.844e8_dp     ! 384 400 km
        ecc   = 0.1_dp         ! eccentricity (non-zero for LRL test)
        mu    = g_newton * (m_a + m_b)

        r_p   = a_orb * (1.0_dp - ecc)             ! periapsis distance
        v_p   = sqrt(mu * (1.0_dp + ecc) / r_p)    ! periapsis speed
        T_orbit = 2.0_dp * 4.0_dp*atan(1.0_dp) * sqrt(a_orb**3 / mu)

        nstep = 1000
        dt    = T_orbit / real(nstep, dp)

        ! Isotropic tiny inertia: gravity-gradient torque = 0, solve3x3 non-singular
        I_eps = 0.0_dp
        I_eps(1,1) = 1.0e-6_dp;  I_eps(2,2) = 1.0e-6_dp;  I_eps(3,3) = 1.0e-6_dp

        ! Initial state: periapsis in x, orbit in xy plane, CM at origin
        ba_0%x     = [-m_b/(m_a+m_b)*r_p, 0.0_dp, 0.0_dp]
        ba_0%v     = [0.0_dp, -m_b/(m_a+m_b)*v_p, 0.0_dp]
        ba_0%q     = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        ba_0%omega = [0.0_dp, 0.0_dp, 0.0_dp]

        bb_0%x     = [+m_a/(m_a+m_b)*r_p, 0.0_dp, 0.0_dp]
        bb_0%v     = [0.0_dp, +m_a/(m_a+m_b)*v_p, 0.0_dp]
        bb_0%q     = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        bb_0%omega = [0.0_dp, 0.0_dp, 0.0_dp]

        state = pack_state(ba_0, bb_0)

        ! Initial mechanical energy
        call total_kinetic_energy(state, m_a, I_eps, m_b, I_eps, E_tr, E_rot, E_ke)
        V_pe = total_potential_energy(state, m_a, m_b)
        E_mech_0 = E_ke + V_pe

        ! Initial LRL magnitude: |LRL| = mu * ecc for Kepler
        r_rel_0 = bb_0%x - ba_0%x
        r_p_0   = sqrt(dot3(r_rel_0, r_rel_0))
        lrl_0   = lrl_mag(ba_0%v, bb_0%v, r_rel_0, mu)

        ! Integrate one orbital period
        do istep = 1, nstep
            state = rk4_step(state, dt, m_a, I_eps, m_b, I_eps)
        end do

        call unpack_state(state, ba_f, bb_f)

        ! Final energy
        call total_kinetic_energy(state, m_a, I_eps, m_b, I_eps, E_tr, E_rot, E_ke)
        V_pe = total_potential_energy(state, m_a, m_b)
        E_mech_f = E_ke + V_pe

        ! Final relative position and LRL
        r_rel_f = bb_f%x - ba_f%x
        lrl_f   = lrl_mag(ba_f%v, bb_f%v, r_rel_f, mu)

        !-------------------------------------------------------------------
        ! T1: Position closure (period check)
        ! After one analytical period the orbit should have closed to < 1e-6
        !-------------------------------------------------------------------
        dr_close = sqrt((r_rel_f(1)-r_rel_0(1))**2 + &
                        (r_rel_f(2)-r_rel_0(2))**2 + &
                        (r_rel_f(3)-r_rel_0(3))**2) / r_p_0
        call report("T1 position closure |dr|/r_p < 1e-6 (period check)", &
                    dr_close < 1.0e-6_dp, dr_close)

        !-------------------------------------------------------------------
        ! T2: Total mechanical energy conserved to 1e-10
        !-------------------------------------------------------------------
        call report("T2 energy |dE/E_0| < 1e-10", &
                    abs(E_mech_f - E_mech_0) / abs(E_mech_0) < 1.0e-10_dp, &
                    abs(E_mech_f - E_mech_0) / abs(E_mech_0))

        !-------------------------------------------------------------------
        ! T3: LRL vector magnitude conserved to 1e-8
        !-------------------------------------------------------------------
        lrl_rdiff = abs(lrl_f - lrl_0) / max(lrl_0, 1.0e-30_dp)
        call report("T3 LRL magnitude |d|LRL||/|LRL| < 1e-8", &
                    lrl_rdiff < 1.0e-8_dp, lrl_rdiff)
    end block

    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_kepler_limit: ", n_pass, " PASS / ", n_fail, " FAIL"
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
    ! lrl_mag: reduced-mass LRL vector magnitude for the two-body system.
    !   A = v_rel x (r_rel x v_rel) - mu * r_hat
    !   |A| = mu * ecc  for Kepler orbit.
    !-----------------------------------------------------------------------
    pure function lrl_mag(v_a, v_b, r_rel, mu_grav) result(lrl_norm)
        real(dp), intent(in) :: v_a(3), v_b(3), r_rel(3), mu_grav
        real(dp)             :: lrl_norm
        real(dp) :: v_rel(3), r_mag, r_hat(3), L_vec(3), A_vec(3)

        v_rel  = v_b - v_a
        r_mag  = sqrt(r_rel(1)**2 + r_rel(2)**2 + r_rel(3)**2)
        r_hat  = r_rel / r_mag
        L_vec  = cross3(r_rel, v_rel)
        A_vec  = cross3(v_rel, L_vec) - mu_grav * r_hat
        lrl_norm = sqrt(A_vec(1)**2 + A_vec(2)**2 + A_vec(3)**2)
    end function lrl_mag

    pure function dot3(a, b) result(d)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: d
        d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
    end function dot3

    pure function cross3(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross3

end program test_kepler_limit
