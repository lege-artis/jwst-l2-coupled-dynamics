! test_orientation_jacobian  (D5)
! --------------------------------
! Tests orientation_response_jacobian from jwst_l2_jacobian.
!
! Background: the canonical-tier §3.7.2 boxed formula has a sign typo —
! the derivation correctly produces J*a = -alpha*[(a×u)×(I*u) + u×I*(a×u)],
! but the box was transcribed with +alpha.  The module was fixed to -alpha
! before this test was authored.  The finite-difference column T1 serves as
! the primary regression guard for that sign.
!
! Tests:
!   T1: FD cross-check.  I=diag(1000,2000,3000), r_rel=[1e8,0,1.5e9], R=I.
!       Central difference at eps=1e-6; max|J_anal-J_fd|/max|J_anal| < 5e-9.
!   T2a: Libration-equilibrium diagonal.  Oblate I=diag(2540,2540,3870),
!       u=e_3.  Assert J(1,1)=J(2,2)=alpha*(I_par-I_perp), all other
!       elements zero, to eps_test_c2 * |J11|.
!   T2b: J(3,3)=0 exactly at equilibrium (symmetry-axis rotation).
!   T3a: Spherical body I=diag(I0,I0,I0): J = 0 identically (exact).
!   T3b: Null-vector property J*u = 0 for generic asymmetric body.
!       u=[1,0,15]/norm; |J*u|/scale < eps_test_c2.
!
! Note: J always has u in its null space (rank <= 2), so no rank-3 check.
!
! Exit 0 on all PASS, exit 1 on any FAIL.
! Canonical reference: Chapter 03 Sec.3.7.2.
! Python oracle: N/A (no C2 gate on this module).
!
! Author: Sonnet Phase A 2026-05-24

program test_orientation_jacobian
    use jwst_l2_constants, only: dp, g_newton, eps_test_c2
    use jwst_l2_jacobian,  only: orientation_response_jacobian
    use jwst_l2_torque,    only: gravity_gradient_torque_body_frame
    implicit none

    integer :: n_pass, n_fail
    n_pass = 0
    n_fail = 0

    !-----------------------------------------------------------------------
    ! T1: Finite-difference cross-check
    ! I = diag(1000, 2000, 3000), r_rel = [1e8, 0, 1.5e9], R_nom = I_3
    !-----------------------------------------------------------------------
    block
        real(dp) :: I_body(3, 3), r_rel(3)
        real(dp) :: R_pm(3, 3)
        real(dp) :: tau_p(3), tau_m(3), J_fd(3, 3), J_an(3, 3)
        real(dp) :: rho_hat(3), r_mag, alpha, m_oth, eps_fd
        real(dp) :: scale_j, max_rdiff

        I_body    = 0.0_dp
        I_body(1,1) = 1000.0_dp
        I_body(2,2) = 2000.0_dp
        I_body(3,3) = 3000.0_dp
        r_rel = [1.0e8_dp, 0.0_dp, 1.5e9_dp]
        m_oth = 5.972e24_dp
        r_mag = sqrt(r_rel(1)**2 + r_rel(2)**2 + r_rel(3)**2)
        rho_hat = r_rel / r_mag
        alpha   = 3.0_dp * g_newton * m_oth / r_mag**3

        J_an = orientation_response_jacobian(I_body, rho_hat, alpha)

        eps_fd = 1.0e-6_dp

        ! Column 1: R = I3 +/- eps * hat(e1), hat(e1)_23=-1, hat(e1)_32=+1
        R_pm = identity3()
        R_pm(2,3) = -eps_fd;  R_pm(3,2) =  eps_fd
        tau_p = gravity_gradient_torque_body_frame(I_body, r_rel, R_pm, m_oth)
        R_pm = identity3()
        R_pm(2,3) =  eps_fd;  R_pm(3,2) = -eps_fd
        tau_m = gravity_gradient_torque_body_frame(I_body, r_rel, R_pm, m_oth)
        J_fd(:,1) = (tau_p - tau_m) / (2.0_dp * eps_fd)

        ! Column 2: R = I3 +/- eps * hat(e2), hat(e2)_13=+1, hat(e2)_31=-1
        R_pm = identity3()
        R_pm(1,3) =  eps_fd;  R_pm(3,1) = -eps_fd
        tau_p = gravity_gradient_torque_body_frame(I_body, r_rel, R_pm, m_oth)
        R_pm = identity3()
        R_pm(1,3) = -eps_fd;  R_pm(3,1) =  eps_fd
        tau_m = gravity_gradient_torque_body_frame(I_body, r_rel, R_pm, m_oth)
        J_fd(:,2) = (tau_p - tau_m) / (2.0_dp * eps_fd)

        ! Column 3: R = I3 +/- eps * hat(e3), hat(e3)_12=-1, hat(e3)_21=+1
        R_pm = identity3()
        R_pm(1,2) = -eps_fd;  R_pm(2,1) =  eps_fd
        tau_p = gravity_gradient_torque_body_frame(I_body, r_rel, R_pm, m_oth)
        R_pm = identity3()
        R_pm(1,2) =  eps_fd;  R_pm(2,1) = -eps_fd
        tau_m = gravity_gradient_torque_body_frame(I_body, r_rel, R_pm, m_oth)
        J_fd(:,3) = (tau_p - tau_m) / (2.0_dp * eps_fd)

        scale_j   = maxval(abs(J_an))
        max_rdiff = maxval(abs(J_an - J_fd)) / max(scale_j, 1.0e-30_dp)
        call report("T1 FD cross-check max|J_an-J_fd|/max|J_an| < 5e-9", &
                    max_rdiff < 5.0e-9_dp, max_rdiff)
    end block

    !-----------------------------------------------------------------------
    ! T2a/T2b: Libration-equilibrium eigenvalues
    ! Oblate I = diag(I_perp, I_perp, I_par), u = e_3
    ! Expected: J = alpha*(I_par-I_perp)*diag(1,1,0)
    !-----------------------------------------------------------------------
    block
        real(dp) :: I_body(3, 3), u_eq(3), J_eq(3, 3)
        real(dp) :: I_perp, I_par, alpha_eq, m_eq, r_eq
        real(dp) :: J11_exp, max_diff, tol

        I_perp = 2540.0_dp
        I_par  = 3870.0000000000005_dp   ! exact fixture value from case 6/7
        m_eq   = 5.972e24_dp
        r_eq   = 1.5e9_dp
        alpha_eq = 3.0_dp * g_newton * m_eq / r_eq**3

        I_body = 0.0_dp
        I_body(1,1) = I_perp;  I_body(2,2) = I_perp;  I_body(3,3) = I_par
        u_eq = [0.0_dp, 0.0_dp, 1.0_dp]

        J_eq = orientation_response_jacobian(I_body, u_eq, alpha_eq)

        J11_exp = alpha_eq * (I_par - I_perp)
        tol = max(eps_test_c2 * abs(J11_exp), eps_test_c2)

        max_diff = max(abs(J_eq(1,1) - J11_exp), &
                       abs(J_eq(2,2) - J11_exp), &
                       abs(J_eq(3,3)),            &
                       abs(J_eq(1,2)), abs(J_eq(1,3)), &
                       abs(J_eq(2,1)), abs(J_eq(2,3)), &
                       abs(J_eq(3,1)), abs(J_eq(3,2)))
        call report("T2a eq J = alpha*(I_par-I_perp)*diag(1,1,0)", &
                    max_diff < tol, max_diff)

        call report("T2b eq J(3,3) = 0 exactly (symmetry-axis rotation)", &
                    abs(J_eq(3,3)) < 1.0e-30_dp, abs(J_eq(3,3)))
    end block

    !-----------------------------------------------------------------------
    ! T3a: Spherical body I = I0*I_3 -> J = 0 identically
    !-----------------------------------------------------------------------
    block
        real(dp) :: I_body(3, 3), u_arb(3), J_sph(3, 3)
        real(dp) :: u_mag, alpha_val, max_abs_j

        I_body = 0.0_dp
        I_body(1,1) = 5000.0_dp;  I_body(2,2) = 5000.0_dp;  I_body(3,3) = 5000.0_dp
        u_arb  = [1.0_dp, 1.0_dp, 1.0_dp]
        u_mag  = sqrt(3.0_dp)
        u_arb  = u_arb / u_mag
        alpha_val = 3.541e-13_dp

        J_sph = orientation_response_jacobian(I_body, u_arb, alpha_val)
        max_abs_j = maxval(abs(J_sph))
        call report("T3a spherical body J = 0 identically", &
                    max_abs_j < 1.0e-30_dp, max_abs_j)
    end block

    !-----------------------------------------------------------------------
    ! T3b: Null-vector property J*u = 0 for asymmetric body
    ! u = [1, 0, 15] / norm (same direction as T1 r_rel, rh2=0 so cancels
    ! exactly in floating point).
    ! Verify |J*u| / (alpha * max_I) < eps_test_c2.
    !-----------------------------------------------------------------------
    block
        real(dp) :: I_body(3, 3), u_nv(3), J_nv(3, 3)
        real(dp) :: r_rel_nv(3), r_mag_nv, alpha_nv, m_nv
        real(dp) :: Ju(3), norm_Ju, scale_jnv

        I_body = 0.0_dp
        I_body(1,1) = 1000.0_dp;  I_body(2,2) = 2000.0_dp;  I_body(3,3) = 3000.0_dp
        r_rel_nv = [1.0e8_dp, 0.0_dp, 1.5e9_dp]
        m_nv  = 5.972e24_dp
        r_mag_nv = sqrt(r_rel_nv(1)**2 + r_rel_nv(2)**2 + r_rel_nv(3)**2)
        u_nv     = r_rel_nv / r_mag_nv
        alpha_nv = 3.0_dp * g_newton * m_nv / r_mag_nv**3

        J_nv = orientation_response_jacobian(I_body, u_nv, alpha_nv)

        ! J * u  (matrix-vector product)
        Ju(1) = J_nv(1,1)*u_nv(1) + J_nv(1,2)*u_nv(2) + J_nv(1,3)*u_nv(3)
        Ju(2) = J_nv(2,1)*u_nv(1) + J_nv(2,2)*u_nv(2) + J_nv(2,3)*u_nv(3)
        Ju(3) = J_nv(3,1)*u_nv(1) + J_nv(3,2)*u_nv(2) + J_nv(3,3)*u_nv(3)

        norm_Ju  = sqrt(Ju(1)**2 + Ju(2)**2 + Ju(3)**2)
        scale_jnv = alpha_nv * 3000.0_dp   ! alpha * max(I_body diagonal)
        call report("T3b null-vector |J*u|/scale < eps_test_c2", &
                    norm_Ju / scale_jnv < eps_test_c2, norm_Ju / scale_jnv)
    end block

    !-----------------------------------------------------------------------
    ! Summary
    !-----------------------------------------------------------------------
    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_orientation_jacobian: ", n_pass, " PASS / ", n_fail, " FAIL"
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

    pure function identity3() result(R)
        real(dp) :: R(3, 3)
        R = 0.0_dp
        R(1,1) = 1.0_dp;  R(2,2) = 1.0_dp;  R(3,3) = 1.0_dp
    end function identity3

end program test_orientation_jacobian
