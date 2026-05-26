! jwst_l2_jacobian
! ----------------
! Orientation-response Jacobian for the gravity-gradient torque on a rigid body.
!
! For a body at equilibrium orientation (or near it), the linearised torque
! response to a small attitude change delta-q is:
!
!   delta_tau_body = J_orient . delta_theta
!
! where J_orient is the 3x3 orientation-response Jacobian.
!
! Closed-form expression (Chapter 03 Sec.3.7.2, boxed equation):
!   J_orient . a = -alpha * [(a x u) x (I . u) + u x I . (a x u)]
! for any 3-vector a, where u = rho_hat_body and alpha = 3 G m_other / r^3.
!
! Sign convention note (OQ-FORT-1, 2026-05-24): the leading factor is -alpha,
! NOT +alpha.  The minus sign arises from delta_u = -delta_theta x u carried
! through the torque linearisation.  Sonnet's D5 T2a test (libration-equilibrium
! eigenvalues for axisymmetric oblate body) is the permanent regression guard:
! at u = (0,0,1) with I = diag(I_perp, I_perp, I_par), J_orient evaluates to
! diag(alpha*(I_par - I_perp), alpha*(I_par - I_perp), 0).  For oblate
! (I_par > I_perp) the eigenvalue must be POSITIVE (destabilising) -- this
! matches -alpha when combined with the on-axis algebraic identity.
!
! The 3x3 matrix is constructed column-by-column by applying the formula to
! each standard basis vector e_1, e_2, e_3.
!
! Canonical reference: Chapter 03 Sec.3.7.2 (boxed equation).
! Engineering reference: Hughes 2004 Sec.5.3 (gravity-gradient linearisation).
! C1 cross-test: tests/test_orientation_jacobian.f90.
!
! Author: Sonnet Phase A 2026-05-24
! License: Apache 2.0 (inherited from parent project)

module jwst_l2_jacobian
    use jwst_l2_constants, only: dp, g_newton
    implicit none
    private

    public :: orientation_response_jacobian

contains

    !-----------------------------------------------------------------------
    ! orientation_response_jacobian
    !   (I_body, rho_hat_body, alpha) -> J_orient(3,3)
    !
    ! Builds the 3x3 Jacobian matrix J such that:
    !   J . a = -alpha * [(a x u) x (I . u) + u x I . (a x u)]
    ! for arbitrary 3-vector a, where u = rho_hat_body.
    ! Note: leading sign is -alpha (NOT +alpha) -- see module-level comment
    ! and the D5 T2a regression guard.
    !
    ! alpha = 3 * G * m_other / r^3 (caller computes this).
    !
    ! Implementation: column j = J . e_j, where e_j is the j-th basis vector.
    ! This gives J as a full 3x3 matrix, inspectable and eigenvalue-extractable.
    !
    ! Canonical reference: Chapter 03 Sec.3.7.2 (boxed formula).
    !-----------------------------------------------------------------------
    pure function orientation_response_jacobian(I_body, rho_hat_body, alpha) result(J_orient)
        real(dp), intent(in) :: I_body(3, 3)
        real(dp), intent(in) :: rho_hat_body(3)
        real(dp), intent(in) :: alpha
        real(dp)             :: J_orient(3, 3)
        real(dp)             :: u(3), Iu(3), ej(3), axu(3), Iaxu(3), axu_x_Iu(3), u_x_Iaxu(3)
        integer              :: jcol

        u  = rho_hat_body
        Iu = matmul(I_body, u)

        do jcol = 1, 3
            ! Build basis vector e_jcol
            ej = 0.0_dp
            ej(jcol) = 1.0_dp

            ! a x u  (where a = e_j)
            axu = cross3(ej, u)
            ! I . (a x u)
            Iaxu = matmul(I_body, axu)
            ! (a x u) x (I . u)
            axu_x_Iu = cross3(axu, Iu)
            ! u x I . (a x u)
            u_x_Iaxu = cross3(u, Iaxu)

            J_orient(:, jcol) = -alpha * (axu_x_Iu + u_x_Iaxu)
        end do
    end function orientation_response_jacobian

    !-----------------------------------------------------------------------
    ! cross3(a, b) -> a x b  (private)
    !-----------------------------------------------------------------------
    pure function cross3(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp)             :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross3

end module jwst_l2_jacobian
