! jwst_l2_torque
! --------------
! Gravity-gradient torque on a rigid body due to a distant point mass.
! Line-for-line Fortran port of dynamics.py::gravity_gradient_torque_body_frame
! (lines 113-131).
!
! Formula (Hughes 2004 Sec.5.3):
!   tau_body = (3 G m_other / r^3) * r_hat_body x (I_body . r_hat_body)
! where r_hat_body = R^T . r_hat_inertial is the separation unit vector
! expressed in the body frame.
!
! Canonical reference: Chapter 03 Sec.3.7 (boxed equation).
! Engineering reference: Hughes 2004 Sec.5.3.
! Python oracle: dynamics.py::gravity_gradient_torque_body_frame (lines 113-131).
! C2 cross-test: tests/test_gravity_gradient_torque.f90.
!
! Conventions (locked per SONNET-BRIEF Sec.2):
!   - r_rel_inertial = x_other - x_self (points FROM self TO other)
!   - R_body_to_inertial maps body vectors to inertial: v_inertial = R . v_body
!   - Result tau_body is in the BODY frame (Euler's equations canonical form)
!
! Author: Sonnet Phase A 2026-05-24
! License: Apache 2.0 (inherited from parent project)

module jwst_l2_torque
    use jwst_l2_constants, only: dp, g_newton
    implicit none
    private

    public :: gravity_gradient_torque_body_frame

contains

    !-----------------------------------------------------------------------
    ! gravity_gradient_torque_body_frame
    !   (I_body, r_rel_inertial, R_body_to_inertial, m_other) -> tau_body(3)
    !
    ! Exact Fortran translation of dynamics.py lines 113-131.
    ! All arithmetic operations are in the same order as Python to maximise
    ! bit-identity with the C2 cross-test oracle at ULP*sqrt(N) tolerance.
    !-----------------------------------------------------------------------
    pure function gravity_gradient_torque_body_frame(I_body, r_rel_inertial, &
                                                      R_body_to_inertial, m_other) &
                                                      result(tau_body)
        real(dp), intent(in) :: I_body(3, 3)
        real(dp), intent(in) :: r_rel_inertial(3)
        real(dp), intent(in) :: R_body_to_inertial(3, 3)
        real(dp), intent(in) :: m_other
        real(dp)             :: tau_body(3)
        real(dp)             :: r_mag, r_hat_inertial(3), r_hat_body(3), I_rhat(3)

        r_mag = sqrt(r_rel_inertial(1)**2 + r_rel_inertial(2)**2 + r_rel_inertial(3)**2)
        if (r_mag < 1.0e-12_dp) then
            tau_body = 0.0_dp
            return
        end if
        r_hat_inertial = r_rel_inertial / r_mag
        r_hat_body = matmul(transpose(R_body_to_inertial), r_hat_inertial)
        I_rhat = matmul(I_body, r_hat_body)
        tau_body = (3.0_dp * g_newton * m_other / r_mag**3) * cross3(r_hat_body, I_rhat)
    end function gravity_gradient_torque_body_frame

    !-----------------------------------------------------------------------
    ! cross3(a, b) -> c = a x b
    ! Private helper: 3-vector cross product.
    !-----------------------------------------------------------------------
    pure function cross3(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp)             :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross3

end module jwst_l2_torque
