! jwst_l2_diagnostics
! --------------------
! Conservation-law diagnostics for the two-body coupled dynamics system.
! Ports dynamics.py diagnostic functions (lines 196-249).
!
! Two-body state layout (26 components, 1-based Fortran indexing):
!   Body A: x(1:3), v(4:6), q(7:10), omega(11:13)
!   Body B: x(14:16), v(17:19), q(20:23), omega(24:26)
!
! Canonical references:
!   - Chapter 02 Sec.2.3 (translational KE + König decomposition)
!   - Chapter 02 Sec.2.5 (inertia tensor in inertial frame: I_inertial = R . I . R^T)
!   - Chapter 02 Sec.2.5.1 (boxed dI_inertial/dt commutator formula)
! Engineering references:
!   - Hughes 2004 Sec.3.2 (angular momentum decomposition)
!   - dynamics.py lines 196-249 (Python oracle)
!
! Author: Sonnet Phase A 2026-05-24
! License: Apache 2.0 (inherited from parent project)

module jwst_l2_diagnostics
    use jwst_l2_constants, only: dp, g_newton
    use jwst_l2_quaternion, only: q_to_matrix
    implicit none
    private

    public :: total_kinetic_energy
    public :: total_potential_energy
    public :: total_angular_momentum
    public :: inertia_inertial_frame
    public :: inertial_inertia_rate

contains

    !-----------------------------------------------------------------------
    ! total_kinetic_energy
    !   (state, m_a, I_a_body, m_b, I_b_body)
    !   -> (E_trans, E_rot, E_total)
    !
    ! Port of dynamics.py::total_kinetic_energy (lines 196-204).
    !   E_trans = 0.5 * m_A * |v_A|^2 + 0.5 * m_B * |v_B|^2
    !   E_rot   = 0.5 * omega_A . I_A . omega_A + 0.5 * omega_B . I_B . omega_B
    !-----------------------------------------------------------------------
    pure subroutine total_kinetic_energy(state, m_a, I_a_body, m_b, I_b_body, &
                                          E_trans, E_rot, E_total)
        real(dp), intent(in)  :: state(26), m_a, m_b
        real(dp), intent(in)  :: I_a_body(3, 3), I_b_body(3, 3)
        real(dp), intent(out) :: E_trans, E_rot, E_total
        real(dp)              :: v_a(3), v_b(3), omega_a(3), omega_b(3)
        real(dp)              :: Iomega_a(3), Iomega_b(3)

        v_a    = state(4:6)
        v_b    = state(17:19)
        omega_a = state(11:13)
        omega_b = state(24:26)

        E_trans = 0.5_dp * m_a * dot3(v_a, v_a) + 0.5_dp * m_b * dot3(v_b, v_b)

        Iomega_a = matmul(I_a_body, omega_a)
        Iomega_b = matmul(I_b_body, omega_b)
        E_rot = 0.5_dp * dot3(omega_a, Iomega_a) + 0.5_dp * dot3(omega_b, Iomega_b)

        E_total = E_trans + E_rot
    end subroutine total_kinetic_energy

    !-----------------------------------------------------------------------
    ! total_potential_energy(state, m_a, m_b) -> V
    !
    ! Newtonian gravitational PE (point-mass monopole leading order).
    ! Port of dynamics.py::total_potential_energy (lines 207-214).
    !-----------------------------------------------------------------------
    pure function total_potential_energy(state, m_a, m_b) result(V)
        real(dp), intent(in) :: state(26), m_a, m_b
        real(dp)             :: V
        real(dp)             :: x_a(3), x_b(3), r_mag

        x_a = state(1:3)
        x_b = state(14:16)
        r_mag = sqrt((x_a(1)-x_b(1))**2 + (x_a(2)-x_b(2))**2 + (x_a(3)-x_b(3))**2)
        if (r_mag < 1.0e-12_dp) then
            V = 0.0_dp
        else
            V = -g_newton * m_a * m_b / r_mag
        end if
    end function total_potential_energy

    !-----------------------------------------------------------------------
    ! total_angular_momentum
    !   (state, m_a, I_a_body, m_b, I_b_body) -> L_total(3)
    !
    ! Total system angular momentum in the inertial frame.
    ! L = L_orbital_A + L_orbital_B + L_spin_A + L_spin_B
    ! L_spin_X = R_X . (I_X_body . omega_X)
    !
    ! Port of dynamics.py::total_angular_momentum (lines 217-233).
    !-----------------------------------------------------------------------
    pure function total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body) result(L_total)
        real(dp), intent(in) :: state(26), m_a, m_b
        real(dp), intent(in) :: I_a_body(3, 3), I_b_body(3, 3)
        real(dp)             :: L_total(3)
        real(dp)             :: x_a(3), v_a(3), q_a(4), omega_a(3)
        real(dp)             :: x_b(3), v_b(3), q_b(4), omega_b(3)
        real(dp)             :: R_a(3, 3), R_b(3, 3)
        real(dp)             :: L_orb_a(3), L_orb_b(3), L_spin_a(3), L_spin_b(3)

        x_a = state(1:3);   v_a = state(4:6)
        q_a = state(7:10);  omega_a = state(11:13)
        x_b = state(14:16); v_b = state(17:19)
        q_b = state(20:23); omega_b = state(24:26)

        L_orb_a = m_a * cross3(x_a, v_a)
        L_orb_b = m_b * cross3(x_b, v_b)

        R_a = q_to_matrix(q_a)
        R_b = q_to_matrix(q_b)
        L_spin_a = matmul(R_a, matmul(I_a_body, omega_a))
        L_spin_b = matmul(R_b, matmul(I_b_body, omega_b))

        L_total = L_orb_a + L_orb_b + L_spin_a + L_spin_b
    end function total_angular_momentum

    !-----------------------------------------------------------------------
    ! inertia_inertial_frame(I_body, q) -> I_inertial(3,3)
    !
    ! Transform inertia tensor from body frame to inertial frame:
    !   I_inertial = R . I_body . R^T
    !
    ! Port of dynamics.py::inertia_inertial_frame (lines 236-240).
    !-----------------------------------------------------------------------
    pure function inertia_inertial_frame(I_body, q) result(I_inertial)
        real(dp), intent(in) :: I_body(3, 3), q(4)
        real(dp)             :: I_inertial(3, 3)
        real(dp)             :: R(3, 3)

        R = q_to_matrix(q)
        I_inertial = matmul(R, matmul(I_body, transpose(R)))
    end function inertia_inertial_frame

    !-----------------------------------------------------------------------
    ! inertial_inertia_rate
    !   (I_body, omega_body, R_body_to_inertial) -> dI_inertial_dt(3,3)
    !
    ! Time derivative of the inertia tensor in the inertial frame via the
    ! matrix commutator (Chapter 02 Sec.2.5.1 boxed equation):
    !
    !   dI_inertial/dt = hat(omega_inertial) . I_inertial - I_inertial . hat(omega_inertial)
    !                  = [hat(omega_inertial), I_inertial]
    !
    ! where omega_inertial = R . omega_body and I_inertial = R . I_body . R^T.
    !
    ! Canonical reference: Chapter 02 Sec.2.5.1 (boxed commutator formula).
    !-----------------------------------------------------------------------
    pure function inertial_inertia_rate(I_body, omega_body, R_body_to_inertial) &
                                        result(dI_dt)
        real(dp), intent(in) :: I_body(3, 3)
        real(dp), intent(in) :: omega_body(3)
        real(dp), intent(in) :: R_body_to_inertial(3, 3)
        real(dp)             :: dI_dt(3, 3)
        real(dp)             :: omega_inertial(3), I_inertial(3, 3)
        real(dp)             :: hat_omega(3, 3)

        omega_inertial = matmul(R_body_to_inertial, omega_body)
        I_inertial = matmul(R_body_to_inertial, matmul(I_body, transpose(R_body_to_inertial)))

        ! hat(omega) = antisymmetric matrix such that hat(w).v = w x v
        hat_omega(1, 1) =  0.0_dp
        hat_omega(1, 2) = -omega_inertial(3)
        hat_omega(1, 3) =  omega_inertial(2)
        hat_omega(2, 1) =  omega_inertial(3)
        hat_omega(2, 2) =  0.0_dp
        hat_omega(2, 3) = -omega_inertial(1)
        hat_omega(3, 1) = -omega_inertial(2)
        hat_omega(3, 2) =  omega_inertial(1)
        hat_omega(3, 3) =  0.0_dp

        ! Commutator [hat_omega, I_inertial]
        dI_dt = matmul(hat_omega, I_inertial) - matmul(I_inertial, hat_omega)
    end function inertial_inertia_rate

    !-----------------------------------------------------------------------
    ! dot3(a, b) -> scalar dot product  (private)
    !-----------------------------------------------------------------------
    pure function dot3(a, b) result(d)
        real(dp), intent(in) :: a(3), b(3)
        real(dp)             :: d
        d = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
    end function dot3

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

end module jwst_l2_diagnostics
