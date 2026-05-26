! jwst_l2_potential
! -----------------
! Multipole gravitational potential between two rigid bodies.
!
! V_mono  = -G * m_A * m_B / |rho|
! V_quad  = V_tidal_A + V_tidal_B
!   V_tidal_X = -(G * m_other / (2 * |rho|^3)) * [tr(I_X) - 3 * u_X . I_X . u_X]
!
! where rho = (x_B - x_A), u_A = R_A^T . rho_hat (direction to B in A body frame),
!       u_B = R_B^T . (-rho_hat) (direction to A in B body frame).
!
! Canonical reference: Chapter 03 Sec.3.4.3 (boxed quadrupole formula).
! Engineering reference: Hughes 2004 Sec.3 + Sec.5.
! Python oracle: dynamics.py::total_potential_energy (monopole), geometries.py
!   (inertia inputs); full multipole is new in this Fortran reference.
! C1 cross-test: tests/test_multipole_potential.f90.
!
! Author: Sonnet Phase A 2026-05-24
! License: Apache 2.0 (inherited from parent project)

module jwst_l2_potential
    use jwst_l2_constants, only: dp, g_newton
    use jwst_l2_quaternion, only: q_to_matrix
    implicit none
    private

    public :: quadrupole_from_inertia
    public :: mutual_potential_mono
    public :: mutual_potential_quad
    public :: total_mutual_potential

contains

    !-----------------------------------------------------------------------
    ! quadrupole_from_inertia(I_body) -> Q_body(3,3)
    !
    ! Trace-free quadrupole tensor: Q = tr(I) * I_3 - 3 * I.
    ! Used in the multipole expansion of the gravitational potential.
    !
    ! Canonical reference: Chapter 03 Sec.3.5 (trace identity).
    !-----------------------------------------------------------------------
    pure function quadrupole_from_inertia(I_body) result(Q_body)
        real(dp), intent(in) :: I_body(3, 3)
        real(dp)             :: Q_body(3, 3)
        real(dp)             :: tr_I
        integer              :: idiag

        tr_I = I_body(1, 1) + I_body(2, 2) + I_body(3, 3)
        Q_body = -3.0_dp * I_body
        do idiag = 1, 3
            Q_body(idiag, idiag) = Q_body(idiag, idiag) + tr_I
        end do
    end function quadrupole_from_inertia

    !-----------------------------------------------------------------------
    ! mutual_potential_mono(m_a, m_b, rho) -> V_mono
    !
    ! Newtonian monopole potential: V = -G * m_A * m_B / |rho|.
    ! rho = (x_B - x_A) is the separation vector in the inertial frame.
    !-----------------------------------------------------------------------
    pure function mutual_potential_mono(m_a, m_b, rho) result(V_mono)
        real(dp), intent(in) :: m_a, m_b, rho(3)
        real(dp)             :: V_mono
        real(dp)             :: r_mag

        r_mag = sqrt(rho(1)**2 + rho(2)**2 + rho(3)**2)
        if (r_mag < 1.0e-12_dp) then
            V_mono = 0.0_dp
        else
            V_mono = -g_newton * m_a * m_b / r_mag
        end if
    end function mutual_potential_mono

    !-----------------------------------------------------------------------
    ! mutual_potential_quad(m_a, m_b, rho, I_a_body, q_a, I_b_body, q_b)
    !                    -> V_quad
    !
    ! Quadrupole (tidal) correction to the mutual potential.
    ! Adds the tidal potential of body A in the field of B and vice versa.
    !
    !   u_A = R_A^T . rho_hat        (direction A -> B in A's body frame)
    !   u_B = R_B^T . (-rho_hat)     (direction B -> A in B's body frame)
    !
    !   V_tidal_A = -(G * m_B / (2 * r^3)) * [tr(I_A) - 3 * u_A . I_A . u_A]
    !   V_tidal_B = -(G * m_A / (2 * r^3)) * [tr(I_B) - 3 * u_B . I_B . u_B]
    !
    ! Canonical reference: Chapter 03 Sec.3.4.3 (boxed formula).
    !-----------------------------------------------------------------------
    pure function mutual_potential_quad(m_a, m_b, rho, &
                                         I_a_body, q_a, I_b_body, q_b) result(V_quad)
        real(dp), intent(in) :: m_a, m_b, rho(3)
        real(dp), intent(in) :: I_a_body(3, 3), q_a(4)
        real(dp), intent(in) :: I_b_body(3, 3), q_b(4)
        real(dp)             :: V_quad
        real(dp)             :: r_mag, rho_hat(3), r3
        real(dp)             :: R_a(3, 3), R_b(3, 3)
        real(dp)             :: u_a(3), u_b(3), Iu_a(3), Iu_b(3)
        real(dp)             :: tr_a, tr_b, quad_a, quad_b
        real(dp)             :: V_tidal_a, V_tidal_b

        r_mag = sqrt(rho(1)**2 + rho(2)**2 + rho(3)**2)
        if (r_mag < 1.0e-12_dp) then
            V_quad = 0.0_dp
            return
        end if

        rho_hat = rho / r_mag
        r3 = r_mag * r_mag * r_mag

        R_a = q_to_matrix(q_a)
        R_b = q_to_matrix(q_b)

        ! Direction A -> B in A's body frame
        u_a = matmul(transpose(R_a), rho_hat)
        ! Direction B -> A in B's body frame
        u_b = -matmul(transpose(R_b), rho_hat)

        Iu_a = matmul(I_a_body, u_a)
        Iu_b = matmul(I_b_body, u_b)

        tr_a   = I_a_body(1, 1) + I_a_body(2, 2) + I_a_body(3, 3)
        tr_b   = I_b_body(1, 1) + I_b_body(2, 2) + I_b_body(3, 3)
        quad_a = u_a(1)*Iu_a(1) + u_a(2)*Iu_a(2) + u_a(3)*Iu_a(3)
        quad_b = u_b(1)*Iu_b(1) + u_b(2)*Iu_b(2) + u_b(3)*Iu_b(3)

        V_tidal_a = -(g_newton * m_b / (2.0_dp * r3)) * (tr_a - 3.0_dp * quad_a)
        V_tidal_b = -(g_newton * m_a / (2.0_dp * r3)) * (tr_b - 3.0_dp * quad_b)

        V_quad = V_tidal_a + V_tidal_b
    end function mutual_potential_quad

    !-----------------------------------------------------------------------
    ! total_mutual_potential(...) -> V_total
    !
    ! Sum of monopole + quadrupole potential.
    !-----------------------------------------------------------------------
    pure function total_mutual_potential(m_a, m_b, rho, &
                                          I_a_body, q_a, I_b_body, q_b) result(V_total)
        real(dp), intent(in) :: m_a, m_b, rho(3)
        real(dp), intent(in) :: I_a_body(3, 3), q_a(4)
        real(dp), intent(in) :: I_b_body(3, 3), q_b(4)
        real(dp)             :: V_total

        V_total = mutual_potential_mono(m_a, m_b, rho) + &
                  mutual_potential_quad(m_a, m_b, rho, I_a_body, q_a, I_b_body, q_b)
    end function total_mutual_potential

end module jwst_l2_potential
