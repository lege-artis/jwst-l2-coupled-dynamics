! jwst_l2_dynamics
! ----------------
! Two-rigid-body coupled dynamics: state representation, equations of motion,
! and RK4 integrator. Line-for-line Fortran port of dynamics.py (lines 73-191).
!
! Two-body state vector (26 components, 1-based):
!   Body A: x(1:3), v(4:6), q(7:10), omega(11:13)
!   Body B: x(14:16), v(17:19), q(20:23), omega(24:26)
!
! Equations of motion:
!   Translation (Newton, inertial frame):
!     dx/dt = v
!     dv/dt = F_grav / m
!   Rotation (Euler, body frame):
!     I . d(omega)/dt = tau - omega x (I . omega)
!     dq/dt = 0.5 * Omega(omega) . q
!
! Canonical references:
!   - Chapter 01 Sec.1.4 (quaternion kinematics)
!   - Chapter 02 Sec.2.7 (Euler equations in body frame)
!   - Press et al. 2007 Numerical Recipes 3rd ed Sec.17.1 (RK4)
! Python oracle: dynamics.py lines 73-191.
! C2 cross-test: tests/test_symmetric_top.f90.
!
! Author: Sonnet Phase A 2026-05-24
! License: Apache 2.0 (inherited from parent project)

module jwst_l2_dynamics
    use jwst_l2_constants, only: dp, g_newton
    use jwst_l2_quaternion, only: q_to_matrix, q_kinematics_matrix, q_normalise
    use jwst_l2_torque,     only: gravity_gradient_torque_body_frame
    implicit none
    private

    ! Public derived type: single rigid-body state
    type, public :: body_state_t
        real(dp) :: x(3)      ! position (inertial frame)
        real(dp) :: v(3)      ! velocity (inertial frame)
        real(dp) :: q(4)      ! attitude quaternion (qw, qx, qy, qz)
        real(dp) :: omega(3)  ! angular velocity (body frame)
    end type body_state_t

    public :: pack_state
    public :: unpack_state
    public :: state_derivative
    public :: rk4_step

contains

    !-----------------------------------------------------------------------
    ! pack_state(body_a, body_b) -> state(26)
    !
    ! Flatten two body states into the 26-component vector.
    ! Port of dynamics.py::pack_state (line 93-94).
    !-----------------------------------------------------------------------
    pure function pack_state(body_a, body_b) result(state)
        type(body_state_t), intent(in) :: body_a, body_b
        real(dp)                       :: state(26)

        state(1:3)   = body_a%x
        state(4:6)   = body_a%v
        state(7:10)  = body_a%q
        state(11:13) = body_a%omega
        state(14:16) = body_b%x
        state(17:19) = body_b%v
        state(20:23) = body_b%q
        state(24:26) = body_b%omega
    end function pack_state

    !-----------------------------------------------------------------------
    ! unpack_state(state, body_a, body_b)
    !
    ! Extract two body states from the 26-component vector.
    ! Port of dynamics.py::unpack_state (lines 97-98).
    !-----------------------------------------------------------------------
    pure subroutine unpack_state(state, body_a, body_b)
        real(dp),           intent(in)  :: state(26)
        type(body_state_t), intent(out) :: body_a, body_b

        body_a%x     = state(1:3)
        body_a%v     = state(4:6)
        body_a%q     = state(7:10)
        body_a%omega = state(11:13)
        body_b%x     = state(14:16)
        body_b%v     = state(17:19)
        body_b%q     = state(20:23)
        body_b%omega = state(24:26)
    end subroutine unpack_state

    !-----------------------------------------------------------------------
    ! state_derivative(state, m_a, I_a_body, m_b, I_b_body) -> dsdt(26)
    !
    ! Compute dstate/dt for the 26-component coupled-body state.
    ! Port of dynamics.py::state_derivative (lines 136-171).
    !-----------------------------------------------------------------------
    pure function state_derivative(state, m_a, I_a_body, m_b, I_b_body) result(dsdt)
        real(dp), intent(in) :: state(26)
        real(dp), intent(in) :: m_a, m_b
        real(dp), intent(in) :: I_a_body(3, 3), I_b_body(3, 3)
        real(dp)             :: dsdt(26)

        type(body_state_t) :: body_a, body_b
        real(dp)           :: F_a(3), F_b(3)
        real(dp)           :: R_a(3, 3), R_b(3, 3)
        real(dp)           :: r_a_to_b(3), r_b_to_a(3)
        real(dp)           :: tau_a(3), tau_b(3)
        real(dp)           :: I_omega_a(3), I_omega_b(3)
        real(dp)           :: domega_a(3), domega_b(3)
        real(dp)           :: dq_a(4), dq_b(4)

        call unpack_state(state, body_a, body_b)

        ! Translational dynamics: F = -G*m_A*m_B*(x_A-x_B)/|r|^3
        F_a = gravitational_force_on_a(body_a%x, body_b%x, m_a, m_b)
        F_b = -F_a   ! Newton's 3rd law

        ! Rotational dynamics
        R_a = q_to_matrix(body_a%q)
        R_b = q_to_matrix(body_b%q)

        r_a_to_b = body_b%x - body_a%x
        r_b_to_a = -r_a_to_b

        tau_a = gravity_gradient_torque_body_frame(I_a_body, r_a_to_b, R_a, m_b)
        tau_b = gravity_gradient_torque_body_frame(I_b_body, r_b_to_a, R_b, m_a)

        ! Euler's equations: I . d(omega)/dt = tau - omega x (I . omega)
        I_omega_a = matmul(I_a_body, body_a%omega)
        I_omega_b = matmul(I_b_body, body_b%omega)
        call solve3x3(I_a_body, tau_a - cross3(body_a%omega, I_omega_a), domega_a)
        call solve3x3(I_b_body, tau_b - cross3(body_b%omega, I_omega_b), domega_b)

        ! Quaternion kinematics: dq/dt = 0.5 * Omega(omega) . q
        dq_a = 0.5_dp * matmul(q_kinematics_matrix(body_a%omega), body_a%q)
        dq_b = 0.5_dp * matmul(q_kinematics_matrix(body_b%omega), body_b%q)

        dsdt(1:3)   = body_a%v
        dsdt(4:6)   = F_a / m_a
        dsdt(7:10)  = dq_a
        dsdt(11:13) = domega_a
        dsdt(14:16) = body_b%v
        dsdt(17:19) = F_b / m_b
        dsdt(20:23) = dq_b
        dsdt(24:26) = domega_b
    end function state_derivative

    !-----------------------------------------------------------------------
    ! rk4_step(state, dt, m_a, I_a_body, m_b, I_b_body) -> new_state(26)
    !
    ! Classical 4th-order Runge-Kutta step. Quaternions renormalised after.
    ! Port of dynamics.py::rk4_step (lines 176-191).
    !
    ! Reference: Press et al. 2007 Numerical Recipes 3rd ed Sec.17.1, eq.17.1.3.
    !-----------------------------------------------------------------------
    pure function rk4_step(state, dt, m_a, I_a_body, m_b, I_b_body) result(new_state)
        real(dp), intent(in) :: state(26), dt, m_a, m_b
        real(dp), intent(in) :: I_a_body(3, 3), I_b_body(3, 3)
        real(dp)             :: new_state(26)
        real(dp)             :: k1(26), k2(26), k3(26), k4(26)

        k1 = state_derivative(state,                        m_a, I_a_body, m_b, I_b_body)
        k2 = state_derivative(state + 0.5_dp*dt*k1,        m_a, I_a_body, m_b, I_b_body)
        k3 = state_derivative(state + 0.5_dp*dt*k2,        m_a, I_a_body, m_b, I_b_body)
        k4 = state_derivative(state + dt*k3,                m_a, I_a_body, m_b, I_b_body)

        new_state = state + (dt / 6.0_dp) * (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)

        ! Renormalise quaternions (counters integration drift)
        new_state(7:10)  = q_normalise(new_state(7:10))
        new_state(20:23) = q_normalise(new_state(20:23))
    end function rk4_step

    !-----------------------------------------------------------------------
    ! gravitational_force_on_a (private)
    ! F_A = -G * m_A * m_B * (x_A - x_B) / |x_A - x_B|^3
    ! Port of dynamics.py::gravitational_force_on_a (lines 103-110).
    !-----------------------------------------------------------------------
    pure function gravitational_force_on_a(x_a, x_b, m_a, m_b) result(F_a)
        real(dp), intent(in) :: x_a(3), x_b(3), m_a, m_b
        real(dp)             :: F_a(3)
        real(dp)             :: r_vec(3), r_mag

        r_vec = x_a - x_b
        r_mag = sqrt(r_vec(1)**2 + r_vec(2)**2 + r_vec(3)**2)
        if (r_mag < 1.0e-12_dp) then
            F_a = 0.0_dp
        else
            F_a = -g_newton * m_a * m_b * r_vec / (r_mag**3)
        end if
    end function gravitational_force_on_a

    !-----------------------------------------------------------------------
    ! solve3x3(A, b, x)
    ! Solve the 3x3 linear system A . x = b by Gaussian elimination
    ! with partial pivoting. Used for the Euler-equation inversion
    ! I_body . domega = rhs.
    !
    ! A must be non-singular (inertia tensor is always SPD, so this holds).
    !-----------------------------------------------------------------------
    pure subroutine solve3x3(A, b, x)
        real(dp), intent(in)  :: A(3, 3), b(3)
        real(dp), intent(out) :: x(3)
        real(dp)              :: M(3, 4)   ! augmented matrix [A | b]
        real(dp)              :: factor, pivot_val, tmp
        integer               :: ipiv, irow, jcol, krow

        ! Build augmented matrix
        M(1, 1:3) = A(1, :); M(1, 4) = b(1)
        M(2, 1:3) = A(2, :); M(2, 4) = b(2)
        M(3, 1:3) = A(3, :); M(3, 4) = b(3)

        ! Forward elimination with partial pivoting
        do krow = 1, 2
            ! Find pivot
            ipiv = krow
            pivot_val = abs(M(krow, krow))
            do irow = krow + 1, 3
                if (abs(M(irow, krow)) > pivot_val) then
                    pivot_val = abs(M(irow, krow))
                    ipiv = irow
                end if
            end do
            ! Swap rows krow and ipiv
            if (ipiv /= krow) then
                do jcol = 1, 4
                    tmp          = M(krow, jcol)
                    M(krow, jcol) = M(ipiv, jcol)
                    M(ipiv, jcol) = tmp
                end do
            end if
            ! Eliminate
            do irow = krow + 1, 3
                factor = M(irow, krow) / M(krow, krow)
                M(irow, krow:4) = M(irow, krow:4) - factor * M(krow, krow:4)
            end do
        end do

        ! Back-substitution
        x(3) = M(3, 4) / M(3, 3)
        x(2) = (M(2, 4) - M(2, 3)*x(3)) / M(2, 2)
        x(1) = (M(1, 4) - M(1, 3)*x(3) - M(1, 2)*x(2)) / M(1, 1)
    end subroutine solve3x3

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

end module jwst_l2_dynamics
