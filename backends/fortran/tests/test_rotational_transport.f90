! test_rotational_transport  (D6)
! ---------------------------------
! Tests inertial_inertia_rate from jwst_l2_diagnostics:
!   dI_inertial/dt = [hat(omega_inertial), I_inertial]
! (Chapter 02 Sec.2.5.1 boxed commutator formula).
!
! Tests:
!   T1: Numerical-differentiation cross-check.
!       Propagate R(t) via Rodrigues formula at +/-dt (dt=1e-4 s).
!       Compute I_inertial(+dt) and I_inertial(-dt) via R*I*R^T similarity.
!       Central-difference dI/dt vs commutator formula at t=0.
!       Body: I=diag(1000,2000,3000), omega=[0.01,0.02,0.03] rad/s, R_0=I_3.
!       Tolerance: max|dI_num-dI_comm|/max|dI_comm| < 1e-10.
!   T2: Off-diagonal rate formula check (Sec.2.5.1 numerical sanity).
!       I=diag(2540,2540,3870), omega=[0.01,0,0] rad/s, R_0=I_3.
!       Assert (dI/dt)_{2,3} = omega_1*(I_2-I_3) = 0.01*(-1330) = -13.3
!       and all other elements within eps_test_c2.
!
! Exit 0 on all PASS, exit 1 on any FAIL.
! Canonical reference: Chapter 02 Sec.2.5.1 (boxed commutator formula).
! Python oracle: dynamics.py::inertial_inertia_rate (lines 236-248).
!
! Author: Sonnet Phase A 2026-05-24

program test_rotational_transport
    use jwst_l2_constants,   only: dp, eps_test_c2
    use jwst_l2_diagnostics, only: inertial_inertia_rate
    implicit none

    integer :: n_pass, n_fail
    n_pass = 0
    n_fail = 0

    !-----------------------------------------------------------------------
    ! T1: Central-difference cross-check
    ! I_body = diag(1000, 2000, 3000), omega_body = [0.01, 0.02, 0.03], R0 = I
    !-----------------------------------------------------------------------
    block
        real(dp) :: I_body(3,3), omega_body(3), R_0(3,3)
        real(dp) :: R_p(3,3), R_m(3,3), I_p(3,3), I_m(3,3)
        real(dp) :: dI_num(3,3), dI_comm(3,3)
        real(dp) :: dt, scale_di, max_rdiff

        I_body    = 0.0_dp
        I_body(1,1) = 1000.0_dp;  I_body(2,2) = 2000.0_dp;  I_body(3,3) = 3000.0_dp
        omega_body  = [0.01_dp, 0.02_dp, 0.03_dp]
        R_0         = eye3()
        dt          = 1.0e-4_dp

        ! Propagate: R(+dt) and R(-dt) via Rodrigues
        R_p = rodrigues(omega_body,  dt)
        R_m = rodrigues(omega_body, -dt)

        ! I_inertial(t) = R(t) * I_body * R(t)^T
        I_p = matmul(R_p, matmul(I_body, transpose(R_p)))
        I_m = matmul(R_m, matmul(I_body, transpose(R_m)))

        ! Numerical dI/dt via central difference
        dI_num = (I_p - I_m) / (2.0_dp * dt)

        ! Commutator formula at t=0 (R_0 = identity)
        dI_comm = inertial_inertia_rate(I_body, omega_body, R_0)

        scale_di  = maxval(abs(dI_comm))
        max_rdiff = maxval(abs(dI_num - dI_comm)) / max(scale_di, 1.0e-30_dp)
        call report("T1 FD vs commutator max|dI_num-dI_comm|/max|dI_comm| < 1e-10", &
                    max_rdiff < 1.0e-10_dp, max_rdiff)
    end block

    !-----------------------------------------------------------------------
    ! T2: Off-diagonal rate formula
    ! I = diag(2540, 2540, 3870), omega = [0.01, 0, 0], R = I
    ! Expected: (dI/dt)_{2,3} = (dI/dt)_{3,2} = omega_1*(I_2-I_3) = -13.3
    !           all other elements = 0
    !-----------------------------------------------------------------------
    block
        real(dp) :: I_body(3,3), omega_body(3), R_0(3,3), dI_dt(3,3)
        real(dp) :: expected_23, tol_off, tol_zero, max_zero_diff

        I_body    = 0.0_dp
        I_body(1,1) = 2540.0_dp;  I_body(2,2) = 2540.0_dp;  I_body(3,3) = 3870.0_dp
        omega_body  = [0.01_dp, 0.0_dp, 0.0_dp]
        R_0         = eye3()

        dI_dt = inertial_inertia_rate(I_body, omega_body, R_0)

        ! (dI/dt)_{2,3} = omega_1 * (I_2 - I_3) = 0.01 * (2540-3870) = -13.3
        expected_23 = omega_body(1) * (I_body(2,2) - I_body(3,3))
        tol_off  = max(eps_test_c2 * abs(expected_23), eps_test_c2)

        call report("T2a (dI/dt)_{2,3} = omega_1*(I_2-I_3) = -13.3", &
                    abs(dI_dt(2,3) - expected_23) < tol_off, &
                    abs(dI_dt(2,3) - expected_23))

        call report("T2b (dI/dt)_{3,2} = (dI/dt)_{2,3} (symmetric)", &
                    abs(dI_dt(3,2) - expected_23) < tol_off, &
                    abs(dI_dt(3,2) - expected_23))

        ! All non-(2,3) and non-(3,2) elements should be zero
        tol_zero = max(eps_test_c2 * abs(expected_23), eps_test_c2)
        max_zero_diff = max(abs(dI_dt(1,1)), abs(dI_dt(1,2)), abs(dI_dt(1,3)), &
                            abs(dI_dt(2,1)), abs(dI_dt(2,2)),                   &
                            abs(dI_dt(3,1)),                  abs(dI_dt(3,3)))
        call report("T2c all other dI/dt elements = 0", &
                    max_zero_diff < tol_zero, max_zero_diff)
    end block

    !-----------------------------------------------------------------------
    ! Summary
    !-----------------------------------------------------------------------
    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_rotational_transport: ", n_pass, " PASS / ", n_fail, " FAIL"
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

    pure function eye3() result(R)
        real(dp) :: R(3,3)
        R = 0.0_dp
        R(1,1) = 1.0_dp;  R(2,2) = 1.0_dp;  R(3,3) = 1.0_dp
    end function eye3

    !-----------------------------------------------------------------------
    ! rodrigues: exact rotation matrix exp(t * hat(omega)) via Rodrigues
    ! formula.  For t < 0, rotation is in the opposite direction.
    !-----------------------------------------------------------------------
    pure function rodrigues(omega_vec, t) result(R)
        real(dp), intent(in) :: omega_vec(3), t
        real(dp)             :: R(3,3)
        real(dp)             :: om_mag, theta, nx, ny, nz, c, s, nc

        om_mag = sqrt(omega_vec(1)**2 + omega_vec(2)**2 + omega_vec(3)**2)
        if (om_mag < 1.0e-30_dp) then
            R = eye3()
            return
        end if
        theta = om_mag * t   ! signed angle
        nx = omega_vec(1) / om_mag
        ny = omega_vec(2) / om_mag
        nz = omega_vec(3) / om_mag
        c  = cos(theta);  s = sin(theta);  nc = 1.0_dp - c

        R(1,1) = c + nc*nx*nx;    R(1,2) = nc*nx*ny - s*nz; R(1,3) = nc*nx*nz + s*ny
        R(2,1) = nc*ny*nx + s*nz; R(2,2) = c + nc*ny*ny;    R(2,3) = nc*ny*nz - s*nx
        R(3,1) = nc*nz*nx - s*ny; R(3,2) = nc*nz*ny + s*nx; R(3,3) = c + nc*nz*nz
    end function rodrigues

end program test_rotational_transport
