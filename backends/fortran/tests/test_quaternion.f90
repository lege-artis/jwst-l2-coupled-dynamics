! test_quaternion
! ---------------
! Sanity tests for jwst_l2_quaternion module. NOT in the D1-D9 list per
! SONNET-BRIEF, but serves as the SMOKE TEST for the Makefile pattern and
! the quaternion primitives that subsequent tests depend on.
!
! Tests:
!   T1. hat_map(a) . b == cross(a, b) for random a, b
!   T2. q_to_matrix(identity_quaternion) == I_3
!   T3. q_to_matrix orthogonality: R . R^T == I_3 for unit quaternion
!   T4. q_kinematics_matrix antisymmetry: Omega^T == -Omega
!
! Exit 0 on all PASS, exit 1 on any FAIL.
!
! Author: Opus prep pass 2026-05-24

program test_quaternion
    use jwst_l2_constants, only: dp, eps_test_default
    use jwst_l2_quaternion
    implicit none

    integer :: n_pass, n_fail
    real(dp) :: a(3), b(3), c1(3), c2(3), maxdiff
    real(dp) :: q_id(4), rmat(3, 3), rmat_rt(3, 3), eye(3, 3)
    real(dp) :: omega(3), omat(4, 4), omat_t(4, 4), omat_anti(4, 4)
    real(dp) :: hat_a(3, 3)
    integer  :: i

    n_pass = 0
    n_fail = 0

    ! --- T1: hat_map(a) . b == cross(a, b) ---
    a = [1.0_dp, 2.0_dp, 3.0_dp]
    b = [-0.5_dp, 0.7_dp, 1.2_dp]
    hat_a = hat_map(a)
    c1 = matmul(hat_a, b)
    ! cross(a, b) by direct formula
    c2(1) = a(2)*b(3) - a(3)*b(2)
    c2(2) = a(3)*b(1) - a(1)*b(3)
    c2(3) = a(1)*b(2) - a(2)*b(1)
    maxdiff = maxval(abs(c1 - c2))
    call report("T1 hat(a) . b == cross(a,b)", maxdiff < eps_test_default, maxdiff)

    ! --- T2: q_to_matrix(identity) == I_3 ---
    q_id = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    rmat = q_to_matrix(q_id)
    eye = 0.0_dp
    do i = 1, 3
        eye(i, i) = 1.0_dp
    end do
    maxdiff = maxval(abs(rmat - eye))
    call report("T2 q_to_matrix(identity) == I_3", maxdiff < eps_test_default, maxdiff)

    ! --- T3: q_to_matrix orthogonality (rotation by 60 deg about (1,1,1)/sqrt(3)) ---
    block
        real(dp) :: angle, axis(3), nrm, qw, sin_h, q_unit(4)
        angle = 60.0_dp * (4.0_dp * atan(1.0_dp)) / 180.0_dp  ! 60 deg in radians
        axis = [1.0_dp, 1.0_dp, 1.0_dp]
        nrm = sqrt(sum(axis*axis))
        axis = axis / nrm
        qw = cos(angle / 2.0_dp)
        sin_h = sin(angle / 2.0_dp)
        q_unit = [qw, sin_h * axis(1), sin_h * axis(2), sin_h * axis(3)]
        rmat = q_to_matrix(q_unit)
        rmat_rt = matmul(rmat, transpose(rmat))
        maxdiff = maxval(abs(rmat_rt - eye))
    end block
    call report("T3 q_to_matrix orthogonality R . R^T = I_3", maxdiff < eps_test_default, maxdiff)

    ! --- T4: q_kinematics_matrix antisymmetry ---
    omega = [0.123_dp, -0.456_dp, 0.789_dp]
    omat = q_kinematics_matrix(omega)
    omat_t = transpose(omat)
    omat_anti = omat + omat_t  ! should be zero
    maxdiff = maxval(abs(omat_anti))
    call report("T4 q_kinematics_matrix antisymmetry Omega^T = -Omega", &
                maxdiff < eps_test_default, maxdiff)

    ! --- summary ---
    write (*, '(a)') ""
    write (*, '(a, i0, a, i0)') "test_quaternion: ", n_pass, " PASS / ", n_fail, " FAIL"
    if (n_fail > 0) then
        stop 1
    end if

contains

    subroutine report(label, passed, residual)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: passed
        real(dp),         intent(in) :: residual
        character(len=8)             :: status_str

        if (passed) then
            status_str = "PASS"
            n_pass = n_pass + 1
        else
            status_str = "FAIL"
            n_fail = n_fail + 1
        end if
        write (*, '(a, a, a, es12.4, a, a)') "  ", trim(label), "  (residual=", residual, ")  ", trim(status_str)
    end subroutine report

end program test_quaternion
