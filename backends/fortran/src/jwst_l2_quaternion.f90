! jwst_l2_quaternion
! ------------------
! Quaternion algebra + hat map (so(3) Lie algebra primitives) for the
! JWST-L2 Lege-Artis Fortran reference. Direct port of dynamics.py
! quaternion utilities, conventions matched line-for-line so the C2
! cross-implementation bit-identity gate passes.
!
! Canonical references:
!   - Chapter 01 Sec.1.4 (hat map and dq/dt = 0.5 * Omega(omega) * q derivation)
!   - Chapter 02 Sec.2.5.1 (the rotational-transport equation that consumes hat_map)
!   - Chapter 02 Sec.2.7.1 (left-invariance discussion of body-frame omega)
!   - dynamics.py lines 41-68 (the Python reference being mirrored)
!
! Conventions per SONNET-BRIEF Sec.2 (locked, NON-NEGOTIABLE):
!   - Quaternion ordering: q = (qw, qx, qy, qz) with qw scalar
!   - Rotation matrix R(q) maps body-frame to inertial-frame:
!         v_inertial = R . v_body
!   - Hat map sign: hat(v) such that hat(v) . b == cross(v, b) for any b
!     Equivalent component form: hat(v)_{ij} = -epsilon_{ijk} v_k
!   - Quaternion kinematics: dq/dt = 0.5 * Omega(omega) * q where omega is
!     in BODY frame (per Euler's equations canonical form)
!
! Author: Opus prep pass 2026-05-24
! License: Apache 2.0 (inherited from parent project)

module jwst_l2_quaternion
    use jwst_l2_constants, only: dp
    implicit none
    private

    public :: q_normalise
    public :: q_to_matrix
    public :: q_kinematics_matrix
    public :: hat_map

contains

    !---------------------------------------------------------------------
    ! q_normalise(q) -> q_unit
    !
    ! Renormalise a quaternion to unit length (counters integration drift).
    ! Port of dynamics.py q_normalise (line 42).
    !---------------------------------------------------------------------
    pure function q_normalise(q) result(q_unit)
        real(dp), intent(in)  :: q(4)
        real(dp)              :: q_unit(4)
        real(dp)              :: nrm
        nrm = sqrt(q(1)**2 + q(2)**2 + q(3)**2 + q(4)**2)
        q_unit = q / nrm
    end function q_normalise

    !---------------------------------------------------------------------
    ! q_to_matrix(q) -> R(3,3)
    !
    ! Rotation matrix from unit quaternion (qw, qx, qy, qz).
    ! R maps body-frame vectors to inertial-frame vectors: v_inertial = R . v_body.
    ! Port of dynamics.py q_to_matrix (lines 47-56) - element-wise identical.
    !
    ! Indexing note: Fortran is 1-based, Python is 0-based. Component layout:
    !   Fortran q(1) = Python qw, q(2) = qx, q(3) = qy, q(4) = qz
    !   Fortran R(1,1) = Python R[0,0], etc.
    !---------------------------------------------------------------------
    pure function q_to_matrix(q) result(rmat)
        real(dp), intent(in)  :: q(4)
        real(dp)              :: rmat(3, 3)
        real(dp)              :: qw, qx, qy, qz

        qw = q(1); qx = q(2); qy = q(3); qz = q(4)

        rmat(1, 1) = 1.0_dp - 2.0_dp * (qy*qy + qz*qz)
        rmat(1, 2) = 2.0_dp * (qx*qy - qz*qw)
        rmat(1, 3) = 2.0_dp * (qx*qz + qy*qw)
        rmat(2, 1) = 2.0_dp * (qx*qy + qz*qw)
        rmat(2, 2) = 1.0_dp - 2.0_dp * (qx*qx + qz*qz)
        rmat(2, 3) = 2.0_dp * (qy*qz - qx*qw)
        rmat(3, 1) = 2.0_dp * (qx*qz - qy*qw)
        rmat(3, 2) = 2.0_dp * (qy*qz + qx*qw)
        rmat(3, 3) = 1.0_dp - 2.0_dp * (qx*qx + qy*qy)
    end function q_to_matrix

    !---------------------------------------------------------------------
    ! q_kinematics_matrix(omega) -> Omega(4,4)
    !
    ! Quaternion rate matrix Omega(omega) such that dq/dt = 0.5 * Omega . q.
    ! omega is in the BODY frame (per dynamics.py line 60-61 convention).
    ! Port of dynamics.py q_kinematics_matrix (lines 59-68) - element-wise identical.
    !---------------------------------------------------------------------
    pure function q_kinematics_matrix(omega) result(omat)
        real(dp), intent(in)  :: omega(3)
        real(dp)              :: omat(4, 4)
        real(dp)              :: wx, wy, wz

        wx = omega(1); wy = omega(2); wz = omega(3)

        ! Row 1
        omat(1, 1) =  0.0_dp; omat(1, 2) = -wx;   omat(1, 3) = -wy;   omat(1, 4) = -wz
        ! Row 2
        omat(2, 1) =  wx;     omat(2, 2) =  0.0_dp; omat(2, 3) =  wz;   omat(2, 4) = -wy
        ! Row 3
        omat(3, 1) =  wy;     omat(3, 2) = -wz;   omat(3, 3) =  0.0_dp; omat(3, 4) =  wx
        ! Row 4
        omat(4, 1) =  wz;     omat(4, 2) =  wy;   omat(4, 3) = -wx;   omat(4, 4) =  0.0_dp
    end function q_kinematics_matrix

    !---------------------------------------------------------------------
    ! hat_map(v) -> hat(3,3)
    !
    ! Hat-map isomorphism so(3) ~= R^3. Returns the antisymmetric 3x3 matrix
    ! such that hat(v) . b = cross(v, b) for any 3-vector b.
    !
    ! Component form: hat(v)_{ij} = -epsilon_{ijk} v_k.
    ! Matrix form:    hat(v) = [[ 0  -v_z  v_y]
    !                            [ v_z  0  -v_x]
    !                            [-v_y v_x   0 ]]
    !
    ! Canonical reference: Chapter 01 Sec.1.4 + Chapter 02 Sec.2.5.1 + Sec.2.7.1.
    ! Python reference: this function is the standalone equivalent of the
    ! hat-map structure implicit inside dynamics.py q_kinematics_matrix and
    ! gravity_gradient_torque_body_frame.
    !
    ! Verification test (in test_quaternion.f90): for arbitrary a, b in R^3,
    ! matmul(hat_map(a), b) must equal cross(a, b) to machine precision.
    !---------------------------------------------------------------------
    pure function hat_map(v) result(hat)
        real(dp), intent(in)  :: v(3)
        real(dp)              :: hat(3, 3)
        real(dp)              :: vx, vy, vz

        vx = v(1); vy = v(2); vz = v(3)

        hat(1, 1) =  0.0_dp; hat(1, 2) = -vz;   hat(1, 3) =  vy
        hat(2, 1) =  vz;     hat(2, 2) =  0.0_dp; hat(2, 3) = -vx
        hat(3, 1) = -vy;     hat(3, 2) =  vx;   hat(3, 3) =  0.0_dp
    end function hat_map

end module jwst_l2_quaternion
