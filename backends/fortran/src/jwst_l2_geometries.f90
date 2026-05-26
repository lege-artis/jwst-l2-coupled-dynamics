! jwst_l2_geometries
! ------------------
! Composite-body inertia primitives for the JWST-L2 Lege-Artis Fortran
! reference. Direct port of geometries.py - formulas locked per OQ-FORT-1
! resolution (Pollard 1966 "Mathematical Introduction to Celestial Mechanics"
! Table 2.1; equivalent in Goldstein 2002 problem 5.16 for the cone case).
!
! Canonical references:
!   - Chapter 02 Sec.2.4 (inertia tensor + parallel-axis theorem derivations)
!   - geometries.py (the Python reference being mirrored)
!   - Hughes 2004 Ch 4 (engineering treatment of composite rigid-body inertia)
!
! Conventions:
!   - All inertia tensors returned in the component's own centroid frame.
!     Callers wanting the tensor at a different reference point apply
!     parallel_axis() with the appropriate offset.
!   - axis_dir is an integer 1, 2, or 3 (corresponding to body x, y, z axis)
!     not the Python string "x"/"y"/"z" - string-handling in Fortran is
!     verbose and the integer index is unambiguous.
!   - All quantities in SI units (kg, m, s).
!
! Composite-body construction:
!   - inertia_at_com_assembly(components, n_comp) computes the composite
!     inertia tensor at the assembly CoM by summing each component's
!     parallel-axis-shifted contribution. This matches the Python
!     CompositeBody.inertia_at_com() method element-for-element.
!
! Author: Opus prep pass 2026-05-24
! License: Apache 2.0 (inherited from parent project)

module jwst_l2_geometries
    use jwst_l2_constants, only: dp
    implicit none
    private

    ! Public elementary inertia primitives
    public :: disk_inertia
    public :: cylinder_inertia
    public :: thin_rod_inertia
    public :: solid_cone_inertia
    public :: parallel_axis

    ! Public composite-body helpers
    public :: rigid_component_t
    public :: composite_inertia_at_com
    public :: composite_center_of_mass
    public :: composite_total_mass

    ! Composite-body component type (matches Python RigidComponent dataclass)
    type :: rigid_component_t
        character(len=32) :: name           ! human-readable label
        real(dp)          :: mass            ! kg
        real(dp)          :: inertia_local(3, 3)  ! in component's centroid frame
        real(dp)          :: position(3)     ! centroid offset from assembly origin
    end type rigid_component_t

    ! Axis-direction enumeration constants (clearer than bare integers at call sites)
    integer, parameter, public :: axis_x = 1
    integer, parameter, public :: axis_y = 2
    integer, parameter, public :: axis_z = 3

contains

    !---------------------------------------------------------------------
    ! disk_inertia(mass, radius, axis_dir) -> I(3,3)
    !
    ! Thin disk: I about symmetry axis = m*R^2/2; I about a diameter = m*R^2/4.
    ! axis_dir = axis_z means the disk plane is perpendicular to z-axis.
    !
    ! Port of geometries.py disk_inertia (lines 25-36).
    !---------------------------------------------------------------------
    pure function disk_inertia(mass, radius, axis_dir) result(inertia)
        real(dp), intent(in)  :: mass, radius
        integer,  intent(in)  :: axis_dir
        real(dp)              :: inertia(3, 3)
        real(dp)              :: i_principal, i_diametral

        i_principal = 0.5_dp * mass * radius * radius   ! along symmetry axis
        i_diametral = 0.25_dp * mass * radius * radius  ! about a diameter

        inertia = 0.0_dp
        select case (axis_dir)
        case (axis_x)
            inertia(1, 1) = i_principal
            inertia(2, 2) = i_diametral
            inertia(3, 3) = i_diametral
        case (axis_y)
            inertia(1, 1) = i_diametral
            inertia(2, 2) = i_principal
            inertia(3, 3) = i_diametral
        case default  ! axis_z (default)
            inertia(1, 1) = i_diametral
            inertia(2, 2) = i_diametral
            inertia(3, 3) = i_principal
        end select
    end function disk_inertia

    !---------------------------------------------------------------------
    ! cylinder_inertia(mass, radius, length, axis_dir) -> I(3,3)
    !
    ! Solid cylinder: I_axial = m*R^2/2;
    ! I_perp = m*(3*R^2 + L^2)/12 about a diameter through centroid.
    !
    ! Port of geometries.py cylinder_inertia (lines 65-76).
    !---------------------------------------------------------------------
    pure function cylinder_inertia(mass, radius, length, axis_dir) result(inertia)
        real(dp), intent(in)  :: mass, radius, length
        integer,  intent(in)  :: axis_dir
        real(dp)              :: inertia(3, 3)
        real(dp)              :: i_axial, i_perp

        i_axial = 0.5_dp * mass * radius * radius
        i_perp  = mass * (3.0_dp * radius * radius + length * length) / 12.0_dp

        inertia = 0.0_dp
        select case (axis_dir)
        case (axis_x)
            inertia(1, 1) = i_axial
            inertia(2, 2) = i_perp
            inertia(3, 3) = i_perp
        case (axis_y)
            inertia(1, 1) = i_perp
            inertia(2, 2) = i_axial
            inertia(3, 3) = i_perp
        case default  ! axis_z
            inertia(1, 1) = i_perp
            inertia(2, 2) = i_perp
            inertia(3, 3) = i_axial
        end select
    end function cylinder_inertia

    !---------------------------------------------------------------------
    ! thin_rod_inertia(mass, length, axis_dir) -> I(3,3)
    !
    ! Thin rod along the given axis: I along its axis = 0;
    ! I about a perpendicular axis through centroid = m*L^2/12.
    !
    ! Port of geometries.py thin_rod_inertia (lines 79-88).
    !---------------------------------------------------------------------
    pure function thin_rod_inertia(mass, length, axis_dir) result(inertia)
        real(dp), intent(in)  :: mass, length
        integer,  intent(in)  :: axis_dir
        real(dp)              :: inertia(3, 3)
        real(dp)              :: i_perp

        i_perp = mass * length * length / 12.0_dp

        inertia = 0.0_dp
        select case (axis_dir)
        case (axis_x)
            inertia(1, 1) = 0.0_dp
            inertia(2, 2) = i_perp
            inertia(3, 3) = i_perp
        case (axis_y)
            inertia(1, 1) = i_perp
            inertia(2, 2) = 0.0_dp
            inertia(3, 3) = i_perp
        case default  ! axis_z
            inertia(1, 1) = i_perp
            inertia(2, 2) = i_perp
            inertia(3, 3) = 0.0_dp
        end select
    end function thin_rod_inertia

    !---------------------------------------------------------------------
    ! solid_cone_inertia(mass, base_radius, height, axis_dir) -> I(3,3)
    !
    ! Solid right circular cone, vertex at origin, axis along +axis_dir,
    ! base at distance `height` from vertex.
    !
    ! Closed-form principal moments about the cone's CENTROID
    ! (Pollard 1966 'Mathematical Introduction to Celestial Mechanics'
    ! Table 2.1; equivalent in Goldstein 2002 problem 5.16):
    !
    !   centroid is at (3/4)*height from vertex along the axis
    !   I_axial  = (3/10) * mass * base_radius^2
    !   I_perp   = (3/80) * mass * (base_radius^2 + 4 * height^2)
    !
    ! Returned tensor is at the cone's centroid in the body frame.
    ! NOTE per geometries.py docstring: callers wanting the tensor at the
    ! assembly origin must apply parallel_axis() with the appropriate offset.
    !
    ! Port of geometries.py solid_cone_inertia (lines 39-62).
    !---------------------------------------------------------------------
    pure function solid_cone_inertia(mass, base_radius, height, axis_dir) result(inertia)
        real(dp), intent(in)  :: mass, base_radius, height
        integer,  intent(in)  :: axis_dir
        real(dp)              :: inertia(3, 3)
        real(dp)              :: i_axial, i_perp

        i_axial = (3.0_dp / 10.0_dp) * mass * base_radius * base_radius
        i_perp  = (3.0_dp / 80.0_dp) * mass * (4.0_dp * height * height + base_radius * base_radius)

        inertia = 0.0_dp
        select case (axis_dir)
        case (axis_x)
            inertia(1, 1) = i_axial
            inertia(2, 2) = i_perp
            inertia(3, 3) = i_perp
        case (axis_y)
            inertia(1, 1) = i_perp
            inertia(2, 2) = i_axial
            inertia(3, 3) = i_perp
        case default  ! axis_z
            inertia(1, 1) = i_perp
            inertia(2, 2) = i_perp
            inertia(3, 3) = i_axial
        end select
    end function solid_cone_inertia

    !---------------------------------------------------------------------
    ! parallel_axis(i_centroid, mass, r_offset) -> I_new(3,3)
    !
    ! Parallel-axis theorem: shift inertia tensor from a component's centroid
    ! to a new reference point at offset r from that centroid.
    !   I_new = I_centroid + mass * (|r|^2 * I_3 - r * r^T)
    !
    ! Canonical reference: Chapter 02 Sec.2.4 (boxed formula).
    ! Port of geometries.py parallel_axis (lines 91-99) - outer product
    ! r * r^T computed element-wise here, matching numpy.outer.
    !---------------------------------------------------------------------
    pure function parallel_axis(i_centroid, mass, r_offset) result(i_new)
        real(dp), intent(in)  :: i_centroid(3, 3)
        real(dp), intent(in)  :: mass
        real(dp), intent(in)  :: r_offset(3)
        real(dp)              :: i_new(3, 3)
        real(dp)              :: rr(3, 3), r_squared
        integer               :: i, j

        r_squared = r_offset(1)**2 + r_offset(2)**2 + r_offset(3)**2

        ! Outer product rr = r_offset * r_offset^T
        do j = 1, 3
            do i = 1, 3
                rr(i, j) = r_offset(i) * r_offset(j)
            end do
        end do

        ! Parallel-axis shift: I_centroid + m * (r^2 * I_3 - rr)
        i_new = i_centroid
        do i = 1, 3
            i_new(i, i) = i_new(i, i) + mass * r_squared
        end do
        i_new = i_new - mass * rr
    end function parallel_axis

    !---------------------------------------------------------------------
    ! composite_total_mass(components, n_comp) -> total_mass
    !
    ! Sum the masses of all components in the composite body.
    !---------------------------------------------------------------------
    pure function composite_total_mass(components, n_comp) result(total_mass)
        type(rigid_component_t), intent(in)  :: components(:)
        integer,                 intent(in)  :: n_comp
        real(dp)                             :: total_mass
        integer                              :: i

        total_mass = 0.0_dp
        do i = 1, n_comp
            total_mass = total_mass + components(i)%mass
        end do
    end function composite_total_mass

    !---------------------------------------------------------------------
    ! composite_center_of_mass(components, n_comp) -> com(3)
    !
    ! Compute the CoM position of the composite body relative to the assembly
    ! origin. CoM = (sum m_i * r_i) / (sum m_i).
    !---------------------------------------------------------------------
    pure function composite_center_of_mass(components, n_comp) result(com)
        type(rigid_component_t), intent(in)  :: components(:)
        integer,                 intent(in)  :: n_comp
        real(dp)                             :: com(3)
        real(dp)                             :: total_mass
        integer                              :: i

        com = 0.0_dp
        total_mass = composite_total_mass(components, n_comp)
        do i = 1, n_comp
            com = com + components(i)%mass * components(i)%position
        end do
        com = com / total_mass
    end function composite_center_of_mass

    !---------------------------------------------------------------------
    ! composite_inertia_at_com(components, n_comp) -> I_total(3,3)
    !
    ! Composite inertia tensor evaluated at the assembly CoM in body frame.
    ! Sum of each component's parallel-axis-translated tensor.
    !
    ! Port of geometries.py CompositeBody.inertia_at_com (lines 132-140).
    !---------------------------------------------------------------------
    pure function composite_inertia_at_com(components, n_comp) result(i_total)
        type(rigid_component_t), intent(in)  :: components(:)
        integer,                 intent(in)  :: n_comp
        real(dp)                             :: i_total(3, 3)
        real(dp)                             :: com(3), r_to_com(3)
        integer                              :: i

        com = composite_center_of_mass(components, n_comp)
        i_total = 0.0_dp
        do i = 1, n_comp
            r_to_com = components(i)%position - com
            i_total = i_total + parallel_axis(components(i)%inertia_local, &
                                              components(i)%mass, &
                                              r_to_com)
        end do
    end function composite_inertia_at_com

end module jwst_l2_geometries
