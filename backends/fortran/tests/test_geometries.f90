! test_geometries
! ---------------
! D1 deliverable: composite-body inertia primitives (Ch 02 Sec.2.4).
!
! Tests:
!   T1. disk_inertia closed form: I_axial = m*R^2/2, I_perp = m*R^2/4
!   T2. cylinder_inertia closed form: I_axial = m*R^2/2, I_perp = m*(3R^2+L^2)/12
!   T3. thin_rod_inertia closed form: I_axial = 0, I_perp = m*L^2/12
!   T4. solid_cone_inertia closed form: I_axial = (3/10)m*r^2, I_perp = (3/80)m*(4h^2+r^2)
!   T5. parallel_axis on a point mass: I_new = m*(|r|^2*I_3 - r*r^T)
!   T6. JWST-like composite reference: I_at_com matches Python first-cut element-wise
!   T7. Probe-like composite reference: I_at_com matches Python first-cut element-wise
!
! Reference values for T6/T7 derived from running geometries.py directly
! (2026-05-24, Opus prep pass). NOTE: these values differ from the
! canonical-tier docs claim of I_jwst = diag(2540, 2540, 3870) - see
! OQ-FORT-6 in PREP-STATUS-2026-05-24.md. The C2 oracle is the Python
! first-cut, so this test asserts against the actual Python output.
!
! Exit 0 on all PASS, exit 1 on any FAIL.
!
! Author: Opus prep pass 2026-05-24

program test_geometries
    use jwst_l2_constants, only: dp, eps_test_default
    use jwst_l2_geometries
    implicit none

    integer :: n_pass, n_fail
    real(dp) :: inertia(3, 3), expected(3, 3), maxdiff
    real(dp) :: r_off(3), i_centroid(3, 3), i_new(3, 3)
    real(dp) :: mass, radius, length, base_radius, height, r_sq
    type(rigid_component_t) :: jwst_components(4), probe_components(2)
    real(dp) :: i_composite(3, 3)

    n_pass = 0
    n_fail = 0

    ! --- T1: disk_inertia (mass=10, R=2, axis_z) ---
    mass = 10.0_dp
    radius = 2.0_dp
    inertia = disk_inertia(mass, radius, axis_z)
    expected = 0.0_dp
    expected(1, 1) = 0.25_dp * mass * radius * radius  ! 10
    expected(2, 2) = 0.25_dp * mass * radius * radius  ! 10
    expected(3, 3) = 0.5_dp * mass * radius * radius   ! 20
    maxdiff = maxval(abs(inertia - expected))
    call report("T1 disk_inertia analytical", maxdiff < eps_test_default, maxdiff)

    ! --- T2: cylinder_inertia (mass=5, R=1, L=3, axis_z) ---
    mass = 5.0_dp; radius = 1.0_dp; length = 3.0_dp
    inertia = cylinder_inertia(mass, radius, length, axis_z)
    expected = 0.0_dp
    expected(1, 1) = mass * (3.0_dp * radius * radius + length * length) / 12.0_dp  ! 5
    expected(2, 2) = mass * (3.0_dp * radius * radius + length * length) / 12.0_dp  ! 5
    expected(3, 3) = 0.5_dp * mass * radius * radius                                  ! 2.5
    maxdiff = maxval(abs(inertia - expected))
    call report("T2 cylinder_inertia analytical", maxdiff < eps_test_default, maxdiff)

    ! --- T3: thin_rod_inertia (mass=2, L=6, axis_z) ---
    mass = 2.0_dp; length = 6.0_dp
    inertia = thin_rod_inertia(mass, length, axis_z)
    expected = 0.0_dp
    expected(1, 1) = mass * length * length / 12.0_dp  ! 6
    expected(2, 2) = mass * length * length / 12.0_dp  ! 6
    expected(3, 3) = 0.0_dp
    maxdiff = maxval(abs(inertia - expected))
    call report("T3 thin_rod_inertia analytical", maxdiff < eps_test_default, maxdiff)

    ! --- T4: solid_cone_inertia (mass=12, r=2, h=4, axis_z) ---
    !   I_axial = (3/10)*12*4 = 14.4
    !   I_perp  = (3/80)*12*(64+4) = (3/80)*12*68 = 30.6
    mass = 12.0_dp; base_radius = 2.0_dp; height = 4.0_dp
    inertia = solid_cone_inertia(mass, base_radius, height, axis_z)
    expected = 0.0_dp
    expected(1, 1) = (3.0_dp / 80.0_dp) * mass * (4.0_dp * height**2 + base_radius**2)
    expected(2, 2) = (3.0_dp / 80.0_dp) * mass * (4.0_dp * height**2 + base_radius**2)
    expected(3, 3) = (3.0_dp / 10.0_dp) * mass * base_radius**2
    maxdiff = maxval(abs(inertia - expected))
    call report("T4 solid_cone_inertia analytical (Pollard 1966)", maxdiff < eps_test_default, maxdiff)

    ! --- T5: parallel_axis on point mass (I_centroid=0, m=10, r=(1,2,3)) ---
    mass = 10.0_dp
    r_off = [1.0_dp, 2.0_dp, 3.0_dp]
    i_centroid = 0.0_dp
    i_new = parallel_axis(i_centroid, mass, r_off)
    ! Expected: m*(|r|^2 * I_3 - r*r^T) = 10*(14*I_3 - rr)
    r_sq = sum(r_off * r_off)  ! = 14
    expected(1, 1) = mass * (r_sq - r_off(1) * r_off(1))  ! 10*(14-1) = 130
    expected(1, 2) = -mass * r_off(1) * r_off(2)          ! -10*2 = -20
    expected(1, 3) = -mass * r_off(1) * r_off(3)          ! -10*3 = -30
    expected(2, 1) = -mass * r_off(2) * r_off(1)          ! -20
    expected(2, 2) = mass * (r_sq - r_off(2) * r_off(2))  ! 10*(14-4) = 100
    expected(2, 3) = -mass * r_off(2) * r_off(3)          ! -60
    expected(3, 1) = -mass * r_off(3) * r_off(1)          ! -30
    expected(3, 2) = -mass * r_off(3) * r_off(2)          ! -60
    expected(3, 3) = mass * (r_sq - r_off(3) * r_off(3))  ! 10*(14-9) = 50
    maxdiff = maxval(abs(i_new - expected))
    call report("T5 parallel_axis on point mass", maxdiff < eps_test_default, maxdiff)

    ! --- T6: JWST-like composite (4 components, must match Python first-cut) ---
    ! Reference values from running geometries.py make_jwst_like().inertia_at_com():
    !   I = diag(2.332262967687075e+04, 2.332262967687075e+04, 1.538437500000000e+04)
    !   total_mass = 2450 kg, CoM = (0, 0, 0.32653061)
    ! NOTE: differs from canonical-tier docs claim (2540, 2540, 3870) - see OQ-FORT-6.
    !
    ! Components per geometries.py make_jwst_like():
    !   sunshield: disk, m=200, R=10, z=-3
    !   bus:       cylinder, m=1500, R=1.5, L=3, z=-1
    !   boom:      thin_rod, m=50, L=4, z=+2
    !   primary:   disk, m=700, R=3.25, z=+4
    jwst_components(1)%name          = "sunshield"
    jwst_components(1)%mass          = 200.0_dp
    jwst_components(1)%inertia_local = disk_inertia(200.0_dp, 10.0_dp, axis_z)
    jwst_components(1)%position      = [0.0_dp, 0.0_dp, -3.0_dp]

    jwst_components(2)%name          = "bus"
    jwst_components(2)%mass          = 1500.0_dp
    jwst_components(2)%inertia_local = cylinder_inertia(1500.0_dp, 1.5_dp, 3.0_dp, axis_z)
    jwst_components(2)%position      = [0.0_dp, 0.0_dp, -1.0_dp]

    jwst_components(3)%name          = "boom"
    jwst_components(3)%mass          = 50.0_dp
    jwst_components(3)%inertia_local = thin_rod_inertia(50.0_dp, 4.0_dp, axis_z)
    jwst_components(3)%position      = [0.0_dp, 0.0_dp, 2.0_dp]

    jwst_components(4)%name          = "primary_mirror"
    jwst_components(4)%mass          = 700.0_dp
    jwst_components(4)%inertia_local = disk_inertia(700.0_dp, 3.25_dp, axis_z)
    jwst_components(4)%position      = [0.0_dp, 0.0_dp, 4.0_dp]

    i_composite = composite_inertia_at_com(jwst_components, 4)
    expected = 0.0_dp
    expected(1, 1) = 2.332262967687075e+04_dp
    expected(2, 2) = 2.332262967687075e+04_dp
    expected(3, 3) = 1.538437500000000e+04_dp
    maxdiff = maxval(abs(i_composite - expected))
    call report("T6 JWST-like composite vs Python first-cut (C2 element-wise)", &
                maxdiff < 1.0e-9_dp, maxdiff)

    ! --- T7: Probe-like composite (2 components, must match Python first-cut) ---
    ! Reference values: I = diag(29.06, 29.06, 5.0), total_mass = 100 kg
    probe_components(1)%name          = "probe_body"
    probe_components(1)%mass          = 80.0_dp
    probe_components(1)%inertia_local = cylinder_inertia(80.0_dp, 0.25_dp, 1.5_dp, axis_z)
    probe_components(1)%position      = [0.0_dp, 0.0_dp, 0.0_dp]

    probe_components(2)%name          = "dish_antenna"
    probe_components(2)%mass          = 20.0_dp
    probe_components(2)%inertia_local = disk_inertia(20.0_dp, 0.5_dp, axis_z)
    probe_components(2)%position      = [0.0_dp, 0.0_dp, 0.85_dp]

    i_composite = composite_inertia_at_com(probe_components, 2)
    expected = 0.0_dp
    expected(1, 1) = 29.06_dp
    expected(2, 2) = 29.06_dp
    expected(3, 3) = 5.00_dp
    maxdiff = maxval(abs(i_composite - expected))
    call report("T7 Probe-like composite vs Python first-cut (C2 element-wise)", &
                maxdiff < 1.0e-12_dp, maxdiff)

    ! --- T8: Oblate reference body (synthetic single cylinder, OQ-FORT-6 Path C resolution) ---
    ! Mirrors geometries.py make_oblate_reference_body() default parameters:
    !   mass = 2450 kg, i_axial = 3870 kg.m^2, i_perp = 2540 kg.m^2
    ! Solve for cylinder dimensions from target moments:
    !   R^2 = 2 * i_axial / mass        =  2 * 3870 / 2450  = 3.1591836...
    !   L^2 = 12 * i_perp / mass - 3*R^2 = 12 * 2540 / 2450 - 3 * 3.1591836... = 2.9632653...
    ! Single cylinder centred at body origin -> no parallel-axis shift.
    ! Expected composite I_at_com = diag(2540, 2540, 3870) to machine precision by construction.
    block
        type(rigid_component_t) :: oblate_components(1)
        real(dp) :: target_mass, target_i_axial, target_i_perp
        real(dp) :: r_cyl, l_cyl

        target_mass    = 2450.0_dp
        target_i_axial = 3870.0_dp
        target_i_perp  = 2540.0_dp

        r_cyl = sqrt(2.0_dp * target_i_axial / target_mass)
        l_cyl = sqrt(12.0_dp * target_i_perp / target_mass - 3.0_dp * (r_cyl * r_cyl))

        oblate_components(1)%name          = "oblate_cylinder"
        oblate_components(1)%mass          = target_mass
        oblate_components(1)%inertia_local = cylinder_inertia(target_mass, r_cyl, l_cyl, axis_z)
        oblate_components(1)%position      = [0.0_dp, 0.0_dp, 0.0_dp]

        i_composite = composite_inertia_at_com(oblate_components, 1)
        expected = 0.0_dp
        expected(1, 1) = target_i_perp
        expected(2, 2) = target_i_perp
        expected(3, 3) = target_i_axial
        maxdiff = maxval(abs(i_composite - expected))
        call report("T8 Oblate reference body diag(2540,2540,3870) by construction", &
                    maxdiff < 1.0e-10_dp, maxdiff)
    end block

    ! --- summary ---
    write (*, '(a)') ""
    write (*, '(a, i0, a, i0)') "test_geometries: ", n_pass, " PASS / ", n_fail, " FAIL"
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

end program test_geometries
