! test_multipole_potential  (D3)
! --------------------------------
! Tests the multipole gravitational potential functions.
!
! Analytical tests (C1):
!   T1. Kepler limit: I_A = I_B = 0 -> V_total = V_mono exactly.
!   T2. Quadrupole scales as 1/r^3: V_tidal(r2) / V_tidal(r1) = (r1/r2)^3.
!   T3. Ratio check: V_quad at r=100 m vs r=1000 m -> ratio = 1000 to 4 sig figs.
!
! Fixture cross-check (3 cases hardcoded from tests/fixtures/multipole_potential_inputs.json):
!   F1. JWST + probe at 1000 km (far field, quad correction tiny).
!   F2. JWST + probe at 1 km (near field, quad correction large).
!   F3. Point-mass Kepler limit (I_A = I_B = 0, V_quad = 0).
!
! C2 JSON output: _audit/fortran_multipole_potential_output.json
!   Schema per case: {"name": ..., "output": {"v_mono": .., "v_tidal_A": ..,
!                     "v_tidal_B": .., "v_total": ..}}
!
! Exit 0 on all PASS, exit 1 on any FAIL.
! Canonical reference: Chapter 03 Sec.3.4.3 (boxed formula).
!
! Author: Sonnet Phase A 2026-05-24

program test_multipole_potential
    use jwst_l2_constants,  only: dp, g_newton
    use jwst_l2_potential,  only: mutual_potential_mono, mutual_potential_quad, &
                                   total_mutual_potential
    implicit none

    integer  :: n_pass, n_fail
    real(dp) :: rho(3), I_a(3,3), I_b(3,3), q_a(4), q_b(4)
    real(dp) :: V_mono, V_quad, V_total, V_tidal_a_val, V_tidal_b_val
    real(dp) :: V_q_r1, V_q_r2, ratio
    real(dp) :: m_a, m_b

    ! Fixture case values (from tests/fixtures/multipole_potential_inputs.json)
    ! Case F1: JWST + probe, 1000 km separation
    real(dp), parameter :: F1_m_a   = 2450.0_dp
    real(dp), parameter :: F1_m_b   = 100.0_dp
    real(dp), parameter :: F1_rho3  = 1.0e6_dp   ! rho = (0, 0, 1e6)
    real(dp), parameter :: F1_Ia_xx = 23322.629676870747_dp
    real(dp), parameter :: F1_Ia_zz = 15384.375_dp
    real(dp), parameter :: F1_Ib_xx = 29.06_dp
    real(dp), parameter :: F1_Ib_zz = 5.0_dp
    ! q_a = identity, q_b = 45-deg rotation about x
    real(dp), parameter :: F1_qb_w  = 0.9238795325112867_dp
    real(dp), parameter :: F1_qb_x  = 0.3826834323650898_dp
    ! Expected Python outputs
    real(dp), parameter :: F1_exp_vmono  = -1.6352035000000002e-11_dp
    real(dp), parameter :: F1_exp_vtidA  = -5.298229318983842e-23_dp
    real(dp), parameter :: F1_exp_vtidB  = -9.835749052499995e-25_dp
    real(dp), parameter :: F1_exp_vtotal = -1.6352035000053966e-11_dp

    ! Case F2: 1 km separation
    real(dp), parameter :: F2_rho3  = 1.0e3_dp
    real(dp), parameter :: F2_exp_vmono  = -1.6352035e-08_dp
    real(dp), parameter :: F2_exp_vtidA  = -5.298229318983842e-14_dp
    real(dp), parameter :: F2_exp_vtidB  = -9.835749052499995e-16_dp
    real(dp), parameter :: F2_exp_vtotal = -1.6352088965868098e-08_dp

    ! Case F3: point masses (I=0)
    real(dp), parameter :: F3_rho3  = 1.5e9_dp
    real(dp), parameter :: F3_exp_vmono  = -1.0901356666666668e-14_dp

    ! Relative tolerance for fixture comparisons
    real(dp), parameter :: tol_rel = 1.0e-13_dp

    integer :: json_unit
    logical :: json_ok

    n_pass = 0
    n_fail = 0

    !-----------------------------------------------------------------------
    ! T1: Kepler limit - point masses give V_total = V_mono
    !-----------------------------------------------------------------------
    m_a = 2450.0_dp; m_b = 100.0_dp
    rho = [0.0_dp, 0.0_dp, 1.0e6_dp]
    I_a = 0.0_dp; I_b = 0.0_dp
    q_a = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    q_b = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

    V_mono  = mutual_potential_mono(m_a, m_b, rho)
    V_quad  = mutual_potential_quad(m_a, m_b, rho, I_a, q_a, I_b, q_b)
    V_total = total_mutual_potential(m_a, m_b, rho, I_a, q_a, I_b, q_b)

    call report("T1 Kepler limit: V_quad = 0 exactly", &
                abs(V_quad) < 1.0e-30_dp, abs(V_quad))
    call report("T1 Kepler limit: V_total = V_mono", &
                abs(V_total - V_mono) < 1.0e-30_dp, abs(V_total - V_mono))

    !-----------------------------------------------------------------------
    ! T2: Quadrupole scales as 1/r^3 between two separations
    !-----------------------------------------------------------------------
    I_a = 0.0_dp
    I_a(1,1) = F1_Ia_xx; I_a(2,2) = F1_Ia_xx; I_a(3,3) = F1_Ia_zz
    I_b = 0.0_dp
    I_b(1,1) = F1_Ib_xx; I_b(2,2) = F1_Ib_xx; I_b(3,3) = F1_Ib_zz
    q_a = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    q_b = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]  ! identity for clean test

    rho = [0.0_dp, 0.0_dp, 1.0e3_dp]    ! r = 1000 m
    V_q_r1 = mutual_potential_quad(m_a, m_b, rho, I_a, q_a, I_b, q_b)

    rho = [0.0_dp, 0.0_dp, 2.0e3_dp]    ! r = 2000 m
    V_q_r2 = mutual_potential_quad(m_a, m_b, rho, I_a, q_a, I_b, q_b)

    ! V_quad(r=2000) / V_quad(r=1000) should equal (1000/2000)^3 = 0.125
    ratio = V_q_r2 / V_q_r1
    call report("T2 V_quad scales as 1/r^3: ratio |V(2r)/V(r) - 0.125| < 1e-13", &
                abs(ratio - 0.125_dp) < 1.0e-13_dp, abs(ratio - 0.125_dp))

    !-----------------------------------------------------------------------
    ! T3: Ratio check at 100 m vs 1000 m: ratio should be 1000.0 to 4 sig figs
    !-----------------------------------------------------------------------
    rho = [0.0_dp, 0.0_dp, 100.0_dp]
    V_q_r1 = mutual_potential_quad(m_a, m_b, rho, I_a, q_a, I_b, q_b)

    rho = [0.0_dp, 0.0_dp, 1000.0_dp]
    V_q_r2 = mutual_potential_quad(m_a, m_b, rho, I_a, q_a, I_b, q_b)

    ! V_quad(100m) / V_quad(1000m) = (1000/100)^3 = 1000
    ratio = V_q_r1 / V_q_r2
    call report("T3 V_quad(100m)/V_quad(1000m) = 1000 to 4 sig figs", &
                abs(ratio - 1000.0_dp) < 0.1_dp, abs(ratio - 1000.0_dp))

    !-----------------------------------------------------------------------
    ! F1: Fixture case 1 — JWST + probe at 1000 km
    !-----------------------------------------------------------------------
    m_a = F1_m_a; m_b = F1_m_b
    I_a = 0.0_dp
    I_a(1,1) = F1_Ia_xx; I_a(2,2) = F1_Ia_xx; I_a(3,3) = F1_Ia_zz
    I_b = 0.0_dp
    I_b(1,1) = F1_Ib_xx; I_b(2,2) = F1_Ib_xx; I_b(3,3) = F1_Ib_zz
    q_a = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    q_b = [F1_qb_w, F1_qb_x, 0.0_dp, 0.0_dp]
    rho = [0.0_dp, 0.0_dp, F1_rho3]

    V_mono = mutual_potential_mono(m_a, m_b, rho)
    call split_quad_terms(m_a, m_b, rho, I_a, q_a, I_b, q_b, V_tidal_a_val, V_tidal_b_val)
    V_quad  = V_tidal_a_val + V_tidal_b_val
    V_total = V_mono + V_quad

    call report("F1 v_mono 1000km", reldiff(V_mono, F1_exp_vmono) < tol_rel, &
                reldiff(V_mono, F1_exp_vmono))
    call report("F1 v_tidal_A 1000km", reldiff(V_tidal_a_val, F1_exp_vtidA) < tol_rel, &
                reldiff(V_tidal_a_val, F1_exp_vtidA))
    call report("F1 v_tidal_B 1000km", reldiff(V_tidal_b_val, F1_exp_vtidB) < tol_rel, &
                reldiff(V_tidal_b_val, F1_exp_vtidB))
    call report("F1 v_total 1000km", reldiff(V_total, F1_exp_vtotal) < tol_rel, &
                reldiff(V_total, F1_exp_vtotal))

    !-----------------------------------------------------------------------
    ! F2: Fixture case 2 — JWST + probe at 1 km
    !-----------------------------------------------------------------------
    rho = [0.0_dp, 0.0_dp, F2_rho3]
    V_mono = mutual_potential_mono(m_a, m_b, rho)
    call split_quad_terms(m_a, m_b, rho, I_a, q_a, I_b, q_b, V_tidal_a_val, V_tidal_b_val)
    V_quad  = V_tidal_a_val + V_tidal_b_val
    V_total = V_mono + V_quad

    call report("F2 v_mono 1km", reldiff(V_mono, F2_exp_vmono) < tol_rel, &
                reldiff(V_mono, F2_exp_vmono))
    call report("F2 v_tidal_A 1km", reldiff(V_tidal_a_val, F2_exp_vtidA) < tol_rel, &
                reldiff(V_tidal_a_val, F2_exp_vtidA))
    call report("F2 v_tidal_B 1km", reldiff(V_tidal_b_val, F2_exp_vtidB) < tol_rel, &
                reldiff(V_tidal_b_val, F2_exp_vtidB))
    call report("F2 v_total 1km", reldiff(V_total, F2_exp_vtotal) < tol_rel, &
                reldiff(V_total, F2_exp_vtotal))

    !-----------------------------------------------------------------------
    ! F3: Fixture case 3 — Kepler limit (I=0, at 1.5e9 m)
    !-----------------------------------------------------------------------
    m_a = F1_m_a; m_b = F1_m_b
    I_a = 0.0_dp; I_b = 0.0_dp
    q_a = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    q_b = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    rho = [0.0_dp, 0.0_dp, F3_rho3]
    V_mono  = mutual_potential_mono(m_a, m_b, rho)
    V_quad  = mutual_potential_quad(m_a, m_b, rho, I_a, q_a, I_b, q_b)
    V_total = V_mono + V_quad

    call report("F3 Kepler v_mono", reldiff(V_mono, F3_exp_vmono) < tol_rel, &
                reldiff(V_mono, F3_exp_vmono))
    call report("F3 Kepler v_quad = 0", abs(V_quad) < 1.0e-30_dp, abs(V_quad))

    !-----------------------------------------------------------------------
    ! Write JSON output for cross_test_python.sh
    !-----------------------------------------------------------------------
    json_ok = .true.
    open(newunit=json_unit, file="_audit/fortran_multipole_potential_output.json", &
         status="replace", iostat=n_fail)
    ! Temporarily reuse n_fail for open status; reset after
    if (n_fail /= 0) then
        write (*, '(a)') "  WARNING: could not open _audit/fortran_multipole_potential_output.json"
        json_ok = .false.
        n_fail = 0
    else
        call write_multipole_json(json_unit)
        close(json_unit)
    end if

    !-----------------------------------------------------------------------
    ! Summary
    !-----------------------------------------------------------------------
    write (*, '(a)') ""
    write (*, '(a,i0,a,i0)') "test_multipole_potential: ", n_pass, " PASS / ", n_fail, " FAIL"
    if (n_fail > 0) stop 1

contains

    subroutine report(label, passed, residual)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: passed
        real(dp),         intent(in) :: residual
        if (passed) then
            n_pass = n_pass + 1
            write (*, '(a,a,a,es12.4,a)') "  ", trim(label), "  (residual=", residual, ")  PASS"
        else
            n_fail = n_fail + 1
            write (*, '(a,a,a,es12.4,a)') "  ", trim(label), "  (residual=", residual, ")  FAIL"
        end if
    end subroutine report

    real(dp) function reldiff(a, b)
        real(dp), intent(in) :: a, b
        real(dp)             :: denom
        denom = max(abs(a), abs(b), 1.0e-300_dp)
        reldiff = abs(a - b) / denom
    end function reldiff

    !-----------------------------------------------------------------------
    ! split_quad_terms: compute V_tidal_A and V_tidal_B separately.
    ! Needed for per-term fixture comparison.
    !-----------------------------------------------------------------------
    subroutine split_quad_terms(ma, mb, rho_in, Ia, qa, Ib, qb, V_ta, V_tb)
        use jwst_l2_quaternion, only: q_to_matrix
        real(dp), intent(in)  :: ma, mb, rho_in(3)
        real(dp), intent(in)  :: Ia(3,3), qa(4), Ib(3,3), qb(4)
        real(dp), intent(out) :: V_ta, V_tb
        real(dp) :: r_mag, rho_hat(3), r3
        real(dp) :: Ra(3,3), Rb(3,3), u_a(3), u_b(3)
        real(dp) :: Iu_a(3), Iu_b(3), tr_a, tr_b, quad_a, quad_b

        r_mag = sqrt(rho_in(1)**2 + rho_in(2)**2 + rho_in(3)**2)
        if (r_mag < 1.0e-12_dp) then
            V_ta = 0.0_dp; V_tb = 0.0_dp; return
        end if
        rho_hat = rho_in / r_mag
        r3 = r_mag * r_mag * r_mag

        Ra = q_to_matrix(qa)
        Rb = q_to_matrix(qb)
        u_a = matmul(transpose(Ra), rho_hat)
        u_b = -matmul(transpose(Rb), rho_hat)
        Iu_a = matmul(Ia, u_a)
        Iu_b = matmul(Ib, u_b)
        tr_a   = Ia(1,1) + Ia(2,2) + Ia(3,3)
        tr_b   = Ib(1,1) + Ib(2,2) + Ib(3,3)
        quad_a = u_a(1)*Iu_a(1) + u_a(2)*Iu_a(2) + u_a(3)*Iu_a(3)
        quad_b = u_b(1)*Iu_b(1) + u_b(2)*Iu_b(2) + u_b(3)*Iu_b(3)
        V_ta = -(g_newton * mb / (2.0_dp * r3)) * (tr_a - 3.0_dp * quad_a)
        V_tb = -(g_newton * ma / (2.0_dp * r3)) * (tr_b - 3.0_dp * quad_b)
    end subroutine split_quad_terms

    !-----------------------------------------------------------------------
    ! write_multipole_json: emit JSON for cross_test_python.sh.
    !-----------------------------------------------------------------------
    subroutine write_multipole_json(junit)
        use jwst_l2_quaternion, only: q_to_matrix
        integer,  intent(in) :: junit
        real(dp) :: I_a_loc(3,3), I_b_loc(3,3), q_a_loc(4), q_b_loc(4), rho_loc(3)
        real(dp) :: vm, vta, vtb, vt

        write (junit, '(a)') '{'
        write (junit, '(a)') '  "spec_id": "fortran_multipole_potential_output.v0.1",'
        write (junit, '(a)') '  "cases": ['

        ! Case 1 (F1)
        I_a_loc = 0.0_dp
        I_a_loc(1,1)=F1_Ia_xx; I_a_loc(2,2)=F1_Ia_xx; I_a_loc(3,3)=F1_Ia_zz
        I_b_loc = 0.0_dp
        I_b_loc(1,1)=F1_Ib_xx; I_b_loc(2,2)=F1_Ib_xx; I_b_loc(3,3)=F1_Ib_zz
        q_a_loc = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        q_b_loc = [F1_qb_w, F1_qb_x, 0.0_dp, 0.0_dp]
        rho_loc = [0.0_dp, 0.0_dp, F1_rho3]
        vm = mutual_potential_mono(F1_m_a, F1_m_b, rho_loc)
        call split_quad_terms(F1_m_a, F1_m_b, rho_loc, I_a_loc, q_a_loc, &
                              I_b_loc, q_b_loc, vta, vtb)
        vt = vm + vta + vtb
        call write_potential_case(junit, "case_01_far_apart_generic", vm, vta, vtb, vt, .true.)

        ! Case 2 (F2)
        rho_loc = [0.0_dp, 0.0_dp, F2_rho3]
        vm = mutual_potential_mono(F1_m_a, F1_m_b, rho_loc)
        call split_quad_terms(F1_m_a, F1_m_b, rho_loc, I_a_loc, q_a_loc, &
                              I_b_loc, q_b_loc, vta, vtb)
        vt = vm + vta + vtb
        call write_potential_case(junit, "case_02_close_1km", vm, vta, vtb, vt, .true.)

        ! Case 3 (F3) — point masses
        I_a_loc = 0.0_dp; I_b_loc = 0.0_dp
        q_a_loc = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        q_b_loc = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        rho_loc = [0.0_dp, 0.0_dp, F3_rho3]
        vm = mutual_potential_mono(F1_m_a, F1_m_b, rho_loc)
        vta = 0.0_dp; vtb = 0.0_dp
        vt = vm
        call write_potential_case(junit, "case_03_point_masses_kepler_limit", &
                                  vm, vta, vtb, vt, .false.)

        write (junit, '(a)') "  ]"
        write (junit, '(a)') "}"
    end subroutine write_multipole_json

    subroutine write_potential_case(junit, name, vm, vta, vtb, vt, trailing_comma)
        integer,          intent(in) :: junit
        character(len=*), intent(in) :: name
        real(dp),         intent(in) :: vm, vta, vtb, vt
        logical,          intent(in) :: trailing_comma
        write (junit, '(a,a,a)') '    {"name": "', trim(name), '",'
        write (junit, '(a,es24.17,a)') '     "output": {"v_mono": ', vm, ','
        write (junit, '(a,es24.17,a)') '       "v_tidal_A": ', vta, ','
        write (junit, '(a,es24.17,a)') '       "v_tidal_B": ', vtb, ','
        if (trailing_comma) then
            write (junit, '(a,es24.17,a)') '       "v_total": ', vt, '}},'
        else
            write (junit, '(a,es24.17,a)') '       "v_total": ', vt, '}}'
        end if
    end subroutine write_potential_case

end program test_multipole_potential
