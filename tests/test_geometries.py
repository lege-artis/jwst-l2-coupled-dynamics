"""
test_geometries.py — unit tests for geometries module.

Tests the closed-form elementary inertia tensors, the parallel-axis
theorem, and the composite-body assembly + principal-axis decomposition.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from geometries import (
    disk_inertia, cylinder_inertia, thin_rod_inertia, solid_cone_inertia,
    parallel_axis, RigidComponent, CompositeBody,
    make_jwst_like, make_probe_like, make_parametric_cone_probe,
)


# ============================== elementary tensors ===============================

class TestDiskInertia:
    def test_symmetry_axis_z_default(self):
        I = disk_inertia(mass=10.0, radius=1.0, axis_dir="z")
        # diagonal in body frame; symmetry axis = z -> I_zz = mR^2/2, I_xx = I_yy = mR^2/4
        assert I[0, 0] == pytest.approx(2.5, rel=1e-12)
        assert I[1, 1] == pytest.approx(2.5, rel=1e-12)
        assert I[2, 2] == pytest.approx(5.0, rel=1e-12)
        # off-diagonals must be exactly zero
        assert np.allclose(I - np.diag(np.diag(I)), 0.0)

    def test_axis_dir_x_permutes(self):
        Iz = disk_inertia(10.0, 1.0, "z")
        Ix = disk_inertia(10.0, 1.0, "x")
        # I along symmetry axis should be the largest moment; under axis_dir = "x"
        # it should now appear in I[0,0]
        assert Ix[0, 0] == pytest.approx(Iz[2, 2], rel=1e-12)

    def test_symmetric_and_positive_definite(self):
        I = disk_inertia(5.0, 2.0, "y")
        assert np.allclose(I, I.T)
        eigvals = np.linalg.eigvalsh(I)
        assert np.all(eigvals > 0)


class TestCylinderInertia:
    def test_axial_and_perpendicular_moments(self):
        # cylinder along z: I_axial = m*R^2/2, I_perp = m*(3*R^2 + L^2)/12
        m, R, L = 12.0, 1.0, 3.0
        I = cylinder_inertia(m, R, L, "z")
        assert I[2, 2] == pytest.approx(m * R * R / 2.0, rel=1e-12)
        expected_perp = m * (3.0 * R * R + L * L) / 12.0
        assert I[0, 0] == pytest.approx(expected_perp, rel=1e-12)
        assert I[1, 1] == pytest.approx(expected_perp, rel=1e-12)


class TestSolidConeInertia:
    def test_axial_moment_closed_form(self):
        # solid cone about axis through centroid: I_axial = 3*m*R^2/10
        m, R, L = 100.0, 0.4, 5.0
        I = solid_cone_inertia(m, R, L, "z")
        assert I[2, 2] == pytest.approx(3.0 * m * R * R / 10.0, rel=1e-12)

    def test_perpendicular_moment_closed_form(self):
        # I_perp at centroid = 3*m*(4*L^2 + R^2)/80
        m, R, L = 100.0, 0.4, 5.0
        I = solid_cone_inertia(m, R, L, "z")
        expected_perp = 3.0 * m * (4.0 * L * L + R * R) / 80.0
        assert I[0, 0] == pytest.approx(expected_perp, rel=1e-12)
        assert I[1, 1] == pytest.approx(expected_perp, rel=1e-12)

    def test_known_anisotropy_at_specific_length(self):
        # docstring claims: L=5, R=0.4, m=100 -> I_axial = 4.8, I_perp = 375.6
        I = solid_cone_inertia(100.0, 0.4, 5.0, "z")
        assert I[2, 2] == pytest.approx(4.8, rel=1e-12)
        assert I[0, 0] == pytest.approx(375.6, rel=1e-12)
        # anisotropy I_perp / I_axial = 78.25
        assert I[0, 0] / I[2, 2] == pytest.approx(78.25, rel=1e-12)


class TestThinRodInertia:
    def test_zero_moment_along_rod_axis(self):
        I = thin_rod_inertia(mass=10.0, length=4.0, axis_dir="z")
        assert I[2, 2] == 0.0  # zero moment about the rod's own axis
        # perpendicular moments: m*L^2/12
        assert I[0, 0] == pytest.approx(10.0 * 16.0 / 12.0, rel=1e-12)


# ============================== parallel axis theorem ===============================

class TestParallelAxis:
    def test_zero_offset_is_identity(self):
        I0 = np.diag([1.0, 2.0, 3.0])
        I_shifted = parallel_axis(I0, mass=5.0, r_offset=np.zeros(3))
        assert np.allclose(I_shifted, I0)

    def test_offset_adds_m_r_squared_to_perpendicular_axes(self):
        # offset along +x by d: I_yy and I_zz gain m*d^2; I_xx unchanged
        m, d = 5.0, 2.0
        I0 = np.diag([1.0, 1.0, 1.0])
        I_shifted = parallel_axis(I0, m, np.array([d, 0.0, 0.0]))
        assert I_shifted[0, 0] == pytest.approx(1.0, rel=1e-12)
        assert I_shifted[1, 1] == pytest.approx(1.0 + m * d * d, rel=1e-12)
        assert I_shifted[2, 2] == pytest.approx(1.0 + m * d * d, rel=1e-12)
        # off-diagonals stay zero for axis-aligned offset
        assert np.allclose(I_shifted - np.diag(np.diag(I_shifted)), 0.0)

    def test_diagonal_offset_introduces_off_diagonals(self):
        # offset by (d, d, 0) -> off-diagonal I_xy = -m*d*d
        m, d = 1.0, 1.0
        I0 = np.zeros((3, 3))
        I_shifted = parallel_axis(I0, m, np.array([d, d, 0.0]))
        # |r|^2 = 2*d^2; rr^T diagonal is (d^2, d^2, 0); off-diag of rr^T at (0,1) = d*d
        # I_new[0,1] = -m * d * d
        assert I_shifted[0, 1] == pytest.approx(-m * d * d, rel=1e-12)
        assert I_shifted[1, 0] == pytest.approx(-m * d * d, rel=1e-12)  # symmetry
        assert I_shifted[0, 0] == pytest.approx(m * d * d, rel=1e-12)
        assert I_shifted[1, 1] == pytest.approx(m * d * d, rel=1e-12)
        assert I_shifted[2, 2] == pytest.approx(2.0 * m * d * d, rel=1e-12)


# ============================== composite body ===============================

class TestCompositeBody:
    def test_total_mass_is_sum(self):
        body = CompositeBody(name="t", components=[
            RigidComponent("a", 10.0, np.eye(3), np.array([0.0, 0.0, 0.0])),
            RigidComponent("b", 20.0, np.eye(3), np.array([1.0, 0.0, 0.0])),
        ])
        assert body.total_mass == pytest.approx(30.0, rel=1e-12)

    def test_center_of_mass_at_balance_point(self):
        body = CompositeBody(name="t", components=[
            RigidComponent("a", 10.0, np.eye(3), np.array([0.0, 0.0, 0.0])),
            RigidComponent("b", 30.0, np.eye(3), np.array([4.0, 0.0, 0.0])),
        ])
        # CoM = (10*0 + 30*4)/40 = 3 along x
        assert body.center_of_mass[0] == pytest.approx(3.0, rel=1e-12)
        assert body.center_of_mass[1] == pytest.approx(0.0, abs=1e-15)
        assert body.center_of_mass[2] == pytest.approx(0.0, abs=1e-15)

    def test_inertia_at_com_is_symmetric_pd(self):
        body = make_jwst_like()
        I = body.inertia_at_com()
        assert np.allclose(I, I.T)
        eigvals = np.linalg.eigvalsh(I)
        assert np.all(eigvals > 0)

    def test_principal_axes_are_orthonormal(self):
        body = make_jwst_like()
        eigvals, eigvecs = body.principal_axes()
        # eigvecs columns are orthonormal: V^T V = I
        product = eigvecs.T @ eigvecs
        assert np.allclose(product, np.eye(3), atol=1e-12)
        # eigvals sorted ascending
        assert eigvals[0] <= eigvals[1] <= eigvals[2]

    def test_jwst_like_anisotropy_matches_docstring(self):
        body = make_jwst_like()
        eigvals = np.sort(np.linalg.eigvalsh(body.inertia_at_com()))
        ratio = eigvals[-1] / eigvals[0]
        # JWST-like docstring: principal moments ~ (15384, 23322, 23322), ratio ~ 1.516
        assert ratio == pytest.approx(1.516, rel=1e-3)

    def test_probe_like_anisotropy_matches_docstring(self):
        body = make_probe_like()
        eigvals = np.sort(np.linalg.eigvalsh(body.inertia_at_com()))
        ratio = eigvals[-1] / eigvals[0]
        # probe-like docstring: principal moments ~ (5, 29.06, 29.06), ratio ~ 5.81
        assert ratio == pytest.approx(5.812, rel=1e-3)


class TestParametricConeProbe:
    @pytest.mark.parametrize("length, expected_axial, expected_perp", [
        (1.0, 4.8, 15.6),
        (5.0, 4.8, 375.6),
        (15.0, 4.8, 3375.6),
    ])
    def test_anisotropy_scales_with_length(self, length, expected_axial, expected_perp):
        # Closed-form check using base_radius = 0.4 default
        body = make_parametric_cone_probe(length=length)
        I = body.inertia_at_com()
        # axial moment is the smallest (prolate cone with R=0.4 << L)
        eigvals = np.sort(np.linalg.eigvalsh(I))
        assert eigvals[0] == pytest.approx(expected_axial, rel=1e-6)
        assert eigvals[-1] == pytest.approx(expected_perp, rel=1e-6)

    def test_centroid_at_body_origin(self):
        body = make_parametric_cone_probe(length=5.0)
        assert np.allclose(body.center_of_mass, np.zeros(3), atol=1e-15)
