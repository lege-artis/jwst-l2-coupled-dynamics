"""
geometries.py — composite-body inertia tensors for the JWST-L2 first-cut example.

Bodies modeled as cones + disks + cylinders + thin rods, decomposed at their
local centroid then translated to the assembly center-of-mass via the
parallel-axis theorem.

References:
- Goldstein, Poole, Safko (2002) "Classical Mechanics" 3rd ed §5 (inertia tensor)
- Hughes (2004) "Spacecraft Attitude Dynamics" §3 (composite rigid body)

Conventions:
- All quantities in SI units (m, kg, s)
- Body frame: spacecraft Z-axis is the symmetry/boom axis (sunshield -> mirror)
- Inertia tensor in body frame; transform to inertial frame via R . I . R^T
- This is a SIMPLIFIED MODEL — illustrative geometry, not exact JWST values
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


# --- elementary inertia tensors (about the component's own centroid) ---

def disk_inertia(mass: float, radius: float, axis_dir: str = "z") -> np.ndarray:
    """Thin disk: I about its symmetry axis = m*R^2 / 2;
       I about a diameter = m*R^2 / 4.
       axis_dir = 'z' means disk plane is perpendicular to z-axis."""
    Ip = 0.5 * mass * radius * radius   # principal (symmetry axis)
    Id = 0.25 * mass * radius * radius  # diametral
    I = np.diag([Id, Id, Ip])           # default: symmetry along z
    if axis_dir == "x":
        I = np.diag([Ip, Id, Id])
    elif axis_dir == "y":
        I = np.diag([Id, Ip, Id])
    return I


def solid_cone_inertia(mass: float, base_radius: float, height: float,
                        axis_dir: str = "z") -> np.ndarray:
    """Solid right circular cone, vertex at origin, axis along +axis_dir,
       base at distance `height` from vertex.

       Closed-form principal moments about the cone's centroid (Pollard 1966
       'Mathematical Introduction to Celestial Mechanics' table 2.1;
       equivalent in Goldstein 2002 problem 5.16):

         centroid is at height * 3/4 from the vertex along the axis
         I_axial   = (3/10) * mass * base_radius^2
         I_perp_at_centroid = (3/80) * mass * (base_radius^2 + 4 * height^2)

       Returns the inertia tensor at the centroid in the body frame.
       NOTE: callers wanting the tensor at the assembly origin must apply
       parallel_axis() with the appropriate centroid offset."""
    Ia = (3.0 / 10.0) * mass * base_radius * base_radius
    Ip = (3.0 / 80.0) * mass * (4.0 * height * height + base_radius * base_radius)
    I = np.diag([Ip, Ip, Ia])
    if axis_dir == "x":
        I = np.diag([Ia, Ip, Ip])
    elif axis_dir == "y":
        I = np.diag([Ip, Ia, Ip])
    return I


def cylinder_inertia(mass: float, radius: float, length: float,
                      axis_dir: str = "z") -> np.ndarray:
    """Solid cylinder: I_axial = m*R^2/2;
       I_perp = m*(3*R^2 + L^2)/12 about a diameter through centroid."""
    Ia = 0.5 * mass * radius * radius
    Ip = mass * (3.0 * radius * radius + length * length) / 12.0
    I = np.diag([Ip, Ip, Ia])
    if axis_dir == "x":
        I = np.diag([Ia, Ip, Ip])
    elif axis_dir == "y":
        I = np.diag([Ip, Ia, Ip])
    return I


def thin_rod_inertia(mass: float, length: float, axis_dir: str = "z") -> np.ndarray:
    """Thin rod along the given axis: I along its axis = 0;
       I about a perpendicular axis through centroid = m*L^2/12."""
    Iperp = mass * length * length / 12.0
    I = np.diag([Iperp, Iperp, 0.0])  # default axis = z
    if axis_dir == "x":
        I = np.diag([0.0, Iperp, Iperp])
    elif axis_dir == "y":
        I = np.diag([Iperp, 0.0, Iperp])
    return I


def parallel_axis(I_centroid: np.ndarray, mass: float,
                  r_offset: np.ndarray) -> np.ndarray:
    """Parallel-axis theorem: shift inertia tensor from a component's centroid
       to a new reference point at offset r from that centroid.
       I_new = I_centroid + m * (|r|^2 * I3 - r r^T)
       where I3 is the 3x3 identity and r r^T is the outer product."""
    rr = np.outer(r_offset, r_offset)
    r2 = float(np.dot(r_offset, r_offset))
    return I_centroid + mass * (r2 * np.eye(3) - rr)


# --- composite body specification ---

@dataclass
class RigidComponent:
    """One geometric component of a composite rigid body."""
    name: str
    mass: float                    # kg
    inertia_local: np.ndarray      # 3x3 in component's own centroid frame
    position: np.ndarray           # 3-vector: centroid offset from assembly origin


@dataclass
class CompositeBody:
    """Composite rigid body: list of components + computed assembly properties."""
    name: str
    components: list[RigidComponent] = field(default_factory=list)

    @property
    def total_mass(self) -> float:
        return sum(c.mass for c in self.components)

    @property
    def center_of_mass(self) -> np.ndarray:
        """CoM position relative to assembly origin."""
        m_total = self.total_mass
        com = np.zeros(3)
        for c in self.components:
            com += c.mass * c.position
        return com / m_total

    def inertia_at_com(self) -> np.ndarray:
        """Composite inertia tensor evaluated at the assembly CoM in body frame.
           Sum of each component's parallel-axis-translated tensor."""
        com = self.center_of_mass
        I_total = np.zeros((3, 3))
        for c in self.components:
            r_to_com = c.position - com
            I_total += parallel_axis(c.inertia_local, c.mass, r_to_com)
        return I_total

    def principal_axes(self) -> tuple[np.ndarray, np.ndarray]:
        """Diagonalise inertia tensor at CoM.
           Returns (eigenvalues, eigenvectors) where eigenvectors[:, k] is
           the k-th principal axis in body frame, eigenvalues sorted ascending."""
        I = self.inertia_at_com()
        eigvals, eigvecs = np.linalg.eigh(I)
        return eigvals, eigvecs


# --- specific bodies: JWST-like + probe-like ---

def make_jwst_like() -> CompositeBody:
    """Simplified JWST-style spacecraft. Geometry along body z-axis:
       sunshield (-z end), spacecraft bus (mid), boom + secondary mirror,
       primary mirror (+z end).

       NOT actual JWST values — illustrative + bounded for a teaching example.
       Real JWST: ~6500 kg total launch mass. We use ~2350 kg simplified."""
    components = []

    # Sunshield: thin disk, plane perpendicular to z-axis, at z = -3 m
    sunshield = RigidComponent(
        name="sunshield",
        mass=200.0,                                   # kg
        inertia_local=disk_inertia(200.0, 10.0, "z"), # R = 10 m
        position=np.array([0.0, 0.0, -3.0]),
    )
    components.append(sunshield)

    # Spacecraft bus: short cylinder along z, centered at z = -1 m
    bus = RigidComponent(
        name="bus",
        mass=1500.0,
        inertia_local=cylinder_inertia(1500.0, 1.5, 3.0, "z"),
        position=np.array([0.0, 0.0, -1.0]),
    )
    components.append(bus)

    # Boom: thin rod along z, length 4 m, centered at z = +2 m
    boom = RigidComponent(
        name="boom",
        mass=50.0,
        inertia_local=thin_rod_inertia(50.0, 4.0, "z"),
        position=np.array([0.0, 0.0, 2.0]),
    )
    components.append(boom)

    # Primary mirror: thin disk, plane perpendicular to z-axis, at z = +4 m
    primary = RigidComponent(
        name="primary_mirror",
        mass=700.0,
        inertia_local=disk_inertia(700.0, 3.25, "z"),
        position=np.array([0.0, 0.0, 4.0]),
    )
    components.append(primary)

    return CompositeBody(name="jwst_like", components=components)


def make_probe_like() -> CompositeBody:
    """Simplified inspection-probe geometry: small cylindrical body with one
       dish antenna. Total mass ~100 kg (CubeSat-class + slightly larger)."""
    components = []

    # Probe main body: cylinder along z, centered at origin
    body = RigidComponent(
        name="probe_body",
        mass=80.0,
        inertia_local=cylinder_inertia(80.0, 0.25, 1.5, "z"),
        position=np.array([0.0, 0.0, 0.0]),
    )
    components.append(body)

    # Dish antenna: thin disk at +z end of cylinder
    dish = RigidComponent(
        name="dish_antenna",
        mass=20.0,
        inertia_local=disk_inertia(20.0, 0.5, "z"),
        position=np.array([0.0, 0.0, 0.85]),
    )
    components.append(dish)

    return CompositeBody(name="probe_like", components=components)


def make_oblate_reference_body(mass: float = 2450.0,
                                i_axial: float = 3870.0,
                                i_perp: float = 2540.0) -> CompositeBody:
    """Synthetic oblate axisymmetric reference body with specified principal moments.

    Constructs a SINGLE uniform cylinder of given mass, dimensioned so that the
    resulting inertia tensor at its centroid matches diag(i_perp, i_perp, i_axial)
    in the body frame (with the symmetry axis along body z). The cylinder is
    centred at the body-frame origin, so no parallel-axis shift is required.

    Solves for R and L from the cylinder primitive formulas:
        I_axial = m R^2 / 2                     ->  R^2 = 2 * i_axial / m
        I_perp  = m (3 R^2 + L^2) / 12          ->  L^2 = 12 * i_perp / m - 3 R^2

    Constraint: requires i_axial > 0 and 4 * i_perp > i_axial (otherwise L^2 < 0
    and the construction is infeasible). The constraint 4 * i_perp > i_axial bounds
    the achievable anisotropy: for I_axial / I_perp > 2 (more anisotropic than a pure
    disk) a single cylinder cannot represent the body and a composite (disk + ring)
    would be required.

    Default values (2450 kg, 3870, 2540 kg*m^2) match the canonical-tier docs
    Chapter 02 Sec.2.4 + Sec.3.7.2 reference for the OBLATE gravity-gradient
    instability worked example. This is the "synthetic oblate reference body" that
    pairs with make_jwst_like() (which models a sunshield-on-boom composite that
    happens to be prolate due to component spread along z) -- see OQ-FORT-6
    investigation 2026-05-24 for the doc-vs-code drift resolution.

    Resolved anisotropy at default values:
        I_axial / I_perp = 3870 / 2540 = 1.524 (oblate)
        max/min eigenvalue ratio = 1.524 (identical because oblate axisymmetric)

    Parameters:
        mass     : total body mass, kg (default 2450 to match the canonical-tier
                   docs claim of make_jwst_like total mass)
        i_axial  : moment of inertia about the symmetry axis (z), kg*m^2 (default 3870)
        i_perp   : moment of inertia about an axis perpendicular to z, kg*m^2
                   (default 2540)
    """
    if i_axial <= 0:
        raise ValueError("i_axial must be positive")
    if 4.0 * i_perp <= i_axial:
        raise ValueError(f"requires 4 * i_perp > i_axial (got i_perp={i_perp}, "
                         f"i_axial={i_axial}); single-cylinder construction infeasible")
    R_squared = 2.0 * i_axial / mass
    L_squared = 12.0 * i_perp / mass - 3.0 * R_squared
    R = float(np.sqrt(R_squared))
    L = float(np.sqrt(L_squared))
    cyl = RigidComponent(
        name="oblate_cylinder",
        mass=mass,
        inertia_local=cylinder_inertia(mass, R, L, "z"),
        position=np.array([0.0, 0.0, 0.0]),
    )
    return CompositeBody(name="oblate_reference_body", components=[cyl])


def make_parametric_cone_probe(length: float,
                                mass: float = 100.0,
                                base_radius: float = 0.4
                                ) -> CompositeBody:
    """Parametric probe shaped as a single solid cone with radial symmetry.

    Geometry: cone of total height `length` along its body z-axis, with the
    centroid placed at the body-frame origin. Base radius is held FIXED
    across the testcase sweep so that varying `length` produces a meaningful
    inertia-anisotropy sweep (NOT just a uniform scaling).

    Inertia anisotropy as a function of length (with base_radius fixed at 0.4 m):
      length=1.0:  I_axial = 4.80,  I_perp =  15.6,  I_perp/I_axial =   3.25
      length=5.0:  I_axial = 4.80,  I_perp = 375.6,  I_perp/I_axial =  78.25
      length=15.0: I_axial = 4.80,  I_perp = 3375.6, I_perp/I_axial = 703.25

    The cone is PROLATE (perpendicular moments larger than axial) at every
    length >= base_radius. As length grows the anisotropy diverges; the cone
    becomes effectively a long thin pencil. This is the "space station
    approaching" parametric shape: simple, radially symmetric, single-
    parameter, with a non-trivial coupling response to gravity-gradient
    torques because the increasing I_perp amplifies the gyroscopic term in
    Euler's equations.

    Parameters:
      length:       total cone height along its symmetry axis (m)
      mass:         total mass (kg); default 100 kg
      base_radius:  cone base radius (m); default 0.4 m (FIXED across the sweep)"""
    components = [
        RigidComponent(
            name=f"cone_L{length}",
            mass=mass,
            inertia_local=solid_cone_inertia(mass, base_radius, length, "z"),
            position=np.array([0.0, 0.0, 0.0]),  # cone centroid at body origin
        )
    ]
    return CompositeBody(name=f"parametric_cone_probe_L{length}", components=components)


# --- quick self-test when invoked directly ---

if __name__ == "__main__":
    jwst = make_jwst_like()
    probe = make_probe_like()
    # also instantiate the parametric cone-probe at three lengths to display
    # the inertia-anisotropy sweep
    cone_short = make_parametric_cone_probe(length=1.0)
    cone_med = make_parametric_cone_probe(length=5.0)
    cone_long = make_parametric_cone_probe(length=15.0)
    for body in (jwst, probe, cone_short, cone_med, cone_long):
        I = body.inertia_at_com()
        eigvals, eigvecs = body.principal_axes()
        print(f"--- {body.name} ---")
        print(f"  total mass:    {body.total_mass:.2f} kg")
        print(f"  CoM (body fr): {body.center_of_mass}")
        print(f"  Inertia tensor (kg.m^2):")
        for row in I:
            print(f"    [{row[0]:14.6e} {row[1]:14.6e} {row[2]:14.6e}]")
        print(f"  Principal moments: {eigvals}")
        print(f"  Inertia ratio max/min: {eigvals[-1] / eigvals[0]:.3f}")
        print()
