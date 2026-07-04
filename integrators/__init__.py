"""
integrators — structure-preserving integrators for two-body SE(3) × SE(3) dynamics.

Provided integrators (Phase 2):
  rk4.rk4_step                     — classical RK4; preserved from dynamics.py (v0.1.x)
  symplectic_se3_variational.se3_variational_step — KDK Störmer-Verlet on SE(3)
                                      (Lee-Leok-McClamroch 2007 variational framework)
  symplectic_se3_splitting.se3_splitting_step     — Yoshida 4th-order splitting
                                      (OQ-PHASE2-1 BOTH backstop; Sitting 2 deliverable)

All integrators share the same 26-component state-vector convention as dynamics.py.
"""
from .rk4 import rk4_step                                   # noqa: F401
from .symplectic_se3_variational import se3_variational_step  # noqa: F401
from .symplectic_se3_splitting import se3_splitting_step      # noqa: F401
