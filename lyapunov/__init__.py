"""
lyapunov — Lyapunov spectrum computation for SE(3) × SE(3) dynamics.

Phase 2 deliverables (Sitting 3):
  cocycle.py   — Euclidean (R^14 quaternion-component) variational-equation
                 propagator  (OQ-PHASE2-2 DEFAULT = Euclidean)
  spectrum.py  — Benettin-Galgani-Giorgilli-Strelcyn 1980 algorithm with
                 pilot-then-production QR-cadence autotuning
  reference.py — known-Lyapunov-spectrum oracles for calibration testbeds

References:
  Benettin, G., Galgani, L., Giorgilli, A., Strelcyn, J.-M. (1980)
    "Lyapunov characteristic exponents for smooth dynamical systems and
    for Hamiltonian systems; a method for computing all of them",
    Meccanica 15, 9–30.
  Skokos, C. (2010) "The Lyapunov Characteristic Exponents and Their
    Computation", Lect. Notes Phys. 790, Springer.
"""
from .spectrum import lyapunov_spectrum   # noqa: F401
