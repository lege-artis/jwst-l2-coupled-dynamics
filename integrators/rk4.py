"""
rk4.py — RK4 integrator re-export for backward compatibility.

Extracted from dynamics.py v0.1.x per Phase 2 §2.1 module-layout contract.
Callers can import from either dynamics or integrators.rk4 — both refer to
the same implementation object.

Reference: Press et al. (2007) Numerical Recipes 3rd ed §17.1, Eq. 17.1.3.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from dynamics import rk4_step, state_derivative  # noqa: F401
