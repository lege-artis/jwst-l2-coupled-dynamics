"""
testbeds — analytical testbed oracles + numerical validation for Phase 2.

Phase 2 deliverables:
  lagrange_top_libration.py   — Ch 06 §6.4.3 oracle + numerical validation
                                 + Δt/2 and Δt/4 convergence-order test
                                 (Sitting 2 deliverable)
  lyapunov_three_scenario.py  — three-scenario Lyapunov calibration:
                                 torque-free integrable + GG-stationary +
                                 GG-Dzhanibekov-IC
                                 (Sitting 3 deliverable)

References per Phase 2 brief §2.5 + §2.6:
  Podolský, J. (2018) Teoretická mechanika, §8.3 (Lagrange top)
  Landau, L.D., Lifshitz, E.M. (1976) Mechanics, §37 (asymmetric top)
  arxiv:2003.13539 — Geometric origin of the tennis racket effect
"""
from .lagrange_top_libration import libration_period_analytical   # noqa: F401
