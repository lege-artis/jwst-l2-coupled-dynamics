"""
conftest.py — pytest configuration for jwst-l2-coupled-dynamics.

Registers custom markers:
  slow — tests that take > 30s (e.g. 10^6-step long-time energy tests).
         Run with: pytest -m slow
         Skip with: pytest -m "not slow"  (default)
"""
import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow (>30s); deselect with -m 'not slow'",
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("-m", default=""):
        # By default, skip slow tests unless explicitly requested
        skip_slow = pytest.mark.skip(reason="slow test: use -m slow to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)
