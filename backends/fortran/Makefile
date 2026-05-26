# JWST-L2 Lege-Artis Fortran Reference — Phase 1 Makefile
# OS-aware build driver per KB-038 (Windows mkdir lacks -p); flags per KB-039
# (no -pedantic to avoid UTF-8 line-length cascade).
#
# Targets:
#   all         (default) build modules + run unit tests
#   build       compile modules only
#   test        run unit tests (depends on build)
#   cross-test  run Python-oracle cross-tests (depends on test + Python+numpy)
#   clean       remove build/ and test executables
#   help        list available targets
#
# Usage:
#   make             # equivalent to `make all`
#   make build       # compile only, no tests
#   make test        # run all unit tests
#   make cross-test  # run cross-tests against Python first-cut
#   make clean       # purge build artifacts

# === Compiler + flags ===
FC      = gfortran
FFLAGS  = -O2 -std=f2018 -Wall -Wextra -fimplicit-none -fcheck=bounds
LDFLAGS =

# === Directory layout ===
SRCDIR    = src
TESTDIR   = tests
TOOLDIR   = tools
BUILDDIR  = build
AUDITDIR  = _audit
FIXTURE_DIR = ../../tests/fixtures

# === OS-aware mkdir + rm + executable extension (KB-038) ===
ifeq ($(OS),Windows_NT)
    MKDIR_BUILD = if not exist $(BUILDDIR) mkdir $(BUILDDIR)
    MKDIR_AUDIT = if not exist $(AUDITDIR) mkdir $(AUDITDIR)
    RM      = del /Q
    EXE     = .exe
    PATHSEP = \\
else
    MKDIR_BUILD = mkdir -p $(BUILDDIR)
    MKDIR_AUDIT = mkdir -p $(AUDITDIR)
    RM      = rm -f
    EXE     =
    PATHSEP = /
endif

# === Module list — dependency-ordered (compile in this order) ===
# Sonnet authoring contract: each module name must match the .f90 file basename
# and the `module ... end module` identifier must match the file basename too.
# Ordering must reflect actual dependency chains (e.g. _quaternion is used by
# _torque, so _quaternion compiles first).
MODULES = jwst_l2_constants \
          jwst_l2_quaternion \
          jwst_l2_geometries \
          jwst_l2_potential \
          jwst_l2_torque \
          jwst_l2_jacobian \
          jwst_l2_diagnostics \
          jwst_l2_dynamics

MODULE_OBJS = $(MODULES:%=$(BUILDDIR)/%.o)
MODULE_MODS = $(MODULES:%=$(BUILDDIR)/%.mod)

# === Test list ===
# Each test_*.f90 in tests/ compiles to a standalone executable in build/.
# Test executables print PASS/FAIL to stdout; exit code 0 on PASS, non-zero on FAIL.
TESTS = test_quaternion \
        test_geometries \
        test_trace_identity \
        test_multipole_potential \
        test_gravity_gradient_torque \
        test_orientation_jacobian \
        test_rotational_transport \
        test_kepler_limit \
        test_symmetric_top \
        test_instability_oblate \
        test_instability_dzhanibekov \
        test_lyapunov_pair \
        test_anisotropy_sweep

TEST_BINS = $(TESTS:%=$(BUILDDIR)/%$(EXE))

# === Phony targets ===
.PHONY: all build test cross-test clean help debug

# === Default: build + test ===
all: build test

help:
	@echo "JWST-L2 Lege-Artis Fortran Reference - Phase 1 Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all         (default) build modules + run unit tests"
	@echo "  build       compile modules only"
	@echo "  test        run unit tests (depends on build)"
	@echo "  cross-test  run cross-tests against Python first-cut"
	@echo "  clean       remove build/ and test executables"
	@echo "  help        this message"
	@echo ""
	@echo "Configuration:"
	@echo "  FC      = $(FC)"
	@echo "  FFLAGS  = $(FFLAGS)"
	@echo "  Build:    $(BUILDDIR)"
	@echo "  Fixtures: $(FIXTURE_DIR)"

# === Build modules ===
build: $(MODULE_OBJS)

# === Compile each module ===
# Sonnet authoring contract: every module file is in src/ and produces both
# build/<name>.o and build/<name>.mod via the same compile step.
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	@$(MKDIR_BUILD)
	$(FC) $(FFLAGS) -J$(BUILDDIR) -c $< -o $@

# === Build all test executables ===
$(BUILDDIR)/%$(EXE): $(TESTDIR)/%.f90 $(MODULE_OBJS)
	@$(MKDIR_BUILD)
	$(FC) $(FFLAGS) -J$(BUILDDIR) -I$(BUILDDIR) $< $(MODULE_OBJS) -o $@ $(LDFLAGS)

# === Run all unit tests + summary ===
test: build $(TEST_BINS)
	@$(MKDIR_AUDIT)
	@echo "=========================================="
	@echo "  JWST-L2 Fortran Phase 1 — Unit test run"
	@echo "=========================================="
	@PASS_COUNT=0; FAIL_COUNT=0; \
	for t in $(TESTS); do \
		printf '%-40s ' "$$t"; \
		if $(BUILDDIR)$(PATHSEP)$$t$(EXE) > $(AUDITDIR)$(PATHSEP)$$t.log 2>&1; then \
			echo "PASS"; \
			PASS_COUNT=$$((PASS_COUNT + 1)); \
		else \
			echo "FAIL (see $(AUDITDIR)/$$t.log)"; \
			FAIL_COUNT=$$((FAIL_COUNT + 1)); \
		fi; \
	done; \
	echo "=========================================="; \
	echo "  Results: $$PASS_COUNT PASS / $$FAIL_COUNT FAIL"; \
	echo "=========================================="; \
	if [ $$FAIL_COUNT -gt 0 ]; then exit 1; fi

# === Cross-tests against Python first-cut ===
cross-test: test
	@$(MKDIR_AUDIT)
	bash $(TOOLDIR)/cross_test_python.sh

# === Clean ===
clean:
	@if [ -d "$(BUILDDIR)" ]; then $(RM) $(BUILDDIR)$(PATHSEP)*.o $(BUILDDIR)$(PATHSEP)*.mod $(BUILDDIR)$(PATHSEP)test_*$(EXE) 2>/dev/null || true; fi

# === Debug help ===
debug:
	@echo "MODULES     = $(MODULES)"
	@echo "MODULE_OBJS = $(MODULE_OBJS)"
	@echo "TESTS       = $(TESTS)"
	@echo "TEST_BINS   = $(TEST_BINS)"
	@echo "OS detected = $(if $(filter Windows_NT,$(OS)),Windows,POSIX)"
	@echo "MKDIR_BUILD = $(MKDIR_BUILD)"
	@echo "EXE         = $(EXE)"
