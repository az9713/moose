# ============================================================
# Case 84: Lossy TEM Cavity Resonator — Q Factor Calculation
#
# References:
#   Pozar, Microwave Engineering, 4th ed. (2012), Sec. 6.7–6.8
#   Collin, Foundations for Microwave Engineering, 2nd ed. (1992), Ch. 7
#   Kong, Electromagnetic Wave Theory, 2nd ed. (1990), Ch. 1
#   MIT 6.635 Advanced Electromagnetism, Spring 2003, Prof. Jin Au Kong
#   OCW: https://ocw.mit.edu/courses/6-635-advanced-electromagnetism-spring-2003/
#
# ============================================================
# PHYSICS — CAVITY RESONANCES IN A LOSSY DIELECTRIC
# ============================================================
#
# A 1D Fabry-Pérot-like resonator is formed by a slab of dielectric
# with COMPLEX permittivity ε_r = ε' − j ε'' enclosed between two
# Perfect Electric Conductor (PEC) walls at x = 0 and x = D.
#
# The convention ε_r = ε' − j ε'' corresponds to a time-harmonic
# convention exp(+jωt); the imaginary part ε'' > 0 represents
# absorption (lossy dielectric).  With the exp(−jωt) convention
# the sign would be reversed; the Kong textbook uses exp(+jωt)
# so we follow that convention here.
#
# Inside the cavity the electric field satisfies the 1D Helmholtz
# eigenvalue equation:
#
#   d²E/dx² + k₀² ε_r E = 0
#
# where k₀ = ω/c is the free-space wavenumber and ε_r is complex.
# This is a SELF-ADJOINT eigenvalue problem (for real ε_r) but
# becomes NON-SELF-ADJOINT with complex ε_r, yielding COMPLEX
# eigenvalues λ = (ω/c)² = k₀² that encode both the resonance
# frequency and the damping rate.
#
# ============================================================
# COUPLED REAL/IMAGINARY EIGENVALUE FORMULATION
# ============================================================
#
# Write E = E_real + j E_imag and λ = λ_real + j λ_imag.
# With ε_r = ε' − j ε'' (where ε', ε'' are REAL positive numbers),
# the complex Helmholtz equation splits into two coupled equations:
#
#   d²E_real/dx²  +  λ ε' E_real  +  λ ε'' E_imag  = 0   ... (1)
#   d²E_imag/dx²  +  λ ε' E_imag  −  λ ε'' E_real  = 0   ... (2)
#
# (After substituting E = E_r + j E_i and λ = λ_r + j λ_i and
#  separating real/imaginary parts, the coupling between equations
#  comes from the ε'' terms.)
#
# This is a GENERALISED EIGENVALUE PROBLEM  A·x = λ·B·x  where:
#   A = stiffness block matrix (from Diffusion kernels)
#   B = ε-weighted mass block matrix (from MatReaction + CoupledForce)
#   x = [E_real; E_imag]  (stacked real/imaginary DOF vectors)
#   λ = (ω/c)²            (complex eigenvalue)
#
# ============================================================
# WEAK FORM AND MOOSE KERNEL MAPPING
# ============================================================
#
# Multiply eq. (1) by test function φ, integrate over [0, D]:
#
#   ∫₀ᴰ (dE_real/dx)(dφ/dx) dx = λ · [-∫₀ᴰ ε' E_real φ dx
#                                        -∫₀ᴰ ε'' E_imag φ dx]
#
# The boundary terms vanish because φ = 0 at the PEC walls
# (homogeneous Dirichlet BC on E_real and E_imag at x = 0, D).
#
# Multiply eq. (2) by test function φ:
#
#   ∫₀ᴰ (dE_imag/dx)(dφ/dx) dx = λ · [-∫₀ᴰ ε' E_imag φ dx
#                                        +∫₀ᴰ ε'' E_real φ dx]
#
# MOOSE kernel contributions (sign conventions):
#
#   Diffusion(E_real):
#     residual += ∫ (dE_real/dx)(dφ/dx) dx     → contributes to A
#
#   MatReaction(E_real, rate = ε'):
#     residual  = -ε' × ∫ E_real φ dx           → tagged 'eigen' → contributes to B
#
#   CoupledForce(variable = E_real, v = E_imag, coef = +ε''):
#     residual  = -ε'' × ∫ E_imag φ dx          → tagged 'eigen' → contributes to B
#
#   Diffusion(E_imag):
#     residual += ∫ (dE_imag/dx)(dφ/dx) dx     → contributes to A
#
#   MatReaction(E_imag, rate = ε'):
#     residual  = -ε' × ∫ E_imag φ dx           → tagged 'eigen' → contributes to B
#
#   CoupledForce(variable = E_imag, v = E_real, coef = -ε''):
#     residual  = +ε'' × ∫ E_real φ dx          → tagged 'eigen' → contributes to B
#
# This assembly produces the generalised eigenproblem:
#   A·[E_r; E_i] = λ · B·[E_r; E_i]
# where the off-diagonal ε'' blocks in B create the lossy coupling.
#
# ============================================================
# BOUNDARY CONDITIONS — PEC WALLS
# ============================================================
#
# At a Perfect Electric Conductor, the tangential electric field
# must vanish: E(0) = E(D) = 0.  In the complex phasor notation:
#
#   E_real(0) = 0,   E_imag(0) = 0     at x = 0
#   E_real(D) = 0,   E_imag(D) = 0     at x = D
#
# Both DirichletBC and EigenDirichletBC are required:
#   DirichletBC      → enforces E = 0 in the stiffness matrix A
#   EigenDirichletBC → enforces E = 0 in the mass matrix B
# Without EigenDirichletBC, the boundary DOFs contribute to B
# with non-zero entries, producing spurious near-zero eigenvalues.
#
# ============================================================
# ANALYTIC PREDICTIONS
# ============================================================
#
# For a lossless cavity (ε'' = 0) the modes are:
#
#   E_m(x) = sin(mπx/D),   m = 1, 2, 3, ...
#   λ_m    = (mπ/D)² / ε'
#
# For small loss ε'' << ε', perturbation theory gives:
#
#   λ_m ≈ (mπ/D)²/ε'  ×  [1 / (1 + j ε''/ε')]
#       ≈ (mπ/D)²/ε'  ×  [1  −  j ε''/ε'  + O((ε''/ε')²)]
#
# So:
#   Re(λ_m) ≈ (mπ/D)²/ε'         (resonant wavenumber squared)
#   Im(λ_m) ≈ -(mπ/D)² ε''/(ε')² (damping rate, negative for lossy)
#
# The Q factor of the m-th mode is:
#
#   Q_m = |Re(λ_m)| / (2 |Im(λ_m)|)   ≈   ε' / (2 ε'')
#
# For ε' = 4.0, ε'' = 0.2: Q ≈ 4.0 / (2 × 0.2) = 10.
#
# Numerical values with D = 1.0:
#   λ₁_real ≈ (π)²/4  = 2.4674
#   λ₁_imag ≈ -(π)² × 0.2/16 ≈ -0.1234
#   Q₁ ≈ 2.4674 / (2 × 0.1234) ≈ 10.0
#
# ============================================================
# NUMERICAL IMPLEMENTATION NOTES
# ============================================================
#
# MatReaction (NON-AD) is used — not ADMatReaction — because:
#   (1) GenericConstantMaterial provides NON-AD material properties.
#   (2) ADMatReaction requires AD properties; mixing causes runtime error.
#   (See MEMORY.md: "AD vs Non-AD Material Properties (Critical)")
#
# CoupledForce is used for the off-diagonal ε'' blocks.
# The CoupledForce residual is: -coef * v * test
# So for the E_real equation cross-term (need -ε'' * E_imag):
#   CoupledForce(v = E_imag, coef = +eps_imag_val) contributes -ε''*E_imag ✓
# For the E_imag equation cross-term (need +ε'' * E_real):
#   CoupledForce(v = E_real, coef = -eps_imag_val) contributes +ε''*E_real ✓
#
# The eigensolver targets shift σ slightly below the expected first
# eigenvalue λ₁ ≈ 2.47 (lossless estimate).  Shift-invert with
# MUMPS as inner solver ensures robust convergence for this small
# 1D problem (200 × 2 = 400 DOFs).
#
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...})
# -----------------------------------------------------------
# D: cavity length [m]
# eps_real: real part of relative permittivity ε' = 4.0 (moderate dielectric)
# eps_imag_val: imaginary part ε'' = 0.2 (10% loss tangent: tan_d = ε''/ε' = 0.05)
# nx: number of elements (200 gives dx = 0.005 m, ~200 elements per standing wave)
#
# Mode 1 predictions:
#   λ₁_lossless = (π/D)² / ε' = 9.8696 / 4.0 = 2.4674
#   Q ≈ ε'/(2ε'') = 4.0/0.4 = 10.0
D            = 1.0    # cavity length [m]
eps_real     = 4.0    # Re(ε_r) — real part of relative permittivity
eps_imag_val = 0.2    # Im(ε_r) — imaginary part (loss term); ε_r = 4 - j0.2
nx           = 200    # number of 1D elements

[Mesh]
  # 1D domain [0, D] — the cavity interior between two PEC walls.
  #
  # GeneratedMeshGenerator in 1D creates a uniform grid of nx EDGE2 elements
  # between xmin = 0 and xmax = D.  The two endpoints are automatically
  # labelled as boundaries 'left' (x=0) and 'right' (x=D).
  #
  # Resolution: nx = 200 → dx = D/nx = 0.005 m
  # The m-th mode has spatial period 2D/m.  For mode 6 (the highest
  # we request): half-period = D/6 ≈ 0.167 m → ~33 elements per half-wave.
  # This is more than adequate for LAGRANGE/FIRST on smooth sine modes.
  [cavity_mesh]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = ${nx}
    xmin = 0
    xmax = ${D}
  []

  # Rename the auto-generated boundary labels to physically descriptive names:
  #   'pec_left'  → x = 0:  Perfect Electric Conductor wall 1
  #   'pec_right' → x = D:  Perfect Electric Conductor wall 2
  [rename_bcs]
    type         = RenameBoundaryGenerator
    input        = cavity_mesh
    old_boundary = 'left right'
    new_boundary = 'pec_left pec_right'
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor.
  # Both E_real and E_imag use first-order Lagrange (C⁰) elements.
  # The PEC BCs (E = 0) at both ends are strongly enforced via DirichletBC.
  # At convergence E_real contains the real part of the m-th eigenmode:
  #   E_real(x) ≈ sin(mπx/D) × cos(phase)
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex electric field phasor.
  # For a lossless cavity (ε'' = 0) the two equations decouple and
  # E_imag = 0 for purely real eigenmodes.  The non-zero ε'' = 0.2 here
  # creates coupling that shifts the eigenvalue into the complex plane
  # and makes E_imag non-trivial (the mode is complex-valued).
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Materials]
  # GenericConstantMaterial: uniform material properties throughout the 1D domain.
  #
  # The cavity is completely filled with the lossy dielectric (no vacuum region).
  # ε_r = ε' − j ε'' = 4.0 − j 0.2.
  #
  # Two properties are needed:
  #   'eps_r_real' = ε' = 4.0  → MatReaction rate for the diagonal mass blocks
  #
  # NOTE: GenericConstantMaterial provides NON-AD properties.
  # MatReaction (not ADMatReaction) must be used — see kernel notes.
  # CoupledForce uses a hard-coded coefficient, not a material property,
  # so no additional material property is needed for the ε'' terms.
  [dielectric_material]
    type       = GenericConstantMaterial
    prop_names  = 'eps_r_real'
    prop_values = '${eps_real}'
  []
[]

[Kernels]
  # ==================================================================
  # STIFFNESS MATRIX A — contributes to the standard residual system
  # ==================================================================
  # These kernels are NOT tagged 'eigen'; they build the A matrix of
  # the generalised eigenproblem A·x = λ·B·x.
  #
  # Diffusion kernel weak form (after integration by parts):
  #   ∫ (dE/dx)(dφ/dx) dx
  # Strong form equivalent: −d²E/dx²
  #
  # The boundary term ∫ (dE/dx) φ |₀ᴰ vanishes because:
  #   - At x = 0: φ = 0 (test function satisfies Dirichlet BC)
  #   - At x = D: φ = 0 (test function satisfies Dirichlet BC)
  # ==================================================================

  # Stiffness contribution for E_real equation (builds A block for E_real row)
  [diff_real]
    type     = Diffusion
    variable = E_real
  []

  # Stiffness contribution for E_imag equation (builds A block for E_imag row)
  [diff_imag]
    type     = Diffusion
    variable = E_imag
  []

  # ==================================================================
  # MASS MATRIX B — tagged 'eigen' → goes to the eigen (mass) system
  # ==================================================================
  # These kernels build the B matrix of A·x = λ·B·x.
  # SLEPc's KRYLOVSCHUR then finds complex pairs (λ, x) satisfying:
  #   A·x = λ · B·x
  #
  # The generalised eigenvalue λ = (ω/c)² encodes the resonance:
  #   Re(λ) ↔ resonant frequency squared
  #   Im(λ) ↔ damping rate (negative for passive lossy cavity)
  # ==================================================================

  # ------------------------------------------------------------------
  # Diagonal block: E_real equation, E_real term
  # Contribution: −ε' × ∫ E_real φ dx  →  B block −ε' × M_rr
  #
  # MatReaction residual = -reaction_rate × E_real × φ
  # With reaction_rate = eps_r_real = 4.0:
  #   B_ij (E_real block) = −∫ ε' φ_j φ_i dx  ✓
  # ------------------------------------------------------------------
  [mass_real_diag]
    type              = MatReaction
    variable          = E_real
    reaction_rate     = eps_r_real
    extra_vector_tags = 'eigen'
  []

  # ------------------------------------------------------------------
  # Off-diagonal block: E_real equation, E_imag term (ε'' coupling)
  # Contribution: −ε'' × ∫ E_imag φ dx  →  B off-diagonal block
  #
  # For the weak form of eq. (1):
  #   LHS of B for E_real eq = -ε' ∫ E_real φ dx − ε'' ∫ E_imag φ dx
  #   (The eigenvalue λ multiplies the entire B term)
  #
  # CoupledForce residual = -coef × v × test
  # With coef = +eps_imag_val = +0.2 and v = E_imag:
  #   residual = -0.2 × E_imag × φ  →  B off-diagonal: -ε'' ∫ E_imag φ dx ✓
  # ------------------------------------------------------------------
  [mass_real_offdiag]
    type              = CoupledForce
    variable          = E_real
    v                 = E_imag
    coef              = ${eps_imag_val}
    extra_vector_tags = 'eigen'
  []

  # ------------------------------------------------------------------
  # Diagonal block: E_imag equation, E_imag term
  # Contribution: −ε' × ∫ E_imag φ dx  →  B block −ε' × M_ii
  #
  # Identical structure to mass_real_diag — both equations share the
  # same diagonal ε' coefficient.
  # ------------------------------------------------------------------
  [mass_imag_diag]
    type              = MatReaction
    variable          = E_imag
    reaction_rate     = eps_r_real
    extra_vector_tags = 'eigen'
  []

  # ------------------------------------------------------------------
  # Off-diagonal block: E_imag equation, E_real term (ε'' coupling)
  # Contribution: +ε'' × ∫ E_real φ dx  →  B off-diagonal block
  #
  # For the weak form of eq. (2), the E_real coupling has opposite sign:
  #   LHS of B for E_imag eq = -ε' ∫ E_imag φ dx + ε'' ∫ E_real φ dx
  #
  # CoupledForce residual = -coef × v × test
  # With coef = -eps_imag_val = -0.2 and v = E_real:
  #   residual = -(-0.2) × E_real × φ = +0.2 × E_real × φ
  #   → B off-diagonal: +ε'' ∫ E_real φ dx ✓
  # The negative coefficient is essential: a positive coef here would
  # give the wrong sign (−ε'' instead of +ε''), breaking the physics.
  # ------------------------------------------------------------------
  [mass_imag_offdiag]
    type              = CoupledForce
    variable          = E_imag
    v                 = E_real
    coef              = ${fparse -eps_imag_val}
    extra_vector_tags = 'eigen'
  []
[]

[BCs]
  # ==================================================================
  # PEC wall at x = 0 (left boundary)
  # ==================================================================
  # A Perfect Electric Conductor forces the tangential E field to zero.
  # The phasor field must vanish at all times: E(0) = 0.
  # This means BOTH components must be zero simultaneously:
  #   E_real(0) = 0   AND   E_imag(0) = 0
  #
  # Two BC types are required for the generalised eigenproblem:
  #   DirichletBC      → zeroes E_real in the stiffness matrix A
  #                       (the standard residual system)
  #   EigenDirichletBC → zeroes E_real in the mass matrix B
  #                       (the 'eigen' tagged system)
  # Without EigenDirichletBC, the boundary nodes in B would generate
  # a contribution -ε' × E × φ at the walls, creating an inconsistency
  # between the constraint in A (E=0) and B (E free), which produces
  # spurious near-zero eigenvalues polluting the spectrum.
  # ==================================================================

  # Left wall — E_real Dirichlet (stiffness system A)
  [pec_left_real]
    type     = DirichletBC
    variable = E_real
    boundary = pec_left
    value    = 0
  []

  # Left wall — E_real EigenDirichlet (mass system B)
  [eigen_pec_left_real]
    type     = EigenDirichletBC
    variable = E_real
    boundary = pec_left
  []

  # Left wall — E_imag Dirichlet (stiffness system A)
  [pec_left_imag]
    type     = DirichletBC
    variable = E_imag
    boundary = pec_left
    value    = 0
  []

  # Left wall — E_imag EigenDirichlet (mass system B)
  [eigen_pec_left_imag]
    type     = EigenDirichletBC
    variable = E_imag
    boundary = pec_left
  []

  # ==================================================================
  # PEC wall at x = D (right boundary)
  # Same physics and reasoning as the left wall.
  # ==================================================================

  # Right wall — E_real Dirichlet (stiffness system A)
  [pec_right_real]
    type     = DirichletBC
    variable = E_real
    boundary = pec_right
    value    = 0
  []

  # Right wall — E_real EigenDirichlet (mass system B)
  [eigen_pec_right_real]
    type     = EigenDirichletBC
    variable = E_real
    boundary = pec_right
  []

  # Right wall — E_imag Dirichlet (stiffness system A)
  [pec_right_imag]
    type     = DirichletBC
    variable = E_imag
    boundary = pec_right
    value    = 0
  []

  # Right wall — E_imag EigenDirichlet (mass system B)
  [eigen_pec_right_imag]
    type     = EigenDirichletBC
    variable = E_imag
    boundary = pec_right
  []
[]

[AuxVariables]
  # E_magnitude_sq = E_real² + E_imag²
  # This is proportional to the time-averaged energy density at each point.
  # For the m-th mode: |E_m(x)|² ≈ sin²(mπx/D) in the lossless limit.
  # With loss the mode shape acquires a small complex component but
  # the magnitude pattern remains essentially sinusoidal.
  [E_magnitude_sq]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Compute |E|² = E_real² + E_imag² at each node post-solve.
  # ParsedAux evaluates the expression pointwise using the nodal DOF values.
  # This is the key diagnostic: peaks at x = D/(2m), 3D/(2m), ...
  # (antinodes of the m-th standing wave) confirm the correct mode shape.
  [E_mag_sq_aux]
    type              = ParsedAux
    variable          = E_magnitude_sq
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[VectorPostprocessors]
  # Reports all n_eigen_pairs eigenvalues after the eigensolver converges.
  # For this problem each eigenvalue is complex: λ = λ_real + j λ_imag.
  # MOOSE/SLEPc reports both parts; the CSV output contains columns for
  # the real and imaginary parts of each eigenvalue.
  #
  # From the eigenvalues the Q factor of each mode can be computed as:
  #   Q_m = |Re(λ_m)| / (2 |Im(λ_m)|)
  # Expected: Q ≈ ε'/(2 ε'') = 4.0/0.4 = 10.0 for all modes.
  # (This result is exact in the large-Q limit ε'' << ε'.)
  [eigenvalues]
    type = Eigenvalues
  []
[]

[Postprocessors]
  # ----------------------------------------------------------------
  # Mode shape diagnostics: sample E_magnitude_sq at the antinodes
  # of the first three modes to verify spatial structure.
  #
  # Mode 1: antinode at x = D/2 = 0.5       → peak of sin(πx/D)
  # Mode 2: antinodes at x = D/4 = 0.25 and x = 3D/4 = 0.75
  # Mode 3: antinode at x = D/6 = 0.1667 and x = D/2 = 0.5
  # ----------------------------------------------------------------

  # At x = D/2 = 0.5 (antinode of mode 1, mode 3, mode 5, ...)
  [E_sq_center]
    type     = PointValue
    variable = E_magnitude_sq
    point    = '0.5 0 0'
  []

  # At x = D/4 = 0.25 (antinode of mode 2, mode 6, ...)
  [E_sq_quarter]
    type     = PointValue
    variable = E_magnitude_sq
    point    = '0.25 0 0'
  []

  # At x = 3D/4 = 0.75 (second antinode of mode 2)
  [E_sq_three_quarter]
    type     = PointValue
    variable = E_magnitude_sq
    point    = '0.75 0 0'
  []

  # PEC wall sanity check: E_magnitude_sq should be ~0 at both ends
  # (enforced by DirichletBC; any non-zero value indicates BC failure)
  [E_sq_at_left_wall]
    type     = PointValue
    variable = E_magnitude_sq
    point    = '0 0 0'
  []

  [E_sq_at_right_wall]
    type     = PointValue
    variable = E_magnitude_sq
    point    = '${D} 0 0'
  []

  # Integral of |E|² over the cavity — proportional to the stored electric
  # field energy in the resonator.  For the m-th normalised mode:
  #   ∫₀ᴰ sin²(mπx/D) dx = D/2
  # The SLEPc normalisation is arbitrary, but this integral gives the
  # mode-energy scale for the eigenvector normalisation used.
  [E_sq_integral]
    type     = ElementIntegralVariablePostprocessor
    variable = E_magnitude_sq
  []
[]

[Executioner]
  type = Eigenvalue

  # KRYLOVSCHUR: SLEPc's Krylov-Schur method for generalised eigenproblems.
  # Required for finding multiple eigenvalues simultaneously (PJFNK power
  # method can only find the single dominant eigenvalue, not the full spectrum).
  solve_type = KRYLOVSCHUR

  # Request 6 eigen pairs to resolve the first six cavity modes:
  #   Mode m: λ_m ≈ (mπ/D)²/ε' = m² × π²/4 ≈ m² × 2.4674
  #   Mode 1: λ ≈  2.47   Mode 2: λ ≈  9.87   Mode 3: λ ≈ 22.21
  #   Mode 4: λ ≈ 39.48   Mode 5: λ ≈ 61.68   Mode 6: λ ≈ 88.82
  # The eigensolver finds modes near the target shift first.
  n_eigen_pairs     = 6
  which_eigen_pairs = TARGET_MAGNITUDE

  # Shift-invert spectral transformation:
  #   eps_target = 2.0: shift σ placed just below the expected first
  #                      eigenvalue λ₁ ≈ 2.467 (lossless estimate).
  #                      Shift-invert maps λ → 1/(λ − σ), so modes near
  #                      σ become the largest in the transformed problem
  #                      and the Krylov-Schur method converges to them first.
  #
  #   st_type = sinvert: activates the shift-invert spectral transformation.
  #
  #   st_ksp_type preonly + st_pc_type lu + mumps: use MUMPS direct sparse
  #     LU solver for the inner linear system (A − σB)⁻¹ at each Arnoldi
  #     step.  For this small 1D problem (400 complex DOFs) MUMPS is very
  #     fast and numerically robust.  The complex shift σ ≈ 2.0 (real) makes
  #     the inner system non-Hermitian but sparse; direct LU handles it exactly.
  petsc_options_iname = '-eps_target -st_type -st_ksp_type -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = '2.0         sinvert  preonly       lu           mumps'
[]

[Outputs]
  # exodus: 1D field profiles — E_real, E_imag, E_magnitude_sq along x.
  #   In ParaView, plot E_magnitude_sq vs x to see the standing-wave
  #   pattern of the computed eigenmode.  For mode 1: a single arch
  #   peaking at x = 0.5.  For mode 2: two arches with a node at x = 0.5.
  #   The SLEPc solver reports the mode with the eigenvalue closest to
  #   the target shift; subsequent modes are stored in higher time steps.
  exodus = true

  # csv: eigenvalue table output.
  #   The Eigenvalues VectorPostprocessor writes λ_real and λ_imag for each
  #   mode to a CSV file.  Post-process with Python:
  #     import numpy as np, pandas as pd
  #     df = pd.read_csv('case84_lossy_tem_cavity_eigenvalues_*.csv')
  #     Q_m = abs(df['real']) / (2 * abs(df['imag']))  # Q factor
  # Expected: Q ≈ 10.0 for all modes (ε'/(2ε'') = 4.0/0.4).
  csv    = true

  # execute_on = FINAL: write output only after the eigensolver converges.
  # Avoids partial output during intermediate SLEPc Krylov iterations.
  execute_on = FINAL
[]
