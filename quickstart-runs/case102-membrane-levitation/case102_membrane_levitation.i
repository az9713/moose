# ============================================================
# Case 102: Electric Field Levitation of Membrane
# MIT 6.641, Lec 18 — Stability / Instabilities
# Prof. Markus Zahn, Spring 2005
#
# A thin elastic membrane (string) of tension S is suspended
# between two fixed supports (x = 0, x = l) at gap distance d
# above a grounded electrode. A voltage V is applied, creating
# an attractive electrostatic force that opposes the tension.
#
# The linearized equation of motion about equilibrium is:
#   ρ ∂²ξ/∂t² = S ∂²ξ/∂x² + (ε₀V²/d³) ξ
#
# The electrostatic term (ε₀V²/d³)ξ acts as a NEGATIVE
# spring constant — it destabilises the membrane.
#
# Eigenvalue analysis: find the natural frequencies ω_n:
#   −S ∂²ξ/∂x² − F_e ξ = ω²ρ ξ    (F_e = ε₀V²/d³)
#
# Eigenvalues: ω²_n = S(nπ/l)² − F_e    (with ρ=1)
#
# Pull-in instability occurs when ω²₁ < 0:
#   F_e > S(π/l)²   →   V > V_crit = d^{3/2}·π/l·√(S/ε₀)
#
# With normalised S=1, l=1, d=1, ε₀=1:
#   V_crit = π ≈ 3.14159
#
# We solve at V = 2 (stable, below V_crit):
#   ω²₁ = π² − 4 ≈ 5.87   (positive → stable oscillation)
#   ω²₂ = 4π² − 4 ≈ 35.48
#   ω²₃ = 9π² − 4 ≈ 84.78
#
# Domain: 1D [0, 1], 200 elements
# BC: ξ(0) = ξ(1) = 0 (clamped ends)
# ============================================================

# Normalised parameters.
# S_tension = 1.0 (membrane tension, used as Diffusion coefficient = 1)
F_e          = 4.0       # destabilising electric force ε₀V²/d³ (V=2)
mass_density = 1.0       # linear mass density ρ [kg/m]
nx_mesh      = 200       # mesh resolution

[Mesh]
  [line]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = ${nx_mesh}
    xmin = 0
    xmax = 1.0
  []
  [rename]
    type = RenameBoundaryGenerator
    input = line
    old_boundary = 'left right'
    new_boundary = 'clamp_left clamp_right'
  []
[]

[Variables]
  # Membrane displacement ξ [m] from equilibrium.
  # The eigensolver finds the normal modes sin(nπx/l)
  # and their associated frequencies ω_n.
  [xi]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # =========================================================
  # STIFFNESS MATRIX A (untagged kernels)
  # =========================================================

  # S ∂²ξ/∂x² → in weak form: S·∫ ∂ξ/∂x · ∂ψ/∂x dV
  # This is the restoring (stabilising) membrane tension.
  [tension_stiffness]
    type     = Diffusion
    variable = xi
  []

  # −F_e·ξ → destabilising electrostatic force.
  # CoefReaction adds coeff·∫ ξ·ψ dV to the residual.
  # With coeff = −F_e = −4, this REDUCES the effective stiffness.
  [electric_destabilise]
    type        = CoefReaction
    variable    = xi
    coefficient = -${F_e}
  []

  # =========================================================
  # MASS MATRIX B (tagged 'eigen')
  # =========================================================

  # ρ·ξ mass term: MatReaction adds −rate·∫ ξ·ψ dV.
  # Tagged 'eigen' so it goes into the B (mass) system.
  # The eigenvalue problem becomes: A·x = ω²·B·x
  # where ω² is the squared natural frequency.
  [mass]
    type              = MatReaction
    variable          = xi
    reaction_rate     = mass_rho
    extra_vector_tags = 'eigen'
  []
[]

[Materials]
  # Mass density for the B matrix.
  [density_mat]
    type        = GenericConstantMaterial
    prop_names  = 'mass_rho'
    prop_values = '${mass_density}'
  []
[]

[BCs]
  # =========================================================
  # Clamped boundary conditions: ξ = 0 at both ends.
  # Both DirichletBC (for A) and EigenDirichletBC (for B)
  # are required to avoid spurious near-zero eigenvalues.
  # =========================================================

  [clamp_left]
    type     = DirichletBC
    variable = xi
    boundary = clamp_left
    value    = 0
  []
  [eigen_clamp_left]
    type     = EigenDirichletBC
    variable = xi
    boundary = clamp_left
  []

  [clamp_right]
    type     = DirichletBC
    variable = xi
    boundary = clamp_right
    value    = 0
  []
  [eigen_clamp_right]
    type     = EigenDirichletBC
    variable = xi
    boundary = clamp_right
  []
[]

[VectorPostprocessors]
  # Reports all eigenvalues ω²_n after the eigensolver converges.
  [eigenvalues]
    type = Eigenvalues
  []
[]

[Postprocessors]
  # Mode shape diagnostics: |ξ|² at key points.

  # Midpoint (x = 0.5): antinode of mode 1, node of mode 2.
  [xi_sq_mid]
    type     = PointValue
    variable = xi
    point    = '0.5 0 0'
  []

  # Quarter point (x = 0.25): antinode of mode 2.
  [xi_sq_quarter]
    type     = PointValue
    variable = xi
    point    = '0.25 0 0'
  []
[]

[Executioner]
  type = Eigenvalue

  # KRYLOVSCHUR: Krylov-Schur eigensolver from SLEPc.
  solve_type = KRYLOVSCHUR

  # Request 4 eigen pairs to see the first four modes.
  # Expected eigenvalues (ω²_n = S(nπ)² − F_e, with S=1, F_e=4):
  #   ω²₁ = π² − 4  ≈  5.87  (mode 1)
  #   ω²₂ = 4π² − 4 ≈ 35.48  (mode 2)
  #   ω²₃ = 9π² − 4 ≈ 84.78  (mode 3)
  #   ω²₄ = 16π²− 4 ≈ 153.91 (mode 4)
  n_eigen_pairs     = 4
  which_eigen_pairs = TARGET_MAGNITUDE

  # Shift-invert targeting near the first eigenvalue ≈ 5.
  petsc_options_iname = '-eps_target -st_type -st_ksp_type -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = '5.0         sinvert  preonly       lu           mumps'
[]

[Outputs]
  exodus     = true
  csv        = true
  execute_on = FINAL
[]
