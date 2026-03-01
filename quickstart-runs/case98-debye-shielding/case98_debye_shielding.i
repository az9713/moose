# ============================================================
# Case 98: Debye Shielding — Linearized Poisson-Boltzmann
# MIT 6.641, Lec 7 — Electrochemistry / Plasma Physics
# Prof. Markus Zahn, Spring 2005
#
# In a plasma or electrolyte, a test charge is screened by
# the surrounding mobile charges. The linearized form of the
# Poisson-Boltzmann equation is:
#
#   ∇²Φ − Φ/λ_D² = 0
#
# where λ_D = √(ε₀ k_B T / (2 n₀ q²)) is the Debye length.
#
# The solution in 2D (cylindrical-like) decays as:
#   Φ(r) ~ K₀(r/λ_D)  (modified Bessel function of 2nd kind)
# For r >> λ_D this approaches ~ (1/√r) exp(−r/λ_D).
#
# This equation has the same form as the screened Poisson
# (Helmholtz) equation and appears in semiconductor physics
# (Thomas-Fermi screening), colloid science (DLVO theory),
# and nuclear physics (Yukawa potential).
#
# Domain: [−2, 2]² square, 40×40 mesh
# BC: Φ = V₀ on a small central circle (electrode), Φ = 0 on outer boundary
# Approximation: electrode represented by a high-value source region
# Parameters: λ_D = 0.2 (normalised Debye length)
# ============================================================

# Inverse Debye length squared: 1/λ_D² = 25.
# λ_D = 0.2 (Debye length, used in comments only)
inv_lambda_D_sq = 25.0  # 1/(0.2²) = 25

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 40
    xmin = -2
    xmax = 2
    ymin = -2
    ymax = 2
  []
[]

[Variables]
  # Electric potential Φ [V].
  # Governed by the screened Poisson equation ∇²Φ − Φ/λ_D² = 0.
  [phi]
    order  = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
[]

[Kernels]
  # ∇²Φ: Laplacian (diffusion) operator.
  [laplacian]
    type     = Diffusion
    variable = phi
  []

  # −Φ/λ_D²: screening term (linear reaction).
  # CoefReaction adds ∫ coeff·Φ·ψ dV with coeff = 1/λ_D² = 25.
  # The positive coefficient acts as a sink, causing the exponential decay.
  [screening]
    type        = CoefReaction
    variable    = phi
    coefficient = ${inv_lambda_D_sq}
  []

  # Localised source: a Gaussian "electrode" at the origin.
  # This mimics a point charge by adding a smooth source term
  # ∫ f(x,y) ψ dV where f is a narrow Gaussian.
  [source]
    type     = BodyForce
    variable = phi
    function = source_function
  []
[]

[Functions]
  # Narrow Gaussian centred at origin with width σ = 0.1.
  # Amplitude chosen to give Φ_max ≈ 1 at origin.
  [source_function]
    type       = ParsedFunction
    expression = '500.0 * exp(-((x*x + y*y)) / (2 * 0.01))'
  []
[]

[BCs]
  # All four outer boundaries: Φ = 0 (far-field, fully screened).
  # The domain extends 10 Debye lengths from the source, so the
  # field has decayed to exp(−10) ≈ 5e-5.
  [outer_zero]
    type     = DirichletBC
    variable = phi
    boundary = 'left right top bottom'
    value    = 0
  []
[]

[Postprocessors]
  # Peak potential at the origin (source location).
  [max_phi]
    type       = ElementExtremeValue
    variable   = phi
    value_type = max
  []

  # Potential at one Debye length from origin.
  [phi_at_1lambda]
    type     = PointValue
    variable = phi
    point    = '0.2 0 0'
  []

  # Potential at two Debye lengths.
  [phi_at_2lambda]
    type     = PointValue
    variable = phi
    point    = '0.4 0 0'
  []

  # Potential at five Debye lengths.
  [phi_at_5lambda]
    type     = PointValue
    variable = phi
    point    = '1.0 0 0'
  []

  # Domain integral of Φ (total electrostatic energy proxy).
  [phi_integral]
    type     = ElementIntegralVariablePostprocessor
    variable = phi
  []
[]

[Executioner]
  type = Steady

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv    = true
[]
