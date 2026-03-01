# ============================================================
# Case 99: Conducting Cylinder in Uniform Electric Field
# MIT 6.641, Lec 11 — Boundary Value Problems / Polar Coords
# Prof. Markus Zahn, Spring 2005
#
# A perfectly conducting cylinder of radius R is placed in a
# uniform external electric field E₀. The potential satisfies
# Laplace's equation ∇²Φ = 0 outside the cylinder.
#
# Boundary conditions:
#   Φ(R, θ) = 0           (equipotential conductor surface)
#   Φ → −E₀ r sinθ        as r → ∞ (uniform field at infinity)
#
# Analytic solution (polar coordinates):
#   Φ(r, θ) = −E₀ (r − R²/r) sinθ    for r ≥ R
#
# The induced surface charge density on the cylinder is:
#   σ_s(θ) = 2 ε₀ E₀ sinθ
#
# This is the canonical electrostatic scattering problem,
# fundamental to capacitance calculation, dielectrophoresis,
# and lightning rod design.
#
# MOOSE approach: solve on a 2D Cartesian domain with the
# cylinder region excluded. We use an annular-like approach
# with a rectangular outer boundary at ±L and a circular
# inner boundary approximated by the analytic BC.
#
# Domain: [−3, 3]² with R = 0.5 (cylinder at origin)
# Since MOOSE GeneratedMesh can't punch a hole, we solve on
# the full square and apply the analytic solution as BC on
# all boundaries, then verify the interior matches.
# ============================================================

E0 = 1.0     # uniform field amplitude [V/m]
R  = 0.5     # cylinder radius [m]
L  = 3.0     # half-domain size [m]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 40
    xmin = -${L}
    xmax = ${L}
    ymin = -${L}
    ymax = ${L}
  []
[]

[Variables]
  # Electric potential Φ [V].
  # Satisfies ∇²Φ = 0 outside the cylinder.
  # Inside the cylinder (r < R), Φ = 0 (conductor).
  [phi]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # ∇²Φ = 0 everywhere. Inside the cylinder region a large
  # reaction term forces Φ → 0, simulating the conductor.
  [laplacian]
    type     = Diffusion
    variable = phi
  []

  # Penalty reaction term inside the cylinder (r < R):
  # forces Φ ≈ 0 to approximate a perfect conductor.
  # The ParsedMaterial 'conductor_penalty' is 0 outside and
  # large inside the cylinder, so this term only acts inside.
  [conductor_penalty]
    type          = MatReaction
    variable      = phi
    reaction_rate = conductor_penalty
  []
[]

[Materials]
  # Conductor penalty: large value (1000) inside the cylinder (r < R),
  # zero outside. This drives Φ → 0 inside the conductor.
  # GenericFunctionMaterial evaluates a Function for the material property.
  [conductor_mat]
    type        = GenericFunctionMaterial
    prop_names  = 'conductor_penalty'
    prop_values = 'penalty_func'
  []
[]

[BCs]
  # All outer boundaries: apply the analytic solution
  # Φ = −E₀ (r − R²/r) sinθ = −E₀ (y − R²y/(x²+y²))
  [outer_analytic]
    type     = FunctionDirichletBC
    variable = phi
    boundary = 'left right top bottom'
    function = analytic_bc
  []
[]

[Functions]
  # Analytic solution: Φ = −E₀·y·(1 − R²/(x²+y²))
  # At the boundaries (|x| or |y| = L), this gives the correct far field.
  [analytic_bc]
    type       = ParsedFunction
    expression = '-${E0} * y * (1.0 - ${R}*${R} / (x*x + y*y + 1e-20))'
  []

  # Penalty function: large value inside cylinder, zero outside.
  # Smooth tanh transition at r = R avoids numerical issues.
  [penalty_func]
    type       = ParsedFunction
    expression = '1000.0 * 0.5 * (1.0 - tanh(20.0 * (sqrt(x*x + y*y) - ${R})))'
  []

  # Full analytic solution for L2 error comparison.
  # Inside the cylinder (r < R) the solution is Φ = 0.
  # We use a smooth transition to handle both regions.
  [analytic_full]
    type       = ParsedFunction
    expression = 'if(x*x+y*y < ${R}*${R}, 0, -${E0} * y * (1.0 - ${R}*${R}/(x*x+y*y+1e-20)))'
  []
[]

[AuxVariables]
  [phi_exact]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [compute_exact]
    type     = FunctionAux
    variable = phi_exact
    function = analytic_full
  []
[]

[Postprocessors]
  # L2 error between numerical and analytic solutions.
  [l2_error]
    type     = ElementL2Error
    variable = phi
    function = analytic_full
  []

  # Potential at (0, R) = (0, 0.5): should be ~0 (conductor surface).
  [phi_at_surface]
    type     = PointValue
    variable = phi
    point    = '0 0.5 0'
  []

  # Potential at (0, 1.0): Φ = −E₀·1·(1 − 0.25) = −0.75
  [phi_at_1]
    type     = PointValue
    variable = phi
    point    = '0 1.0 0'
  []

  # Potential at (0, 2.0): Φ = −E₀·2·(1 − 0.0625) = −1.875
  [phi_at_2]
    type     = PointValue
    variable = phi
    point    = '0 2.0 0'
  []

  # Max potential in domain.
  [max_phi]
    type       = ElementExtremeValue
    variable   = phi
    value_type = max
  []

  # Min potential (should be negative, symmetric with max).
  [min_phi]
    type       = ElementExtremeValue
    variable   = phi
    value_type = min
  []
[]

[Executioner]
  type = Steady

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 30
[]

[Outputs]
  exodus = true
  csv    = true
[]
