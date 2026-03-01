# ============================================================
# Case 94: Spatially Periodic Potential Sheet
# MIT 6.641, Lec 10 — Electrostatics / Separation of Variables
# Prof. Markus Zahn, Spring 2005
#
# A sheet of surface charge with sinusoidal distribution
# Φ(x, y=0) = V₀ sin(ax) is placed at y = 0.
# The potential satisfies Laplace's equation ∇²Φ = 0
# in the upper half-plane, with the analytic solution:
#
#   Φ(x, y) = V₀ sin(ax) exp(−a|y|)
#
# The key physics: fields decay exponentially away from
# a periodic charge sheet, with decay length 1/a = λ/(2π).
#
# Domain: [0, 2π/a] × [0, L_y], 40×40 mesh
# BC: Φ = V₀ sin(ax) at y = 0 (bottom), Φ = 0 at y = L_y (top)
#     Periodic in x (or zero Neumann — natural BC)
# Parameters: a = 2π (wavelength λ = 1), V₀ = 1.0, L_y = 2.0
# ============================================================

# Wavenumber a = 2π corresponds to wavelength λ = 1.
# Decay length = 1/a ≈ 0.159 m.
# At y = 1 the field has decayed to exp(−2π) ≈ 0.002.
V0   = 1.0
a    = 6.2831853   # 2π
L_y  = 2.0

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 40
    xmin = 0
    xmax = 1.0      # one full wavelength λ = 2π/a = 1
    ymin = 0
    ymax = ${L_y}
  []
[]

[Variables]
  # Electric potential Φ [V].
  # Governed by Laplace's equation ∇²Φ = 0.
  [phi]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # ∇²Φ = 0 → the Diffusion kernel provides ∫ ∇Φ · ∇ψ dV = 0
  # (weak form of Laplace's equation with no source term).
  [laplacian]
    type     = Diffusion
    variable = phi
  []
[]

[BCs]
  # Bottom boundary (y = 0): sinusoidal potential sheet
  # Φ(x, 0) = V₀ sin(ax) = sin(2πx) for one wavelength.
  [bottom_sinusoidal]
    type     = FunctionDirichletBC
    variable = phi
    boundary = bottom
    function = bottom_potential
  []

  # Top boundary (y = L_y): far-field condition Φ → 0.
  # At y = 2.0 the analytic solution gives Φ/V₀ = exp(−2π·2) ≈ 3.5e-6,
  # so Φ = 0 is an excellent approximation.
  [top_zero]
    type     = DirichletBC
    variable = phi
    boundary = top
    value    = 0
  []

  # Left and right boundaries: natural (zero-flux) Neumann BC.
  # Since sin(2πx) = 0 at x = 0 and x = 1 and the solution is periodic,
  # the natural BC ∂Φ/∂n = 0 is consistent at these boundaries.
  # No explicit BC needed — MOOSE applies zero-flux by default.
[]

[Functions]
  # Sinusoidal potential at the bottom boundary.
  # ParsedFunction evaluates V₀ sin(ax) = sin(2πx) at each node.
  [bottom_potential]
    type       = ParsedFunction
    expression = '${V0} * sin(${a} * x)'
  []

  # Analytic solution for validation: Φ = V₀ sin(ax) exp(−ay)
  [analytic_solution]
    type       = ParsedFunction
    expression = '${V0} * sin(${a} * x) * exp(-${a} * y)'
  []
[]

[AuxVariables]
  # Analytic solution for pointwise comparison.
  [phi_exact]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [compute_exact]
    type     = FunctionAux
    variable = phi_exact
    function = analytic_solution
  []
[]

[Postprocessors]
  # L2 error between numerical and analytic solutions.
  [l2_error]
    type     = ElementL2Error
    variable = phi
    function = analytic_solution
  []

  # Maximum potential in the domain (should be ~V₀ = 1.0 at bottom).
  [max_phi]
    type       = ElementExtremeValue
    variable   = phi
    value_type = max
  []

  # Potential at (0.25, 0): should be V₀ sin(π/2) = 1.0
  [phi_at_quarter]
    type     = PointValue
    variable = phi
    point    = '0.25 0 0'
  []

  # Potential at (0.25, 0.159): should be V₀ sin(π/2) exp(−1) ≈ 0.368
  # This is one decay length away from the sheet.
  [phi_at_one_decay]
    type     = PointValue
    variable = phi
    point    = '0.25 0.159 0'
  []

  # Potential at (0.25, 0.5): should be V₀ sin(π/2) exp(−π) ≈ 0.043
  [phi_at_half]
    type     = PointValue
    variable = phi
    point    = '0.25 0.5 0'
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
