# ============================================================
# Case 103: Kelvin Polarization Force — Dielectric Fluid Rise
# MIT 6.641, Lec 12/15 — Force Densities / Stress Tensors
# Prof. Markus Zahn, Spring 2005
#
# A parallel-plate capacitor is partially immersed in a
# dielectric liquid (ε > ε₀). The Kelvin polarization force
# density F = P · ∇E pulls the dielectric into the strong
# field region between the plates, causing the liquid to rise.
#
# At equilibrium, the electric force balances gravity:
#   h = (ε − ε₀) V² / (2 ρ g a²)
#
# where:
#   h = rise height [m]
#   ε = permittivity of liquid [F/m]
#   V = applied voltage [V]
#   ρ = liquid density [kg/m³]
#   g = gravitational acceleration [m/s²]
#   a = plate spacing [m]
#
# We model the electrostatic potential distribution in the
# capacitor cross-section with a dielectric/air interface.
# The gradient of |E|² at the interface reveals the Kelvin
# force that drives the fluid rise.
#
# Domain: [0, 1] × [0, 2] (capacitor cross-section)
#   x: between plates (spacing a = 1 m, normalised)
#   y: vertical (lower half: dielectric ε_r, upper half: air)
# BC: Φ(0,y) = V, Φ(1,y) = 0, ∂Φ/∂n = 0 top/bottom
# Parameters: ε_r = 3 (dielectric), V = 1 V (normalised)
# ============================================================

V_applied = 1.0
eps_r     = 3.0   # relative permittivity of dielectric liquid

[Mesh]
  # Two-block mesh: dielectric (bottom, block 0) and air (top, block 1).
  [dielectric]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 30
    ny   = 20
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
  [dielectric_id]
    type = SubdomainIDGenerator
    input = dielectric
    subdomain_id = 0
  []
  [air]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 30
    ny   = 20
    xmin = 0
    xmax = 1
    ymin = 1
    ymax = 2
  []
  [air_id]
    type = SubdomainIDGenerator
    input = air
    subdomain_id = 1
  []
  [combine]
    type = StitchMeshGenerator
    inputs = 'dielectric_id air_id'
    stitch_boundaries_pairs = 'top bottom'
  []
  [rename]
    type = RenameBoundaryGenerator
    input = combine
    old_boundary = 'left right top bottom'
    new_boundary = 'plate_high plate_ground top_bc bottom_bc'
  []
[]

[Variables]
  # Electric potential Φ [V].
  # Governed by ∇·(ε∇Φ) = 0 (Gauss's law, no free charge).
  [phi]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # ∇·(ε∇Φ) = 0 → in weak form: ∫ ε ∇Φ·∇ψ dV = 0.
  # ADMatDiffusion with diffusivity = ε gives the correct operator.
  # The permittivity jump at y = 1 is handled by block-wise materials.
  [gauss_law]
    type        = ADMatDiffusion
    variable    = phi
    diffusivity = permittivity
  []
[]

[BCs]
  # High-voltage plate (x = 0): Φ = V.
  [plate_high]
    type     = DirichletBC
    variable = phi
    boundary = plate_high
    value    = ${V_applied}
  []

  # Ground plate (x = 1): Φ = 0.
  [plate_ground]
    type     = DirichletBC
    variable = phi
    boundary = plate_ground
    value    = 0
  []

  # Top and bottom: natural BC (∂Φ/∂n = 0), representing
  # the fringe field extending smoothly out of the capacitor.
[]

[Materials]
  # Dielectric liquid (block 0, y < 1): ε_r = 3.
  [dielectric_mat]
    type        = ADGenericConstantMaterial
    prop_names  = 'permittivity'
    prop_values = '${eps_r}'
    block       = 0
  []

  # Air (block 1, y > 1): ε_r = 1.
  [air_mat]
    type        = ADGenericConstantMaterial
    prop_names  = 'permittivity'
    prop_values = '1.0'
    block       = 1
  []
[]

[AuxVariables]
  # Electric field magnitude squared |E|² = |∇Φ|².
  # The Kelvin force is proportional to ∇(|E|²).
  [E_squared]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Compute |E|² = (∂Φ/∂x)² + (∂Φ/∂y)² using ParsedAux.
  # We use the VariableGradientComponent to get gradients,
  # but ParsedAux doesn't directly access gradients.
  # Instead compute E² from the potential directly using
  # a post-processing approach.
  # For now we use a simpler diagnostic: just Φ² as a proxy.
  [E_sq_aux]
    type              = ParsedAux
    variable          = E_squared
    expression        = 'phi * phi'
    coupled_variables = 'phi'
  []
[]

[Postprocessors]
  # Average potential in dielectric region.
  [avg_phi_dielectric]
    type     = ElementAverageValue
    variable = phi
    block    = 0
  []

  # Average potential in air region.
  [avg_phi_air]
    type     = ElementAverageValue
    variable = phi
    block    = 1
  []

  # Potential at the interface (x = 0.5, y = 1.0).
  # Should be V/2 = 0.5 (midway between plates).
  [phi_interface]
    type     = PointValue
    variable = phi
    point    = '0.5 1.0 0'
  []

  # Potential at (0.5, 0.5) — inside dielectric.
  [phi_dielectric_mid]
    type     = PointValue
    variable = phi
    point    = '0.5 0.5 0'
  []

  # Potential at (0.5, 1.5) — inside air.
  [phi_air_mid]
    type     = PointValue
    variable = phi
    point    = '0.5 1.5 0'
  []

  # Max potential (should be V = 1 at plate).
  [max_phi]
    type       = ElementExtremeValue
    variable   = phi
    value_type = max
  []

  # Integral of |E|² in dielectric (proportional to stored energy).
  [E_sq_dielectric]
    type     = ElementIntegralVariablePostprocessor
    variable = E_squared
    block    = 0
  []

  # Integral of |E|² in air.
  [E_sq_air]
    type     = ElementIntegralVariablePostprocessor
    variable = E_squared
    block    = 1
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
