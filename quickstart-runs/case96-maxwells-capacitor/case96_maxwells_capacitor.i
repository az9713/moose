# ============================================================
# Case 96: Maxwell's Capacitor — Two-Layer Dielectric
# MIT 6.641, Lec 7 — Interfacial Polarization
# Prof. Markus Zahn, Spring 2005
#
# A step voltage V is applied across two dielectric layers
# with different permittivities (ε_a, ε_b) and conductivities
# (σ_a, σ_b). Interfacial charge builds up at the boundary
# between the two layers:
#
#   σ_s(t) = V(ε_b σ_a − ε_a σ_b) / (σ_a d_b + σ_b d_a)
#            × (1 − exp(−t/τ))
#
# where τ = (ε_a d_b + ε_b d_a) / (σ_a d_b + σ_b d_a)
#
# This is critical for high-voltage insulation design where
# mismatched ε,σ layers cause interfacial charge buildup.
#
# We model this as a 1D transient problem for the electric
# potential Φ, with ohmic conduction:
#   ∂(ε ∂Φ/∂x)/∂t + ∂(σ ∂Φ/∂x)/∂x = 0  (continuity)
#
# Simplified approach: solve ∂Φ/∂t = (σ/ε) ∂²Φ/∂x²
# in each layer with matched flux at the interface.
#
# Domain: 1D [0, d_a+d_b], two blocks
# BC: Φ(0) = V, Φ(d_a+d_b) = 0
# Parameters: d_a = d_b = 0.5
#   Layer a: ε_a = 2, σ_a = 1e-3  (low conductivity)
#   Layer b: ε_b = 5, σ_b = 1e-2  (higher conductivity)
#   τ = (ε_a·d_b + ε_b·d_a)/(σ_a·d_b + σ_b·d_a)
#     = (2×0.5 + 5×0.5)/(1e-3×0.5 + 1e-2×0.5)
#     = 3.5 / 0.0055 = 636.4 s
# ============================================================

# We rescale for a fast simulation by using σ/ε ratios directly.
# Layer a: diffusivity = σ_a/ε_a = 0.5, Layer b: σ_b/ε_b = 2.0
# This gives τ ≈ 3.5/5.5 ≈ 0.636 s (with normalised σ,ε).
V_applied = 1.0

[Mesh]
  # Two-block 1D mesh: layer_a from 0 to 0.5, layer_b from 0.5 to 1.0.
  [layer_a]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = 50
    xmin = 0
    xmax = 0.5
  []
  [layer_a_id]
    type = SubdomainIDGenerator
    input = layer_a
    subdomain_id = 0
  []
  [layer_b]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = 50
    xmin = 0.5
    xmax = 1.0
  []
  [layer_b_id]
    type = SubdomainIDGenerator
    input = layer_b
    subdomain_id = 1
  []
  [combine]
    type = StitchMeshGenerator
    inputs = 'layer_a_id layer_b_id'
    stitch_boundaries_pairs = 'right left'
  []
  # Rename boundaries for clarity
  [rename]
    type = RenameBoundaryGenerator
    input = combine
    old_boundary = 'left right'
    new_boundary = 'anode cathode'
  []
[]

[Variables]
  # Electric potential Φ [V].
  # Transient: ε ∂²Φ/∂x∂t = −∂J/∂x where J = σ E = −σ ∂Φ/∂x
  # Combined: this reduces to solving the conduction equation
  # with material-dependent σ/ε ratio.
  [phi]
    order  = FIRST
    family = LAGRANGE
    initial_condition = 0
  []
[]

[Kernels]
  # Time derivative: ε ∂Φ/∂t contributes via material-weighted time derivative.
  # We use ADTimeDerivative with a coefficient from the material system.
  [phi_time]
    type     = ADTimeDerivative
    variable = phi
  []

  # Conduction term: (σ/ε) ∂²Φ/∂x².
  # ADMatDiffusion with diffusivity = σ/ε gives the correct spatial operator.
  [phi_conduction]
    type        = ADMatDiffusion
    variable    = phi
    diffusivity = sigma_over_eps
  []
[]

[BCs]
  # Anode (x = 0): step voltage V = 1.0 applied at t = 0.
  [anode]
    type     = DirichletBC
    variable = phi
    boundary = anode
    value    = ${V_applied}
  []

  # Cathode (x = 1.0): grounded.
  [cathode]
    type     = DirichletBC
    variable = phi
    boundary = cathode
    value    = 0
  []
[]

[Materials]
  # Layer a (block 0): σ_a/ε_a = 1/2 = 0.5
  [layer_a_props]
    type        = ADGenericConstantMaterial
    prop_names  = 'sigma_over_eps'
    prop_values = '0.5'
    block       = 0
  []

  # Layer b (block 1): σ_b/ε_b = 10/5 = 2.0
  [layer_b_props]
    type        = ADGenericConstantMaterial
    prop_names  = 'sigma_over_eps'
    prop_values = '2.0'
    block       = 1
  []
[]

[Postprocessors]
  # Potential at the interface (x = 0.5): tracks the transient
  # redistribution of voltage between the two layers.
  # Initially Φ(0.5) follows the capacitive divider: ε_b/(ε_a+ε_b) = 5/7 ≈ 0.714
  # At steady state it follows the resistive divider: σ_b/(σ_a+σ_b) → σ ratio
  [phi_interface]
    type     = PointValue
    variable = phi
    point    = '0.5 0 0'
  []

  # Average potential in layer a.
  [avg_phi_a]
    type     = ElementAverageValue
    variable = phi
    block    = 0
  []

  # Average potential in layer b.
  [avg_phi_b]
    type     = ElementAverageValue
    variable = phi
    block    = 1
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  # Run for 5τ ≈ 3.2 s to see full transient and approach steady state.
  # τ ≈ 0.636 s, so we run to 3.2 s with dt = 0.02 s (160 steps).
  dt       = 0.02
  end_time = 3.2

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
