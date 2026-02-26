# ============================================================
# Case 35: Dispersive Pulse Broadening in an Optical Fiber
# Haus, Electromagnetic Noise and Quantum Optical Measurements (2000), Ch. 4
#
# A Gaussian pulse envelope A(x,t) propagating along a fiber broadens
# due to group velocity dispersion (GVD).  The governing equation for
# the complex envelope in the co-moving frame is:
#
#   ∂A/∂t + v_g · ∂A/∂x = D_gvd · ∂²A/∂x²
#
# where v_g is the group velocity and D_gvd is the GVD coefficient
# (= β₂/2 in standard notation).  Here we work with real amplitude
# to capture the broadening while keeping the MOOSE model simple.
#
# Analytical pulse width evolution:
#   w(z) = w₀ · √(1 + (z / z_d)²)
#
# Dispersion length (distance for √2 broadening):
#   z_d = w₀² / (2 · D_gvd) = 0.04 / (2 · 0.01) = 2.0
#
# So a pulse centred at x = 2 traveling at v_g = 1 should show
# 63 % broadening after traveling one dispersion length (z = z_d = 2,
# i.e., at t = 2), and √2 ≈ 41 % extra width at that point.
#
# MOOSE mapping:
#   ∂A/∂t          → ADTimeDerivative
#   v_g · ∂A/∂x   → ADConservativeAdvection  (velocity_material = vel_vec)
#   D_gvd · ∂²A/∂x² → ADMatDiffusion        (diffusivity = D_coeff)
#
# Initial condition: Gaussian centered at x = 2, σ² = 0.04 (σ ≈ 0.28):
#   A(x, 0) = exp(-(x-2)² / (2·0.04)) = exp(-(x-2)² / 0.08)
#
# Boundary conditions:
#   left  : Dirichlet A = 0  (pulse moves away from inlet)
#   right : zero-flux Neumann (natural, no block required)
#
# Domain: [0, 10] × [0, 0.2], quasi-1D with 100 × 2 quad elements
# Run time: t = 0 → 6.0  (pulse travels 6 units; passes z_d at t = 2)
# ============================================================

v_g   = 1.0    # group velocity
D_gvd = 0.01   # GVD coefficient  (β₂/2, units: length²/time)

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 100
  ny   = 2
  xmin = 0
  xmax = 10
  ymin = 0
  ymax = 0.2
[]

[Variables]
  # A — real pulse envelope amplitude
  [A]
  []
[]

[ICs]
  # Gaussian pulse centred at x = 2, half-width σ² = 0.04.
  # Analytic form: A(x,0) = exp(-(x-2)² / 0.08)
  [pulse]
    type     = FunctionIC
    variable = A
    function = 'exp(-(x-2)^2/0.08)'
  []
[]

[Kernels]
  # ∂A/∂t — transient accumulation term
  [time_deriv]
    type     = ADTimeDerivative
    variable = A
  []

  # v_g · ∂A/∂x — advection of the envelope at the group velocity.
  # ADConservativeAdvection implements the conservative (divergence) form
  #   -∫ ∇φ_i · (v A) dV
  # with the velocity supplied as an AD vector material property.
  [advection]
    type              = ADConservativeAdvection
    variable          = A
    velocity_material = vel_vec   # vector material property (v_g, 0, 0)
    upwinding_type    = full      # full upwinding suppresses Gibbs oscillations
  []

  # D_gvd · ∂²A/∂x² — GVD broadening term (diffusion in space)
  [diffusion]
    type        = ADMatDiffusion
    variable    = A
    diffusivity = D_coeff
  []
[]

[BCs]
  # Fix A = 0 at the inlet (left wall).  The pulse starts at x = 2
  # and moves rightward, so this boundary stays near zero throughout.
  [left]
    type     = ADDirichletBC
    variable = A
    boundary = left
    value    = 0
  []
  # Right boundary: zero-flux Neumann (natural BC, no block required).
  # Pulse exits the domain at late times without spurious reflection.
[]

[Materials]
  # Group-velocity vector supplied as an AD material property.
  # ADConservativeAdvection reads this via velocity_material = 'vel_vec'.
  [velocity_mat]
    type        = ADGenericConstantVectorMaterial
    prop_names  = 'vel_vec'
    prop_values = '${v_g} 0 0'   # v = (v_g, 0, 0)
  []

  # GVD diffusivity scalar supplied as an AD material property.
  # ADMatDiffusion reads this via diffusivity = 'D_coeff'.
  [diffusivity_mat]
    type        = ADGenericConstantMaterial
    prop_names  = 'D_coeff'
    prop_values = '${D_gvd}'
  []
[]

[Postprocessors]
  # Peak amplitude: decreases as 1/√(1 + (t/z_d)²) due to broadening
  # while area is conserved.
  [max_A]
    type       = ElementExtremeValue
    variable   = A
    value_type = max
  []

  # Total (integrated) amplitude: approximately conserved until the pulse
  # reaches the right boundary (no outflow loss for 0 ≤ t ≤ ~5).
  [total_A]
    type     = ElementIntegralVariablePostprocessor
    variable = A
  []
[]

[Executioner]
  type = Transient

  # Newton with direct LU: robust for the coupled advection-diffusion
  # system; avoids the instability of PJFNK for stiff transient advection.
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  dt       = 0.02
  end_time = 6.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true   # full field output for ParaView visualisation
  csv    = true   # postprocessor time-series (max_A, total_A)
[]
