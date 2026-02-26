# ============================================================
# Case 24: Charge Drift-Diffusion Between Parallel Plates
# Melcher, Continuum Electromechanics (MIT, 1981), Ch. 5 §5.5-5.7
#
# Unipolar charge injection from an electrode into a dielectric:
#   ∂ρ_e/∂t + div(μ_i · E · ρ_e) = D_i · ∇²ρ_e  (drift-diffusion)
#   -div(ε · grad(φ)) = ρ_e                         (Poisson)
#   E = -grad(φ)
#
# Drift velocity approximation: v = μ_i · E_applied (uniform prescribed field)
# Valid when space charge is weak relative to the applied field —
# the standard starting point for space-charge-limited conduction analysis.
#
# Domain: [0,1] × [0,0.067] quasi-1D slab, 60×4 mesh
# IC: ρ_e = 0 (charge-free gap at t=0), φ = 0
# BC: ρ_e = 1 at left (injection electrode), ρ_e = 0 at right (collecting electrode)
#     φ = 1 at left (anode), φ = 0 at right (cathode)
# Parameters: μ_i = 1.0 m²/(V·s), D_i = 0.01 m²/s, ε = 1.0, V_applied = 1.0 V
#
# Physical picture: positive ions are injected at the left electrode and drift
# rightward under the applied electric field. Diffusion broadens the charge
# front. Poisson's equation is solved simultaneously to capture the self-consistent
# potential distortion caused by the accumulating space charge.
# ============================================================

# --- Top-level parameters (modify here to explore parameter space) ---
mu_ion = 1.0     # ion mobility [m²/(V·s)]
D_ion  = 0.01    # ion diffusivity [m²/s]
eps    = 1.0     # permittivity [C²/(N·m²)] (dimensionless normalised units)

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 60      # 60 elements along x for good drift-front resolution
    ny   = 4       # thin slab: physics independent of y
    xmin = 0
    xmax = 1       # inter-electrode gap d = 1 m (normalised)
    ymin = 0
    ymax = 0.067   # quasi-1D strip width
  []
[]

[Variables]
  # Free charge density ρ_e [C/m³] — the drifting and diffusing quantity.
  # Zero initial condition: no charge in the gap before injection begins.
  [rho_e]
    order  = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []

  # Electric potential φ [V] — driven by ρ_e via Poisson's equation.
  # Solved simultaneously with ρ_e to capture the coupled field distortion.
  [phi]
    order  = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
[]

[Kernels]
  # -------------------------------------------------------
  # Drift-diffusion equation for ρ_e:
  #   ∂ρ_e/∂t + div(v · ρ_e) - D_i · ∇²ρ_e = 0
  # where v = -μ_i · grad(φ) ≈ μ_i · E_applied (uniform, prescribed)
  # -------------------------------------------------------

  # ∂ρ_e/∂t: accumulation term. AD-compatible time derivative.
  [rho_time]
    type     = ADTimeDerivative
    variable = rho_e
  []

  # D_i · ∇²ρ_e: Fickian diffusion of charge.
  # ADMatDiffusion reads the scalar material property 'D_ion'.
  # Without diffusion the charge front would be a perfect step; diffusion
  # rounds and broadens it, which is physically correct and improves stability.
  [rho_diffusion]
    type        = ADMatDiffusion
    variable    = rho_e
    diffusivity = D_ion   # material property defined in [Materials]
  []

  # div(v · ρ_e): conservative advection with upwinding.
  # ConservativeAdvection with velocity_as_variable_gradient uses the gradient
  # of phi as the advection velocity: v = coef * grad(φ).
  # The ion drift velocity is v = μ_i · E = -μ_i · grad(φ), so coef = -μ_i.
  # velocity_scalar_coef references a material property (neg_mobility = -1.0).
  # Full upwinding adds numerical diffusion to prevent oscillations at the
  # sharp charge front without losing conservation.
  [rho_advection]
    type                          = ConservativeAdvection
    variable                      = rho_e
    velocity_as_variable_gradient = phi
    velocity_scalar_coef          = neg_mobility
    upwinding_type                = full
  []

  # -------------------------------------------------------
  # Poisson equation for φ:
  #   -div(ε · grad(φ)) = ρ_e
  # -------------------------------------------------------

  # -div(ε · grad(φ)): Laplacian operator weighted by permittivity.
  # ADHeatConduction implements ∫ k · grad(T) · grad(φ_i) dV, which is the
  # weak form of -div(k · grad(T)). Here k → ε and T → φ, so the mapping
  # is exact. The property name 'permittivity' is passed via thermal_conductivity.
  [phi_laplacian]
    type                 = ADHeatConduction
    variable             = phi
    thermal_conductivity = permittivity   # maps ε to the diffusion coefficient
  []

  # ρ_e source on the right-hand side of Poisson's equation.
  # ADCoupledForce adds -∫ ρ_e · φ_i dV to the residual for φ, which places
  # the charge density on the RHS of the Poisson equation as a source term.
  # The space charge bends the potential away from the linear Laplace solution.
  [phi_source]
    type     = ADCoupledForce
    variable = phi
    v        = rho_e
  []
[]

[BCs]
  # --- Charge density boundary conditions ---

  # Left electrode (x=0): continuous injection of charge.
  # ρ_e = 1 is the prescribed injection density — the source of all charge
  # that enters the gap. This models a strong-injection (ohmic-contact) electrode.
  [rho_inject]
    type     = ADDirichletBC
    variable = rho_e
    boundary = left
    value    = 1.0
  []

  # Right electrode (x=1): collecting electrode absorbs all arriving charge.
  # ρ_e = 0 is the sink condition — the electrode immediately neutralises
  # any charge that reaches it, keeping the surface charge-free.
  [rho_collect]
    type     = ADDirichletBC
    variable = rho_e
    boundary = right
    value    = 0.0
  []

  # --- Electric potential boundary conditions ---

  # Left electrode: anode at high potential (drives ion injection).
  [phi_anode]
    type     = ADDirichletBC
    variable = phi
    boundary = left
    value    = 1.0   # V_applied = 1 V (normalised)
  []

  # Right electrode: cathode grounded (reference potential).
  [phi_cathode]
    type     = ADDirichletBC
    variable = phi
    boundary = right
    value    = 0.0
  []
[]

[Materials]
  # Ion diffusivity D_i [m²/s].
  # ADGenericConstantMaterial is the AD-compatible version of GenericConstantMaterial.
  # The property name 'D_ion' is consumed by ADMatDiffusion (rho_diffusion kernel).
  [ion_diffusivity]
    type        = ADGenericConstantMaterial
    prop_names  = 'D_ion'
    prop_values = '${D_ion}'   # 0.01 m²/s — diffusion broadens the charge front
  []

  # Permittivity ε [C²/(N·m²)] — consumed by ADHeatConduction (phi_laplacian kernel).
  # Here ε = 1 (normalised). In physical units, ε₀ = 8.85e-12 F/m.
  [permittivity_mat]
    type        = ADGenericConstantMaterial
    prop_names  = 'permittivity'
    prop_values = '${eps}'   # 1.0 (normalised)
  []

  # Negative ion mobility: scalar coefficient for the drift velocity.
  # velocity_as_variable_gradient uses v = coef * grad(φ).
  # Ion drift is v = -μ_i * grad(φ), so coef = -μ_i = -1.0.
  # This is now SELF-CONSISTENT: the actual gradient of the solved
  # potential φ is used (not a prescribed uniform field), so the
  # drift velocity responds to space-charge distortion of the field.
  [neg_mobility_mat]
    type        = GenericConstantMaterial
    prop_names  = 'neg_mobility'
    prop_values = '-${mu_ion}'   # -μ_i = -1.0 m²/(V·s)
  []
[]

[Postprocessors]
  # Domain-average charge density: rises from 0 as the charge front crosses the gap.
  # At steady state, the average equals ρ_inject / 2 ≈ 0.5 for uniform drift.
  [avg_rho]
    type     = ElementAverageValue
    variable = rho_e
  []

  # Maximum charge density: should stay near 1.0 (injection value) once the
  # injecting-electrode BC establishes itself. Confirms injection is active.
  [max_rho]
    type       = ElementExtremeValue
    variable   = rho_e
    value_type = max
  []

  # Maximum electric potential: should remain near 1.0 (anode value).
  # Distortion from 1.0 indicates strong space-charge feedback on the field.
  [max_phi]
    type       = ElementExtremeValue
    variable   = phi
    value_type = max
  []

  # Minimum electric potential: monitors the cathode side — should stay near 0.0.
  [min_phi]
    type       = ElementExtremeValue
    variable   = phi
    value_type = min
  []
[]

[Executioner]
  type = Transient

  # Newton with exact AD Jacobians.
  # The full SMP preconditioner (below) ensures the rho_e–phi coupling is
  # captured in the Jacobian for tight convergence.
  solve_type          = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # dt = 0.02 s — CFL-like criterion: dt < dx/v_drift = (1/60)/1 ≈ 0.017 s.
  # Using dt = 0.02 with the implicit BDF1 scheme is stable for any Courant
  # number (implicit methods are unconditionally stable), but the BDF1 temporal
  # truncation error grows with dt. 0.02 is a good balance of speed and accuracy.
  dt       = 0.02
  end_time = 1.0    # one convective transit time: t_transit = d/v = 1/1 = 1 s

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20
[]

# Full single-matrix preconditioning.
# The SMP block with full = true assembles off-diagonal Jacobian blocks that
# couple rho_e and phi. Without this, each variable is preconditioned in
# isolation and the rho_e–phi coupling is handled only approximately,
# leading to more Newton iterations or even divergence.
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Outputs]
  exodus = true   # Full spatial fields (rho_e, phi) at every timestep for ParaView
  csv    = true   # Postprocessor time history (avg_rho, max_rho, max_phi, min_phi vs. t)
[]
