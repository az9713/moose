# ============================================================
# Case 22: Charge Relaxation in an Ohmic Medium
# Melcher, Continuum Electromechanics (MIT, 1981), Ch. 5 §5.9
#
# Free charge in a conducting medium decays exponentially:
#   ∂ρ_e/∂t + (σ/ε)·ρ_e = 0    (relaxation, τ = ε/σ)
#   -div(ε·grad(φ)) = ρ_e       (Poisson equation)
#
# Domain: unit square [0,1]², 30×30 mesh
# IC: Gaussian charge blob at center, φ = 0
# BC: φ = 0 on all walls (grounded enclosure)
# Parameters: ε = 1, σ = 10  →  τ_e = 0.1 s
# ============================================================

# HIT top-level variable: σ/ε ratio.
# With ε = 1 and σ = 10 the charge relaxation time is τ = ε/σ = 0.1 s.
# Changing this single value re-scales the entire decay rate.
sigma_over_eps = 10.0

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 30
    ny   = 30
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[Variables]
  # Free charge density ρ_e [C/m²] — the primary evolution variable.
  # Governed by the ODE ∂ρ_e/∂t + (σ/ε)·ρ_e = 0 at each point.
  [rho_e]
    order  = FIRST
    family = LAGRANGE
    [InitialCondition]
      # Gaussian blob centred at (0.5, 0.5) with width √0.01 = 0.1.
      # Peak value of 1 provides a non-trivial initial charge distribution.
      type        = FunctionIC
      function    = 'exp(-((x-0.5)^2+(y-0.5)^2)/0.01)'
    []
  []

  # Electric potential φ [V] — solved quasi-statically via Poisson's equation.
  # The potential responds instantly to the current charge distribution.
  [phi]
    order  = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
[]

[Kernels]
  # -------------------------------------------------------
  # Charge relaxation equation for ρ_e:
  #   ∂ρ_e/∂t + (σ/ε)·ρ_e = 0
  # -------------------------------------------------------

  # ∂ρ_e/∂t: standard first-order time derivative.
  # This is the accumulation term — charges cannot teleport.
  [rho_e_time]
    type     = ADTimeDerivative
    variable = rho_e
  []

  # (σ/ε)·ρ_e: linear decay (Ohmic conduction drains the free charge).
  # ADReaction adds  ∫ rate·ρ_e·φ_i dV  to the residual, which represents
  # the volumetric loss of charge due to the finite conductivity of the medium.
  [rho_e_decay]
    type     = ADReaction
    variable = rho_e
    rate     = ${sigma_over_eps}   # σ/ε = 10 → τ = 0.1 s
  []

  # -------------------------------------------------------
  # Poisson equation for φ:
  #   -div(ε·grad(φ)) = ρ_e
  # -------------------------------------------------------

  # -div(ε·grad(φ)): diffusion-like operator with the permittivity ε playing
  # the role of thermal conductivity.  ADHeatConduction adds
  #   ∫ ε·grad(φ)·grad(φ_i) dV
  # which is the weak form of −div(ε·grad(φ)).
  # The material property name 'permittivity' is wired in via the
  # thermal_conductivity parameter so no rename is needed in [Materials].
  [phi_laplacian]
    type                 = ADHeatConduction
    variable             = phi
    thermal_conductivity = permittivity
  []

  # ρ_e source on the right-hand side of Poisson's equation.
  # ADCoupledForce adds  -∫ ρ_e·φ_i dV  to the residual for φ,
  # which represents the charge density driving the potential.
  [phi_source]
    type     = ADCoupledForce
    variable = phi
    v        = rho_e
  []
[]

[BCs]
  # φ = 0 on all four walls: grounded conducting enclosure.
  # Free charge inside cannot escape through the boundaries
  # (zero-flux Neumann on ρ_e is the natural BC — no explicit entry needed).
  [phi_ground]
    type     = ADDirichletBC
    variable = phi
    boundary = 'left right top bottom'
    value    = 0.0
  []
[]

[Materials]
  # Permittivity ε = 1 [C²/(N·m²)] — used by ADHeatConduction (phi_laplacian).
  # The property name 'permittivity' is referenced directly by the kernel via
  # thermal_conductivity = permittivity, so this is the only material needed.
  [permittivity]
    type        = ADGenericConstantMaterial
    prop_names  = 'permittivity'
    prop_values = '1.0'
  []
[]

[Postprocessors]
  # Domain-average free charge density: should decay as exp(−t/τ) = exp(−10t).
  # Monitoring this confirms the correct exponential relaxation rate.
  [avg_rho_e]
    type     = ElementAverageValue
    variable = rho_e
  []

  # Peak free charge density: tracks the maximum of the Gaussian blob.
  # The spatial shape remains Gaussian (linear PDE), so this also decays as exp(−10t).
  [max_rho_e]
    type       = ElementExtremeValue
    variable   = rho_e
    value_type = max
  []

  # Maximum electric potential: the potential bump is driven solely by ρ_e,
  # so it collapses to zero at the same rate the charge dissipates.
  [max_phi]
    type       = ElementExtremeValue
    variable   = phi
    value_type = max
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # Fixed timestep dt = 0.01 s over 0.5 s total (50 steps).
  # dt/τ = 0.01/0.1 = 0.1 — well-resolved temporal resolution.
  # By t = 0.5 s = 5τ the charge has decayed to exp(−5) ≈ 0.007.
  dt       = 0.01
  end_time = 0.5

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20
[]

# Full single-matrix preconditioning captures the rho_e–phi coupling in the
# Jacobian, which is essential for tight convergence of the coupled system.
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Outputs]
  exodus = true   # Full spatial fields (rho_e and phi) at every timestep
  csv    = true   # Postprocessor history (avg_rho_e, max_rho_e, max_phi vs. time)
[]
