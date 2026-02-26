# ============================================================
# Case 36: Soliton Pulse Propagation — Nonlinear Balance
# Haus, Electromagnetic Noise and Quantum Optical Measurements
#   (Springer, 2000), Ch. 10
#
# The nonlinear Schrodinger equation for the real pulse envelope A(x,t)
# combines group-velocity dispersion with Kerr (self-phase-modulation)
# nonlinearity:
#
#   ∂A/∂t + v_g · ∂A/∂x = D · ∂²A/∂x² − α · A³
#
# At the soliton balance condition the dispersive broadening is exactly
# cancelled by the nonlinear self-steepening, and the pulse propagates
# without changing shape.
#
# Soliton balance condition:
#   α · A₀² · w₀² = 2 · D
#
# With A₀ = 1 (peak amplitude), w₀ = 1 (half-width), D = 0.05:
#   α_soliton = 2 · D / (A₀² · w₀²) = 2 · 0.05 / (1 · 1) = 0.1
#
# The exact soliton initial condition is a sech profile:
#   A(x, 0) = A₀ / cosh((x − x₀) / w₀) = 1 / cosh(x − 5)
#
# Expected behaviour with α = 0.1 (soliton balance):
#   - Peak amplitude max_A stays near 1.0 throughout t ∈ [0, 10]
#   - Total integrated amplitude total_A is approximately conserved
#   - Pulse centroid advances at v_g = 1 (one unit per unit time)
#
# Regime comparison (change alpha at the top of this file):
#   alpha = 0.0  → pure GVD dispersion (reproduces Case 35 with D = 0.05)
#                  peak decays as 1/√(1 + t/z_d),  z_d = w₀²/(2D) = 10
#   alpha = 0.1  → soliton balance: peak stays constant (this case)
#   alpha = 0.3  → over-nonlinear: pulse compresses then develops
#                  Akhmediev-breather-like oscillations around a narrower core
#
# MOOSE kernel mapping:
#   ∂A/∂t               → ADTimeDerivative
#   v_g · ∂A/∂x         → ADConservativeAdvection  (velocity_material = vel_vec)
#   D · ∂²A/∂x²         → ADMatDiffusion           (diffusivity = D_coeff)
#   −α · A³              → ADMatReaction            (reaction_rate = nl_rate)
#
# ADMatReaction sign convention:
#   Residual contribution = −reaction_rate · A · φ_i
#   With nl_rate = −α · A²:  residual = +α · A³ · φ_i
#   Strong-form PDE term:    −α · A³   (correct self-focusing sign)
#
# DerivativeParsedMaterial JIT workaround (Docker compatibility):
#   MOOSE's expression parser can invoke an LLVM JIT compiler that may be
#   unavailable or crash inside Docker containers without a full toolchain.
#   Setting disable_fpoptimizer = true and enable_jit = false forces the
#   interpreter path, which is slower but fully portable.
#
# Domain:     [0, 20] × [0, 0.2],  quasi-1D,  100 × 2 QUAD4 elements
# Run time:   t = 0 → 10  (pulse travels 10 units at v_g = 1)
# ============================================================

# HIT top-level variables — all three must be referenced below or
# MOOSE will error on "unused parameter".
v_g   = 1.0    # group velocity  [length/time]
D     = 0.05   # GVD coefficient [length²/time]  (dispersion broadening)
alpha = 0.1    # Kerr nonlinearity coefficient   (soliton balance at 0.1)

[Mesh]
  type = GeneratedMesh
  dim  = 2
  # 100 elements along x give h = 0.2 per element; the sech half-width
  # w₀ = 1 spans ~5 elements — adequate resolution for a smooth sech profile.
  # 2 elements in y make this a quasi-1D domain; only x-direction physics matter.
  nx   = 100
  ny   = 2
  xmin = 0
  xmax = 20
  ymin = 0
  ymax = 0.2
[]

[Variables]
  # A — real pulse envelope amplitude (sech-shaped soliton).
  # Initialised by [ICs] below to A(x,0) = 1/cosh(x−5).
  [A]
  []
[]

[ICs]
  # Exact soliton initial condition: A₀/cosh((x − x₀)/w₀)
  # with A₀ = 1, x₀ = 5, w₀ = 1.
  # The pulse is centred far enough from the left boundary (x = 5)
  # that the Dirichlet BC at x = 0 has negligible effect initially.
  [soliton]
    type     = FunctionIC
    variable = A
    function = 'cosh_inv'
  []
[]

[Functions]
  # sech profile: 1/cosh((x−5)/1) = 2/(exp(x−5) + exp(−(x−5)))
  # ParsedFunction evaluates this analytically at every node.
  [cosh_inv]
    type       = ParsedFunction
    expression = '1.0/cosh((x-5.0)/1.0)'
  []
[]

[Kernels]
  # ∂A/∂t — transient accumulation term.
  # ADTimeDerivative uses automatic differentiation for the exact Jacobian;
  # no separate Jacobian assembly is needed.
  [time_deriv]
    type     = ADTimeDerivative
    variable = A
  []

  # v_g · ∂A/∂x — advection of the soliton envelope at the group velocity.
  # ADConservativeAdvection implements the conservative (divergence) form:
  #   ∫ A · (v · ∇φ_i) dV  (integrated by parts from ∫ ∇·(vA)·φ_i dV)
  # Full upwinding eliminates spurious Gibbs oscillations at the pulse flanks;
  # the small amount of numerical diffusion added is negligible compared to D.
  [advection]
    type              = ADConservativeAdvection
    variable          = A
    velocity_material = vel_vec    # vector material property (v_g, 0, 0)
    upwinding_type    = full       # full upwinding suppresses Gibbs oscillations
  []

  # D · ∂²A/∂x² — group-velocity dispersion (GVD) broadening term.
  # Without nonlinearity (alpha = 0) this spreads the pulse as a Gaussian;
  # the Kerr term below opposes that spreading to form a soliton at balance.
  [diffusion]
    type        = ADMatDiffusion
    variable    = A
    diffusivity = D_coeff   # scalar AD material property from [Materials]
  []

  # −α · A³ — Kerr self-phase-modulation nonlinearity (self-focusing term).
  # MatReaction residual convention: F_i = −nl_rate · A · φ_i
  # With nl_rate = −α·A² from DerivativeParsedMaterial:
  #   F_i = +α · A³ · φ_i   →   strong form adds −α·A³ on the left side.
  # This is the correct sign for self-focusing that opposes GVD broadening.
  # NOTE: MatReaction (non-AD) is used because DerivativeParsedMaterial
  # produces non-AD material properties.  The AD/non-AD property mismatch
  # causes a runtime error with ADMatReaction.
  [nonlinear]
    type          = MatReaction
    variable      = A
    reaction_rate = nl_rate   # nonlinear rate material property from [Materials]
  []
[]

[BCs]
  # A = 0 at the left inlet.  The soliton starts at x = 5 and moves rightward,
  # so the left boundary sees only the exponentially decaying sech tail (~e⁻⁵).
  # Dirichlet BC here is effectively zero throughout the run.
  [left]
    type     = ADDirichletBC
    variable = A
    boundary = left
    value    = 0
  []
  # Right boundary: zero-flux Neumann (natural BC, no block required).
  # The soliton exits the domain near t = 15 and should not reflect.
[]

[Materials]
  # Group-velocity vector supplied as an AD vector material property.
  # ADConservativeAdvection reads this via velocity_material = 'vel_vec'.
  # The y- and z-components are zero; the domain is quasi-1D in x.
  [velocity_mat]
    type        = ADGenericConstantVectorMaterial
    prop_names  = 'vel_vec'
    prop_values = '${v_g} 0 0'   # velocity vector (v_g, 0, 0)
  []

  # GVD diffusivity scalar supplied as an AD material property.
  # ADMatDiffusion reads this via diffusivity = 'D_coeff'.
  [diffusivity_mat]
    type        = ADGenericConstantMaterial
    prop_names  = 'D_coeff'
    prop_values = '${D}'   # GVD coefficient from top-level variable
  []

  # Nonlinear Kerr rate: nl_rate = −α · A²
  # Combined with MatReaction (F_i = −rate · A · φ_i) this gives:
  #   F_i = −(−α·A²) · A · φ_i = +α · A³ · φ_i
  # which is the weak form of the −α·A³ strong-form source term.
  #
  # DerivativeParsedMaterial computes nl_rate and its derivative ∂(nl_rate)/∂A
  # automatically; the derivative is used by ADMatReaction for Jacobian assembly.
  #
  # Docker JIT workaround: disable_fpoptimizer = true, enable_jit = false
  # ensures the expression is evaluated through the interpreter (not LLVM JIT),
  # which is fully portable inside containers without a C++ toolchain.
  [nl_coeff]
    type              = DerivativeParsedMaterial
    property_name     = nl_rate
    coupled_variables = 'A'
    expression        = '-${alpha}*A*A'   # = -α·A²  (nl_rate for MatReaction)
    derivative_order  = 1
    disable_fpoptimizer = true
    enable_jit          = false   # JIT requires mpicxx; disable for Docker portability
  []
[]

[Postprocessors]
  # Peak amplitude anywhere in the domain.
  # At soliton balance (alpha = 0.1) this stays near 1.0 throughout.
  # With alpha = 0 (pure GVD) it decays as 1/√(1 + t/z_d), z_d = w₀²/(2D) = 10.
  # With alpha = 0.3 (over-nonlinear) it exceeds 1 during compression phases.
  [max_A]
    type       = ElementExtremeValue
    variable   = A
    value_type = max
  []

  # Total integrated amplitude ∫A dV over the domain.
  # For a perfect soliton this is conserved: ∫sech(ξ) dξ = π (times domain width).
  # Gradual decrease indicates numerical dissipation or outflow from the right
  # boundary once the pulse centre passes x ≈ 18.
  [total_A]
    type     = ElementIntegralVariablePostprocessor
    variable = A
  []
[]

[Executioner]
  type = Transient

  # Newton with direct LU preconditioner.
  # The nonlinear Kerr term makes PJFNK less reliable for stiff transients;
  # Newton with AD Jacobians gives quadratic convergence per timestep.
  solve_type          = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  dt       = 0.02
  end_time = 10.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20

  # Adaptive timestepping: grow dt when Newton converges easily (few iterations),
  # cut it when convergence is slow.  The soliton is smooth and well-behaved so
  # growth_factor = 1.5 typically allows the step to reach 0.1 after a few steps.
  [TimeStepper]
    type               = IterationAdaptiveDT
    dt                 = 0.02
    growth_factor      = 1.5
    cutback_factor     = 0.5
    optimal_iterations = 8
  []
[]

[Outputs]
  exodus = true   # full A(x,y,t) field at every timestep for ParaView animation
  csv    = true   # max_A, total_A time-series for soliton stability verification
[]
