# ============================================================
# Case 25: Induction Heating — Magnetic Diffusion + Joule Heat
# Melcher, Continuum Electromechanics (MIT, 1981), Ch. 6 §6.7-6.8
#
# An oscillating magnetic field diffuses into a conducting slab.
# Eddy currents induced by the changing field generate Joule
# heating concentrated within the electromagnetic skin depth:
#
#   ∂B/∂t = D_m · ∇²B                    (magnetic diffusion)
#   ρcp·∂T/∂t = k·∇²T + Q_eddy          (heat equation)
#   Q_eddy = Q_coeff · (∂B/∂x)²          (eddy-current heating)
#
# The heating coefficient absorbs the physical constants:
#   Q_coeff = D_m / μ₀  (non-dimensionalised here as 0.5)
#
# Skin depth:  δ ~ √(2·D_m / ω)
# With D_m = 0.005, ω = 2π/τ = 4π:  δ ≈ √(0.01/4π) ≈ 0.028
# Most heating is therefore confined to x ≲ 0.03.
#
# Domain: [0,1] × [0,0.04] quasi-1D slab, 50×2 mesh
# IC:  B = 0,  T = 300 K
# BC:  B(0,t) = sin(2π·t/τ),  B(1,t) = 0  (far-field)
#      T = 300 K at both ends (thermal sinks)
# ============================================================

D_m     = 0.005   # magnetic diffusivity D_m = 1/(μ₀·σ) [m²/s, normalised]
k_th    = 1.0     # thermal conductivity k [W/(m·K)]
rho_cp  = 1.0     # volumetric heat capacity ρ·cp [J/(m³·K)]
Q_coeff = 0.5     # eddy-heating coefficient D_m/μ₀ (non-dimensional)
tau     = 0.5     # period of the applied oscillating field [s]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    # 50 elements along x resolves the skin depth (δ ≈ 0.028 → ~1.4 elements/δ).
    # 2 elements along y keeps the problem quasi-1D while satisfying MOOSE's 2D requirement.
    nx   = 50
    ny   = 2
    xmin = 0
    xmax = 1       # slab thickness (normalised)
    ymin = 0
    ymax = 0.04    # thin strip; physics is y-independent
  []
[]

[Variables]
  # Magnetic flux density B [T, normalised].
  # IC = 0: the conductor is field-free before the oscillating source is switched on.
  [B]
    initial_condition = 0.0
  []

  # Temperature T [K].
  # IC = 300 K: uniform reference temperature at t = 0.
  [T]
    initial_condition = 300.0
  []
[]

[AuxVariables]
  # x-component of the magnetic field gradient, ∂B/∂x.
  # Stored as an element-average (CONSTANT MONOMIAL) because
  # VariableGradientComponent computes gradients at quadrature points.
  # Physically this quantity is proportional to the eddy-current density J = (1/μ₀)·∂B/∂x.
  [dBdx]
    order  = CONSTANT
    family = MONOMIAL
  []

  # Eddy-current heating power density Q_eddy = Q_coeff · (∂B/∂x)².
  # Also element-average (CONSTANT MONOMIAL) to be consistent with dBdx.
  # This auxiliary variable is the source term for the heat equation.
  [eddy_heat]
    order  = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  # Compute ∂B/∂x at element centres.
  # VariableGradientComponent extracts one spatial component of the gradient
  # of a LAGRANGE variable and stores it in a MONOMIAL AuxVariable.
  [calc_dBdx]
    type              = VariableGradientComponent
    variable          = dBdx
    gradient_variable = B
    component         = x
    execute_on        = 'TIMESTEP_END'
  []

  # Compute Q_eddy = Q_coeff · dBdx² from the stored gradient.
  # ParsedAux evaluates an arbitrary algebraic expression involving coupled
  # AuxVariables, avoiding the need for a custom C++ material or kernel.
  [calc_eddy_heat]
    type              = ParsedAux
    variable          = eddy_heat
    coupled_variables = 'dBdx'
    expression        = '${Q_coeff} * dBdx * dBdx'
    execute_on        = 'TIMESTEP_END'
  []
[]

[Kernels]
  # ------------------------------------------------------------------
  # Magnetic diffusion equation:  ∂B/∂t = D_m · ∇²B
  # ------------------------------------------------------------------

  # ∂B/∂t — accumulation of magnetic flux density in the conductor.
  [B_time]
    type     = ADTimeDerivative
    variable = B
  []

  # D_m · ∇²B — magnetic diffusion driven by field gradients.
  # ADMatDiffusion reads the diffusivity from the material property named
  # by the 'diffusivity' parameter (here: 'mag_diffusivity').
  [B_diff]
    type        = ADMatDiffusion
    variable    = B
    diffusivity = mag_diffusivity
  []

  # ------------------------------------------------------------------
  # Heat equation:  ρcp·∂T/∂t = k·∇²T + Q_eddy
  # ------------------------------------------------------------------

  # ρcp·∂T/∂t — thermal inertia.
  # ADHeatConductionTimeDerivative reads 'density' and 'specific_heat'
  # material properties and forms ρ·cp·∂T/∂t automatically.
  [T_time]
    type     = ADHeatConductionTimeDerivative
    variable = T
  []

  # k·∇²T — Fourier heat conduction.
  # ADHeatConduction reads the 'thermal_conductivity' material property.
  [T_diff]
    type     = ADHeatConduction
    variable = T
  []

  # Q_eddy source — eddy-current heating drives the temperature field.
  # CoupledForce adds +∫ v·ψ_i dV to the T residual, where v = eddy_heat.
  # Because eddy_heat is updated at TIMESTEP_END and kernels are evaluated
  # at the start of the next step, this coupling is lagged by one timestep.
  # For the small dt used here the lag error is negligible.
  [T_source]
    type     = CoupledForce
    variable = T
    v        = eddy_heat
  []
[]

[BCs]
  # Left boundary (x = 0): oscillating applied field.
  # The function sin(2π·t/τ) represents the surface value of B imposed by an
  # external coil.  The amplitude is 1 (normalised).
  [B_left]
    type     = ADFunctionDirichletBC
    variable = B
    boundary = left
    function = 'sin(2*pi*t/${tau})'
  []

  # Right boundary (x = 1): far-field condition B = 0.
  # The slab is thick enough that the oscillating field does not penetrate
  # to the right edge (skin depth δ ≈ 0.028 ≪ 1).
  [B_right]
    type     = ADDirichletBC
    variable = B
    boundary = right
    value    = 0.0
  []

  # Thermal BCs: both ends held at 300 K (electrode / heat-sink surfaces).
  [T_left]
    type     = ADDirichletBC
    variable = T
    boundary = left
    value    = 300.0
  []
  [T_right]
    type     = ADDirichletBC
    variable = T
    boundary = right
    value    = 300.0
  []
[]

[Materials]
  # Magnetic diffusivity D_m = 1/(μ₀·σ).
  # The property name 'mag_diffusivity' is passed to ADMatDiffusion via
  # the 'diffusivity' parameter.
  [mag_diff]
    type        = ADGenericConstantMaterial
    prop_names  = 'mag_diffusivity'
    prop_values = '${D_m}'
  []

  # Thermal properties used by ADHeatConduction and ADHeatConductionTimeDerivative.
  # 'specific_heat' and 'density' are combined as ρ·cp = 1·1 = 1.0.
  [thermal]
    type        = ADGenericConstantMaterial
    prop_names  = 'thermal_conductivity specific_heat density'
    prop_values = '${k_th}            1.0          ${rho_cp}'
  []
[]

[Postprocessors]
  # Maximum temperature — rises as eddy-current heating accumulates near
  # the left surface.  Should exceed 300 K after the first few oscillations.
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []

  # Domain-average temperature — rises more slowly; indicates total power deposition.
  [avg_T]
    type     = ElementAverageValue
    variable = T
  []

  # Maximum B anywhere in the domain — confirms the oscillating BC is active.
  [max_B]
    type       = ElementExtremeValue
    variable   = B
    value_type = max
  []

  # Peak eddy-current heating power density — localised near the skin depth.
  [max_eddy_heat]
    type       = ElementExtremeValue
    variable   = eddy_heat
    value_type = max
  []
[]

[Executioner]
  type = Transient

  # NEWTON with exact AD Jacobians.  Even though the B–T coupling is lagged
  # through the AuxVariable, the within-timestep convergence for each field
  # is quadratic.
  solve_type          = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # dt = 0.01 s gives 50 steps per oscillation period τ = 0.5 s.
  # 500 steps to t = 5 s = 10 oscillation periods.
  dt       = 0.01
  end_time = 5.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20
[]

# Full single-matrix preconditioning captures the off-diagonal coupling
# between B and T in the Jacobian for robust convergence.
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Outputs]
  exodus = true   # full B, T, dBdx, eddy_heat fields at every timestep
  csv    = true   # max_T, avg_T, max_B, max_eddy_heat vs. time
[]
