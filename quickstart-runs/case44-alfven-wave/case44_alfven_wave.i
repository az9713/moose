# ============================================================
# Case 44: Alfven Wave Propagation — MHD Elsasser Variables
# Rieutord, Fluid Dynamics (Springer, 2015), Ch 10, Sec 10.4
#
# Transverse MHD wave in a conducting fluid with background
# field B_0. Using Elsasser variables:
#   d+ = vy + by/sqrt(mu_0*rho)
#   d- = vy - by/sqrt(mu_0*rho)
# the MHD equations decouple into two advection-diffusion
# equations with wave speeds +/- v_A (Alfven speed):
#
#   d(d+)/dt + v_A * d(d+)/dx = D_eff * d²(d+)/dx²
#   d(d-)/dt - v_A * d(d-)/dx = D_eff * d²(d-)/dx²
#
# A Gaussian initial pulse in d+ propagates rightward at v_A
# while d- remains zero (no leftward-propagating wave).
#
# Parameters: v_A=1.0, D_eff=0.01
# IC: d+(x,0) = exp(-(x-3)^2/0.25), d-(x,0) = 0
# Gaussian half-width w = 0.5
#
# Validation at t=5:
#   d+ peak at x = 3 + v_A*t = 8
#   Peak amplitude decays as 1/sqrt(1 + 4*D*t/w^2)
#     = 1/sqrt(1 + 4*0.01*5/0.25) = 1/sqrt(1.8) ~ 0.745
#   d- remains ~ 0 everywhere
#
# Domain: [0, 12] x [0, 0.12], 120x2 quad elements (quasi-1D)
# ============================================================

v_A = 1.0
D_eff = 0.01

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 12
  ymin = 0
  ymax = 0.12
  nx = 120
  ny = 2
[]

[Variables]
  [d_plus]
  []
  [d_minus]
  []
[]

[ICs]
  # Gaussian pulse in d+, centered at x=3
  [d_plus_ic]
    type = FunctionIC
    variable = d_plus
    function = 'exp(-((x-3)^2)/0.25)'
  []
  # d- starts at zero (no leftward wave)
  [d_minus_ic]
    type = FunctionIC
    variable = d_minus
    function = '0'
  []
[]

[Kernels]
  # --- d+ equation: d(d+)/dt + v_A*d(d+)/dx = D*d²(d+)/dx² ---
  [dp_time]
    type = ADTimeDerivative
    variable = d_plus
  []
  [dp_advection]
    type = ADConservativeAdvection
    variable = d_plus
    velocity_material = vel_plus
    upwinding_type = full
  []
  [dp_diffusion]
    type = ADMatDiffusion
    variable = d_plus
    diffusivity = D_coeff
  []

  # --- d- equation: d(d-)/dt - v_A*d(d-)/dx = D*d²(d-)/dx² ---
  [dm_time]
    type = ADTimeDerivative
    variable = d_minus
  []
  [dm_advection]
    type = ADConservativeAdvection
    variable = d_minus
    velocity_material = vel_minus
    upwinding_type = full
  []
  [dm_diffusion]
    type = ADMatDiffusion
    variable = d_minus
    diffusivity = D_coeff
  []
[]

[BCs]
  # d+ = 0 at both boundaries (pulse is well within domain)
  [dp_left]
    type = ADDirichletBC
    variable = d_plus
    boundary = left
    value = 0
  []
  [dp_right]
    type = ADDirichletBC
    variable = d_plus
    boundary = right
    value = 0
  []
  # d- = 0 at both boundaries
  [dm_left]
    type = ADDirichletBC
    variable = d_minus
    boundary = left
    value = 0
  []
  [dm_right]
    type = ADDirichletBC
    variable = d_minus
    boundary = right
    value = 0
  []
[]

[Materials]
  # Rightward Alfven velocity for d+
  [vel_plus_mat]
    type = ADGenericConstantVectorMaterial
    prop_names = 'vel_plus'
    prop_values = '${v_A} 0 0'
  []
  # Leftward Alfven velocity for d- (negative v_A)
  [vel_minus_mat]
    type = ADGenericConstantVectorMaterial
    prop_names = 'vel_minus'
    prop_values = '${fparse -v_A} 0 0'
  []
  # Diffusivity (resistive + viscous dissipation combined)
  [diff_mat]
    type = ADGenericConstantMaterial
    prop_names = 'D_coeff'
    prop_values = '${D_eff}'
  []
[]

[Postprocessors]
  # Peak amplitude of d+ (should decay as Gaussian broadens)
  [max_d_plus]
    type = ElementExtremeValue
    variable = d_plus
    value_type = max
  []
  # Peak of d- (should stay near zero)
  [max_d_minus]
    type = ElementExtremeValue
    variable = d_minus
    value_type = max
  []
  # Total energy proxy: integral of d+^2 + d-^2
  [total_d_plus]
    type = ElementIntegralVariablePostprocessor
    variable = d_plus
  []
  [total_d_minus]
    type = ElementIntegralVariablePostprocessor
    variable = d_minus
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  dt = 0.05
  end_time = 6.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 15

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 6
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
