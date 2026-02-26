# ============================================================
# Case 42: Sod Shock Tube — 1D Riemann Problem
# Rieutord, Fluid Dynamics (Springer, 2015), Ch 5, Sec 5.5
#
# A membrane at x=0.5 separates high-pressure gas (left) from
# low-pressure gas (right). At t=0 the membrane bursts, producing
# three distinct waves: a leftward rarefaction fan, a rightward
# contact discontinuity, and a rightward shock wave.
#
# Rankine-Hugoniot jump conditions govern the shock speed and
# post-shock state. The exact solution involves all three wave
# families of the Euler equations.
#
# Left state:  rho=1.0,   p=1.0,   u=0
# Right state: rho=0.125, p=0.1,   u=0
# gamma = 1.4 (ideal diatomic gas)
#
# Validation at t=0.2:
#   Shock at x ~ 0.85
#   Contact discontinuity at x ~ 0.69
#   Rarefaction fan from x ~ 0.26 to x ~ 0.49
#   Post-shock pressure p* ~ 0.303
#
# Domain: [0, 1], 200 FV cells
# MOOSE: CNS HLLC flux splitting, ExplicitSSPRungeKutta order 2
# ============================================================

# Conservative variable initial values
# E = p / ((gamma-1) * rho) for a gas at rest (kinetic energy = 0)
# Left:  E_left  = 1.0 / (0.4 * 1.0)   = 2.5
# Right: E_right = 0.1 / (0.4 * 0.125) = 2.0
rho_left = 1.0
E_left = 2.5
u_left = 1e-15

rho_right = 0.125
E_right = 2.0
u_right = 1e-15

middle = 0.5

[GlobalParams]
  fp = fp
[]

[Mesh]
  [cartesian]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = 1
    nx = 200
  []
[]

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
  []
[]

[Variables]
  [rho]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [rho_u]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [rho_E]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[FVKernels]
  [mass_time]
    type = FVTimeKernel
    variable = rho
  []
  [mass_advection]
    type = CNSFVMassHLLC
    variable = rho
  []

  [momentum_time]
    type = FVTimeKernel
    variable = rho_u
  []
  [momentum_advection]
    type = CNSFVMomentumHLLC
    variable = rho_u
    momentum_component = x
  []

  [fluid_energy_time]
    type = FVTimeKernel
    variable = rho_E
  []
  [fluid_energy_advection]
    type = CNSFVFluidEnergyHLLC
    variable = rho_E
  []
[]

[FVBCs]
  [mass_implicit]
    type = CNSFVHLLCMassImplicitBC
    variable = rho
    fp = fp
    boundary = 'left right'
  []
  [mom_implicit]
    type = CNSFVHLLCMomentumImplicitBC
    variable = rho_u
    momentum_component = x
    fp = fp
    boundary = 'left right'
  []
  [fluid_energy_implicit]
    type = CNSFVHLLCFluidEnergyImplicitBC
    variable = rho_E
    fp = fp
    boundary = 'left right'
  []
[]

[ICs]
  [rho_ic]
    type = FunctionIC
    variable = rho
    function = 'if (x < ${middle}, ${rho_left}, ${rho_right})'
  []
  [rho_u_ic]
    type = FunctionIC
    variable = rho_u
    function = 'if (x < ${middle}, ${fparse rho_left * u_left}, ${fparse rho_right * u_right})'
  []
  [rho_E_ic]
    type = FunctionIC
    variable = rho_E
    function = 'if (x < ${middle}, ${fparse E_left * rho_left}, ${fparse E_right * rho_right})'
  []
[]

[Materials]
  [var_mat]
    type = ConservedVarValuesMaterial
    rho = rho
    rhou = rho_u
    rho_et = rho_E
    fp = fp
  []
[]

[Postprocessors]
  # Track density extremes to confirm shock and rarefaction formation
  [max_rho]
    type = ElementExtremeValue
    variable = rho
    value_type = max
  []
  [min_rho]
    type = ElementExtremeValue
    variable = rho
    value_type = min
  []
  # Track momentum to confirm wave propagation
  [max_rho_u]
    type = ElementExtremeValue
    variable = rho_u
    value_type = max
  []
  # Track total energy (should be conserved)
  [total_rho_E]
    type = ElementIntegralVariablePostprocessor
    variable = rho_E
  []
  # Total mass (should be conserved)
  [total_mass]
    type = ElementIntegralVariablePostprocessor
    variable = rho
  []
[]

[Executioner]
  type = Transient
  [TimeIntegrator]
    type = ExplicitSSPRungeKutta
    order = 2
  []
  l_tol = 1e-8

  start_time = 0.0
  dt = 5e-4
  end_time = 0.2
  abort_on_solve_fail = true
[]

[Outputs]
  exodus = true
  csv = true
[]
