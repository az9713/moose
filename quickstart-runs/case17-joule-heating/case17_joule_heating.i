# ============================================================
# Case 17: Joule Heating — Electric Current Generates Heat
# Coupled electrostatic + thermal problem:
#   -div(sigma * grad(V)) = 0          (Laplace equation for V)
#   rho*cp * dT/dt = div(k*grad(T)) + sigma*|grad(V)|^2   (heat equation with Joule source)
#
# Domain: 2D rectangle, x in [0,2], y in [0,1]
# V = 10 V on left, V = 0 V on right
# T = 300 K on left and right (electrodes as heat sinks)
# Top and bottom: insulated (natural Neumann, zero flux)
# sigma = 1e6 S/m, k = 50 W/(m K), rho = 8000 kg/m^3, cp = 500 J/(kg K)
# Solved with combined-opt (heat_transfer module)
# ============================================================

[Mesh]
  [gen]
    type  = GeneratedMeshGenerator
    dim   = 2
    nx    = 40
    ny    = 20
    xmin  = 0
    xmax  = 2
    ymin  = 0
    ymax  = 1
  []
[]

[Variables]
  # Electric potential [V].  Solved first (Laplace equation, no time dependence).
  [V]
    initial_condition = 0.0
  []

  # Temperature [K].  Driven transient by Joule heat source.
  [T]
    initial_condition = 300.0
  []
[]

[Kernels]
  # -------------------------------------------------------
  # Electric potential equation: -div(sigma * grad(V)) = 0
  # ADHeatConduction reads the material property named
  # 'thermal_conductivity' by default.  We alias that
  # property to 'electrical_conductivity' via the parameter.
  # -------------------------------------------------------
  [V_diff]
    type                 = ADHeatConduction
    variable             = V
    thermal_conductivity = electrical_conductivity
  []

  # -------------------------------------------------------
  # Heat equation: rho*cp * dT/dt = div(k*grad(T)) + Q_joule
  # -------------------------------------------------------
  [T_time]
    type     = ADHeatConductionTimeDerivative
    variable = T
  []

  [T_diff]
    type     = ADHeatConduction
    variable = T
  []

  # Joule heating source: Q = sigma * |grad(V)|^2
  # ADJouleHeatingSource reads the Joule heating value from
  # the 'electric_field_heating' material property computed
  # by ElectromagneticHeatingMaterial below.
  [T_joule]
    type         = ADJouleHeatingSource
    variable     = T
    heating_term = electric_field_heating
  []
[]

[BCs]
  # Electric potential: 10 V at left electrode, 0 V at right electrode.
  # Top and bottom have no BC — natural zero-flux Neumann (insulated).
  [V_left]
    type     = ADDirichletBC
    variable = V
    boundary = left
    value    = 10.0
  []
  [V_right]
    type     = ADDirichletBC
    variable = V
    boundary = right
    value    = 0.0
  []

  # Temperature: 300 K at both electrodes (act as heat sinks).
  # Top and bottom have no BC — natural zero-flux Neumann (insulated walls).
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
  # Thermal properties (used by ADHeatConduction and ADHeatConductionTimeDerivative)
  [thermal]
    type        = ADGenericConstantMaterial
    prop_names  = 'thermal_conductivity specific_heat density'
    prop_values = '50.0              500.0        8000.0'
  []

  # Electrical conductivity [S/m]
  [electrical]
    type        = ADGenericConstantMaterial
    prop_names  = 'electrical_conductivity'
    prop_values = '1e6'
  []

  # ElectromagneticHeatingMaterial couples V into a material property
  # 'electric_field_heating' = sigma * |grad(V)|^2.
  # The formulation is 'time' (not frequency domain) and the solver
  # is 'electrostatic' (scalar potential, not a vector EM field).
  [joule_material]
    type                       = ElectromagneticHeatingMaterial
    electric_field             = V
    electric_field_heating_name = electric_field_heating
    electrical_conductivity    = electrical_conductivity
    formulation                = time
    solver                     = electrostatic
  []
[]

[Postprocessors]
  # Track maximum temperature to see how fast Joule heating warms the conductor.
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []

  # Domain-average temperature gives a single number for each timestep.
  [avg_T]
    type     = ElementAverageValue
    variable = T
  []

  # Verify the electric potential field is steady and correct.
  [max_V]
    type       = ElementExtremeValue
    variable   = V
    value_type = max
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # Timestep: 0.25 s, run for 5 s total (20 steps).
  dt       = 0.25
  end_time = 5.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20
[]

# Full single-matrix preconditioning to capture V-T coupling in the Jacobian.
[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Outputs]
  exodus = true   # Full spatial fields (V and T) at every timestep
  csv    = true   # Postprocessor history (max_T, avg_T, max_V vs. time)
[]
