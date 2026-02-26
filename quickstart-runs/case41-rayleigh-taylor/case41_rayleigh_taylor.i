# ============================================================
# Case 41: Rayleigh-Taylor Instability
# Rieutord, Fluid Dynamics (Springer, 2015), Ch 6, Sec 6.3
#
# Heavy fluid atop light fluid in a gravitational field.
# A small sinusoidal perturbation at the density interface
# grows into mushroom-shaped fingers as the heavy fluid
# sinks through the light fluid.
#
# Boussinesq approximation with T as density marker:
#   rho_eff = rho * (1 - alpha*(T - T_ref))
# With alpha=1, T_ref=0.5, g=(0,-1,0):
#   T < T_ref → heavier (sinks)
#   T > T_ref → lighter (rises)
#
# IC: T = 0.5*(1 - tanh((y-1)/0.02))
#   y >> 1: T → 0 (cold, heavy on top → UNSTABLE)
#   y << 1: T → 1 (hot, light on bottom)
# Perturbation: + 0.05*cos(2*pi*x)*exp(-(y-1)^2/0.01)
#
# Parameters: rho=1, mu=0.001 (Re~1000), kappa=0.001
# Domain: [0, 1] x [0, 2] (tall), 25x50
#
# Linear growth rate: sigma ~ sqrt(A*g*k)
#   A = Atwood ~ 0.5, k = 2*pi
#   sigma ~ sqrt(0.5 * 1 * 2*pi) ~ 1.77
# ============================================================

mu_val = 0.001
k_val = 0.001

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 2
    nx = 25
    ny = 50
  []
[]

[Modules]
  [NavierStokesFV]
    compressibility = 'incompressible'
    porous_medium_treatment = false
    add_energy_equation = true

    density = 'rho'
    dynamic_viscosity = 'mu'
    thermal_conductivity = 'k'
    specific_heat = 'cp'

    initial_velocity = '1e-15 1e-15 0'
    initial_pressure = 0.0
    initial_temperature = T_init_func

    # All four walls: no-slip + insulated
    wall_boundaries = 'left right top bottom'
    momentum_wall_types = 'noslip noslip noslip noslip'
    energy_wall_types = 'heatflux heatflux heatflux heatflux'
    energy_wall_functors = '0 0 0 0'

    # Boussinesq buoyancy
    boussinesq_approximation = true
    gravity = '0 -1 0'
    ref_temperature = 0.5
    thermal_expansion = 'alpha'

    # Pin pressure (closed box)
    pin_pressure = true
    pinned_pressure_type = average
    pinned_pressure_value = 0

    mass_advection_interpolation = 'average'
    momentum_advection_interpolation = 'upwind'
    energy_advection_interpolation = 'upwind'
  []
[]

[FunctorMaterials]
  [const]
    type = ADGenericFunctorMaterial
    prop_names = 'rho mu k cp alpha'
    prop_values = '1.0 ${mu_val} ${k_val} 1.0 1.0'
  []
[]

# IC: heavy fluid (T=0) on top, light fluid (T=1) on bottom
# Smooth interface via tanh + sinusoidal perturbation
[Functions]
  [T_init_func]
    type = ParsedFunction
    expression = '0.5*(1.0 - tanh((y-1.0)/0.02)) + 0.05*cos(2*pi*x)*exp(-((y-1.0)^2)/0.01)'
  []
[]

[Postprocessors]
  [max_vel_x]
    type = ADElementExtremeFunctorValue
    functor = vel_x
  []
  [max_vel_y]
    type = ADElementExtremeFunctorValue
    functor = vel_y
  []
  [min_vel_y]
    type = ADElementExtremeFunctorValue
    functor = vel_y
    value_type = min
  []
  [avg_T]
    type = ElementAverageValue
    variable = T_fluid
  []
  [max_T]
    type = ElementExtremeValue
    variable = T_fluid
    value_type = max
  []
  [min_T]
    type = ElementExtremeValue
    variable = T_fluid
    value_type = min
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20

  dt = 0.05
  end_time = 3.0

  automatic_scaling = true

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    growth_factor = 1.2
    cutback_factor = 0.5
    optimal_iterations = 8
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
