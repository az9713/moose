# ============================================================
# Case 37: Rayleigh-Benard Convection Onset
# Rieutord, Fluid Dynamics (Springer, 2015), Ch 7, Sec 7.5
#
# Fluid heated from below between two rigid plates. Above the
# critical Rayleigh number Ra_c = 1708, the conduction state
# becomes unstable and convective rolls form. The Nusselt
# number Nu > 1 indicates convective heat transport exceeding
# pure conduction.
#
# Governing equations (Boussinesq incompressible NS + energy):
#   div(u) = 0
#   rho*(du/dt + u.grad(u)) = -grad(p) + mu*lap(u) - rho*alpha*(T-T_ref)*g
#   dT/dt + u.grad(T) = kappa*lap(T)
#
# Non-dimensional parameters:
#   Ra = g*alpha*dT*H^3/(nu*kappa) = 2000 (above Ra_c=1708)
#   Pr = nu/kappa = 0.71 (air)
#
# Material properties (non-dimensional, rho=1, cp=1):
#   nu = mu = sqrt(Pr/Ra) = sqrt(0.71/2000) = 0.01884
#   kappa = k = nu/Pr = 0.01884/0.71 = 0.02653
#   alpha = 1.0, g = (0,-1,0), T_ref = 0.5
#
# Domain: [0, 2] x [0, 1], 40x20 (aspect ratio 2:1)
# BCs: bottom T=1 (hot), top T=0 (cold), sides insulated
#       all walls no-slip
# IC: T = 1-y + 0.01*sin(pi*x)*sin(pi*y) perturbation
#
# Differentiation from Case 16:
#   Case 16: hot LEFT wall, Ra=10000, steady state
#   Case 37: hot BOTTOM, Ra=2000 near onset, transient dynamics
#
# Validation: Nu > 1 at steady state, convective rolls visible
# ============================================================

nu    = 0.01884
kappa = 0.02653

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 1
    nx = 40
    ny = 20
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

    # All four walls are no-slip
    wall_boundaries = 'left right top bottom'
    momentum_wall_types = 'noslip noslip noslip noslip'

    # Energy BCs (order matches wall_boundaries: left right top bottom):
    #   left:   heatflux 0  (insulated side wall)
    #   right:  heatflux 0  (insulated side wall)
    #   top:    fixed-temperature 0  (cold plate)
    #   bottom: fixed-temperature 1  (hot plate)
    energy_wall_types = 'heatflux heatflux fixed-temperature fixed-temperature'
    energy_wall_functors = '0 0 0 1'

    # Boussinesq buoyancy
    boussinesq_approximation = true
    gravity = '0 -1 0'
    ref_temperature = 0.5
    thermal_expansion = 'alpha'

    # Pin pressure (closed cavity)
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
    prop_values = '1.0 ${nu} ${kappa} 1.0 1.0'
  []
[]

# Seed the instability with a sinusoidal temperature perturbation
# T = 1 - y + 0.01*sin(pi*x)*sin(pi*y)
[Functions]
  [T_init_func]
    type = ParsedFunction
    expression = '1 - y + 0.01*sin(pi*x)*sin(pi*y)'
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
  [avg_T]
    type = ElementAverageValue
    variable = T_fluid
  []
  # Nusselt number proxy: average T on hot wall
  # For pure conduction Nu=1, for convection Nu>1
  [T_hot_wall]
    type = SideAverageFunctorPostprocessor
    functor = T_fluid
    boundary = bottom
  []
  [T_cold_wall]
    type = SideAverageFunctorPostprocessor
    functor = T_fluid
    boundary = top
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

  dt = 0.5
  end_time = 20.0

  automatic_scaling = true

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    growth_factor = 1.5
    cutback_factor = 0.5
    optimal_iterations = 8
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
