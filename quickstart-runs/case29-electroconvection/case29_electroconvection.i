# ============================================================
# Case 29: Electroconvection — EHD-Enhanced Natural Convection
# Melcher, Continuum Electromechanics (MIT, 1981), Ch. 9 §9.12
#
# Natural convection in a differentially heated cavity (Case 16)
# with an additional electric body force that can enhance or
# suppress the buoyancy-driven flow:
#
#   ρ(v·∇)v = -∇p + μ∇²v + ρα(T-T_ref)g + f_EHD
#   f_EHD_y = Fe · (T_fluid - 0.5)
#
# The EHD force mimics dielectrophoresis in a non-uniform field
# where permittivity depends on temperature.
#
# Fe = 0:   pure natural convection (reproduces Case 16)
# Fe > 0:   EHD enhances convection  (stronger circulation)
# Fe < 0:   EHD suppresses convection (weaker circulation)
#
# Domain: [0,1]², 25×25 mesh
# BCs: same as Case 16 (hot left, cold right, insulated top/bottom)
# Ra = 10000, Pr = 0.71
# Fe = 5.0 (enhancement mode)
# ============================================================

# Non-dimensional kinematic viscosity and thermal diffusivity
# Pr = nu/kappa = 0.71,   Ra = g*alpha*dT*L^3/(nu*kappa) = 10000
# nu = sqrt(Pr/Ra) = sqrt(0.71/10000) = 0.008426
# kappa = nu/Pr = 0.008426/0.71 = 0.011867
# Check: Ra = 1*1*1*1/(0.008426*0.011867) = 1/0.0001 = 10000  (approx)
nu    = 0.008426   # = mu  since rho = 1
kappa = 0.011867   # = k   since rho*cp = 1

# EHD force parameter
# Fe = 0:   pure natural convection (reproduces Case 16, Nu ~ 2.24)
# Fe > 0:   EHD force reinforces buoyancy near hot wall (T > 0.5),
#           enhancing the circulation and raising Nu above 2.24
# Fe < 0:   EHD force opposes buoyancy, suppressing circulation
# The EHD force f_y = Fe*(T-0.5) has the same form as Boussinesq buoyancy
# f_y = alpha*(T-T_ref)*|g|.  By combining them into an effective thermal
# expansion coefficient alpha_eff = alpha_buoy + Fe = 1.0 + 5.0 = 6.0,
# we capture both effects within the built-in Boussinesq mechanism.
# Fe = 5.0 (enhancement mode)
alpha_eff = 6.0    # = 1.0 (buoyancy) + 5.0 (EHD enhancement)

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx   = 25
    ny   = 25
  []
[]

# Note: [Modules/NavierStokesFV] is a deprecated but still-functional syntax.
# The modern equivalent is [Physics/NavierStokes/Flow/...] and
# [Physics/NavierStokes/FluidHeatTransfer/...].
# This block constructs all FV variables (vel_x, vel_y, pressure, T_fluid),
# all momentum/continuity/energy kernels, Rhie-Chow interpolation, and all
# wall BCs automatically.  The EHD body force is appended separately below
# in [FVKernels].
[Modules]
  [NavierStokesFV]
    compressibility          = 'incompressible'
    porous_medium_treatment  = false
    add_energy_equation      = true

    # Material property references (functor names from [FunctorMaterials])
    density              = 'rho'
    dynamic_viscosity    = 'mu'
    thermal_conductivity = 'k'
    specific_heat        = 'cp'

    # Initial conditions: fluid at rest, mid-range temperature
    initial_velocity    = '1e-15 1e-15 0'
    initial_pressure    = 0.0
    initial_temperature = 0.5

    # All four walls are no-slip for velocity
    wall_boundaries       = 'left right top bottom'
    momentum_wall_types   = 'noslip noslip noslip noslip'

    # Energy BCs:
    #   left  -> fixed T = 1 (hot wall)
    #   right -> fixed T = 0 (cold wall)
    #   top   -> zero heat flux (insulated)
    #   bottom-> zero heat flux (insulated)
    energy_wall_types    = 'fixed-temperature fixed-temperature heatflux heatflux'
    energy_wall_functors = '1 0 0 0'

    # Boussinesq buoyancy: rho = rho_ref * (1 - alpha*(T - T_ref))
    boussinesq_approximation = true
    gravity                  = '0 -1 0'
    ref_temperature          = 0.5
    thermal_expansion        = 'alpha'

    # Pin the pressure to zero average (closed cavity with no outlets)
    pin_pressure          = true
    pinned_pressure_type  = average
    pinned_pressure_value = 0

    # Upwind for momentum and energy: more stable than central differencing
    # when convection dominates in the thermal boundary layers
    mass_advection_interpolation     = 'average'
    momentum_advection_interpolation = 'upwind'
    energy_advection_interpolation   = 'upwind'
  []
[]

[FunctorMaterials]
  # Fluid properties (non-dimensional): rho = 1, cp = 1, so
  #   nu = mu/rho = mu  and  kappa = k/(rho*cp) = k
  # Values are derived from Ra = 10000, Pr = 0.71 (see header).
  #
  # The effective thermal expansion coefficient alpha_eff combines
  # natural buoyancy (α = 1.0) with the EHD force (Fe = 5.0):
  #   f_total_y = alpha_eff * (T - T_ref) * |g| = 6.0 * (T - 0.5)
  # This exploits the mathematical equivalence between Boussinesq
  # buoyancy and the dielectrophoretic EHD body force.
  [const]
    type        = ADGenericFunctorMaterial
    prop_names  = 'rho  mu       k         cp  alpha'
    prop_values = '1.0  ${nu}   ${kappa}  1.0  ${alpha_eff}'
  []
[]

[Postprocessors]
  # Maximum velocity components: larger magnitudes than Case 16 (Fe=0)
  # confirm that EHD is augmenting the buoyancy-driven circulation.
  [max_vel_x]
    type    = ADElementExtremeFunctorValue
    functor = vel_x
  []
  [max_vel_y]
    type    = ADElementExtremeFunctorValue
    functor = vel_y
  []

  # Average temperature: should remain near 0.5 by symmetry regardless of Fe.
  [avg_T]
    type     = ElementAverageValue
    variable = T_fluid
  []

  # Side average of T on the hot wall (left boundary).
  # With the fixed-temperature BC, this reports T = 1 and is mainly a
  # diagnostic sanity check that the wall BC is applied correctly.
  [Nu_hot_wall]
    type     = SideAverageFunctorPostprocessor
    functor  = 'T_fluid'
    boundary = 'left'
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'NEWTON'

  # Direct LU factorization: robust for the coupled flow-heat-EHD system.
  # NONZERO shift prevents zero-pivot failures in the pressure null space.
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  # Tight nonlinear tolerance to capture the correct EHD-modified flow structure.
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30

  # Automatic scaling balances velocity, pressure, and temperature DOFs,
  # which can differ by orders of magnitude in the coupled system.
  automatic_scaling = true
[]

[Outputs]
  exodus = true   # Full field output (vel_x, vel_y, pressure, T_fluid) for ParaView
  csv    = true   # Postprocessor values at convergence
[]
