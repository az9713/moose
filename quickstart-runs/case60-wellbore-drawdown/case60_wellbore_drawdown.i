# ============================================================
# Case 60: Wellbore Drawdown — Radial Flow to a Production Well
# Simulates pressure drawdown around a pumping well in a
# confined aquifer using axisymmetric (RZ) coordinates.
#
# Theis solution (1935) for radial flow:
#   p(r,t) = p0 - Q*mu/(4*pi*k*H) * W(r^2*S/(4*k*t/mu))
# where W(u) is the well function (exponential integral).
#
# Parameters:
#   r_well = 0.1 m (wellbore radius)
#   r_outer = 100 m (outer boundary)
#   H = 10 m (aquifer thickness, modeled as unit thickness in RZ)
#   k = 1e-12 m^2 (~1 Darcy)
#   mu = 0.001 Pa.s
#   porosity = 0.2
#   Initial pressure p0 = 10 MPa (1000 m depth)
#   Well flow rate Q = 0.001 m^3/s (via pressure BC)
#   Drawdown ~ Q*mu/(4*pi*k*H) * ln(r_outer/r_well) ~ 70 kPa
#
# Domain: r in [0.1, 100] m (log-spaced would be ideal but
#   GeneratedMesh with bias_x gives reasonable refinement)
# coord_type = RZ (x = r, y = z)
# Run for 10000 s to approach steady-state drawdown cone
# ============================================================

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 50
    ny = 2
    xmin = 0.1       # wellbore radius
    xmax = 100.0     # outer boundary radius
    ymin = 0
    ymax = 1.0       # unit thickness
    bias_x = 1.15    # refine near wellbore
  []
  coord_type = RZ     # x = r, y = z
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
    initial_condition = 1e7   # 10 MPa initial pressure
  []
[]

[PorousFlowBasicTHM]
  porepressure    = porepressure
  coupling_type   = Hydro
  gravity         = '0 0 0'     # horizontal aquifer, ignore gravity
  fp              = water
  multiply_by_density = true
[]

[FluidProperties]
  [water]
    type                 = SimpleFluidProperties
    density0             = 1000
    viscosity            = 0.001
    thermal_expansion    = 0
    cp                   = 4186
    cv                   = 4186
    thermal_conductivity = 0.6
  []
[]

[Materials]
  [porosity]
    type          = PorousFlowPorosity
    porosity_zero = 0.2
    mechanical    = false
    thermal       = false
    fluid         = false
  []
  [biot_modulus]
    type                  = PorousFlowConstantBiotModulus
    biot_coefficient      = 1.0
    solid_bulk_compliance = 1e-10
    fluid_bulk_modulus    = 2e9
  []
  [thermal_expansion]
    type                 = PorousFlowConstantThermalExpansionCoefficient
    biot_coefficient     = 1.0
    drained_coefficient  = 0.0
    fluid_coefficient    = 0.0
  []
  [permeability]
    type         = PorousFlowPermeabilityConst
    permeability = '1e-12 0 0  0 1e-12 0  0 0 1e-12'
  []
  # Note: thermal materials removed — Hydro-only coupling
[]

[BCs]
  # Well boundary (left = inner radius): reduced pressure (pumping)
  [well_pressure]
    type     = DirichletBC
    variable = porepressure
    boundary = left
    value    = 9.9e6     # 9.9 MPa (100 kPa drawdown at well)
  []
  # Outer boundary: constant pressure (infinite aquifer approx)
  [outer_pressure]
    type     = DirichletBC
    variable = porepressure
    boundary = right
    value    = 1e7       # 10 MPa (undisturbed)
  []
[]

[Postprocessors]
  # Pressure at wellbore
  [p_well]
    type     = PointValue
    variable = porepressure
    point    = '0.1 0.5 0'
  []
  # Pressure at 1 m
  [p_1m]
    type     = PointValue
    variable = porepressure
    point    = '1.0 0.5 0'
  []
  # Pressure at 10 m
  [p_10m]
    type     = PointValue
    variable = porepressure
    point    = '10.0 0.5 0'
  []
  # Pressure at 50 m
  [p_50m]
    type     = PointValue
    variable = porepressure
    point    = '50.0 0.5 0'
  []
  # Average pressure
  [p_avg]
    type     = ElementAverageValue
    variable = porepressure
  []
[]

[Executioner]
  type       = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  dt       = 10        # 10 s time steps (capture early transient)
  end_time = 500       # 500 s total

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
