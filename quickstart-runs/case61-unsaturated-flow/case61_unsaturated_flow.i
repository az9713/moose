# ============================================================
# Case 61: Unsaturated Flow — Wetting Front in Soil Column
# Simulates infiltration of water into an initially dry soil
# column using Richards' equation (single-phase unsaturated).
#
# Physics: Richards' equation for variably saturated flow
#   d(phi*S)/dt + div(rho*k*kr/mu * (grad(P) - rho*g)) = 0
# where S = S(P) is the saturation-pressure relationship
# and kr = kr(S) is the relative permeability.
#
# Van Genuchten-Mualem model:
#   S_eff = (1 + (alpha*|P|)^n)^(-m), m = 1 - 1/n
#   kr = sqrt(S_eff) * (1 - (1 - S_eff^(1/m))^m)^2
#
# Parameters (typical sandy soil):
#   k = 1e-11 m^2 (saturated permeability)
#   porosity = 0.35
#   alpha_vG = 3.5e-4 /Pa (Van Genuchten alpha)
#   n_vG = 2.0 (Van Genuchten n)
#   Initial pressure = -50 kPa (dry, below atmospheric)
#   Top BC: P = 0 (saturated, water ponded on surface)
#   Bottom: free drainage (zero flux)
#
# Domain: vertical column [0, 2] m, 100 elements
# Run for 3600 s (1 hour) to see wetting front advance
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 2
  ny   = 100
  xmin = 0
  xmax = 0.5
  ymin = 0
  ymax = 2.0        # 2 m column height
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [porepressure]
    initial_condition = -5e4   # -50 kPa (initially unsaturated)
  []
[]

[PorousFlowUnsaturated]
  porepressure    = porepressure
  coupling_type   = Hydro
  gravity         = '0 -9.81 0'
  fp              = water
  relative_permeability_type = COREY
  relative_permeability_exponent = 2
  van_genuchten_alpha = 3.5e-4
  van_genuchten_m = 0.5
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
    porosity_zero = 0.35
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
    permeability = '1e-11 0 0  0 1e-11 0  0 0 1e-11'
  []
  [rock_heat]
    type                   = PorousFlowMatrixInternalEnergy
    density                = 2000
    specific_heat_capacity = 800
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.0 0 0  0 2.0 0  0 0 2.0'
  []
[]

[BCs]
  # Top: saturated (water ponded on surface)
  [top_saturated]
    type     = DirichletBC
    variable = porepressure
    boundary = top
    value    = 0      # P = 0 (atmospheric = saturated)
  []
  # Bottom: free drainage (natural zero-flux BC)
[]

[Postprocessors]
  # Pressure at base
  [p_base]
    type     = PointValue
    variable = porepressure
    point    = '0.25 0 0'
  []
  # Pressure at 0.5 m height
  [p_0p5m]
    type     = PointValue
    variable = porepressure
    point    = '0.25 0.5 0'
  []
  # Pressure at 1.0 m height
  [p_1p0m]
    type     = PointValue
    variable = porepressure
    point    = '0.25 1.0 0'
  []
  # Pressure at 1.5 m height
  [p_1p5m]
    type     = PointValue
    variable = porepressure
    point    = '0.25 1.5 0'
  []
  # Average pressure (tracks wetting progress)
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

  dt       = 60        # 60 s time steps
  end_time = 3600      # 1 hour total

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20
[]

[Outputs]
  exodus = true
  csv    = true
[]
