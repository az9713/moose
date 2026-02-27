# ============================================================
# Case 56: Fuel Pin Heat Transfer — Radial Temperature Profile
# Axisymmetric (RZ) steady-state heat conduction in a nuclear
# fuel pin with:
#   - UO2 fuel pellet (r = 0 to r_fuel = 4.1 mm)
#   - Gap (modeled as thin region, r_fuel to r_clad_in = 4.2 mm)
#   - Zircaloy cladding (r_clad_in to r_clad_out = 4.75 mm)
#
# Volumetric heat generation q''' = 4e8 W/m^3 in fuel only.
# Coolant at T_cool = 573 K with heat transfer coeff h = 3e4 W/m^2/K.
#
# Analytical centerline temperature (fuel only, no gap):
#   T_center = T_surface + q'''*r_f^2/(4*k_fuel)
#   ~ 573 + 4e8*(4.1e-3)^2 / (4*3.0) ~ 573 + 561 = 1134 K
#
# With gap resistance, T_center will be higher (~1200-1400 K).
#
# Domain: r in [0, 4.75e-3] m, z in [0, 0.01] m
# Two blocks: fuel (block 0, r<4.1e-3) and clad (block 1, r>4.2e-3)
# ============================================================

[Mesh]
  # Simple single-block mesh for fuel + cladding
  # (Gap conductance modeled via material property jump)
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40       # radial direction
    ny = 2        # axial (short segment)
    xmin = 0.0
    xmax = 4.75e-3   # outer cladding radius (m)
    ymin = 0.0
    ymax = 0.01      # 1 cm axial segment
  []
  coord_type = RZ   # x = r, y = z
[]

[Variables]
  [T]
    initial_condition = 600.0   # K
  []
[]

[Kernels]
  # Heat conduction: -div(k*grad(T)) = q'''
  [conduction]
    type = ADHeatConduction
    variable = T
  []
  # Volumetric heat source in fuel region only (r < 4.1e-3 m)
  [heat_source]
    type = BodyForce
    variable = T
    function = heat_gen
  []
[]

[Functions]
  # Volumetric heat generation: uniform in fuel (r < 4.1 mm), zero in gap/clad
  [heat_gen]
    type = ParsedFunction
    expression = 'if(x < 4.1e-3, 4.0e8, 0.0)'   # W/m^3
  []
  # Spatially varying thermal conductivity
  [k_func]
    type = ParsedFunction
    # Fuel (r < 4.1 mm): k = 3.0 W/m/K (UO2)
    # Gap (4.1-4.2 mm): k = 0.5 W/m/K (He gas)
    # Clad (> 4.2 mm): k = 16.0 W/m/K (Zircaloy)
    expression = 'if(x < 4.1e-3, 3.0, if(x < 4.2e-3, 0.5, 16.0))'
  []
[]

[Materials]
  [thermal]
    type = ADGenericFunctionMaterial
    prop_names  = 'thermal_conductivity'
    prop_values = 'k_func'
  []
  # Density * specific heat (not needed for steady-state but required by material)
  [density_cp]
    type = ADGenericConstantMaterial
    prop_names = 'density specific_heat'
    prop_values = '10000.0 300.0'
  []
[]

[BCs]
  # Outer surface: convective cooling T_cool = 573 K, h = 30000 W/m^2/K
  [coolant]
    type = ConvectiveHeatFluxBC
    variable = T
    boundary = right
    T_infinity = 573.0   # K — coolant bulk temperature
    heat_transfer_coefficient = 30000.0   # W/m^2/K
  []
  # Top/bottom: insulated (natural BC)
  # Left (r=0): symmetry axis — natural BC (zero flux)
[]

[Postprocessors]
  # Fuel centerline temperature (should be ~1200-1400 K)
  [T_centerline]
    type = PointValue
    variable = T
    point = '0 0.005 0'
  []
  # Fuel surface temperature (at r = 4.1 mm)
  [T_fuel_surface]
    type = PointValue
    variable = T
    point = '4.1e-3 0.005 0'
  []
  # Clad inner temperature
  [T_clad_inner]
    type = PointValue
    variable = T
    point = '4.2e-3 0.005 0'
  []
  # Clad outer temperature
  [T_clad_outer]
    type = PointValue
    variable = T
    point = '4.75e-3 0.005 0'
  []
  # Average fuel temperature
  [T_fuel_avg]
    type = ElementAverageValue
    variable = T
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv    = true
[]
