# ============================================================
# Case 51: Power-Law Creep — Column Under Sustained Compression
# A 2D column under sustained compressive pressure exhibits
# time-dependent creep deformation.  Norton power-law:
#   d(eps_cr)/dt = A * sigma^n * exp(-Q/(R*T))
# At constant temperature T=1000K, the Arrhenius factor is
# exp(-3e5/(8.314*1000)) ~ 2e-16, giving observable creep
# strain accumulation over ~100 seconds.
#
# All units in SI (Pa, m, s, K) to match the MOOSE test suite.
# Parameters: A = 1e-15 /s/Pa^n, n = 4, Q = 3e5 J/mol
# E = 2e11 Pa, nu = 0.3, T = 1000 K constant
# sigma_applied = 10 MPa = 1e7 Pa via pressure BC
#
# Domain: [0,1] x [0,2], 5x10 elements (plane strain column)
# Loading: uniform pressure on top face
# ============================================================

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 10
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 2
  []
[]

# Temperature held constant (needed for Arrhenius factor in creep law)
[Variables]
  [temp]
    initial_condition = 1000.0   # K
  []
[]

# Dummy thermal kernel to keep the temperature field alive
[Kernels]
  [heat_diff]
    type = Diffusion
    variable = temp
  []
[]

# QuasiStatic action with finite strain (required for creep return-mapping)
[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    generate_output = 'stress_yy stress_xx vonmises_stress creep_strain_yy elastic_strain_yy'
  []
[]

[Functions]
  # Ramp the load over the first time step, then hold constant
  [load_ramp]
    type = PiecewiseLinear
    x = '0     1     100'
    y = '0     1     1'
  []
[]

[BCs]
  # Bottom face: fully fixed
  [fix_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0.0
  []
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []

  # Left face: symmetry (no horizontal displacement)
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []

  # Uniform downward pressure on top surface (10 MPa = 1e7 Pa)
  [Pressure]
    [top_pressure]
      boundary = top
      factor   = -1.0e7     # 10 MPa compressive (Pa)
      function = load_ramp
    []
  []

  # Temperature: hold constant everywhere
  [temp_fix]
    type = DirichletBC
    variable = temp
    boundary = 'left right top bottom'
    value = 1000.0
  []
[]

[Materials]
  # Isotropic elasticity: steel (SI units, Pa)
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.0e11    # 200 GPa in Pa
    poissons_ratio = 0.3
  []
  # Power-law creep model (Norton)
  [power_law_creep]
    type = PowerLawCreepStressUpdate
    coefficient = 1.0e-15      # A in d(eps)/dt = A * sigma^n * exp(-Q/RT)
    n_exponent = 4              # stress exponent
    activation_energy = 3.0e5   # J/mol
    temperature = temp
  []
  # Stress computation with creep
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'power_law_creep'
    tangent_operator = elastic
  []
[]

[Postprocessors]
  # Average axial stress (Pa)
  [avg_stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  []
  # Max von Mises stress
  [max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
  []
  # Average creep strain in y (should grow over time)
  [avg_creep_strain_yy]
    type = ElementAverageValue
    variable = creep_strain_yy
  []
  # Max vertical displacement (column shortening)
  [max_disp_y]
    type = ElementExtremeValue
    variable = disp_y
  []
  # Average elastic strain
  [avg_elastic_strain_yy]
    type = ElementAverageValue
    variable = elastic_strain_yy
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  line_search = 'none'

  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 30
  l_max_its  = 50
  l_tol      = 1e-5

  # 10 time steps over 100 seconds of creep
  dt       = 10.0
  end_time = 100.0
[]

[Outputs]
  exodus = true
  csv    = true
[]
