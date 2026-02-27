# ============================================================
# Case 49: J2 Plasticity — Uniaxial Tension with Isotropic Hardening
# A 2D plane-strain bar is pulled in the y-direction.
# Material yields at sigma_y = 250 MPa with linear hardening
# (H = 1 GPa).  Validates stress-strain curve against the
# analytical bilinear response:
#   elastic:  sigma = E * epsilon
#   plastic:  sigma = sigma_y + H * epsilon_p
#
# Domain: [0,1] x [0,5] mm, 5x25 elements (plane strain)
# Loading: displacement-controlled, disp_y up to 0.05 mm on top
# ============================================================

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 5
    ny = 25
    xmin = 0
    xmax = 1     # 1 mm wide bar
    ymin = 0
    ymax = 5     # 5 mm tall bar
  []
[]

# Ramp the top displacement linearly in time
[Functions]
  [top_pull]
    type = ParsedFunction
    expression = '0.05 * t'   # 0.05 mm max displacement at t=1
  []
  # Piecewise linear hardening function: yield stress vs plastic strain
  [hardening_func]
    type = PiecewiseLinear
    x = '0.0    0.10'    # plastic strain
    y = '250.0  350.0'   # yield stress (MPa) — H = (350-250)/0.10 = 1000 MPa
  []
[]

# QuasiStatic action sets up strain, stress divergence kernels, and variables
[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL                 # small-strain formulation
    incremental = true             # incremental strain needed for plasticity return-map
    add_variables = true           # auto-create disp_x, disp_y
    generate_output = 'stress_yy stress_xx vonmises_stress plastic_strain_yy elastic_strain_yy'
  []
[]

[BCs]
  # Bottom face: fixed in both directions (symmetry + support)
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []
  [fix_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0.0
  []

  # Top face: prescribed displacement (tension)
  [pull_top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = top_pull
  []

  # Constrain lateral expansion on the left edge (symmetry)
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []
[]

[Materials]
  # Isotropic elasticity: E = 200 GPa, nu = 0.3
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 200e3     # MPa
    poissons_ratio = 0.3
  []
  # J2 isotropic plasticity with piecewise-linear hardening
  [plasticity]
    type = IsotropicPlasticityStressUpdate
    yield_stress = 250.0       # MPa — initial yield
    hardening_function = hardening_func
  []
  # Stress computation: radial return with inelastic model
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'plasticity'
    tangent_operator = elastic
  []
[]

[Postprocessors]
  # Average axial stress (should match analytical bilinear curve)
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
  # Average plastic strain in y-direction
  [avg_plastic_strain_yy]
    type = ElementAverageValue
    variable = plastic_strain_yy
  []
  # Average elastic strain in y-direction
  [avg_elastic_strain_yy]
    type = ElementAverageValue
    variable = elastic_strain_yy
  []
  # Top displacement (loading parameter)
  [top_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  nl_max_its = 30
  l_max_its  = 50

  # 50 time steps from t=0 to t=1
  dt       = 0.02
  end_time = 1.0
[]

[Outputs]
  exodus = true
  csv    = true
[]
