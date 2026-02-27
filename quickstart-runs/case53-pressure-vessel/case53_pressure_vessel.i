# ============================================================
# Case 53: Pressure Vessel — Thick-Walled Cylinder (Lamé Solution)
# Axisymmetric (RZ) analysis of a thick-walled cylinder under
# internal pressure.  Validates against the Lamé analytical
# solution:
#   sigma_r(r) = A - B/r^2
#   sigma_theta(r) = A + B/r^2
# where A = p_i * r_i^2 / (r_o^2 - r_i^2)
#       B = p_i * r_i^2 * r_o^2 / (r_o^2 - r_i^2)
#
# Parameters: r_i = 1 m, r_o = 2 m, p_i = 100 MPa
#   A = 100 * 1 / (4-1) = 33.33 MPa
#   B = 100 * 1 * 4 / (4-1) = 133.33 MPa*m^2
#   sigma_r(r_i)     = 33.33 - 133.33 = -100 MPa (= -p_i, check!)
#   sigma_r(r_o)     = 33.33 - 33.33  =    0 MPa (traction-free, check!)
#   sigma_theta(r_i) = 33.33 + 133.33 = 166.67 MPa (max hoop stress)
#   sigma_theta(r_o) = 33.33 + 33.33  =  66.67 MPa
#
# Domain: r in [1,2], z in [0,1] (2D axisymmetric, coord_type = RZ)
#   x = r, y = z in MOOSE convention
# Mesh: 20x4 elements
# ============================================================

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20       # radial direction
    ny = 4        # axial direction
    xmin = 1.0    # inner radius r_i = 1 m
    xmax = 2.0    # outer radius r_o = 2 m
    ymin = 0.0
    ymax = 1.0    # cylinder height (arbitrary for plane-strain-like condition)
  []
  # RZ coordinate system: x -> r, y -> z, theta is the hoop direction
  coord_type = RZ
[]

# SolidMechanics QuasiStatic action for axisymmetric
[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz vonmises_stress'
    # stress_xx = sigma_r (radial)
    # stress_yy = sigma_z (axial)
    # stress_zz = sigma_theta (hoop, out-of-plane in RZ)
  []
[]

[BCs]
  # Internal pressure on the inner surface (left boundary, r = r_i)
  [Pressure]
    [internal_pressure]
      boundary = left
      factor = 100.0     # 100 MPa internal pressure
    []
  []

  # Outer surface (right): traction-free (natural BC, no explicit BC needed)

  # Bottom: fix axial displacement (plane-strain-like constraint)
  [fix_z_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  []

  # Top: fix axial displacement (generalized plane strain)
  [fix_z_top]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0.0
  []
[]

[Materials]
  # Steel: E = 200 GPa, nu = 0.3
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 200e3   # MPa
    poissons_ratio = 0.3
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

# Analytical solution for validation
[Functions]
  [sigma_r_exact]
    type = ParsedFunction
    # A - B/r^2 where r = x in RZ, A = 33.333, B = 133.333
    expression = '33.3333 - 133.3333 / (x * x)'
  []
  [sigma_theta_exact]
    type = ParsedFunction
    # A + B/r^2
    expression = '33.3333 + 133.3333 / (x * x)'
  []
[]

[Postprocessors]
  # Radial stress at inner surface (should be -100 MPa)
  [sigma_r_inner]
    type = PointValue
    variable = stress_xx
    point = '1.0 0.5 0'
  []
  # Radial stress at outer surface (should be 0 MPa)
  [sigma_r_outer]
    type = PointValue
    variable = stress_xx
    point = '2.0 0.5 0'
  []
  # Hoop stress at inner surface (should be 166.67 MPa)
  [sigma_theta_inner]
    type = PointValue
    variable = stress_zz
    point = '1.0 0.5 0'
  []
  # Hoop stress at outer surface (should be 66.67 MPa)
  [sigma_theta_outer]
    type = PointValue
    variable = stress_zz
    point = '2.0 0.5 0'
  []
  # Max von Mises stress
  [max_vonmises]
    type = ElementExtremeValue
    variable = vonmises_stress
    value_type = max
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
