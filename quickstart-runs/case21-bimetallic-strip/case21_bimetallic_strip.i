# ============================================================
# Case 21: Thermo-Mechanical Bimetallic Strip
# Two bonded metals heated uniformly — differential thermal
# expansion causes the strip to bend (classic bimetallic effect).
# Steel (bottom, block 0): alpha=12e-6 /K, E=200 GPa, nu=0.3
# Aluminum (top, block 1): alpha=23e-6 /K, E=70 GPa, nu=0.33
# T_ref = 300 K  ->  T_final = 500 K  (deltaT = 200 K)
# Pinned at left end, free at right end.
# Requires: combined-opt (SolidMechanics module)
# ============================================================

[Mesh]
  # Stage 1: create a 40x8 strip, 10 m long, 1 m tall.
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 8
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 1
  []

  # Stage 2: reassign the top half (y in [0.5, 1.0]) to block 1 (aluminum).
  # Elements NOT reassigned remain in block 0 (steel).
  # ny=8 gives element edges at y = 0, 0.125, 0.25, 0.375, 0.5, ...
  # The interface at y=0.5 coincides exactly with element boundaries.
  [top_block]
    type        = SubdomainBoundingBoxGenerator
    input       = gen
    bottom_left = '0 0.5 0'
    top_right   = '10 1.0 0'
    block_id    = 1
  []
[]

# Declare that disp_x and disp_y are the displacement variables used
# throughout the solid mechanics system.
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

# The QuasiStatic action creates the displacement variables, the
# TotalLagrangian (or small-strain) kernels, and the stress output.
# eigenstrain_names links the thermal expansion eigenstrain into the
# constitutive update for every element on every block.
[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain              = SMALL
        add_variables       = true
        generate_output     = 'vonmises_stress stress_xx stress_yy'
        eigenstrain_names   = 'thermal_eigenstrain'
      []
    []
  []
[]

# T is prescribed (not solved) — it represents the uniform final
# temperature after heating.  initial_condition sets T = 500 K
# throughout the domain.  The stress-free reference is 300 K,
# so the effective deltaT = 200 K.
[AuxVariables]
  [T]
    initial_condition = 500
  []
[]

[Materials]
  # ---- Bottom strip: Steel (block 0) ----
  [elasticity_steel]
    type           = ComputeIsotropicElasticityTensor
    block          = 0
    youngs_modulus = 200e9
    poissons_ratio = 0.3
  []
  [thermal_expansion_steel]
    type                    = ComputeThermalExpansionEigenstrain
    block                   = 0
    eigenstrain_name        = thermal_eigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 12e-6
    temperature             = T
  []
  [stress_steel]
    type  = ComputeLinearElasticStress
    block = 0
  []

  # ---- Top strip: Aluminum (block 1) ----
  [elasticity_aluminum]
    type           = ComputeIsotropicElasticityTensor
    block          = 1
    youngs_modulus = 70e9
    poissons_ratio = 0.33
  []
  [thermal_expansion_aluminum]
    type                    = ComputeThermalExpansionEigenstrain
    block                   = 1
    eigenstrain_name        = thermal_eigenstrain
    stress_free_temperature = 300
    thermal_expansion_coeff = 23e-6
    temperature             = T
  []
  [stress_aluminum]
    type  = ComputeLinearElasticStress
    block = 1
  []
[]

[BCs]
  # Pin the entire left edge: no horizontal or vertical displacement.
  # This removes all rigid-body modes while allowing the strip to
  # bend freely to the right.
  [fix_x_left]
    type     = DirichletBC
    variable = disp_x
    boundary = left
    value    = 0
  []
  [fix_y_left]
    type     = DirichletBC
    variable = disp_y
    boundary = left
    value    = 0
  []
[]

[Postprocessors]
  # Maximum upward displacement anywhere in the strip.
  [max_disp_y]
    type       = NodalExtremeValue
    variable   = disp_y
    value_type = max
  []
  # Maximum downward displacement (bending direction).
  [min_disp_y]
    type       = NodalExtremeValue
    variable   = disp_y
    value_type = min
  []
  # Vertical displacement at the free tip (midline of right edge).
  # Aluminum expands more -> strip bends toward steel -> tip moves down.
  [tip_disp_y]
    type     = PointValue
    variable = disp_y
    point    = '10 0.5 0'
  []
  # Peak von Mises stress (expected at or near the bonded interface).
  [max_vonmises]
    type     = ElementExtremeValue
    variable = vonmises_stress
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv    = true
[]
