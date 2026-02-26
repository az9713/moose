# ============================================================
# Case 14: Thermoelasticity — Heated Plate with Thermal Stress
# Steady-state heat conduction drives thermal expansion in a
# 2D elastic solid.  Hot left (T=500K), cold right (T=300K).
# One-way coupling: T field -> eigenstrain -> displacement/stress
# Requires: combined-opt  (heat_transfer + solid_mechanics modules)
# ============================================================

# -------------------------------------------------------------
# Material constants for structural steel
# -------------------------------------------------------------
E     = 200e9   # Young's modulus, Pa
nu    = 0.3     # Poisson's ratio, dimensionless
alpha = 12e-6   # coefficient of thermal expansion, 1/K
k_th  = 50      # thermal conductivity, W/(m K)
cp    = 500     # specific heat, J/(kg K)
T_ref = 300     # stress-free (reference) temperature, K

[GlobalParams]
  # Displacement variable names used by both the SolidMechanics
  # action and any kernel/BC/postprocessor that needs them.
  displacements = 'disp_x disp_y'
[]

[Mesh]
  # 20x20 structured quad mesh on a 1 m x 1 m plate.
  # Named boundaries on a 2D GeneratedMesh:
  #   left   (x=0), right  (x=1)
  #   bottom (y=0), top    (y=1)
  type = GeneratedMesh
  dim  = 2
  nx   = 20
  ny   = 20
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

# -------------------------------------------------------------
# Temperature variable (heat conduction problem)
# Displacement variables are created automatically by the
# Physics/SolidMechanics/QuasiStatic action (add_variables=true)
# -------------------------------------------------------------
[Variables]
  [T]
    # Initial guess: linear ramp from 300 to 500 K helps convergence
    [InitialCondition]
      type  = FunctionIC
      function = '300 + 200*(1-x)'
    []
  []
[]

# -------------------------------------------------------------
# Kernels
# -------------------------------------------------------------
[Kernels]
  # Steady-state heat conduction: -div( k * grad(T) ) = 0
  # ADHeatConduction reads thermal_conductivity from the material.
  [heat_conduction]
    type     = ADHeatConduction
    variable = T
  []
[]

# -------------------------------------------------------------
# SolidMechanics QuasiStatic action
# This is the modern replacement for [Modules/TensorMechanics/Master].
# It automatically:
#   - Creates disp_x and disp_y variables (add_variables = true)
#   - Adds ADStressDivergenceTensors kernels for each component
#   - Adds ADComputeSmallStrain material
#   - Creates an AuxVariable and AuxKernel for vonmises_stress
# -------------------------------------------------------------
[Physics/SolidMechanics/QuasiStatic]
  [solid]
    strain             = SMALL                # small-strain (linear) formulation
    add_variables      = true                 # generate disp_x, disp_y automatically
    eigenstrain_names  = 'thermal_eigenstrain' # must match eigenstrain_name in Material
    generate_output    = 'vonmises_stress'    # creates vonmises_stress AuxVariable
    use_automatic_differentiation = true      # use AD for exact Jacobian
  []
[]

# -------------------------------------------------------------
# Boundary conditions
# -------------------------------------------------------------
[BCs]
  # Temperature: hot left wall, cold right wall
  [T_hot]
    type     = DirichletBC
    variable = T
    boundary = left
    value    = 500   # K
  []
  [T_cold]
    type     = DirichletBC
    variable = T
    boundary = right
    value    = 300   # K
  []

  # Displacement: pin the bottom edge (y=0) to prevent rigid-body motion.
  # disp_x = 0 prevents horizontal sliding.
  # disp_y = 0 prevents vertical lifting.
  [pin_bottom_x]
    type     = DirichletBC
    variable = disp_x
    boundary = bottom
    value    = 0
  []
  [pin_bottom_y]
    type     = DirichletBC
    variable = disp_y
    boundary = bottom
    value    = 0
  []
[]

# -------------------------------------------------------------
# Materials
# -------------------------------------------------------------
[Materials]
  # Thermal material: provides 'thermal_conductivity' and 'specific_heat'
  # properties consumed by ADHeatConduction kernel.
  [thermal_props]
    type                = ADHeatConductionMaterial
    thermal_conductivity = ${k_th}   # W/(m K)
    specific_heat        = ${cp}      # J/(kg K)
  []

  # Isotropic linear elasticity tensor from E and nu.
  # Provides the 'elasticity_tensor' material property.
  [elasticity_tensor]
    type          = ADComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  []

  # Thermal eigenstrain: maps temperature change to a stress-free
  # volumetric strain.  The eigenstrain name must match the list
  # declared in the QuasiStatic action above.
  [thermal_eigenstrain]
    type                    = ADComputeThermalExpansionEigenstrain
    temperature             = T               # coupled temperature variable
    thermal_expansion_coeff = ${alpha}        # 1/K
    stress_free_temperature  = ${T_ref}        # K — strain is zero at this temperature
    eigenstrain_name        = thermal_eigenstrain
  []

  # Linear elastic stress from total mechanical strain
  # (total strain minus eigenstrain).
  [stress]
    type = ADComputeLinearElasticStress
  []
[]

# -------------------------------------------------------------
# Postprocessors
# -------------------------------------------------------------
[Postprocessors]
  # --- Temperature diagnostics ---
  [max_temperature]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
  [avg_temperature]
    type     = ElementAverageValue
    variable = T
  []

  # --- Displacement diagnostics ---
  # Horizontal (x) displacement driven by thermal gradient
  [max_disp_x]
    type       = ElementExtremeValue
    variable   = disp_x
    value_type = max
  []
  # Vertical (y) displacement due to Poisson effect and constraint
  [max_disp_y]
    type       = ElementExtremeValue
    variable   = disp_y
    value_type = max
  []

  # --- Stress diagnostic ---
  [max_vonmises]
    type       = ElementExtremeValue
    variable   = vonmises_stress
    value_type = max
  []
[]

# -------------------------------------------------------------
# Executioner
# -------------------------------------------------------------
[Executioner]
  # Steady solves the coupled T + (disp_x, disp_y) system once.
  type = Steady

  # NEWTON uses the exact (AD-assembled) Jacobian — efficient for
  # this fully-coupled linear problem.
  solve_type = 'NEWTON'

  # LU direct factorization: ideal for small-to-medium 2D problems
  # where the system fits comfortably in memory.
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

# -------------------------------------------------------------
# Outputs
# -------------------------------------------------------------
[Outputs]
  # Exodus file contains: T, disp_x, disp_y, vonmises_stress fields
  exodus = true

  # CSV file contains: postprocessor scalar values
  csv = true
[]
