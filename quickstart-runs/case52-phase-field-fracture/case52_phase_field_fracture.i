# ============================================================
# Case 52: Phase-Field Fracture — Notched Specimen Under Tension
# A 2D domain with a pre-existing edge notch (from left side to
# center) is pulled vertically.  The phase-field damage variable
# c evolves from 0 (intact) to 1 (fully cracked) based on the
# competition between elastic strain energy and fracture energy.
#
# Physics:
#   - Elasticity with degradation: sigma = g(c) * C : epsilon
#   - Damage evolution: gc/(2l) * c - gc*l * laplacian(c) = -g'(c)*psi_e
#   - g(c) = (1-c)^2 + eta  (quadratic degradation)
#   - Spectral decomposition prevents damage in compression
#
# Parameters: E=210 GPa, nu=0.3, gc=2.7e-3 kN/mm, l=0.04 mm
# Domain: [0,1] x [0,0.5], notch from (0,0) to (0.5,0) on bottom
# ============================================================

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 20
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 0.5
  []
  # Create a node set for the "noncrack" portion of the bottom boundary.
  # The notch is on the left half of the bottom (x=0 to x=0.5).
  # The noncrack region is the right half (x=0.5 to x=1.0).
  [noncrack]
    type = BoundingBoxNodeSetGenerator
    new_boundary = noncrack
    bottom_left = '0.5 0 0'
    top_right = '1 0 0'
    input = gen
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

# Phase-field damage variable 'c' is solved by the Nonconserved action.
# This creates: TimeDerivative, ACBulk (Allen-Cahn), ACInterface kernels.
[Modules]
  [PhaseField]
    [Nonconserved]
      [c]
        free_energy = F
        kappa = kappa_op
        mobility = L
      []
    []
  []
[]

# Solid mechanics action sets up displacement variables and stress divergence
[Physics/SolidMechanics/QuasiStatic]
  [mech]
    add_variables = true
    strain = SMALL
    additional_generate_output = 'stress_yy vonmises_stress'
    save_in = 'resid_x resid_y'
  []
[]

# Residual auxiliary variables (for reaction force postprocessors)
[AuxVariables]
  [resid_x]
  []
  [resid_y]
  []
[]

# Off-diagonal coupling between displacement and damage
[Kernels]
  [solid_x]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_x
    component = 0
    c = c
  []
  [solid_y]
    type = PhaseFieldFractureMechanicsOffDiag
    variable = disp_y
    component = 1
    c = c
  []
[]

[BCs]
  # Pull the top face upward
  [ydisp_top]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = '1e-3 * t'   # slow loading to resolve crack path
  []
  # Fix horizontal displacement on top to prevent shearing
  [xfix_top]
    type = DirichletBC
    variable = disp_x
    boundary = top
    value = 0
  []
  # Fix y-displacement on the noncrack portion of bottom (right half)
  [yfix_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = noncrack
    value = 0
  []
[]

[Materials]
  # Phase-field parameters
  [pf_constants]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco'
    prop_values = '2.7e-3 0.04 1e-4'   # gc (kN/mm), regularization length, viscosity
  []
  # Mobility L = 1/(gc * visco)
  [define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    property_name = L
    expression = '1.0/(gc_prop * visco)'
  []
  # Gradient energy coefficient kappa = gc * l
  [define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    property_name = kappa_op
    expression = 'gc_prop * l'
  []
  # Elasticity tensor (steel)
  [elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '120.0 80.0'
    fill_method = symmetric_isotropic   # lambda=80 GPa, mu=120 GPa → E~210 GPa, nu~0.3
  []
  # Stress with phase-field degradation and spectral decomposition
  [damage_stress]
    type = ComputeLinearElasticPFFractureStress
    c = c
    E_name = 'elastic_energy'
    D_name = 'degradation'
    F_name = 'local_fracture_energy'
    decomposition_type = strain_spectral
  []
  # Quadratic degradation function: g(c) = (1-c)^2 + eta
  [degradation]
    type = DerivativeParsedMaterial
    property_name = degradation
    coupled_variables = 'c'
    expression = '(1.0-c)^2*(1.0 - eta) + eta'
    constant_names       = 'eta'
    constant_expressions = '1e-6'
    derivative_order = 2
    disable_fpoptimizer = true
    enable_jit = false
  []
  # Local fracture energy density: gc/(2l) * c^2
  [local_fracture_energy]
    type = DerivativeParsedMaterial
    property_name = local_fracture_energy
    coupled_variables = 'c'
    material_property_names = 'gc_prop l'
    expression = 'c^2 * gc_prop / 2 / l'
    derivative_order = 2
    disable_fpoptimizer = true
    enable_jit = false
  []
  # Total free energy: elastic_energy + local_fracture_energy
  [fracture_driving_energy]
    type = DerivativeSumMaterial
    coupled_variables = c
    sum_materials = 'elastic_energy local_fracture_energy'
    derivative_order = 2
    property_name = F
  []
[]

[Postprocessors]
  # Reaction force on top boundary (load-displacement curve)
  [reaction_force_y]
    type = NodalSum
    variable = resid_y
    boundary = top
  []
  # Max damage
  [max_damage]
    type = ElementExtremeValue
    variable = c
    value_type = max
  []
  # Displacement of top surface
  [top_disp_y]
    type = SideAverageValue
    variable = disp_y
    boundary = top
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30
  l_max_its  = 50

  dt = 1.0
  end_time = 10.0
[]

[Outputs]
  exodus = true
  csv    = true
[]
