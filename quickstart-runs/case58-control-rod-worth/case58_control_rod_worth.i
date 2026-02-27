# ============================================================
# Case 58: Control Rod Worth — Eigenvalue Shift from Absorber
# Compares k_eff of a bare slab vs. a slab with a central
# absorber region, demonstrating reactivity control.
#
# Physics: 1-group neutron diffusion eigenvalue problem
#   -D * laplacian(phi) + Sigma_a(x) * phi = (1/k) * nu_Sf * phi
#
# Two configurations computed in a single run using a
# spatially varying absorption cross section:
#   Sigma_a(x) = Sigma_a0 + Delta_Sigma_a * rod_function(x)
#
# The control rod is modeled as a region of enhanced absorption
# in the center of the slab (x in [8, 12] out of [0, 20]).
#
# Parameters:
#   D = 1.0 cm
#   Sigma_a0 = 0.08 /cm (base absorption)
#   nu_Sigma_f = 0.12 /cm (fission source)
#   Delta_Sigma_a = 0.05 /cm (rod absorption penalty)
#   L = 20 cm slab
#
# Without rod: k_inf = nu_Sf/Sigma_a0 = 0.12/0.08 = 1.50
#   B^2 = (nu_Sf - Sigma_a0)/D = 0.04
#   k_eff = nu_Sf / (Sigma_a0 + D*B^2) = nu_Sf / (Sigma_a0 + 0.04)
#   With vacuum BCs: B = pi/L = 0.1571, B^2 = 0.02467
#   k_eff ~ 0.12 / (0.08 + 0.02467) = 1.146
#
# With rod: enhanced absorption in center depresses flux there,
#   reducing k_eff. The reactivity worth rho = (k_rod - k_no_rod)/k_rod.
#
# Domain: [0, 20] cm, 100 elements, quasi-1D (100x2)
# ============================================================

L = 20.0

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 2
    xmin = 0
    xmax = ${L}
    ymin = 0
    ymax = 0.5
  []
[]

[Variables]
  [phi]
    initial_condition = 1.0
  []
[]

[Functions]
  # Spatially varying absorption: base + rod in center [8, 12]
  # Negated because ADMatReaction residual = -mob*u*test,
  # but we need +Sigma_a*phi*test (positive = loss)
  [sigma_a_func]
    type = ParsedFunction
    expression = 'if(x > 8.0 & x < 12.0, -0.13, -0.08)'
    # -0.13 = -(0.08 + 0.05) in rod region, -0.08 elsewhere
  []
[]

[Kernels]
  # Diffusion: -D * laplacian(phi)
  [diffusion]
    type = MatDiffusion
    variable = phi
    diffusivity = D
  []
  # Spatially varying absorption: +Sigma_a(x) * phi
  # ADMatReaction residual = -reaction_rate * u * test
  # We need +Sigma_a*phi*test, so reaction_rate = -Sigma_a (negated in function)
  [absorption]
    type = ADMatReaction
    variable = phi
    reaction_rate = sigma_a_ad
  []
  # Fission source: -(1/k) * nu_Sf * phi (tagged as eigen)
  [fission]
    type = Reaction
    variable = phi
    rate = -0.12    # -nu_Sigma_f
    extra_vector_tags = 'eigen'
  []
[]

[Materials]
  [diffusion_coeff]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '1.0'
  []
  # Position-dependent absorption cross section via function
  [sigma_a_mat]
    type = ADGenericFunctionMaterial
    prop_names = 'sigma_a_ad'
    prop_values = 'sigma_a_func'
  []
[]

[BCs]
  # Vacuum BCs: phi = 0 at both ends
  [left_vacuum]
    type = DirichletBC
    variable = phi
    boundary = left
    value = 0
  []
  [right_vacuum]
    type = DirichletBC
    variable = phi
    boundary = right
    value = 0
  []
  # EigenDirichletBC for eigenvalue system
  [left_eigen]
    type = EigenDirichletBC
    variable = phi
    boundary = left
  []
  [right_eigen]
    type = EigenDirichletBC
    variable = phi
    boundary = right
  []
[]

[VectorPostprocessors]
  # Sample flux profile along the slab
  [flux_profile]
    type = LineValueSampler
    variable = phi
    start_point = '0 0.25 0'
    end_point = '${L} 0.25 0'
    num_points = 101
    sort_by = x
  []
[]

[Postprocessors]
  # Flux at center of slab (in rod region)
  [phi_center]
    type = PointValue
    variable = phi
    point = '10 0.25 0'
  []
  # Flux at quarter point (outside rod region)
  [phi_quarter]
    type = PointValue
    variable = phi
    point = '5 0.25 0'
  []
[]

[Executioner]
  type = Eigenvalue
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
  csv    = true
[]
