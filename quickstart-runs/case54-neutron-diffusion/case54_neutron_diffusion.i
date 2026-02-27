# ============================================================
# Case 54: 1-Group Neutron Diffusion — Bare Slab Criticality
# Solves the one-speed neutron diffusion eigenvalue problem:
#   -D * laplacian(phi) + Sigma_a * phi = (1/k) * nu_Sigma_f * phi
# on a 1D slab [0, L] with vacuum BCs (phi = 0 at boundaries).
#
# Analytical critical half-thickness:
#   B^2 = (nu_Sigma_f - Sigma_a) / D = (pi/L)^2
# For D = 1.0 cm, Sigma_a = 0.1 /cm, nu_Sigma_f = 0.15 /cm:
#   B^2 = (0.15 - 0.10) / 1.0 = 0.05
#   L_crit = pi / sqrt(0.05) = 14.05 cm
#
# We use L = 14.05 cm so k_eff ~ 1.0 (critical).
# The eigenvalue solver finds k_eff and the fundamental mode
# phi(x) ~ cos(pi*x/L) (shifted to match [0,L] domain).
#
# Domain: [0, 14.05] cm, 100 elements (quasi-1D: 100x2)
# ============================================================

L = 14.05   # slab thickness (cm) — near-critical

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
  [phi]   # scalar neutron flux
  []
[]

[Kernels]
  # Diffusion term: -D * laplacian(phi) → integral(D * grad(phi) . grad(test))
  [diffusion]
    type = MatDiffusion
    variable = phi
    diffusivity = D
  []
  # Absorption term: +Sigma_a * phi
  [absorption]
    type = Reaction
    variable = phi
    rate = 0.1    # Sigma_a = 0.1 /cm
  []
  # Fission source: -(1/k) * nu_Sigma_f * phi  (goes to eigen tag)
  [fission]
    type = Reaction
    variable = phi
    rate = -0.15   # -nu_Sigma_f = -0.15 /cm (negative => source)
    extra_vector_tags = 'eigen'
  []
[]

[Materials]
  [diffusion_coeff]
    type = GenericConstantMaterial
    prop_names = 'D'
    prop_values = '1.0'   # diffusion coefficient (cm)
  []
[]

[BCs]
  # Vacuum BCs: phi = 0 at both ends (extrapolated boundary)
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
  # Eigenvalue system also needs EigenDirichletBC
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

[Executioner]
  type = Eigenvalue
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[VectorPostprocessors]
  # Sample flux profile along the slab centerline
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
  # k-effective (the eigenvalue)
  [k_eff]
    type = VectorPostprocessorComponent
    vectorpostprocessor = flux_profile
    vector_name = phi
    index = 50   # midpoint flux (for normalization check)
  []
[]

[Outputs]
  exodus = true
  csv    = true
[]
