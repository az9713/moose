# ============================================================
# Case 55: 2-Group Neutron Diffusion — Fast/Thermal with Fission
# Solves the two-group neutron diffusion equations:
#   Group 1 (fast):   -D1*laplacian(phi1) + Sigma_r1*phi1 = (1/k)*chi1*[nu_Sf1*phi1 + nu_Sf2*phi2]
#   Group 2 (thermal): -D2*laplacian(phi2) + Sigma_a2*phi2 = Sigma_s12*phi1
#
# where Sigma_r1 = Sigma_a1 + Sigma_s12 (removal from group 1)
#
# Parameters (typical water-moderated reactor):
#   D1 = 1.5 cm, D2 = 0.4 cm
#   Sigma_a1 = 0.01 /cm, Sigma_a2 = 0.08 /cm
#   Sigma_s12 = 0.02 /cm (downscatter fast→thermal)
#   Sigma_r1 = Sigma_a1 + Sigma_s12 = 0.03 /cm
#   nu_Sigma_f1 = 0.005 /cm, nu_Sigma_f2 = 0.10 /cm
#   chi1 = 1.0, chi2 = 0.0 (all fission neutrons born fast)
#
# Domain: [0, 40] cm slab, 80x2 elements, vacuum BCs
# Expected k_eff ~ 1.0 for this slab width
# ============================================================

L = 40.0

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 80
    ny = 2
    xmin = 0
    xmax = ${L}
    ymin = 0
    ymax = 0.5
  []
[]

[Variables]
  [phi1]   # fast group flux
  []
  [phi2]   # thermal group flux
  []
[]

[Kernels]
  # --- Group 1 (fast) ---
  # Diffusion: -D1*laplacian(phi1)
  [diff1]
    type = MatDiffusion
    variable = phi1
    diffusivity = D1
  []
  # Removal: +Sigma_r1*phi1  (absorption + downscatter out)
  [removal1]
    type = Reaction
    variable = phi1
    rate = 0.03   # Sigma_r1 = 0.03 /cm
  []
  # Fission source from group 1: -(1/k)*chi1*nu_Sf1*phi1
  [fission11]
    type = Reaction
    variable = phi1
    rate = -0.005   # -nu_Sigma_f1
    extra_vector_tags = 'eigen'
  []
  # Fission source from group 2 into group 1: -(1/k)*chi1*nu_Sf2*phi2
  [fission21]
    type = CoupledForce
    variable = phi1
    v = phi2
    coef = 0.10   # nu_Sigma_f2 (positive because CoupledForce has -coef*v*test)
    extra_vector_tags = 'eigen'
  []

  # --- Group 2 (thermal) ---
  # Diffusion: -D2*laplacian(phi2)
  [diff2]
    type = MatDiffusion
    variable = phi2
    diffusivity = D2
  []
  # Absorption: +Sigma_a2*phi2
  [absorption2]
    type = Reaction
    variable = phi2
    rate = 0.08   # Sigma_a2 = 0.08 /cm
  []
  # Downscatter source from group 1: -Sigma_s12*phi1
  [scatter12]
    type = CoupledForce
    variable = phi2
    v = phi1
    coef = 0.02   # Sigma_s12 (positive, because CoupledForce residual = -coef*v*test)
  []
[]

[Materials]
  [group_diffusion]
    type = GenericConstantMaterial
    prop_names  = 'D1   D2'
    prop_values = '1.5  0.4'
  []
[]

[BCs]
  # Vacuum BCs for both groups
  [left_phi1]
    type = DirichletBC
    variable = phi1
    boundary = left
    value = 0
  []
  [right_phi1]
    type = DirichletBC
    variable = phi1
    boundary = right
    value = 0
  []
  [left_phi2]
    type = DirichletBC
    variable = phi2
    boundary = left
    value = 0
  []
  [right_phi2]
    type = DirichletBC
    variable = phi2
    boundary = right
    value = 0
  []
  # Eigen BCs
  [left_phi1_eigen]
    type = EigenDirichletBC
    variable = phi1
    boundary = left
  []
  [right_phi1_eigen]
    type = EigenDirichletBC
    variable = phi1
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
  [flux_profiles]
    type = LineValueSampler
    variable = 'phi1 phi2'
    start_point = '0 0.25 0'
    end_point = '${L} 0.25 0'
    num_points = 81
    sort_by = x
  []
[]

[Postprocessors]
  # Fast flux at center
  [phi1_center]
    type = PointValue
    variable = phi1
    point = '20 0.25 0'
  []
  # Thermal flux at center
  [phi2_center]
    type = PointValue
    variable = phi2
    point = '20 0.25 0'
  []
[]

[Outputs]
  exodus = true
  csv    = true
[]
