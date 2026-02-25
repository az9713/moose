# ============================================================
# Case 6: Two-Region Domain with Different Conductivities
# Left half:  k=1.0, right half: k=5.0
# u=0 on left boundary, u=1 on right boundary
# ============================================================

[Mesh]
  # We build the mesh in stages using MeshGenerators.
  # First create a uniform mesh, then carve out a second subdomain.
  [gen]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 20
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []

  # SubdomainBoundingBoxGenerator reassigns all elements whose centroid
  # falls inside the given box to a new block id.
  [right_half]
    type        = SubdomainBoundingBoxGenerator
    input       = gen           # takes the mesh from [gen] above
    bottom_left = '0.5 0.0 0'  # box corners (z ignored in 2-D)
    top_right   = '1.0 1.0 0'
    block_id    = 2             # new subdomain id for the right half
  []
  # Elements NOT reassigned remain in block 0 (the default).
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [diffusion]
    type        = MatDiffusion
    variable    = u
    diffusivity = k   # material property 'k', different per block
  []
[]

[BCs]
  [left]
    type     = DirichletBC
    variable = u
    boundary = left
    value    = 0
  []
  [right]
    type     = DirichletBC
    variable = u
    boundary = right
    value    = 1
  []
[]

[Materials]
  # Left half (block 0, default) — lower conductivity.
  [mat_left]
    type        = GenericConstantMaterial
    block       = 0          # applies only to elements in block 0 (default)
    prop_names  = 'k'
    prop_values = '1.0'
  []

  # Right half (block 2) — five times higher conductivity.
  [mat_right]
    type        = GenericConstantMaterial
    block       = 2          # applies only to elements in block 2
    prop_names  = 'k'
    prop_values = '5.0'
  []
[]

[Postprocessors]
  # Compute average u in each subdomain separately.
  [avg_u_left]
    type     = ElementAverageValue
    variable = u
    block    = 0
  []
  [avg_u_right]
    type     = ElementAverageValue
    variable = u
    block    = 2
  []
[]

[Executioner]
  type       = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
  csv    = true
[]
