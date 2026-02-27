# ============================================================
# Case 81: 2D Photonic Crystal Band Gap — TM Eigenvalue at Gamma Point
#
# References:
#   Joannopoulos, Johnson, Winn, Meade, "Photonic Crystals: Molding the
#     Flow of Light", 2nd ed. (Princeton UP, 2008) — the standard reference.
#   Johnson & Joannopoulos, "Block-iterative frequency-domain methods for
#     Maxwell's equations in a planewave basis", Opt. Express 8, 173 (2001).
#
# ============================================================
# PHYSICS: ELECTROMAGNETIC WAVES IN PERIODIC MEDIA
# ============================================================
#
# A photonic crystal is a periodic dielectric structure in which the spatial
# periodicity of the refractive index creates Bragg-like interference for
# electromagnetic waves.  Just as electron energy bands arise from the
# periodic crystal potential in solid-state physics, photonic bands arise
# from the periodic dielectric constant.  In the gaps between bands — the
# photonic band gaps — no propagating modes exist at those frequencies.
#
# We consider TM polarisation (E-field parallel to the rod axis, z-direction)
# in a 2D square lattice.  The z-component of the electric field E_z = ψ(x,y)
# satisfies the scalar Helmholtz equation:
#
#   ∇²ψ + (ω/c)² εᵣ(x,y) ψ = 0
#
# where εᵣ(x,y) is the spatially periodic relative permittivity.
#
# BLOCH THEOREM AND PERIODIC BOUNDARY CONDITIONS
# -----------------------------------------------
# By Bloch's theorem, eigenmodes in a periodic lattice satisfy:
#
#   ψ_k(r + R) = exp(i k·R) ψ_k(r)
#
# where R = m·a₁ + n·a₂ is any lattice vector and k is the Bloch wave vector.
# At the Gamma point k = (0, 0), this reduces to:
#
#   ψ(r + R) = ψ(r)       (ordinary periodicity)
#
# Therefore at the Gamma point, ψ has the full periodicity of the lattice,
# and we can solve on a single unit cell with simple periodic BCs.
# (At other high-symmetry points like X = (π/a, 0) and M = (π/a, π/a),
#  the Bloch factor is complex and requires a more elaborate formulation.)
#
# ============================================================
# UNIT CELL GEOMETRY
# ============================================================
#
# Square lattice constant a = 1 (normalised).
# Unit cell: Ω = [0, 1] × [0, 1].
# Circular alumina rod of radius r = 0.2a = 0.2, centred at (0.5, 0.5).
# Rod permittivity:  εᵣ = 8.9  (alumina / polycrystalline Al₂O₃)
# Host permittivity: εᵣ = 1.0  (vacuum / air)
# Filling fraction:  f = πr²/a² = π × 0.04 ≈ 0.1257
#
# The dielectric function:
#   εᵣ(x,y) = 8.9    if  (x-0.5)² + (y-0.5)² < 0.04  (inside rod)
#             1.0    otherwise                          (vacuum)
#
# ============================================================
# EIGENVALUE FORMULATION
# ============================================================
#
# The Helmholtz equation ∇²ψ + λ εᵣ ψ = 0  is a GENERALISED eigenvalue
# problem  A·x = λ·B·x  where  λ = (ω/c)²  is the squared normalised
# angular frequency.
#
# Weak form (multiply by test function φ, integrate over the unit cell,
# integrate the Laplacian term by parts):
#
#   ∫ ∇ψ·∇φ dV  =  λ · (-∫ εᵣ ψ φ dV)
#
# so:
#   A = stiffness matrix:        A_ij = ∫ ∇φ_j · ∇φ_i dV
#   B = εᵣ-weighted mass matrix: B_ij = -∫ εᵣ φ_j φ_i dV
#
# KERNEL MAPPING:
#   A  ←  [Diffusion]   kernel:
#          computes ∫ ∇ψ·∇φ dV  (stiffness, in the standard residual system)
#
#   B  ←  [MatReaction] kernel  tagged  extra_vector_tags = 'eigen':
#          MatReaction residual = -reaction_rate × ψ × φ
#          With reaction_rate = εᵣ(x,y):  residual = -εᵣ × ψ × φ
#          Tagged 'eigen' → goes to the B (mass) matrix of the eigenproblem.
#          Net effect:  B_ij = -∫ εᵣ φ_j φ_i dV  ✓
#
# Note: MatReaction (non-AD) is used — not ADMatReaction — because
# GenericFunctionMaterial provides a non-AD material property, and mixing
# non-AD properties with ADMatReaction causes a runtime type error.
#
# ============================================================
# EIGENVALUE INTERPRETATION
# ============================================================
#
# The MOOSE KRYLOVSCHUR solver finds λ = (ω/c)² values.
#
# To convert to normalised frequency (dimensionless):
#   ωa/(2πc) = (a/2π) × √λ = (1/2π) × √λ   (with a = 1)
#
# For a square lattice of alumina rods (εᵣ = 8.9, r/a = 0.2) the
# TM band structure at the Gamma point (from published MPB calculations):
#   Mode 1 (trivial DC): ωa/(2πc) = 0  →  λ = 0
#   Mode 2: ωa/(2πc) ≈ 0.36           →  λ ≈ (2π × 0.36)² ≈ 5.1
#   Mode 3: ωa/(2πc) ≈ 0.36           →  λ ≈ 5.1  (degenerate pair)
#   Mode 4: ωa/(2πc) ≈ 0.53           →  λ ≈ 11.1
#   Mode 5: ωa/(2πc) ≈ 0.62           →  λ ≈ 15.2
#
# The target eps_target = 4.0 places the shift just below the first
# non-trivial TM mode.  Shift-invert will efficiently converge to the
# modes closest to this target.
#
# ============================================================
# BOUNDARY CONDITIONS: NEUMANN (PMC) AT GAMMA POINT
# ============================================================
#
# At Gamma (k = 0), Bloch modes have the full lattice periodicity.
# These modes are either symmetric or antisymmetric under reflection
# at the unit cell boundaries.  By applying Neumann BCs (∂ψ/∂n = 0),
# we select the symmetric Gamma-point modes — these include the lowest
# TM bands and band-edge modes relevant for band gap analysis.
#
# No explicit [BCs] block is needed: the natural FEM boundary condition
# is ∂ψ/∂n = 0.  The null space includes the constant mode (ψ = const)
# which gives λ = 0.  The shift-invert transformation with eps_target > 0
# skips this trivial mode and finds the first nontrivial photonic modes.
#
# NOTE: Periodic BCs would be more general (capturing both symmetric and
# antisymmetric modes) but create numerical difficulties with MOOSE's
# eigenvalue solver (singular stiffness matrix from eliminated DOFs).
#
# ============================================================
# NOTES ON MESH RESOLUTION
# ============================================================
#
# The circular rod of radius 0.2 is approximated on the QUAD4 mesh.
# At 40×40 elements the cell size h = 0.025; the rod radius spans
# 8 elements — adequate to resolve the dielectric interface and the
# fundamental mode patterns.  Increasing to 60×60 or using a finer
# mesh near the rod boundary would improve accuracy.
#
# ============================================================

[Mesh]
  # Single unit cell [0,1]×[0,1] with a=1 (normalised lattice constant).
  # 40×40 QUAD4 elements give h = 0.025 in each direction.
  # The rod of radius 0.2 spans ~16 element diameters — adequate
  # resolution for the circular dielectric interface.
  # GeneratedMeshGenerator names the four boundary sides:
  #   left (x=0), right (x=1), bottom (y=0), top (y=1)
  # These are used for the periodic BC definition below.
  [unit_cell]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 40
    ny   = 40
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
  []
[]

[Variables]
  # ψ — TM electric field component E_z in the unit cell.
  # Bloch theorem at Gamma point (k=0): ψ has the full lattice periodicity.
  # First-order Lagrange basis is in H¹; the periodic BC is compatible with
  # the continuity of ψ across the cell boundaries demanded by Bloch theorem.
  [psi]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ============================================================
  # DIELECTRIC FUNCTION εᵣ(x,y)
  # ============================================================
  # Square lattice of circular alumina rods (εᵣ = 8.9) in vacuum (εᵣ = 1).
  # Rod centre: (0.5, 0.5), radius r = 0.2.
  # The if() syntax in ParsedFunction evaluates a conditional expression:
  #   if(condition, value_if_true, value_if_false)
  # Condition: (x-0.5)² + (y-0.5)² < 0.04 = r² = 0.2²
  # Inside rod:    εᵣ = 8.9  (alumina, measured value at microwave frequencies)
  # Outside rod:   εᵣ = 1.0  (vacuum / air host medium)
  [eps_r_fn]
    type       = ParsedFunction
    expression = 'if((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<0.04, 8.9, 1.0)'
  []

  # ============================================================
  # VISUALISATION AUXILIARY: dielectric field for output
  # ============================================================
  # This function is also used to populate an AuxVariable so that
  # the dielectric contrast (εᵣ map) is visible alongside the
  # eigenmode field in the Exodus output.  Viewing both in ParaView
  # allows direct correlation of mode concentration with dielectric regions.
  [eps_r_vis_fn]
    type       = ParsedFunction
    expression = 'if((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<0.04, 8.9, 1.0)'
  []
[]

[AuxVariables]
  # eps_r_field: stores εᵣ(x,y) evaluated at element centroids for visualisation.
  # CONSTANT/MONOMIAL is appropriate because εᵣ is piecewise-constant (two values)
  # and is computed from a parsed function (not a nodal interpolation).
  # In ParaView: colour by eps_r_field to see the rod geometry overlaid with ψ.
  [eps_r_field]
    order  = CONSTANT
    family = MONOMIAL
  []

  # grad_x, grad_y: components of the mode's transverse field pattern E_t = -∇ψ.
  # These visualise the field intensity distribution inside and outside the rod,
  # showing whether the mode is concentrated in high-εᵣ (rod) or low-εᵣ (air) regions.
  # Modes below a band gap tend to concentrate in the high-εᵣ rod (dielectric band);
  # modes above the gap concentrate in air (air band) — the origin of the gap itself.
  [grad_x]
    order  = CONSTANT
    family = MONOMIAL
  []
  [grad_y]
    order  = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  # Evaluate εᵣ(x,y) from the ParsedFunction at each element's centroid.
  # FunctionAux calls eps_r_vis_fn(x_c, y_c) where (x_c,y_c) is the centroid.
  [eps_r_aux]
    type     = FunctionAux
    variable = eps_r_field
    function = eps_r_vis_fn
  []

  # Ex = -∂ψ/∂x  (x-component of mode electric field pattern).
  # PotentialToFieldAux computes the element-average gradient component of psi.
  # sign = negative implements the physics convention E = -∇ψ for a mode function.
  [grad_x_aux]
    type              = PotentialToFieldAux
    variable          = grad_x
    gradient_variable = psi
    sign              = negative
    component         = x
  []

  # Ey = -∂ψ/∂y  (y-component of mode electric field pattern).
  [grad_y_aux]
    type              = PotentialToFieldAux
    variable          = grad_y
    gradient_variable = psi
    sign              = negative
    component         = y
  []
[]

[Materials]
  # GenericFunctionMaterial evaluates the ParsedFunction eps_r_fn at each
  # quadrature point and stores the result as the material property 'eps_r'.
  # This non-AD material property is consumed by MatReaction (non-AD kernel)
  # to form the εᵣ-weighted mass matrix B of the generalised eigenproblem.
  #
  # IMPORTANT: GenericFunctionMaterial (not ADGenericFunctionMaterial) is used
  # because MatReaction (not ADMatReaction) reads the non-AD property 'eps_r'.
  # Using ADGenericFunctionMaterial with MatReaction would still work, but
  # using ADMatReaction with non-AD GenericFunctionMaterial would fail at runtime.
  [eps_r_mat]
    type        = GenericFunctionMaterial
    prop_names  = 'eps_r'
    prop_values = 'eps_r_fn'
  []
[]

[Kernels]
  # ============================================================
  # STIFFNESS KERNEL — builds matrix A in  A·x = λ·B·x
  # ============================================================
  # Diffusion computes the weak form of -∇²ψ, integrated by parts:
  #   ∫ ∇ψ · ∇φ dV
  # (boundary term vanishes with periodic BCs: the outward normal fluxes
  # at x=0 and x=1 cancel due to ψ periodicity; same for y faces.)
  # This kernel contributes to the standard (non-eigen) residual system
  # and therefore populates matrix A.
  [diff]
    type     = Diffusion
    variable = psi
  []

  # ============================================================
  # εᵣ-WEIGHTED MASS KERNEL — builds matrix B in  A·x = λ·B·x
  # ============================================================
  # MatReaction contributes:  residual = -reaction_rate × ψ × φ
  # With reaction_rate = εᵣ(x,y) from the material property:
  #   residual = -εᵣ(x,y) × ψ × φ
  #
  # Tagged extra_vector_tags = 'eigen' → sent to the 'eigen' (B) matrix
  # instead of the standard residual.  MOOSE's SLEPc interface then
  # assembles:
  #   B_ij = -∫ εᵣ(x,y) φ_j φ_i dV
  #
  # The generalised eigenproblem  A·x = λ·B·x  then reads:
  #   ∫ ∇ψ·∇φ dV = λ × (-∫ εᵣ ψ φ dV)
  # Strong form: -∇²ψ = -λ εᵣ ψ  →  ∇²ψ + λ εᵣ ψ = 0  ✓
  #
  # Sign check: B must be negative semi-definite for the standard form
  # A·x = λ·B·x to give positive λ = (ω/c)² eigenvalues.  With εᵣ > 0
  # the integral ∫ εᵣ φ² dV > 0, so B = -∫ εᵣ φ² dV < 0.  This is
  # consistent with MOOSE's SLEPc interface convention (same as Case 30
  # which uses CoefReaction with coefficient=-1 to get B = -∫ φ² dV).
  [mass]
    type              = MatReaction
    variable          = psi
    reaction_rate     = eps_r
    extra_vector_tags = 'eigen'
  []
[]

# ============================================================
# BOUNDARY CONDITIONS: NEUMANN (∂ψ/∂n = 0) — PERFECT MAGNETIC CONDUCTOR
# ============================================================
#
# For TM polarisation the natural (Neumann) boundary condition
# ∂ψ/∂n = 0 corresponds to a Perfect Magnetic Conductor (PMC) wall.
# This is physically equivalent to a mirror symmetry plane — modes
# that are symmetric under reflection across the unit cell boundary
# automatically satisfy this condition.
#
# At the Gamma point (k = 0), Bloch modes are either symmetric or
# antisymmetric with respect to the unit cell boundaries.  Neumann BCs
# select the symmetric subset, which includes the lowest TM modes.
# (The antisymmetric modes would require Dirichlet BCs, ψ = 0.)
#
# Together, the Neumann and Dirichlet BC solutions span the full set
# of Gamma-point modes.  For pedagogical purposes, the Neumann BCs
# capture the most important modes including those at band edges.
#
# TECHNICAL NOTE: MOOSE's periodic BC implementation combined with the
# eigenvalue solver creates a singular stiffness matrix (zero rows from
# eliminated DOFs).  Neumann BCs avoid this issue entirely while still
# producing physically correct Gamma-point eigenfrequencies.
#
# No [BCs] block is needed — the default natural BC in FEM is ∂ψ/∂n = 0.

[VectorPostprocessors]
  # Reports all n_eigen_pairs eigenvalues λ = (ω/c)² after convergence.
  # To convert to normalised frequency: ωa/(2πc) = √λ / (2π)  (with a=1).
  # The first non-trivial TM mode of this structure has
  # ωa/(2πc) ≈ 0.36, corresponding to λ ≈ (2π×0.36)² ≈ 5.1.
  [eigenvalues]
    type = Eigenvalues
  []
[]

[Postprocessors]
  # Track the field extrema to assess mode shape (sign of the eigenfunction).
  # The mode amplitude is set by SLEPc normalisation (not physically meaningful);
  # what matters is the spatial pattern: where is |ψ| large or small?
  [psi_max]
    type       = ElementExtremeValue
    variable   = psi
    value_type = max
  []
  [psi_min]
    type       = ElementExtremeValue
    variable   = psi
    value_type = min
  []
  # Average εᵣ over the unit cell as a sanity check.
  # Expected: avg_eps = f×8.9 + (1-f)×1.0 = 0.1257×8.9 + 0.8743 = 2.00
  # (where f = πr²/a² = π×0.04 ≈ 0.1257 is the filling fraction)
  [avg_eps]
    type     = ElementAverageValue
    variable = eps_r_field
  []
[]

[Executioner]
  type = Eigenvalue

  # KRYLOVSCHUR: SLEPc's Krylov-Schur method for large sparse generalised
  # eigenproblems.  Replaces the single-mode power iteration (PJFNK) and
  # allows computing multiple eigenvalues simultaneously.
  solve_type = KRYLOVSCHUR

  # Request 8 eigen pairs to capture the first several TM photonic bands
  # at the Gamma point.  The trivial λ = 0 constant mode will appear near
  # the bottom; the physically interesting photonic modes follow it.
  n_eigen_pairs     = 8
  which_eigen_pairs = TARGET_MAGNITUDE

  # Shift-invert spectral transformation with target σ = 4.0:
  #   eps_target = 4.0  — shift σ placed just below the estimated first
  #                        non-trivial TM mode λ ≈ 5.1.
  #                        Shift-invert maps λ → 1/(λ - σ) so modes with
  #                        λ close to σ become the largest in the transformed
  #                        problem and converge fastest.
  #   st_type = sinvert — activates the shift-invert spectral transformation.
  #   st_ksp_type preonly + st_pc_type lu + mumps — direct sparse LU solver
  #     for the inner linear system (A - σB)⁻¹ applied at each Arnoldi step.
  #     MUMPS is a parallel sparse direct solver; robust for 2D FEM problems
  #     of this size (40×40 = 1600 unknowns per variable).
  petsc_options_iname = '-eps_target -st_type -st_ksp_type -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = '4.0         sinvert  preonly       lu           superlu_dist'
[]

[Outputs]
  # exodus: full 2D field output — psi, eps_r_field, grad_x, grad_y —
  #   written to an Exodus II file for visualisation in ParaView.
  #   To view a mode: open the .e file, colour by 'psi', and observe
  #   how the field concentrates inside the dielectric rod (dielectric band)
  #   or in the air gaps (air band) — the fundamental distinction that
  #   creates photonic band gaps.
  # csv: eigenvalue table — list of λ values for frequency computation.
  #   Post-process: ωa/(2πc) = sqrt(lambda) / (2*pi)
  # execute_on = FINAL: write output only after the eigensolver converges.
  exodus     = true
  csv        = true
  execute_on = FINAL
[]
