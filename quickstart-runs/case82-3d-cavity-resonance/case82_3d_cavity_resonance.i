# ============================================================
# Case 82: 3D Rectangular Cavity Resonant Modes (TM-like scalar modes)
# Extends Case 30 (2D waveguide cutoff) to a full 3D closed cavity.
#
# Reference: Pozar, Microwave Engineering, 4th ed. (2012), Sec. 6.7
#            Collin, Foundations for Microwave Engineering, 2nd ed. (1992), Ch. 7
#
# PHYSICS — THE 3D CAVITY EIGENVALUE PROBLEM
# -------------------------------------------
# A perfectly conducting (PEC) rectangular cavity occupies the volume
#   Ω = [0,a] × [0,b] × [0,d]
# with a = 3, b = 2, d = 1 (all dimensions in metres).
#
# Inside a PEC cavity each Cartesian field component ψ satisfies the
# scalar Helmholtz eigenvalue equation:
#
#   ∇²ψ + k²ψ = 0   in Ω
#   ψ = 0            on ∂Ω  (PEC wall condition: tangential E = 0)
#
# Separating variables in Cartesian coordinates gives:
#   ψ(x,y,z) = sin(mπx/a) · sin(nπy/b) · sin(pπz/d)
#
# The resonant eigenvalues are:
#   k²(m,n,p) = (mπ/a)² + (nπ/b)² + (pπ/d)²,   m,n,p = 1, 2, 3, …
#
# Note: This problem is separable and the modes are products of 1D sine
# solutions along each axis — an elegant generalisation of the 2D
# waveguide case (Case 30 had ψ = 0 only on four faces; here all six
# faces carry the PEC condition).
#
# ANALYTICAL EIGENVALUES FOR a=3, b=2, d=1
# ------------------------------------------
# Constants:
#   (π/3)²  =  1.0966   (π/2)²  =  2.4674   (π/1)²  =  9.8696
#   (2π/3)² =  4.3864   (2π/2)² =  9.8696   (2π/1)² = 39.478
#   (3π/3)² =  9.8696
#
# Mode    m  n  p   (mπ/a)²   (nπ/b)²   (pπ/d)²   k²_exact
# ---------------------------------------------------------------
# TM₁₁₁  1  1  1   1.0966    2.4674    9.8696    13.434   ← target
# TM₂₁₁  2  1  1   4.3864    2.4674    9.8696    16.723
# TM₁₂₁  1  2  1   1.0966    9.8696    9.8696    20.836
# TM₃₁₁  3  1  1   9.8696    2.4674    9.8696    22.207
# TM₂₂₁  2  2  1   4.3864    9.8696    9.8696    24.126
# TM₁₁₂  1  1  2   1.0966    2.4674   39.478     43.042
# ---------------------------------------------------------------
#
# HOW THIS EXTENDS CASE 30 (2D WAVEGUIDE)
# ----------------------------------------
# Case 30: 2D domain [0,a]×[0,b] — TE/TM cutoff.
#   • 4 PEC faces (left, right, top, bottom).
#   • 2 index integers (m, n).
#   • mesh: QUAD4 elements.
#
# Case 82: 3D domain [0,a]×[0,b]×[0,d] — full resonant cavity.
#   • 6 PEC faces (left, right, top, bottom, front, back).
#   • 3 index integers (m, n, p) — extra mode family in z.
#   • mesh: HEX8 elements (dim=3 in GeneratedMeshGenerator).
#   • PotentialToFieldAux now computes all three gradient components
#     (Ex, Ey, Ez) instead of just two.
#
# The MOOSE formulation is otherwise identical:
#   [Diffusion]   → stiffness matrix A = ∫ ∇ψ·∇φ dV
#   [CoefReaction] (extra_vector_tags = 'eigen') → mass matrix B = ∫ ψ·φ dV
#   [DirichletBC] + [EigenDirichletBC]  →  homogeneous BC on A and B
#   Eigenvalue executioner with KRYLOVSCHUR + shift-invert at σ = 13.0
#   (just below k²(1,1,1) = 13.434 so the solver converges to the
#    cavity modes ascending from the fundamental).
#
# Expected MOOSE eigenvalues (15×10×5 HEX8 mesh, FEM discretisation error
# ~ O(h²) for LAGRANGE/FIRST):
#   λ₁ ≈ 13.5,  λ₂ ≈ 16.8,  λ₃ ≈ 20.9,  λ₄ ≈ 22.3,  λ₅ ≈ 24.2
# ============================================================

[Mesh]
  # GeneratedMeshGenerator with dim=3 produces a structured HEX8 mesh
  # on the rectangular parallelepiped [0,3]×[0,2]×[0,1].
  # Resolution 15×10×5 = 750 elements — sufficient for the first 5 modes
  # and small enough to run inside Docker without memory pressure.
  # The six boundary names created automatically by GeneratedMeshGenerator
  # for a 3D box are: left, right, bottom, top, back, front
  #   left   → x = 0    right  → x = 3
  #   bottom → y = 0    top    → y = 2
  #   back   → z = 0    front  → z = 1
  [gmg]
    type = GeneratedMeshGenerator
    dim  = 3
    nx   = 15
    ny   = 10
    nz   = 5
    xmin = 0
    xmax = 3
    ymin = 0
    ymax = 2
    zmin = 0
    zmax = 1
  []
[]

[Variables]
  # ψ — scalar potential whose nodal values form the eigenvector.
  # Physically ψ represents the dominant Cartesian field component of
  # the TM-like cavity mode (e.g. Ez for TM modes excited along z).
  # First-order Lagrange basis is in H¹, satisfying the continuity
  # requirement for the Helmholtz operator; the Dirichlet condition
  # ψ = 0 on ∂Ω is imposed strongly.
  [psi]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  # Transverse (and longitudinal) electric field components derived from
  #   E = -∇ψ
  # Stored as CONSTANT/MONOMIAL (piecewise-constant per element) because
  # they are computed as the element-average gradient of the nodal
  # LAGRANGE variable — exactly what PotentialToFieldAux produces.
  # These are useful for visualising the mode field patterns in ParaView.
  [Ex]
    order  = CONSTANT
    family = MONOMIAL
  []
  [Ey]
    order  = CONSTANT
    family = MONOMIAL
  []
  # Ez is the third gradient component — new relative to Case 30 (2D).
  # For a TM mode with ψ = Ez, this component gives the axial field
  # itself (as a piecewise element average of −∂ψ/∂z).
  [Ez]
    order  = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  # -----------------------------------------------------------------------
  # Stiffness kernel — weak form of −∇²ψ.
  # Integration by parts: ∫ ∇ψ·∇φ dV  (boundary term vanishes because
  # ψ = 0 on ∂Ω).  This kernel contributes to the standard residual
  # system and therefore builds matrix A in the generalised eigenproblem
  #   A·x = λ·B·x
  # -----------------------------------------------------------------------
  [diff]
    type     = Diffusion
    variable = psi
  []

  # -----------------------------------------------------------------------
  # Mass kernel — weak form of −k²ψ.
  # CoefReaction contributes −coefficient · ∫ ψ·φ dV to the system
  # tagged 'eigen'.  Setting coefficient = -1 means the 'eigen' matrix
  # receives  +∫ ψ·φ dV, which is the positive-definite mass matrix B.
  # SLEPc then seeks the smallest λ = k² satisfying A·x = λ·B·x,
  # i.e. the resonant eigenvalues of the cavity in ascending order.
  # -----------------------------------------------------------------------
  [coeff]
    type              = CoefReaction
    variable          = psi
    coefficient       = -1
    extra_vector_tags = 'eigen'
  []
[]

[AuxKernels]
  # -----------------------------------------------------------------------
  # Electric field reconstruction.
  # PotentialToFieldAux evaluates the element-average gradient of psi
  # and stores the requested Cartesian component in the aux variable.
  # sign = negative implements E = -∇ψ (gradient of a scalar potential
  # with a minus sign, standard in EM theory).
  # -----------------------------------------------------------------------

  # Ex = -∂ψ/∂x  (x-component of the transverse electric field pattern)
  [Ex_aux]
    type              = PotentialToFieldAux
    variable          = Ex
    gradient_variable = psi
    sign              = negative
    component         = x
  []

  # Ey = -∂ψ/∂y  (y-component of the transverse electric field pattern)
  [Ey_aux]
    type              = PotentialToFieldAux
    variable          = Ey
    gradient_variable = psi
    sign              = negative
    component         = y
  []

  # Ez = -∂ψ/∂z  (z-component; new in 3D — absent in Case 30).
  # For a pure TM₁₁₁ mode ψ ∝ sin(πx/3)sin(πy/2)sin(πz/1), this
  # component is largest near z = 0.5 where cos(πz/1) is maximal.
  [Ez_aux]
    type              = PotentialToFieldAux
    variable          = Ez
    gradient_variable = psi
    sign              = negative
    component         = z
  []
[]

[BCs]
  # -----------------------------------------------------------------------
  # PEC boundary conditions — ψ = 0 on all six cavity walls.
  #
  # In a 3D GeneratedMesh the boundary set names are:
  #   left, right   → faces perpendicular to x
  #   bottom, top   → faces perpendicular to y
  #   back, front   → faces perpendicular to z
  #
  # Two BC objects are required (as in Case 30):
  #   DirichletBC       → enforces ψ = 0 in the stiffness system (A)
  #   EigenDirichletBC  → enforces ψ = 0 in the mass system (B)
  # Without EigenDirichletBC the boundary nodes would contribute to B
  # with non-zero entries, producing spurious modes near λ = 0.
  # -----------------------------------------------------------------------

  # Standard Dirichlet — enforces zero field on all six PEC faces
  # in the primary (stiffness) system.
  [pec_walls]
    type     = DirichletBC
    variable = psi
    boundary = 'left right top bottom front back'
    value    = 0
  []

  # Eigen Dirichlet — simultaneously enforces the same constraint
  # in the eigen (mass) system so that A·x = λ·B·x sees zero BCs
  # on both sides of the generalised eigenproblem.
  [eigen_pec_walls]
    type     = EigenDirichletBC
    variable = psi
    boundary = 'left right top bottom front back'
  []
[]

[VectorPostprocessors]
  # Reports all computed eigenvalues (k² values) as a CSV column after
  # the eigensolver converges.  The first entry is the fundamental mode
  # TM₁₁₁ with k²_exact = 13.434.  Compare with analytical values in
  # the header table to assess discretisation accuracy.
  [eigenvalues]
    type = Eigenvalues
  []
[]

[Executioner]
  type = Eigenvalue

  # KRYLOVSCHUR: SLEPc's Krylov-Schur algorithm — the recommended method
  # for computing multiple eigenvalues of a large sparse eigenproblem.
  # (Default PJFNK would use a nonlinear power-iteration scheme capable
  # of finding only the single dominant eigenvalue.)
  solve_type = KRYLOVSCHUR

  # Request 6 eigen pairs so that the five physically distinct cavity modes
  # listed in the header are all captured:
  #   TM₁₁₁ (k²=13.434), TM₂₁₁ (16.723), TM₁₂₁ (20.836),
  #   TM₃₁₁ (22.207), TM₂₂₁ (24.126) — plus one extra for robustness.
  n_eigen_pairs     = 6
  which_eigen_pairs = TARGET_MAGNITUDE

  # Shift-invert spectral transformation:
  #   eps_target = 13.0  — shift σ placed just below the fundamental
  #                         cavity eigenvalue k²(1,1,1) = 13.434.
  #                         Shift-invert maps λ → 1/(λ-σ) so the
  #                         eigenvalues closest to σ become the largest
  #                         in the transformed problem and converge first.
  #   st_type = sinvert  — activates the shift-invert transformation.
  #   st_ksp_type/st_pc_type/st_pc_factor_mat_solver_type — configure the
  #     inner linear solver used for the (A - σB)⁻¹ action: direct LU
  #     factorisation via MUMPS (robust for 3D FEM problems).
  petsc_options_iname = '-eps_target -st_type -st_ksp_type -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = '13.0        sinvert  preonly       lu           mumps'
[]

[Outputs]
  # exodus: full 3D field output — ψ, Ex, Ey, Ez written to an Exodus II
  #   file for 3D mode-shape visualisation in ParaView (volume render or
  #   isosurface at ψ = ±0.5 to see the standing-wave nodal planes).
  # csv: eigenvalue table for direct comparison with analytical k² values.
  # execute_on = FINAL: write output only after the eigensolver converges
  #   (avoids partial output during SLEPc iterations).
  exodus     = true
  csv        = true
  execute_on = FINAL
[]
