# ============================================================
# Case 89: Dielectric Slab Waveguide — TE Guided Mode Eigenvalues
# Saleh & Teich, "Fundamentals of Photonics", 3rd Ed. (2019), Ch. 8
# Pozar, "Microwave Engineering", 4th Ed. (2012), Sec. 3.1
# Yariv & Yeh, "Photonics: Optical Electronics in Modern
#   Communications", 6th Ed. (2007), Ch. 3
#
# MIT 6.635 connection: Lectures 3-5 (guided waves, dielectric waveguides,
#   dispersion relations, eigenmode analysis)
#
# ============================================================
# PHYSICAL BACKGROUND
# ============================================================
#
# A symmetric planar dielectric waveguide consists of three layers:
#
#   Cladding:  n2 (lower refractive index)  for |x| > d/2
#   Core:      n1 (higher refractive index) for |x| < d/2
#   Cladding:  n2                            for |x| > d/2
#
# Light is guided by total internal reflection (TIR) when the
# propagation constant beta satisfies:
#
#   n2 k0 < beta < n1 k0   (guidance condition)
#
# Below n2*k0: radiation modes (light leaks into cladding)
# Above n1*k0: physically impossible (exceeds plane-wave limit in core)
# Between:     guided (bound) modes — evanescent in cladding, oscillatory in core
#
# ============================================================
# TRANSVERSE EIGENVALUE EQUATION
# ============================================================
#
# For TE modes (electric field polarised in y-direction, propagation in z),
# the transverse field psi(x) satisfies:
#
#   d^2 psi/dx^2 + (n(x)^2 k0^2 - beta^2) psi = 0
#
# where n(x) = n1 for |x| < d/2 (core), n(x) = n2 otherwise (cladding).
#
# This is a Sturm-Liouville eigenvalue problem. Defining lambda = beta^2:
#
#   d^2 psi/dx^2 + n^2(x) k0^2 psi = lambda psi
#
# With boundary conditions:
#   psi(x) -> 0 as x -> +/- infinity (evanescent decay in cladding)
#
# MOOSE solves this as a GENERALISED eigenvalue problem A psi = lambda B psi:
#
#   A_ij = integral( d phi_j/dx * d phi_i/dx dx )    (Diffusion, stiffness)
#         - integral( n^2(x) k0^2 * phi_j * phi_i dx ) (MatReaction with rate = n^2 k0^2)
#
# Wait — careful! Let's derive properly.
#
# Weak form: multiply d^2psi/dx^2 + n^2 k0^2 psi = lambda psi by test phi,
# integrate by parts (boundary term at x = +/-W/2 vanishes since psi = 0):
#
#   -integral( dpsi/dx * dphi/dx dx ) + integral( n^2 k0^2 psi phi dx )
#   = lambda * integral( psi phi dx )
#
# MOOSE kernel residuals:
#   Diffusion:    +integral( dpsi/dx * dphi/dx dx )  [strong form: -d^2psi/dx^2]
#   MatReaction(rate): -rate * integral( psi phi dx ) [strong form: -rate * psi]
#
# To get the A matrix = LHS of weak form = -integral(grad grad) + integral(n^2 k0^2):
#
#   Standard residual = Diffusion + MatReaction(rate = -n^2 k0^2):
#     = +integral(grad grad) - (-n^2 k0^2) * integral(psi phi)
#     = +integral(grad grad) + n^2 k0^2 * integral(psi phi)
#     This gives A representing: -d^2psi/dx^2 - n^2 k0^2 psi = RHS
#     Strong form: d^2psi/dx^2 + n^2 k0^2 psi = -RHS ... not quite right.
#
# Better: use MatReaction with rate = +n^2 k0^2 for the A matrix term:
#   Residual_A = Diffusion + MatReaction(rate = +n^2 k0^2)
#              = +integral(grad grad) - n^2 k0^2 * integral(psi phi)
#   Strong form from residual_A = 0:
#     -nabla^2 psi - n^2 k0^2 psi = 0
#     <=>  nabla^2 psi + n^2 k0^2 psi = 0  ... still Helmholtz, not eigenvalue.
#
# For the generalised eigenvalue A psi = lambda B psi, we need to SEPARATE
# the lambda term onto the B (eigen) matrix side. The correct formulation is:
#
#   A matrix (standard residual): Diffusion - MatReaction(rate = n^2 k0^2)
#     A_ij = integral(grad phi_j . grad phi_i) - n^2 k0^2 integral(phi_j phi_i)
#     Strong form of A: -d^2psi/dx^2 - n^2 k0^2 psi
#
#   B matrix (eigen-tagged): CoefReaction with coefficient = -1
#     B_ij = -(-1) * integral(phi_j phi_i) = +integral(phi_j phi_i)
#     Wait: CoefReaction residual = +coef * integral(psi phi)
#     With coef = -1: residual_B = -integral(psi phi)
#     So B_ij = -integral(phi_j phi_i)  (negative mass matrix)
#
#   Eigenproblem: A psi = lambda B psi
#   [+integral(grad grad) - n^2 k0^2 integral(psi phi)] psi
#   = lambda * [-integral(psi phi)] psi
#
#   Strong form:
#   (-d^2psi/dx^2 - n^2 k0^2 psi) = lambda * (-psi)
#   => d^2psi/dx^2 + n^2 k0^2 psi = lambda psi   checkmark
#
# Therefore with eigenvalue lambda = beta^2:
#   The A matrix contains:
#     [diff]     Diffusion:         +integral(grad psi . grad phi)
#     [n2k0sq]   MatReaction:       -n^2 k0^2 * integral(psi phi)   [NOT eigen-tagged]
#   The B matrix contains:
#     [mass]     CoefReaction(coefficient=-1, tagged 'eigen'):
#                -1 * integral(psi phi) => B_ij = -integral(phi_j phi_i)
#
# Note: Using MatReaction (non-AD) paired with GenericFunctionMaterial (non-AD),
# following the established Case 81 pattern for eigenvalue problems.
# ADMatReaction + non-AD material property would fail at runtime.
#
# ============================================================
# PARAMETERS AND EXPECTED EIGENVALUES
# ============================================================
#
# Core:     n1 = 1.5,  d = 1.0 m  (core thickness)
# Cladding: n2 = 1.0
# lambda_0 = 1.0 m,  k0 = 2pi rad/m
#
# V-number (normalised frequency):
#   V = k0 * (d/2) * sqrt(n1^2 - n2^2)
#     = 2pi * 0.5 * sqrt(1.5^2 - 1.0^2)
#     = pi * sqrt(2.25 - 1.0)
#     = pi * sqrt(1.25)
#     = pi * 1.11803...
#     = 3.5124...
#
# Number of guided modes: N = ceil(V/pi) = ceil(1.118) = 2
# The waveguide supports exactly 2 TE guided modes (TE0 and TE1).
# (Each additional mode requires V to increase by pi, so V > pi for mode 1,
#  V > 2pi for mode 2, etc. Here 2pi < V < 3pi, so modes TE0 and TE1 exist.)
#
# Wait, more precisely: modes exist for V > m*pi/2, m = 1, 2, 3...
# Cut-off of TE_m: V_c(m) = m*pi/2
# TE0: V > 0 (always guided) -> TE0 exists
# TE1: V > pi/2 ~ 1.571 -> V = 3.512 > 1.571 -> TE1 exists
# TE2: V > pi ~ 3.141 -> V = 3.512 > pi -> TE2 also exists!
# TE3: V > 3pi/2 ~ 4.712 -> V = 3.512 < 4.712 -> TE3 does not exist
#
# The waveguide actually supports 3 guided TE modes (TE0, TE1, TE2).
#
# Numerical ranges:
#   k0^2 = (2pi)^2 = 39.4784 m^-2
#   n1^2 k0^2 = (1.5)^2 * (2pi)^2 = 2.25 * 39.4784 = 88.826 m^-2
#   n2^2 k0^2 = (1.0)^2 * (2pi)^2 = 1.0  * 39.4784 = 39.4784 m^-2
#
# Guidance condition in terms of lambda = beta^2:
#   n2^2 k0^2 < lambda < n1^2 k0^2
#   39.4784 < lambda < 88.826
#
# Characteristic equations for TE modes (symmetric slab, derived from
# matching conditions at x = d/2):
#   Even modes (TE0, TE2, ...): kappa * tan(kappa * d/2) = gamma
#   Odd modes  (TE1, TE3, ...): kappa * cot(kappa * d/2) = -gamma  (or -kappa*tan = gamma)
# where:
#   kappa = sqrt(n1^2 k0^2 - beta^2) = sqrt(n1^2 k0^2 - lambda)  (oscillatory in core)
#   gamma = sqrt(beta^2 - n2^2 k0^2) = sqrt(lambda - n2^2 k0^2)  (decaying in cladding)
#
# The FEM eigenvalue MOOSE computes converges to these analytically derived beta^2 values.
# Approximate estimates from graphical solution of the characteristic equations:
#   TE0 (even, fundamental):  lambda ~ 80  (kappa ~ 2.98, gamma ~ 6.54)
#   TE1 (odd,  first):        lambda ~ 60  (kappa ~ 5.29, gamma ~ 4.61)
#   TE2 (even, second):       lambda ~ 43  (kappa ~ 6.97, gamma ~ 1.84)
#
# These are rough estimates; the FEM solution gives exact values.
# Target for shift-invert: eps_target = 39.5 (just above n2^2 k0^2 = 39.48)
# This ensures the solver finds eigenvalues just above the radiation continuum,
# starting from the highest-order guided mode (TE2) and working up to TE0.
#
# ============================================================
# DOMAIN AND BOUNDARY CONDITIONS
# ============================================================
#
# Domain: 1D, [-W/2, W/2] = [-5, 5] m
# The evanescent decay length in the cladding is:
#   delta = 1/gamma = 1/sqrt(lambda - n2^2 k0^2)
# For TE0 (lambda ~ 80): delta ~ 1/sqrt(80-39.48) ~ 1/6.37 ~ 0.157 m
# For TE2 (lambda ~ 43): delta ~ 1/sqrt(43-39.48) ~ 1/1.88 ~ 0.53 m
#
# At x = +/-5 m (domain boundary), the evanescent field for all modes:
#   psi ~ exp(-gamma*(5 - d/2)) = exp(-gamma*4.5)
# For TE2 (weakest decay): exp(-1.88*4.5) = exp(-8.46) ~ 2e-4 (very small)
# Therefore the Dirichlet BC psi=0 at x=+/-5 is an excellent approximation.
#
# Boundary conditions: DirichletBC psi=0 at x=-5 and x=+5.
# In the standard Dirichlet form: psi = 0 on left and right boundaries.
# EigenDirichletBC: also enforces psi=0 in the B (mass) matrix to avoid
# spurious modes from the boundary nodes.
#
# ============================================================

[Mesh]
  # 1D domain from x = -5 m to x = +5 m (total width W = 10 m).
  # Core occupies |x| < d/2 = 0.5 m; cladding occupies 0.5 < |x| < 5 m.
  #
  # Resolution:
  #   nx = 400: dx = 10/400 = 0.025 m
  #   Core (|x| < 0.5): 40 elements spanning the core half-width
  #   Cladding half-space: (5-0.5)/0.025 = 180 elements for evanescent decay
  #
  # The free-space wavelength lambda_0 = 1.0 m contains 40 elements —
  # adequate for the smooth mode shape (no singularity). The mode field
  # psi(x) in the core is a cosine/sine function with period ~= 2*pi/kappa;
  # for TE0 kappa ~ 2.98 rad/m -> period ~ 2.1 m -> ~84 elements per period.
  # This is more than sufficient for FIRST-order Lagrange elements.
  #
  # GeneratedMeshGenerator in 1D names two boundaries: left (x=-5) and right (x=5).
  [core_cladding]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = 400
    xmin = -5
    xmax = 5
  []
[]

[Variables]
  # psi — transverse mode field for TE polarisation.
  # In a planar waveguide with propagation in z, psi(x) represents:
  #   E_y(x) for TE modes (electric field y-component as function of x)
  #
  # First-order Lagrange basis: C^0 continuous across the core-cladding
  # interface at x = +/-0.5 m, which is the correct interface condition
  # for TE modes (E_y is continuous across dielectric interfaces).
  #
  # The mode field psi is normalised by SLEPc to unit norm internally;
  # the sign of the eigenfunction is arbitrary.
  [psi]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # n2k0sq_fn: n^2(x) * k0^2 — the position-dependent coefficient
  # for the MatReaction kernel in the A matrix.
  #
  # n(x) = n1 = 1.5 for |x| < d/2 = 0.5 m  (core)
  # n(x) = n2 = 1.0 for |x| >= 0.5 m        (cladding)
  #
  # Numerical values:
  #   Core:     n1^2 k0^2 = (1.5)^2 * (2pi)^2 = 2.25 * 39.47841760 = 88.82643961 m^-2
  #   Cladding: n2^2 k0^2 = (1.0)^2 * (2pi)^2 = 1.00 * 39.47841760 = 39.47841760 m^-2
  #
  # The if() expression evaluates the piecewise-constant profile:
  #   If |x| < 0.5: return core value (88.826...)
  #   Otherwise:    return cladding value (39.478...)
  #
  # This is the refractive-index step function that creates the waveguiding
  # mechanism. The discontinuity at x = +/-0.5 m models the abrupt
  # core-cladding interface (step-index waveguide).
  # ------------------------------------------------------------------
  [n2k0sq_fn]
    type       = ParsedFunction
    expression = 'if(abs(x) < 0.5, 88.82643961, 39.47841760)'
  []

  # ------------------------------------------------------------------
  # eps_r_vis_fn: n^2(x) for visualisation only.
  # Stored as an AuxVariable (CONSTANT MONOMIAL) so ParaView can display
  # the refractive-index profile alongside the mode field psi.
  # Shows the step from n=1 in cladding to n=1.5 in core as a clear
  # two-level discontinuous field.
  # ------------------------------------------------------------------
  [eps_r_vis_fn]
    type       = ParsedFunction
    expression = 'if(abs(x) < 0.5, 2.25, 1.0)'
  []
[]

[AuxVariables]
  # n_sq_field: stores n^2(x) evaluated at element centroids for visualisation.
  # CONSTANT/MONOMIAL: appropriate for a piecewise-constant profile —
  # represents the average refractive index in each element.
  # In ParaView, colour by n_sq_field to see the two-layer structure
  # (n^2 = 2.25 inside core, n^2 = 1.0 in cladding).
  [n_sq_field]
    order  = CONSTANT
    family = MONOMIAL
  []

  # dpsi_dx: transverse gradient of the mode field, d psi/dx.
  # For TE modes, this is proportional to the magnetic field component:
  #   H_z ~ dpsi/dx (in the appropriate normalisation).
  # Mode shape: in the core dpsi/dx varies as a sine/cosine; in the
  # cladding dpsi/dx decays exponentially — showing the guided nature.
  [dpsi_dx]
    order  = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  # Evaluate n^2(x) from the ParsedFunction at each element's centroid.
  # FunctionAux calls eps_r_vis_fn(x_c) where x_c is the element centroid.
  # The piecewise-constant nature is correctly captured: core elements
  # (centroid inside |x| < 0.5) get n^2 = 2.25; cladding elements get 1.0.
  [n_sq_aux]
    type     = FunctionAux
    variable = n_sq_field
    function = eps_r_vis_fn
  []

  # d psi/dx: transverse derivative of the mode field.
  # PotentialToFieldAux computes the element-average gradient component.
  # sign = positive: returns +d psi/dx (not the negative E = -grad psi convention,
  # since psi here IS the mode field, not a potential).
  # For TE0 (even, symmetric): dpsi/dx is an odd function of x —
  #   positive slope for x < 0, negative for x > 0.
  # For TE1 (odd, antisymmetric): dpsi/dx is an even function of x.
  [dpsi_dx_aux]
    type              = PotentialToFieldAux
    variable          = dpsi_dx
    gradient_variable = psi
    sign              = positive
    component         = x
  []
[]

[Materials]
  # GenericFunctionMaterial (non-AD): wraps n2k0sq_fn as the non-AD material
  # property 'n2k0sq'. This is consumed by MatReaction (non-AD kernel).
  #
  # IMPORTANT: Use GenericFunctionMaterial (not ADGenericFunctionMaterial) because:
  #   - MatReaction reads non-AD properties.
  #   - For eigenvalue problems, AD is not required (the Jacobian is assembled
  #     exactly from the linear structure; no Newton iteration is needed for psi).
  #   - This is the same pattern as Case 81 (photonic crystal eigenvalue):
  #     GenericFunctionMaterial + MatReaction (non-AD), NOT ADMatReaction.
  #   - Mixing: ADMatReaction + non-AD material would cause a runtime type error.
  [n2k0sq_mat]
    type        = GenericFunctionMaterial
    prop_names  = 'n2k0sq'
    prop_values = 'n2k0sq_fn'
  []
[]

[Kernels]
  # ==================================================================
  # A MATRIX — Standard (non-eigen-tagged) residual
  # ==================================================================
  #
  # Weak form of the operator (d^2/dx^2 + n^2 k0^2) on the A side:
  #   A_ij = -integral( dphi_j/dx * dphi_i/dx dx ) + integral( n^2 k0^2 phi_j phi_i dx )
  #        = -stiffness + n^2 k0^2 mass
  #
  # The full eigenproblem is:
  #   A psi = lambda B psi
  #   [-stiffness + n^2 k0^2 mass] psi = lambda [-mass] psi
  #   => stiffness psi - n^2 k0^2 mass psi = -lambda mass psi
  #   => stiffness psi = (n^2 k0^2 - lambda) mass psi
  #   Strong form: -d^2 psi/dx^2 = (n^2 k0^2 - lambda) psi
  #   <=> d^2 psi/dx^2 + n^2 k0^2 psi = lambda psi   checkmark
  #
  # ==================================================================

  # Stiffness kernel — Diffusion builds +integral( dpsi/dx * dphi/dx dx ).
  # In standard (non-eigen) residual; contributes to A matrix.
  # Strong form: -d^2 psi/dx^2 appears on the left side of A psi = lambda B psi.
  [diff]
    type     = Diffusion
    variable = psi
  []

  # n^2 k0^2 mass kernel — MatReaction with reaction_rate = n^2(x) k0^2.
  # MatReaction residual = -reaction_rate * integral(psi phi) = -n^2 k0^2 integral(psi phi).
  # Combined with Diffusion in the standard residual:
  #   residual_A = +integral(dpsi/dx dphi/dx dx) - n^2 k0^2 integral(psi phi dx)
  # This represents A psi with:
  #   A_ij = integral(dphi_j/dx dphi_i/dx) - n^2 k0^2 integral(phi_j phi_i)
  #
  # Note: NOT eigen-tagged — goes to the A (standard) system, not the B (mass) system.
  # Note: MatReaction (non-AD) used, consistent with GenericFunctionMaterial (non-AD).
  # Note: reaction_rate parameter (current MOOSE version name, not 'mob_name').
  [n2k0sq_kernel]
    type          = MatReaction
    variable      = psi
    reaction_rate = n2k0sq
  []

  # ==================================================================
  # B MATRIX — Eigen-tagged residual
  # ==================================================================
  #
  # The B (mass) matrix for the generalised eigenproblem A psi = lambda B psi.
  # We need B_ij = -integral(phi_j phi_i dx) (negative mass matrix).
  #
  # CoefReaction with coefficient = -1:
  #   residual_B = (-1) * integral(psi phi dx) [CoefReaction sign is +coef]
  #   Wait: CoefReaction residual = +coefficient * integral(psi phi)
  #   With coefficient = -1: residual_B = -integral(psi phi)
  #   => B_ij = -integral(phi_j phi_i)
  #
  # MOOSE eigensolver: A psi = lambda B psi
  # With our B: A psi = lambda * (-mass) psi
  # => A psi + lambda * mass psi = 0
  # => [stiffness - n^2 k0^2 mass] psi = -lambda * mass psi
  # => d^2 psi/dx^2 + n^2 k0^2 psi = lambda psi   checkmark
  #
  # The negative sign on the B matrix means the physically relevant
  # eigenvalues lambda = beta^2 are POSITIVE. SLEPc with TARGET_MAGNITUDE
  # and shift-invert near lambda ~ 60 will find these positive guided-mode
  # eigenvalues efficiently.
  # ==================================================================
  [mass]
    type              = CoefReaction
    variable          = psi
    coefficient       = -1
    extra_vector_tags = 'eigen'
  []
[]

[BCs]
  # ==================================================================
  # Dirichlet BCs at the domain boundaries x = +/-5 m.
  #
  # Physical justification:
  # The guided mode field decays exponentially in the cladding:
  #   psi(x) ~ exp(-gamma * (|x| - d/2))  for |x| > d/2
  # At x = +/-5 m (4.5 m into the cladding from the core edge):
  #   Most weakly decaying mode (TE2, gamma ~ 1.88): psi ~ exp(-8.46) ~ 2e-4
  # This is negligible, so psi = 0 is an excellent approximation.
  #
  # Two BC objects required for the eigenvalue problem:
  #   DirichletBC      -> enforces psi = 0 in the A (stiffness) matrix.
  #   EigenDirichletBC -> enforces psi = 0 in the B (mass) matrix.
  # Without EigenDirichletBC, boundary nodes contribute non-trivially
  # to the B matrix, producing spurious eigenvalues near the BC nodes.
  # ==================================================================

  # Standard Dirichlet — psi = 0 on both domain boundaries.
  # Applied to the standard (A) system.
  [pec_walls]
    type     = DirichletBC
    variable = psi
    boundary = 'left right'
    value    = 0
  []

  # Eigen Dirichlet — simultaneously enforces psi = 0 in the B system.
  # This is essential for the KRYLOVSCHUR eigensolver to avoid spurious
  # modes associated with the boundary nodes.
  # Same pattern as Case 82 (3D cavity resonance).
  [eigen_pec_walls]
    type     = EigenDirichletBC
    variable = psi
    boundary = 'left right'
  []
[]

[VectorPostprocessors]
  # Reports all n_eigen_pairs computed eigenvalues lambda = beta^2.
  # To convert to propagation constant: beta = sqrt(lambda) [rad/m]
  # To convert to effective index: n_eff = beta/k0 = sqrt(lambda)/(2pi)
  #
  # Expected eigenvalues (guided modes, n2^2 k0^2 < lambda < n1^2 k0^2):
  #   39.478 < lambda < 88.826
  # Approximate values from characteristic equations:
  #   TE0 (fundamental, even):   lambda ~ 80  (n_eff ~ 1.42)
  #   TE1 (first order, odd):    lambda ~ 60  (n_eff ~ 1.23)
  #   TE2 (second order, even):  lambda ~ 43  (n_eff ~ 1.04)
  # These values are rough estimates; the FEM solution gives precise values.
  [eigenvalues]
    type = Eigenvalues
  []
[]

[Postprocessors]
  # ==================================================================
  # MODE SHAPE DIAGNOSTICS
  # ==================================================================
  # PointValue postprocessors sample psi at key locations to assess
  # the mode shape and verify correct mode identification.
  #
  # Mode identification from psi values:
  #   TE0 (even/symmetric):
  #     psi(0) = max (peak at core center)
  #     psi(0.5) > 0 (same sign as center, smooth decay across interface)
  #     psi(2) ~ 0 (evanescent decay in cladding)
  #
  #   TE1 (odd/antisymmetric):
  #     psi(0) = 0 (node at center by antisymmetry)
  #     psi(0.25) > 0, psi(-0.25) < 0 (or vice versa)
  #     psi(0.5) = interface value
  #
  #   TE2 (even with one interior node):
  #     psi(0) = local maximum OR node (depending on normalisation)
  #     Has one additional node inside the core
  # ==================================================================

  # At the core center (x=0): maximum for even modes, zero for odd modes.
  [psi_at_center]
    type     = PointValue
    variable = psi
    point    = '0 0 0'
  []

  # At the core-cladding interface (x=0.5): continuity of psi here.
  # The mode field is continuous across the dielectric step.
  [psi_at_interface]
    type     = PointValue
    variable = psi
    point    = '0.5 0 0'
  []

  # Inside the core at x=0.25: between center and edge.
  # Helps distinguish TE0 (monotone from center to edge) from
  # TE2 (has a local extremum inside the core).
  [psi_at_core_mid]
    type     = PointValue
    variable = psi
    point    = '0.25 0 0'
  []

  # In the cladding at x=1.0 (0.5 m from interface):
  # For guided modes: |psi(1.0)| << |psi(0.5)|
  # Decay should be approximately exponential.
  [psi_at_cladding_near]
    type     = PointValue
    variable = psi
    point    = '1 0 0'
  []

  # Far cladding at x=3.0: should be ~ 0 for all guided modes.
  [psi_at_cladding_far]
    type     = PointValue
    variable = psi
    point    = '3 0 0'
  []

  # Mode field extrema: maximum and minimum values of psi in the domain.
  # The ratio psi_min/psi_max identifies the mode symmetry:
  #   |psi_min/psi_max| ~ 1 (same magnitude, opposite sign) -> odd mode (TE1)
  #   psi_min ~ 0                                            -> even mode (TE0, TE2)
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

  # Average n^2 in the domain: sanity check.
  # Expected: avg_n2 = (core fraction) * n1^2 + (cladding fraction) * n2^2
  #         = (1.0/10.0) * 2.25 + (9.0/10.0) * 1.0 = 0.225 + 0.9 = 1.125
  # (Core width 1.0 m out of total domain width 10.0 m)
  [avg_n2]
    type     = ElementAverageValue
    variable = n_sq_field
  []
[]

[Executioner]
  type = Eigenvalue

  # KRYLOVSCHUR: SLEPc's Krylov-Schur algorithm for the generalised
  # eigenvalue problem A psi = lambda B psi. This algorithm is superior
  # to power iteration (PJFNK) for finding multiple eigenvalues and is
  # robust for the symmetric, well-conditioned 1D FEM stiffness system.
  solve_type = KRYLOVSCHUR

  # Request 6 eigen pairs to capture all guided modes plus a few
  # near-continuum modes for diagnostic purposes.
  # Guided modes expected: TE0, TE1, TE2 (3 modes for this waveguide).
  # Extra pairs (4-6) capture modes just below the cut-off.
  n_eigen_pairs     = 6
  which_eigen_pairs = TARGET_MAGNITUDE

  # Shift-invert spectral transformation:
  #   eps_target = 39.5 — shift sigma placed just ABOVE n2^2 k0^2 = 39.478
  #
  # Why target just above the cladding cutoff?
  # Shift-invert maps lambda -> 1/(lambda - sigma). With sigma = 39.5:
  #   - Radiation modes (lambda < 39.478): negative lambda-sigma, converge slowly
  #   - TE2 (lambda ~ 43, closest to sigma): |lambda-sigma|=3.5, largest 1/(lambda-sigma) -> converges first
  #   - TE1 (lambda ~ 60): |lambda-sigma|=20.5 -> converges second
  #   - TE0 (lambda ~ 80): |lambda-sigma|=40.5 -> converges third
  #
  # This ordering (highest mode first, fundamental last) is fine — SLEPc
  # finds all requested pairs simultaneously.
  #
  # st_type = sinvert: activates shift-invert spectral transformation.
  # LU inner solver (MUMPS): robust direct solver for the 400-DOF 1D system.
  # The (A - sigma B) matrix is 400x400 — trivially small for direct methods.
  petsc_options_iname = '-eps_target -st_type -st_ksp_type -st_pc_type -st_pc_factor_mat_solver_type'
  petsc_options_value = '39.5        sinvert  preonly       lu           mumps'
[]

[Outputs]
  # exodus: 1D mode field output — psi, n_sq_field, dpsi_dx along x.
  # Visualise in ParaView:
  #   - Plot psi vs. x to see the transverse mode profile:
  #       * TE0: smooth Gaussian-like peak centred at x=0 with exponential tails
  #       * TE1: antisymmetric profile with a zero crossing at x=0
  #       * TE2: symmetric with one interior zero crossing inside the core
  #   - Overlay n_sq_field to see the core (n^2=2.25) and cladding (n^2=1) regions.
  #   - Use "Plot Over Line" to extract psi(x) and compare with analytic modes.
  # execute_on = FINAL: write after eigensolver convergence (not during iterations).
  exodus     = true
  csv        = true
  execute_on = FINAL
[]
