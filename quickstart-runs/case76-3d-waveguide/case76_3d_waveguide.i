# ============================================================
# Case 76: 3D Rectangular Waveguide — TE10 Mode Propagation
# MIT 6.635 Advanced Electromagnetism, Lectures 3–4
#
# Reference: Pozar, Microwave Engineering, 4th ed. (2012), Sec. 3.3
#            Griffiths, Introduction to Electrodynamics, 4th ed., Sec. 9.5
#
# PHYSICS — THE RECTANGULAR WAVEGUIDE TE10 MODE
# ----------------------------------------------
# A rectangular metallic waveguide with perfectly conducting (PEC)
# walls occupies the volume:
#   Ω = [0, a] × [0, b] × [0, L]
#   a = 3 m  (x-direction, wide dimension)
#   b = 2 m  (y-direction, narrow dimension)
#   L = 10 m (z-direction, propagation axis)
#
# PEC walls on 4 longitudinal faces:  x=0, x=a, y=0, y=b
# Port injection at  z = 0  (boundary 'back')
# Absorbing (matched) BC at z = L  (boundary 'front')
#
# The TE10 mode has only a y-component of electric field:
#   E_y(x, z) = sin(pi*x/a) * exp(-j*beta_10*z)
#   E_x = 0,  E_z = 0
#
# The TE10 cutoff wavenumber and propagation constant:
#   k_c10 = pi/a = pi/3 ≈ 1.0472 m^-1
#   beta_10 = sqrt(k0^2 - (pi/a)^2)
#
# Operating conditions:
#   k0 = 1.5 rad/m
#   beta_10 = sqrt(1.5^2 - (pi/3)^2)
#           = sqrt(2.25 - 1.0966)
#           = sqrt(1.1534)
#           ≈ 1.0740 rad/m
#
# Single-mode condition:  only TE10 propagates.
# The next mode (TE20) has cutoff k_c20 = 2*pi/a = 2.094 m^-1 > k0 = 1.5.
# All higher-order modes are evanescent. TE10 is the ONLY propagating mode.
#
# GOVERNING EQUATION — VECTOR HELMHOLTZ
# ---------------------------------------
# In the frequency domain, with exp(+j*omega*t) convention removed,
# the electric field satisfies the vector Helmholtz equation:
#
#   curl(curl(E)) - k0^2 * E = 0   in Ω
#
# with PEC boundary condition (tangential E = 0) on the four walls.
#
# Separating into real and imaginary parts (E = E_real + j*E_imag):
#   curl(curl(E_real)) - k0^2 * E_real = 0
#   curl(curl(E_imag)) - k0^2 * E_imag = 0
#
# The two components are decoupled in the bulk equation (lossless vacuum,
# no imaginary part of k0^2). They couple ONLY through the port Robin BC.
#
# MOOSE KERNEL RESIDUALS
# -----------------------
# CurlCurlField:      adds +∫ (curl E)·(curl v) dV
#                     (weak form of curl-curl, integration by parts)
#
# VectorFunctionReaction (sign=negative, function=k0^2):
#                     adds -k0^2 * ∫ E·v dV
#                     (the -k0^2*E term in the Helmholtz equation)
#
# Total residual:  ∫ (curl E)·(curl v) dV - k0^2 * ∫ E·v dV = BC terms
# which enforces:  curl(curl(E)) - k0^2 * E = 0  ✓
#
# NEDELEC_ONE ELEMENTS
# ---------------------
# Nedelec edge elements of the first kind are REQUIRED for the curl-curl
# operator. They live in H(curl) — tangential components are continuous
# across element faces, while normal components can be discontinuous.
# This is exactly the right function space for the electric field in
# Maxwell's equations (tangential E continuous at interfaces).
#
# FIRST-order NEDELEC_ONE on HEX20 (second-order hexahedral):
# Each edge carries one DOF per HEX20 element.
# The mesh 10 × 8 × 20 = 1600 elements gives 20 elements per free-space
# wavelength lambda_0 = 2*pi/1.5 ≈ 4.19 m → element size ≈ 0.5 m,
# which is λ/8 — sufficient for accurate wave propagation.
#
# PORT BOUNDARY CONDITION (z = 0, boundary 'back')
# --------------------------------------------------
# VectorEMRobinBC in port mode injects the TE10 mode profile.
#
# The first-order Robin port condition (Jin, FEM in EM, Eq. 13.46):
#   n × curl(E) + j*beta * n × (n × E) = source_term   at z = 0
#
# where n = -z_hat (outward normal points in -z direction at back face),
# and the source term encodes the incoming TE10 field:
#   source = n × curl(E_inc) + j*beta * n × (n × E_inc)
#          = [2 * j*beta * n × (n × E_inc)] + n × curl(E_inc)
#
# TE10 incoming field at z = 0 (unit amplitude):
#   E_inc = (0, sin(pi*x/a), 0)     <- purely real
#
# Curl of E_inc = sin(pi*x/a) * exp(-j*beta*z) evaluated at z = 0:
#   curl(E_inc)_x = dEz/dy - dEy/dz = -(-j*beta)*sin(pi*x/a) = j*beta*sin(pi*x/a)
#   curl(E_inc)_y = dEx/dz - dEz/dx = 0
#   curl(E_inc)_z = dEy/dx - dEx/dy = (pi/a)*cos(pi*x/a)
#
# Real part of curl:   (0,          0, (pi/a)*cos(pi*x/a))
# Imaginary part of curl: (beta*sin(pi*x/a), 0, 0)
#
# These are provided to VectorEMRobinBC via real_incoming and imag_incoming
# ParsedVectorFunction objects (which must supply the curl via curl_x/y/z).
#
# ABSORBING BOUNDARY CONDITION (z = L, boundary 'front')
# --------------------------------------------------------
# VectorEMRobinBC in absorbing mode implements the Silver-Muller ABC:
#   n × curl(E) + j*beta * n × (n × E) = 0   at z = L
#
# This absorbs the outgoing TE10 wave without reflection. The absorbing
# BC uses the same beta (TE10 propagation constant) as the port injection.
#
# PEC WALLS (x=0, x=a, y=0, y=b)
# --------------------------------
# VectorCurlPenaltyDirichletBC enforces tangential E = 0 weakly through
# a penalty method:
#   penalty * ∫ (n × E)·(n × v) dS = 0   on PEC walls
#
# The penalty value 1e8 is chosen large enough to enforce the BC strongly
# while keeping the condition number manageable. For NEDELEC_ONE elements,
# this is the standard approach (strong Dirichlet is not available for
# vector Nedelec DOFs in the same way as for Lagrange DOFs).
#
# COMPUTED QUANTITIES
# --------------------
# AuxVariables: E_y_real, E_y_imag — y-component of the complex E-field.
# |E_y|^2 = E_y_real^2 + E_y_imag^2 — field intensity at each point.
# These are visualised in ParaView as cross-section slices at z = L/4, L/2, 3L/4.
#
# Postprocessors: max |E_y| and point values at z-midplane center
# to verify the sinusoidal sin(pi*x/a) profile and exp(-j*beta*z) propagation.
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables
# -----------------------------------------------------------
# Waveguide geometry
a      = 3          # width in x-direction [m]
b      = 2          # height in y-direction [m]
L      = 10         # length in z-direction [m]

# Wave numbers
k0     = 1.5        # free-space wavenumber [rad/m]

# TE10 mode:
#   k_c10 = pi/a = pi/3 = 1.04720 rad/m
#   beta_10 = sqrt(k0^2 - (pi/a)^2)
#           = sqrt(2.25 - (pi/3)^2)
#           = sqrt(2.25 - 1.09662)
#           = sqrt(1.15338)
#           = 1.07395 rad/m
beta10 = 1.07395    # TE10 propagation constant [rad/m]

# pi/a for the transverse profile function: pi/3 = 1.04720
pi_over_a = 1.04720   # pi/a [rad/m]

[Mesh]
  # 3D structured hexahedral mesh on [0,a] x [0,b] x [0,L].
  # Resolution: 10 x 8 x 20 = 1600 HEX elements.
  # With HEX20 (second-order hexahedral), each element has 20 nodes.
  # Element size: dx = 0.3 m, dy = 0.25 m, dz = 0.5 m.
  # At k0 = 1.5 rad/m, lambda_0 = 2*pi/1.5 = 4.19 m, giving about
  # 8 elements per free-space wavelength in z — adequate for accuracy.
  #
  # GeneratedMeshGenerator boundary names for a 3D box:
  #   left   -> x = 0      right  -> x = a
  #   bottom -> y = 0      top    -> y = b
  #   back   -> z = 0      front  -> z = L
  #
  # PEC walls: left, right, bottom, top  (4 longitudinal faces)
  # Port:      back   (injection at z = 0)
  # Absorber:  front  (matched termination at z = L)
  [gmg]
    type  = GeneratedMeshGenerator
    dim   = 3
    nx    = 10
    ny    = 8
    nz    = 20
    xmin  = 0
    xmax  = ${a}
    ymin  = 0
    ymax  = ${b}
    zmin  = 0
    zmax  = ${L}
    elem_type = HEX20
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor.
  # E = E_real + j * E_imag with exp(+j*omega*t) suppressed.
  # NEDELEC_ONE (first-order Nedelec edge elements) is required for
  # the curl-curl operator to be well-posed. These elements carry
  # tangential DOFs along each mesh edge, ensuring continuity of the
  # tangential field component across element boundaries — the correct
  # physical condition for the electric field in free space.
  [E_real]
    order  = FIRST
    family = NEDELEC_ONE
  []

  # E_imag — imaginary part of the complex electric field phasor.
  # Same element type as E_real. The two components are decoupled in
  # the bulk Helmholtz equation for lossless (real k0^2) media,
  # but are coupled at the port boundary through the Robin condition.
  [E_imag]
    order  = FIRST
    family = NEDELEC_ONE
  []
[]

[AuxVariables]
  # Y-component of E_real, extracted from the Nedelec vector variable.
  # NEDELEC_ONE variables are edge-based (elemental), so AuxVariables
  # that couple to them must also be elemental (CONSTANT MONOMIAL).
  # Used to observe the sin(pi*x/a) transverse profile and propagation.
  [E_y_real]
    order  = CONSTANT
    family = MONOMIAL
  []

  # Y-component of E_imag. Together with E_y_real, we get the full
  # complex phasor amplitude: E_y = E_y_real + j * E_y_imag.
  # The phase shift along z: E_y_real ~ cos(beta*z), E_y_imag ~ -sin(beta*z)
  # for a forward-propagating wave injected with unit real amplitude.
  [E_y_imag]
    order  = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # k0sq_fn: wavenumber squared k0^2 = 1.5^2 = 2.25 [rad^2/m^2]
  # Used by VectorFunctionReaction to implement the -k0^2 * E term.
  # VectorFunctionReaction with sign=negative adds: -k0sq * E * test
  # which is the correct discretisation of -k0^2 * E in the
  # curl-curl Helmholtz equation.
  # Computed as ${k0}^2 so that k0 is the single source of truth.
  # ------------------------------------------------------------------
  [k0sq_fn]
    type       = ParsedFunction
    expression = '${k0} * ${k0}'
  []

  # ------------------------------------------------------------------
  # beta10_fn: TE10 propagation constant as a ParsedFunction.
  # VectorEMRobinBC takes beta as a FunctionName, so even a constant
  # must be wrapped in a ParsedFunction.
  # beta_10 = sqrt(k0^2 - (pi/a)^2) = sqrt(2.25 - 1.0966) = 1.07395
  # ------------------------------------------------------------------
  [beta10_fn]
    type       = ParsedFunction
    expression = '${beta10}'
  []

  # ------------------------------------------------------------------
  # TE10 incoming field (real part) at the port z = 0.
  #
  # The TE10 mode profile (unit amplitude):
  #   E_inc_y(x) = sin(pi*x/a) = sin(1.04720 * x)
  #   E_inc_x = 0,  E_inc_z = 0
  #
  # For a wave travelling in the +z direction with exp(-j*beta*z):
  #   E_inc = (0, sin(pi*x/a) * exp(-j*beta*z), 0)
  #
  # At z = 0, the real part of the field is:
  #   E_inc_real = (0, sin(pi*x/a), 0)
  #
  # Curl of E_inc * exp(-j*beta*z) evaluated at z = 0:
  #   curl_x = dEz/dy - dEy/dz = 0 - (-j*beta)*sin(pi*x/a) = +j*beta*sin(pi*x/a)
  #   curl_y = dEx/dz - dEz/dx = 0
  #   curl_z = dEy/dx - dEx/dy = (pi/a)*cos(pi*x/a)
  #
  # Real part of curl at z = 0:
  #   curl_x_real = 0  (comes from j*beta — purely imaginary)
  #   curl_y_real = 0
  #   curl_z_real = (pi/a)*cos(pi*x/a)  [from the spatial derivative in z-direction]
  # ------------------------------------------------------------------
  [inc_real]
    type         = ParsedVectorFunction
    expression_x = '0'
    expression_y = 'sin(${pi_over_a} * x)'
    expression_z = '0'
    curl_x       = '0'
    curl_y       = '0'
    curl_z       = '${pi_over_a} * cos(${pi_over_a} * x)'
  []

  # ------------------------------------------------------------------
  # TE10 incoming field (imaginary part) at the port z = 0.
  #
  # At z = 0, the imaginary part of the field is zero:
  #   E_inc_imag = (0, 0, 0)
  #
  # The imaginary part of the curl at z = 0:
  #   curl_x_imag = +beta * sin(pi*x/a)  [from j*beta term]
  #   curl_y_imag = 0
  #   curl_z_imag = 0
  # ------------------------------------------------------------------
  [inc_imag]
    type         = ParsedVectorFunction
    expression_x = '0'
    expression_y = '0'
    expression_z = '0'
    curl_x       = '${beta10} * sin(${pi_over_a} * x)'
    curl_y       = '0'
    curl_z       = '0'
  []
[]

[Kernels]
  # ==================================================================
  # REAL COMPONENT EQUATION:
  #   curl(curl(E_real)) - k0^2 * E_real = 0
  #
  # Weak form after integration by parts:
  #   ∫ (curl E_real)·(curl v) dV - k0^2 * ∫ E_real·v dV = BC terms
  # ==================================================================

  # CurlCurlField: contributes +∫ (curl E_real)·(curl v) dV
  # This is the weak form of curl(curl(E)) with the sign convention
  # that the natural BC (boundary term from integration by parts)
  # is: -∫ (n × curl E)·v dS, absorbed into the Robin BC.
  [curl_curl_real]
    type     = CurlCurlField
    variable = E_real
  []

  # VectorFunctionReaction (sign=negative): contributes -k0^2 * ∫ E_real·v dV
  # With sign=negative and function=k0sq_fn=2.25, this adds:
  #   residual += -2.25 * ∫ E_real·v dV
  # which is the -k0^2 * E_real term needed for the Helmholtz equation.
  [wave_reaction_real]
    type     = VectorFunctionReaction
    variable = E_real
    function = k0sq_fn
    sign     = negative
  []

  # ==================================================================
  # IMAGINARY COMPONENT EQUATION:
  #   curl(curl(E_imag)) - k0^2 * E_imag = 0
  #
  # Identical structure to the real equation — lossless medium means
  # bulk equations for real and imaginary parts are decoupled.
  # Coupling enters only through the port Robin boundary condition.
  # ==================================================================

  [curl_curl_imag]
    type     = CurlCurlField
    variable = E_imag
  []

  [wave_reaction_imag]
    type     = VectorFunctionReaction
    variable = E_imag
    function = k0sq_fn
    sign     = negative
  []
[]

[AuxKernels]
  # Extract the y-component of the NEDELEC_ONE vector variable E_real
  # into a scalar Lagrange variable for easy postprocessing/visualisation.
  # VectorVariableComponentAux picks component 1 (y) from a vector variable.
  [E_y_real_aux]
    type      = VectorVariableComponentAux
    variable  = E_y_real
    vector_variable = E_real
    component = y
  []

  [E_y_imag_aux]
    type      = VectorVariableComponentAux
    variable  = E_y_imag
    vector_variable = E_imag
    component = y
  []
[]

[BCs]
  # ==================================================================
  # PEC WALLS: VectorCurlPenaltyDirichletBC on four longitudinal faces
  # Enforces: n × E = 0  (tangential electric field = 0)
  # Applied to: left (x=0), right (x=a), bottom (y=0), top (y=b)
  #
  # The penalty method adds:
  #   penalty * ∫ (n × E)·(n × v) dS = penalty * ∫ (n × E_prescribed)·(n × v) dS
  # with E_prescribed = 0 (zero tangential field on PEC).
  #
  # Penalty = 1e8 is large enough to enforce the BC strongly without
  # excessive ill-conditioning. The exact value is not critical as long
  # as penalty >> max(k0^2 * h^2) ~ 1.5^2 * 0.25 ~ 0.56.
  # ==================================================================

  # Real component — PEC walls (tangential E_real = 0)
  [pec_walls_real]
    type     = VectorCurlPenaltyDirichletBC
    variable = E_real
    boundary = 'left right bottom top'
    penalty  = 1e8
  []

  # Imaginary component — PEC walls (tangential E_imag = 0)
  [pec_walls_imag]
    type     = VectorCurlPenaltyDirichletBC
    variable = E_imag
    boundary = 'left right bottom top'
    penalty  = 1e8
  []

  # ==================================================================
  # PORT INJECTION (z = 0, boundary 'back')
  # VectorEMRobinBC in port mode simultaneously:
  #   (1) Injects the TE10 mode profile sin(pi*x/a)
  #   (2) Absorbs any reflected wave that returns to the port
  #
  # The port condition (Jin, FEM in EM, Ch. 13):
  #   n × curl(E) + j*beta * n × (n × E) = S_port
  # where S_port = n × curl(E_inc) + j*beta * n × (n × E_inc)
  #
  # This is implemented by VectorEMRobinBC with:
  #   mode = port
  #   real_incoming = inc_real  (real part of incoming vector field)
  #   imag_incoming = inc_imag  (imaginary part of incoming vector field)
  #
  # Note: both real_incoming and imag_incoming MUST be set when mode=port.
  # VectorEMRobinBC calls vectorValue() for the field and curl() for the
  # curl contribution — hence the curl_x/y/z parameters in the functions
  # above are critical.
  #
  # The 'beta' parameter is the propagation constant of the mode being
  # injected — beta_10 = 1.07395 for the TE10 mode at k0 = 1.5.
  # ==================================================================

  # Real-component port BC at z = 0
  [port_real]
    type          = VectorEMRobinBC
    variable      = E_real
    coupled_field = E_imag
    boundary      = back
    component     = real
    mode          = port
    beta          = beta10_fn
    real_incoming = inc_real
    imag_incoming = inc_imag
  []

  # Imaginary-component port BC at z = 0
  [port_imag]
    type          = VectorEMRobinBC
    variable      = E_imag
    coupled_field = E_real
    boundary      = back
    component     = imaginary
    mode          = port
    beta          = beta10_fn
    real_incoming = inc_real
    imag_incoming = inc_imag
  []

  # ==================================================================
  # ABSORBING BC (z = L, boundary 'front')
  # VectorEMRobinBC in absorbing mode implements the Silver-Muller ABC:
  #   n × curl(E) + j*beta * n × (n × E) = 0   at z = L
  #
  # This is a first-order absorbing condition that is exact for a plane
  # wave propagating at angle theta to the normal. For the TE10 mode
  # (normally incident on the flat termination), this condition is
  # effectively perfect — all outgoing TE10 power is absorbed.
  #
  # The beta here is again the TE10 propagation constant, ensuring that
  # the characteristic impedance matched to the injected mode is used.
  # ==================================================================

  # Real-component absorbing BC at z = L
  [absorb_real]
    type          = VectorEMRobinBC
    variable      = E_real
    coupled_field = E_imag
    boundary      = front
    component     = real
    mode          = absorbing
    beta          = beta10_fn
  []

  # Imaginary-component absorbing BC at z = L
  [absorb_imag]
    type          = VectorEMRobinBC
    variable      = E_imag
    coupled_field = E_real
    boundary      = front
    component     = imaginary
    mode          = absorbing
    beta          = beta10_fn
  []
[]

[Postprocessors]
  # ------------------------------------------------------------------
  # Maximum |E_y| over the entire domain.
  # For a perfectly transmitted TE10 wave with unit amplitude,
  # max|E_y| should be close to 1 (some overshoot is possible near
  # the port due to the mode matching; values of 1.0 ± 0.05 indicate
  # a healthy simulation).
  # ------------------------------------------------------------------
  [max_E_y_real]
    type     = ElementExtremeValue
    variable = E_y_real
    value_type = max
  []

  [max_E_y_imag]
    type     = ElementExtremeValue
    variable = E_y_imag
    value_type = max
  []

  # ------------------------------------------------------------------
  # Point values at the waveguide centre cross-section (z = L/2 = 5).
  # At x = a/2 = 1.5, y = b/2 = 1.0, z = 5.0:
  #   E_y_real(L/2) = sin(pi*1.5/3) * cos(beta10*5)
  #                 = sin(pi/2) * cos(1.07395*5)
  #                 = 1.0 * cos(5.3698)
  #                 = 1.0 * cos(5.37)
  #                 ≈ 0.776
  #   E_y_imag(L/2) = -sin(pi*1.5/3) * sin(beta10*5)
  #                 = -1.0 * sin(5.37)
  #                 ≈ -0.631
  # These predicted values verify correct phase propagation.
  # ------------------------------------------------------------------
  [E_y_real_center]
    type     = PointValue
    variable = E_y_real
    point    = '1.5 1.0 5.0'
  []

  [E_y_imag_center]
    type     = PointValue
    variable = E_y_imag
    point    = '1.5 1.0 5.0'
  []

  # ------------------------------------------------------------------
  # Point value near the port (z = 0.25 m, quarter-element spacing).
  # At x = a/2 = 1.5 (field maximum of sin profile):
  #   E_y_real(0.25) = cos(beta10*0.25) ≈ cos(0.268) ≈ 0.964
  #   E_y_imag(0.25) = -sin(beta10*0.25) ≈ -sin(0.268) ≈ -0.265
  # ------------------------------------------------------------------
  [E_y_real_near_port]
    type     = PointValue
    variable = E_y_real
    point    = '1.5 1.0 0.25'
  []

  [E_y_imag_near_port]
    type     = PointValue
    variable = E_y_imag
    point    = '1.5 1.0 0.25'
  []

  # ------------------------------------------------------------------
  # Point value near the absorbing end (z = L - 0.25 = 9.75 m).
  # Checks that the wave is propagating cleanly to the termination.
  # ------------------------------------------------------------------
  [E_y_real_near_end]
    type     = PointValue
    variable = E_y_real
    point    = '1.5 1.0 9.75'
  []

  [E_y_imag_near_end]
    type     = PointValue
    variable = E_y_imag
    point    = '1.5 1.0 9.75'
  []

  # ------------------------------------------------------------------
  # Off-centre points to verify the sin(pi*x/a) transverse profile.
  # At x = a/4 = 0.75, sin(pi*0.75/3) = sin(pi/4) = 0.707.
  # At x = 3*a/4 = 2.25, sin(pi*2.25/3) = sin(3*pi/4) = 0.707.
  # Both should give the same amplitude by symmetry.
  # ------------------------------------------------------------------
  [E_y_real_quarter_x]
    type     = PointValue
    variable = E_y_real
    point    = '0.75 1.0 5.0'
  []

  [E_y_real_threequarter_x]
    type     = PointValue
    variable = E_y_real
    point    = '2.25 1.0 5.0'
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner is REQUIRED for the coupled
  # (E_real, E_imag) system. The VectorEMRobinBC Robin condition introduces
  # off-diagonal Jacobian blocks between E_real and E_imag. Without
  # full=true the SMP preconditioner would omit these cross-coupling blocks,
  # leading to poor convergence of the linear solver near the port.
  #
  # Direct LU factorisation (pc_type = lu) is used for robustness.
  # The curl-curl system is indefinite (negative eigenvalues exist below
  # the lowest eigenfrequency), so iterative methods like CG would fail.
  # LU handles the indefinite system exactly at the cost of more memory.
  # For 1600 HEX20 elements with 2 NEDELEC_ONE variables, the system
  # has O(10,000) DOFs — well within the range where direct LU is fast.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU: exact factorisation, robust for indefinite wave problems.
  # MUMPS is preferred for 3D problems due to its fill-reducing ordering
  # and parallel capability, but sequential superlu_dist also works.
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  # Tight tolerances for the wave propagation problem.
  # The linear system from the Helmholtz equation is well-conditioned
  # with direct LU, so these tolerances should be easily achieved.
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # Exodus output: full 3D vector field (E_real, E_imag) and scalar
  # auxiliary variables (E_y_real, E_y_imag) for ParaView visualisation.
  # Use ParaView's "Slice" filter to view cross-sections at multiple z-planes.
  # Use "Warp by Scalar" or "Glyph" on E_y_real to see the standing wave.
  # Use "Calculator" to compute |E_y|^2 = E_y_real^2 + E_y_imag^2.
  exodus = true

  # CSV output: all postprocessor values (point values and extrema)
  # for comparison with the analytical TE10 solution.
  csv    = true

  # Write output only after the steady-state solve converges.
  # For transient problems this would write every time step, but here
  # we only want the final converged solution.
  execute_on = FINAL
[]
