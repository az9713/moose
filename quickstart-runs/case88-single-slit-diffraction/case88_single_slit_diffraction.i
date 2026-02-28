# ============================================================
# Case 88: Single-Slit Diffraction — Fraunhofer Limit and Resolution
# Hecht, "Optics", 5th Ed. (2016), Ch. 10
# Born & Wolf, "Principles of Optics", 7th Ed. (1999), Ch. 8
# Jin, "The Finite Element Method in Electromagnetics", 3rd Ed. (2014), Ch. 9
#
# MIT 6.635 connection: Lecture 3-4 (diffraction theory, Huygens principle)
#
# ============================================================
# HISTORICAL AND PHYSICAL BACKGROUND
# ============================================================
#
# Single-slit diffraction is the canonical demonstration of wave nature
# of light, first analysed quantitatively by Augustin-Jean Fresnel (1818)
# and later simplified in the far-field (Fraunhofer) limit. It underlies
# the Rayleigh criterion for optical resolution:
#
#   sin(theta_min) = 1.22 lambda / D   (circular aperture)
#   sin(theta_null) = lambda / a        (single slit, first null)
#
# where a is the slit width. The single-slit case is the building block
# for understanding gratings, arrays, and the diffraction limit of imaging
# systems.
#
# ============================================================
# GEOMETRY AND SETUP
# ============================================================
#
#   x = 0                                     x = 15
#   Screen/Slit plane                         Observation / ABC port
#
#        |y = +1.5                                    |
#        |  slit: E = 1 (plane wave amplitude)         |
#        |y = 0                                        |
#        |  slit: E = 1                                |
#        |y = -1.5                                     |
#        |  screen: E = 0 (PEC)                        |
#        |y = -8                                       |
#   _____|_______________________________________________|
#   left boundary                              right boundary
#   (FunctionDirichletBC)                     (EMRobinBC absorbing)
#
# Domain:   [0, 15] x [-8, 8]   (x is the propagation direction)
# Slit:     |y| < a/2 = 1.5 m   (width a = 3lambda_0 with lambda_0 = 1)
# Screen:   |y| >= 1.5 m on x = 0 (Dirichlet E = 0, PEC condition)
#
# ============================================================
# PHYSICS: FREQUENCY-DOMAIN HELMHOLTZ EQUATION
# ============================================================
#
# For TE polarisation (E_z out of the 2D plane, propagation in x), the
# scalar Helmholtz equation in free space is:
#
#   nabla^2 E_z + k0^2 E_z = 0
#
# where k0 = 2pi/lambda_0 = 2pi (with lambda_0 = 1 m).
#
# The complex phasor field E_z = E_real + j E_imag splits into two
# real-valued equations:
#
#   nabla^2 E_real + k0^2 E_real = 0
#   nabla^2 E_imag + k0^2 E_imag = 0
#
# These are coupled only through the Robin absorbing BCs on the open
# boundaries (right, top, bottom). In the bulk (lossless, purely real
# wavenumber) the two equations are identical and decoupled.
#
# ============================================================
# BOUNDARY CONDITIONS
# ============================================================
#
# LEFT BOUNDARY (x = 0) — Screen and Slit Plane:
# ------------------------------------------------
# This is the key boundary for the diffraction problem. We use
# FunctionDirichletBC with a piecewise function:
#
#   E_real(0, y) = if(|y| < a/2, 1.0, 0.0)
#   E_imag(0, y) = 0.0   (plane wave arrives with real amplitude at x=0)
#
# Physical interpretation:
#   - At the slit (|y| < 1.5): the incoming plane wave E_inc = exp(jk0 x)
#     at x=0 has amplitude E_real=1, E_imag=0. We prescribe this as
#     a Dirichlet condition, equivalent to injecting the Huygens sources
#     at the aperture (Kirchhoff diffraction theory).
#   - At the screen (|y| >= 1.5): PEC wall, E_total = 0, so E_real = E_imag = 0.
#
# This direct "aperture field" formulation is mathematically equivalent
# to the Kirchhoff boundary condition used in scalar diffraction theory:
#   - Inside the aperture: field is the same as if no screen were present
#   - Outside (on the screen): field is exactly zero
#
# Note: This is the Kirchhoff approximation — it ignores edge diffraction
# (the Sommerfeld correction) and is accurate when a >> lambda_0.
# Here a = 3 lambda_0, which is in the regime where the Fraunhofer
# pattern is well-developed with several visible fringes.
#
# RIGHT BOUNDARY (x = 15, observation plane):
# --------------------------------------------
# EMRobinBC absorbing — absorbs the outgoing diffracted wave without
# significant spurious reflections. The domain is 15 lambda_0 long,
# which is well into the Fraunhofer (far-field) regime:
#
#   Far-field condition: z >> a^2 / lambda_0 = 9 / 1 = 9 m
#   Observation distance: 15 m > 9 m  -> Fraunhofer regime
#
# TOP AND BOTTOM BOUNDARIES (y = +/-8):
# --------------------------------------
# EMRobinBC absorbing only. The domain half-height 8 m > a/2 = 1.5 m
# by a factor of 5, which is sufficient for the first several diffraction
# maxima to develop without boundary truncation effects.
#
# ============================================================
# ANALYTICAL DIFFRACTION PATTERN (FRAUNHOFER)
# ============================================================
#
# In the Fraunhofer limit (far field), the intensity pattern from a
# single slit of width a illuminated by a plane wave is:
#
#   I(theta) = I_0 * sinc^2(beta)
#
# where:
#   beta = (pi a / lambda_0) * sin(theta)
#   sinc(beta) = sin(beta) / beta
#   theta = angle from the optical axis (x-direction)
#
# Nulls (zeros of I) occur at:
#   sin(theta_m) = m * lambda_0 / a  for m = +/-1, +/-2, ...
#
# With a = 3 lambda_0 = 3 m, the first nulls are at:
#   sin(theta_1) = 1/3   ->   theta_1 = arcsin(1/3) = 19.47 degrees
#   At observation distance x = 10:
#     y_null = 10 * tan(theta_1) = 10 * sin(theta_1)/cos(theta_1)
#            = 10 * (1/3) / sqrt(1 - 1/9)
#            = 10 * (1/3) / sqrt(8/9)
#            = 10 / (3 * sqrt(8/9))
#            ~ 3.536 m  (y-position of first null at x=10)
#
# Secondary maxima occur at approximately:
#   sin(theta) ~ (2m+1) * lambda_0 / (2a)   for m = 1, 2, ...
#   m=1: sin(theta) = 3/(6) = 0.5 -> theta = 30 deg, y ~ 5.77 m at x=10
#   (secondary maximum intensity ~ 4/(9 pi^2) * I_0 ~ 0.045 I_0)
#
# ============================================================
# NUMERICAL VALUES
# ============================================================
#
# lambda_0 = 1.0 m
# k0 = 2 pi / lambda_0 = 6.283185307... rad/m
# k0^2 = 4 pi^2 = 39.4784176... m^-2
# Slit half-width: a/2 = 1.5 m
# Domain: [0, 15] x [-8, 8] = 15 m x 16 m
# Mesh: 150 x 160 quadrilateral elements
#   dx = 15/150 = 0.1 m = lambda_0/10 (10 elements per wavelength)
#   dy = 16/160 = 0.1 m = lambda_0/10
# Resolution is adequate for the smooth diffraction pattern.
#
# ============================================================
# KERNEL SIGN CONVENTIONS (verified against upstream cases)
# ============================================================
#
# Diffusion kernel (MOOSE standard):
#   Residual: +integral( nabla E . nabla v dV )
#   Strong form equivalent: -nabla^2 E
#
# ADMatReaction kernel:
#   Residual: -reaction_rate * integral( E * v dV )
#   Strong form equivalent: -reaction_rate * E
#
# Combined (Diffusion + ADMatReaction with rate = k0^2):
#   Strong form: -nabla^2 E - k0^2 E = 0
#   Equivalent:  nabla^2 E + k0^2 E = 0  (Helmholtz)  checkmark
#
# EMRobinBC (Jin, "FEM in Electromagnetics", Eq. 9.60):
#   Implements:  dn E + j k0 E = 0  (absorbing BC, no injection)
#   sign = negative: consistent with Diffusion kernel convention.
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...})
# -----------------------------------------------------------
# k0 = 2pi/lambda_0 = 2pi with lambda_0 = 1 m
# E0 = 0: no incident wave injected from the Robin port boundaries
#         (the excitation comes from the Dirichlet BC on the left face)
# theta = 0: normal incidence angle placeholder for EMRobinBC
k0    = 6.283185307179586    # free-space wavenumber [rad/m], lambda_0 = 1 m
E0    = 0                    # no incident wave from Robin ports
theta = 0                    # incidence angle [degrees] (placeholder)

[Mesh]
  # 2D rectangular domain [0, 15] x [-8, 8].
  # Propagation direction: x (left = screen/slit at x=0; right = far field at x=15)
  # Transverse direction:  y (slit at |y| < 1.5 m; screen at |y| >= 1.5 m)
  #
  # Element counts:
  #   nx = 150: 15 m / 150 = 0.1 m per element = lambda_0/10 (10 per wavelength)
  #   ny = 160: 16 m / 160 = 0.1 m per element = lambda_0/10
  #
  # Resolution requirement: >= 10 elements per wavelength for FIRST-order
  # Lagrange elements on a smooth Helmholtz problem. The diffraction pattern
  # is smooth away from the slit edges, so lambda_0/10 resolution is adequate.
  # The slit edges at y = +/-1.5 m are aligned with mesh nodes (1.5/0.1 = 15
  # elements from y=0 to y=+/-1.5), ensuring exact representation of the
  # step discontinuity in the Dirichlet BC.
  #
  # GeneratedMeshGenerator in 2D names the four boundaries:
  #   left   -> x = 0  (screen/slit plane)
  #   right  -> x = 15 (far-field observation boundary)
  #   bottom -> y = -8  (absorbing boundary)
  #   top    -> y = +8  (absorbing boundary)
  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 150
    ny   = 160
    xmin = 0
    xmax = 15
    ymin = -8
    ymax = 8
  []
[]

[Variables]
  # E_real — real part of the complex phasor field E_z = E_real + j E_imag.
  # Both components use first-order Lagrange elements, which enforce C^0
  # continuity across elements and can represent the smooth far-field
  # diffraction pattern. On the left boundary the Dirichlet BC prescribes
  # a step function (1 in the slit, 0 on the screen), which Lagrange elements
  # can represent exactly at the nodes.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the phasor.
  # The incoming plane wave at x=0 has E_inc = exp(jk0*0) = 1 + 0j, so
  # the imaginary part of the aperture field is zero. E_imag is driven only
  # by the coupling through the Robin absorbing BCs on the open boundaries.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # slit_fn: aperture field distribution at x = 0 (left boundary)
  #
  # Kirchhoff boundary condition for the single-slit diffraction:
  #   E(0, y) = 1  if |y| < a/2 = 1.5 m  (open aperture)
  #   E(0, y) = 0  if |y| >= 1.5 m        (PEC screen)
  #
  # This step function is the exact Kirchhoff boundary condition used in
  # scalar diffraction theory (see Goodman, "Introduction to Fourier
  # Optics", Ch. 3). It assumes:
  #   (1) The incident wave has uniform amplitude E_0 = 1 across the slit.
  #   (2) The screen blocks all field instantaneously (zero tangential E).
  #   (3) Edge diffraction effects at y = +/-1.5 are neglected (valid when
  #       a >> lambda_0, here a = 3 lambda_0 which is marginal).
  #
  # The 'abs(y) < 1.5' condition creates a symmetric aperture about y=0,
  # producing a symmetric (even) diffraction pattern.
  # ------------------------------------------------------------------
  [slit_fn]
    type       = ParsedFunction
    expression = 'if(abs(y) < 1.5, 1.0, 0.0)'
  []

  # cosTheta = cos(0 degrees) = 1.0 for normal incidence.
  # The EMRobinBC uses k_eff = k0 * cos(theta) as the effective wavenumber
  # along the boundary normal. At normal incidence (theta=0), k_eff = k0.
  # For the absorbing-only BCs (profile_func_real = 0), this just sets
  # the impedance of the ABC layer.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # Wrap k0^2 as an AD material property consumed by ADMatReaction.
  # ADGenericConstantMaterial provides a spatially uniform ADReal material
  # property. For free space (uniform k0^2 everywhere in the domain),
  # this is simpler and more efficient than ADGenericFunctionMaterial.
  #
  # Both E_real and E_imag share the same coefficient k0^2, consistent
  # with the uniform vacuum propagation medium.
  #
  # k0^2 = (2pi)^2 = 4 pi^2 = 39.47841760435743 m^-2
  [k0sq_mat]
    type        = ADGenericConstantMaterial
    prop_names  = 'k0sq'
    prop_values = '39.47841760435743'
  []
[]

[Kernels]
  # ==================================================================
  # E_real EQUATION:
  #   nabla^2 E_real + k0^2 E_real = 0
  #
  # Weak form (multiply by test function v, integrate by parts):
  #   integral( nabla E_real . nabla v dV ) - integral( k0^2 E_real v dV ) = 0
  # (boundary term absorbed by the Robin BCs on the open faces)
  # ==================================================================

  # Diffusion kernel: adds +integral( nabla E_real . nabla v dV ) to residual.
  # Represents the strong-form operator -nabla^2 E_real.
  [diff_real]
    type     = Diffusion
    variable = E_real
  []

  # ADMatReaction: adds -k0^2 * integral( E_real * v dV ) to residual.
  # Note the NEGATIVE sign convention: the kernel subtracts rate*E from the
  # residual, so with rate = +k0^2 the combined strong form is:
  #   -nabla^2 E_real - k0^2 E_real = 0   <=>   nabla^2 E_real + k0^2 E_real = 0  checkmark
  [react_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq
  []

  # ==================================================================
  # E_imag EQUATION:
  #   nabla^2 E_imag + k0^2 E_imag = 0
  #
  # Identical bulk equation as E_real (lossless, uniform medium).
  # The two components are coupled only through the EMRobinBC terms
  # on the right, top, and bottom boundaries (the j factor in the
  # Sommerfeld ABC mixes real and imaginary parts).
  # ==================================================================

  [diff_imag]
    type     = Diffusion
    variable = E_imag
  []

  [react_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq
  []
[]

[BCs]
  # ==================================================================
  # LEFT BOUNDARY (x = 0): Screen/Slit Plane — Kirchhoff BC
  # ==================================================================
  # FunctionDirichletBC prescribes E = slit_fn at every left-boundary node.
  # The function slit_fn returns:
  #   1.0 for |y| < 1.5 m  (slit aperture: incoming plane wave amplitude)
  #   0.0 for |y| >= 1.5 m (conducting screen: PEC condition)
  #
  # This implements the Kirchhoff diffraction approximation:
  #   - The aperture is illuminated by an incident plane wave with E=1 at x=0.
  #   - The perfectly conducting screen zeroes the tangential E field outside
  #     the aperture, which is the physical PEC boundary condition.
  #
  # E_imag = 0 on the entire left boundary: the incoming plane wave
  # exp(jk0 x) at x=0 evaluates to exp(0) = 1 + 0j, so the imaginary
  # part is zero at the aperture plane.
  [slit_real]
    type     = FunctionDirichletBC
    variable = E_real
    boundary = left
    function = slit_fn
  []

  [slit_imag]
    type     = DirichletBC
    variable = E_imag
    boundary = left
    value    = 0
  []

  # ==================================================================
  # RIGHT BOUNDARY (x = 15 m): Far-Field Observation / Absorbing Port
  # ==================================================================
  # EMRobinBC absorbing only (profile_func_real = 0):
  #   dn E + j k0 E = 0   (no incident wave injected, pure ABC)
  #
  # This is a first-order Sommerfeld absorbing boundary condition. It
  # is appropriate here because the diffracted wave at x=15 is nearly
  # a plane wave propagating in the +x direction (at least for the
  # central maximum and first few sidelobes within +/-30 degrees).
  # The ABC absorbs outgoing waves with low spurious reflection for
  # near-normally incident waves.
  #
  # Observation: At x=15, we sample the diffraction pattern as a
  # function of y. This is the "Fraunhofer observation screen."
  [abc_real_right]
    type              = EMRobinBC
    variable          = E_real
    boundary          = right
    component         = real
    coeff_real        = ${k0}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [abc_imag_right]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = right
    component         = imaginary
    coeff_real        = ${k0}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  # ==================================================================
  # TOP BOUNDARY (y = +8 m): Absorbing BC
  # ==================================================================
  # The diffracted wave hits the top boundary at large angles from
  # normal. The first-order ABC (EMRobinBC) has increasing reflection
  # error at glancing angles (theta -> 90 deg), but since the field
  # intensity decreases rapidly away from the slit (the sinc^2 envelope
  # falls off as y^-2 in the Fraunhofer limit), the boundary truncation
  # error is small for y >= 8 m = 5.33 * (a/2).
  [abc_real_top]
    type              = EMRobinBC
    variable          = E_real
    boundary          = top
    component         = real
    coeff_real        = ${k0}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [abc_imag_top]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = top
    component         = imaginary
    coeff_real        = ${k0}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  # ==================================================================
  # BOTTOM BOUNDARY (y = -8 m): Absorbing BC
  # ==================================================================
  # Symmetric to the top boundary by the symmetry of the aperture.
  # The field is even in y (slit symmetric about y=0), so the pattern
  # at y=-8 m mirrors that at y=+8 m.
  [abc_real_bottom]
    type              = EMRobinBC
    variable          = E_real
    boundary          = bottom
    component         = real
    coeff_real        = ${k0}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [abc_imag_bottom]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = bottom
    component         = imaginary
    coeff_real        = ${k0}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []
[]

[AuxVariables]
  # E_intensity = |E_z|^2 = E_real^2 + E_imag^2
  # This is the time-averaged field intensity proportional to the
  # cycle-averaged Poynting flux (power per unit transverse length).
  # Plotting E_intensity in ParaView reveals the classic sinc^2
  # diffraction envelope: a bright central maximum at y=0 flanked
  # by symmetrically placed secondary maxima and nulls.
  [E_intensity]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # ParsedAux evaluates |E_z|^2 = E_real^2 + E_imag^2 at every mesh node
  # using the converged primary field values. This is a purely algebraic
  # auxiliary computation — no PDE solve is required.
  [intensity_aux]
    type              = ParsedAux
    variable          = E_intensity
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # DIFFRACTION PATTERN SAMPLING AT x = 10 m (Fraunhofer observation)
  # ==================================================================
  # We sample E_real, E_imag, and E_intensity at 13 y-positions along
  # the vertical line x = 10 m, which is well into the Fraunhofer
  # regime (10 m > a^2/lambda_0 = 9 m).
  #
  # The observation points are chosen to resolve:
  #   - Central maximum at y = 0
  #   - First nulls at y ~ +/-3.536 m (sin(theta)=1/3, theta=19.47 deg)
  #   - Secondary maxima near y ~ +/-5.77 m (sin(theta)=0.5, theta=30 deg)
  #
  # Analytical E_intensity (unnormalised sinc^2 in far field):
  #   beta(y) = pi * a / lambda_0 * y / sqrt(x^2 + y^2)
  #           = 3pi * y / sqrt(100 + y^2)   at x=10
  #   I(y) = I_0 * (sin(beta)/beta)^2
  #
  # Central maximum (y=0): I = I_0 (maximum by L'Hopital, sinc(0) = 1)
  # First null (y~3.536):  I = 0
  # Secondary max (y~5.77): I ~ 0.045 I_0  (much weaker)
  # ==================================================================

  # -- On-axis: central maximum --
  [E_intensity_y0]
    type     = PointValue
    variable = E_intensity
    point    = '10 0 0'
  []
  [E_real_y0]
    type     = PointValue
    variable = E_real
    point    = '10 0 0'
  []
  [E_imag_y0]
    type     = PointValue
    variable = E_imag
    point    = '10 0 0'
  []

  # -- y = 1.0 m: inside the central lobe --
  [E_intensity_y1p0]
    type     = PointValue
    variable = E_intensity
    point    = '10 1 0'
  []

  # -- y = 2.0 m: approaching the first null --
  [E_intensity_y2p0]
    type     = PointValue
    variable = E_intensity
    point    = '10 2 0'
  []

  # -- y = 3.0 m: near but inside first null --
  # Expected: I still above zero but declining rapidly
  [E_intensity_y3p0]
    type     = PointValue
    variable = E_intensity
    point    = '10 3 0'
  []

  # -- y = 3.5 m: near the first null location --
  # Analytic first null: y_null = 10*tan(arcsin(1/3)) ~ 3.536 m
  # Expected: I ~ 0 (null of sinc^2 pattern)
  [E_intensity_y3p5]
    type     = PointValue
    variable = E_intensity
    point    = '10 3.5 0'
  []

  # -- y = 4.0 m: just past the first null (entering secondary lobe) --
  [E_intensity_y4p0]
    type     = PointValue
    variable = E_intensity
    point    = '10 4 0'
  []

  # -- y = 5.0 m: rising toward secondary maximum --
  [E_intensity_y5p0]
    type     = PointValue
    variable = E_intensity
    point    = '10 5 0'
  []

  # -- y = 5.8 m: near secondary maximum --
  # Analytic: sin(theta) ~ 1/2 -> theta ~ 30 deg, y ~ 10*tan(30deg) = 5.77 m
  # Expected: I ~ 0.045 * I_central (small secondary lobe)
  [E_intensity_y5p8]
    type     = PointValue
    variable = E_intensity
    point    = '10 5.8 0'
  []

  # -- y = 7.0 m: second null region --
  # Second null: sin(theta) = 2/3 -> theta ~ 41.8 deg, y ~ 10*tan(41.8) ~ 8.9 m
  # Near-null region
  [E_intensity_y7p0]
    type     = PointValue
    variable = E_intensity
    point    = '10 7 0'
  []

  # ==================================================================
  # NEAR-FIELD DIAGNOSTICS at x = 1 m (one wavelength from slit)
  # ==================================================================
  # Near the slit (Fresnel regime: x ~ a^2/lambda_0 = 9 m, so x=1 is
  # far in the Fresnel near field). The field here is the direct superposition
  # of Huygens secondary wavelets from all points of the slit, not yet
  # collapsed into the sinc^2 envelope. Sampling at x=1 illustrates the
  # transition from near to far field.
  [E_intensity_near_y0]
    type     = PointValue
    variable = E_intensity
    point    = '1 0 0'
  []
  [E_intensity_near_y1p5]
    type     = PointValue
    variable = E_intensity
    point    = '1 1.5 0'
  []

  # ==================================================================
  # MID-FIELD DIAGNOSTICS at x = 5 m
  # ==================================================================
  # x = 5 m: intermediate regime (Fresnel, transitioning to Fraunhofer).
  # The sinc^2 envelope begins to form but is not yet fully developed.
  [E_intensity_mid_y0]
    type     = PointValue
    variable = E_intensity
    point    = '5 0 0'
  []
  [E_intensity_mid_y3p5]
    type     = PointValue
    variable = E_intensity
    point    = '5 3.5 0'
  []

  # ==================================================================
  # INTEGRATED POWER
  # ==================================================================
  # Total field energy in the domain: proportional to integral |E|^2 dA.
  # This gives a global measure of field level and can be used to
  # verify that the ABC is working (power injected = power absorbed).
  [total_field_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = E_intensity
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner.
  # The E_real and E_imag variables are coupled through the EMRobinBC
  # terms on the right, top, and bottom boundaries. These BCs introduce
  # off-diagonal Jacobian blocks between E_real and E_imag (from the
  # j k0 factor in the Sommerfeld ABC: the term j k0 E couples the
  # imaginary part of E_real equation to E_imag and vice versa).
  #
  # Without full=true, the SMP ignores these off-diagonal blocks and
  # Newton convergence degrades. With full=true, the exact Jacobian
  # is assembled and Newton converges in 1-2 iterations (the system
  # is linear, so Newton converges in exactly 1 step with a correct
  # Jacobian and LU factorisation).
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation via PETSc.
  # The 2D Helmholtz system at this resolution has:
  #   150 x 160 QUAD4 elements = 24,000 elements
  #   (151 x 161) x 2 variables = 48,722 DOFs (FIRST-order Lagrange nodes)
  # This is comfortably within the range where LU is efficient in serial.
  #
  # LU is preferred over iterative solvers because:
  #   (1) The Helmholtz operator can be indefinite (k0^2 term changes sign
  #       of eigenvalue spectrum), which stalls GMRES without good preconditioning.
  #   (2) The system is linear, so Newton converges in exactly 1 step with LU.
  #   (3) The off-diagonal Robin BC coupling makes the system non-symmetric,
  #       ruling out standard conjugate-gradient methods.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: full 2D field — E_real, E_imag, E_intensity on the 150x160 mesh.
  # Visualise in ParaView:
  #   - Colour by E_intensity: reveals the classic single-slit diffraction
  #     pattern — a bright central stripe at y=0 flanked by weaker fringes.
  #   - The slit at x=0 is visible as a bright stripe of width a=3m at y=0.
  #   - The shadow regions (|y|>1.5 at x near 0) gradually acquire field
  #     amplitude due to diffraction as x increases.
  #   - At x=10-15, the sinc^2 envelope is clearly visible with nulls at
  #     y ~ +/-3.5 m and secondary maxima at y ~ +/-5.8 m.
  #   - Use the "Plot Over Line" filter in ParaView along x=10 to extract
  #     the intensity vs. y profile and compare with the analytic sinc^2 curve.
  exodus = true

  # csv: all postprocessor values — intensity samples at 13 y-positions
  # along x=10, near-field probes at x=1 and x=5, and total field energy.
  # The key result: E_intensity_y0 >> E_intensity_y3p5 (null) and
  # E_intensity_y5p8 / E_intensity_y0 ~ 0.045 (secondary lobe ratio).
  csv = true
[]
