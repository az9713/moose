# ============================================================
# Case 91: Radar Cross Section of a Finite Flat Conducting Plate
# Balanis, "Advanced Engineering Electromagnetics", 2nd Ed. (2012), Ch. 11
# Harrington, "Time-Harmonic Electromagnetic Fields" (1961), Ch. 3
# Jin, "The Finite Element Method in Electromagnetics", 3rd Ed. (2014), Ch. 9
# Knott, Shaeffer & Tuley, "Radar Cross Section", 2nd Ed. (2004), Ch. 5
#
# MIT 6.635 Advanced Electromagnetism (Prof. Jin Au Kong, Spring 2003)
# Lecture themes: scattering, EFIE, aperture diffraction, RCS
#
# ============================================================
# PHYSICAL BACKGROUND — RADAR CROSS SECTION
# ============================================================
#
# Radar Cross Section (RCS) σ is the effective area that a target
# presents to an incoming radar wave. It is defined in terms of the
# ratio of scattered power density to incident power density:
#
#   σ = lim_{R→∞} 4πR² |E_scat|² / |E_inc|²   (3D definition)
#
# In the 2D problem (cylinder of infinite extent, unit cross-section),
# the scattering width W replaces σ:
#
#   W = lim_{ρ→∞} 2πρ |E_scat|² / |E_inc|²     (2D scattering width)
#
# For a FLAT PLATE (strip) of width L = 3λ₀ at normal incidence, the
# physical optics prediction for the monostatic (backscatter) RCS is:
#
#   W_monostatic ≈ k₀ L²                  (normal incidence)
#                = 2π × (3λ₀)²/λ₀
#                = 2π × 9 λ₀ = 18π λ₀
#
# For L = 3.0, λ₀ = 1.0: W_mono ≈ 2π × k₀ × L²/4 = π × 3² = 28.27 m²/m
# (exact physical optics; FEM captures diffraction corrections)
#
# The bistatic RCS pattern (scattered power as function of observation
# angle φ) has a main lobe in the specular direction and side lobes
# from diffraction at the plate edges (end-fire directions).
#
# ============================================================
# SCATTERED-FIELD FORMULATION
# ============================================================
#
# The total field is decomposed as:
#   E_z = E_inc + E_scat
#
# Incident plane wave (normal incidence, +x direction):
#   E_inc = E₀ exp(+j k₀ x)
#   E_inc_real(x,y) = E₀ cos(k₀ x)
#   E_inc_imag(x,y) = E₀ sin(k₀ x)
#
# The plate (PEC, conducting strip) is modelled as a thin high-loss
# dielectric strip (same approach as Case 90):
#   |x| < 0.1  AND  |y| < L/2 = 1.5
#   ε_r = 1 + j × 100  (very high loss → skin depth ≈ 0.016 λ₀)
#
# Scattered field equation:
#   ∇²E_scat + k₀² ε_r(x,y) E_scat = 0  (free space surrounding the plate)
#
# In the total-field approach, we solve for E_z directly (total field),
# inject the incident wave from the west boundary, and absorb on east,
# north, south. Then compute the scattered field analytically:
#   E_scat_real(x,y) = E_real(x,y) − cos(k₀ x)
#   E_scat_imag(x,y) = E_imag(x,y) − sin(k₀ x)
#   |E_scat|² = E_scat_real² + E_scat_imag²
#
# ============================================================
# GEOMETRY
# ============================================================
#
#   Domain: [-8, 8] × [-8, 8]  (16 × 16 units; λ₀ = 1.0 unit)
#   Plate:  |x| < 0.1, |y| < 1.5  (width L = 3.0 = 3λ₀, thin strip)
#   Incident wave: plane wave from west (E₀ = 1, normal incidence)
#
#   y
#   ^
#   |   8
#   |   .  .  .  .  .  .  .  .
#   |   .                    .
#   |   .       ┃ ← plate    .
#   |   .       ┃            .  → E_inc = exp(+jk₀x)
#   |   .       ┃            .
#   |   .  .  .  .  .  .  .  .
#   |  -8           0        8 → x
#
#   Observation circle at R = 5 (well within the domain boundary at 8).
#   Angles φ = 0°, 30°, 60°, 90°, 120°, 150°, 180° on the circle:
#     φ = 0°:   (5, 0)         forward scatter
#     φ = 90°:  (0, 5)         broadside left
#     φ = 180°: (-5, 0)        backscatter
#
# ============================================================
# GOVERNING EQUATIONS — COMPLEX HELMHOLTZ WITH LOSSY PLATE
# ============================================================
#
# For TE polarisation (E = E_z ẑ):
#   ∇²E_z + k₀² ε_r(x,y) E_z = 0
#
# The complex ε_r = ε_r' + j ε_r'' leads to coupled equations for
# E_real and E_imag. Following the same sign convention as Case 90:
#
#   E_real:  ∇²E_r + k₀²ε_r' E_r − k₀²ε_r'' E_i = 0
#   E_imag:  ∇²E_i + k₀²ε_r' E_i + k₀²ε_r'' E_r = 0
#
# Numerical values:
#   k₀ = 2π = 6.2831853...  (λ₀ = 1.0)
#   k₀² = 4π² = 39.4784176...
#
#   Free space: ε_r' = 1.0,  ε_r'' = 0.0  → k₀² ε_r' = 39.478,  k₀² ε_r'' = 0
#   Plate:      ε_r' = 1.0,  ε_r'' = 100  → k₀² ε_r' = 39.478,  k₀² ε_r'' = 3947.84
#
# ============================================================
# SCATTERED FIELD AND RCS EXTRACTION
# ============================================================
#
# After solving for the total field E = E_real + j E_imag, the
# scattered field is computed analytically as:
#
#   E_scat_real = E_real − cos(k₀ x)
#   E_scat_imag = E_imag − sin(k₀ x)
#   |E_scat|² = E_scat_real² + E_scat_imag²
#
# The 2D scattering width (RCS proxy) at observation point (x_obs, y_obs):
#   W(φ) ∝ 2π R × |E_scat(R, φ)|² / |E_inc|²
#         = 2π R × |E_scat|²  (since |E_inc| = 1)
#
# For R = 5 at seven angles, the normalised scattered intensity map:
#   |E_scat(5, φ)|²  for φ ∈ {0°, 30°, 60°, 90°, 120°, 150°, 180°}
#
# Physical optics prediction (normal incidence, plate width L = 3):
#   Main lobe at φ = 0° (forward scatter) and φ = 180° (backscatter)
#   First nulls at sin(φ) = ±1/3 → φ ≈ ±19.5° from forward direction
#   Side-lobe level ≈ −13 dB below the main lobe (sinc² pattern)
#
# ============================================================
# KERNEL AND SIGN CONVENTIONS
# ============================================================
#
# (Same conventions as Case 90 — see detailed derivation there)
#
#   E_real equation:
#     Diffusion(E_r)                  → −∇²E_r  (strong form)
#     ADMatReaction(k0sq_eps_pr)      → −k₀²ε_r' E_r
#     ADMatCoupledForce(E_i, neg_pp)  → +k₀²ε_r'' E_i  (note: neg_pp < 0)
#
#   E_imag equation:
#     Diffusion(E_i)                  → −∇²E_i
#     ADMatReaction(k0sq_eps_pr)      → −k₀²ε_r' E_i
#     ADMatCoupledForce(E_r, pos_pp)  → −k₀²ε_r'' E_r  (note: pos_pp > 0)
#
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...} throughout)
# -----------------------------------------------------------
# k = k₀ = 2π rad/unit  (λ₀ = 1.0 unit)
# L = plate half-width = 1.5 (full width = 3.0 = 3λ₀)
# E0 = incident field amplitude = 1 V/m
# theta = 0 (normal incidence; cosTheta = cos(0) = 1)
k     = 6.283185307179586    # k₀ = 2π [rad/unit], λ₀ = 1.0 unit
E0    = 1                    # incident field amplitude [V/m]
theta = 0                    # incidence angle [degrees]

[Mesh]
  # 2D square domain [-8, 8] × [-8, 8].
  # Mesh: 160 × 160 = 25,600 quad elements.
  # Element size: 16/160 = 0.1 unit = 0.1 λ₀.
  # Elements per wavelength: 10 (minimum for FIRST-order Lagrange).
  # Domain: 16λ₀ × 16λ₀ gives ≥ 5λ₀ clearance from plate edges to the
  # ABC boundaries, minimising spurious reflections from the first-order ABC.
  #
  # The plate at |x| < 0.1 is only 2 elements thick in x — this is
  # marginal but sufficient for the skin-depth model (ε_r'' = 100 gives
  # a very short skin depth δ_s ≈ 1/(k₀√50) ≈ 0.023 λ₀). The field
  # barely penetrates the plate, so resolution inside the plate is not
  # critical to the far-field scattered field accuracy.
  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 160
    ny   = 160
    xmin = -8
    xmax = 8
    ymin = -8
    ymax = 8
  []

  # Rename boundaries to compass-direction names.
  [rename]
    type         = RenameBoundaryGenerator
    input        = domain
    old_boundary = 'left right bottom top'
    new_boundary = 'west east south north'
  []
[]

[Variables]
  # E_real — real part of the complex total electric field phasor E_z.
  # Includes the incident wave (exp(+j k₀ x) real part = cos(k₀ x))
  # plus the scattered field from the conducting plate.
  # First-order Lagrange elements enforce C⁰ continuity at the plate
  # boundary, which is correct for TE polarisation (tangential E_z
  # is continuous across the dielectric-to-conductor interface).
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex total field phasor.
  # The total field at the plate surface is the superposition of the
  # incident wave and the induced scattering; the high loss forces the
  # net field inside the plate toward zero, mimicking PEC behaviour.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # k0sq_eps_pr_fn: k₀² × Re(ε_r(x,y)) — reaction coefficient.
  #
  # ε_r' = 1.0 everywhere (free space and in the conducting plate layer,
  # since the PEC effect is achieved entirely through the imaginary part).
  # k₀² × 1.0 = (2π)² = 4π² = 39.47841760435743 unit⁻²
  #
  # Uniform coefficient — no spatial variation in the real permittivity.
  # ------------------------------------------------------------------
  [k0sq_eps_pr_fn]
    type       = ParsedFunction
    expression = '39.47841760435743'
  []

  # ------------------------------------------------------------------
  # k0sq_eps_pp_fn: +k₀² × Im(ε_r(x,y)) — positive cross-coupling.
  #
  # The flat conducting plate (PEC model) occupies:
  #   |x| < 0.1  AND  |y| < 1.5
  # (thickness 0.2 units = 0.2 λ₀ in x; length L = 3.0 = 3λ₀ in y)
  #
  # Inside the plate: ε_r'' = 100
  #   k₀² × 100 = 39.4784 × 100 = 3947.84177... unit⁻²
  # This very large imaginary permittivity gives:
  #   skin depth δ_s = 1/(k₀ √(ε_r''/2)) = 1/(2π √50) ≈ 0.023 unit
  # The field decays e-fold in 0.023 λ₀, so the plate (thickness 0.2λ₀)
  # is over 8 skin depths thick — effectively an opaque PEC body.
  #
  # Outside the plate: ε_r'' = 0 → no cross-coupling.
  # ------------------------------------------------------------------
  [k0sq_eps_pp_fn]
    type       = ParsedFunction
    expression = 'if(abs(x) < 0.1 & abs(y) < 1.5, 3947.841760435743, 0.0)'
  []

  # ------------------------------------------------------------------
  # neg_k0sq_eps_pp_fn: −k₀² × Im(ε_r(x,y)) — negative cross-coupling.
  #
  # Used in the E_real equation (sign convention: see header derivation).
  # Returns −3947.84 inside the plate, 0 outside.
  # ------------------------------------------------------------------
  [neg_k0sq_eps_pp_fn]
    type       = ParsedFunction
    expression = 'if(abs(x) < 0.1 & abs(y) < 1.5, -3947.841760435743, 0.0)'
  []

  # ------------------------------------------------------------------
  # cosTheta: cos(theta) for EMRobinBC.
  # Normal incidence: theta = 0 → cos(0) = 1.0.
  # ------------------------------------------------------------------
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # ------------------------------------------------------------------
  # Wrap the three ParsedFunctions as AD material properties.
  # ADGenericFunctionMaterial evaluates each function pointwise at each
  # quadrature point. The large jump in k₀²ε_r'' at the plate boundary
  # (0 → 3947.84) is handled correctly by the weak-form integration;
  # automatic differentiation provides the exact Jacobian without any
  # manual derivation.
  # ------------------------------------------------------------------

  # Uniform real reaction coefficient k₀²ε_r' = 39.478
  [k0sq_eps_pr_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'k0sq_eps_pr'
    prop_values = 'k0sq_eps_pr_fn'
  []

  # Positive cross-coupling: +k₀²ε_r'' (into E_imag equation from E_real)
  [k0sq_eps_pp_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'k0sq_eps_pp'
    prop_values = 'k0sq_eps_pp_fn'
  []

  # Negative cross-coupling: −k₀²ε_r'' (into E_real equation from E_imag)
  [neg_k0sq_eps_pp_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'neg_k0sq_eps_pp'
    prop_values = 'neg_k0sq_eps_pp_fn'
  []
[]

[Kernels]
  # ==================================================================
  # E_real EQUATION:
  #   ∇²E_r + k₀²ε_r' E_r − k₀²ε_r'' E_i = 0
  # ==================================================================

  # Laplacian: contributes −∇²E_r to the strong-form residual.
  [diff_real]
    type     = Diffusion
    variable = E_real
  []

  # Reaction: −k₀²ε_r' E_r (uniform coefficient everywhere).
  [helm_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq_eps_pr
  []

  # Cross-coupling: ADMatCoupledForce with neg_k0sq_eps_pp.
  # Residual contribution: −(−k₀²ε_r'') × E_i × test = +k₀²ε_r'' E_i
  # in strong form. Inside the plate this is large and positive, adding
  # the effective "loss" that forces E_real to near-zero inside the conductor.
  [cross_real]
    type          = ADMatCoupledForce
    variable      = E_real
    v             = E_imag
    mat_prop_coef = neg_k0sq_eps_pp
  []

  # ==================================================================
  # E_imag EQUATION:
  #   ∇²E_i + k₀²ε_r' E_i + k₀²ε_r'' E_r = 0
  # ==================================================================

  # Laplacian: contributes −∇²E_i to the strong form.
  [diff_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Reaction: −k₀²ε_r' E_i.
  [helm_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq_eps_pr
  []

  # Cross-coupling: −k₀²ε_r'' × E_r in the strong form.
  # The ADMatCoupledForce with k0sq_eps_pp contributes −k₀²ε_r''×E_r.
  # Inside the plate, this large term forces E_imag to track E_real
  # in the opposite sense, jointly driving both to zero (PEC effect).
  [cross_imag]
    type          = ADMatCoupledForce
    variable      = E_imag
    v             = E_real
    mat_prop_coef = k0sq_eps_pp
  []
[]

[BCs]
  # ==================================================================
  # WEST BOUNDARY (x = -8): WAVE INJECTION PORT
  # ==================================================================
  # The incident plane wave E_inc = exp(+j k₀ x) travels in the +x
  # direction and enters from the left (west). The EMRobinBC in port
  # mode simultaneously injects the incident wave and absorbs any
  # backscattered field leaving the domain through the west face.
  #
  # Jin port condition at x = -8 (outward normal = -x̂):
  #   ∂E/∂n + j k₀ E = 2 j k₀ E₀ exp(-j k₀ × 8)
  # This injects the plane wave E₀ = 1 from the left.
  # profile_func_real = E0 = 1 sets the amplitude.
  [port_real]
    type              = EMRobinBC
    variable          = E_real
    boundary          = west
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = port
  []

  [port_imag]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = west
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = port
  []

  # ==================================================================
  # EAST BOUNDARY (x = +8): ABSORBING BOUNDARY CONDITION
  # ==================================================================
  # The transmitted (forward-scattered) field exits through the east face.
  # mode = absorbing implements ∂E/∂n + j k₀ E = 0.
  # The incident wave also arrives at x = +8 (it propagates in +x),
  # but since we are solving for the TOTAL field (not just scattered field),
  # the east boundary absorbs the total outgoing field.
  [abc_real_east]
    type              = EMRobinBC
    variable          = E_real
    boundary          = east
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = absorbing
  []

  [abc_imag_east]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = east
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = absorbing
  []

  # ==================================================================
  # NORTH BOUNDARY (y = +8): ABSORBING BOUNDARY CONDITION
  # ==================================================================
  # The diffracted field from the plate edges (at y = ±1.5) radiates
  # upward and downward. The north ABC absorbs this edge-diffracted
  # scattered field. The 6.5λ₀ clearance from plate edge to boundary
  # (|y| = 1.5 to y = 8) is adequate for the first-order ABC.
  [abc_real_north]
    type              = EMRobinBC
    variable          = E_real
    boundary          = north
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = absorbing
  []

  [abc_imag_north]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = north
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = absorbing
  []

  # ==================================================================
  # SOUTH BOUNDARY (y = -8): ABSORBING BOUNDARY CONDITION
  # ==================================================================
  # Symmetric to the north boundary. The problem and domain are symmetric
  # about y = 0, so the south and north ABCs see identical field patterns.
  [abc_real_south]
    type              = EMRobinBC
    variable          = E_real
    boundary          = south
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = absorbing
  []

  [abc_imag_south]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = south
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = absorbing
  []
[]

[AuxVariables]
  # ------------------------------------------------------------------
  # E_intensity: |E_total|² = E_real² + E_imag²
  # Total field intensity — the complete (incident + scattered) field.
  # Inside the conducting plate, this should be near zero (PEC condition).
  # Outside, it shows the interference pattern of incident + scattered.
  # ------------------------------------------------------------------
  [E_intensity]
    order  = FIRST
    family = LAGRANGE
  []

  # ------------------------------------------------------------------
  # Es_real: E_scat_real = E_real − cos(k₀ x)
  # Real part of the scattered field, computed analytically by
  # subtracting the known incident wave from the total field.
  # E_inc_real(x,y) = cos(k₀ x) = cos(2π x)
  # This variable isolates the field produced purely by the plate.
  # ------------------------------------------------------------------
  [Es_real]
    order  = FIRST
    family = LAGRANGE
  []

  # ------------------------------------------------------------------
  # Es_imag: E_scat_imag = E_imag − sin(k₀ x)
  # Imaginary part of the scattered field.
  # E_inc_imag(x,y) = sin(k₀ x) = sin(2π x)
  # ------------------------------------------------------------------
  [Es_imag]
    order  = FIRST
    family = LAGRANGE
  []

  # ------------------------------------------------------------------
  # Es_sq: |E_scat|² = Es_real² + Es_imag²
  # Scattered field intensity — the quantity directly related to RCS.
  # The 2D scattering width at observation point (x_obs, y_obs) at
  # distance R from the plate is:
  #   W(φ) = 2π R × Es_sq(x_obs, y_obs)   [valid for R ≫ λ₀]
  # With R = 5 and λ₀ = 1.0, R/λ₀ = 5 (moderate far field; true far
  # field requires R → ∞, but R = 5 gives a good indicator of the
  # scattering pattern).
  # ------------------------------------------------------------------
  [Es_sq]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Total field intensity |E_total|² — all nodes, computed after solve.
  [intensity_aux]
    type              = ParsedAux
    variable          = E_intensity
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []

  # Scattered field real part: subtract E_inc_real = cos(k₀ x) from total.
  # k₀ = 2π = 6.283185307179586. We use the numeric value directly.
  # use_xyzt = true: enables the spatial coordinate x in the expression.
  [Es_real_aux]
    type              = ParsedAux
    variable          = Es_real
    expression        = 'E_real - cos(6.283185307179586 * x)'
    coupled_variables = 'E_real'
    use_xyzt          = true
  []

  # Scattered field imaginary part: subtract E_inc_imag = sin(k₀ x).
  [Es_imag_aux]
    type              = ParsedAux
    variable          = Es_imag
    expression        = 'E_imag - sin(6.283185307179586 * x)'
    coupled_variables = 'E_imag'
    use_xyzt          = true
  []

  # Scattered field intensity |E_scat|².
  [Es_sq_aux]
    type              = ParsedAux
    variable          = Es_sq
    expression        = 'Es_real * Es_real + Es_imag * Es_imag'
    coupled_variables = 'Es_real Es_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # VERIFICATION: incident field on-axis without plate
  # ==================================================================
  # At points far from the plate and in the unperturbed region, the
  # total field should equal the incident wave. At (5, 5) — diagonal
  # from the plate, distance ≈ 7 λ₀ — the plate has negligible effect.
  # |E_total|² should be ≈ 1 (incident wave amplitude).

  [E_intensity_check]
    type     = PointValue
    variable = E_intensity
    point    = '5 5 0'
  []

  # ==================================================================
  # PLATE VERIFICATION: field inside conductor
  # ==================================================================
  # At the plate centre (0, 0), the high loss should drive |E_total|²
  # close to zero (PEC condition: tangential E → 0 on a perfect conductor).
  # The residual |E|² inside the plate is a measure of how well the
  # lossy-dielectric model approximates PEC (should be ≪ 1).

  [E_intensity_plate_centre]
    type     = PointValue
    variable = E_intensity
    point    = '0 0 0'
  []

  # ==================================================================
  # RCS PATTERN: scattered field at R = 5 on observation circle
  # ==================================================================
  # Seven observation angles φ mapping the bistatic RCS pattern.
  # All points are at radius R = 5 from the plate centre (origin).
  # Recall the geometry: plate at x = 0, incident wave in +x direction.
  # - φ = 0° (forward scatter at (5, 0)): transmitted-side main lobe
  # - φ = 90° (broadside left, (0, 5)): side scattering from plate edges
  # - φ = 180° (backscatter at (-5, 0)): specular reflection main lobe

  # Forward scatter φ = 0°: (5, 0)
  [Es_sq_phi_000]
    type     = PointValue
    variable = Es_sq
    point    = '5 0 0'
  []

  # φ = 30°: (5 cos30°, 5 sin30°) = (4.330, 2.500)
  [Es_sq_phi_030]
    type     = PointValue
    variable = Es_sq
    point    = '4.330127 2.5 0'
  []

  # φ = 60°: (5 cos60°, 5 sin60°) = (2.500, 4.330)
  [Es_sq_phi_060]
    type     = PointValue
    variable = Es_sq
    point    = '2.5 4.330127 0'
  []

  # φ = 90° (broadside): (0, 5) — edge-diffracted field
  [Es_sq_phi_090]
    type     = PointValue
    variable = Es_sq
    point    = '0 5 0'
  []

  # φ = 120°: (-2.500, 4.330)
  [Es_sq_phi_120]
    type     = PointValue
    variable = Es_sq
    point    = '-2.5 4.330127 0'
  []

  # φ = 150°: (-4.330, 2.500)
  [Es_sq_phi_150]
    type     = PointValue
    variable = Es_sq
    point    = '-4.330127 2.5 0'
  []

  # Backscatter φ = 180°: (-5, 0) — monostatic RCS direction
  # Physical optics predicts maximum backscatter for a flat plate at
  # normal incidence. The scattering width:
  #   W_180° ≈ k₀ L² / (4π) = 2π × 9 / (4π) = 4.5 unit
  # At R = 5: |E_scat|² ≈ W / (2π R) = 4.5 / (2π × 5) ≈ 0.143
  [Es_sq_phi_180]
    type     = PointValue
    variable = Es_sq
    point    = '-5 0 0'
  []

  # ==================================================================
  # SYMMETRIC ANGLE CHECKS (negative-y hemisphere)
  # ==================================================================
  # The problem is symmetric about y = 0 (symmetric plate, symmetric
  # incident wave). These negative-y probes verify the symmetry:
  # Es_sq at (x, -y) should equal Es_sq at (x, +y).

  # φ = -30° (= 330°): (4.330, -2.500)
  [Es_sq_phi_330]
    type     = PointValue
    variable = Es_sq
    point    = '4.330127 -2.5 0'
  []

  # φ = -90° (= 270°): (0, -5)
  [Es_sq_phi_270]
    type     = PointValue
    variable = Es_sq
    point    = '0 -5 0'
  []

  # ==================================================================
  # TOTAL SCATTERED POWER (integrated over the domain)
  # ==================================================================
  # ∫ |E_scat|² dA — global measure of total scattered energy.
  # For a fixed incident power (|E_inc| = 1), this integral captures
  # the total scattering cross-section (over all angles and both
  # hemispheres). Comparing across frequency or plate width validates
  # the power balance (should scale as L²/λ₀ for flat plate).
  [total_scattered_power]
    type     = ElementIntegralVariablePostprocessor
    variable = Es_sq
  []

  # Total field intensity integral — monitors overall field level.
  [total_field_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = E_intensity
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner — required for this coupled system.
  #
  # Off-diagonal coupling sources:
  #   (1) EMRobinBC at four boundaries: j k₀ coupling between E_real
  #       and E_imag at all boundary DOFs.
  #   (2) ADMatCoupledForce kernels in the plate region (|x|<0.1, |y|<1.5):
  #       the enormous coefficient k₀²ε_r'' = 3947.84 dominates the local
  #       residual. Without full=true the diagonal-only preconditioner
  #       sees only the Helmholtz diagonal (≈39.5) and completely misses
  #       the plate coupling (≈3948), leading to Newton divergence.
  #
  # With full=true, the complete off-diagonal Jacobian is assembled and
  # passed to the LU solver, which handles the high-contrast system exactly.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU via MUMPS.
  # The 2D total-field system has 160×160×2 = 51,200 DOFs.
  # LU is chosen because:
  #   (1) The Helmholtz operator is indefinite (k₀² > 0, mixed-sign eigenvalues).
  #   (2) The plate region creates extremely large off-diagonal Jacobian entries
  #       (k₀²ε_r'' = 3947.84), making iterative Krylov solvers unreliable
  #       without specialised preconditioning (shifted Laplacian, etc.).
  #   (3) The system is linear — Newton converges in 1 iteration. LU cost
  #       is paid once and amortised over all applications.
  # MUMPS efficiently handles the sparse 51k×51k system.
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: 2D field output (E_real, E_imag, E_intensity, Es_real, Es_imag, Es_sq)
  # on the full mesh. Visualise in ParaView:
  #   - E_intensity: total field pattern — shadow region behind the plate
  #     (reduced intensity at x > 0.1, |y| < 1.5), standing-wave fringes
  #     on the incident side (x < -0.1) from incident + backscattered superposition
  #   - Es_sq: scattered field intensity — shows the bistatic RCS pattern.
  #     Main lobes in forward (φ=0°) and backward (φ=180°) directions,
  #     edge-diffraction lobes at broadside (φ=±90°)
  #   - Use a Threshold filter on E_intensity < 0.1 to visualise the shadow
  #     and the plate interior (near-zero field inside the conductor)
  #   - Plot Es_sq along a circle of radius R=5 to directly visualise the RCS
  #     angular pattern; compare the seven postprocessor values
  #   - The plate edges (y = ±1.5) are visible as bright spots in Es_sq
  #     due to edge diffraction — the Keller GTD (Geometrical Theory of
  #     Diffraction) contribution to the RCS pattern
  exodus = true

  # csv: RCS postprocessor values at seven observation angles.
  # The key RCS results:
  #   Es_sq_phi_180: monostatic (backscatter) RCS proxy — largest for normal incidence
  #   Es_sq_phi_000: forward scatter (bistatic RCS at 0°)
  #   Es_sq_phi_090 and Es_sq_phi_270: broadside edge diffraction (should be equal)
  # Physical optics predicts the backscatter lobe width (first null) at:
  #   φ_null = pi - arcsin(lambda_0 / L) = pi - arcsin(1/3) ≈ pi - 19.5° = 160.5°
  csv    = true
[]
