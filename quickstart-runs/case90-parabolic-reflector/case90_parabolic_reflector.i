# ============================================================
# Case 90: Parabolic Reflector — Paraxial Focusing of a Plane Wave
# Born & Wolf, Principles of Optics, 7th Ed. (1999), Ch. 5
# Balanis, "Advanced Engineering Electromagnetics", 2nd Ed. (2012), Ch. 12
# Jin, "The Finite Element Method in Electromagnetics", 3rd Ed. (2014), Ch. 9
#
# MIT 6.635 Advanced Electromagnetism (Prof. Jin Au Kong, Spring 2003)
# Lecture themes: scattering, aperture antennas, focusing systems
#
# ============================================================
# HISTORICAL AND PHYSICAL BACKGROUND
# ============================================================
#
# A parabolic reflector antenna is the canonical device for converting
# a plane wave into a focused beam (or equivalently, for converting a
# point source at the focus into a collimated beam). The parabola
# y² = 4 f x (in 2D) has the geometric property that every ray
# arriving parallel to the axis of symmetry is reflected through
# the focal point at (f, 0).
#
# In electromagnetic terms (2D TE-polarised Helmholtz):
#
#   An incident plane wave E_inc = exp(-j k₀ x)  (travelling in the
#   -x direction) illuminates a parabolic conducting surface. The
#   surface reflects the wave, and constructive interference of the
#   reflected rays produces a high-intensity focus at x = f.
#
# The physics is captured by two ingredients:
#
#   (1) GEOMETRY: the parabola x = y²/(4f) maps every incoming ray
#       (in the -x direction) to a reflected ray that passes through
#       the focal point (f, 0).
#
#   (2) PHASE COHERENCE: all ray paths from the plane wavefront to
#       the focus have EQUAL total path length (this is the definition
#       of a parabola from its optical-path-length property). Equal
#       path lengths → in-phase arrival → constructive interference.
#
# ============================================================
# IMPLEMENTATION STRATEGY — ABSORBING REFLECTOR MATERIAL
# ============================================================
#
# The cleanest FEM approach for a curved-surface reflector on a
# structured rectangular mesh is to embed the reflector as a LOSSY
# DIELECTRIC LAYER that effectively behaves as a PEC (Perfect Electric
# Conductor) for the wave:
#
#   ε_r = ε_r' + j ε_r''   with ε_r'' ≫ 1  (high-loss material)
#
# A large imaginary permittivity ε_r'' creates a very short skin depth
#   δ_s = 1/(k₀ √(ε_r''/2)) → 0  as ε_r'' → ∞
# The wave is absorbed/reflected at the surface without penetrating
# deeply into the reflector material.
#
# For ε_r'' = 50 and k₀ = 2π ≈ 6.283:
#   δ_s ≈ 1/(k₀ √(ε_r''/2)) = 1/(6.283 × √25) ≈ 1/(6.283×5) ≈ 0.032
# which is much smaller than the wavelength λ₀ = 1.0, so the absorber
# effectively mimics a PEC reflector surface.
#
# ============================================================
# GEOMETRY
# ============================================================
#
#   Domain: [-2, 10] × [-5, 5]   (12 × 10 units; λ₀ = 1 unit)
#
#   Parabola:  x_parabola(y) = y²/(4f)  with f = 2.0 (focal length)
#   Aperture:  |y| < 3.0  (reflector half-width = 3λ₀)
#   At |y| = 3: x_parabola = 9/8 = 1.125
#
#   Reflector layer: x < y²/(4f)  AND  |y| < 3  AND  x < 1.5
#   (the condition x < 1.5 caps the region on the right to avoid
#    filling the entire left-half domain; the parabola only extends
#    to x ≈ 1.125 at |y| = 3)
#
#   Incident wave:  E_inc = exp(-j k₀ x) — plane wave in -x direction.
#   Injected via EMRobinBC at x = 10 (right boundary).
#
#   Focal point:  (f, 0) = (2, 0)
#
#   Absorbing BCs on left (x = -2), top (y = 5), bottom (y = -5).
#   Note: the wave travels in the -x direction so the RIGHT boundary
#   is the wave injection port and the LEFT boundary is the absorber
#   for the focused reflected wave.
#
#        x=-2       x=0        x=2        x=10
#          |                   * ← focus   |
#          |       /|          |            |← plane wave
#          |      / |          |           E_inc = e^{-jk₀x}
#          |     /  |          |            |
#          |    /   | ← parabolic           |
#  ABC     |   | reflector  |              Port (injection)
#          |    \   |          |            |
#          |     \  |          |            |
#          |      \ |          |            |
#          |       \|          |            |
#
# ============================================================
# GOVERNING EQUATION — 2D SCALAR HELMHOLTZ
# ============================================================
#
# For TE polarisation (E = E_z ẑ, ∂/∂z = 0):
#
#   ∇²E_z + k₀² ε_r(x,y) E_z = 0
#
# where ε_r(x,y) is piecewise:
#   ε_r = 1            in free space
#   ε_r = 1 + j × 50  inside the parabolic reflector layer
#
# The complex permittivity in the lossy reflector is handled by
# splitting E_z = E_real + j E_imag and the equation into real/imaginary
# components:
#
#   E_real:  ∇²E_r + k₀²ε_r' E_r − k₀²ε_r'' E_i = 0
#   E_imag:  ∇²E_i + k₀²ε_r' E_i + k₀²ε_r'' E_r = 0
#
# (Signs from splitting ε_r = ε_r' + j ε_r'' and matching real/imag:
#  k₀²ε_r E_z = k₀²(ε_r' + j ε_r'')(E_r + j E_i)
#             = k₀²(ε_r' E_r − ε_r'' E_i) + j k₀²(ε_r' E_i + ε_r'' E_r))
#
# ============================================================
# WAVE INJECTION — RIGHT BOUNDARY
# ============================================================
#
# The plane wave E_inc = exp(-j k₀ x) travels in the -x direction.
# At the right boundary (x = 10), the EMRobinBC port condition injects
# this wave. The EMRobinBC injects exp(+j k₀ x) in the +x direction by
# default (port mode). To inject a wave in the -x direction we place
# the injection port on the RIGHT boundary: at x = 10, the outward
# normal is +x̂, so a "port wave" exp(+j k₀ × 10) propagates inward
# (leftward) into the domain.
#
# Effectively: the right-boundary EMRobinBC injects exp(-j k₀ x) from
# the right, propagating toward decreasing x. This is the standard
# configuration for a parabolic reflector illuminated from the far right.
#
# ============================================================
# NUMERICAL VALUES
# ============================================================
#
# λ₀ = 1.0 unit  →  k₀ = 2π = 6.2831853... rad/unit
# k₀² = 4π² = 39.4784176... unit⁻²
#
# Free space (no reflector):
#   ε_r'  = 1.0
#   ε_r'' = 0.0
#   k₀²ε_r' = 39.4784176  (reaction coefficient)
#   k₀²ε_r'' = 0.0         (no cross-coupling)
#
# Inside parabolic reflector (x < y²/8  AND  |y| < 3  AND  x < 1.5):
#   ε_r'  = 1.0
#   ε_r'' = 50.0
#   k₀²ε_r'  = 39.4784176   (same real part as free space)
#   k₀²ε_r'' = 50 × 39.4784176 / 39.4784176 × ... = k₀² × 50
#             = 39.4784176 × 50 = 1973.9208...  (large cross-coupling)
#
# Focal length:  f = 2.0 (focal point at (2, 0))
# 4f = 8.0  (parabola:  x = y²/8)
#
# ============================================================
# WEAK FORM AND KERNEL SIGN CONVENTIONS
# ============================================================
#
# MOOSE kernel sign conventions:
#
#   Diffusion(variable = E):
#     Residual: +∫ ∇E · ∇v dV
#     Strong form equivalent: −∇²E
#
#   ADMatReaction(reaction_rate = r):
#     Residual: −r × ∫ E v dV
#     Strong form equivalent: −r × E
#
#   ADMatCoupledForce(v = u, mat_prop_coef = c):
#     Residual: −c × ∫ u × test dV
#     Strong form equivalent: −c × u
#
# Assembling E_real equation (strong form: ∇²E_r + k₀²ε_r'E_r − k₀²ε_r''E_i = 0):
#   Diffusion(E_r)                           → −∇²E_r
#   ADMatReaction(k0sq_eps_pr)               → −k₀²ε_r' E_r
#   ADMatCoupledForce(E_i, neg_k0sq_eps_pp)  → −(−k₀²ε_r'') E_i = +k₀²ε_r'' E_i
#   Residual = 0 → strong form: −∇²E_r − k₀²ε_r'E_r + k₀²ε_r''E_i = 0
#   ↔ ∇²E_r + k₀²ε_r'E_r − k₀²ε_r''E_i = 0  ✓
#
# Assembling E_imag equation (strong form: ∇²E_i + k₀²ε_r'E_i + k₀²ε_r''E_r = 0):
#   Diffusion(E_i)                           → −∇²E_i
#   ADMatReaction(k0sq_eps_pr)               → −k₀²ε_r' E_i
#   ADMatCoupledForce(E_r, k0sq_eps_pp)      → −k₀²ε_r'' E_r
#   Residual = 0 → strong form: −∇²E_i − k₀²ε_r'E_i − k₀²ε_r''E_r = 0
#   ↔ ∇²E_i + k₀²ε_r'E_i + k₀²ε_r''E_r = 0  ✓
#
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...} throughout)
# -----------------------------------------------------------
# k = k₀ = 2π rad/unit  (λ₀ = 1.0 unit)
# f = focal length = 2.0 units
# E0 = incident field amplitude = 1 V/m
# theta = 0 (normal incidence on the right port boundary)
k     = 6.283185307179586    # k₀ = 2π [rad/unit], λ₀ = 1.0 unit
f     = 2.0                  # focal length [units]; focus at (f, 0) = (2, 0)
E0    = 1                    # incident wave amplitude [V/m]
theta = 0                    # incidence angle [degrees] for cosTheta function

[Mesh]
  # 2D rectangular domain [-2, 10] × [-5, 5].
  # Total size: 12 × 10 = 120 square units.
  # λ₀ = 1.0 unit; mesh resolution: 120 × 100 = 12,000 quad elements.
  # Element size: 12/120 = 0.1 unit in x, 10/100 = 0.1 unit in y.
  # Elements per wavelength: λ₀/dx = 1.0/0.1 = 10 per λ₀.
  # This is the minimum for first-order Lagrange; fine enough for the
  # smooth field away from the reflector, adequate for the PEC-like layer.
  #
  # The reflector occupies the region x < y²/8, |y| < 3, x < 1.5.
  # At the parabola tip (y=0): x = 0.
  # At the aperture edges (|y|=3): x = 9/8 = 1.125.
  # The reflector is roughly 10–15 elements thick (0.1–0.15 λ₀ → ~1 skin depth).
  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 120
    ny   = 100
    xmin = -2
    xmax = 10
    ymin = -5
    ymax = 5
  []

  # Rename boundaries to compass directions for clarity.
  # GeneratedMeshGenerator produces: left, right, bottom, top.
  [rename]
    type         = RenameBoundaryGenerator
    input        = domain
    old_boundary = 'left right bottom top'
    new_boundary = 'west east south north'
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor E_z = E_real + j E_imag.
  # The total field is the solution to the 2D Helmholtz equation with the
  # parabolic absorbing reflector embedded as a high-loss dielectric region.
  # The incident plane wave enters from the right (east boundary) via EMRobinBC.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex phasor.
  # In free space (away from the reflector), E_real and E_imag are decoupled
  # in the bulk and interact only through the boundary conditions. Inside the
  # lossy reflector layer, they are coupled via the ADMatCoupledForce kernels.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # k0sq_eps_pr_fn: k₀² × Re(ε_r(x,y)) — reaction coefficient for both
  # E_real and E_imag equations.
  #
  # ε_r' = 1.0 everywhere (free space and in the reflector layer, since
  # the reflector uses only a large imaginary part).
  # k₀² × ε_r' = (2π)² × 1.0 = 4π² = 39.47841760435743 unit⁻²
  #
  # The reaction coefficient is spatially uniform for the REAL part of ε_r.
  # ------------------------------------------------------------------
  [k0sq_eps_pr_fn]
    type       = ParsedFunction
    expression = '39.47841760435743'
  []

  # ------------------------------------------------------------------
  # k0sq_eps_pp_fn: +k₀² × Im(ε_r(x,y)) — cross-coupling from E_real
  # into E_imag (positive sign in the E_imag equation).
  #
  # The lossy reflector layer is defined by:
  #   x < y²/(4f) = y²/8.0   (interior side of parabola)
  #   AND |y| < 3.0            (within aperture)
  #   AND x < 1.5              (cap on the right to exclude far-right domain)
  #   AND x > -1.5             (cap on the left to confine to near the vertex)
  #
  # The condition x < y²/8 combined with x < 1.5 naturally captures the
  # entire parabolic surface region (since the parabola reaches x = 1.125
  # at |y| = 3, the right cap x < 1.5 has no effect except preventing the
  # region from spilling past the parabola tip at large |y| outside aperture).
  #
  # Inside the reflector: ε_r'' = 50 → k₀² × 50 = 1973.920880...
  # Outside:              ε_r'' = 0  → coefficient = 0
  # ------------------------------------------------------------------
  [k0sq_eps_pp_fn]
    type       = ParsedFunction
    expression = 'if(x < y*y/(4.0*${f}) & abs(y) < 3.0 & x < 1.5 & x > -1.5, 1973.9208802178714, 0.0)'
  []

  # ------------------------------------------------------------------
  # neg_k0sq_eps_pp_fn: −k₀² × Im(ε_r(x,y)) — cross-coupling from E_imag
  # into E_real (negative sign in the E_real equation).
  #
  # The E_real equation has the NEGATIVE cross-coupling term (see sign
  # convention derivation above). This function returns the negative of
  # k0sq_eps_pp_fn and is used in the ADMatCoupledForce kernel for E_real.
  # ------------------------------------------------------------------
  [neg_k0sq_eps_pp_fn]
    type       = ParsedFunction
    expression = 'if(x < y*y/(4.0*${f}) & abs(y) < 3.0 & x < 1.5 & x > -1.5, -1973.9208802178714, 0.0)'
  []

  # ------------------------------------------------------------------
  # cosTheta: cos(theta) for the EMRobinBC boundary condition.
  # At theta = 0 (normal incidence), cosTheta = cos(0) = 1.0.
  # This factor enters the effective wavenumber k_eff = k₀ cos(θ) in
  # the first-order absorbing/port boundary condition.
  # ------------------------------------------------------------------
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # ------------------------------------------------------------------
  # Wrap the three ParsedFunctions as AD material properties.
  # ADGenericFunctionMaterial evaluates each function at each quadrature
  # point and returns an ADReal, enabling automatic differentiation for
  # the Newton Jacobian. This is essential for the spatially discontinuous
  # material properties at the reflector boundary.
  # ------------------------------------------------------------------

  # k₀² Re(ε_r): uniform reaction coefficient (same everywhere)
  [k0sq_eps_pr_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'k0sq_eps_pr'
    prop_values = 'k0sq_eps_pr_fn'
  []

  # +k₀² Im(ε_r): cross-coupling E_real → E_imag (nonzero inside reflector)
  [k0sq_eps_pp_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'k0sq_eps_pp'
    prop_values = 'k0sq_eps_pp_fn'
  []

  # −k₀² Im(ε_r): cross-coupling E_imag → E_real (negative, inside reflector)
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
  #
  # Three kernels: Diffusion + ADMatReaction + ADMatCoupledForce
  # The ADMatCoupledForce uses NEGATIVE k₀²ε_r'' to subtract the
  # cross-coupling contribution (see sign convention derivation above).
  # ==================================================================

  # Laplacian term: contributes −∇²E_r to the strong-form residual.
  # With unit diffusivity (= 1), this is the standard scalar Laplacian.
  [diff_real]
    type     = Diffusion
    variable = E_real
  []

  # Reaction term: contributes −k₀²ε_r' E_r to the strong form.
  # Combined with Diffusion: strong form = −∇²E_r − k₀²ε_r' E_r = 0
  # → ∇²E_r + k₀²ε_r' E_r = 0 (correct Helmholtz sign convention).
  [helm_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq_eps_pr
  []

  # Cross-coupling from E_imag: −(−k₀²ε_r'') × E_i = +k₀²ε_r'' E_i
  # in the strong form. Inside the reflector this creates a large
  # damping-like coupling that dissipates the field — the FEM equivalent
  # of the wave being absorbed/reflected at a PEC surface.
  # Outside the reflector (neg_k0sq_eps_pp = 0): no contribution.
  [cross_real]
    type          = ADMatCoupledForce
    variable      = E_real
    v             = E_imag
    mat_prop_coef = neg_k0sq_eps_pp
  []

  # ==================================================================
  # E_imag EQUATION:
  #   ∇²E_i + k₀²ε_r' E_i + k₀²ε_r'' E_r = 0
  #
  # Three kernels: Diffusion + ADMatReaction + ADMatCoupledForce
  # The ADMatCoupledForce uses POSITIVE k₀²ε_r'' (cross-coupling E_real
  # into E_imag with positive sign, consistent with the derivation).
  # ==================================================================

  # Laplacian term for E_imag.
  [diff_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Reaction term for E_imag: same coefficient k₀²ε_r' as E_real.
  [helm_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq_eps_pr
  []

  # Cross-coupling from E_real: −k₀²ε_r'' × E_r in the strong form.
  # In the reflector this adds a large positive coupling from E_real
  # into E_imag, completing the dissipative high-loss behaviour.
  [cross_imag]
    type          = ADMatCoupledForce
    variable      = E_imag
    v             = E_real
    mat_prop_coef = k0sq_eps_pp
  []
[]

[BCs]
  # ==================================================================
  # EAST BOUNDARY (x = 10): WAVE INJECTION PORT
  # ==================================================================
  # The plane wave E_inc = exp(-j k₀ x) travels in the -x direction
  # and enters the domain from the right.
  #
  # The EMRobinBC in port mode at the east boundary (outward normal = +x̂)
  # implements the condition:
  #   ∂E/∂x + j k₀ E = 2 j k₀ E₀ exp(+j k₀ x)
  #
  # At x = 10, the incident wave has accumulated phase:
  #   E_inc(10) = exp(-j k₀ × 10)
  # The port condition correctly injects this wave propagating in the
  # -x direction. profile_func_real = E0 = 1 sets the amplitude.
  #
  # Sign convention: sign = negative, consistent with all previous
  # EM cases (77, 79, 83). The 'negative' sign ensures the correct
  # orientation of the Robin BC residual contribution.
  [port_real]
    type              = EMRobinBC
    variable          = E_real
    boundary          = east
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
    boundary          = east
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
  # WEST BOUNDARY (x = -2): ABSORBING BOUNDARY CONDITION
  # ==================================================================
  # The focused reflected wave propagates in the +x direction (away from
  # the parabola toward the focus at x = 2 and continuing leftward past
  # the focus). The west boundary at x = -2 is 4 units (4λ₀) from the
  # reflector vertex and 4 units past the focus. The first-order ABC
  # ∂E/∂n + j k₀ E = 0 absorbs the outgoing reflected/focused wave.
  [abc_real_west]
    type              = EMRobinBC
    variable          = E_real
    boundary          = west
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = absorbing
  []

  [abc_imag_west]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = west
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
    mode              = absorbing
  []

  # ==================================================================
  # NORTH BOUNDARY (y = +5): ABSORBING BOUNDARY CONDITION
  # ==================================================================
  # The incident wave has its wavefronts perpendicular to x (plane wave
  # in -x direction), so it hits the top and bottom boundaries at zero
  # angle. However, scattered field from the reflector edges and the
  # focused beam diverging past the focus may hit these boundaries at
  # oblique angles. The first-order ABC is adequate for this domain size.
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
  # SOUTH BOUNDARY (y = -5): ABSORBING BOUNDARY CONDITION
  # ==================================================================
  # Symmetric to the north boundary. The domain is symmetric about y = 0
  # (both the incident plane wave and the parabola are symmetric). The
  # first-order ABC on the south boundary absorbs the same field patterns.
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
  # E_intensity: |E_z|² = E_real² + E_imag² — time-averaged field intensity.
  #
  # This is proportional to the cycle-averaged power density (W/m²).
  # Visualising E_intensity in ParaView reveals:
  #   - Incident plane wave: uniform intensity |E_inc|² = 1 from x = 10 left
  #   - Reflector shadow: reduced intensity behind the parabola
  #   - Focal spot at (2, 0): bright intensity peak — constructive interference
  #     of all reflected rays converging to the focal point
  #   - Diverging beam past the focus: intensity decreasing as 1/(x-f) away
  # The ratio of focal peak intensity to incident intensity (|E_focus|²/1)
  # measures the focusing gain: for a 6λ aperture, geometric gain ≈ 4πA/λ²
  # in 3D, or proportional to aperture/λ in 2D (≈ 2×3 = 6 times for this case).
  # ------------------------------------------------------------------
  [E_intensity]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Compute |E_z|² = E_real² + E_imag² nodewise after the Newton solve.
  # ParsedAux is evaluated at each mesh node using the primary field values.
  [intensity_aux]
    type              = ParsedAux
    variable          = E_intensity
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # Field diagnostics to verify parabolic focusing
  #
  # The key observable: intensity at the focal point (2, 0) should be
  # significantly HIGHER than the incident wave intensity (= 1), because
  # the reflector collects power from the full aperture (width 6λ₀) and
  # focuses it to a diffraction-limited spot (width ≈ λ₀).
  # Expected 2D focusing gain ≈ aperture / (π × λ₀/2) ≈ 6/1.57 ≈ 4.
  #
  # Additional probes along the optical axis (y = 0) map the axial
  # intensity distribution and confirm a local maximum at x = f = 2.
  # ==================================================================

  # --- FOCAL POINT (2, 0) ---
  # Primary measurement: intensity peak at the predicted focus.
  # For a perfect parabolic PEC reflector, all reflected rays pass
  # through this point in phase → maximum intensity.
  [E_real_focus]
    type     = PointValue
    variable = E_real
    point    = '2 0 0'
  []

  [E_imag_focus]
    type     = PointValue
    variable = E_imag
    point    = '2 0 0'
  []

  [E_intensity_focus]
    type     = PointValue
    variable = E_intensity
    point    = '2 0 0'
  []

  # --- INCIDENT WAVE CHECK: far from reflector (8, 0) ---
  # At x = 8, the plane wave has not yet hit the reflector.
  # Intensity should equal |E_inc|² = 1.0 (unperturbed incident wave).
  # This verifies that the port injection and absorbing conditions
  # are functioning correctly.
  [E_intensity_incident]
    type     = PointValue
    variable = E_intensity
    point    = '8 0 0'
  []

  # --- PRE-FOCUS (1, 0): converging reflected beam ---
  # At x = 1 (one unit before the focus), the reflected beam is still
  # converging. Intensity here should be below the focal peak but above
  # the incident intensity, indicating convergence.
  [E_intensity_prefocus]
    type     = PointValue
    variable = E_intensity
    point    = '1 0 0'
  []

  # --- POST-FOCUS (-1, 0): diverging beam past the focus ---
  # At x = -1 (one unit past the focus in the -x direction), the beam
  # has passed through the focus and is diverging. Intensity should be
  # falling from the peak, confirming a genuine focal maximum at (2, 0).
  [E_intensity_postfocus]
    type     = PointValue
    variable = E_intensity
    point    = '-1 0 0'
  []

  # --- NEAR-FOCUS OFF-AXIS (2, 1): first dark ring ---
  # The focused spot has a diffraction-limited width ≈ f λ₀ / D where
  # D = 6 is the aperture diameter. At 1 unit off axis the intensity
  # should be significantly less than the peak, characterising the
  # focal spot size.
  [E_intensity_focus_offaxis]
    type     = PointValue
    variable = E_intensity
    point    = '2 1 0'
  []

  # --- REFLECTOR REGION CHECK (0.5, 0): inside absorber ---
  # Near the parabola vertex (y = 0, x ≈ 0), we are inside the
  # absorbing reflector layer. Intensity should be near zero —
  # the high-loss material prevents the field from penetrating.
  [E_intensity_inside_reflector]
    type     = PointValue
    variable = E_intensity
    point    = '0.5 0 0'
  []

  # --- ABOVE APERTURE EDGE (5, 3.5): unilluminated region ---
  # At y = 3.5 (above the aperture limit of |y| = 3), the parabola
  # is not present and the field should be close to the incident
  # plane wave (intensity ≈ 1), unaffected by the reflector.
  [E_intensity_above_aperture]
    type     = PointValue
    variable = E_intensity
    point    = '5 3.5 0'
  []

  # --- TOTAL FIELD ENERGY in the domain ---
  # Global energy diagnostic. The reflector concentrates but does
  # not add energy; this integral is dominated by the incident wave
  # uniformly filling the 12×10 domain + the focal intensity peak.
  [total_field_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = E_intensity
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner — essential for this coupled system.
  #
  # Two sources of off-diagonal Jacobian coupling require full=true:
  #   (1) EMRobinBC at all four boundaries: the j k₀ term couples
  #       E_real and E_imag at every boundary DOF.
  #   (2) ADMatCoupledForce kernels in the reflector region: the
  #       large cross-coupling coefficient k₀²×50 = 1973.9 creates
  #       strong off-diagonal blocks throughout the reflector volume.
  #
  # Without full=true the preconditioner ignores these cross-blocks,
  # and the Newton iteration diverges in the high-loss reflector region
  # where the off-diagonal terms dominate the diagonal.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation via MUMPS.
  # The 2D coupled Helmholtz system has 120×100×2 = 24,000 DOFs.
  # LU is preferred over iterative solvers because:
  #   (1) The Helmholtz operator is indefinite (k₀² term breaks coercivity).
  #   (2) The reflector introduces off-diagonal blocks with magnitude ~1974,
  #       which are orders of magnitude larger than the diagonal blocks
  #       (k₀²ε_r' ≈ 39). This large condition number defeats ILU and GMRES
  #       without careful preconditioning.
  #   (3) The system is linear (E_z appears linearly in all kernels and BCs),
  #       so Newton converges in exactly 1 iteration — LU is used once.
  # MUMPS handles the non-symmetric system (E_real and E_imag are coupled
  # asymmetrically through the cross terms) efficiently.
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: 2D field output (E_real, E_imag, E_intensity) on the full mesh.
  # Visualise in ParaView:
  #   - E_intensity: shows the focusing pattern — bright focal spot at (2,0),
  #     uniform incident wave at x > 3, shadow behind the parabola at x < 0
  #   - E_real, E_imag: phase pattern of the reflected/converging wave
  #   - Use a clip plane at y = 0 (optical axis) and plot E_intensity vs x
  #     to see the axial intensity profile with the focal peak at x = 2
  #   - The parabolic shape of the high-intensity boundary of the reflector
  #     is visible in the E_intensity field as a sharp drop to near-zero
  #     inside the absorbing layer
  exodus = true

  # csv: postprocessor values — intensity at focal point, incident check,
  # pre/post-focus, off-axis, and total domain energy.
  # The key figure of merit: E_intensity_focus / E_intensity_incident
  # gives the focusing gain of the parabolic reflector. For a lossless
  # parabola with aperture D = 6λ₀ and focal length f = 2λ₀, the 2D
  # geometric focusing gain scales as D/λ₀ = 6.
  csv    = true
[]
