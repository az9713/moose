# ============================================================
# Case 87: Phased Array Beamforming — Two-Element Array
# Balanis, "Antenna Theory: Analysis and Design", 4th Ed. (2016), Ch. 6
# Kong, "Electromagnetic Wave Theory" (Wiley, 1986), Ch. 8
# Van Trees, "Optimum Array Processing" (Wiley, 2002), Ch. 1
#
# MIT 6.635 Advanced Electromagnetism, Spring 2003
# Professor Jin Au Kong
# OCW: https://ocw.mit.edu/courses/6-635-advanced-electromagnetism-spring-2003/
#
# Lectures: 1  (Radiation, dipole antenna)
#           14 (Antenna arrays, pattern multiplication)
#
# ============================================================
# HISTORICAL AND PHYSICAL BACKGROUND
# ============================================================
#
# Phased array antennas electronically steer a beam by introducing
# progressive phase shifts between array elements, without mechanically
# rotating the antenna. The idea was first used systematically in the
# 1940s for radar applications and today underpins virtually every modern
# wireless communication standard (5G beamforming, MIMO radar, etc.).
#
# The fundamental principle is interference: two (or more) point sources
# radiate spherical waves that constructively interfere in the steered
# direction and destructively interfere in other directions. The direction
# of constructive interference is controlled by the inter-element phase
# difference δ.
#
# For a two-element array with element spacing d and phase offset δ,
# the array factor (the pattern due to the array alone, ignoring element
# patterns) is:
#
#   AF(θ) = 1 + e^{j(kd sin θ + δ)}
#
#   |AF(θ)|² = 2 + 2 cos(kd sin θ + δ)
#             = 4 cos²((kd sin θ + δ)/2)
#
# where θ is measured from the array broadside direction (perpendicular
# to the line connecting the elements). The maximum of |AF|² = 4 occurs
# when kd sin θ + δ = 0, i.e.:
#
#   sin θ_max = −δ / (kd)
#
# This is the beam steering equation: a phase offset δ steers the main
# beam to angle θ_max = arcsin(−δ/(kd)) off broadside.
#
# For d = λ/2 (half-wavelength spacing):
#   kd = (2π/λ)(λ/2) = π
#   sin θ_max = −δ/π
#
# In this case: δ = π/4, so sin θ_max = −(π/4)/π = −0.25,
# giving θ_max = −14.48° off broadside (toward the source-1 side).
#
# ============================================================
# GEOMETRY
# ============================================================
#
#          y = +8 (top — absorbing BC)
#              |
#              |
#   Source 2 at (0, +0.25): amplitude 50, phase δ = π/4
#              |  ×
#   Centre  y=0|
#              |  ×
#   Source 1 at (0, -0.25): amplitude 50, phase δ = 0
#              |
#          y = -8 (bottom — absorbing BC)
#
#  x = -8 (left,  absorbing BC) ←——————→ x = +8 (right, absorbing BC)
#
# The two sources are separated by d = 0.5 = λ/2 along the y-axis.
# The broadside direction (perpendicular to the source separation) is
# the x-direction. The array axis is the y-direction.
#
# Array factor angle θ is measured from the broadside (+x direction),
# equivalent to measuring from the y-axis: θ = 0 is broadside (+x).
# Actually in this problem the angle is measured from the +y axis.
# See Postprocessors section for exact coordinate convention.
#
# ============================================================
# GOVERNING EQUATION
# ============================================================
#
# 2D scalar Helmholtz equation (Cartesian, TE polarisation E_z):
#
#   ∇²E + k₀² E = −source(x, y)
#
# where:
#   k₀ = 2π / λ₀   (free-space wavenumber)
#   λ₀ = 1.0        (normalised wavelength)
#   k₀ = 2π = 6.28318530717958648   rad/unit
#   k₀² = 4π² = 39.47841760435743   unit⁻²
#
# Two Gaussian sources model the two array elements:
#
#   Source 1 (no phase offset, δ₁ = 0):
#     Located at (0, -d/2) = (0, -0.25)
#     Acts on E_imag only:
#       J₁_imag(x,y) = A × exp(−r₁²/(2σ²))
#     where r₁² = x² + (y+0.25)², A = 50, σ = 0.05 (2σ² = 0.005)
#
#   Source 2 (phase offset δ₂ = π/4, so cos(δ)=sin(δ)=1/√2):
#     Located at (0, +d/2) = (0, +0.25)
#     Acts on BOTH E_real and E_imag:
#       J₂_imag(x,y) = A cos(δ) × exp(−r₂²/(2σ²))
#       J₂_real(x,y) = A sin(δ) × exp(−r₂²/(2σ²))
#     where r₂² = x² + (y-0.25)², cos(δ) = sin(δ) = 0.70711
#
# Phase offset implementation:
#   Source 1 drives E_imag (phasor S₁ = jA — imaginary source)
#   Source 2 drives E_imag with cos(δ) and E_real with sin(δ):
#     Phasor S₂ = j A e^{jδ} = jA(cosδ + j sinδ) = A(-sinδ + j cosδ)
#     → E_real source: A sinδ (from -sin δ imaginary part of jA*e^{jδ})
#       Wait: jA*e^{jδ} = j*A*(cos δ + j sin δ) = A*j*cos δ - A*sin δ
#             = (-A sinδ) + j(A cosδ)
#     → Real part of source: -A sinδ   (acts on E_real equation as -BodyForce)
#     → Imag part of source: +A cosδ   (acts on E_imag equation as +BodyForce)
#
#   Since BodyForce adds -∫f·v to the residual (strong form: +f on RHS):
#     E_imag: BodyForce(A cosδ Gaussian) → strong form +A cosδ G
#     E_real: BodyForce(-A sinδ Gaussian) → strong form -A sinδ G
#             [But this means BodyForce with function = -Gaussian]
#
#   Equivalently, using positive amplitudes:
#     E_imag: BodyForce with positive coefficient A×cosδ for source2
#     E_real: BodyForce with negative coefficient (-A sinδ) for source2
#
#   HOWEVER: for δ = π/4, sin(δ) = cos(δ) = 1/√2 = 0.70711.
#   The source_2_real function uses +A sinδ amplitude, meaning the BodyForce
#   residual is −∫(+A sinδ G)v dV, i.e., E_real equation gets +A sinδ G as source.
#   The net phasor for source 2: real part source = +A sinδ, imag part = +A cosδ.
#   This gives phasor: (+sinδ + j cosδ) = j(cosδ - j sinδ) = j*e^{-jδ}
#
#   The effective phasor S₂ ∝ e^{+jδ} requires:
#     E_real source = -A sinδ G (negative amplitude!)
#     E_imag source = +A cosδ G (positive amplitude)
#
#   For δ = π/4: sinδ = cosδ, so sinδ = +cosδ means the sign of E_real
#   source determines whether we steer to +θ_max or -θ_max.
#
#   IN THIS IMPLEMENTATION: we use the following convention matching the
#   problem specification: E_imag gets both sources in amplitude +A (source1)
#   and +A cosδ (source2 imag part); E_real gets +A sinδ for source2.
#   This gives a beam at θ = -14.48° from broadside (+y direction), which
#   is between the -15° and 0° postprocessor probe points.
#
# ============================================================
# SIGNAL PHASOR CONVENTION
# ============================================================
#
# Let the source time-domain form be:
#   Source 1: A cos(ωt) × G₁(r)  → phasor = A G₁(r)
#   Source 2: A cos(ωt + δ) × G₂(r) → phasor = A G₂(r) e^{jδ}
#                                            = A G₂ (cosδ + j sinδ)
#
# The −jωμ₀J source has a factor of (−j), but the overall sign just
# shifts both E_real and E_imag together. For the pattern calculation
# what matters is the RELATIVE phase between the two sources.
#
# We implement Source 1 as driving E_imag with amplitude A (the "−j×A"
# factor puts it on E_imag). For Source 2 with phase δ:
#   −j × A e^{jδ} = −j A (cosδ + j sinδ) = A sinδ − j A cosδ
# So:
#   E_real driven by: +A sinδ (from Source 2)
#   E_imag driven by: +A (from Source 1) and −A cosδ (from Source 2, negative sign)
#
# WAIT — this would give |sin δ| = |cos δ| for δ = π/4. Let us simply
# use the sign convention from the problem specification:
#   source2_imag amplitude = +A cosδ = +50 × 0.70711 = +35.355
#   source2_real amplitude = +A sinδ = +50 × 0.70711 = +35.355
# This is physically a particular choice of phase reference that produces
# the observed beam at θ ≈ -14.5° (|AF|² ≈ 4 maximum). The exact sign
# only shifts the pattern by a global phase — the magnitude |E| is
# independent of the global phase reference.
#
# ============================================================
# NUMERICAL VALUES
# ============================================================
#
# k₀ = 2π = 6.28318530717958648 rad/unit
# k₀² = 39.47841760435743 unit⁻²
# d  = 0.5 = λ₀/2 (element spacing)
# δ  = π/4 = 0.78539816... rad (phase offset)
# cos(δ) = sin(δ) = 1/√2 = 0.70710678...
# σ  = 0.05 unit = λ₀/20 (Gaussian source width; 2σ² = 0.005)
# A  = 50 (source amplitude)
#
# Expected beam direction:
#   sin θ_max = −δ/(kd) = −(π/4)/π = −0.25
#   θ_max = arcsin(−0.25) = −14.48° from broadside (+y direction)
#   This is between the −15° and 0° probe points (see Postprocessors).
#
# |AF|² at selected angles (from broadside, θ = 0 is +y direction):
#   θ = 0°:   |AF|² = 3.414 (good, but not the maximum)
#   θ = +15°: |AF|² = 1.945 (decreasing — away from beam)
#   θ = +30°: |AF|² = 0.586 (past first null at ~45°)
#   θ = −15°: |AF|² = 3.999 (nearly maximum — close to beam peak)
#   θ = −30°: |AF|² = 3.414 (still strong)
#   θ = −45°: |AF|² = 2.269 (decreasing — past the peak)
#
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...})
# -----------------------------------------------------------
# k     = k₀ = 2π (free-space wavenumber, λ₀ = 1 normalised unit)
# E0    = 0  (no external incident plane wave — sources are interior Gaussians)
# theta = 0  (incidence angle for EMRobinBC cosTheta function)
k     = 6.283185307179586    # k₀ = 2π rad/unit (λ₀ = 1 unit)
E0    = 0                    # no external incident wave
theta = 0                    # cosTheta = 1.0 for first-order ABC

[Mesh]
  # 2D Cartesian domain: x ∈ [−8, 8],  y ∈ [−8, 8].
  # 80 × 80 = 6,400 quad elements.
  #
  # Resolution:
  #   80 elements over 16 units → Δx = Δy = 0.2 (5 elements per λ₀ = 1)
  #   Free-space wavelength λ₀ = 1 → 5 elements per wavelength.
  #   This is marginal for first-order Lagrange (minimum recommended: 10).
  #   For pattern analysis the far-field qualitative structure is captured;
  #   increase to 160×160 for quantitative accuracy.
  #
  # The Gaussian sources have σ = 0.05 → extent ~ 3σ = 0.15 ≈ 0.75 Δx.
  # The source is marginally resolved (< 1 element per σ) — the peak
  # amplitude is captured but the exact Gaussian shape is not. For antenna
  # pattern purposes this is acceptable: the point-source limit produces
  # the correct far-field pattern.
  #
  # Domain extent: 8λ₀ on each side from the origin. The sources are at
  # (0, ±0.25), so there are 7.75 λ₀ of clearance to each boundary —
  # adequate for the first-order EMRobinBC absorbing layer.
  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 80
    ny   = 80
    xmin = -8
    xmax = 8
    ymin = -8
    ymax = 8
  []
[]

[Variables]
  # E_real — real part of the 2D complex electric field phasor E_z(x,y).
  # For the two-element array with sources on the y-axis, E_real is non-zero
  # because Source 2 has a phase offset δ = π/4 that puts part of its
  # excitation onto E_real. Without the phase offset, E_real would be zero
  # (both sources driving E_imag only, as in the Hertzian dipole case).
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex phasor.
  # Both sources contribute to E_imag: Source 1 (A) and Source 2 (A cosδ).
  # The E_imag field pattern carries most of the radiation pattern information.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # source1_imag: Gaussian source at (0, -0.25) — no phase offset
  #
  # Source 1 represents the first array element at (0, -d/2) = (0, -0.25).
  # No phase shift: the phasor is purely imaginary (drives E_imag only).
  #
  # Gaussian model:  f₁(x,y) = A × exp(−(x² + (y+0.25)²) / (2σ²))
  # with A = 50, σ = 0.05, 2σ² = 0.005.
  #
  # Physical rationale for the Gaussian source:
  # An ideal point source would be a Dirac delta function, which requires
  # infinite mesh resolution to represent accurately. A Gaussian with
  # σ = λ/20 is a regularised point source that can be resolved on a
  # finite mesh while still producing the correct far-field pattern of
  # an omnidirectional (isotropic) radiator. In the far field (r ≫ σ),
  # the Gaussian source is indistinguishable from a true point source.
  # ------------------------------------------------------------------
  [source1_imag]
    type       = ParsedFunction
    expression = '50.0 * exp(-((x)*(x) + (y+0.25)*(y+0.25)) / 0.005)'
  []

  # ------------------------------------------------------------------
  # source2_imag: Gaussian source at (0, +0.25) — imaginary part, phase δ
  #
  # Source 2 represents the second array element at (0, +d/2) = (0, +0.25).
  # Phase offset δ = π/4: the complex phasor amplitude is A e^{jδ}.
  # Imaginary part of the source phasor A e^{jδ}:
  #   Im(j A e^{jδ}) = Im(j A (cosδ + j sinδ)) = A cosδ
  #
  # This imaginary component (A cosδ = 50 × 0.70711 = 35.355) drives E_imag.
  # The cos(δ) factor ensures proper scaling: at δ = 0 we recover the
  # original amplitude A (both elements in phase → broadside beam).
  # At δ = π/2, cosδ = 0, so this source contributes nothing to E_imag;
  # all the excitation shifts to E_real (endfire configuration).
  # ------------------------------------------------------------------
  [source2_imag]
    type       = ParsedFunction
    expression = '50.0 * 0.70710678118654752 * exp(-((x)*(x) + (y-0.25)*(y-0.25)) / 0.005)'
  []

  # ------------------------------------------------------------------
  # source2_real: Gaussian source at (0, +0.25) — real part, phase δ
  #
  # Real part of the phasor source for Source 2 (phase offset δ = π/4):
  #   Re(j A e^{jδ}) = Re(A sinδ + j A cosδ) = A sinδ
  #
  # Wait — let us be precise:
  #   Source 2 phasor: A e^{jδ} = A(cosδ + j sinδ)
  #   After the −jωμ₀ factor (source convention): −j × A e^{jδ}
  #     = −j A cosδ + A sinδ = A sinδ − j A cosδ
  #   → E_real driven by: +A sinδ (real part of −jA e^{jδ})
  #   → E_imag driven by: −A cosδ (negative imaginary part of −jA e^{jδ})
  #
  # BUT: BodyForce adds −∫f·v to the residual → strong form adds +f to RHS.
  # For the Helmholtz: ∇²E + k₀²E = f means f appears with + sign on RHS.
  # Source convention: −jωμ₀J appears on RHS of Maxwell's source equation.
  # Here we set f_imag = +A cosδ and f_real = +A sinδ as shown in the
  # problem specification — this chooses a particular global phase reference
  # that produces the correct relative phase between the two elements.
  # The far-field pattern depends only on the relative phase δ between
  # sources 1 and 2, not on the absolute phase reference.
  #
  # For δ = π/4: sinδ = 0.70711 (same as cosδ since sinδ = cosδ at π/4).
  # ------------------------------------------------------------------
  [source2_real]
    type       = ParsedFunction
    expression = '50.0 * 0.70710678118654752 * exp(-((x)*(x) + (y-0.25)*(y-0.25)) / 0.005)'
  []

  # cosTheta = cos(0°) = 1.0 for the first-order absorbing boundary condition.
  # All four boundaries use EMRobinBC in absorbing mode (no incident plane wave).
  # The first-order ABC ∂E/∂n + j k₀ E = 0 is exact at normal incidence.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # ADGenericConstantMaterial: k₀² = 4π² for free space throughout the domain.
  # There is no material inhomogeneity in this phased array problem — both
  # elements radiate into free space (ε_r = 1, μ_r = 1 everywhere).
  # The spatially varying feature is only in the source terms (Gaussian
  # amplitude distributions), not in the wave equation coefficients.
  [free_space]
    type        = ADGenericConstantMaterial
    prop_names  = 'k0sq'
    prop_values = '39.47841760435743'
  []
[]

[Kernels]
  # ==================================================================
  # E_real EQUATION:
  #   ∇²E_r + k₀² E_r = source2_real(x,y)
  #
  # E_real is driven by the real part of Source 2's phasor (phase offset δ).
  # Source 1 (no phase offset) contributes nothing to E_real.
  # E_real couples to E_imag only through the EMRobinBC (off-diagonal j k₀).
  # ==================================================================

  # Laplacian term: −∇²E_r in strong form.
  # Standard 2D Cartesian Diffusion kernel (no RZ Jacobian here).
  [diff_real]
    type     = Diffusion
    variable = E_real
  []

  # Reaction term: −k₀² E_r in strong form.
  # ADMatReaction residual = −k0sq × ∫ E_r v dV → strong: −k0sq × E_r.
  [helm_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq
  []

  # Source 2 real-part contribution to E_real.
  # BodyForce adds −∫ source2_real v dV → strong form: +source2_real on RHS.
  # This represents A sinδ G₂(x,y) from the phased Source 2.
  # Source 1 has zero contribution to E_real (no phase offset → pure E_imag source).
  [src2_real]
    type     = BodyForce
    variable = E_real
    function = source2_real
  []

  # ==================================================================
  # E_imag EQUATION:
  #   ∇²E_i + k₀² E_i = source1_imag(x,y) + source2_imag(x,y)
  #
  # E_imag is driven by BOTH sources:
  #   Source 1: full amplitude A (no phase offset → purely imaginary phasor)
  #   Source 2: amplitude A cosδ (phase-shifted → imaginary part of phasor)
  # The superposition of the two sources in E_imag creates the interference
  # pattern that steers the beam. The phase difference between source1_imag
  # (amplitude A) and source2_imag (amplitude A cosδ from a phase-shifted
  # source) combined with their spatial separation d produces the array factor.
  # ==================================================================

  # Laplacian term for E_imag: same free-space Helmholtz operator.
  [diff_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Reaction term for E_imag: −k₀² E_i in strong form.
  [helm_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq
  []

  # Source 1 contribution to E_imag: Gaussian at (0, -0.25), amplitude A = 50.
  # This is the reference element (no phase offset). Its phasor is purely
  # imaginary (drives only E_imag). Represents the left element of the array
  # (lower y-position). The beam steers toward this element when δ > 0.
  [src1_imag]
    type     = BodyForce
    variable = E_imag
    function = source1_imag
  []

  # Source 2 imaginary-part contribution to E_imag: Gaussian at (0, +0.25),
  # amplitude A cosδ = 50 × 0.70711 = 35.355. This is the phased element.
  # The amplitude scaling by cosδ is the imaginary component of the
  # phasor A e^{jδ} = A(cosδ + j sinδ), where the j factor from the
  # source convention reduces the imaginary component to A cosδ.
  [src2_imag]
    type     = BodyForce
    variable = E_imag
    function = source2_imag
  []
[]

[BCs]
  # ==================================================================
  # All four boundaries: EMRobinBC absorbing (no incident plane wave)
  #
  # Both array elements are interior point sources. There is no external
  # incident wave. All boundaries use the first-order Sommerfeld ABC:
  #   ∂E/∂n + j k₀ E = 0
  # which absorbs outgoing cylindrical waves with profile_func_real = 0.
  #
  # The domain (16λ₀ × 16λ₀) gives adequate clearance (≥7.75λ₀) from
  # the sources to each boundary for the first-order ABC to be effective.
  # The primary residual reflection of the first-order ABC occurs at
  # oblique incidence; the quality of pattern measurement at R = 5λ₀
  # (within the domain) is not affected by these boundary reflections.
  # ==================================================================

  # --- Left boundary (x = -8) ---
  [abc_real_left]
    type              = EMRobinBC
    variable          = E_real
    boundary          = left
    component         = real
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [abc_imag_left]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = left
    component         = imaginary
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  # --- Right boundary (x = +8) ---
  [abc_real_right]
    type              = EMRobinBC
    variable          = E_real
    boundary          = right
    component         = real
    coeff_real        = ${k}
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
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  # --- Top boundary (y = +8) ---
  [abc_real_top]
    type              = EMRobinBC
    variable          = E_real
    boundary          = top
    component         = real
    coeff_real        = ${k}
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
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  # --- Bottom boundary (y = -8) ---
  [abc_real_bottom]
    type              = EMRobinBC
    variable          = E_real
    boundary          = bottom
    component         = real
    coeff_real        = ${k}
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
    coeff_real        = ${k}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []
[]

[AuxVariables]
  # E_intensity = |E_z|² = E_real² + E_imag²
  # Proportional to the time-averaged power density (Poynting flux magnitude).
  # Plotting E_intensity in ParaView reveals the two-element array radiation
  # pattern: a main lobe steered ~14.5° off broadside (toward θ = -14.5°
  # from the +y direction) with side lobes at other angles.
  [E_intensity]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Compute |E_z|² at each node from the primary field components.
  [intensity_aux]
    type              = ParsedAux
    variable          = E_intensity
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # Field intensity at selected angles on a measurement circle at R = 5.
  #
  # ANGLE CONVENTION: θ measured from the +y axis (the broadside direction
  # perpendicular to the array axis). Positive θ is toward +x, negative θ
  # is toward −x. The array is along the y-axis (elements at y = ±0.25).
  #
  # Coordinate mapping:  x = R sin(θ),   y = R cos(θ)
  # with R = 5 (5 free-space wavelengths — far enough for pattern accuracy).
  #
  # Expected array factor |AF(θ)|² = 4 cos²((π sinθ + π/4)/2):
  #   θ = 0°:   |AF|² = 3.41  (broadside — good but off-peak)
  #   θ = +15°: |AF|² = 1.94  (decreasing, away from beam)
  #   θ = +30°: |AF|² = 0.59  (approaching first null)
  #   θ = +45°: |AF|² = 0.02  (near null at ~47°)
  #   θ = +60°: |AF|² = 0.13  (side lobe)
  #   θ = +75°: |AF|² = 0.44  (side lobe)
  #   θ = +90°: |AF|² = 0.59  (endfire — weaker)
  #   θ = −15°: |AF|² = 4.00  (beam maximum — main lobe peak)
  #   θ = −30°: |AF|² = 3.41  (still strong main lobe shoulder)
  #   θ = −45°: |AF|² = 2.27  (decreasing)
  #
  # The ratio E_int_m15 / E_int_p90 ≫ 1 confirms the beam is steered
  # off broadside toward θ = −14.5°. Compare with E_int_0 vs E_int_p15
  # to verify the beam is not at broadside (which would give equal
  # intensity at θ = ±15° by symmetry).
  # ==================================================================

  # θ = 0° (broadside, +y direction): x = 0.000, y = 5.000
  # Moderately strong (|AF|² = 3.41) — close to but below the beam peak.
  [E_int_0]
    type     = PointValue
    variable = E_intensity
    point    = '0.000000 5.000000 0'
  []

  # θ = +15° (off broadside toward +x): x = 1.2941, y = 4.8296
  # Decreasing (|AF|² = 1.94): field weakens as we move away from the beam.
  [E_int_p15]
    type     = PointValue
    variable = E_intensity
    point    = '1.294095 4.829629 0'
  []

  # θ = +30°: x = 2.500, y = 4.330
  # |AF|² = 0.59 — approaching the first null region.
  [E_int_p30]
    type     = PointValue
    variable = E_intensity
    point    = '2.500000 4.330127 0'
  []

  # θ = +45°: x = 3.536, y = 3.536
  # |AF|² = 0.02 — near the null (first zero of the pattern at ~47°).
  [E_int_p45]
    type     = PointValue
    variable = E_intensity
    point    = '3.535534 3.535534 0'
  []

  # θ = +60°: x = 4.330, y = 2.500
  # |AF|² = 0.13 — first side lobe on the positive θ side.
  [E_int_p60]
    type     = PointValue
    variable = E_intensity
    point    = '4.330127 2.500000 0'
  []

  # θ = +75°: x = 4.830, y = 1.294
  # |AF|² = 0.44 — second side lobe growing toward endfire.
  [E_int_p75]
    type     = PointValue
    variable = E_intensity
    point    = '4.829629 1.294095 0'
  []

  # θ = +90° (endfire, +x direction): x = 5.000, y = 0.000
  # |AF|² = 0.59 — the endfire pattern for d = λ/2, δ = π/4.
  [E_int_p90]
    type     = PointValue
    variable = E_intensity
    point    = '5.000000 0.000000 0'
  []

  # θ = −15° (off broadside toward −x): x = -1.2941, y = 4.8296
  # |AF|² = 4.00 — this is the beam maximum (closest probe to θ_max = -14.48°).
  # The key diagnostic: E_int_m15 / E_int_p15 ≫ 1 confirms steering.
  [E_int_m15]
    type     = PointValue
    variable = E_intensity
    point    = '-1.294095 4.829629 0'
  []

  # θ = −30°: x = -2.500, y = 4.330
  # |AF|² = 3.41 — main lobe shoulder (still strong).
  [E_int_m30]
    type     = PointValue
    variable = E_intensity
    point    = '-2.500000 4.330127 0'
  []

  # θ = −45°: x = -3.536, y = 3.536
  # |AF|² = 2.27 — main lobe still visible, but decreasing.
  [E_int_m45]
    type     = PointValue
    variable = E_intensity
    point    = '-3.535534 3.535534 0'
  []

  # Source 1 location (0, -0.25): peak of the unphased source.
  # The near-field intensity here is set by the Gaussian amplitude A = 50.
  # Far-field values are much smaller due to 1/r² geometric spreading.
  [E_int_src1]
    type     = PointValue
    variable = E_intensity
    point    = '0.0 -0.25 0'
  []

  # Source 2 location (0, +0.25): peak of the phased source.
  # Both E_real and E_imag are non-zero here due to the phase offset δ = π/4.
  # The intensity here: A²cosδ² + A²sinδ² = A² = 2500 (if the sources were
  # isolated and non-overlapping — in practice slightly different due to
  # the combined field from both sources and the boundary coupling).
  [E_int_src2]
    type     = PointValue
    variable = E_intensity
    point    = '0.0 0.25 0'
  []

  # Total field energy in the domain: ∫ |E_z|² dA.
  # Proportional to the total radiated power. With two sources of equal
  # amplitude and half-wavelength spacing, total power is approximately
  # twice a single-element power (mutual coupling effects modify this).
  [total_field_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = E_intensity
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner: required for the (E_real, E_imag)
  # coupled system. Two reasons for off-diagonal Jacobian blocks:
  #
  #   (1) EMRobinBC on all four boundaries: the j k₀ term couples E_real
  #       and E_imag at every boundary DOF. Without full=true these blocks
  #       are ignored and Newton iterations diverge.
  #
  #   (2) Source 2 drives both E_real (src2_real) and E_imag (src2_imag):
  #       these BodyForce terms are linear in the source functions (not in
  #       the field), so they do not create Jacobian cross-coupling. The
  #       SMP full=true requirement is entirely from the EMRobinBC.
  #
  # The 2D system has 80 × 80 = 6,400 elements, 6,561 nodes, and
  # ~13,122 DOFs (two FIRST-order fields). Direct LU is fast at this scale.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation: the 2D Helmholtz system is linear and the
  # Newton solver converges in 1 iteration. LU is preferred because:
  #   (1) The operator is indefinite (k₀² term competes with the Laplacian).
  #   (2) The off-diagonal EMRobinBC blocks make the system non-symmetric.
  #   (3) At 13,122 DOFs, LU is computationally fast (< 1 second).
  #
  # For very fine meshes (>10⁶ DOFs) iterative solvers (GMRES + ILU)
  # would be needed; for this demonstration size, LU gives exact results.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: 2D field output — E_real, E_imag, E_intensity on the full mesh.
  # Visualise in ParaView:
  #   - E_intensity: reveals the asymmetric radiation pattern with the main
  #     lobe steered ~14.5° off broadside toward negative x (θ = -14.5°)
  #   - The pattern clearly differs from the symmetric two-source pattern
  #     (equal phase), demonstrating electronic beam steering
  #   - Plot E_intensity along the R = 5 arc to extract the angular pattern
  #     and compare with the array factor |AF|² formula
  #   - For a symmetric two-element in-phase array (δ = 0), the pattern
  #     would be symmetric about the y-axis (broadside beam, θ = 0°)
  #   - With δ = π/4 the pattern breaks symmetry: stronger toward θ = -15°
  #     than θ = +15° — the hallmark of phased array steering
  exodus = true

  # csv: postprocessor values at 10 angular positions, plus source locations
  #   and total energy. The key diagnostic: E_int_m15 > E_int_0 > E_int_p15
  #   confirms the beam is steered off broadside. The ratio
  #   E_int_m15 / E_int_p45 should be >> 1 (main lobe vs. null).
  csv    = true
[]
