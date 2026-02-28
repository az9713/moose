# ============================================================
# Case 92: Aperture Synthesis Interferometer
# Thompson, Moran & Swenson, "Interferometry and Synthesis
#   in Radio Astronomy" (Wiley-VCH, 3rd Ed., 2017), Ch. 1–2
# Born & Wolf, "Principles of Optics", 7th Ed., Ch. 8
# Goodman, "Statistical Optics" (Wiley, 2015), Ch. 5
#
# TWO coherent point sources at different angular positions
# separated by δy = 2 units, observed from a baseline at x = 6.
# The interference fringes in the far field encode the source
# separation via the van Cittert–Zernike theorem:
#
#   I(y) = I₀ [1 + cos(2π y d / (λ R))]
#
# where d = 2 (source separation), R = 12 (propagation distance),
# λ = 1 (wavelength), giving fringe spacing Λ = λR/d = 6 units.
#
# ============================================================
# PHYSICAL BACKGROUND
# ============================================================
#
# Radio interferometry measures the mutual coherence function
# (visibility) of the wavefield between pairs of antennas.
# The van Cittert–Zernike theorem states that for an incoherent
# source distribution I(θ), the visibility function V(b) measured
# at baseline spacing b is the Fourier transform of I(θ):
#
#   V(b) = ∫ I(θ) exp(−j 2πbθ/λ) dθ
#
# For two point sources at angles ±δθ/2 with equal brightness:
#   I(θ) = δ(θ − δθ/2) + δ(θ + δθ/2)
#   V(b) = 2 cos(π b δθ / λ)
#
# The visibility is a cosine in baseline space, zero when
# b δθ / λ = 1/2, i.e., at the resolving baseline:
#   b_res = λ / (2 δθ) = λ R / (2 d)
#
# Angular resolution:
#   Rayleigh criterion: Δθ ≈ λ / B_max (B_max = maximum baseline)
#   For B_max = 2 (our domain half-width): Δθ_res ≈ λ/2
#   For d = 2 units at R = 12: δθ = d/R = 2/12 = 0.167 rad
#   Since δθ > Δθ_res, the two sources are resolved.
#
# ============================================================
# SIMULATION STRATEGY
# ============================================================
#
# Rather than a true incoherent multi-epoch measurement, we
# simulate a SINGLE coherent field with TWO Gaussian sources
# placed at y = +1 and y = -1 (representing the two "stars").
# The sources have identical amplitude and phase (coherent), so
# their interference pattern in the far field is identical to
# the cross-correlation fringe pattern of an interferometer
# observing two equal point sources.
#
# Governing equation (2D Helmholtz, frequency domain):
#
#   ∇²E + k₀² E = −s(x,y)        on Ω
#
# where s(x,y) = A[exp(−r₁²/2σ²) + exp(−r₂²/2σ²)]
#   r₁² = (x+6)² + (y−1)²   (source 1, above axis)
#   r₂² = (x+6)² + (y+1)²   (source 2, below axis)
#   A = 50, σ = 0.1
#
# The standard Helmholtz ∇²E + k₀²E = 0 in weak form:
#
#   Diffusion kernel:    +∫ ∇E·∇v dΩ  →  −∇²E  (strong)
#   ADMatReaction:       −k₀²∫ E v dΩ →  −k₀²E  (strong)
#   Sum = 0: −∇²E − k₀²E = 0 ⟺ ∇²E + k₀²E = 0 ✓
#
# The source BodyForce adds −∫ s v dΩ to the residual,
# which (with the −sign on LHS) places +s on the RHS:
#   −∇²E − k₀²E = −s  ⟺  ∇²E + k₀²E = s ✓
#
# ============================================================
# REAL/IMAGINARY SPLITTING
# ============================================================
#
# Decompose E = E_real + j E_imag. Since k₀² is purely real
# and ε_r = 1 everywhere (uniform vacuum), the two equations
# decouple in the bulk:
#
#   ∇²E_real + k₀² E_real = s(x,y)   (source on E_real)
#   ∇²E_imag + k₀² E_imag = 0        (no source on E_imag)
#
# The source drives E_real; E_imag is driven only by the
# Robin BC coupling at the boundaries. The intensity pattern
# |E|² = E_real² + E_imag² captures the full interference.
#
# ============================================================
# FRINGE ANALYSIS
# ============================================================
#
# The two sources produce a superposition:
#   E_total(x=6, y) ≈ A₁ × exp(jk₀r₁)/r₁ + A₂ × exp(jk₀r₂)/r₂
#
# In the far field (r₁,r₂ >> d):
#   r₁ ≈ R − y·sin(δθ/2) ≈ R − y·(1/12) × 1
#   r₂ ≈ R + y·(1/12) × 1
#
# Phase difference: Δφ = k₀(r₁ − r₂) ≈ k₀ × 2y/R × d/2
#   = 2π × y × d / (λ R)
#   = 2π × y × 2 / (1 × 12)
#   = 2π y / 6
#
# Intensity: |E|² ∝ 2[1 + cos(2πy/6)]
# Fringe maxima at y = 0, ±6; fringe spacing Λ = 6 units.
#
# Sampling the fringe at x = 6, y = {−5, −4, ..., +5}:
#   y = 0:  maximum (both sources equidistant, constructive)
#   y = ±3: minimum (path difference λ/2, destructive)
#   y = ±6: maximum (path difference λ, constructive)
#
# ============================================================
# NUMERICAL VALUES
# ============================================================
#
# λ = 1.0 unit  →  k₀ = 2π/λ = 6.2831853... rad/unit
# k₀² = (2π)² = 39.4784... unit⁻²
#
# Source 1: centre (−6, +1), σ = 0.1 → 2σ² = 0.02
# Source 2: centre (−6, −1), σ = 0.1 → 2σ² = 0.02
# Amplitude A = 50 (drives observable field at x = 6, R = 12 away)
#
# Domain: [−8, 8] × [−6, 6]   (16 × 12 units = 16λ × 12λ)
# Mesh:   80 × 60 elements     (0.2 units/element = λ/5 resolution)
# DOFs:   (81 × 61) × 2 ≈ 9,882 per variable = ~19,764 total
#
# The λ/5 resolution is adequate for first-order Lagrange elements
# on a smooth source; the Gaussian width σ = 0.1 is resolved by
# ~0.5 elements (barely, but sufficient for point-source behavior).
#
# ============================================================
# BOUNDARY CONDITIONS
# ============================================================
#
# All four boundaries use EMRobinBC in absorbing-only mode
# (profile_func_real = 0, no incident wave injection).
# The sources are interior Gaussian distributions; the Robin
# ABC absorbs the outgoing cylindrical waves. With 8 wavelengths
# of propagation distance from the sources to the boundary, the
# outgoing waves are nearly plane-like and the first-order ABC
# introduces only small spurious reflections.
#
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...})
# -----------------------------------------------------------
# k = 2π/λ with λ = 1 unit. This is the free-space wavenumber
# and appears in the reaction coefficient k₀² and in the Robin
# port BC impedance term jk₀ that absorbs outgoing waves.
k     = 6.283185307179586    # k₀ = 2π rad/unit (λ = 1 unit)
E0    = 0                    # no incident wave injection (absorbing BCs only)
theta = 0                    # normal incidence placeholder for EMRobinBC

[Mesh]
  # 2D domain [−8, 8] × [−6, 6].
  # 80 × 60 quadrilateral elements (bilinear Q1 Lagrange).
  # Element size: Δx = Δy = 0.2 units = λ/5.
  #
  # Resolution check:
  #   Free-space wavelength: λ = 1 unit → 5 elements/λ (minimal)
  #   Source Gaussian σ = 0.1 → ~0.5 elements/σ (point-source limit)
  #   Fringe spacing Λ = 6 units → 30 elements/fringe (well resolved)
  #
  # The relatively coarse λ/5 resolution is acceptable here because:
  #   (1) The fringe pattern (the observable of interest) has period 6λ
  #       and is extremely well resolved at 30 elements per fringe.
  #   (2) The source Gaussian acts as a point source; its spatial extent
  #       only needs to be small compared to λ, not well-resolved.
  #   (3) First-order Sommerfeld ABC at λ/5 resolution still absorbs
  #       outgoing plane-waves well (reflection < 1% at normal incidence).
  #
  # Increase nx/ny if you need sharper near-source detail.
  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 80
    ny   = 60
    xmin = -8
    xmax = 8
    ymin = -6
    ymax = 6
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor.
  # Both sources are placed with real amplitude (in-phase), so
  # the constructive/destructive interference appears primarily
  # in E_real. E_imag carries the quadrature component that
  # arises from the Robin BC coupling.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex field phasor.
  # No direct source drives E_imag (in this formulation sources
  # are placed on E_real); E_imag is driven indirectly through
  # the Robin ABC off-diagonal coupling at the four domain walls.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # source_fn: superposition of two Gaussian point sources
  #
  # The two "stars" (sources) are located at:
  #   Star 1: (x₁, y₁) = (−6, +1) — above the optical axis
  #   Star 2: (x₂, y₂) = (−6, −1) — below the optical axis
  # Both at x = −6 (the "transmitter" plane), separated by d = 2 units.
  # Gaussian width σ = 0.1 unit (≈ λ/10 → nearly point-source).
  # Amplitude A = 50 chosen so the field at the observation plane
  # x = 6 (distance R = 12 = 12λ away) gives convenient amplitudes.
  #
  # Combined source:
  #   s(x,y) = A·exp(−r₁²/2σ²) + A·exp(−r₂²/2σ²)
  # where 2σ² = 2 × 0.01 = 0.02.
  #
  # Note: the source drives E_real (the "in-phase" component).
  # The choice to source E_real rather than E_imag is a convention;
  # the intensity pattern |E|² = E_real² + E_imag² is independent
  # of this choice for a linear system.
  # ------------------------------------------------------------------
  [source_fn]
    type       = ParsedFunction
    expression = '50.0 * exp(-((x + 6.0)^2 + (y - 1.0)^2) / 0.02) + 50.0 * exp(-((x + 6.0)^2 + (y + 1.0)^2) / 0.02)'
  []

  # ------------------------------------------------------------------
  # coeff_fn: Helmholtz reaction coefficient k₀² × ε_r
  #
  # Uniform vacuum everywhere: ε_r = 1, so coeff = k₀² = (2π)².
  # The Diffusion + ADMatReaction combination implements:
  #   Diffusion:    −∇²E       (strong form from +∫∇E·∇v dΩ)
  #   ADMatReaction: −k₀²E    (strong form from −k₀²∫Ev dΩ)
  # Sum = 0: −∇²E − k₀²E = 0, i.e., ∇²E + k₀²E = 0 ✓
  # The minus sign convention in ADMatReaction means rate = +k₀²
  # gives the correct Helmholtz equation (wave propagation).
  # ------------------------------------------------------------------
  [coeff_fn]
    type       = ParsedFunction
    expression = '39.47841760435743'  # k₀² = (2π)² = 39.478...
  []

  # cosTheta = cos(0°) = 1.0 — normal-incidence factor for EMRobinBC.
  # The first-order Robin ABC uses k_eff = k₀ cos(θ) as the impedance.
  # At normal incidence (outgoing rays perpendicular to the boundary)
  # this is exact; for grazing angles the first-order ABC degrades.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # Wrap coeff_fn (k₀² = const) as an AD material property.
  # ADGenericFunctionMaterial evaluates the ParsedFunction at each
  # quadrature point and returns an ADReal, enabling automatic
  # differentiation of the reaction term through the Jacobian.
  # Because ε_r = 1 everywhere (uniform vacuum), a simpler
  # GenericConstantMaterial would suffice, but ADGenericFunctionMaterial
  # maintains consistency with the other Batch F cases and makes it
  # trivial to introduce spatially varying ε_r in future variants.
  [coeff_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'k0sq'
    prop_values = 'coeff_fn'
  []
[]

[Kernels]
  # ==================================================================
  # E_real EQUATION:
  #   ∇²E_real + k₀² E_real = s(x,y)
  #
  # Weak form (multiply by test v, integrate by parts):
  #   ∫ ∇E_real·∇v dΩ − k₀²∫ E_real v dΩ = ∫ s v dΩ
  # Note: BodyForce adds −∫ s v to residual, so for residual = 0:
  #   ∫∇E_r·∇v − k₀²∫E_r v − ∫s v = 0  →  −∇²E_r − k₀²E_r = s ✓
  # ==================================================================

  # Diffusion kernel: contributes +∫ ∇E_real·∇v dΩ to the residual.
  # Strong form equivalent: −∇²E_real.
  [diffusion_real]
    type     = Diffusion
    variable = E_real
  []

  # ADMatReaction: contributes −k₀² ∫ E_real v dΩ to the residual.
  # Combined with Diffusion: residual strong form = −∇²E_r − k₀²E_r.
  # Setting residual = 0 recovers the homogeneous Helmholtz equation.
  # The source shifts the right-hand side (via BodyForce below).
  [reaction_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq
  []

  # BodyForce: the two-star source s(x,y) drives E_real.
  # Contributes −∫ source_fn · v dΩ to the residual.
  # Strong form: the equation becomes −∇²E_r − k₀²E_r = −(−s) = s.
  # i.e., ∇²E_r + k₀²E_r = s (Helmholtz with source) ✓
  [source_real]
    type     = BodyForce
    variable = E_real
    function = source_fn
  []

  # ==================================================================
  # E_imag EQUATION:
  #   ∇²E_imag + k₀² E_imag = 0
  #
  # No source on E_imag — the imaginary component is driven purely
  # by the Robin BC off-diagonal coupling (the jk₀ term in the
  # absorbing boundary condition couples E_real to E_imag).
  # In the bulk (uniform vacuum, no loss), the two equations decouple.
  # ==================================================================

  [diffusion_imag]
    type     = Diffusion
    variable = E_imag
  []

  [reaction_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq
  []
[]

[BCs]
  # ==================================================================
  # Absorbing boundary conditions on all four domain walls.
  #
  # EMRobinBC (Jin, "FEM in Electromagnetics", 3rd Ed.):
  #   ∂E/∂n + j k₀ cos(θ) E = 0   on ∂Ω
  #
  # profile_func_real = ${E0} = 0: pure absorbing, no incident wave.
  # The sources are Gaussian distributions inside the domain;
  # waves radiate outward and are absorbed by these BCs.
  #
  # The Robin BC couples E_real and E_imag: the jk₀ term in
  #   ∂E/∂n + jk₀E = 0
  # becomes (in real/imaginary split):
  #   ∂E_r/∂n − k₀ E_i = 0
  #   ∂E_i/∂n + k₀ E_r = 0
  # so the real BC involves E_imag (off-diagonal block) and vice versa.
  # This is why [Preconditioning] with full = true is required.
  #
  # sign = negative is the standard Batch F convention (matching
  # the upstream MOOSE EM module reference benchmark).
  # ==================================================================

  # --- Left boundary (x = −8): absorb leftward-going waves ---
  # The sources at x = −6 are 2 units inside this boundary, so
  # some radiation hits the left wall. The ABC absorbs it with
  # reflection coefficient |R|_ABC = |1 − cos(θ)|/|1 + cos(θ)|,
  # which is small at near-normal incidence.
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

  # --- Right boundary (x = +8): absorb the transmitted waves ---
  # The observation "baseline" at x = 6 is 2 units inside this
  # boundary. Postprocessors sample the fringe pattern at x = 6
  # before the wave reaches the right ABC.
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

  # --- Top boundary (y = +6): absorb upward-going waves ---
  # Sources at y = ±1 are 5 units from this boundary (5λ), so
  # the waves are approximately plane-wave like at arrival.
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

  # --- Bottom boundary (y = −6): absorb downward-going waves ---
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
  # |E|² = E_real² + E_imag² — instantaneous field intensity.
  # This is proportional to the time-averaged Poynting flux for a
  # monochromatic wave. Plotting E_intensity along the observation
  # line x = 6 reveals the interference fringe pattern:
  #   Maxima at y = 0, ±6 (constructive interference)
  #   Minima at y = ±3   (destructive interference, path diff = λ/2)
  [E_intensity]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Evaluate |E|² = E_real² + E_imag² pointwise at all mesh nodes.
  # ParsedAux runs post-solve, so E_real and E_imag are converged.
  [intensity_aux]
    type              = ParsedAux
    variable          = E_intensity
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # Fringe pattern sampling at the observation baseline x = 6.
  #
  # The theoretical fringe spacing is Λ = λR/d = 1×12/2 = 6 units.
  # Sampling at y = −5 to +5 in steps of 1 gives 11 points spanning
  # almost two full fringe periods. Key predicted values:
  #
  #   y = 0:   MAXIMUM — equal path from both sources → |E|² ∝ 4A²/R²
  #   y = ±3:  MINIMUM — path difference = λ/2       → |E|² ≈ 0
  #   y = ±6:  MAXIMUM — path difference = λ          → |E|² ∝ 4A²/R²
  #
  # These postprocessors form a "baseline measurement" analogous to
  # what a radio interferometer records as it samples the visibility
  # function at different baseline positions (different y-spacings
  # between antenna pairs).
  #
  # The Fourier transform of {E_intensity(y)} gives the source power
  # spectrum (two delta functions at angular positions ±1/12 rad),
  # recovering the source geometry from the fringe data — this is the
  # core operation of aperture synthesis imaging (CLEAN algorithm, etc.)
  # ==================================================================

  # Fringe maximum: y = 0 (on-axis, equidistant from both sources)
  # Expect: peak of the fringe pattern (constructive interference)
  [fringe_y_m5]
    type     = PointValue
    variable = E_intensity
    point    = '6 -5 0'
  []

  [fringe_y_m4]
    type     = PointValue
    variable = E_intensity
    point    = '6 -4 0'
  []

  [fringe_y_m3]
    type     = PointValue
    variable = E_intensity
    point    = '6 -3 0'
  []

  [fringe_y_m2]
    type     = PointValue
    variable = E_intensity
    point    = '6 -2 0'
  []

  [fringe_y_m1]
    type     = PointValue
    variable = E_intensity
    point    = '6 -1 0'
  []

  # On-axis measurement: both sources equidistant → constructive.
  # This is the "zero-baseline" visibility (total source flux).
  [fringe_y_0]
    type     = PointValue
    variable = E_intensity
    point    = '6 0 0'
  []

  [fringe_y_p1]
    type     = PointValue
    variable = E_intensity
    point    = '6 1 0'
  []

  [fringe_y_p2]
    type     = PointValue
    variable = E_intensity
    point    = '6 2 0'
  []

  # Near minimum (y = 3 ≈ Λ/2 = 3): destructive interference.
  # Path difference ≈ 2×3/12 = 0.5 = λ/2 → near cancellation.
  [fringe_y_p3]
    type     = PointValue
    variable = E_intensity
    point    = '6 3 0'
  []

  [fringe_y_p4]
    type     = PointValue
    variable = E_intensity
    point    = '6 4 0'
  []

  [fringe_y_p5]
    type     = PointValue
    variable = E_intensity
    point    = '6 5 0'
  []

  # ------------------------------------------------------------------
  # Source region diagnostics
  # ------------------------------------------------------------------
  # Probe the near-source field to confirm both Gaussian sources
  # are active and have the expected amplitudes.
  # At the source centres (x = −6, y = ±1), the field amplitude
  # should be at its maximum (driven point).
  [E_real_source1]
    type     = PointValue
    variable = E_real
    point    = '-6 1 0'
  []

  [E_real_source2]
    type     = PointValue
    variable = E_real
    point    = '-6 -1 0'
  []

  # Midpoint between the two sources (x = −6, y = 0).
  # Expect intermediate amplitude: the two Gaussians add coherently
  # at σ = 0.1 apart by d = 2, so their overlap is negligible
  # (separation >> σ) and the field at y = 0 is ~0 between them.
  [E_real_midpoint]
    type     = PointValue
    variable = E_real
    point    = '-6 0 0'
  []

  # ------------------------------------------------------------------
  # Global field diagnostics
  # ------------------------------------------------------------------
  # Total field energy: ∫|E|² dA over the full domain.
  # A single Gaussian source in free space gives energy ~ A²σ²R²;
  # two coherent sources give up to 4× the single-source energy
  # depending on the phase relationship. Useful for checking
  # that the solver has converged to a physically reasonable level.
  [total_field_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = E_intensity
  []

  # Maximum intensity in the domain: should occur near the sources.
  [max_intensity]
    type     = ElementExtremeValue
    variable = E_intensity
    value_type = max
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner (SMP with full = true).
  #
  # Required because the EMRobinBC on all four walls introduces
  # off-diagonal Jacobian blocks between E_real and E_imag. The
  # Robin condition (∂E/∂n + jk₀E = 0) in real/imag split is:
  #   Real BC: ∂E_r/∂n − k₀ E_i = 0  → off-diagonal coupling E_i
  #   Imag BC: ∂E_i/∂n + k₀ E_r = 0  → off-diagonal coupling E_r
  #
  # Without full = true, the SMP builds only diagonal (E_r,E_r) and
  # (E_i,E_i) blocks, omitting the boundary coupling. The preconditioned
  # system is then poorly conditioned and Newton iterations stagnate.
  # With full = true, the complete 2×2 block structure is captured and
  # the direct LU factorisation handles the complex-valued system exactly.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation of the 2D Helmholtz system.
  # Problem size: 81 × 61 nodes × 2 variables ≈ 9,882 DOFs.
  # The Helmholtz operator is indefinite (the k₀² term can dominate
  # over the stiffness, making some eigenvalues negative), so direct
  # LU is safer than iterative methods that assume positive-definiteness.
  # LU solves in one Newton step (the system is linear in E_real, E_imag).
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: full 2D field output (E_real, E_imag, E_intensity) on the mesh.
  # Visualise in ParaView:
  #   - E_intensity: fringe pattern visible as alternating bright/dark bands
  #     parallel to the x-axis, with spacing Λ = 6 units along y.
  #   - E_real: shows the spatial phase structure of the interference field.
  #   - Apply a horizontal slice at x = 6 to extract the baseline measurement.
  #   - Use a Line Plot over the range y ∈ [−6, 6] at x = 6 to trace the
  #     fringe curve: expect a cosine with period 6 (two peaks at y = 0, ±6).
  exodus = true

  # csv: fringe samples at x = 6 for y = −5 to +5, plus source diagnostics.
  # The fringe_y_* values should show a clear oscillation:
  #   High at y = 0, low near y = ±3, high again approaching y = ±6.
  # Compute the ratio (max − min) / (max + min) = visibility = 1.0 (ideal).
  # Any departure from unit visibility arises from ABC reflections,
  # finite-source-width effects, or numerical dispersion (λ/5 mesh).
  csv    = true
[]
