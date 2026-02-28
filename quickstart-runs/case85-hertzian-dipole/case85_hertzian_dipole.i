# ============================================================
# Case 85: Hertzian Dipole Radiation Pattern
# 2D Axisymmetric Frequency-Domain Helmholtz — RZ Coordinates
#
# References:
#   Kong, J.A., Electromagnetic Wave Theory, 2nd ed. (EMW, 1990),
#     Ch. 2 (radiation from electric dipoles)
#   Pozar, D.M., Microwave Engineering, 4th ed. (Wiley, 2012),
#     Sec. 9.2–9.3 (infinitesimal dipole radiation pattern)
#   Balanis, C.A., Antenna Theory, 3rd ed. (Wiley, 2005), Ch. 4
#   MIT 6.635 Advanced Electromagnetism, Spring 2003, Prof. Jin Au Kong
#   OCW: https://ocw.mit.edu/courses/6-635-advanced-electromagnetism-spring-2003/
#
# ============================================================
# PHYSICAL BACKGROUND: THE HERTZIAN DIPOLE
# ============================================================
#
# A Hertzian (infinitesimal) dipole is a current element of
# infinitesimal length Δl oriented along the z-axis, carrying
# a current I·exp(+jωt).  It is the fundamental building block
# of antenna theory: any real antenna can be decomposed into a
# superposition of Hertzian dipoles via the Fourier integral.
#
# The dipole radiates electromagnetic waves with a characteristic
# doughnut-shaped radiation pattern:
#
#   U(θ) ∝ sin²(θ)
#
# where θ is the polar angle measured from the dipole axis (z-axis):
#   - Zero radiation along the dipole axis (θ = 0°, 180°)
#   - Maximum radiation in the equatorial plane (θ = 90°)
#
# This sin²θ pattern arises from the azimuthal component of the
# electric field:
#   E_θ ∝ sin(θ)   in the far field (kr >> 1)
# so the power density |E|² ∝ sin²(θ).
#
# ============================================================
# COORDINATE SYSTEM — RZ AXISYMMETRIC FORMULATION
# ============================================================
#
# In MOOSE's RZ coordinate system (coord_type = RZ):
#   - The horizontal axis (x in MOOSE) represents the RADIAL direction r
#   - The vertical axis   (y in MOOSE) represents the AXIAL direction z
#   - Azimuthal symmetry: ∂/∂φ = 0 for all fields
#   - The domain is the half-plane r ≥ 0 (with r = 0 on the left boundary)
#
# The dipole is oriented along the z-axis (y in MOOSE) and located
# at the origin (r=0, z=0).  By azimuthal symmetry, the only
# non-zero field component is E_φ (azimuthal electric field), which
# in RZ coordinates is a scalar function E_φ(r, z).
#
# IMPORTANT: The polar angle θ in 3D spherical coordinates is related
# to the RZ coordinates by:
#   r = R · sin(θ),   z = R · cos(θ)   (where R = sqrt(r² + z²))
# So:
#   θ = 0°:  the z-axis (r = 0, z > 0)     → zero radiation
#   θ = 90°: the equatorial plane (z = 0)   → maximum radiation
#
# ============================================================
# GOVERNING EQUATION IN RZ COORDINATES
# ============================================================
#
# For TE_φ polarisation (E = E_φ(r,z) · φ̂) in cylindrical coordinates,
# the Helmholtz equation for E_φ is:
#
#   (1/r) ∂/∂r (r ∂E_φ/∂r) + ∂²E_φ/∂z² − E_φ/r² + k₀² E_φ = −J_source
#
# The −E_φ/r² term (from the vector Laplacian in cylindrical coordinates)
# distinguishes the azimuthal field from the scalar (TE_z or TM_z) case.
#
# HOWEVER: For a z-directed Hertzian dipole source, the appropriate
# field to solve is the z-component of the vector potential A_z, which
# satisfies the SCALAR Helmholtz equation:
#
#   ∇²A_z + k₀² A_z = −μ₀ J_z(r, z)
#
# where ∇² in cylindrical coordinates with ∂/∂φ = 0 is:
#
#   ∇²A_z = (1/r)∂/∂r(r ∂A_z/∂r) + ∂²A_z/∂z²
#
# In MOOSE with coord_type = RZ, the standard MOOSE Diffusion kernel
# automatically handles the (1/r)∂/∂r(r·...) term — it computes the
# cylindrical-coordinate Laplacian correctly.  This is because MOOSE
# modifies the integration measure dV = r dr dφ dz and the gradient
# operator to account for the RZ coordinate system.
#
# So the MOOSE input solves:
#
#   ∇²_cyl A_z + k₀² A_z = −source(r,z)
#
# where ∇²_cyl is automatically the cylindrical Laplacian.
# The fields are then related to A_z by:
#   E_z ∝ ∂²A_z/∂z² + k₀² A_z   (far field: E_z ∝ sin(θ)/R × exp(-jk₀R))
#   E_r ∝ ∂²A_z/∂r∂z             (far field: E_r ∝ cos(θ)/R × exp(-jk₀R))
#
# The radiation pattern is dominated by E_θ in the far field, and
# |E_θ|² ∝ sin²(θ), giving the characteristic dipole pattern.
#
# We label the field variable E_real + j E_imag for continuity with
# the rest of the Batch F cases, noting that it physically represents
# the vector potential A_z (which has the same sinusoidal pattern as E_θ).
#
# ============================================================
# DOMAIN AND SOURCE
# ============================================================
#
# λ₀ = 1.0 (normalised), k₀ = 2π/λ₀ = 2π rad/unit
# k₀² = (2π)² ≈ 39.478 rad²/unit²
#
# Domain:
#   r ∈ [0, R_max] = [0, 6.0]   (6 wavelengths in radial direction)
#   z ∈ [−Z_max, Z_max] = [−6, 6]  (±6 wavelengths in axial direction)
#
# Mesh: 60 elements in r × 120 elements in z = 7200 QUAD4 elements.
# Resolution: 60/6 = 10 elements per wavelength in r,
#             120/12 = 10 elements per wavelength in z.
# This is the minimum for the first-order ABC to work reasonably well;
# increasing to 20 el/λ would improve absorbing BC performance.
#
# Source: Gaussian approximation to a point dipole at (r=0, z=0):
#   J(r,z) = A · exp(−(r² + z²) / (2σ²))
#   A = 100,  σ = 0.05 (≈ λ₀/20 — narrow compared to wavelength)
#
# The source is applied as a BodyForce on E_imag only.  This follows
# the convention that the source −jωμ₀J (with the +jωt convention)
# is purely imaginary; its effect appears in the imaginary field equation.
# (Same convention as Cases 77, 79, 83.)
#
# ============================================================
# BOUNDARY CONDITIONS
# ============================================================
#
# Left boundary (r = 0, the symmetry axis):
#   NATURAL NEUMANN BC — no explicit BC block needed.
#   For the scalar Helmholtz (vector potential A_z), the natural
#   condition ∂A_z/∂r = 0 at r = 0 is physically correct: by symmetry,
#   A_z is an even function of r (A_z(−r,z) = A_z(+r,z)), so ∂A_z/∂r = 0
#   at r = 0.  MOOSE's RZ formulation handles this automatically.
#
# Right boundary (r = R_max = 6):
# Top boundary   (z = +Z_max = +6):
# Bottom boundary (z = −Z_max = −6):
#   EMRobinBC — first-order absorbing (Sommerfeld radiation condition).
#   ∂E/∂n + j k₀ E = 0  absorbs outgoing cylindrical waves.
#   With profile_func_real = 0 (no incident wave injection), these
#   are pure absorbing boundaries.
#
# NOTE: The first-order ABC is approximate; at 6λ from the source,
# the outgoing waves are nearly plane-wave-like and the ABC is
# reasonably effective.  A perfectly matched layer (PML) would give
# much better absorption, but is beyond the scope of this case.
#
# ============================================================
# RADIATION PATTERN EXTRACTION
# ============================================================
#
# The radiation pattern U(θ) ∝ |E(R,θ)|² is sampled at a circle of
# radius R = 4λ₀ = 4.0 units.  Points on this circle at angles θ
# from the z-axis have coordinates:
#
#   r = R · sin(θ),   z = R · cos(θ)
#
# Angles probed:
#   θ = 0°:   r = 0.000, z = 4.000  → expected |E|² ≈ 0 (on z-axis)
#   θ = 30°:  r = 2.000, z = 3.464  → expected |E|² ∝ sin²(30°) = 0.25
#   θ = 45°:  r = 2.828, z = 2.828  → expected |E|² ∝ sin²(45°) = 0.50
#   θ = 60°:  r = 3.464, z = 2.000  → expected |E|² ∝ sin²(60°) = 0.75
#   θ = 90°:  r = 4.000, z = 0.000  → expected |E|² ∝ sin²(90°) = 1.00 (max)
#
# The ratios E_intensity(θ) / E_intensity(90°) should reproduce the
# sin²(θ) pattern when the domain is electrically large enough.
# With 6λ domain and a first-order ABC, expect some deviation from
# the ideal pattern (±5–10%) due to ABC reflections and near-field effects.
#
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...})
# -----------------------------------------------------------
# k:    free-space wavenumber k₀ = 2π/λ₀ = 2π (normalised: λ₀ = 1)
# E0:   incident field amplitude = 0 (pure absorbing BC, no incident wave)
# theta: incidence angle = 0 (normal incidence placeholder for EMRobinBC)
# R_max: radial domain extent = 6 λ₀
# Z_max: axial half-extent = 6 λ₀ (domain z ∈ [−Z_max, +Z_max])
k     = 6.283185307179586  # k₀ = 2π rad/unit (λ₀ = 1.0 unit)
E0    = 0                  # no incident plane wave (source is interior Gaussian)
theta = 0                  # placeholder for EMRobinBC (normal incidence)
R_max = 6.0                # radial domain extent [λ₀]
Z_max = 6.0                # axial half-domain extent [λ₀]

[Mesh]
  # 2D RZ mesh: x-axis = r (radial), y-axis = z (axial).
  # Domain: r ∈ [0, 6] × z ∈ [−6, 6].
  #
  # GeneratedMeshGenerator in 2D names the four boundary sides automatically:
  #   'left'   → x = 0 (r = 0):    symmetry axis (natural Neumann, no explicit BC)
  #   'right'  → x = R_max (r=6):  outer radial absorbing boundary (EMRobinBC)
  #   'bottom' → y = −Z_max:       lower axial absorbing boundary (EMRobinBC)
  #   'top'    → y = +Z_max:       upper axial absorbing boundary (EMRobinBC)
  #
  # Mesh resolution:
  #   Radial:   nx = 60 elements over 6 λ₀ → 10 elements/λ₀
  #   Axial:    ny = 120 elements over 12 λ₀ → 10 elements/λ₀
  #   Total:    60 × 120 = 7200 QUAD4 elements (14,762 nodes for LAGRANGE/FIRST)
  #
  # The Gaussian source has σ = 0.05 λ₀; with dx = dz ≈ 0.1, there are
  # approximately 2 elements across the source radius — marginal but
  # sufficient to capture the source amplitude and drive the far-field pattern.
  # A finer mesh (e.g., 120×240) would give a slightly sharper near-field
  # but the far-field radiation pattern is insensitive to source resolution.
  # coord_type = RZ tells MOOSE to use cylindrical coordinates.
  # This modifies:
  #   (1) The volume integration measure: dV = r × dx × dy (instead of dx × dy)
  #   (2) The divergence operator in Diffusion: computes (1/r)∂/∂r(r ∂/∂r) + ∂²/∂z²
  #       automatically — the cylindrical Laplacian without any special kernels.
  #   x → r (radial, r ≥ 0),  y → z (axial, any sign)
  coord_type = RZ

  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 60    # radial divisions (x-direction in MOOSE = r in physics)
    ny   = 120   # axial divisions  (y-direction in MOOSE = z in physics)
    xmin = 0
    xmax = ${R_max}
    ymin = ${fparse -Z_max}
    ymax = ${Z_max}
  []
[]

[Variables]
  # E_real — real part of the complex scalar field A_z = E_real + j E_imag.
  # In the RZ axisymmetric context this represents the real part of the
  # z-component of the magnetic vector potential (or equivalently, the
  # azimuthal electric field for TE_φ modes).
  # FIRST LAGRANGE basis: C⁰ continuous, appropriate for the scalar Helmholtz.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex field.
  # The BodyForce source acts directly on E_imag (imaginary source convention
  # from the +jωt time dependence of the harmonic source).
  # In the far field: E_real and E_imag encode the cos(k₀R)/R and sin(k₀R)/R
  # components of the cylindrical outgoing wave respectively.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # source_fn: Gaussian approximation to a z-directed point dipole at origin.
  #
  # Physical interpretation:
  #   J_z(r, z) ≈ I·Δl·δ(r)·δ(z)/(2πr)   [ideal point dipole in RZ]
  #
  # The δ-function singularity is regularised by a Gaussian of width σ:
  #   J_z(r, z) ≈ A × exp(−(r² + z²) / (2σ²))
  # with A = 100, σ = 0.05 λ₀.
  #
  # The source is localised within ~3σ = 0.15 λ₀ of the origin — much
  # smaller than λ₀ = 1.0, so it acts as a true Hertzian dipole for
  # field points in the far field (kr >> 1, i.e., r >> 0.16 units).
  #
  # Source magnitude: 2σ² = 2 × 0.05² = 0.005
  # At r = z = 0: J = 100 × exp(0) = 100 (peak amplitude)
  # At r = σ:     J = 100 × exp(−0.5) ≈ 60.65 (−4.3 dB point)
  # At r = 2σ:    J = 100 × exp(−2.0) ≈ 13.53 (−17.4 dB point)
  # Effectively zero for r > 0.2.
  # ------------------------------------------------------------------
  [source_fn]
    type       = ParsedFunction
    expression = '100.0 * exp(-(x * x + y * y) / 0.005)'
  []

  # ------------------------------------------------------------------
  # cosTheta: placeholder for EMRobinBC incidence angle factor.
  # All boundaries are absorbing only (profile_func_real = ${E0} = 0),
  # so the exact value of cosTheta merely sets the ABC wavenumber:
  #   k_eff = k₀ × cosTheta = k₀ × cos(0) = k₀  (normal incidence)
  # For a curved wave front from a point source, the true incidence
  # angle varies across the boundary; normal incidence (θ = 0) is the
  # best first-order approximation for a distant source.
  # ------------------------------------------------------------------
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # ------------------------------------------------------------------
  # k₀² material property for the ADMatReaction Helmholtz term.
  #
  # The Helmholtz equation in cylindrical RZ coordinates:
  #   ∇²_cyl E + k₀² E = −source
  # requires the spatially constant coefficient k₀² = (2π)² ≈ 39.478.
  #
  # ADGenericConstantMaterial provides an AD (automatic differentiation)
  # material property that ADMatReaction can consume.  Using AD is
  # important here because the Newton Jacobian of ADMatReaction depends
  # on the reaction rate; with AD this is computed exactly.
  #
  # NOTE: The entire domain has k₀² constant (uniform free space); there
  # is no spatially varying permittivity in this case (the source is a
  # localised current, not a dielectric contrast).
  # ------------------------------------------------------------------
  [k0sq_material]
    type        = ADGenericConstantMaterial
    prop_names  = 'k0sq'
    prop_values = '${fparse k*k}'
  []
[]

[Kernels]
  # ==================================================================
  # E_real EQUATION:
  #   ∇²_cyl E_real + k₀² E_real = 0
  #
  # (No source term — the source is purely imaginary in phasor form.)
  # Weak form: ∫ ∇E_real · ∇φ dV − ∫ k₀² E_real φ dV = 0
  # (with the cylindrical volume element dV = r dr dz, handled by MOOSE)
  # ==================================================================

  # Cylindrical Laplacian for E_real.
  # In RZ mode, MOOSE Diffusion computes:
  #   ∫ (∂E/∂r ∂φ/∂r + ∂E/∂z ∂φ/∂z) r dr dz
  # which corresponds to the strong form −(1/r)∂/∂r(r ∂E/∂r) − ∂²E/∂z²
  # — the cylindrical vector Laplacian for a scalar field.
  [diff_real]
    type     = Diffusion
    variable = E_real
  []

  # Helmholtz k₀² term for E_real.
  # ADMatReaction residual = −rate × ∫ E_real φ dV = −k₀² ∫ E_real φ r dr dz
  # Strong form equivalent: −k₀² E_real.
  # Combined with Diffusion: −∇²E_real − k₀² E_real = 0
  # ↔ ∇²E_real + k₀² E_real = 0 ✓
  [helm_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq
  []

  # ==================================================================
  # E_imag EQUATION:
  #   ∇²_cyl E_imag + k₀² E_imag = −source(r, z)
  #
  # The Gaussian source drives the E_imag equation (imaginary source
  # from the +jωt convention).  E_real is driven indirectly through
  # the absorbing BCs which couple E_real ↔ E_imag at the boundaries.
  # ==================================================================

  # Cylindrical Laplacian for E_imag.
  [diff_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Helmholtz k₀² term for E_imag.
  # Same structure as for E_real; both equations share the same k₀².
  [helm_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq
  []

  # Gaussian point source.
  # BodyForce residual: −∫ source_fn × φ dV (with RZ volume element r dr dz)
  # Strong form: E_imag equation RHS = +source_fn
  # This drives the imaginary part of the azimuthal vector potential,
  # which in the far field produces both E_imag (cos-like phase) and
  # E_real (sin-like phase) components due to the ABC coupling at boundaries.
  [source_imag]
    type     = BodyForce
    variable = E_imag
    function = source_fn
  []
[]

[BCs]
  # ==================================================================
  # Left boundary (r = 0): SYMMETRY AXIS — NO EXPLICIT BC
  # ==================================================================
  # The natural Neumann BC ∂E/∂r = 0 at r = 0 is physically correct
  # for the scalar field A_z(r,z).  By azimuthal symmetry A_z must be
  # an even function of r: A_z(−r,z) = A_z(+r,z), implying ∂A_z/∂r|_{r=0} = 0.
  # MOOSE applies this automatically (no [BCs] block entry needed for r=0).
  #
  # Do NOT apply a DirichletBC at r = 0 — that would force E = 0 on the
  # symmetry axis, which is only correct for higher-order azimuthal modes
  # (m ≥ 1), not for the m = 0 (z-directed dipole) case solved here.
  # ==================================================================

  # ==================================================================
  # Right boundary (r = R_max = 6): Absorbing Robin BC
  # ==================================================================
  # EMRobinBC absorbs outgoing cylindrical waves from the source/scatterer.
  # The boundary condition in strong form (on the total field E = E_r + jE_i):
  #   ∂E/∂n + j k₀ E = 0
  # where n is the outward normal (pointing in the +r direction here).
  # With profile_func_real = 0 (E0 = 0): pure absorption, no injection.
  #
  # Sign = negative: consistent with the Diffusion kernel convention used in
  # Cases 74–83.  EMRobinBC adds to the residual:
  #   R_real += +k₀ × ∫ E_imag φ dS  (from j k₀ E = j k₀(E_r+jE_i) → real part: −k₀ E_i)
  #   R_imag += −k₀ × ∫ E_real φ dS  (imaginary part: +k₀ E_r)
  # The sign=negative convention reverses these to match the Diffusion sign.
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

  # ==================================================================
  # Top boundary (z = +Z_max = +6): Absorbing Robin BC
  # Outward normal points in the +z direction on this face.
  # Same absorbing structure as the right boundary.
  # ==================================================================
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

  # ==================================================================
  # Bottom boundary (z = −Z_max = −6): Absorbing Robin BC
  # Outward normal points in the −z direction on this face.
  # Same absorbing structure as the top boundary.
  # ==================================================================
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
  # E_intensity = |E|² = E_real² + E_imag²
  # Proportional to the cycle-averaged radiated power density (Poynting magnitude).
  # In the far field (r >> λ₀):
  #   E_intensity(r,z) ≈ C/R² × sin²(θ)    where R = sqrt(r²+z²), sin(θ) = r/R
  # This is the characteristic sin²θ dipole radiation pattern.
  # Plotting E_intensity in ParaView with the RZ solution reveals the
  # doughnut-shaped radiation lobe in 2D (which corresponds to a torus in 3D).
  [E_intensity]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Compute E_intensity = E_real² + E_imag² nodewise using ParsedAux.
  # This auxiliary variable is used both for visualisation and for
  # the PointValue postprocessors that extract the radiation pattern.
  [intensity_aux]
    type              = ParsedAux
    variable          = E_intensity
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # Radiation pattern: sample E_intensity at 5 angles on a circle of
  # radius R = 4 λ₀.
  #
  # Coordinates: (r, z) = (R·sin(θ), R·cos(θ)) with R = 4.0 units.
  # Expected: E_intensity(θ) / E_intensity(90°) ≈ sin²(θ)
  #
  # Note: the MOOSE PointValue postprocessor uses (x, y, z) coordinates
  # where for RZ: x = r, y = z (axial), z = 0 (azimuthal = 0).
  # ==================================================================

  # θ = 0° (on the z-axis): r = R·sin(0°) = 0, z = R·cos(0°) = 4.0
  # Expected: sin²(0°) = 0.0 — zero radiation along the dipole axis.
  # Due to the finite source size and ABC reflections, expect a small
  # but non-zero residual value (perhaps 1–5% of the equatorial intensity).
  [E_intensity_theta0]
    type     = PointValue
    variable = E_intensity
    point    = '0 4.0 0'
  []

  # θ = 30°: r = 4·sin(30°) = 2.0, z = 4·cos(30°) = 3.4641
  # Expected: sin²(30°) = 0.25 → E_intensity ≈ 0.25 × I_max
  [E_intensity_theta30]
    type     = PointValue
    variable = E_intensity
    point    = '2.0 3.4641 0'
  []

  # θ = 45°: r = 4·sin(45°) = 2.8284, z = 4·cos(45°) = 2.8284
  # Expected: sin²(45°) = 0.50 → E_intensity ≈ 0.50 × I_max
  [E_intensity_theta45]
    type     = PointValue
    variable = E_intensity
    point    = '2.8284 2.8284 0'
  []

  # θ = 60°: r = 4·sin(60°) = 3.4641, z = 4·cos(60°) = 2.0
  # Expected: sin²(60°) = 0.75 → E_intensity ≈ 0.75 × I_max
  [E_intensity_theta60]
    type     = PointValue
    variable = E_intensity
    point    = '3.4641 2.0 0'
  []

  # θ = 90°: r = 4·sin(90°) = 4.0, z = 4·cos(90°) = 0.0  (equatorial plane)
  # Expected: sin²(90°) = 1.0 → maximum radiation.
  # This is the point of maximum field intensity; the ratios
  # E_intensity(θ) / E_intensity(90°) give the normalised radiation pattern.
  [E_intensity_theta90]
    type     = PointValue
    variable = E_intensity
    point    = '4.0 0.0 0'
  []

  # ------------------------------------------------------------------
  # Additional diagnostics
  # ------------------------------------------------------------------

  # Near-field intensity at the source location (r=0, z=0).
  # This is the peak driving field; much larger than the far-field values.
  # Not part of the radiation pattern (near-field), but useful for
  # checking that the source is driving the field correctly.
  [E_intensity_source]
    type     = PointValue
    variable = E_intensity
    point    = '0 0 0'
  []

  # Intensity at θ = 120° (= 180° − 60°): by symmetry should equal θ = 60°.
  # r = 4·sin(120°) = 3.4641, z = 4·cos(120°) = −2.0
  # Checks the axial symmetry of the computed pattern.
  [E_intensity_theta120]
    type     = PointValue
    variable = E_intensity
    point    = '3.4641 -2.0 0'
  []

  # Intensity at θ = 150° (= 180° − 30°): by symmetry should equal θ = 30°.
  # r = 4·sin(150°) = 2.0, z = 4·cos(150°) = −3.4641
  [E_intensity_theta150]
    type     = PointValue
    variable = E_intensity
    point    = '2.0 -3.4641 0'
  []

  # Total radiated power: ∫ E_intensity dV (over the RZ half-domain).
  # In MOOSE RZ, ElementIntegralVariablePostprocessor integrates with the
  # cylindrical measure dV = r dr dz (the 2πr factor from azimuthal
  # integration is NOT included by default — it would need a separate
  # functor weight).  This gives the 2D integral ∫∫ |E|² r dr dz.
  # Useful for monitoring overall field level and checking energy balance.
  [total_radiated_power]
    type     = ElementIntegralVariablePostprocessor
    variable = E_intensity
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner: required because the EMRobinBC on
  # the three outer boundaries introduces off-diagonal Jacobian blocks
  # coupling E_real and E_imag at the boundary DOFs.
  #
  # Without full=true:
  #   - The SMP builds only diagonal blocks (E_real × E_real, E_imag × E_imag)
  #   - Off-diagonal blocks from EMRobinBC are ignored
  #   - The preconditioner is poor near boundaries → slow Newton convergence
  #
  # With full=true:
  #   - All Jacobian blocks including off-diagonal are included in the SMP
  #   - Direct LU solves the coupled E_real/E_imag system exactly
  #   - For a linear Helmholtz problem, Newton converges in 1 iteration
  #
  # The Helmholtz system is linear (E appears at most linearly in all
  # kernels and BCs), so the Newton iteration converges exactly in one
  # step regardless — the preconditioner choice mainly affects the
  # initial linear solve quality within that step.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation via MUMPS.
  # The 2D axisymmetric Helmholtz system has 60 × 120 × 2 ≈ 14,762 × 2
  # DOFs. MUMPS handles this size efficiently on a single processor.
  #
  # Direct LU is strongly preferred over iterative methods here because:
  #   (1) The Helmholtz operator ∇² + k₀² is indefinite: its eigenvalues
  #       span positive and negative values, making GMRES poorly conditioned.
  #   (2) The off-diagonal Robin BC coupling makes the system non-Hermitian,
  #       ruling out standard symmetric Krylov methods (CG, MINRES).
  #   (3) For a linear problem this size, direct LU takes ~1 second and
  #       requires no iteration — far more reliable than iterative solvers.
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       mumps'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: 2D (RZ half-plane) field output — E_real, E_imag, E_intensity.
  # Visualise in ParaView:
  #   - E_intensity: shows the radiation lobe pattern in the r-z plane.
  #     The doughnut shape (maximum at z=0, zero at r=0) is the sin²θ pattern.
  #   - Use "Warp by Scalar" or "Plot Over Line" along the arc r²+z²=16
  #     to extract the angular radiation pattern.
  #   - Apply "Reflect" filter with normal [−1,0,0] to show the full 2D cross
  #     section (both positive and negative r), revealing the symmetric lobe.
  #   - For 3D visualisation: use "Rotational Extrusion" around the z-axis
  #     to create the full torus-shaped radiation pattern.
  exodus = true

  # csv: postprocessor values — E_intensity at the five radiation-pattern
  # sample points plus source, symmetry checks, and total power.
  # Key result: ratios E_intensity(θ) / E_intensity(90°) should match sin²(θ).
  csv    = true
[]
