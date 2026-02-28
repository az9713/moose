# ============================================================
# Case 86: Half-Wave Dipole Antenna — Radiation Pattern
# Balanis, "Antenna Theory: Analysis and Design", 4th Ed. (2016), Ch. 4
# Kong, "Electromagnetic Wave Theory" (Wiley, 1986), Ch. 7
#
# MIT 6.635 Advanced Electromagnetism, Spring 2003
# Professor Jin Au Kong
# OCW: https://ocw.mit.edu/courses/6-635-advanced-electromagnetism-spring-2003/
#
# Lectures: 1  (Radiation, dipole antenna)
#           14 (Antenna arrays, linear dipole)
#
# ============================================================
# HISTORICAL AND PHYSICAL BACKGROUND
# ============================================================
#
# The half-wave dipole antenna is the canonical radiating element in
# antenna theory. Proposed in its modern resonant form by Heinrich Hertz
# in 1888, the half-wave dipole (total length L = λ/2) draws its name
# from the fact that its physical length equals one half of the operating
# wavelength. At this length the antenna presents a nearly-real input
# impedance (~73 Ω) and achieves a directivity of 1.64 dBi (Balanis,
# 4th Ed., Eq. 4-64).
#
# The defining characteristic of the half-wave dipole is its current
# distribution:
#
#   I(z) = I₀ cos(πz/L)    for |z| ≤ L/2 = λ/4
#
# This cosine distribution is the natural eigenfunction of the thin-wire
# integral equation and differs from the uniform (Hertzian) distribution
# of the infinitesimal dipole (Case 85). The cosine distribution produces
# the far-field radiation pattern:
#
#   E_θ(θ) ∝ [cos((π/2) cos θ)] / sin θ
#
# and the power gain pattern:
#
#   G(θ) ∝ cos²((π/2) cos θ) / sin²θ
#
# with maximum at θ = 90° (broadside, perpendicular to the dipole axis)
# and nulls at θ = 0° and θ = 180° (along the dipole axis).
#
# The half-wave dipole's pattern is slightly more directive than the
# Hertzian dipole (sin²θ pattern), with a narrower main lobe and the
# same two nulls at the poles.
#
# ============================================================
# AXISYMMETRIC (RZ) FORMULATION
# ============================================================
#
# The half-wave dipole is rotationally symmetric about the z-axis.
# All field quantities are independent of the azimuthal angle φ.
# MOOSE's coord_type = RZ reduces the 3D problem to 2D:
#
#   In RZ coordinates:  x ≡ r  (radial, r ≥ 0)
#                        y ≡ z  (axial, −∞ < z < +∞)
#
# The 3D Laplacian in cylindrical coordinates reduces to:
#
#   ∇²E_φ = (1/r)(∂/∂r)(r ∂E_φ/∂r) + ∂²E_φ/∂z²
#
# where E_φ is the sole non-zero field component for an axisymmetric
# source along the z-axis (for azimuthal mode m = 0).
#
# In the simplified scalar Helmholtz model used here (far-field pattern
# only), we treat E as the z-component of the electric field driven by
# the cosine current distribution, analogous to Case 85 but with the
# spatially extended source:
#
#   ∇²E + k₀² E = −source(r, z)
#
# where the source is zero except in a thin strip near r = 0,
# representing the distributed current along the dipole axis.
#
# ============================================================
# GEOMETRY
# ============================================================
#
#   r = 0               r = 6
#   |                    |
#   z = +6 ─────────────── (top, absorbing BC)
#   |                    |
#   | dipole arm (+z):   |
#   | J(r,z) = I₀cos(πz/L)  for r < 0.05, 0 < z < 0.25
#   |           ^axis       |
#   z = 0  ─── feed point ──
#   |                    |
#   | dipole arm (−z):   |
#   | J(r,z) = I₀cos(πz/L)  for r < 0.05, -0.25 < z < 0
#   |                    |
#   z = -6 ─────────────── (bottom, absorbing BC)
#   |                    |
#   left: r=0 axis       right: r=6 (absorbing BC)
#   (natural Neumann)
#
# Domain (RZ coords): r ∈ [0, 6],  z ∈ [−6, 6]
# Dipole: total length L = λ/2 = 0.5, current strip r < 0.05, |z| < 0.25
# Measurement arc: r = 4 (points at selected θ from z-axis)
#
# ============================================================
# GOVERNING EQUATION
# ============================================================
#
# Frequency-domain Helmholtz equation (scalar 2D axisymmetric):
#
#   ∇²E + k₀² E = −source(r, z)
#
# where:
#   k₀ = 2π / λ₀  (free-space wavenumber)
#   λ₀ = 1.0       (normalised wavelength)
#   k₀ = 2π = 6.2831853...   rad/unit
#   k₀² = 4π² = 39.4784176...  unit⁻²
#
# The source represents the impressed current density of the half-wave dipole:
#
#   source(r,z) = I₀ × cos(πz/L)   for r < r_wire, |z| < L/2
#                 0                  otherwise
#
# where L = λ/2 = 0.5, r_wire = 0.05 (thin strip near r = 0),
# and I₀ = 100 (amplitude chosen to give a clear radiation pattern).
#
# ============================================================
# REAL/IMAGINARY SPLITTING
# ============================================================
#
# The complex phasor field E = E_real + j E_imag satisfies (with
# real coefficient k₀² and the source on the imaginary part):
#
#   E_real:  ∇²E_r + k₀² E_r = 0          (no source)
#   E_imag:  ∇²E_i + k₀² E_i = −source    (source drives imaginary part)
#
# The source is placed on E_imag because the phasor of a real current
# J_z(t) = I₀ cos(ωt) is −jωμ₀I₀, which has only an imaginary component.
# The EMRobinBC couples E_real and E_imag through the boundary integrals.
#
# KERNEL SIGN CONVENTION (MOOSE):
#   Diffusion:     residual = +∫ ∇E · ∇v dV  → strong: −∇²E
#   ADMatReaction: residual = −rate × ∫ E v dV → strong: −rate × E
#   BodyForce:     residual = −∫ f v dV        → strong: −f (adds +f to LHS)
#
# Assembling E_imag equation strong form from MOOSE kernels:
#   −∇²E_i − k₀²E_i = source(r,z)
#   ⟺  ∇²E_i + k₀²E_i = −source(r,z)   ✓
#
# ============================================================
# BOUNDARY CONDITIONS
# ============================================================
#
# Left boundary (r = 0):
#   Natural Neumann (no explicit BC). In RZ coordinates, r = 0 is the
#   symmetry axis. The natural boundary condition ∂E/∂r = 0 at r = 0
#   is automatically satisfied by the axisymmetric formulation for
#   m = 0 (azimuthal mode). No explicit BC block needed.
#
# Right, top, bottom boundaries:
#   EMRobinBC absorbing (profile_func_real = 0): absorbs outgoing
#   cylindrical/spherical waves without injecting any incident field.
#   The first-order ABC: ∂E/∂n + j k₀ E = 0.
#
# ============================================================
# NUMERICAL VALUES
# ============================================================
#
# λ₀ = 1.0 unit  →  k₀ = 2π = 6.28318530717958648 rad/unit
# k₀² = 4π² = 39.47841760435743 unit⁻²
# L    = λ₀/2 = 0.5 (dipole half-length = L/2 = 0.25)
# σ    = source strip width: r < r_wire = 0.05
# π/L  = 2π   (so cos(πz/L) = cos(2πz) — full cosine over λ/2)
#
# Source amplitude: 100.0 × cos(πz/L) = 100.0 × cos(2π z / 0.5)
#   = 100.0 × cos(4π z)   in the strip r < 0.05, |z| < 0.25
#
# Wait — note: π/L = π/0.5 = 2π, so cos(πz/L) = cos(2πz).
# For |z| < L/2 = 0.25 this covers exactly one quarter period of cos(2πz)
# on each arm, giving the positive lobe of the cosine.
# The ParsedFunction below uses the exact expression π/L × z = (π/0.5)×z.
#
# Measurement arc at r_meas = 4 units (4 λ₀ from source):
#   θ = 10° from z-axis:  r = 4 sin(10°) = 0.6946,  z = 4 cos(10°) = 3.9392
#   θ = 30° from z-axis:  r = 4 sin(30°) = 2.0000,  z = 4 cos(30°) = 3.4641
#   θ = 45° from z-axis:  r = 4 sin(45°) = 2.8284,  z = 4 cos(45°) = 2.8284
#   θ = 60° from z-axis:  r = 4 sin(60°) = 3.4641,  z = 4 cos(60°) = 2.0000
#   θ = 75° from z-axis:  r = 4 sin(75°) = 3.8637,  z = 4 cos(75°) = 1.0353
#   θ = 90° from z-axis:  r = 4 sin(90°) = 4.0000,  z = 4 cos(90°) = 0.0000
#
# Expected normalised gain pattern (Balanis, Eq. 4-58a):
#   G(θ) ∝ cos²((π/2) cos θ) / sin²θ
#   G(10°) = 0.0189  (very small, close to axis null)
#   G(30°) = 0.1746
#   G(45°) = 0.3943
#   G(60°) = 0.6667
#   G(75°) = 0.9042
#   G(90°) = 1.0000  (maximum — broadside to dipole axis)
#
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables (referenced with ${...})
# -----------------------------------------------------------
# k   = k₀ = 2π (free-space wavenumber, λ₀ = 1 normalised unit)
# E0  = 0  (no incident wave from port — source is interior distributed current)
# theta = 0  (incidence angle placeholder for EMRobinBC, normal incidence)
k     = 6.283185307179586    # k₀ = 2π rad/unit (λ₀ = 1 unit)
E0    = 0                    # no external incident wave (interior source only)
theta = 0                    # cosTheta = 1 for first-order ABC

[Mesh]
  # 2D axisymmetric domain: r ∈ [0, 6],  z ∈ [−6, 6].
  # coord_type = RZ converts this to cylindrical coordinates.
  # In MOOSE RZ: x ≡ r (horizontal axis), y ≡ z (vertical axis).
  #
  # Resolution:
  #   60 elements over r ∈ [0, 6] → Δr = 0.1 (10 elements per λ₀ = 1)
  #   120 elements over z ∈ [−6, 6] → Δz = 0.1 (10 elements per λ₀)
  #   Free-space wavelength λ₀ = 1 → 10 elements per wavelength (adequate for
  #   FIRST-order Lagrange; 15–20 preferred for higher accuracy).
  #   Source strip r < 0.05 = 0.5 Δr — marginal, but sufficient to
  #   capture the cosine distribution in the far field where λ₀/r ≪ 1.
  #
  # The domain extends 6λ₀ in both r and z, giving ~6 free-space
  # wavelengths of clearance from the dipole (at r ≈ 0, z = 0) to the
  # absorbing boundary. This is adequate for the first-order Robin ABC.
  # coord_type = RZ activates the 2D axisymmetric (cylindrical) coordinate
  # formulation. MOOSE modifies the weak form integrals to include the
  # Jacobian factor r (the radial distance), replacing Cartesian area
  # elements dA = dx dy with cylindrical area elements dA = r dr dz.
  coord_type = RZ

  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 60
    ny   = 120
    xmin = 0
    xmax = 6
    ymin = -6
    ymax = 6
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor.
  # In the RZ axisymmetric model, E represents E_φ (azimuthal component
  # for m=0 mode) or equivalently the z-component of the far field.
  # FIRST/LAGRANGE: standard continuous Galerkin with linear basis functions.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex phasor E = E_real + j E_imag.
  # The distributed current source (cosine distribution along z) drives
  # E_imag directly through the BodyForce kernel. E_real is driven by
  # coupling through the absorbing boundary conditions.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # source_fn: Distributed cosine current source — the half-wave dipole
  #
  # Models the impressed current density J_z(r, z) of a half-wave dipole:
  #
  #   J(r, z) = I₀ × cos(πz/L)    if r < r_wire AND |z| < L/2
  #             0                   otherwise
  #
  # where:
  #   I₀    = 100     (normalised amplitude)
  #   L     = 0.5     (total dipole length = λ₀/2)
  #   L/2   = 0.25    (half-dipole length = λ₀/4)
  #   r_wire = 0.05   (source strip width, approximately 5% of λ₀)
  #   π/L   = π/0.5 = 2π  (spatial frequency of cosine)
  #
  # In RZ coordinates: x = r, y = z.
  # The ParsedFunction expression uses 'x' for r and 'y' for z.
  #
  # Condition: x < 0.05 (within strip) AND abs(y) < 0.25 (within dipole arm)
  # Using the boolean short-circuit: if(A & B, val, 0)
  # Note: MOOSE ParsedFunction uses '&' for logical AND.
  #
  # Physical rationale for the cosine distribution:
  # The current along a centre-fed thin dipole satisfies the wire antenna
  # integral equation. For a resonant half-wave dipole the sinusoidal
  # standing wave I(z) = I₀ sin(k₀(L/2 − |z|)) = I₀ cos(k₀z) (for L = λ/2)
  # emerges naturally. The cosine form correctly gives I = 0 at the tips
  # (z = ±L/2) and I = I₀ at the feed (z = 0). This is physically
  # required because current cannot flow past the open endpoints.
  # ------------------------------------------------------------------
  [source_fn]
    type       = ParsedFunction
    expression = 'if(x < 0.05 & abs(y) < 0.25, 100.0 * cos(3.14159265358979323846 * y / 0.5), 0.0)'
  []

  # cosTheta = cos(θ) for the EMRobinBC first-order ABC.
  # At normal incidence (θ = 0), cos(θ) = 1 and the ABC is:
  #   ∂E/∂n + j k₀ E = 0
  # For an interior point source producing outgoing cylindrical waves,
  # the first-order ABC is most effective at near-normal incidence.
  # The 6λ₀ domain ensures the waves are nearly plane-wave-like at the
  # boundary, minimising reflection from the ABC.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # ADGenericConstantMaterial: wrap the uniform free-space k₀² as an AD
  # material property, required by ADMatReaction.
  #
  # The half-wave dipole radiates in free space (ε_r = 1, μ_r = 1)
  # throughout the entire RZ domain. There is no material inhomogeneity
  # (unlike Cases 77, 79, 83). The spatially-varying factor is only in
  # the source term, not in the wave equation coefficients.
  #
  # k0sq = k₀² = (2π)² = 4π² = 39.47841760435743 unit⁻²
  [free_space]
    type        = ADGenericConstantMaterial
    prop_names  = 'k0sq'
    prop_values = '39.47841760435743'
  []
[]

[Kernels]
  # ==================================================================
  # E_real EQUATION:
  #   ∇²E_r + k₀² E_r = 0
  #
  # In MOOSE axisymmetric RZ, the Diffusion kernel automatically
  # includes the cylindrical Jacobian factor r in the weak form:
  #   ∫ r ∇E_r · ∇v dA  − k₀² ∫ r E_r v dA = 0
  # which corresponds to the strong form:
  #   (1/r)∂(r ∂E_r/∂r)/∂r + ∂²E_r/∂z² + k₀² E_r = 0
  # ==================================================================

  # Laplacian term for E_real: −∇²E_r in the RZ cylindrical sense.
  # MOOSE applies the coord_type = RZ modification to the Diffusion kernel
  # automatically, so no special cylindrical kernel is needed.
  [diff_real]
    type     = Diffusion
    variable = E_real
  []

  # Reaction term: −k₀² E_r in strong form.
  # ADMatReaction residual = −k0sq × ∫ E_r v dA → strong: −k0sq × E_r.
  # Combined with Diffusion: −∇²E_r − k₀²E_r = 0 ↔ ∇²E_r + k₀²E_r = 0 ✓
  [helm_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq
  []

  # No BodyForce for E_real: the current source −jωμ₀J has only an
  # imaginary phasor component, so the source drives E_imag only.
  # E_real is excited only through the off-diagonal coupling in the
  # EMRobinBC (the j k₀ term couples real ↔ imaginary at the boundary).

  # ==================================================================
  # E_imag EQUATION:
  #   ∇²E_i + k₀² E_i = −source(r, z)
  #
  # The distributed cosine source drives E_imag. BodyForce adds
  # −∫ source_fn × v dA, so in strong form the equation has
  # source_fn on the RHS: ∇²E_i + k₀² E_i = source_fn.
  # This corresponds to the −jωμ₀ J_z source (imaginary part).
  # ==================================================================

  # Laplacian term for E_imag: same cylindrical Helmholtz operator as E_real.
  [diff_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Reaction term for E_imag: −k₀² E_i in strong form.
  # Same k0sq material property as for E_real — free space throughout.
  [helm_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq
  []

  # Distributed cosine source for E_imag — the half-wave dipole current.
  # BodyForce adds −∫ source_fn × v dA to the residual.
  # Strong form: E_imag equation has +source_fn on the RHS.
  # This source is non-zero only in the thin strip near r = 0 where the
  # dipole current flows, simulating the antenna wire with a finite width
  # approximation. The strip width r < 0.05 ≈ λ₀/20, small enough that
  # the far-field pattern is well-approximated by a line source.
  [source_imag]
    type     = BodyForce
    variable = E_imag
    function = source_fn
  []
[]

[BCs]
  # ==================================================================
  # LEFT BOUNDARY (r = 0): SYMMETRY AXIS — NO EXPLICIT BC
  # ==================================================================
  # In RZ coordinates the left boundary (x = 0, i.e. r = 0) is the
  # symmetry axis. For the m = 0 axisymmetric mode the natural boundary
  # condition ∂E/∂r = 0 applies automatically. No explicit BC block is
  # needed: MOOSE's weak form with the RZ Jacobian naturally enforces
  # the axis condition through the Neumann term which vanishes at r = 0
  # (the cylindrical Jacobian r = 0 makes the boundary integral zero).
  #
  # MOOSE's RZ formulation: weak form has ∫ r ∇E · ∇v dA. On the axis
  # (r = 0) the boundary flux r × ∂E/∂n = 0 × ∂E/∂r = 0 identically.
  # This is the correct physical condition for m = 0 symmetry.

  # ==================================================================
  # RIGHT BOUNDARY (r = 6): ABSORBING BC
  # ==================================================================
  # Absorbs outgoing cylindrical waves from the dipole.
  # The first-order ABC ∂E/∂r + j k₀ E = 0 is appropriate since at
  # r = 6 = 6λ₀ the outgoing waves are nearly plane-wave-like.
  # No incident wave injection (profile_func_real = 0).
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
  # TOP BOUNDARY (z = +6): ABSORBING BC
  # ==================================================================
  # Absorbs outgoing waves from the upper dipole arm region.
  # The domain extends 6λ₀ above the dipole tip (at z = +0.25), giving
  # ample clearance for the first-order ABC to be effective.
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
  # BOTTOM BOUNDARY (z = −6): ABSORBING BC
  # ==================================================================
  # Absorbs outgoing waves from the lower dipole arm region.
  # Symmetric to the top boundary — the dipole is symmetric about z = 0.
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
  # Proportional to the time-averaged radiated power density.
  # In the far field (r ≫ λ₀), E_intensity varies as 1/r² multiplied
  # by the gain pattern G(θ). Plotting E_intensity × r² in ParaView
  # (using a Calculator filter) extracts the angular gain pattern.
  [E_intensity]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Compute |E|² nodewise from the primary field components.
  # This auxiliary variable captures the radiated field intensity
  # at each point in the RZ domain.
  [intensity_aux]
    type              = ParsedAux
    variable          = E_intensity
    expression        = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # Field intensity at selected points on a measurement arc at r_meas = 4
  # in RZ coordinates (x = r, y = z).
  #
  # The angle θ is measured from the +z axis (dipole axis). The expected
  # normalised gain pattern is G(θ) ∝ cos²((π/2)cosθ) / sin²θ, which:
  #   - Peaks at θ = 90° (broadside, equatorial plane z = 0)
  #   - Decreases toward θ = 0°, 180° (along dipole axis)
  #   - Goes to zero at θ = 0° and θ = 180° (axial nulls)
  #
  # All points are in the y ≥ 0 (z ≥ 0) half, sufficient to map the
  # pattern due to symmetry. The pattern in the z < 0 hemisphere is
  # identical by symmetry of the cosine source distribution.
  #
  # Mapping: x = r = R_meas × sin(θ),  y = z = R_meas × cos(θ)
  # R_meas = 4 (4 free-space wavelengths from the dipole — far field)
  # ==================================================================

  # θ = 10° from z-axis: r = 4 sin(10°) = 0.6946, z = 4 cos(10°) = 3.9392
  # Near the dipole axis — G(10°) ≈ 0.019 (very small)
  [E_int_theta10]
    type     = PointValue
    variable = E_intensity
    point    = '0.694593 3.939231 0'
  []

  # θ = 30°: r = 4 sin(30°) = 2.000, z = 4 cos(30°) = 3.4641
  # G(30°) ≈ 0.175 — pattern is building toward broadside
  [E_int_theta30]
    type     = PointValue
    variable = E_intensity
    point    = '2.000000 3.464102 0'
  []

  # θ = 45°: r = 4 sin(45°) = 2.8284, z = 4 cos(45°) = 2.8284
  # G(45°) ≈ 0.394 — half-way between axis and broadside
  [E_int_theta45]
    type     = PointValue
    variable = E_intensity
    point    = '2.828427 2.828427 0'
  []

  # θ = 60°: r = 4 sin(60°) = 3.4641, z = 4 cos(60°) = 2.000
  # G(60°) ≈ 0.667 — two-thirds of maximum
  [E_int_theta60]
    type     = PointValue
    variable = E_intensity
    point    = '3.464102 2.000000 0'
  []

  # θ = 75°: r = 4 sin(75°) = 3.8637, z = 4 cos(75°) = 1.0353
  # G(75°) ≈ 0.904 — near maximum
  [E_int_theta75]
    type     = PointValue
    variable = E_intensity
    point    = '3.863703 1.035276 0'
  []

  # θ = 90°: r = 4 sin(90°) = 4.000, z = 4 cos(90°) = 0.000
  # G(90°) = 1.000 — maximum, the equatorial plane (broadside)
  # The half-wave dipole radiates most strongly in the plane perpendicular
  # to the dipole axis (the z = 0 plane). This is the reference point
  # for normalising the gain pattern.
  [E_int_theta90]
    type     = PointValue
    variable = E_intensity
    point    = '4.000000 0.000000 0'
  []

  # On-axis probe (r = 0.1, z = 1.0): should be near-zero due to the
  # null in the radiation pattern along the dipole axis.
  # Placed at r = 0.1 (slightly off-axis to avoid r = 0 singularity).
  [E_int_axis]
    type     = PointValue
    variable = E_intensity
    point    = '0.1 1.0 0'
  []

  # Feed point field (r = 0.1, z = 0): near the dipole centre.
  # This is the maximum current location — highest source amplitude.
  # The near-field intensity here should be much larger than the far-field
  # values at the measurement arc (near-to-far-field ratio).
  [E_int_feed]
    type     = PointValue
    variable = E_intensity
    point    = '0.1 0.0 0'
  []

  # Total radiated field energy: ∫ |E|² dA (cylindrical area element in RZ).
  # Provides a measure of total radiated power (proportional to ∫ r|E|² dr dz).
  # Useful for monitoring convergence and normalising the gain pattern.
  [total_radiated_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = E_intensity
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner: required for the coupled
  # (E_real, E_imag) system. The EMRobinBC on the right, top, and bottom
  # boundaries introduces off-diagonal Jacobian blocks between E_real and
  # E_imag through the j k₀ coupling term:
  #   j k₀ E = j k₀ (E_real + j E_imag) → couples E_real ↔ E_imag
  # Without full = true, these off-diagonal blocks are ignored and the
  # Newton solver has incorrect Jacobian structure, leading to divergence.
  #
  # The 2D RZ system has 60 × 120 = 7,200 elements, 7,381 nodes, and
  # ~14,762 DOFs (two variables on FIRST-order mesh). Direct LU is
  # fast at this scale (< 1 second on modern hardware).
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation: the axisymmetric Helmholtz system is linear
  # (E appears only linearly in all kernels and BCs). Newton converges in
  # 1 iteration. Direct LU is preferred over iterative solvers because:
  #   (1) The Helmholtz operator is indefinite (k₀² term can dominate the
  #       positive-definite Laplacian at lower mesh resolutions).
  #   (2) The off-diagonal EMRobinBC blocks make the system non-symmetric.
  #   (3) At 14,762 DOFs, LU is computationally inexpensive.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: 2D axisymmetric field output — E_real, E_imag, E_intensity
  #   on the full RZ mesh. Visualise in ParaView:
  #   - E_intensity in the r-z plane reveals the dipole radiation pattern
  #     as rings of decreasing intensity in the equatorial plane (z = 0)
  #     and near-zero intensity along the z-axis (dipole null)
  #   - The cosine source distribution is visible as a concentration of
  #     high intensity near r ≈ 0, |z| < 0.25 (the dipole wire region)
  #   - Use the Revolve filter in ParaView to visualise the full 3D
  #     toroidal radiation pattern from the 2D RZ solution
  #   - Plot E_intensity vs. angle θ along a circle at r = 4 to extract
  #     the gain pattern and compare with the analytic formula
  exodus = true

  # csv: postprocessor values — E_intensity at θ = 10°, 30°, 45°, 60°,
  #   75°, 90° on the R = 4 measurement arc, plus on-axis null check
  #   and total radiated energy. The ratios E_int_thetaN / E_int_theta90
  #   give the normalised gain pattern for comparison with theory.
  csv    = true
[]
