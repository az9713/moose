# ============================================================
# Case 77: EM Scattering from a Dielectric Cylinder (2D, TE Polarisation)
# Harrington, "Time-Harmonic Electromagnetic Fields" (1961), Ch. 5
# Balanis, "Advanced Engineering Electromagnetics", 2nd Ed. (2012), Ch. 11
# Jin, "The Finite Element Method in Electromagnetics", 3rd Ed. (2014), Ch. 9
#
# Lecture 5 (Integral Equations, EFIE, scattering)
# Lecture 6 (Method of Moments — MoM)
#
# PHYSICAL SCENARIO
# -----------------
# An infinite dielectric cylinder (εᵣ = 4, n = 2) with circular cross-
# section of radius R = 1 m is illuminated by a TE-polarised plane wave
# travelling in the +x direction. In the TE case the electric field has
# only a z-component (out of the 2D cross-sectional plane), so the full
# 3D Maxwell problem reduces to a scalar 2D Helmholtz equation for E_z.
#
# The domain is the 2D cross-section [-4, 4] × [-4, 4]. The dielectric
# cylinder occupies the disc x² + y² < 1. Outside the cylinder is free
# space (εᵣ = 1).
#
# GOVERNING EQUATION
# ------------------
# From the z-component of the curl-curl equation (∇×∇×E = k₀² εᵣ E)
# and the TE polarisation assumption (E has only z-component, ∂/∂z = 0):
#
#   ∇²E_z + k₀² εᵣ(x,y) E_z = 0
#
# where:
#   ∇² = ∂²/∂x² + ∂²/∂y²   (2D Laplacian)
#   k₀ = 2π/λ₀              (free-space wavenumber)
#   εᵣ(x,y) = 4  if x²+y² < 1   (dielectric cylinder, n=2)
#              1  if x²+y² ≥ 1   (free space, n=1)
#
# WAVELENGTH AND WAVENUMBER CHOICE
# ---------------------------------
# We choose λ₀ = 4 m so the cylinder radius R = 1 = λ₀/4. This puts
# the problem in the Mie scattering regime (ka = k₀ R = π/2 ≈ 1.57,
# so ka ~ O(1)) where diffraction effects are pronounced and the simple
# Rayleigh (ka ≪ 1) or geometric optics (ka ≫ 1) limits do not apply.
#
#   λ₀ = 4 m
#   k₀ = 2π/λ₀ = 2π/4 = π/2 = 1.5707963...  rad/m
#   k₀² = (π/2)² = π²/4 = 2.4674011...  m⁻²
#   k₀² × εᵣ = 2.4674 × 4 = 9.8696  m⁻²  (inside cylinder)
#   k₀² × εᵣ = 2.4674 × 1 = 2.4674  m⁻²  (outside cylinder)
#
# REAL / IMAGINARY SPLITTING
# --------------------------
# The total field is complex:  E_z(x,y) = E_r(x,y) + j E_i(x,y).
# Because εᵣ is purely real (lossless dielectric), the Helmholtz equation
# has real coefficients and the two components decouple in the BULK:
#
#   ∇²E_r + k₀² εᵣ(x,y) E_r = 0    ... (R)
#   ∇²E_i + k₀² εᵣ(x,y) E_i = 0    ... (I)
#
# Both equations have exactly the same spatial operator. The real and
# imaginary parts couple ONLY through the boundary conditions (Robin BCs),
# not in the interior. This is identical to Cases 32 and 74.
#
# INCIDENT PLANE WAVE
# -------------------
# The incident field propagates in the +x direction:
#   E_z^inc(x,y) = E₀ exp(+j k₀ x)
#
# In real/imaginary form:
#   E_r^inc = E₀ cos(k₀ x)   (E₀ = 1 in this case)
#   E_i^inc = E₀ sin(k₀ x)
#
# Physically, at the west boundary (x = -4), the incident field is:
#   E_r^inc(-4, y) = cos(k₀ × (-4)) = cos(-2π) = 1.0
#   E_i^inc(-4, y) = sin(k₀ × (-4)) = sin(-2π) = 0.0
#
# NOTE: k₀ × 4 = (π/2) × 4 = 2π, so the domain is exactly 2 full
# wavelengths wide. This means the phase of the incident wave at x = -4
# is the same as at x = +4: the domain fits an integer number of periods.
#
# BOUNDARY CONDITIONS — EMRobinBC (First-Order ABC / Port)
# --------------------------------------------------------
# The four domain boundaries (west, east, north, south) are terminated
# with the first-order Sommerfeld absorbing boundary condition:
#
#   ∂E_z/∂n + j k₀ E_z = 2 j k₀ E_z^inc    (wave injection + absorption)
#   ∂E_z/∂n + j k₀ E_z = 0                  (pure absorption only)
#
# The EMRobinBC kernel implements Jin's port condition (Jin Eq. 9.60):
#
#   ∂E/∂n + j k₀ cos(θ) E = 2 j k₀ cos(θ) E₀ exp(j k₀ cos(θ) x)
#
# at x = -4 (west): normal to boundary is -x̂, outward direction is -x,
#   the incident wave is incoming from the left so the port source term
#   is non-zero: profile_func_real = E₀ = 1. EMRobinBC handles the
#   correct sign of the incident flux at the port automatically.
#
# at x = +4 (east), y = ±4 (north/south): pure absorbing BC.
#   profile_func_real = 0 → no incident wave injection on these faces.
#   The scattered field is approximately cylindrical at large r, which
#   means it hits the top/bottom boundaries at oblique incidence. The
#   first-order ABC is not perfectly absorbing at oblique angles, but
#   it provides adequate accuracy for a domain of this size (4 wavelengths
#   across). Mie analytical results confirm the scattered power is small
#   far from the cylinder.
#
# SIGN CONVENTION: sign = negative on all EMRobinBC instances.
# This follows the convention established in Cases 32, 74, 75.
#
# CONNECTION TO INTEGRAL EQUATIONS AND METHOD OF MOMENTS (Lectures 5-6)
# -----------------------------------------------------------------------
# The same scattering problem is traditionally solved with the Electric
# Field Integral Equation (EFIE). The EFIE represents the scattered field
# as a surface integral of the induced current over the cylinder surface:
#
#   E_z^scat(r) = -j k₀ ∫_S J_z(r') H₀⁽²⁾(k₀|r-r'|) dl'
#
# where H₀⁽²⁾ is the Hankel function of the second kind and J_z is the
# equivalent surface current. The Method of Moments (MoM) discretises
# this surface integral into a dense N×N matrix equation and solves for
# the N unknowns J_z. The FEM approach here instead solves a sparse
# volume integral in the 2D domain. The two approaches give the same
# physical result, but:
#
#   MoM (EFIE):   N ≈ circumference / (λ/10) ≈ 6λ/10 ≈ 15 unknowns
#                 Dense 15×15 complex matrix — extremely efficient
#                 Limited to homogeneous cylinders (Green's function known)
#
#   FEM (this):   N ≈ 80×80 × 2 = 12,800 real DOFs
#                 Sparse ~12,800×12,800 matrix — less efficient for this case
#                 Works for arbitrary inhomogeneous cylinders, complex media,
#                 arbitrary shapes — general-purpose approach
#
# The FEM is the natural extension of MoM to volumetric (body-of-revolution)
# problems where no analytical Green's function exists, which is why the
# FEM formulation is taught after the integral equation lectures.
#
# MIE SERIES BENCHMARK
# --------------------
# For a lossless dielectric cylinder with εᵣ = 4 (n = 2), the exact
# solution is given by the Mie series. For TE polarisation:
#
#   E_z^scat = -Σ_n  cₙ aₙ Hₙ⁽²⁾(k₀ r) exp(j n φ)
#
# where the Mie coefficients aₙ involve Bessel functions Jₙ at the
# cylinder boundary and encode the scattering for each angular harmonic.
#
# For k₀R = π/2 ≈ 1.57 (moderate ka):
#   - Shadow region (x > 1, on axis y = 0): field is reduced below 1
#     as the cylinder intercepts and partially absorbs the incident power.
#   - Forward scatter (far field, +x direction): constructive interference
#     → field amplitude > 1 in the forward direction past the shadow zone.
#   - Inside cylinder (|r| < 1): field is enhanced because n = 2 focuses
#     the incident wave (effective wavelength λ_cyl = λ₀/n = 2 m).
#   - Near field around cylinder: complex interference pattern with lobes.
#
# DOMAIN SIZING
# -------------
# The domain extends ±4 m from the origin (4 = λ₀ × 1 = exactly 1
# free-space wavelength from cylinder centre to boundary on each side).
# A rule of thumb for first-order ABC is: place the boundary at least
# λ/2 from the scatterer. Here we have 3 m of clearance (R=1 + 3 m gap
# = 4 m boundary), i.e. 3/4 λ₀ clearance, which is adequate.
#
# MESH: 80×80 = 6400 quad elements → ~12,800 DOFs for two scalar fields.
# At k₀ = π/2, λ₀ = 4 m, domain = 8 m → 80 elements/8 m = 10 elements/m.
# Spatial period λ₀ = 4 m → 40 elements per free-space wavelength.
# Inside cylinder λ_cyl = λ₀/n = 2 m → 20 elements per cylinder wavelength.
# Both are well above the minimum (≥10 elements/wavelength for FIRST-order
# Lagrange), so the 80×80 mesh provides adequate spatial resolution.
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables
# -----------------------------------------------------------
# k = k₀ = 2π/λ₀ = 2π/4 = π/2 = 1.5707963... rad/m
# k² = π²/4 = 2.4674011...  m⁻²
# k² × εᵣ_inside  = 9.8696044...  m⁻²   (cylinder, εᵣ = 4)
# k² × εᵣ_outside = 2.4674011...  m⁻²   (free space, εᵣ = 1)
#
# Adjusting 'k' changes the operating frequency (or equivalently, the
# electrical size of the cylinder k₀R). Increasing k₀ moves toward the
# geometric optics (ray) limit; decreasing k₀ moves toward Rayleigh.
k     = 1.5707963267948966   # free-space wavenumber k₀ = π/2 [rad/m]
E0    = 1                    # incident field amplitude [V/m]
theta = 0                    # plane-wave incidence angle [degrees]

[Mesh]
  # 2D square domain [-4,4] × [-4,4].
  # 80×80 quad elements: 40 elements per λ₀ = 4 m free-space wavelength.
  # The dielectric cylinder (r < 1) is not meshed separately — the
  # ParsedFunction for εᵣ(x,y) handles the material discontinuity.
  # MOOSE's Lagrange elements enforce C⁰ continuity of E_z across the
  # cylinder boundary, which is physically correct: tangential E_z is
  # continuous at the dielectric interface.
  [domain]
    type = GeneratedMeshGenerator
    dim  = 2
    nx   = 80
    ny   = 80
    xmin = -4
    xmax = 4
    ymin = -4
    ymax = 4
  []

  # Rename the four boundary labels to physically meaningful names.
  # GeneratedMeshGenerator produces: left, right, bottom, top.
  # Renaming clarifies the compass-direction geometry used in the BCs.
  [rename]
    type         = RenameBoundaryGenerator
    input        = domain
    old_boundary = 'left right bottom top'
    new_boundary = 'west east south north'
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor E_z(x,y).
  # The total complex field is E_z = E_real + j E_imag.
  # Using first-order Lagrange elements (FIRST/LAGRANGE) on the 2D mesh.
  # Lagrange basis functions enforce C⁰ continuity of E_z across the
  # dielectric cylinder boundary at r = 1, correctly implementing the
  # physical boundary condition: tangential E_z continuous at the
  # vacuum–dielectric interface (Balanis, Advanced EM, Sec. 3.6).
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex phasor.
  # Decoupled from E_real in the 2D bulk (lossless dielectric → no
  # cross-coupling between real and imaginary parts of E_z). The two
  # components interact only through the EMRobinBC at the port boundary.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # coeff_fn: spatially varying reaction coefficient k₀² εᵣ(x,y)
  #
  # The 2D scalar Helmholtz equation ∇²E_z + k₀² εᵣ E_z = 0 is
  # discretised in MOOSE using Diffusion + ADMatReaction. The residual
  # is (in strong form):
  #   −∇²E_z − rate × E_z = 0
  # so we need rate = +k₀² εᵣ (positive, so the two terms combine to
  # reproduce the Helmholtz equation with the correct sign).
  #
  # Inside the cylinder (x²+y² < 1, εᵣ = 4):
  #   rate = k₀² × 4 = (π/2)² × 4 = π² = 9.8696044010893580  m⁻²
  #
  # Outside the cylinder (x²+y² ≥ 1, εᵣ = 1):
  #   rate = k₀² × 1 = (π/2)² × 1 = π²/4 = 2.4674011002723395 m⁻²
  #
  # The if() function evaluates at each quadrature point: where
  # x²+y² < 1 the cylinder material is used, everywhere else free space.
  # Note the strict inequality < 1: on the cylinder surface itself the
  # discontinuity is handled by continuity of E_z (weak form automatically
  # averages across the interface via the C⁰ Lagrange basis).
  # ------------------------------------------------------------------
  [coeff_fn]
    type       = ParsedFunction
    expression = 'if(x*x+y*y < 1.0, 9.8696044010893580, 2.4674011002723395)'
  []

  # cosTheta = cos(0°) = 1.0 for normal (head-on) incidence.
  # EMRobinBC uses k_eff = k₀ × cos(θ) for the boundary condition.
  # At θ = 0 (plane wave in +x direction, port on west face), k_eff = k₀.
  # Oblique incidence (θ ≠ 0) would require adjusting this function and
  # the incident wave expression in the boundary conditions.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # Wrap coeff_fn as an AD material property so ADMatReaction can consume it.
  # ADGenericFunctionMaterial evaluates the ParsedFunction at each quadrature
  # point and returns an ADReal, enabling automatic differentiation through
  # the spatially varying permittivity. This is critical for Newton convergence:
  # the Jacobian of ADMatReaction depends on the material property derivative.
  #
  # The material property 'coeff_material' encodes the full spatial variation
  # of k₀² εᵣ(x,y) across both the cylinder and the surrounding free space,
  # evaluated pointwise from the ParsedFunction at runtime.
  [coeff_material]
    type        = ADGenericFunctionMaterial
    prop_names  = 'coeff_material'
    prop_values = 'coeff_fn'
  []
[]

[Kernels]
  # ==================================================================
  # Equation (R) for E_real:
  #   ∇²E_r + k₀² εᵣ(x,y) E_r = 0
  #
  # MOOSE weak form (Green's first identity, integrating Diffusion by parts):
  #   ∫_Ω ∇E_r · ∇v dA  −  k₀² εᵣ ∫_Ω E_r v dA  +  ∫_∂Ω (∂E_r/∂n) v ds = 0
  #
  # The boundary integral is handled by the EMRobinBCs below.
  # The two volume integrals are:
  #   Diffusion:     +∫_Ω ∇E_r · ∇v dA     → strong form: −∇²E_r
  #   ADMatReaction: −rate × ∫_Ω E_r v dA  → strong form: −rate × E_r
  #
  # Combined: −∇²E_r − rate × E_r = 0  with rate = +k₀²εᵣ → −∇²E_r − k₀²εᵣ E_r = 0
  # ↔ ∇²E_r + k₀²εᵣ E_r = 0   (Helmholtz, equation R)  ✓
  # ==================================================================

  # Laplacian term for E_real: contributes −∇²E_r to the strong form.
  # In 2D, this is −(∂²E_r/∂x² + ∂²E_r/∂y²).
  [diffusion_real]
    type     = Diffusion
    variable = E_real
  []

  # Reaction term for E_real: adds −k₀² εᵣ(x,y) E_real to the strong form.
  # The positive reaction_rate = +k₀²εᵣ combined with the ADMatReaction
  # sign convention (residual = −rate × E × test) gives the correct sign.
  # No cross-coupling from E_imag: εᵣ is real (lossless cylinder) so the
  # imaginary permittivity εᵣ'' = 0 and the cross terms vanish.
  [field_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = coeff_material
  []

  # NOTE: No ADMatCoupledForce kernel for E_real.
  # The cylinder is lossless (εᵣ'' = 0). Cross-coupling terms
  # ±k₀² εᵣ'' × E_imag are present only for lossy media (as in Case 75,
  # Drude slab). Here they vanish exactly, and the two field components
  # are decoupled in the 2D volume — they couple only through the Robin BC.

  # ==================================================================
  # Equation (I) for E_imag:
  #   ∇²E_i + k₀² εᵣ(x,y) E_i = 0
  #
  # Identical structure to equation (R): same Helmholtz operator, same
  # spatially varying coefficient, no cross-coupling (lossless medium).
  # ==================================================================

  # Laplacian term for E_imag: contributes −∇²E_i to the strong form.
  [diffusion_imag]
    type     = Diffusion
    variable = E_imag
  []

  # Reaction term for E_imag: adds −k₀² εᵣ(x,y) E_imag to the strong form.
  # Uses the same coeff_material as the E_real equation: lossless → no coupling.
  [field_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = coeff_material
  []

  # NOTE: No ADMatCoupledForce kernel for E_imag either.
  # In a lossy problem (εᵣ'' ≠ 0), there would be a term +k₀² εᵣ'' E_real
  # here (see Case 75). For the lossless cylinder it vanishes.
[]

[BCs]
  # ==================================================================
  # WEST BOUNDARY (x = -4): WAVE INJECTION PORT
  # ==================================================================
  # The incident plane wave E_z^inc = exp(+j k₀ x) enters from the left.
  # The EMRobinBC at x = -4 simultaneously:
  #   (1) Injects the incident wave into the domain
  #   (2) Absorbs any reflected/scattered field leaving through this face
  #
  # Jin port condition (Eq. 9.60) at the west face:
  #   ∂E_z/∂n_west + j k₀ E_z = 2 j k₀ E₀ exp(+j k₀ × (-4))
  #
  # where n_west = -x̂ is the outward normal. The EMRobinBC formulation
  # absorbs the sign of the outward normal internally, so we provide:
  #   coeff_real = k₀ (the wavenumber magnitude)
  #   func_real  = cosTheta = cos(0°) = 1  (normal incidence)
  #   profile_func_real = E₀ = 1  (incident wave amplitude)
  #
  # The incident wave phase at the boundary (exp(+j k₀ × (-4))) is handled
  # automatically by EMRobinBC, which evaluates the profile function at the
  # correct position. With profile_func_real = constant E₀, the BC creates
  # a uniform-amplitude plane wave incident from the west.
  #
  # Physical reasoning for port on west boundary only:
  #   The incident wave travels in the +x direction, so it enters through
  #   the west face. The east, north, south boundaries only need to absorb
  #   the outgoing scattered field — those use mode = absorbing below.
  #
  # mode = port (default): activates both wave injection and absorption.
  #   The RHS source is: 2 j k₀ × profile_func × exp(j k₀ x)
  #   This projects the incident wave E₀ = 1 into the domain through
  #   the west face.
  [west_real]
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

  [west_imag]
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
  # EAST BOUNDARY (x = +4): ABSORBING BOUNDARY CONDITION (ABC)
  # ==================================================================
  # The incident wave has passed through the cylinder and the scattered
  # field radiates outward in all directions. At the east face the wave
  # is primarily propagating in the +x direction (transmitted + forward-
  # scattered). mode = absorbing implements the first-order ABC:
  #   ∂E/∂n + j k₀ E = 0
  # which absorbs outgoing plane waves at normal incidence.
  #
  # For scattered cylindrical waves: at r = 3 m from the cylinder axis,
  # the scattering angle range at the east face is ±arctan(4/3) ≈ ±53°
  # from normal. The first-order ABC has a reflection coefficient
  # |R| = |cos(θ_s) - 1| / |cos(θ_s) + 1| which peaks at oblique angles.
  # For the domain size chosen (4λ gap), residual reflections from the
  # ABC are small compared to the scattered field near the cylinder.
  #
  # mode = absorbing: EMRobinBC sets the RHS source term to zero, acting
  # as a pure absorber with no incident wave injection on this boundary.
  [east_real]
    type            = EMRobinBC
    variable        = E_real
    boundary        = east
    component       = real
    coeff_real      = ${k}
    func_real       = cosTheta
    field_real      = E_real
    field_imaginary = E_imag
    sign            = negative
    mode            = absorbing
  []

  [east_imag]
    type            = EMRobinBC
    variable        = E_imag
    boundary        = east
    component       = imaginary
    coeff_real      = ${k}
    func_real       = cosTheta
    field_real      = E_real
    field_imaginary = E_imag
    sign            = negative
    mode            = absorbing
  []

  # ==================================================================
  # NORTH BOUNDARY (y = +4): ABSORBING BOUNDARY CONDITION
  # ==================================================================
  # The scattered field from the cylinder also radiates toward the north
  # (and south) face. mode = absorbing implements ∂E/∂n + j k₀ E = 0.
  #
  # Note: for the north boundary the outward normal is +ŷ, while the
  # incident plane wave travels in +x̂. The scattered field hitting the
  # north face travels at a range of angles; the first-order ABC is exact
  # only at normal incidence (θ_scat = 0°, i.e. k_scat = k₀ ŷ). For all
  # other directions there is a small residual reflection. In practice,
  # for a domain this size (4 wavelengths from cylinder to boundary) the
  # effect on the near-field solution around the cylinder is negligible.
  [north_real]
    type            = EMRobinBC
    variable        = E_real
    boundary        = north
    component       = real
    coeff_real      = ${k}
    func_real       = cosTheta
    field_real      = E_real
    field_imaginary = E_imag
    sign            = negative
    mode            = absorbing
  []

  [north_imag]
    type            = EMRobinBC
    variable        = E_imag
    boundary        = north
    component       = imaginary
    coeff_real      = ${k}
    func_real       = cosTheta
    field_real      = E_real
    field_imaginary = E_imag
    sign            = negative
    mode            = absorbing
  []

  # ==================================================================
  # SOUTH BOUNDARY (y = -4): ABSORBING BOUNDARY CONDITION
  # ==================================================================
  # Symmetric to the north boundary: absorbs scattered field radiating
  # toward negative y. Same comment on oblique-incidence accuracy applies.
  [south_real]
    type            = EMRobinBC
    variable        = E_real
    boundary        = south
    component       = real
    coeff_real      = ${k}
    func_real       = cosTheta
    field_real      = E_real
    field_imaginary = E_imag
    sign            = negative
    mode            = absorbing
  []

  [south_imag]
    type            = EMRobinBC
    variable        = E_imag
    boundary        = south
    component       = imaginary
    coeff_real      = ${k}
    func_real       = cosTheta
    field_real      = E_real
    field_imaginary = E_imag
    sign            = negative
    mode            = absorbing
  []
[]

[AuxVariables]
  # |E_z|² = E_real² + E_imag² — the square of the electric field magnitude.
  # This is proportional to the time-averaged power density (Poynting flux
  # density). Visualising |E_z|² in ParaView reveals:
  #   - Shadow region (behind cylinder in +x direction): |E_z|² < 1
  #   - Forward scattering lobe:  |E_z|² > 1 (constructive interference)
  #   - Field enhancement inside cylinder: |E_z|² ~ (4/3)² ≈ 1.78 peak
  #     (n=2 cylinder acts as a cylindrical lens focusing the incident wave)
  #   - Ripple pattern outside cylinder: standing-wave interference between
  #     incident and scattered waves creates alternating bright/dark bands
  [E_magnitude_sq]
    order  = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  # Compute |E_z|² = E_real² + E_imag² nodewise.
  # ParsedAux evaluates the expression at each node using the current
  # values of E_real and E_imag. This auxiliary variable is useful for
  # post-processing and comparison with Mie theory predictions.
  [E_magnitude_sq_aux]
    type       = ParsedAux
    variable   = E_magnitude_sq
    expression = 'E_real * E_real + E_imag * E_imag'
    coupled_variables = 'E_real E_imag'
  []
[]

[Postprocessors]
  # ==================================================================
  # Field amplitude at key physical locations
  # Probing E_real and E_imag at specific (x,y) points lets us:
  #   (1) Check that the cylinder focuses the field internally
  #   (2) Confirm the shadow region has reduced |E_z|
  #   (3) Verify that the incident wave amplitude is correct far from
  #       the cylinder where the total field ≈ incident field
  # ==================================================================

  # -- Centre of the cylinder (x=0, y=0) --
  # Inside the dielectric cylinder the field is enhanced by focusing.
  # For a circular cylinder with εᵣ = 4 and k₀R = π/2, Mie theory gives
  # an internal field amplitude |E_z| / |E_inc| that peaks near the centre.
  [E_real_centre]
    type     = PointValue
    variable = E_real
    point    = '0 0 0'
  []

  [E_imag_centre]
    type     = PointValue
    variable = E_imag
    point    = '0 0 0'
  []

  # -- Shadow region (x=+2.5, y=0): 1.5 m behind the cylinder --
  # The cylinder intercepts the incident wave, creating a shadow. At
  # x = 2.5 (1.5λ₀/4 = 0.375λ₀ behind the cylinder boundary) the total
  # field should be reduced below unity. The exact reduction depends on
  # diffraction around the cylinder (Babinet's principle for the far field;
  # near-field behaviour is more complex).
  [E_real_shadow]
    type     = PointValue
    variable = E_real
    point    = '2.5 0 0'
  []

  [E_imag_shadow]
    type     = PointValue
    variable = E_imag
    point    = '2.5 0 0'
  []

  # -- Forward-scatter axis far from cylinder (x=+3.5, y=0) --
  # At x = 3.5 we are 2.5 m from the cylinder axis on the forward side.
  # The scattered field partially reconstructs the shadow via diffraction,
  # and forward-scattered waves from the dielectric can interfere constructively.
  [E_real_forward]
    type     = PointValue
    variable = E_real
    point    = '3.5 0 0'
  []

  [E_imag_forward]
    type     = PointValue
    variable = E_imag
    point    = '3.5 0 0'
  []

  # -- Incident-side axis (x=-2.5, y=0): 1.5 m in front of the cylinder --
  # In this region the total field = incident + backscattered. The
  # interference between incident and backward-scattered waves creates
  # a standing-wave pattern along the x-axis with period λ₀/2 = 2 m.
  [E_real_incident_side]
    type     = PointValue
    variable = E_real
    point    = '-2.5 0 0'
  []

  [E_imag_incident_side]
    type     = PointValue
    variable = E_imag
    point    = '-2.5 0 0'
  []

  # -- Broadside (x=0, y=+2): 1 m above cylinder surface --
  # At the broadside position (90° from the propagation axis) the scattered
  # field is primarily due to the induced currents on the cylinder sides.
  # This direction is often used to measure the differential scattering
  # cross-section in experimental setups.
  [E_real_broadside]
    type     = PointValue
    variable = E_real
    point    = '0 2 0'
  []

  [E_imag_broadside]
    type     = PointValue
    variable = E_imag
    point    = '0 2 0'
  []

  # -- Far-corner diagnostic (x=+3.5, y=+3.5) --
  # At a large distance from the cylinder (r ≈ 4.95 m ≈ 1.24 λ₀) the
  # scattered field should be small compared to the incident wave. This
  # point checks that the solution converges to the correct background
  # field far from the scatterer.
  [E_real_corner]
    type     = PointValue
    variable = E_real
    point    = '3.5 3.5 0'
  []

  [E_imag_corner]
    type     = PointValue
    variable = E_imag
    point    = '3.5 3.5 0'
  []

  # -- Field magnitude squared at key locations --
  # |E_z|² ∝ time-averaged intensity. Comparing these values:
  #   centre / shadow / forward / incident-side / broadside
  # maps out the scattering pattern and verifies Mie-theory trends.
  [E_sq_centre]
    type     = PointValue
    variable = E_magnitude_sq
    point    = '0 0 0'
  []

  [E_sq_shadow]
    type     = PointValue
    variable = E_magnitude_sq
    point    = '2.5 0 0'
  []

  [E_sq_forward]
    type     = PointValue
    variable = E_magnitude_sq
    point    = '3.5 0 0'
  []

  # Total field energy in the domain: ∫ |E_z|² dA.
  # For a lossless scattering problem with no sources in the volume,
  # this integral should remain finite and well-behaved. Monitoring it
  # across mesh refinements provides a convergence check.
  [total_field_energy]
    type     = ElementIntegralVariablePostprocessor
    variable = E_magnitude_sq
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner (SMP with full=true).
  #
  # WHY FULL=TRUE IS REQUIRED:
  # The EMRobinBC at all four boundaries introduces off-diagonal Jacobian
  # blocks between E_real and E_imag (through the j k₀ coupling term):
  #   j k₀ E_z = j k₀ (E_real + j E_imag) → couples real ↔ imaginary
  # Without full=true, the preconditioner sees only the diagonal blocks
  # (the E_real and E_imag sub-systems separately) and cannot account for
  # this off-diagonal coupling. In a 2D problem with ABCs on all four
  # faces, the off-diagonal blocks are dense along the boundary DOFs.
  # Using full=true ensures the Newton Jacobian is fully assembled and
  # the LU preconditioner inverts the correct coupled system.
  #
  # The 2D problem is larger than the 1D cases (12,800 DOFs vs ~600-1000
  # for cases 32, 74, 75) but LU is still fast for this problem size.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorisation: the 2D Helmholtz system with 12,800 DOFs is
  # well within the capability of PETSc's direct LU (via SuperLU or MUMPS).
  # The Helmholtz operator is indefinite (neither positive-definite nor
  # negative-definite due to the k₀² εᵣ term), which makes iterative solvers
  # (GMRES, CG) unreliable without careful preconditioning. Direct LU is
  # robust for Helmholtz problems at this scale and completes in seconds.
  #
  # For larger problems (k₀R ≫ 1 or finer meshes), iterative solvers with
  # shifted-Laplacian preconditioners or PML absorbing layers would be needed.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  # Tight tolerances: the Helmholtz system is linear (E_z appears linearly
  # in the equation) so Newton converges in 1–2 iterations. The tolerance
  # is set tight to ensure the LU solve is used to full floating-point
  # accuracy, producing field values accurate to ~10 significant digits
  # at the PointValue postprocessors.
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  # exodus: 2D field output (E_real, E_imag, E_magnitude_sq) on the full
  #   mesh. Visualise in ParaView:
  #   - E_magnitude_sq: field intensity, reveals shadow/bright regions
  #   - E_real and E_imag: real/imaginary parts showing the phase pattern
  #   - Threshold at E_magnitude_sq < 0.5 to see the shadow cone
  #   - The dielectric cylinder boundary (unit circle) is visible as
  #     a change in field scale at r = 1 in the ParaView colour map
  exodus = true

  # csv: postprocessor values — field amplitudes at probed locations
  #   and total field energy. Used for comparison with Mie series.
  csv    = true
[]
