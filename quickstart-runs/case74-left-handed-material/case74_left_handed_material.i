# ============================================================
# Case 74: EM Wave in a Left-Handed Material (LHM) Slab
# Veselago, Soviet Physics Uspekhi 10 (1968) 509–514
#
# A plane wave at f = 20 MHz propagates from a vacuum port
# (x = 120) through vacuum, then through an LHM slab (40–80),
# then through vacuum again, and finally reflects off a PEC
# wall at x = 0.
#
# Physical layout:
#   x = 0       x = 40      x = 80      x = 120
#   |<-- vac -->|<-- LHM -->|<-- vac -->|
#   PEC         ε=-1.5     ε=1         port (Robin BC)
#   (E=0)       μ=-1.5      μ=1         E_inc injected
#
# Left-handed materials (LHM):
# ---------------------------------------------------------------
# A medium with BOTH ε_r < 0 and μ_r < 0 simultaneously supports
# electromagnetic propagation, but with a reversed phase velocity.
# The refractive index n = -sqrt(ε_r μ_r) is NEGATIVE:
#
#   n_LHM = -sqrt((-1.5)×(-1.5)) = -1.5
#
# while the energy (group) velocity remains forward (into the slab).
# This is the hallmark of left-handed or double-negative media,
# first proposed theoretically by Veselago (1968) and realised
# experimentally via split-ring resonator metamaterials (Smith et al.,
# Science 292, 2001).
#
# Governing equation (frequency domain, 1D, normal incidence):
# ---------------------------------------------------------------
# For anisotropic media with spatially varying μ_r(x) and ε_r(x),
# the 1-D Helmholtz equation for the transverse electric field is:
#
#   d/dx [ (1/μ_r) dE/dx ] + k₀² ε_r(x) E = 0
#
# where k₀ = 2πf/c = 0.41888 rad/m.
#
# This is NOT the same as the scalar Laplacian form d²E/dx² + k₀² n² E = 0,
# because the interface condition from the 1/μ_r coefficient differs:
#
#   Continuity: E continuous,  (1/μ_r) dE/dx continuous
#
# For vacuum→LHM:   (1/1)×dE/dx|vac = (1/(-1.5))×dE/dx|LHM
#   ⟹  dE/dx|LHM = -1.5 × dE/dx|vac   (SIGN REVERSAL at interface)
#
# This sign reversal is what distinguishes n = -1.5 from n = +1.5.
# The bulk equation d²E/dx² + k₀² n² E = 0 is identical for both signs
# of n — only the interface condition differs.
#
# MOOSE enforces this interface condition automatically through the
# weak form of ADMatDiffusion: the diffusivity jump in (1/μ_r) is
# naturally built into the test-function formulation, producing the
# correct boundary flux matching at the vacuum–LHM interfaces.
#
# Weak form and residual contributions:
# ---------------------------------------------------------------
# Expanding the strong form:
#   (1/μ_r) d²E/dx² + k₀² ε_r E = 0
# rearranged (dividing through by +1 is valid even when 1/μ_r < 0):
#   −d/dx((1/μ_r) dE/dx) − k₀² ε_r E = 0  ← weak form source
#
# MOOSE residual contributions:
#   ADMatDiffusion(diffusivity = 1/μ_r):  adds  +∫(1/μ_r)∇E·∇v dV
#     strong form equivalent: −d/dx((1/μ_r) dE/dx)
#
#   ADMatReaction(reaction_rate = k₀² ε_r):  adds  −k₀² ε_r × ∫ E v dV
#     strong form equivalent: −k₀² ε_r × E
#
# Total residual strong-form:  −d/dx((1/μ_r) dE/dx) − k₀² ε_r E = 0
# which is exactly the rearranged governing equation. ✓
#
# Numerical values:
#   Vacuum (x<40 and x>80):  1/μ_r = 1.0,     k₀² ε_r = 0.17546
#   LHM slab (40≤x≤80):      1/μ_r = -0.66667, k₀² ε_r = -0.26320
#
# Note that the LHM reaction coefficient k₀² ε_r is NEGATIVE.
# ADMatReaction adds −(−0.26320)×E = +0.26320×E, which combined with
# the negative-diffusivity ADMatDiffusion term gives a well-posed system.
#
# Real/imaginary splitting (lossless case):
# ---------------------------------------------------------------
# ε_r and μ_r are purely real (lossless LHM), so the complex Helmholtz
# equation splits into two decoupled real-valued equations:
#
#   d/dx[(1/μ_r) dE_real/dx] + k₀² ε_r E_real = 0
#   d/dx[(1/μ_r) dE_imag/dx] + k₀² ε_r E_imag = 0
#
# The two components couple only through the Robin port BC.
#
# Robin port boundary condition (right boundary, x = 120):
# ---------------------------------------------------------------
# EMRobinBC (Jin, "FEM in Electromagnetics", 3rd Ed., Eq. 9.60):
#
#   ∂E/∂n + j k₀ E = 2 j k₀ E0   at x = L
#
# This simultaneously injects the incident wave and absorbs reflections.
#
# PEC boundary condition (left boundary, x = 0):
# ---------------------------------------------------------------
# DirichletBC: E_real = E_imag = 0.
# The PEC reflection interacts with the LHM slab to produce a
# standing wave whose phase structure reveals the negative-n effect.
#
# Comparison with RHM (positive-n) slab:
# ---------------------------------------------------------------
# For an RHM slab (ε_r = +2.25, μ_r = +1), n = +1.5 and n² = 2.25.
# That case uses plain Diffusion + ADMatReaction (case 32 pattern).
# The LHM (ε_r = -1.5, μ_r = -1.5) has n² = 2.25 too, yet the field
# profile DIFFERS because of the interface condition sign reversal.
# Running both cases and overlaying E_real(x) demonstrates this.
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables
# -----------------------------------------------------------
# k0 = 2π × 20 MHz / c = 2π × 20e6 / 3e8 = 0.41888 rad/m
# L  = 120 m = 8 × λ₀ (λ₀ = 15 m at 20 MHz)
# The LHM slab occupies 40 ≤ x ≤ 80 (thickness = 40 m ≈ 2.67 λ₀)
k0    = 0.41887902047863906   # free-space wavenumber [rad/m]
L     = 120                   # total domain length [m]
E0    = 1                     # incident field amplitude [V/m]
theta = 0                     # incidence angle [degrees] — normal incidence

[Mesh]
  # 1D domain from x = 0 (PEC wall) to x = L (port).
  # 600 elements gives ~5.33 elements per free-space wavelength (λ₀ = 15 m)
  # and ~8 per LHM wavelength (λ_LHM = λ₀/|n| = 15/1.5 = 10 m).
  # FIRST-order Lagrange resolves both the fast and slow oscillations.
  [domain]
    type = GeneratedMeshGenerator
    dim  = 1
    nx   = 600
    xmin = 0
    xmax = ${L}
  []

  # Rename auto-generated labels to physically meaningful names.
  # 'metal' → x = 0:   Perfect Electric Conductor (PEC), E = 0
  # 'port'  → x = L:   open port — injects incident wave, absorbs reflections
  [rename]
    type         = RenameBoundaryGenerator
    input        = domain
    old_boundary = 'left right'
    new_boundary = 'metal port'
  []
[]

[Variables]
  # E_real — real part of the complex electric field phasor E(x) = E_real + j E_imag.
  # Lagrange elements enforce C⁰ continuity across all material interfaces,
  # implementing the physical condition that tangential E is continuous.
  [E_real]
    order  = FIRST
    family = LAGRANGE
  []

  # E_imag — imaginary part of the complex phasor.
  # Decoupled from E_real in the bulk (lossless) but coupled at the port.
  [E_imag]
    order  = FIRST
    family = LAGRANGE
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # inv_mu_fn: diffusion coefficient = 1/μ_r(x)
  #
  # This is the key function that encodes the LHM physics.
  # In the vacuum regions, 1/μ_r = 1/1 = 1.0 (normal diffusion).
  # In the LHM slab, 1/μ_r = 1/(-1.5) = -2/3 ≈ -0.66667.
  # The NEGATIVE diffusivity in the slab makes ADMatDiffusion subtract
  # the gradient term there, which physically corresponds to the reversed
  # relationship between the field gradient and the displacement current.
  #
  # The interface condition (1/μ_r)dE/dx = continuous is automatically
  # enforced by the weak form: no special interface kernels are needed.
  # ------------------------------------------------------------------
  [inv_mu_fn]
    type       = ParsedFunction
    expression = 'if(x>=40 & x<=80, -0.66667, 1.0)'
  []

  # ------------------------------------------------------------------
  # k0sq_eps_fn: reaction coefficient = k₀² ε_r(x)
  #
  # Used by ADMatReaction to implement the −k₀² ε_r E term.
  # Recall ADMatReaction residual sign: adds −rate × E × test.
  # So reaction_rate = +k₀² ε_r gives strong-form contribution −k₀² ε_r E,
  # which combines with the diffusion term to give the Helmholtz equation.
  #
  # Numerical values:
  #   Vacuum: k₀² × (+1)    = 0.41888² × 1    = +0.17546 m⁻²
  #   LHM:    k₀² × (-1.5)  = 0.41888² × (-1.5) = -0.26320 m⁻²
  #
  # The NEGATIVE reaction rate in the LHM means ADMatReaction adds
  # −(−0.26320)×E = +0.26320×E to the residual there. This, together
  # with the negative-diffusivity diffusion term, enforces the wave
  # equation inside the double-negative medium.
  # ------------------------------------------------------------------
  [k0sq_eps_fn]
    type       = ParsedFunction
    expression = 'if(x>=40 & x<=80, -0.26320, 0.17546)'
  []

  # cosTheta = cos(0°) = 1 for normal incidence.
  # EMRobinBC uses this to compute k_eff = k₀ cos(θ) for oblique incidence.
  [cosTheta]
    type       = ParsedFunction
    expression = 'cos(${theta})'
  []
[]

[Materials]
  # Wrap inv_mu_fn as an AD material property for ADMatDiffusion.
  # ADGenericFunctionMaterial evaluates the ParsedFunction at each quadrature
  # point and returns an ADReal, enabling automatic differentiation.
  # The NEGATIVE value in the LHM slab is handled transparently by the kernel.
  [inv_mu_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'inv_mu_r'
    prop_values = 'inv_mu_fn'
  []

  # Wrap k0sq_eps_fn as an AD material property for ADMatReaction.
  [k0sq_eps_mat]
    type        = ADGenericFunctionMaterial
    prop_names  = 'k0sq_eps_r'
    prop_values = 'k0sq_eps_fn'
  []
[]

[Kernels]
  # ==================================================================
  # Equation for E_real:
  #   d/dx[(1/μ_r) dE_real/dx] + k₀² ε_r(x) E_real = 0
  # Weak form:
  #   ∫ (1/μ_r)(dE_real/dx)(dv/dx) dx − k₀² ε_r ∫ E_real v dx = 0
  # ==================================================================

  # ADMatDiffusion kernel: adds +∫(1/μ_r)∇E_real·∇v dV.
  # When diffusivity = inv_mu_r = 1/μ_r, this is the correct discretisation
  # of the 1/μ_r-weighted Laplacian from the Maxwell curl equations.
  # In the LHM slab (diffusivity = -0.66667) the sign of this integral flips,
  # which together with the negative reaction rate keeps the system balanced.
  [diffusion_real]
    type        = ADMatDiffusion
    variable    = E_real
    diffusivity = inv_mu_r
  []

  # ADMatReaction: adds −k₀² ε_r × ∫ E_real v dx.
  # reaction_rate = k₀² ε_r (which is −0.26320 in the LHM slab).
  # ADMatReaction residual sign: −rate × E × test = −(−0.26320) × E = +0.26320 E.
  # Combined with the negative diffusion term in the slab, the total residual
  # reproduces the Helmholtz equation with the correct propagation constant.
  [field_real]
    type          = ADMatReaction
    variable      = E_real
    reaction_rate = k0sq_eps_r
  []

  # NOTE: No cross-coupling kernel for E_real.
  # The medium is lossless: Im(ε_r) = Im(μ_r) = 0.
  # Cross-coupling terms (proportional to ε_r'' or μ_r'') vanish.

  # ==================================================================
  # Equation for E_imag:
  #   d/dx[(1/μ_r) dE_imag/dx] + k₀² ε_r(x) E_imag = 0
  # Identical structure to E_real — lossless medium decouples the two.
  # ==================================================================

  [diffusion_imag]
    type        = ADMatDiffusion
    variable    = E_imag
    diffusivity = inv_mu_r
  []

  [field_imag]
    type          = ADMatReaction
    variable      = E_imag
    reaction_rate = k0sq_eps_r
  []

  # NOTE: No cross-coupling kernel for E_imag either.
[]

[BCs]
  # ------------------------------------------------------------------
  # PEC wall at x = 0 (boundary 'metal')
  # ------------------------------------------------------------------
  # Perfect Electric Conductor: tangential E = 0.
  # In phasor form:  E(0) = 0  ⟹  E_real(0) = 0  AND  E_imag(0) = 0.
  # The PEC reflection, combined with the LHM slab, creates a standing
  # wave whose spatial phase pattern reveals the negative-index effect.
  [pec_real]
    type     = DirichletBC
    variable = E_real
    boundary = metal
    value    = 0
  []

  [pec_imag]
    type     = DirichletBC
    variable = E_imag
    boundary = metal
    value    = 0
  []

  # ------------------------------------------------------------------
  # Robin port condition at x = L (boundary 'port')
  # ------------------------------------------------------------------
  # EMRobinBC (Jin, "FEM in Electromagnetics", Eq. 9.60):
  #
  #   ∂E/∂n + j k₀ cos(θ) E = 2 j k₀ cos(θ) E0 exp(j k₀ cos(θ) x)
  #
  # For normal incidence (θ = 0, cos(θ) = 1, E0 = 1 V/m):
  #   ∂E/∂x + j k₀ E = 2 j k₀ exp(j k₀ x)  at x = L
  #
  # The Robin BC couples E_real and E_imag at the port boundary only.
  # It simultaneously:
  #   (1) Injects the plane-wave incident field  E_inc = E0 exp(j k₀ x)
  #   (2) Absorbs the outgoing reflected wave transparently
  #
  # Parameters match case 32 exactly (same frequency, same port structure).
  [port_real]
    type              = EMRobinBC
    variable          = E_real
    boundary          = port
    component         = real
    coeff_real        = ${k0}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []

  [port_imag]
    type              = EMRobinBC
    variable          = E_imag
    boundary          = port
    component         = imaginary
    coeff_real        = ${k0}
    func_real         = cosTheta
    profile_func_real = ${E0}
    field_real        = E_real
    field_imaginary   = E_imag
    sign              = negative
  []
[]

[Postprocessors]
  # Power reflection coefficient |R|² at the port.
  # ReflectionCoefficient reconstructs the reflected amplitude:
  #   E_total(L) = E_inc(L) + R_complex × E_ref(L)
  #   |R|² = |R_complex|²
  # For the LHM slab-on-PEC geometry, the phase reversal at both
  # vacuum–LHM interfaces and the PEC reflection interact to produce
  # a Fabry-Pérot pattern that differs from the RHM n=+1.5 case.
  [reflection_coefficient]
    type                     = ReflectionCoefficient
    k                        = ${k0}
    theta                    = ${theta}
    length                   = ${L}
    incoming_field_magnitude = ${E0}
    field_real               = E_real
    field_imag               = E_imag
    boundary                 = port
  []

  # Field values at key x-locations for diagnostics and comparison.
  # Vacuum left of slab (x = 20): standing wave in left vacuum region.
  [E_real_x20]
    type     = PointValue
    variable = E_real
    point    = '20 0 0'
  []

  [E_imag_x20]
    type     = PointValue
    variable = E_imag
    point    = '20 0 0'
  []

  # LHM slab centre (x = 60): phase reversal is most visible here.
  [E_real_x60]
    type     = PointValue
    variable = E_real
    point    = '60 0 0'
  []

  [E_imag_x60]
    type     = PointValue
    variable = E_imag
    point    = '60 0 0'
  []

  # Vacuum right of slab (x = 100): incident + reflected wave superposition.
  [E_real_x100]
    type     = PointValue
    variable = E_real
    point    = '100 0 0'
  []

  [E_imag_x100]
    type     = PointValue
    variable = E_imag
    point    = '100 0 0'
  []

  # Port (x = L = 120): total field at the injection boundary.
  [E_real_at_port]
    type     = PointValue
    variable = E_real
    point    = '${L} 0 0'
  []

  [E_imag_at_port]
    type     = PointValue
    variable = E_imag
    point    = '${L} 0 0'
  []
[]

[Preconditioning]
  # Full Single Matrix Preconditioner: required for the coupled (E_real, E_imag)
  # system. The Robin BC introduces off-diagonal Jacobian blocks between the
  # two field components. Without full=true the preconditioner omits these
  # couplings and convergence degrades near the port.
  #
  # Note on negative diffusivity: the system is non-coercive in the LHM slab
  # (the bilinear form is not positive-definite there). This is physically
  # correct for a Helmholtz-type wave problem. Direct LU handles it exactly;
  # iterative solvers would need careful preconditioning for large problems.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Steady
  solve_type = NEWTON

  # Direct LU factorization: the 1D system (1200 DOFs for 600 elements with
  # two variables) is well within the range where direct LU is fast and robust.
  # The non-coercive (indefinite) system arising from the LHM negative diffusivity
  # and negative reaction rate makes iterative solvers unreliable near resonance;
  # LU avoids all convergence issues.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true   # Full 1D field (E_real, E_imag) along x — visualise in ParaView or matplotlib
  csv    = true   # reflection_coefficient and point values
[]
