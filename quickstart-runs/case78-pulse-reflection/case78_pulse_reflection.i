# ============================================================
# Case 78: Time-Domain EM Pulse Reflection from a Dielectric Slab (2D)
# MIT 6.635 Advanced Electromagnetism — Lecture 7
#
# A narrow-band Gaussian EM pulse propagates in the +x direction,
# hits a dielectric slab (eps_r = 4, n = 2), and splits into a
# reflected pulse traveling back in -x and a transmitted pulse
# continuing in +x at reduced speed c/n = c/2.
#
# Physical layout (x-axis):
#   x = 0            x = x1         x = x2        x = Lx
#   |<---- vacuum --->|<-- slab (n=2) -->|<-- vacuum --->|
#   absorbing BC     interface 1    interface 2    source + absorbing BC
#
# The slab occupies x in [x1, x2] = [0.3, 0.5] m (slab_thick = 0.2 m).
# Domain is Lx = 0.8 m wide, Ly = 0.08 m tall (effectively 1D physics).
# A QUAD9 mesh with ~1800 elements is used.
#
# ============================================================
# Governing Equation
# ============================================================
#
# Maxwell's curl equations in a non-magnetic (mu_r = 1) medium:
#
#   curl(curl(E)) + (eps_r / c^2) * d^2E/dt^2 = 0
#
# In component form for TE polarization (E_y only, no variation in y):
#
#   -d^2 E_y / dx^2 + (eps_r / c^2) * d^2 E_y / dt^2 = 0
#
# This is the 1D (or 2D TE) wave equation for the y-component.
# The curl-curl form is general and handles both the propagating
# vacuum regions and the dielectric slab naturally.
#
# ============================================================
# Analytic Signal Representation
# ============================================================
#
# Instead of a single real-valued time-domain variable, MOOSE's
# electromagnetics module uses a two-variable analytic signal approach:
#
#   E(x, y, t) = E_real(x, y, t) + j * E_imag(x, y, t)
#
# where E_real and E_imag are the in-phase (cosine) and
# quadrature (sine) components of the complex analytic signal.
# The physical electric field is the real part: E_real(x, y, t).
#
# This formulation enables the VectorTransientAbsorbingBC to implement
# the correct first-order Mur absorbing boundary condition, which
# requires knowledge of both components of the analytic signal to
# correctly absorb outgoing waves without reflection.
#
# Both variables satisfy the same PDE independently in the bulk;
# they are coupled only through the absorbing BCs.
#
# ============================================================
# Gaussian Pulse Source
# ============================================================
#
# The incident pulse is injected at the right boundary (x = Lx)
# via a VectorCurlPenaltyDirichletBC with a modulated Gaussian
# envelope function:
#
#   E_inc(t) = exp(-(t - t0)^2 / (2 * sigma_t^2)) * cos(omega_0 * t)
#
# Parameters:
#   f0 = 3e9 Hz       (carrier frequency, 3 GHz)
#   omega_0 = 2*pi*f0 = 18.85e9 rad/s
#   lambda0 = c/f0 = 0.1 m     (free-space wavelength)
#   sigma_t = 4 / omega_0      (pulse width: ~4 cycles FWHM)
#   t0 = 6 * sigma_t           (pulse center: starts near t=0)
#   c = 3e8 m/s
#
# The corresponding quadrature component uses sin(omega_0 * t).
#
# ============================================================
# Fresnel Reflection/Transmission Coefficients
# ============================================================
#
# For normal incidence from vacuum (n1=1) to dielectric (n2=2):
#
#   R = (n1 - n2) / (n1 + n2) = (1 - 2) / (1 + 2) = -1/3
#   T = 2*n1 / (n1 + n2) = 2/3
#
# Reflected power: |R|^2 = 1/9  ≈ 11.1%
# Transmitted power: n2/n1 * |T|^2 = 2 * 4/9 = 8/9 ≈ 88.9%
# (Energy conservation: |R|^2 + n2/n1*|T|^2 = 1 ✓)
#
# The reflected pulse has NEGATIVE amplitude (phase shift of π).
# The transmitted pulse propagates at c/n2 = c/2, so it arrives
# at x = Lx/2 at time t ≈ 2*Lx/(2*c) — twice the vacuum travel time.
#
# ============================================================
# Time-Stepping: NewmarkBeta
# ============================================================
#
# NewmarkBeta (beta=0.25, gamma=0.5) is the implicit trapezoidal
# rule for second-order ODEs. It is unconditionally stable for
# linear problems and second-order accurate in time.
#
# For numerical accuracy, the time step should resolve the carrier:
#   dt < 1/(20 * f0) = 1/(20 * 3e9) ≈ 1.67e-11 s
# We use dt = 1e-11 s (0.01 ns) and run for 200 steps (2 ns total),
# which is long enough for the pulse to traverse the domain (~2.67 ns).
#
# ============================================================
# References
# ============================================================
#
# Griffiths, "Introduction to Electrodynamics", 4th Ed., Sec. 9.3
# Jackson, "Classical Electrodynamics", 3rd Ed., Ch. 7
# Jin, "The Finite Element Method in Electromagnetics", 3rd Ed.
# Mur, "Absorbing boundary conditions...", IEEE Trans. EMC 23 (1981)
# ============================================================

# -----------------------------------------------------------
# HIT top-level variables
# -----------------------------------------------------------
c       = 3e8          # speed of light in vacuum [m/s]
# f0 = 3e9 Hz (carrier frequency), lambda0 = c/f0 = 0.1 m
omega0  = 18849555921.538759   # 2*pi*3e9 [rad/s]

# Domain geometry
Lx      = 0.8          # domain length in x [m] = 8 wavelengths
Ly      = 0.08         # domain height in y [m] = 0.8 wavelengths (thin strip)
nx      = 200          # elements in x  (0.004 m per element ≈ 25 per wavelength)
ny      = 5            # elements in y  (thin strip, minimal transverse resolution needed)

# Dielectric slab
x1      = 0.3          # slab left interface [m]
x2      = 0.5          # slab right interface [m]
eps_slab = 4.0         # relative permittivity of slab (n = sqrt(4) = 2)

# Gaussian pulse parameters
sigma_t = 2.122e-10    # pulse width ≈ 4/omega0 [s] (~4 carrier cycles)
t0      = 1.273e-9     # pulse center = 6*sigma_t [s] (starts near zero at t=0)
# Admittance of free space: Y0 = 1/(mu_0*c) = eps_0*c
# VectorTransientAbsorbingBC default = 1/(mu_0*c) = 2653.07 S/m -- free space admittance

# Time stepping
dt      = 1e-11        # time step [s] = 0.01 ns (resolves carrier at 20 pts/cycle)
num_steps = 200        # total time = 2 ns (enough for pulse to cross 0.6 m at c)

[Mesh]
  # 2D rectangular mesh: x in [0, Lx], y in [0, Ly].
  # We use QUAD9 (9-node serendipity quads) as required for NEDELEC_ONE edge elements.
  # NEDELEC_ONE on QUAD9 gives the correct curl-conforming interpolation for E-fields.
  #
  # Physical regions:
  #   "slab"    — dielectric with eps_r = 4 (x in [x1, x2])
  #   "vacuum"  — free space with eps_r = 1 (x < x1 or x > x2)
  # Boundaries:
  #   "left"    — x = 0     (absorbing BC for E_real + E_imag)
  #   "right"   — x = Lx    (pulse source + absorbing BC)
  #   "top"     — y = Ly    (PEC: tangential E = 0, enforced by symmetry)
  #   "bottom"  — y = 0     (PEC: tangential E = 0, enforced by symmetry)
  [full_domain]
    type   = GeneratedMeshGenerator
    dim    = 2
    elem_type = QUAD9
    nx     = ${nx}
    ny     = ${ny}
    xmin   = 0
    xmax   = ${Lx}
    ymin   = 0
    ymax   = ${Ly}
  []

  # Rename boundaries to physically meaningful names
  [rename_bounds]
    type         = RenameBoundaryGenerator
    input        = full_domain
    old_boundary = 'left right bottom top'
    new_boundary = 'absorb_left source absorb_bottom absorb_top'
  []
[]

[Variables]
  # E_real — in-phase (cosine) component of the analytic signal.
  # This is the physically observable electric field.
  # NEDELEC_ONE (curl-conforming edge elements) enforces tangential
  # continuity of E across all material interfaces automatically.
  # The y-component of E (E_y) is the dominant TE-polarized component.
  [E_real]
    family = NEDELEC_ONE
    order  = FIRST
  []

  # E_imag — quadrature (sine) component of the analytic signal.
  # Required by VectorTransientAbsorbingBC to implement the Mur
  # first-order absorbing BC (which needs the complex analytic signal).
  # In the bulk, E_imag satisfies the identical wave equation as E_real.
  [E_imag]
    family = NEDELEC_ONE
    order  = FIRST
  []
[]

[Functions]
  # ------------------------------------------------------------------
  # eps_r_fn: spatially varying relative permittivity
  #
  # eps_r = 4.0 inside the slab [x1, x2], 1.0 elsewhere (vacuum).
  # This is the refractive-index-squared that determines wave speed.
  # ------------------------------------------------------------------
  [eps_r_fn]
    type       = ParsedFunction
    expression = 'if(x >= ${x1} & x <= ${x2}, ${eps_slab}, 1.0)'
  []

  # ------------------------------------------------------------------
  # coeff_fn: second-time-derivative coefficient = eps_r / c^2
  #
  # This is the coefficient multiplying d^2E/dt^2 in the wave equation:
  #   curl(curl(E)) + (eps_r/c^2) * d^2E/dt^2 = 0
  #
  # In vacuum: eps_r/c^2 = 1/(3e8)^2 = 1.111e-17 s^2/m^2
  # In slab:   eps_r/c^2 = 4/(3e8)^2 = 4.444e-17 s^2/m^2
  # ------------------------------------------------------------------
  [coeff_fn]
    type       = ParsedFunction
    expression = 'if(x >= ${x1} & x <= ${x2}, ${eps_slab}, 1.0) / (${c} * ${c})'
  []

  # ------------------------------------------------------------------
  # Gaussian pulse source functions (applied at x = Lx boundary)
  #
  # The modulated Gaussian pulse is a common test signal in FDTD/FEM:
  #   E_inc(t) = A * exp(-(t - t0)^2 / (2 * sigma_t^2)) * cos(omega_0 * t)
  #
  # A = 1.0 V/m (unit amplitude)
  # sigma_t = 4/omega_0 ≈ 2.12e-10 s  (pulse width, ~4 carrier cycles FWHM)
  # t0 = 6*sigma_t ≈ 1.27e-9 s        (delay so pulse starts near zero at t=0)
  #
  # E_y component only; E_x = 0 (TE polarization, E perpendicular to x-hat).
  #
  # The quadrature component E_imag uses sin(omega_0*t) so that the
  # analytic signal pair correctly represents the causal modulated pulse.
  # ------------------------------------------------------------------
  [pulse_real_fn]
    type       = ParsedVectorFunction
    # E_x = 0 (propagation direction), E_y = Gaussian envelope × cos
    expression_x = '0'
    expression_y = 'exp(-(t - ${t0})^2 / (2 * ${sigma_t}^2)) * cos(${omega0} * t)'
  []

  [pulse_imag_fn]
    type       = ParsedVectorFunction
    # Quadrature component: same envelope, but sin instead of cos
    expression_x = '0'
    expression_y = 'exp(-(t - ${t0})^2 / (2 * ${sigma_t}^2)) * sin(${omega0} * t)'
  []
[]

[Kernels]
  # ==================================================================
  # Wave equation for E_real:
  #   curl(curl(E_real)) + (eps_r/c^2) * d^2E_real/dt^2 = 0
  #
  # Weak form (integrate against test function v, integrate curl by parts):
  #   ∫ (curl E_real)·(curl v) dV + ∫ (eps_r/c^2) (d^2E_real/dt^2)·v dV
  #     - ∮ (n × curl E_real)·v dS = 0
  #
  # The boundary integrals are handled by the BCs (absorbing + source).
  # ==================================================================

  # CurlCurlField: implements ∫ (curl E)·(curl v) dV
  # This is the spatial part of the curl-curl operator.
  # coeff = 1.0 (vacuum permeability mu_r = 1 everywhere).
  [curl_curl_real]
    type     = CurlCurlField
    variable = E_real
    coeff    = 1.0
  []

  # VectorSecondTimeDerivative: implements ∫ coefficient * (d^2E/dt^2)·v dV
  # coefficient function = eps_r(x)/c^2 (spatially varying).
  # Uses NewmarkBeta second-derivative estimate d^2E/dt^2.
  [time_deriv_real]
    type        = VectorSecondTimeDerivative
    variable    = E_real
    coefficient = coeff_fn
  []

  # ==================================================================
  # Wave equation for E_imag (identical structure to E_real):
  #   curl(curl(E_imag)) + (eps_r/c^2) * d^2E_imag/dt^2 = 0
  #
  # In the bulk both components decouple; coupling occurs only at the
  # absorbing BCs through the complex admittance condition.
  # ==================================================================

  [curl_curl_imag]
    type     = CurlCurlField
    variable = E_imag
    coeff    = 1.0
  []

  [time_deriv_imag]
    type        = VectorSecondTimeDerivative
    variable    = E_imag
    coefficient = coeff_fn
  []
[]

[BCs]
  # ==================================================================
  # Gaussian pulse source at x = Lx (right boundary "source")
  # ==================================================================
  #
  # VectorCurlPenaltyDirichletBC enforces the incoming pulse in a
  # weak (penalty) sense by adding to the residual:
  #
  #   penalty * ∫ (E - E_inc) × n · (v × n) dS
  #
  # The penalty must be large enough to enforce the BC strongly, but
  # not so large as to ill-condition the system. 1e6 is standard for
  # NEDELEC_ONE problems at this scale.
  #
  # The pulse is E_y-polarized; E_x = 0 (TE wave propagating in +x).
  # Applying the source only at the right ensures a right-to-left
  # injected wave — but physically we want a left-to-right pulse.
  # Therefore we inject at "source" (right boundary) and the absorbing
  # BC at "absorb_left" prevents reflection back into the domain.
  # The pulse travels LEFT from x=Lx toward the slab at x=[x1,x2],
  # hits it, and creates reflected (going right back to x=Lx) and
  # transmitted (going left past x=x1) pulses.
  #
  # NOTE: For a left-to-right pulse, we inject at the LEFT boundary.
  # Here the source is at the RIGHT boundary "source" (x = Lx = 0.8 m).
  # The pulse propagates in -x direction, hits the slab from the right,
  # reflects rightward and transmits leftward. This is equivalent by
  # symmetry to the case stated in the problem description.
  #
  [source_real]
    type       = VectorCurlPenaltyDirichletBC
    variable   = E_real
    boundary   = source
    function   = pulse_real_fn
    penalty    = 1e6
  []

  [source_imag]
    type       = VectorCurlPenaltyDirichletBC
    variable   = E_imag
    boundary   = source
    function   = pulse_imag_fn
    penalty    = 1e6
  []

  # ==================================================================
  # First-order Mur absorbing BC at x = 0 (left boundary "absorb_left")
  # ==================================================================
  #
  # VectorTransientAbsorbingBC implements the first-order Mur absorbing BC:
  #
  #   n × curl(E) + (admittance / mu_0) * n × (n × dE/dt) = 0
  #
  # where admittance = 1/(mu_0 * c) for vacuum = free-space admittance.
  #
  # For a wave hitting this boundary from the right: the BC absorbs the
  # normally incident component with zero reflection for plane waves.
  # Oblique incidence introduces small errors (first-order Mur).
  #
  # The BC requires both E_real and E_imag (coupled_field) because the
  # absorbing condition involves the time derivative of the complex
  # analytic signal.
  #
  [absorb_left_real]
    type          = VectorTransientAbsorbingBC
    variable      = E_real
    coupled_field = E_imag
    boundary      = absorb_left
    component     = real
    # admittance defaults to 1/(mu_0*c) = free-space admittance; vacuum BC
  []

  [absorb_left_imag]
    type          = VectorTransientAbsorbingBC
    variable      = E_imag
    coupled_field = E_real
    boundary      = absorb_left
    component     = imaginary
  []

  # ==================================================================
  # Absorbing BC at top and bottom (y = 0 and y = Ly)
  # ==================================================================
  #
  # For TE polarization (E in y-direction) with a uniform source at x=Lx,
  # the field has no y-variation in the interior. The top/bottom walls
  # are periodic by symmetry; we use absorbing BCs to avoid reflections
  # from these artificial boundaries.
  #
  # Note: For E_y-polarized fields on a boundary with normal in y-hat,
  # n × E = 0 for the E_x component (tangential to y-boundary). The
  # absorbing BC still correctly handles the E_y component which is
  # tangential to the y-boundaries.
  #
  [absorb_top_real]
    type          = VectorTransientAbsorbingBC
    variable      = E_real
    coupled_field = E_imag
    boundary      = absorb_top
    component     = real
  []

  [absorb_top_imag]
    type          = VectorTransientAbsorbingBC
    variable      = E_imag
    coupled_field = E_real
    boundary      = absorb_top
    component     = imaginary
  []

  [absorb_bottom_real]
    type          = VectorTransientAbsorbingBC
    variable      = E_real
    coupled_field = E_imag
    boundary      = absorb_bottom
    component     = real
  []

  [absorb_bottom_imag]
    type          = VectorTransientAbsorbingBC
    variable      = E_imag
    coupled_field = E_real
    boundary      = absorb_bottom
    component     = imaginary
  []
[]

[AuxVariables]
  # NEDELEC_ONE is an edge-element family — its DOFs live on edges, not nodes.
  # Standard scalar postprocessors (NodalMax, ElementIntegral) cannot operate
  # directly on NEDELEC_ONE variables. We therefore extract the E_y component
  # into a MONOMIAL (element-piecewise-constant) AuxVariable at each time step,
  # then apply standard postprocessors to that scalar field.
  #
  # VectorVariableComponentAux samples the NEDELEC_ONE E_real/E_imag at
  # quadrature points and stores the y-component (component = y).

  [E_y_real]
    # y-component of E_real, stored as piecewise-constant over elements
    family = MONOMIAL
    order  = CONSTANT
  []

  [E_y_imag]
    # y-component of E_imag
    family = MONOMIAL
    order  = CONSTANT
  []

  [E_amplitude]
    # Total amplitude sqrt(E_real_y^2 + E_imag_y^2) — computed in AuxKernels
    family = MONOMIAL
    order  = CONSTANT
  []
[]

[AuxKernels]
  # Extract y-component of E_real (the physically observable TE field).
  # This runs after each time step (execute_on = TIMESTEP_END).
  [extract_Ey_real]
    type            = VectorVariableComponentAux
    variable        = E_y_real
    vector_variable = E_real
    component       = y
    execute_on      = TIMESTEP_END
  []

  # Extract y-component of E_imag (quadrature component of the analytic signal).
  [extract_Ey_imag]
    type            = VectorVariableComponentAux
    variable        = E_y_imag
    vector_variable = E_imag
    component       = y
    execute_on      = TIMESTEP_END
  []

  # Compute the instantaneous amplitude envelope = sqrt(E_y_real^2 + E_y_imag^2).
  # For a narrowband pulse, this is the amplitude of the analytic signal.
  [compute_amplitude]
    type       = ParsedAux
    variable   = E_amplitude
    expression = 'sqrt(E_y_real * E_y_real + E_y_imag * E_y_imag)'
    coupled_variables = 'E_y_real E_y_imag'
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  # ============================================================
  # Point probes of E_y_real at key x-locations along centerline
  # (y = Ly/2 = 0.04 m).  These track the time history of the
  # in-phase field component — oscillations + envelope are visible.
  #
  # NOTE: Using MONOMIAL AuxVariable E_y_real (not NEDELEC_ONE E_real)
  # so that PointValue interpolates a standard scalar field.
  # ============================================================

  # Left vacuum (x = 0.1 m) — watches the transmitted pulse arrive
  [Ey_x01]
    type     = PointValue
    variable = E_y_real
    point    = '0.1 0.04 0'
  []

  # Left vacuum (x = 0.2 m) — transmitted pulse monitoring
  [Ey_x02]
    type     = PointValue
    variable = E_y_real
    point    = '0.2 0.04 0'
  []

  # Left interface of slab (x = x1 = 0.3 m) — wave entering slab
  [Ey_x03]
    type     = PointValue
    variable = E_y_real
    point    = '0.3 0.04 0'
  []

  # Slab center (x = 0.4 m) — wave propagating at half-speed inside slab
  [Ey_slab]
    type     = PointValue
    variable = E_y_real
    point    = '0.4 0.04 0'
  []

  # Right interface of slab (x = x2 = 0.5 m) — pulse exiting slab
  [Ey_x05]
    type     = PointValue
    variable = E_y_real
    point    = '0.5 0.04 0'
  []

  # Right vacuum (x = 0.6 m) — incident + reflected superposition
  [Ey_x06]
    type     = PointValue
    variable = E_y_real
    point    = '0.6 0.04 0'
  []

  # Right vacuum (x = 0.7 m) — close to source, sees both incident and reflected
  [Ey_x07]
    type     = PointValue
    variable = E_y_real
    point    = '0.7 0.04 0'
  []

  # At source boundary (x = Lx = 0.8 m) — monitors reflected pulse returning
  [Ey_source]
    type     = PointValue
    variable = E_y_real
    point    = '0.799 0.04 0'
  []

  # ============================================================
  # Amplitude envelope probes (slow-varying, shows pulse envelope)
  # ============================================================

  [Amp_x02]
    type     = PointValue
    variable = E_amplitude
    point    = '0.2 0.04 0'
  []

  [Amp_x04]
    type     = PointValue
    variable = E_amplitude
    point    = '0.4 0.04 0'
  []

  [Amp_x06]
    type     = PointValue
    variable = E_amplitude
    point    = '0.6 0.04 0'
  []

  # ============================================================
  # Global energy-like measure: integral of E_y_real^2 over domain
  # Monitors total field energy (should be ~conserved after pulse passes)
  # ============================================================
  [E_energy]
    type               = ElementIntegralVariablePostprocessor
    variable           = E_y_real
    # This integrates E_y (scalar MONOMIAL) — a proxy for field energy density
  []

  # Peak field amplitude in the domain at each time step
  [E_peak]
    type     = ElementExtremeValue
    variable = E_amplitude
    value_type = max
  []
[]

[Preconditioning]
  # Full SMP required because VectorTransientAbsorbingBC creates
  # off-diagonal Jacobian blocks coupling E_real to E_imag (and vice versa).
  # Without full = true the preconditioner omits these and may fail to converge.
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type       = Transient
  solve_type = NEWTON

  # LU direct solver is robust for the indefinite wave-equation systems.
  # For larger problems, GMRES with hypre or ilu preconditioner is needed.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  # NewmarkBeta time integrator (beta=0.25, gamma=0.5 = implicit trapezoidal rule).
  # Unconditionally stable for linear problems; second-order accurate in time.
  # Uses dotDot() (second time derivative) required by VectorSecondTimeDerivative.
  [TimeIntegrator]
    type  = NewmarkBeta
    beta  = 0.25
    gamma = 0.5
  []

  dt        = ${dt}
  num_steps = ${num_steps}

  # Nonlinear solver tolerances. The curl-curl + NewmarkBeta system requires
  # several Newton iterations per step due to the time integrator coupling.
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  nl_max_its = 15

  # Output solution every 10 steps for efficiency
  [TimeStepper]
    type = ConstantDT
    dt   = ${dt}
  []
[]

[Outputs]
  # Full time-history in Exodus for ParaView visualization:
  # load case78_pulse_reflection_out.e and animate E_real to see
  # the pulse approaching, splitting, and separating.
  [exodus_out]
    type               = Exodus
    file_base          = case78_pulse_reflection_out
    time_step_interval = 10   # Save every 10th step (20 frames total at 200 steps)
  []

  # CSV for postprocessor time series (point values, energy, peak)
  [csv_out]
    type               = CSV
    file_base          = case78_pulse_reflection
    time_step_interval = 1    # Save every step for detailed time traces
  []
[]
