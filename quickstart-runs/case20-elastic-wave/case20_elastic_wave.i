# ============================================================
# Case 20: Elastic Wave Propagation in a Bar
# Dynamic solid mechanics — Newmark-beta time integration
#
# A thin 2D strip (100 x 5 elements) approximates a 1D bar.
# A short pressure pulse is applied to the left face at t=0.
# The compressive stress wave travels rightward at the
# longitudinal wave speed c = sqrt(E/rho), reflects off the
# free right end (sign reversal: compression -> tension),
# and travels back toward the left.
#
# Domain:    x in [0, 10] m,  y in [0, 0.5] m
# Material:  Steel — E = 200 GPa,  nu = 0.0,  rho = 7800 kg/m^3
# Wave speed: c = sqrt(E/rho) = sqrt(200e9 / 7800) ~ 5064 m/s
# Transit time (one way): L/c = 10 / 5064 ~ 0.001975 s
#
# Newmark-beta parameters (unconditionally stable implicit scheme):
#   beta = 0.25, gamma = 0.5  (trapezoidal rule — no numerical damping)
#
# Time step: dt = 2e-5 s
#   CFL check: c * dt / dx = 5064 * 2e-5 / 0.1 ~ 1.01
#   (implicit Newmark can run at CFL > 1; this is fine)
#
# Run with:
#   combined-opt -i case20_elastic_wave.i
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 100
  ny   = 5
  xmin = 0.0
  xmax = 10.0
  ymin = 0.0
  ymax = 0.5
[]

# Provide displacements list once for all SolidMechanics objects.
[GlobalParams]
  displacements = 'disp_x disp_y'
[]

# Physics/SolidMechanics/Dynamic sets up:
#   - disp_x, disp_y variables (add_variables = true)
#   - vel_x, vel_y, accel_x, accel_y AuxVariables
#   - StressDivergenceTensors kernels for each displacement direction
#   - InertialForce kernels with Newmark-beta integration
#   - NewmarkAccelAux / NewmarkVelAux AuxKernels
# Specifying density here wires it to the inertial force kernels.
[Physics]
  [SolidMechanics]
    [Dynamic]
      [all]
        add_variables  = true
        newmark_beta   = 0.25
        newmark_gamma  = 0.5
        strain         = SMALL
        density        = 7800       # kg/m^3 — steel
        generate_output = 'stress_xx stress_yy vonmises_stress'
      []
    []
  []
[]

[BCs]
  # ---- Fixed end (left): prescribed displacement for the pulse ----
  # We drive a displacement pulse rather than a Neumann traction,
  # which avoids the coupling between the Pressure BC and the
  # Newmark predictor.  The amplitude is chosen so that the
  # equivalent stress sigma = E * du/dx ~ 1 MPa for a 1 mm pulse.
  #
  # Actually we use a Pressure (traction) BC nested inside [BCs],
  # following the MOOSE Dynamic action convention.
  # The [Pressure] sub-block applies a normal traction to the boundary.
  # A positive value pushes inward (compression in +x direction).
  [Pressure]
    [pulse_left]
      boundary = left
      function = pressure_pulse
      # factor scales the function value; keep at 1 since the
      # function already returns Pascals.
      factor = 1
    []
  []

  # ---- Constrain y-displacement on the bottom edge ----
  # Prevents rigid-body rotation and keeps the bar in plane.
  # y is also fixed on top to enforce plane-strain-like symmetry
  # (the strip is very thin so transverse modes are negligible).
  [fix_y_bottom]
    type     = DirichletBC
    variable = disp_y
    boundary = bottom
    value    = 0.0
  []
  [fix_y_top]
    type     = DirichletBC
    variable = disp_y
    boundary = top
    value    = 0.0
  []

  # ---- Left end: allow only axial motion ----
  # The pressure pulse drives disp_x; pin disp_y on left face
  # so the left end does not drift laterally.
  [fix_y_left]
    type     = DirichletBC
    variable = disp_y
    boundary = left
    value    = 0.0
  []

  # Right end is free (no BC on disp_x at right) — this is where
  # the wave reflects.  At a free end the stress goes to zero,
  # which doubles the displacement and reverses the wave sign.
[]

# Pressure pulse: 1 MPa ramp-up over 0.2 ms, hold for 0.2 ms,
# ramp-down over 0.2 ms.  Total pulse width = 0.4 ms.
# This is short compared to the 2 ms transit time so the pulse
# travels as a spatially compact packet.
[Functions]
  [pressure_pulse]
    type = PiecewiseLinear
    x = '0.0    0.0002  0.0004  0.0006'
    y = '0.0    1.0e6   1.0e6   0.0'
  []
[]

[Materials]
  # Isotropic elasticity tensor for steel.
  # Poisson ratio = 0.0 decouples x and y stress/strain, giving a
  # clean 1D wave without lateral constraint effects.
  [elasticity]
    type           = ComputeIsotropicElasticityTensor
    youngs_modulus = 200.0e9    # 200 GPa — steel
    poissons_ratio = 0.0        # pure 1D wave propagation
  []
  # Linear elastic stress from small-strain tensor.
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Postprocessors]
  # Axial displacement at the free right end (mid-height).
  # Shows the arriving wave front and the doubled displacement
  # at the moment of reflection.
  [disp_x_right]
    type  = PointValue
    variable = disp_x
    point = '10.0 0.25 0'
  []

  # Axial displacement at the driven left end (mid-height).
  # Shows the initial displacement from the pulse and then the
  # returning wave arriving back from the right.
  [disp_x_left]
    type  = PointValue
    variable = disp_x
    point = '0.0 0.25 0'
  []

  # Domain-average axial stress — useful for checking overall
  # energy conservation (should return to near zero after the
  # pulse has fully traversed and exited both ends).
  [avg_stress_xx]
    type     = ElementAverageValue
    variable = stress_xx
  []
[]

[Executioner]
  type = Transient

  # Implicit solve — suitable for the implicit Newmark-beta scheme.
  # LU factorization is robust for this small 2D problem.
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  # ---- Time control ----
  # Wave transit time one way: ~0.002 s
  # Run for 0.006 s to observe:
  #   t ~ 0.002 s : pulse arrives at right end (free reflection)
  #   t ~ 0.004 s : reflected wave arrives back at left end
  #   t ~ 0.006 s : second reflection from left (fixed via load release)
  dt       = 2.0e-5    # 20 microseconds per step — resolves the wave
  end_time = 0.006     # 300 time steps total

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-12
  nl_max_its = 20

  l_tol     = 1e-8
  l_max_its = 150
[]

[Outputs]
  csv = true
  [exodus]
    type     = Exodus
    # Write every 5th step — keeps file size manageable while still
    # providing smooth animation of the traveling wave front.
    time_step_interval = 5
  []
[]
