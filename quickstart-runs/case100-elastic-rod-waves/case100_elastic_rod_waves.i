# ============================================================
# Case 100: Elastic Wave Propagation on Thin Rod
# MIT 6.641, Lec 16 — Elastodynamics / Waves
# Prof. Markus Zahn, Spring 2005
#
# Longitudinal elastic waves in a thin rod satisfy:
#
#   ∂²ξ/∂t² = (E/ρ) ∂²ξ/∂x²
#
# where ξ is the axial displacement, E is Young's modulus,
# ρ is the density, and the phase velocity is v_p = √(E/ρ).
#
# A sinusoidal displacement is applied at one end:
#   ξ(0, t) = A sin(ωt)
#
# The other end is free (zero stress), creating standing
# wave patterns at resonant frequencies:
#   f_n = n v_p / (2L),   n = 1, 2, 3, ...
#
# Applications: ultrasonic sensors, non-destructive testing,
# musical instruments, pile driving dynamics.
#
# This uses the solid_mechanics module with dynamic Newmark-β
# time integration, identical to Case 20 but with sinusoidal
# driving for resonance analysis.
#
# Domain: [0, 10] × [0, 0.5], 100×5 elements (thin strip)
# Material: Aluminium — E = 70 GPa, ν = 0, ρ = 2700 kg/m³
# Wave speed: v_p = √(70e9/2700) ≈ 5092 m/s
# Transit time: L/v_p = 10/5092 ≈ 0.00196 s
# Fundamental resonance: f₁ = v_p/(2L) ≈ 254.6 Hz
# Driving frequency: ω = 2π·250 ≈ 1571 rad/s (near resonance)
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

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

# Dynamic solid mechanics with Newmark-β integration.
# Sets up stress divergence + inertial force + velocity/accel aux variables.
[Physics]
  [SolidMechanics]
    [Dynamic]
      [all]
        add_variables  = true
        newmark_beta   = 0.25
        newmark_gamma  = 0.5
        strain         = SMALL
        density        = 2700     # kg/m³ — aluminium
        generate_output = 'stress_xx vonmises_stress'
      []
    []
  []
[]

[BCs]
  # Left end: sinusoidal displacement driving at ~250 Hz.
  # ξ(0,t) = A sin(ωt), A = 1e-4 m (0.1 mm amplitude).
  [drive_left]
    type     = FunctionDirichletBC
    variable = disp_x
    boundary = left
    function = sinusoidal_drive
  []

  # Constrain y-displacement everywhere to enforce 1D wave motion.
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
  [fix_y_left]
    type     = DirichletBC
    variable = disp_y
    boundary = left
    value    = 0.0
  []
  [fix_y_right]
    type     = DirichletBC
    variable = disp_y
    boundary = right
    value    = 0.0
  []

  # Right end: free (zero stress). Natural Neumann BC.
  # The wave reflects here with sign reversal of stress.
[]

[Functions]
  # Sinusoidal driving displacement at ~250 Hz.
  # ω = 2π × 250 = 1570.8 rad/s
  # Amplitude = 1e-4 m for small-strain regime.
  [sinusoidal_drive]
    type       = ParsedFunction
    expression = '1.0e-4 * sin(1570.8 * t)'
  []
[]

[Materials]
  [elasticity]
    type           = ComputeIsotropicElasticityTensor
    youngs_modulus  = 70.0e9    # 70 GPa — aluminium
    poissons_ratio  = 0.0      # pure 1D wave
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Postprocessors]
  # Displacement at the free right end: shows the standing wave buildup.
  [disp_x_right]
    type     = PointValue
    variable = disp_x
    point    = '10.0 0.25 0'
  []

  # Displacement at the driven left end.
  [disp_x_left]
    type     = PointValue
    variable = disp_x
    point    = '0.0 0.25 0'
  []

  # Displacement at rod midpoint: node for mode 1, antinode for mode 2.
  [disp_x_mid]
    type     = PointValue
    variable = disp_x
    point    = '5.0 0.25 0'
  []

  # Average axial stress for energy monitoring.
  [avg_stress_xx]
    type     = ElementAverageValue
    variable = stress_xx
  []
[]

[Executioner]
  type = Transient

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  # Run for 10 ms ≈ 5 transit times, enough for standing wave buildup.
  # dt = 2e-5 s gives ~50 steps per period at 250 Hz (T = 4ms).
  dt       = 2.0e-5
  end_time = 0.01

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-12
  nl_max_its = 20

  l_tol     = 1e-8
  l_max_its = 150
[]

[Outputs]
  csv = true
  [exodus]
    type = Exodus
    time_step_interval = 5
  []
[]
