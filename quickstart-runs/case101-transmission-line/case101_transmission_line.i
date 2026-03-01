# ============================================================
# Case 101: Transmission Line — RC Cable Transient Response
# MIT 6.641, Lec 17 — Transmission Lines
# Prof. Markus Zahn, Spring 2005
#
# For a lossy transmission line with series resistance R and
# shunt capacitance C per unit length (the RC limit where
# inductance L is negligible), the telegrapher's equation
# reduces to the diffusion equation:
#
#   ∂v/∂t = D · ∂²v/∂x²     where D = 1/(R·C)
#
# This describes signal propagation on submarine telegraph
# cables, VLSI interconnects, and biological axons (cable eq.).
#
# A step voltage V₀ is applied at x = 0 and the other end
# (x = L) is open-circuit (zero flux). The analytic solution
# for the semi-infinite line is:
#
#   v(x,t) = V₀ · erfc(x / (2√(Dt)))
#
# The voltage "diffuses" along the cable with characteristic
# delay time t_delay = x²·R·C / 4 for each position x.
#
# Domain: [0, 1] × [0, 0.04] quasi-1D, 100×2 mesh
# IC: v = 0 (uncharged cable)
# BC: v(0,t) = V₀ = 1 (step voltage), ∂v/∂x(L,t) = 0 (open end)
# Parameters: D = 0.1 (normalised RC diffusivity)
# ============================================================

D_rc = 0.1   # diffusivity = 1/(R·C) [m²/s]
V0   = 1.0   # step voltage [V]

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 100     # fine resolution along cable
  ny   = 2       # quasi-1D
  xmin = 0
  xmax = 1.0     # cable length L = 1 m
  ymin = 0
  ymax = 0.04
[]

[Variables]
  # Voltage v along the cable [V].
  # Initially zero (uncharged), driven by step source at x = 0.
  [v]
    initial_condition = 0
  []
[]

[Kernels]
  # ∂v/∂t: rate of voltage change at each point.
  [v_time]
    type     = ADTimeDerivative
    variable = v
  []

  # D · ∂²v/∂x²: voltage diffusion along the cable.
  # ADMatDiffusion with diffusivity = D = 1/(RC).
  [v_diffusion]
    type        = ADMatDiffusion
    variable    = v
    diffusivity = diffusivity
  []
[]

[BCs]
  # Source end (x = 0): step voltage V₀ = 1.0 applied at t = 0.
  # A small ramp (0.01 s) avoids numerical shock.
  [source_step]
    type     = FunctionDirichletBC
    variable = v
    boundary = left
    function = step_voltage
  []

  # Load end (x = 1): open circuit → zero flux (natural Neumann).
  # No explicit BC needed; MOOSE applies dv/dx = 0 by default.
[]

[Functions]
  # Step voltage with smooth ramp-up over 0.01 s.
  [step_voltage]
    type = PiecewiseLinear
    x    = '0      0.01    100'
    y    = '0      ${V0}   ${V0}'
  []
[]

[Materials]
  [rc_cable]
    type        = ADGenericConstantMaterial
    prop_names  = 'diffusivity'
    prop_values = '${D_rc}'
  []
[]

[Postprocessors]
  # Voltage at x = 0.25 (quarter length): tracks the diffusing front.
  [v_quarter]
    type     = PointValue
    variable = v
    point    = '0.25 0.02 0'
  []

  # Voltage at x = 0.50 (midpoint).
  [v_mid]
    type     = PointValue
    variable = v
    point    = '0.50 0.02 0'
  []

  # Voltage at x = 0.75 (three-quarter length).
  [v_three_quarter]
    type     = PointValue
    variable = v
    point    = '0.75 0.02 0'
  []

  # Voltage at x = 1.0 (load end): shows the latest arrival.
  [v_load]
    type     = PointValue
    variable = v
    point    = '1.0 0.02 0'
  []

  # Average voltage along cable: approaches V₀ at steady state.
  [avg_v]
    type     = ElementAverageValue
    variable = v
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  # Diffusion timescale to x = 1: t ~ L²/(4D) = 1/(4·0.1) = 2.5 s.
  # Run to 5 s to see full charging.
  dt       = 0.05
  end_time = 5.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv    = true
[]
