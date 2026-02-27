# ============================================================
# Case 72: THM Pipe Flow — 1D Single-Phase Compressible
# Uses the thermal_hydraulics module's component-based DSL
# to model 1D compressible gas flow through a pipe with
# friction.
#
# An ideal gas (gamma = 1.4) flows through a pipe driven by
# an inlet mass flow rate and constrained by an outlet
# pressure. The high friction factor (f = 10) creates a
# significant pressure drop along the pipe.
#
# Components:
#   FlowChannel1Phase: L = 1 m, A = 1 m^2, 50 elements
#   Inlet: m_dot = 2 kg/s, T = 500 K
#   Outlet: p = 2e5 Pa
#
# Initial conditions: T = 300 K, p = 1e5 Pa, vel = 1 m/s
# These are far from steady state — the transient shows the
# system adjusting to the boundary conditions.
#
# Expected: Pressure rises to accommodate the outlet BC;
# temperature propagates from inlet (500 K) toward outlet;
# velocity adjusts to satisfy mass flow constraint.
#
# Time: t in [0, 0.5], dt = 0.025 (20 steps)
# ============================================================

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
    gamma = 1.4
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Components]
  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'pipe:in'
    m_dot = 2
    T = 500
  []

  [pipe]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    gravity_vector = '0 0 0'
    length = 1.0
    n_elems = 50
    A = 1.0

    initial_T = 300
    initial_p = 1e5
    initial_vel = 1

    f = 10.0
    closures = simple_closures
    fp = fp

    # Critical for conditioning: scale energy equation
    scaling_factor_1phase = '1 1 1e-5'
  []

  [outlet]
    type = Outlet1Phase
    input = 'pipe:out'
    p = 2e5
  []
[]

[Postprocessors]
  [T_near_inlet]
    type = PointValue
    variable = T
    point = '0.01 0 0'
    execute_on = 'initial timestep_end'
  []
  [T_mid]
    type = PointValue
    variable = T
    point = '0.5 0 0'
    execute_on = 'initial timestep_end'
  []
  [T_near_outlet]
    type = PointValue
    variable = T
    point = '0.99 0 0'
    execute_on = 'initial timestep_end'
  []
  [p_near_inlet]
    type = PointValue
    variable = p
    point = '0.01 0 0'
    execute_on = 'initial timestep_end'
  []
  [p_near_outlet]
    type = PointValue
    variable = p
    point = '0.99 0 0'
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 0.5
  dt = 0.025
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  nl_max_its = 15
  l_tol = 1e-3
  l_max_its = 10
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  [Quadrature]
    type = GAUSS
    order = SECOND
  []
[]

[Outputs]
  exodus = true
  csv = true
[]
