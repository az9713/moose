# ============================================================
# Case 28: Two-Way Joule Heating — Temperature-Dependent σ(T)
# Melcher, Continuum Electromechanics (MIT, 1981), Ch. 10 §10.1
#
# Extends Case 17 with metallic conductivity that decreases
# as temperature rises (negative feedback):
#   -div(σ(T)·grad(V)) = 0                (current conservation)
#   ρcp·∂T/∂t = div(k·∇T) + σ(T)·|∇V|²   (Joule heating)
#   σ(T) = σ₀ / (1 + α·(T - T_ref))       (metallic: σ drops with T)
#
# Domain: [0,2] × [0,1], 40×20 mesh (same as Case 17)
# BC: V = 10 V at left, V = 0 at right
#     T = 300 K at left and right (electrode heat sinks)
# Parameters: σ₀ = 1e6 S/m, α = 0.004 /K, T_ref = 300 K
#             k = 50 W/(m·K), ρ = 8000 kg/m³, cp = 500 J/(kg·K)
# ============================================================

# σ(T) = σ₀ / (1 + α·(T - T_ref))
# σ₀ = 1e6, α = 0.004, T_ref = 300 → tabulated via ADPiecewiseLinearInterpolationMaterial

[Mesh]
  [gen]
    type  = GeneratedMeshGenerator
    dim   = 2
    nx    = 40
    ny    = 20
    xmin  = 0
    xmax  = 2
    ymin  = 0
    ymax  = 1
  []
[]

[Variables]
  # Electric potential [V]
  [V]
    initial_condition = 0.0
  []
  # Temperature [K]
  [T]
    initial_condition = 300.0
  []
[]

[Kernels]
  # --- Electric potential: -div(σ(T)·grad(V)) = 0 ---
  # ADHeatConduction reads 'thermal_conductivity' by default.
  # We point it to 'electrical_conductivity' (computed from T).
  [V_diff]
    type                 = ADHeatConduction
    variable             = V
    thermal_conductivity = electrical_conductivity
  []

  # --- Heat equation: ρcp·∂T/∂t = div(k·grad(T)) + Q_joule ---
  [T_time]
    type     = ADHeatConductionTimeDerivative
    variable = T
  []
  [T_diff]
    type     = ADHeatConduction
    variable = T
  []
  # Joule heating source: Q = σ(T)·|grad(V)|²
  [T_joule]
    type         = ADJouleHeatingSource
    variable     = T
    heating_term = electric_field_heating
  []
[]

[BCs]
  # Voltage: 10 V at left electrode, 0 V at right
  [V_left]
    type     = ADDirichletBC
    variable = V
    boundary = left
    value    = 10.0
  []
  [V_right]
    type     = ADDirichletBC
    variable = V
    boundary = right
    value    = 0.0
  []
  # Temperature: 300 K at both electrodes (heat sinks)
  [T_left]
    type     = ADDirichletBC
    variable = T
    boundary = left
    value    = 300.0
  []
  [T_right]
    type     = ADDirichletBC
    variable = T
    boundary = right
    value    = 300.0
  []
[]

[Materials]
  # Thermal properties (used by ADHeatConduction and ADHeatConductionTimeDerivative)
  [thermal]
    type        = ADGenericConstantMaterial
    prop_names  = 'thermal_conductivity specific_heat density'
    prop_values = '50.0              500.0        8000.0'
  []

  # Temperature-dependent electrical conductivity:
  # σ(T) = σ₀ / (1 + α·(T - T_ref))
  # This is the standard metallic model: resistivity increases linearly
  # with temperature, so conductivity decreases hyperbolically.
  # ADPiecewiseLinearInterpolationMaterial provides AD-compatible tabulated
  # data without needing JIT compilation (which is unavailable in Docker).
  # Tabulated from σ(T) = 1e6 / (1 + 0.004*(T-300)):
  #   T=300 → σ=1e6,  T=400 → σ=714286,  T=500 → σ=555556,  T=600 → σ=454545
  [sigma_of_T]
    type         = ADPiecewiseLinearInterpolationMaterial
    property     = electrical_conductivity
    variable     = T
    x            = '250  300   350     400     450     500     550     600     700'
    y            = '1.25e6  1e6 833333  714286  625000  555556  500000  454545  384615'
    extrapolation = true
  []

  # ElectromagneticHeatingMaterial computes Q = σ·|grad(V)|²
  # It reads 'electrical_conductivity' (now T-dependent) and the
  # gradient of V to produce the 'electric_field_heating' property.
  [joule_material]
    type                        = ElectromagneticHeatingMaterial
    electric_field              = V
    electric_field_heating_name = electric_field_heating
    electrical_conductivity     = electrical_conductivity
    formulation                 = time
    solver                      = electrostatic
  []
[]

[Postprocessors]
  [max_T]
    type       = ElementExtremeValue
    variable   = T
    value_type = max
  []
  [avg_T]
    type     = ElementAverageValue
    variable = T
  []
  [max_V]
    type       = ElementExtremeValue
    variable   = V
    value_type = max
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
  dt       = 0.25
  end_time = 5.0
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 20
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Outputs]
  exodus = true
  csv    = true
[]
