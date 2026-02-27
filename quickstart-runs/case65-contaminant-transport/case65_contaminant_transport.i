# ============================================================
# Case 65: Contaminant Transport — Advection-Dispersion-Reaction
# A dissolved contaminant is transported through a 1D porous
# medium by Darcy flow while undergoing first-order decay.
# Uses the chemical_reactions module's transport kernels.
#
# Governing equation:
#   phi * dc/dt + div(c*q) - D*∇²c + phi*k*c = 0
# where q = -K*grad(P) is the Darcy flux, D = diffusivity,
# k = first-order decay rate, phi = porosity.
#
# Setup:
#   - Prescribed linear pressure field P(x) = 2 - x (left=2, right=0)
#   - Darcy flux q = -K * dP/dx = K (left to right)
#   - Contaminant pulse at left end (BoundingBoxIC)
#   - ChemicalOutFlowBC at right (advective outflow)
#   - First-order decay k = 0.01 /s
#
# Parameters:
#   phi = 0.3, D = 1e-3 m²/s, K = 1e-3 m/s
#   Domain: [0, 1] x [0, 0.1], 50x2 elements
#   Run for 100 s to see advection + dispersion + decay
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim  = 2
  nx   = 50
  ny   = 2
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 0.1
[]

[Variables]
  [c]
    order = FIRST
    family = LAGRANGE
    [InitialCondition]
      type = BoundingBoxIC
      x1 = 0.0
      y1 = 0.0
      x2 = 0.05
      y2 = 0.1
      inside = 1.0
      outside = 0.0
      variable = c
    []
  []
[]

[AuxVariables]
  [pressure]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [pressure_ic]
    type = FunctionIC
    variable = pressure
    function = '2-x'    # linear gradient, left=2, right=0
  []
[]

[Kernels]
  # Porosity-weighted time derivative
  [c_time]
    type = PrimaryTimeDerivative
    variable = c
  []
  # Molecular diffusion
  [c_diff]
    type = PrimaryDiffusion
    variable = c
  []
  # Darcy convection (uses pressure gradient)
  [c_conv]
    type = PrimaryConvection
    variable = c
    p = pressure
  []
  # First-order decay: +k*c*test in residual
  [c_decay]
    type = CoefReaction
    variable = c
    coefficient = 0.01    # k = 0.01 /s decay rate
  []
[]

[BCs]
  # Left: fixed concentration source
  [c_left]
    type = DirichletBC
    variable = c
    boundary = left
    value = 1.0
  []
  # Right: advective outflow
  [c_right]
    type = ChemicalOutFlowBC
    variable = c
    boundary = right
  []
[]

[Materials]
  [porous]
    type = GenericConstantMaterial
    prop_names  = 'diffusivity conductivity porosity'
    prop_values = '1e-3        1e-3         0.3'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  # Concentration at x=0.25
  [c_025]
    type = PointValue
    variable = c
    point = '0.25 0.05 0'
  []
  # Concentration at x=0.5
  [c_050]
    type = PointValue
    variable = c
    point = '0.5 0.05 0'
  []
  # Concentration at x=0.75
  [c_075]
    type = PointValue
    variable = c
    point = '0.75 0.05 0'
  []
  # Average concentration (should decrease from decay)
  [c_avg]
    type = ElementAverageValue
    variable = c
  []
  # Total mass (integral)
  [c_total]
    type = ElementIntegralVariablePostprocessor
    variable = c
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  dt       = 2
  end_time = 100

  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-8
[]

[Outputs]
  exodus = true
  csv    = true
[]
