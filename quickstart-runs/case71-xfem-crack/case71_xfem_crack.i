# ============================================================
# Case 71: XFEM — Heat Conduction with Stationary Crack
# Transient heat conduction in a square plate containing a
# stationary edge crack. The crack runs vertically from the
# bottom edge to the center, acting as a perfect insulator.
#
# XFEM (eXtended Finite Element Method) enriches the standard
# FE space to represent the discontinuity without remeshing.
# Elements cut by the crack have enriched DOFs on each side.
#
# Physics:
#   dT/dt = nabla^2 T  in [0,1]^2
#   Crack: insulated line from (0.5, 0) to (0.5, 0.5)
#   T(x=0) = 1 (hot), T(x=1) = 0 (cold)
#   dT/dn = 0 on top and bottom (insulated)
#   T(t=0) = 0
#
# Expected behavior:
#   Heat enters from the left. In the upper half (above the
#   crack tip), heat flows freely left-to-right. In the lower
#   half, the crack blocks heat flow at x=0.5:
#   - Lower-left quadrant: temperature rises toward 1
#   - Lower-right quadrant: temperature stays near 0
#   - Sharp discontinuity across the crack face
#
# Domain: [0, 1]^2, 20x20 QUAD4 elements
# Time: t in [0, 0.5], dt = 0.02 (25 steps)
# ============================================================

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  elem_type = QUAD4
[]

[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]

[UserObjects]
  # Vertical crack from bottom edge to center
  [crack]
    type = LineSegmentCutUserObject
    cut_data = '0.5 0.0 0.5 0.5'
  []
[]

[Variables]
  [T]
  []
[]

[Kernels]
  [time]
    type = TimeDerivative
    variable = T
  []
  [diff]
    type = Diffusion
    variable = T
  []
[]

[BCs]
  [left_hot]
    type = DirichletBC
    variable = T
    boundary = left
    value = 1.0
  []
  [right_cold]
    type = DirichletBC
    variable = T
    boundary = right
    value = 0.0
  []
  # Top and bottom: natural BC (insulated) — no block needed
[]

[Postprocessors]
  [T_avg]
    type = ElementAverageValue
    variable = T
    execute_on = 'initial timestep_end'
  []
  [T_max]
    type = ElementExtremeValue
    variable = T
    value_type = max
    execute_on = 'initial timestep_end'
  []
  [T_min]
    type = ElementExtremeValue
    variable = T
    value_type = min
    execute_on = 'initial timestep_end'
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  dt = 0.02
  end_time = 0.5
  l_tol = 1e-6
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
[]
