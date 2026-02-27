# Case 69: MultiApp Coupled Diffusion — Bidirectional Picard Iteration

## Overview

Many real-world physics problems involve distinct subsystems that are tightly coupled but most naturally described by separate sets of equations. MOOSE's MultiApp system provides a principled way to compose these subsystems: each physics lives in its own "app" with its own mesh, variables, and solvers, and data is exchanged between apps via Transfer objects. The coupling loop — solve main, transfer to sub, solve sub, transfer back, repeat — is the Picard (fixed-point) iteration algorithm, which converges when all transferred fields are self-consistent across the boundary between apps.

This case models a bidirectional coupling between a heat conduction problem and a reaction-diffusion source problem. The main app solves dT/dt = nabla^2 T + q(x,y), where q is a heat source produced by the sub-app. The sub-app solves dq/dt = nabla^2 q - q + S(x,y) - 0.5*T, where S is a fixed Gaussian forcing and the -0.5*T term provides negative feedback from the main app's temperature. As T rises, q is suppressed; as q falls, T stops rising. This feedback loop ensures the coupled system converges to a stable steady state rather than growing without bound.

Key concepts demonstrated:

- `TransientMultiApp` for running a sub-application synchronized with the main transient time stepper
- `MultiAppNearestNodeTransfer` for interpolating field data between apps at shared time steps
- Picard (fixed-point) iteration controlled by `fixed_point_max_its`, `fixed_point_rel_tol`, and `fixed_point_abs_tol`
- `CoupledForce` for adding a field-dependent source term, and `Reaction` for a linear decay term
- Two-way transfer: main sends T to sub, sub sends q back to main

---

## The Physics

### Main App — Heat Conduction with External Source

The main app governs the temperature field T(x,y,t) on the unit square [0,1]^2:

```
dT/dt = nabla^2 T + q(x,y,t)
```

Boundary conditions: T = 0 on all four sides (Dirichlet). Initial condition: T = 0.

The source q is not solved by the main app — it is received via transfer from the sub-app. The main app stores q in an AuxVariable `source_from_sub` and applies it as a body force via the `CoupledForce` kernel. At the start of each Picard iteration, q holds the value from the previous iteration (or the previous time step for the first iteration).

### Sub-App — Reaction-Diffusion Source with Temperature Feedback

The sub-app governs the source field q(x,y,t) on the same unit square:

```
dq/dt = nabla^2 q - q + S(x,y) - 0.5 * T_from_main
```

where:
- `S(x,y) = 10 * exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.05)` — Gaussian forcing, amplitude 10, centered at (0.5, 0.5), half-width 0.1
- `-q` — first-order linear decay (enforced by the `Reaction` kernel)
- `-0.5 * T_from_main` — negative feedback: higher temperature suppresses the source

Boundary conditions: q = 0 on all four sides. Initial condition: q = 0.

The temperature T_from_main is an AuxVariable in the sub-app that receives the main app's T field via transfer. The feedback coefficient 0.5 was chosen so that the steady-state source is significantly reduced from its zero-feedback value, making the coupling effect visible.

### Steady-State Analysis

At steady state (all time derivatives vanish), the coupled system satisfies:

```
nabla^2 T + q = 0
nabla^2 q - q + S - 0.5*T = 0
```

The feedback loop means: larger S drives larger q, which drives larger T, which suppresses q. The balance is:

```
q_ss ~ S / (1 + 0.5 * chi)
```

where chi is the effective Green's function relating T to q (a positive quantity). The result is that the steady-state source is smaller than S alone would predict. The CSV data confirms: maximum T reaches only ~0.0123 despite the Gaussian amplitude of 10.

### Picard Convergence

At each time step, MOOSE performs Picard iterations:
1. Solve main app (with q from previous Picard iteration)
2. Transfer T from main to sub
3. Solve sub-app (with T from step 2)
4. Transfer q from sub to main
5. Check convergence: if |q_new - q_old| / |q_old| < `fixed_point_rel_tol`, stop; else go to step 1

Convergence is fast here — typically 2 iterations per step — because the coupling is moderate (coefficient 0.5) and the time step is not too large.

### Domain

- Square [0,1]^2, 20x20 QUAD elements (both apps share the same mesh)
- Time: t in [0, 2.0], dt = 0.2 (10 time steps)
- The steady state is effectively reached by t = 1.0

---

## Input File Walkthrough

The input consists of two files: `case69_multiapp_coupling.i` (main app) and `case69_sub.i` (sub-app).

### Main App: `[Mesh]`

```
type = GeneratedMesh
dim = 2
nx = 20
ny = 20
```

A 20x20 structured quadrilateral mesh on [0,1]^2. Both the main app and sub-app use an identical mesh, which makes the `MultiAppNearestNodeTransfer` trivially exact — each node maps to itself.

### Main App: `[Variables]` and `[AuxVariables]`

```
[Variables]
  [T]
    initial_condition = 0
  []
[]
[AuxVariables]
  [source_from_sub]
  []
[]
```

`T` is the solved PDE variable. `source_from_sub` is an AuxVariable that holds the q field received from the sub-app. AuxVariables are not solved — they are updated by AuxKernels or Transfers.

### Main App: `[Kernels]`

```
[time]
  type = TimeDerivative
  variable = T
[]
[diffusion]
  type = Diffusion
  variable = T
[]
[source]
  type = CoupledForce
  variable = T
  v = source_from_sub
[]
```

`CoupledForce` adds -v * test to the residual, which places +v on the right-hand side of the PDE. This implements the `+ q` source term in the heat equation.

### Main App: `[MultiApps]`

```
[sub]
  type = TransientMultiApp
  positions = '0 0 0'
  input_files = case69_sub.i
  execute_on = 'timestep_end'
[]
```

`TransientMultiApp` runs the sub-app synchronously with the main app's time stepper. `positions = '0 0 0'` specifies the origin of the sub-app's coordinate system relative to the main app (no offset here since both domains coincide). `execute_on = timestep_end` means the sub-app is advanced after the main app completes a time step.

### Main App: `[Transfers]`

```
[source_from_sub]
  type = MultiAppNearestNodeTransfer
  from_multi_app = sub
  source_variable = q
  variable = source_from_sub
[]
[T_to_sub]
  type = MultiAppNearestNodeTransfer
  to_multi_app = sub
  source_variable = T
  variable = T_from_main
[]
```

Two transfers operate at each Picard iteration. `from_multi_app` pulls data from the sub-app into the main app; `to_multi_app` pushes data from the main app into the sub-app. `MultiAppNearestNodeTransfer` maps values by finding the nearest node in the target mesh — exact for coincident meshes.

### Main App: `[Executioner]`

```
type = Transient
dt = 0.2
end_time = 2.0
fixed_point_max_its = 10
fixed_point_rel_tol = 1e-6
fixed_point_abs_tol = 1e-8
```

The Picard loop parameters control how many iterations are allowed (`fixed_point_max_its = 10`) and what convergence criterion is used. If the relative change in the transferred variables between Picard iterations drops below 1e-6, the loop exits early. In practice, this problem converges in 2-3 iterations.

### Sub-App: `[Kernels]`

```
[decay]
  type = Reaction
  variable = q
[]
[source]
  type = BodyForce
  variable = q
  function = source_func
[]
[feedback]
  type = CoupledForce
  variable = q
  v = T_from_main
  coef = -0.5
[]
```

`Reaction` adds +lambda*u*test to the residual (with lambda=1 by default), implementing the -q decay term. `BodyForce` adds the Gaussian function S(x,y). `CoupledForce` with `coef = -0.5` adds +0.5*T to the residual, which implements -0.5*T on the right-hand side of the PDE (the negative feedback).

---

## Running the Simulation

The main app automatically spawns the sub-app; only the main input file is needed on the command line.

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case69-multiapp-coupling \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case69_multiapp_coupling.i 2>&1 | tail -30'
```

Output files:
- `case69_multiapp_coupling_out.e` — Main app Exodus with T and source_from_sub fields
- `case69_multiapp_coupling_out.csv` — Time series of T_max, T_avg, source_avg
- `case69_multiapp_coupling_out_sub0.e` — Sub-app Exodus with q and T_from_main fields

---

## Expected Results

### Time Evolution

| Time | T_max | T_avg | source_avg |
|------|-------|-------|------------|
| 0.0  | 0.000 | 0.000 | 0.000 |
| 0.2  | 0.00803 | 0.00307 | 0.0728 |
| 0.4  | 0.01109 | 0.00430 | 0.0875 |
| 0.6  | 0.01199 | 0.00466 | 0.0904 |
| 0.8  | 0.01223 | 0.00476 | 0.0910 |
| 1.0  | 0.01229 | 0.00478 | 0.0911 |
| 2.0  | 0.01231 | 0.00479 | 0.0912 |

The system reaches a near-steady state by t = 1.0. The temperature maximum is at the center of the domain (x = y = 0.5) where the Gaussian source is strongest.

### Feedback Attenuation

Without the -0.5*T feedback term, the sub-app source would reach a higher steady state controlled only by the S - q balance. The feedback suppresses q from its uncoupled value and consequently limits T. The source_avg at steady state (~0.091) is notably lower than the peak Gaussian value S_max = 10 * exp(0) = 10, because the diffusion and Dirichlet BCs spread and damp both fields significantly.

### Spatial Structure

The T field peaks at the center and decays to zero at the boundaries (Dirichlet T = 0). The q field has a similar Gaussian spatial structure, modified by the temperature feedback. Because the feedback is proportional to T, and T peaks at the center, the source suppression is also strongest at the center.

---

## Key Takeaways

- The `TransientMultiApp` + `MultiAppNearestNodeTransfer` combination is the standard pattern for operator-splitting multi-physics in MOOSE. Each physics can use different element types, mesh resolutions, or even different MOOSE apps built from different module combinations — the Transfer system handles interpolation between non-matching meshes.
- Picard (fixed-point) iteration is controlled through the Executioner block via `fixed_point_max_its` and `fixed_point_rel_tol`. The iteration automatically stops when the fields transferred between apps change by less than the specified tolerance between successive iterations.
- `CoupledForce` is a general-purpose kernel for adding field-dependent source terms. With `coef = -0.5`, it contributes -0.5*v to the right-hand side. Negative coefficients implement suppression/feedback mechanisms; positive coefficients implement forcing.
- The `Reaction` kernel adds a first-order linear reaction term (+lambda*u) to the residual. Combined with a negative forcing term, it provides the decay needed to prevent unbounded growth of a variable driven by external sources.
- `execute_on = timestep_end` in the MultiApps block means the sub-app advances after the main app's Newton solve completes. The order of execution (main first, then sub) combined with the Picard loop ensures both apps are in equilibrium with each other at the end of each time step.
- AuxVariables serve as the communication buffers between apps. The main app stores q received from the sub in `source_from_sub`; the sub stores T received from the main in `T_from_main`. These are read-only from the perspective of the PDE solvers but writable by Transfers.
- Bidirectional coupling with negative feedback is a stable pattern for coupled multi-physics. One-way coupling (no feedback) or positive feedback can lead to either uncoupled behavior or numerical instability. The feedback coefficient 0.5 here ensures the coupled eigenvalue is real and negative.
