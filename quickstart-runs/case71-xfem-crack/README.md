# Case 71: XFEM — Heat Conduction with Stationary Crack

## Overview

The extended finite element method (XFEM) is a powerful technique for representing discontinuities — cracks, material interfaces, phase boundaries — within a standard finite element mesh without requiring the mesh to conform to the geometry of the discontinuity. In classical FEM, a crack must either lie along element boundaries (requiring mesh generation aligned to the crack) or be represented by very fine mesh refinement. As a crack propagates, the mesh must be updated — an expensive and difficult operation. XFEM avoids this entirely by enriching the standard polynomial basis functions with discontinuous enrichment functions that capture the crack behavior, allowing the crack to pass through elements at arbitrary angles without any remeshing.

This case demonstrates XFEM for transient heat conduction in a unit square plate containing a stationary vertical crack. The crack runs from the bottom edge at x = 0.5 up to the center of the plate at y = 0.5, acting as a perfect thermal insulator (zero heat flux across the crack face). Heat flows from left (T = 1) to right (T = 0), but the crack blocks flow in the lower half of the domain. The result is a striking temperature field: the lower-left quadrant heats up while the lower-right quadrant stays cold, with a sharp discontinuity across the crack. Above the crack tip at y = 0.5, heat flows freely and the temperature varies smoothly.

Key concepts demonstrated:

- `[XFEM]` block with `qrule = volfrac` for cut-element quadrature
- `LineSegmentCutUserObject` to define the crack geometry as a line segment
- `output_cut_plane = true` to write the crack plane to the Exodus output for visualization
- No Constraint block needed for an insulating crack — the XFEM cut naturally creates a zero-flux interface
- Sharp temperature discontinuities across the crack captured without any mesh modification

---

## The Physics

### Governing Equation

Transient heat conduction with uniform thermal diffusivity (set to 1 in non-dimensional units):

```
dT/dt = nabla^2 T     in Omega = [0,1]^2
```

Subject to:
- T = 1 on the left boundary (x = 0)
- T = 0 on the right boundary (x = 1)
- dT/dn = 0 on the top and bottom boundaries (natural BC — insulated)
- T = 0 at t = 0 everywhere (cold initial condition)
- Zero heat flux across both crack faces (insulating crack)

### XFEM Enrichment

In regions away from the crack, the standard bilinear (Q1) basis functions are used. Elements that are cut by the crack geometry receive enriched degrees of freedom. The enrichment uses the generalized Heaviside function:

```
phi_H(x) = H(f(x)) * phi_std(x)
```

where H is the Heaviside step function (1 on one side of the crack, 0 on the other) and f(x) is the level set function that defines the crack geometry. The enrichment allows elements to carry two distinct sets of DOFs — one for each side of the crack — with no coupling between them across the crack face.

For an insulating crack, the natural boundary condition (dT/dn = 0) is automatically satisfied because the enrichment creates a zero-flux interface: the two sides of a cut element evolve independently with no heat exchange. This is the key advantage of XFEM for insulating cracks — no additional constraint equations or penalty terms are needed.

### Crack Geometry

The crack is a vertical line segment from (0.5, 0.0) to (0.5, 0.5). This is an edge crack (it starts at the bottom boundary) that penetrates halfway through the domain. Elements in the lower half of the mesh that span x = 0.5 are cut by the crack; elements in the upper half are unaffected.

### Expected Temperature Field

By t = 0.5 (near steady state), the expected temperature distribution is:

- **Upper half (y > 0.5)**: Nearly 1D conduction from left to right, T(x) ~ 1 - x, smooth and undisturbed
- **Lower-left quadrant (x < 0.5, y < 0.5)**: Hot region, T ~ 1, heat enters from the left but cannot exit right across the crack
- **Lower-right quadrant (x > 0.5, y < 0.5)**: Cold region, T ~ 0, no heat source because the crack blocks the supply from the left
- **Crack tip region**: Heat can flow around the tip (from lower-left to upper half to lower-right), creating a gradual transition above the tip

The average temperature at steady state converges to approximately 0.496 — slightly below 0.5 because the lower-right quadrant is trapped at near-zero temperature.

---

## Input File Walkthrough

The input file is `case71_xfem_crack.i`.

### `[Mesh]`

```
type = GeneratedMesh
dim = 2
nx = 20
ny = 20
elem_type = QUAD4
```

A 20x20 structured QUAD4 mesh on [0,1]^2. XFEM requires QUAD4 (bilinear quadrilateral) or TRI3 elements — higher-order elements are not generally supported by the XFEM implementation in MOOSE. The mesh does not need to conform to the crack path; elements at x = 0.5 will be cut at run time.

### `[XFEM]`

```
[XFEM]
  qrule = volfrac
  output_cut_plane = true
[]
```

This block activates the XFEM system. `qrule = volfrac` selects the volume fraction quadrature rule for cut elements: the standard Gaussian quadrature points are used on each sub-domain of the cut element, weighted by the volume fraction of each sub-domain. This is the simplest and most robust quadrature approach for scalar problems. `output_cut_plane = true` causes the crack geometry to be written to the Exodus output file as an additional block, which allows visualization of the crack position in ParaView.

### `[UserObjects]`

```
[crack]
  type = LineSegmentCutUserObject
  cut_data = '0.5 0.0 0.5 0.5'
[]
```

`LineSegmentCutUserObject` defines the crack as a straight line segment specified by its two endpoints: (x1=0.5, y1=0.0, x2=0.5, y2=0.5). The XFEM system uses this object to identify which elements are cut, compute the enrichment functions, and determine the integration sub-domains in cut elements. For stationary cracks, the UserObject is computed once at the start and never updated.

The `cut_data` format is `'x1 y1 x2 y2'` for 2D cracks. In 3D, `LineSegmentCutUserObject` becomes a plane defined by a line segment swept along the normal direction.

### `[Variables]`

```
[T]
[]
```

A single temperature variable with no explicit initial condition (defaults to zero). XFEM automatically adds enriched DOFs to elements that are cut by the crack; these extra DOFs are also initialized to zero.

### `[Kernels]`

```
[time]
  type = TimeDerivative
  variable = T
[]
[diff]
  type = Diffusion
  variable = T
[]
```

Standard time derivative and Laplacian. XFEM modifies these kernels internally to operate correctly on cut elements — the volume integration is split across the two sub-domains, and no flux is computed across the crack face. No special XFEM-specific kernels are required for the standard heat equation.

### `[BCs]`

```
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
```

Only two Dirichlet BCs are needed. The top and bottom boundaries are naturally insulated (no Neumann block needed — zero flux is the default). The crack faces themselves have no BCs; the XFEM enrichment automatically enforces the insulating condition by decoupling the two sides.

### `[Postprocessors]`

```
[T_avg]  type = ElementAverageValue
[T_max]  type = ElementExtremeValue  value_type = max
[T_min]  type = ElementExtremeValue  value_type = min
```

These postprocessors operate on the standard variable `T`. In cut elements, XFEM provides separate values on each side of the crack; the postprocessors see the average value within each element sub-domain as a single value. This means T_avg reflects the true volume-averaged temperature, correctly accounting for the temperature jump across the crack.

### `[Executioner]`

```
type = Transient
dt = 0.02
end_time = 0.5
nl_rel_tol = 1e-10
nl_abs_tol = 1e-10
```

Twenty-five time steps of dt = 0.02 s. The tight nonlinear tolerances are appropriate for the enriched system where small errors in the cut-element quadrature can propagate. The BDF1 (backward Euler) time integrator is the default and is unconditionally stable for the heat equation.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case71-xfem-crack \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case71_xfem_crack.i 2>&1 | tail -30'
```

Output files:
- `case71_xfem_crack_out.e` — Exodus with temperature field and crack cut plane geometry
- `case71_xfem_crack_out.csv` — Time series of T_avg, T_max, T_min

---

## Expected Results

### Temperature Evolution

| Time (s) | T_avg | T_max | T_min |
|----------|-------|-------|-------|
| 0.00     | 0.000 | 0.000 | 0.000 |
| 0.02     | 0.142 | 0.937 | 0.000 |
| 0.06     | 0.262 | 0.973 | 0.000 |
| 0.10     | 0.335 | 0.982 | 0.001 |
| 0.20     | 0.433 | 0.991 | 0.002 |
| 0.30     | 0.473 | 0.994 | 0.003 |
| 0.40     | 0.489 | 0.995 | 0.004 |
| 0.50     | 0.496 | 0.995 | 0.004 |

T_max rapidly reaches ~0.99 (heat penetrates the left side quickly), while T_min stays near zero (the lower-right quadrant is trapped behind the crack). T_avg converges toward ~0.5 but approaches from below because the crack effectively makes the lower-right quadrant inaccessible to heat.

### Temperature Discontinuity

The most striking feature visible in the Exodus output is the sharp temperature jump across the crack. Elements cut by the crack carry two temperature values — one on each side — with no constraint relating them. Visualizing T in ParaView will show a distinct color change at x = 0.5 for y < 0.5, with no such discontinuity for y > 0.5 (above the crack tip). The `output_cut_plane = true` setting adds the crack geometry as a visible surface in the Exodus file.

### Comparison to Uncracked Case

Without the crack (standard Case 03 heat conduction setup), the temperature field would be symmetric about y = 0.5 and converge to T(x) = 1 - x at steady state, giving T_avg = 0.5 exactly. The crack reduces the accessible area for heat transport, pushing T_avg below 0.5 and creating the asymmetric pattern.

---

## Key Takeaways

- The `[XFEM]` block activates the extended finite element method for the entire simulation. The crack geometry is defined separately by a `UserObject` (`LineSegmentCutUserObject` for 2D straight cracks). No modifications to the mesh, variables, or kernels are required — XFEM operates transparently on top of the standard FEM framework.
- `qrule = volfrac` is the recommended quadrature rule for scalar XFEM problems. It uses volume fractions to weight the standard Gauss points, which is computationally efficient and sufficiently accurate for smooth solution components within each cut sub-domain. More sophisticated quadrature rules (like `moment_fitting`) provide higher accuracy but greater complexity.
- For an insulating crack, no Constraint block or penalty method is needed. The XFEM enrichment automatically decouples the two sides of cut elements, making the crack face a natural (zero-flux) boundary. This is a significant simplification compared to explicitly meshing the crack — no crack-face mesh, no boundary conditions, no interface elements.
- `output_cut_plane = true` writes the crack geometry to the Exodus file as an auxiliary element block. In ParaView, this appears as a separate surface that can be displayed or hidden independently of the solution fields. For crack growth problems, this allows visualization of the evolving crack front at each time step.
- XFEM in MOOSE supports both stationary cracks (this case) and propagating cracks (using `GeometricCutUserObject` subclasses that update the crack geometry based on stress intensity factors or other fracture criteria). Stationary crack cases like this one are useful for verification and for problems where the crack geometry is prescribed rather than predicted.
- The temperature field near the crack tip has a weak singularity — the flux is finite at the tip of an insulating crack in heat conduction, unlike the stress singularity (r^(-0.5)) at a traction-free crack tip in elasticity. This means standard polynomial enrichment is adequate for the thermal problem without adding crack-tip enrichment functions.
- XFEM is particularly valuable for problems where the crack path is not known a priori, where cracks propagate or branch, or where remeshing would be prohibitively expensive (3D crack growth). The ability to define a crack with a single `cut_data` line — rather than designing a conforming mesh — dramatically reduces the preprocessing effort.
