# Case 76: 3D Rectangular Waveguide — TE10 Mode Propagation

**MIT 6.635 Lectures:** 3–4 (Green's Functions for Planarly Layered Media)

## Overview

This case simulates the TE10 mode propagating through a 3D rectangular metallic waveguide using the **vector Helmholtz equation** with **Nedelec edge elements** (NEDELEC_ONE). A TE10 mode is injected at one end of the waveguide via a port Robin boundary condition, propagates along the waveguide axis, and is absorbed at the far end by a matched absorbing boundary.

This is the first case in the series to use **vector finite elements** — the natural discretisation for Maxwell's curl-curl equations. Unlike scalar Lagrange elements, Nedelec edge elements live in H(curl) space, with tangential continuity across element faces and allowed normal discontinuities, exactly matching the physics of the electric field.

## The Physics

**Governing equation** (vector Helmholtz, frequency domain):

    curl(curl(E)) - k0^2 * E = 0

**Waveguide parameters:**
- Cross-section: a = 3 m (x) × b = 2 m (y), length L = 10 m (z)
- k0 = 1.5 rad/m → only TE10 propagates (k_c10 = π/3 = 1.047 < k0 < k_c20 = 2π/3 = 2.094)
- Propagation constant: β10 = sqrt(k0² - (π/a)²) = 1.074 rad/m

**TE10 mode profile:** E_y = sin(πx/a) · exp(-jβ10·z), with E_x = E_z = 0.

**Boundary conditions:**
- PEC (perfect electric conductor) on 4 longitudinal walls: tangential E = 0
- Port injection at z = 0: VectorEMRobinBC with TE10 profile
- Absorbing BC at z = L: VectorEMRobinBC (matched termination)

## MOOSE Objects Used

- **CurlCurlField**: curl-curl stiffness ∫ (curl E)·(curl v) dV
- **VectorFunctionReaction** (sign=negative): -k0² mass term -k0² ∫ E·v dV
- **VectorEMRobinBC**: port injection (z=0) and absorbing (z=L) Robin conditions
- **VectorCurlPenaltyDirichletBC**: PEC walls (tangential E = 0 via penalty)
- **VectorVariableComponentAux**: extracts E_y from vector variable for visualisation
- **NEDELEC_ONE** on **HEX20**: first-order edge elements on second-order hexahedra

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case76-3d-waveguide \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case76_3d_waveguide.i 2>&1 | tail -30'
```

## Expected Results

- Converges in 1 Newton iteration (linear system with direct LU solver)
- max |E_y_real| ≈ 0.94, max |E_y_imag| ≈ 0.93 (close to unit amplitude)
- Symmetric transverse profile: E_y at x=a/4 equals E_y at x=3a/4
- Phase propagation along z visible as cos(β·z) in E_y_real and -sin(β·z) in E_y_imag

## Key Takeaways

- **Vector elements are essential** for curl-curl: Nedelec edge elements naturally enforce tangential continuity of E across element boundaries
- **Single-mode operation**: choosing k0 between TE10 and TE20 cutoffs ensures only one propagating mode
- **Port + absorbing BCs**: VectorEMRobinBC can both inject and absorb — same object, different source profile
- **3D visualisation**: isosurfaces of |E_y| in ParaView show the sin(πx/a) transverse profile propagating along z
