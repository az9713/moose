# Case 96: Maxwell's Capacitor — Two-Layer Dielectric with Interfacial Polarization

## Overview

Maxwell's capacitor is the canonical model for interfacial polarization: what happens when a step voltage is applied across two dielectric layers that have different permittivities and different conductivities. At the moment the voltage is applied (t = 0+), the system behaves as two capacitors in series — charge distributes according to the dielectric constants. At steady state (t >> tau) it behaves as two resistors in series — charge distributes according to the conductivities. During the transient, free charge accumulates at the interface between the layers, building from zero up to a maximum that drives the final resistive partition.

This is the physical origin of "interfacial polarization" (also called Maxwell-Wagner-Sillars polarization), a loss mechanism that appears in high-voltage cable insulation, multilayer ceramic capacitors, biological cell membranes, and polymer composites. MIT 6.641 Lecture 7 derives the exact analytic expression for the interfacial surface charge density and the characteristic time constant tau = (epsilon_a d_b + epsilon_b d_a) / (sigma_a d_b + sigma_b d_a).

The MOOSE implementation uses a two-block 1D mesh with material properties assigned per block. `StitchMeshGenerator` connects the two independently generated half-domains at the shared interface without requiring any explicit interface conditions — continuity of flux is enforced automatically by the Galerkin finite-element method through the shared nodes at x = 0.5. The problem exercises the `ADMatDiffusion` kernel with a spatially varying diffusivity (sigma/epsilon) drawn from the material system.

---

## The Physics

**Governing equation**

In each layer, the current-continuity equation combined with Gauss's law gives:

```
epsilon * d^2 Phi / (dx dt) + sigma * d^2 Phi / dx^2 = 0
```

This can be rewritten as a diffusion equation for the potential:

```
dPhi/dt = (sigma / epsilon) * d^2 Phi / dx^2
```

with layer-dependent diffusivity D = sigma / epsilon.

**Analytic time constant**

```
tau = (epsilon_a * d_b + epsilon_b * d_a) / (sigma_a * d_b + sigma_b * d_a)
```

With the normalized parameters used here (d_a = d_b = 0.5, epsilon_a = 2, epsilon_b = 5, sigma_a = 1, sigma_b = 10):

```
tau = (2 * 0.5 + 5 * 0.5) / (1 * 0.5 + 10 * 0.5) = 3.5 / 5.5 ≈ 0.636 s
```

**Interfacial surface charge**

At steady state the surface charge at the interface is:

```
sigma_s(inf) = epsilon_0 * V * (epsilon_b * sigma_a - epsilon_a * sigma_b) / (sigma_a * d_b + sigma_b * d_a)
```

The sign of sigma_s depends on whether epsilon_b * sigma_a exceeds epsilon_a * sigma_b — that is, whether the more permittive layer is also the less conductive one.

**Boundary conditions**

| Boundary | Condition |
|----------|-----------|
| Anode (x = 0) | Phi = 1.0 V (step voltage applied at t = 0) |
| Cathode (x = 1.0) | Phi = 0 V (grounded) |

**Initial condition**

Phi = 0 everywhere. The step voltage is implemented by applying the Dirichlet BC at t = 0 with initial_condition = 0 in the variable block — MOOSE applies the BC from the first step onward.

**Domain and discretization**

- 1D domain [0, 1]: layer_a from x = 0 to 0.5, layer_b from x = 0.5 to 1.0
- 50 elements per layer, 100 total

**Material properties**

| Layer | Block | sigma / epsilon (D) | Physical interpretation |
|-------|-------|----------------------|------------------------|
| a | 0 | 0.5 s^-1 | Low conductivity, moderate permittivity |
| b | 1 | 2.0 s^-1 | Higher conductivity, higher permittivity |

---

## Input File Walkthrough

**[Mesh]**

The mesh is assembled from three generators:

1. `GeneratedMeshGenerator` creates layer_a (0 to 0.5) and layer_b (0.5 to 1.0) separately.
2. `SubdomainIDGenerator` assigns block IDs 0 and 1 to the two halves.
3. `StitchMeshGenerator` merges them at the shared boundary (`right` of layer_a stitched to `left` of layer_b).
4. `RenameBoundaryGenerator` gives the outer boundaries readable names: `anode` and `cathode`.

The stitch operation creates coincident nodes at x = 0.5 that are merged into a single shared node, automatically enforcing continuity of Phi and normal flux across the interface.

**[Variables]**

`phi` is the electric potential. It starts at zero (initial_condition = 0) and is driven to the step-voltage profile by the Dirichlet BCs.

**[Kernels]**

```
[phi_time]       type = ADTimeDerivative   — dPhi/dt
[phi_conduction] type = ADMatDiffusion     — (sigma/epsilon) * d^2 Phi/dx^2
                        diffusivity = sigma_over_eps
```

`ADMatDiffusion` looks up the material property named `sigma_over_eps` from the `[Materials]` system at each quadrature point, so it automatically uses 0.5 in layer_a and 2.0 in layer_b without any conditional logic in the kernel.

**[Materials]**

```
[layer_a_props]  ADGenericConstantMaterial  sigma_over_eps = 0.5   block = 0
[layer_b_props]  ADGenericConstantMaterial  sigma_over_eps = 2.0   block = 1
```

Block-restricted materials are the standard MOOSE pattern for heterogeneous domains. The `block = N` parameter confines each material object to its assigned subdomain.

**[Postprocessors]**

| Name | Measures |
|------|---------|
| phi_interface | Phi at x = 0.5 — tracks the transient voltage redistribution |
| avg_phi_a | Mean potential in layer a |
| avg_phi_b | Mean potential in layer b |

The interface potential starts near the capacitive divider value (epsilon_b/(epsilon_a + epsilon_b) = 5/7 ≈ 0.714) and relaxes toward the resistive divider value (sigma_b/(sigma_a + sigma_b) = 10/11 ≈ 0.909).

**[Executioner]**

Transient NEWTON with LU direct factorization. Time step dt = 0.02 s, running to t = 3.2 s (five relaxation times, 5 * 0.636 = 3.18 s), which is enough to see the transient reach 99.3% of steady state.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case96-maxwells-capacitor \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case96_maxwells_capacitor.i 2>&1 | tail -30'
```

The 160-step transient runs in a few seconds. The Exodus output shows the 1D potential profile evolving from a straight line (immediate ohmic response) toward the resistive partition. The CSV tracks the three postprocessors versus time.

---

## Expected Results

At t = 0+ the potential profile is set by the initial condition; the step voltage builds up over the first time step. The interface potential phi_interface should relax exponentially:

```
phi_interface(t) = phi_ss - (phi_ss - phi_0) * exp(-t / tau)
```

where phi_ss ≈ 0.909 (resistive divider), phi_0 ≈ 0.714 (capacitive divider), and tau ≈ 0.636 s.

| Time (s) | phi_interface (approx.) |
|----------|------------------------|
| 0.0 | 0.714 |
| 0.636 | 0.818 |
| 1.272 | 0.872 |
| 3.2 | 0.907 |

In the Exodus output, plotting the potential profile at several times shows a kink at x = 0.5 whose slope changes as the interfacial charge builds up — the electric field (slope) is different on the two sides, and the discontinuity grows over the relaxation time.

---

## Key Takeaways

- Maxwell's capacitor is the archetype for interfacial polarization: mismatched epsilon and sigma cause free charge to accumulate at internal interfaces under applied voltage.
- The transient time constant tau = (epsilon_a d_b + epsilon_b d_a) / (sigma_a d_b + sigma_b d_a) is a weighted average of the two layer relaxation times.
- `StitchMeshGenerator` joins independently constructed mesh blocks at a shared boundary, with coincident-node merging automatically enforcing interface continuity.
- `SubdomainIDGenerator` + block-restricted `Materials` is the standard MOOSE pattern for heterogeneous material domains.
- `ADMatDiffusion` with a material-system diffusivity automatically interpolates the correct sigma/epsilon in each element, with no conditional logic needed in the kernel.
- The capacitive-to-resistive voltage redistribution is a practical engineering concern: high-voltage insulation must be designed so that the steady-state electric field (resistive partition) does not exceed the breakdown threshold in either layer.
