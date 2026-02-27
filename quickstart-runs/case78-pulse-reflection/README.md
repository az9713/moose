# Case 78: Time-Domain EM Pulse Reflection from a Dielectric Slab

**MIT 6.635 Lecture:** 7 (Time Domain Method of Moments)

## Overview

This case simulates a **time-domain electromagnetic pulse** propagating through vacuum, hitting a dielectric slab (n = 2), and splitting into reflected and transmitted pulses. It demonstrates Fresnel reflection/transmission, pulse speed reduction inside a dielectric, and time-resolved wave dynamics using the **NewmarkBeta** time integrator with **Nedelec edge elements**.

This is the only time-domain case in the EM series. All other cases use frequency-domain (steady-state) formulations. The time-domain approach directly shows causality — the reflected pulse arrives back at the source only after the forward pulse reaches the slab interface.

## The Physics

**Governing equation** (vector wave equation):

    curl(curl(E)) + (eps_r / c^2) * d^2E/dt^2 = 0

**Domain:** [0, 0.8] x [0, 0.08] m — effectively 1D with a thin 2D strip.

**Dielectric slab:** x in [0.3, 0.5] m with eps_r = 4 (n = 2).

**Gaussian pulse source:** Carrier f0 = 3 GHz, modulated Gaussian envelope with ~4 carrier cycles, injected at x = 0.8 m propagating in -x direction.

**Fresnel coefficients** (normal incidence, n1=1 to n2=2):
- R = -1/3 (reflected amplitude), |R|^2 = 11.1%
- T = 2/3 (transmitted amplitude), transmitted power = 88.9%

**Time integration:** NewmarkBeta (beta=0.25, gamma=0.5) — implicit trapezoidal rule for second-order-in-time PDEs. Unconditionally stable, second-order accurate.

## MOOSE Objects Used

- **CurlCurlField**: curl-curl spatial operator (stiffness)
- **VectorSecondTimeDerivative**: (eps_r/c^2) d^2E/dt^2 mass term
- **VectorTransientAbsorbingBC**: Mur first-order absorbing BC at domain boundaries
- **VectorCurlPenaltyDirichletBC**: Gaussian pulse injection at right boundary
- **NewmarkBeta** time integrator: handles second-order time derivatives
- **NEDELEC_ONE** on **QUAD9**: edge elements for the curl-curl formulation

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case78-pulse-reflection \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case78_pulse_reflection.i 2>&1 | tail -30'
```

## Expected Results

- 200 timesteps with dt = 1e-11 s (total 2 ns simulation)
- Each timestep converges in ~6 Newton iterations (curl-curl + NewmarkBeta coupling)
- Pulse propagates at c in vacuum, slows to c/2 inside the slab
- Reflected pulse amplitude ≈ 1/3 of incident (with π phase flip)
- Transmitted pulse amplitude ≈ 2/3 of incident
- Time snapshots show: approach → interface interaction → split → separation

## Key Takeaways

- **NewmarkBeta** enables implicit time integration of second-order wave equations without reducing to first order
- **Nedelec elements work for time-domain** problems — same spatial discretisation as frequency-domain cases
- **Fresnel coefficients** are automatically captured by the eps_r discontinuity in the FEM formulation
- **Absorbing BCs** (Mur/Robin) prevent artificial reflections from domain boundaries
- **Pulse speed reduction** inside the dielectric is a direct observable consequence of eps_r > 1
