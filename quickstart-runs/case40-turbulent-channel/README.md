# Case 40: Turbulent Channel Flow — RANS k-epsilon Model

**Reference:** Rieutord, *Fluid Dynamics* (Springer, 2015), Ch. 9, Sec. 9.8

---

## Overview

This case simulates fully-developed turbulent flow in a 2D plane channel using
Reynolds-Averaged Navier-Stokes (RANS) with the standard k-epsilon turbulence
closure. It is a canonical benchmark in computational fluid dynamics: the
time-averaged velocity profile follows the well-known **log-law of the wall**,
and the k-epsilon model is specifically designed to reproduce this behavior.

The channel has half-width H = 1, total width 2H = 2, and length L = 30. Flow
is driven by specifying a bulk inlet velocity of U_bulk = 1. At Re = 13700 the
flow is firmly in the turbulent regime and the viscous sublayer is thin enough
that wall functions are both necessary and accurate.

---

## Physics: Reynolds-Averaged Navier-Stokes

Turbulent flows contain a vast range of length and time scales. Direct numerical
simulation resolves every eddy but is prohibitively expensive except at low Re.
RANS instead decomposes every flow quantity into a time-mean part and a
fluctuating part:

```
u_i = <u_i> + u_i'
```

Substituting into the Navier-Stokes equations and averaging eliminates the
fluctuating terms except for the **Reynolds stress tensor**:

```
tau_ij^R = -rho * <u_i' * u_j'>
```

This term represents the momentum transfer by turbulent eddies. It must be
modeled to close the system.

---

## k-epsilon Turbulence Model

The standard k-epsilon model (Launder & Spalding, 1974) uses the **Boussinesq
eddy-viscosity hypothesis**: the Reynolds stresses are proportional to the
mean strain rate, just like laminar viscous stresses, but with an enhanced
**turbulent viscosity** mu_t:

```
tau_ij^R = 2 * mu_t * S_ij  -  (2/3) * rho * k * delta_ij
```

where S_ij = (1/2)(du_i/dx_j + du_j/dx_i) is the mean strain rate tensor and
k is the turbulent kinetic energy (TKE).

The turbulent viscosity is:

```
mu_t = rho * C_mu * k^2 / epsilon
```

Two additional transport equations close the model.

### TKE Equation (k)

```
d(rho*k)/dt + div(rho*u*k) = div[(mu + mu_t/sigma_k) * grad(k)] + P_k - rho*epsilon
```

where P_k = mu_t * |grad(u)|^2 is the turbulence production rate (energy
extracted from the mean flow by Reynolds stresses).

### Dissipation Rate Equation (epsilon)

```
d(rho*eps)/dt + div(rho*u*eps) = div[(mu + mu_t/sigma_eps) * grad(eps)]
                                 + C1_eps * (eps/k) * P_k
                                 - C2_eps * rho * eps^2 / k
```

The first source term represents production of epsilon proportional to turbulence
production, and the second represents the cascade and dissipation of epsilon.

### Model Constants

| Constant | Value | Role |
|----------|-------|------|
| C_mu | 0.09 | Turbulent viscosity coefficient |
| C1_eps | 1.44 | Epsilon production coefficient |
| C2_eps | 1.92 | Epsilon destruction coefficient |
| sigma_k | 1.0 | Turbulent Prandtl number for k |
| sigma_eps | 1.3 | Turbulent Prandtl number for epsilon |

These constants are calibrated against a wide range of turbulent shear flows.

---

## Wall Functions and the Log-Law

### The Law of the Wall

Near a solid wall, turbulent flow develops a universal structure in non-dimensional
(wall) coordinates:

```
y+ = y * u_tau / nu        (wall-normal distance in viscous units)
u+ = u / u_tau             (velocity in friction-velocity units)
```

where u_tau = sqrt(tau_w / rho) is the friction velocity and tau_w is the
wall shear stress.

Three layers exist:

| Region | y+ range | Profile |
|--------|----------|---------|
| Viscous sublayer | 0 < y+ < 5 | u+ = y+ (linear) |
| Buffer layer | 5 < y+ < 30 | transition |
| Log-law region | 30 < y+ < 300 | u+ = (1/kappa)*ln(y+) + B |

The log-law constants are kappa = 0.41 (von Karman constant) and B = 5.2.

### Why Wall Functions Are Needed

The k-epsilon model is derived under the assumption of fully turbulent flow. In
the viscous sublayer the turbulence is damped by viscosity and the model breaks
down. Rather than refine the mesh to resolve y+ ~ 1 (which would require
millions of cells at Re = 13700), **wall functions** bridge the solution from
the first interior cell to the wall using the known log-law profile.

This case uses `wall_treatment = 'eq_newton'`: the equilibrium Newton wall
function, which iteratively solves for u_tau using the log-law, then uses
u_tau to set consistent boundary conditions on k and epsilon:

```
k_wall = u_tau^2 / sqrt(C_mu)
eps_wall = u_tau^3 / (kappa * y_P)
```

where y_P is the distance from the wall to the first cell center.

---

## SIMPLE Segregated Solver

The SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm
solves the coupled velocity-pressure system by a sequence of linear solves
rather than one monolithic nonlinear solve:

```
1. Solve momentum equations for u*, v* (ignoring pressure correction)
2. Form the pressure-correction equation (Poisson-like)
3. Solve for pressure correction p'
4. Update velocities: u = u* - (1/a_P) * grad(p')
5. Update pressure: p = p* + alpha_p * p'
6. Solve turbulence equations (k, epsilon)
7. Update mu_t from k and epsilon
8. Repeat until convergence
```

Under-relaxation is essential for stability:
- Momentum relaxation: 0.7 (typical range 0.5-0.8)
- Pressure relaxation: 0.3 (typical range 0.2-0.4)
- Turbulence relaxation: 0.2 (conservative; k-eps can be stiff near walls)

MOOSE implements SIMPLE via the `LinearFV` infrastructure:
`MooseLinearVariableFVReal` variables are solved with `LinearFVKernels`,
which assemble standard sparse linear systems. Each physical quantity gets
its own named solver system (`u_system`, `v_system`, `pressure_system`,
`TKE_system`, `TKED_system`), allowing different preconditioners per equation.

---

## Problem Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Re | 13700 | Based on bulk velocity and channel full width 2H |
| rho | 1.0 | Fluid density |
| U_bulk | 1.0 | Bulk inlet velocity |
| mu | 1.4599e-4 | = rho * U_bulk * 2H / Re |
| H | 1.0 | Channel half-width |
| L | 30.0 | Channel length (ensures fully developed flow) |
| Mesh | 10 x 10 cells | Two blocks stitched, biased toward walls |
| Wall bias | 0.7 | Cell size ratio between successive layers |

### Initial Conditions

The turbulence quantities are initialized from the turbulence intensity
I = 0.16 * Re^(-1/8), which gives:

```
k_init = 1.5 * (I * U_bulk)^2
eps_init = C_mu^(3/4) * k_init^(3/2) / (2*H)
mu_t_init = rho * C_mu * k_init^2 / eps_init
```

---

## Mesh Design

The mesh uses two `GeneratedMeshGenerator` blocks stitched together with
`StitchMeshGenerator`. Block 1 covers y in [0, H] with `bias_y = 0.7`
(cells shrink toward y = 0, i.e., the centerline). Block 2 covers y in
[-H, 0] with `bias_y = 1/0.7` (cells shrink toward y = 0 from below).

After stitching, the bottom/top boundaries of the two blocks merge into the
internal interface, and the original bottom/top of the combined domain become
the wall boundaries. The bias concentrates cells near both walls where
velocity gradients are steep and wall functions need accurate y_P values.

---

## Running the Case

### Docker (recommended)

```bash
docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs/case40-turbulent-channel:/work" \
  -w /work \
  idaholab/moose:latest \
  combined-opt -i case40_turbulent_channel.i
```

### Expected Runtime

Approximately 2-5 minutes for 500 SIMPLE iterations on a coarse 10x10 mesh.

---

## Expected Results

### Velocity Profile

After convergence the x-velocity profile at the outlet cross-section should
show the characteristic turbulent channel shape:

- **Near the wall (y+ < 5):** thin viscous sublayer, steep gradient
- **Log-law region (30 < y+ < 300):** u+ = (1/0.41)*ln(y+) + 5.2
- **Centerline:** u_centerline ~ 1.2 * U_bulk = 1.2

The turbulent profile is much flatter than the parabolic laminar Poiseuille
profile (where u_centerline = 1.5 * U_bulk), because turbulent mixing
redistributes momentum more uniformly across the channel.

### Postprocessor Output

| Quantity | Expected Value | Physical Meaning |
|----------|---------------|------------------|
| vel_x_centerline | ~1.18 - 1.22 | Peak velocity at x=27, y=0 |
| avg_vel_x | ~1.0 | Bulk velocity conservation |
| max_TKE | ~0.003 - 0.01 | Peak turbulent kinetic energy |
| max_mu_t | ~0.01 - 0.05 | Peak turbulent viscosity |

### Turbulence Fields

- **TKE** peaks near the wall (y+ ~ 15) where production P_k = mu_t * |grad(u)|^2
  is maximum due to the steep mean velocity gradient.
- **Epsilon** peaks very close to the wall and decays rapidly into the bulk.
- **mu_t** peaks in the log-law region and is roughly uniform in the core;
  mu_t / mu ~ 10-100 confirms the turbulence-dominated momentum transport.

---

## Validation Against the Log-Law

To verify the simulation against the log-law, extract the wall-normal velocity
profile at x = 27 (90% of channel length, where flow is fully developed):

1. Compute friction velocity: u_tau = sqrt(tau_w / rho)
   where tau_w = mu * du/dy|_wall
2. Compute y+ = y_P * u_tau / nu for each cell center
3. Compute u+ = u / u_tau
4. Plot u+ vs. log(y+) and compare to the line u+ = (1/0.41)*ln(y+) + 5.2

For this coarse mesh (10x10 cells), the first wall cell center sits at
y_P ~ 0.15 (after bias), giving y+ ~ 100-200 — solidly in the log-law region.
The wall function should place the near-wall velocity directly on the log-law.

---

## Key Takeaways

### RANS Modeling Philosophy

RANS models do not simulate turbulence; they model its *statistical effect* on
the mean flow. The k-epsilon model captures the essential physics — production,
transport, and dissipation of turbulent kinetic energy — in two equations.
This makes it orders of magnitude cheaper than DNS or LES while retaining
engineering accuracy for attached, fully-developed flows.

### Limitations of k-epsilon

The standard k-epsilon model works well for:
- Fully turbulent, attached shear flows (pipes, channels, boundary layers)
- Flows far from separation

It struggles with:
- Strongly curved streamlines (centrifugal effects)
- Flows with adverse pressure gradients near separation
- Buoyancy-driven flows (requires additional buoyancy terms)
- Transitional flows (laminar-to-turbulent)

More advanced models (k-omega SST, RSM, LES) address some of these limitations
at higher computational cost.

### Wall Function Sensitivity

The accuracy of wall functions depends on y+. The log-law region (30 < y+ < 300)
is where standard wall functions are valid. If cells are too fine (y+ < 5), the
wall function is applied in the viscous sublayer where it is incorrect. If cells
are too coarse (y+ > 300), the cell center is outside the log-law region.
Checking y+ values is a critical step in any RANS simulation quality assessment.

### SIMPLE Convergence

The SIMPLE algorithm is conditionally stable; convergence depends on under-
relaxation factors. The values used here (0.7/0.3/0.2) are conservative choices
that sacrifice convergence speed for robustness. In a production setting you
might use 0.8/0.4/0.5 on a well-resolved mesh after verifying stability.

The `continue_on_max_its = true` flag allows the run to complete even if the
500-iteration limit is reached before the residual tolerance is met, making it
easier to examine intermediate results during development.

---

## Connection to Prior Cases

| Case | Topic | Connection |
|------|-------|------------|
| Case 07 | Laminar Poiseuille flow | Same geometry; laminar limit Re << 2300 |
| Case 08 | Navier-Stokes lid-driven cavity | Incompressible NS, SIMPLE |
| Cases 11-12 | FV Navier-Stokes | LinearFV infrastructure, RhieChow |
| Case 40 | **Turbulent channel (RANS k-eps)** | Full turbulence closure |

The progression from laminar Poiseuille (parabolic profile, exact analytical
solution) to turbulent channel (log-law profile, RANS modeling) illustrates
both the richness of turbulent physics and the engineering pragmatism of
closure modeling.
