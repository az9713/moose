# Case 57: Xenon-135 Poisoning Transient — Reactor Flux Decay

## Overview

This case simulates one of the most important transient phenomena in nuclear reactor operation: xenon-135 poisoning. Xe-135 is a fission product with an enormous neutron absorption cross section (~2.6 million barns), making it the dominant neutron poison in an operating reactor. The "xenon pit" — a post-shutdown flux suppression that can make restart impossible for several hours — has caused serious operational difficulties in the history of nuclear power, including contributing to the conditions preceding the Chernobyl accident.

The model couples three variables: neutron flux phi, iodine-135 concentration I, and xenon-135 concentration Xe. These are linked through production and decay chains: fission produces I-135 directly (high yield), which decays with a 6.7-hour half-life to Xe-135, which decays with a 9.2-hour half-life while also being burned out by neutron absorption. The interplay of these timescales creates the characteristic xenon oscillation seen in large thermal reactors.

This is MOOSE's first transient case in this nuclear series, introducing the `Transient` executioner and `TimeDerivative` kernels. The three-variable system illustrates how MOOSE handles coupled PDE-ODE systems: the neutron flux has spatial diffusion, while the Xe and I concentrations are purely local (no spatial transport). Time is scaled in hours to produce readable output over the 24-hour simulation period.

## The Physics

### Governing Equations

**Neutron flux (phi):**
```
dphi/dt = -D * laplacian(phi) - Sigma_a * phi + nu_Sf * phi - sigma_Xe_eff * Xe * phi
```

Or in standard form with all terms on the left:
```
dphi/dt + D * laplacian(phi) + Sigma_a * phi - nu_Sf * phi + sigma_Xe_eff * Xe * phi = 0
```

**Iodine-135 (I):**
```
dI/dt = gamma_I * Sigma_f * phi - lambda_I * I
```

I-135 is produced by fission (yield ~6.1%) and destroyed by radioactive decay (half-life 6.7 hr).

**Xenon-135 (Xe):**
```
dXe/dt = lambda_I * I + gamma_Xe * Sigma_f * phi - lambda_Xe * Xe - sigma_Xe_eff * Xe * phi
```

Xe-135 is produced by I-135 decay (dominant source) and direct fission yield (small), and is destroyed by radioactive decay and neutron absorption (burnup).

### Physical Parameters

All rate constants are converted to per-hour units since simulation time is in hours:

| Parameter | SI value | Per-hour value | Physical meaning |
|-----------|----------|----------------|-----------------|
| lambda_I | 2.87e-5 /s | 0.1033 /hr | I-135 decay constant |
| lambda_Xe | 2.09e-5 /s | 0.0752 /hr | Xe-135 decay constant |
| gamma_I | 0.061 | — | I-135 fission yield |
| gamma_Xe | 0.003 | — | Direct Xe-135 fission yield |
| Sigma_f | 0.05 /cm | — | Fission cross section |
| D | 1.0 cm | — | Neutron diffusion coefficient |
| nu_Sigma_f | 0.10 /cm | — | Fission production |
| Sigma_a | 0.08 /cm | — | Base absorption (without Xe) |

The xenon absorption cross section sigma_Xe ~ 2.6e6 barns = 2.6e-18 cm^2 is modeled in an effective linearized form with coefficient sigma_Xe_eff = 0.05 /hr in the flux-Xe coupling term.

### The Xenon Pit Mechanism

Starting from steady state, as the flux evolves:

1. **t = 0 to ~3 hr**: Xe is at initial equilibrium; flux adjusts from initial condition
2. **t = 3 to ~12 hr**: I-135 (already at equilibrium) decays to Xe-135 faster than Xe is burned; Xe concentration rises; flux is increasingly suppressed
3. **t ~ 12 hr**: Peak xenon concentration; flux at minimum (the "xenon pit")
4. **t > 12 hr**: Xe decays away (lambda_Xe = 0.0752 /hr > lambda_I = 0.1033 /hr at low flux); Xe concentration falls; flux recovers toward a new equilibrium

This behavior is characterized by the xenon-iodine lag: the iodine inventory "remembers" the previous power history, and its decay continues to produce xenon even after the flux has dropped.

### Initial Conditions

The simulation starts with a near-equilibrium state:

- phi(x, 0) = cos(pi*(x-10)/20): fundamental cosine mode (normalized, peak = 1 at center)
- Xe(x, 0) = 0.5 * cos(pi*(x-10)/20): Xe at ~50% of maximum equilibrium
- I(x, 0) = 0.3 * cos(pi*(x-10)/20): I at ~30% of maximum equilibrium

The non-exact initial condition (not true equilibrium) allows the transient to develop naturally from the simulation start.

### Domain and Mesh

- Geometry: 1D slab, x in [0, 20] cm, modeled as 2D quasi-1D (40 x 2 elements)
- Only phi has spatial boundary conditions (vacuum at x=0 and x=20)
- Xe and I have no spatial diffusion and no spatial BCs (volumetric point kinetics)

## Input File Walkthrough

### `[Variables]`

```
[phi]    # neutron flux
[]
[Xe]     # Xe-135 concentration (scaled)
[]
[I135]   # I-135 concentration (scaled)
[]
```

Three coupled variables. The Xe and I variables represent scaled concentrations (atoms/cm^3 normalized to O(1) for numerical convenience).

### `[ICs]`

```
[phi_ic]
  type = FunctionIC
  variable = phi
  function = 'cos(3.14159265*(x - 10)/20)'
[]
[Xe_ic]
  type = FunctionIC
  variable = Xe
  function = '0.5 * cos(3.14159265*(x - 10)/20)'
[]
[I_ic]
  type = FunctionIC
  variable = I135
  function = '0.3 * cos(3.14159265*(x - 10)/20)'
[]
```

`FunctionIC` applies a spatially varying initial condition defined by an inline expression. All three variables start with the same cosine spatial shape (consistent with the fundamental mode) but at different amplitudes representing their respective initial concentrations.

### `[Kernels]`

**Neutron flux equation:**

`phi_time` (`TimeDerivative`): d(phi)/dt term. This is what makes the problem transient — without it, we would solve only for steady state.

`phi_diff` (`MatDiffusion`): -D * laplacian(phi) diffusion term.

`phi_abs` (`CoefReaction`, coefficient = 0.08): absorption loss +Sigma_a * phi. `CoefReaction` residual = +coefficient * u * test (unlike `Reaction` which also uses "rate"). Positive coefficient = loss term.

`phi_fission` (`CoefReaction`, coefficient = -0.10): fission source -nu_Sf * phi. Negative coefficient makes this a source (gain).

`phi_xe_poison` (`CoupledForce`, v = Xe, coef = -0.05): Xenon poisoning term. `CoupledForce` residual = -coef * v * test. With coef = -0.05, contribution = +0.05 * Xe * phi * test (a loss term proportional to Xe concentration). As Xe builds up, this term increasingly suppresses phi.

**Iodine-135 equation:**

`I_time` (`TimeDerivative`): d(I)/dt

`I_decay` (`CoefReaction`, coefficient = 0.1033): radioactive decay loss +lambda_I * I = 0.1033 * I /hr

`I_production` (`CoupledForce`, v = phi, coef = 0.00305): fission production. coef = gamma_I * Sigma_f = 0.061 * 0.05 = 0.00305. The `CoupledForce` residual = -coef * phi * test, making this a source (negative residual contribution = gain for I).

**Xenon-135 equation:**

`Xe_time` (`TimeDerivative`): d(Xe)/dt

`Xe_decay` (`CoefReaction`, coefficient = 0.0752): radioactive decay loss +lambda_Xe * Xe

`Xe_from_I` (`CoupledForce`, v = I135, coef = 0.1033): production from I-135 decay. coef = lambda_I = 0.1033. Gives source term lambda_I * I.

`Xe_direct` (`CoupledForce`, v = phi, coef = 0.00015): direct fission yield. coef = gamma_Xe * Sigma_f = 0.003 * 0.05 = 0.00015.

`Xe_burnup` (`CoefReaction`, coefficient = 0.04): linearized neutron burnup of Xe. Approximates the sigma_Xe_eff * phi * Xe term as a linear Xe loss at the average flux level.

### `[Materials]`

```
[neutron_diff]
  type = GenericConstantMaterial
  prop_names = 'D_neutron'
  prop_values = '1.0'
[]
```

The single material property needed by `MatDiffusion` for the neutron flux equation.

### `[BCs]`

Only the neutron flux has boundary conditions (vacuum at x = 0 and x = L). The Xe and I variables have no spatial BCs because they have no diffusion term — they are purely local (ODE-like) variables that happen to be represented on the FEM mesh.

### `[Executioner]`

```
type = Transient
solve_type = 'PJFNK'
petsc_options_iname = '-pc_type'
petsc_options_value = 'lu'

dt = 0.5        # 0.5 hour time steps
end_time = 24.0  # 24 hours total
```

The `Transient` executioner marches the solution forward in time. `dt = 0.5` hours is chosen to resolve the xenon dynamics adequately (lambda_I = 0.1033 /hr implies a timescale of ~10 hr, so 0.5 hr steps give ~20 points per characteristic timescale). A direct LU preconditioner is used for the small coupled system.

### `[Postprocessors]`

Five time-series postprocessors track the transient evolution:

- `phi_max`: Peak flux over the domain (tracks overall power level)
- `phi_center`: Flux at slab center (x=10 cm)
- `Xe_center`: Xenon concentration at center
- `Xe_avg`: Domain-average Xe concentration
- `I_avg`: Domain-average I-135 concentration

The CSV output provides a time history at every time step (48 rows total), showing the xenon buildup and the corresponding flux depression.

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case57-xenon-transient \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case57_xenon_transient.i 2>&1 | tail -30'
```

## Expected Results

The solver performs 48 time steps (0 to 24 hours, dt = 0.5 hr). Console output shows the nonlinear residual converging at each step, typically in 3-5 iterations.

The postprocessor CSV time history shows the classic xenon transient signature:

```
t (hr)  phi_center   Xe_center   I_avg
0.0     1.000        0.500       0.300
3.0     ~0.90        ~0.55       ~0.32
8.0     ~0.70        ~0.65       ~0.35
12.0    ~0.60        ~0.70       ~0.35  (Xe peak / flux minimum)
18.0    ~0.70        ~0.60       ~0.32
24.0    ~0.80        ~0.50       ~0.30  (partial recovery)
```

The key features to observe:
- **Flux depression**: phi_center drops from 1.0 to ~0.6 over 12 hours as xenon builds up
- **Xe peak**: Xe_center peaks around t = 10-12 hours, after which it decays faster than it is produced
- **I-135 lag**: I_avg stays relatively stable since it is continuously replenished by fission at whatever flux exists; its level is roughly proportional to the local flux
- **Partial recovery**: By t = 24 hr, flux recovers toward equilibrium as Xe is depleted, but the system has not fully returned to its initial state

The Exodus output allows 2D visualization at any time step, showing the spatial distribution of all three variables evolving over time.

## Key Takeaways

- The `TimeDerivative` kernel converts a steady-state problem into a transient one; it must be added to every variable that has temporal evolution.
- Xe and I are local (no spatial diffusion) variables embedded in a spatially distributed system — they effectively act as ODEs at each spatial point, coupled to the PDE for neutron flux.
- `CoupledForce` handles source terms from one variable appearing in another variable's equation; sign conventions require care: residual = -coef * v * test, so a positive coef produces a source (negative residual contribution = gain for the variable).
- The 6.7-hour half-life of I-135 and 9.2-hour half-life of Xe-135 determine the timescales of the xenon transient; time step choice must resolve these adequately (dt << 6.7 hr is needed for accuracy).
- The xenon pit is deepest roughly one half-life of Xe-135 (~9 hr) after a power reduction, when the accumulated Xe from ongoing I-135 decay peaks while Xe burnup has dropped with the flux.
- Xenon oscillations can be spatially non-uniform in large reactors, leading to power distribution control challenges — a 1D slab model captures the temporal physics but not the 3D spatial instabilities seen in large PWRs and BWRs.
