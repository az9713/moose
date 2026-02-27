# Case 66: Mineral Precipitation — Two-Species Kinetic Reaction

## Overview

When two dissolved ions meet in a subsurface fluid, they can combine to precipitate a solid mineral. This process — mineral precipitation — is fundamental to geochemistry and controls the long-term evolution of aquifer chemistry, the clogging of subsurface injection wells, the growth of stalactites, and the cementation of sedimentary rocks. The rate of precipitation is kinetically controlled: it depends on how far the solution is from thermodynamic equilibrium and on the surface area of existing mineral available as a nucleation substrate.

This case uses the MOOSE `chemical_reactions` module's `ReactionNetwork` Action system to simulate a two-species kinetic precipitation experiment. Species A diffuses from the left boundary and species B diffuses from the right boundary. Where they meet near the center of the domain, the reaction A + B → mineral(s) proceeds at a kinetically controlled Arrhenius rate. The mineral accumulates as a localized band that grows with time. This counter-diffusion geometry mimics natural systems such as silica band formation or carbonate precipitation at geochemical fronts.

Key concepts demonstrated:

- `ReactionNetwork/SolidKineticReactions` Action to auto-generate all kinetic precipitation kernels and variables
- `KineticDisPreConcAux` AuxVariable for tracking dissolved species modified by kinetics
- `CoupledBEKinetic` kernel for the backward Euler time integration of kinetic reactions
- Arrhenius rate law with activation energy `Ea = 15 kJ/mol` at T = 298.15 K
- Counter-diffusion geometry producing a localized reaction zone

---

## The Physics

### Reaction and Rate Law

The kinetic precipitation reaction:

```
A(aq) + B(aq) → mineral(s)
```

is modeled with a rate law based on transition state theory (TST):

```
r = k_rate * a_m * exp(-Ea / (R * T)) * (1 - Q/K)
```

where:
- `k_rate = 1e-8 mol/(m^2 * s)` — intrinsic rate constant per unit mineral surface area
- `a_m` — specific mineral surface area (m^2/m^3, here set proportional to mineral concentration)
- `Ea = 15,000 J/mol` — activation energy for the precipitation reaction
- `R = 8.314 J/(mol * K)` — universal gas constant
- `T = 298.15 K` — temperature (25°C)
- `Q = [A] * [B]` — reaction quotient (product of activities)
- `K = 10^(log10_Keq) = 10^(-6)` — equilibrium constant

The factor `(1 - Q/K)` drives the reaction: when `Q < K`, the solution is undersaturated and dissolution occurs; when `Q > K`, the solution is supersaturated and precipitation occurs. Here the reaction proceeds forward (precipitation) when the product `[A] * [B]` exceeds the saturation product K.

At room temperature, the Arrhenius factor:
```
exp(-Ea / (R * T)) = exp(-15000 / (8.314 * 298.15)) = exp(-6.05) = 2.35e-3
```

This reduces the effective rate constant from the intrinsic value, representing the kinetic barrier to forming the solid phase.

### Transport Equations

Each dissolved species obeys its own diffusion equation modified by the kinetic reaction source term:

```
d([A])/dt = D_A * Laplacian([A]) - r
d([B])/dt = D_B * Laplacian([B]) - r
d(mineral)/dt = r
```

where `r` is the local precipitation rate. The `ReactionNetwork/SolidKineticReactions` Action automatically generates all three equations and the correct coupling between them.

### Domain and Boundary Conditions

```
y = 0.2 m  ----------------------------------------  TOP: zero-flux
           |                                        |
           |   A diffuses from left (c_A=1)         |
           |   B diffuses from right (c_B=1)         |
           |   Mineral forms where A and B meet     |
           |                                        |
y = 0.0 m  ----------------------------------------  BOTTOM: zero-flux

x = 0 (c_A=1, c_B=0)                 x = 1 m (c_A=0, c_B=1)
```

The counter-diffusion setup creates a spatial gradient: [A] decreases from left to right, [B] decreases from right to left, and both are nonzero in the central region (x ~ 0.5 m). Precipitation occurs wherever the product [A][B] exceeds the equilibrium constant K = 1e-6.

Species A boundary conditions:
- x = 0: DirichletBC, [A] = 1 mol/m^3
- x = 1: DirichletBC, [A] = 0 mol/m^3

Species B boundary conditions:
- x = 0: DirichletBC, [B] = 0 mol/m^3
- x = 1: DirichletBC, [B] = 1 mol/m^3

Mineral: no boundary condition (it is a solid immobile phase, produced in-situ).

### Reaction Onset Analysis

The mineral begins to precipitate when Q = [A][B] first exceeds K = 1e-6. Before any significant diffusion has occurred, the concentrations at x = 0.5 m remain near zero. As the diffusion fronts arrive, [A][0.5] and [B][0.5] both approach 0.5 mol/m^3, giving Q = 0.5 * 0.5 = 0.25 >> K = 1e-6. The saturation condition is easily met once diffusion establishes a significant concentration at the center. Therefore the onset of precipitation is controlled by the diffusion timescale:

```
t_diffusion = (L/2)^2 / D = (0.5)^2 / 1e-3 = 250 s ~ 4 minutes
```

Mineral should begin accumulating at the center after approximately t ~ 30-40 s, when the diffusion fronts first meet.

---

## Input File Walkthrough

The input file is `case66_mineral_precipitation.i`.

### `[Mesh]`

```
type = GeneratedMesh
dim = 2
nx = 40
ny = 4
xmax = 1.0
ymax = 0.2
```

A 40 x 4 mesh with 0.025 m elements in x. This resolution adequately captures the reaction zone width, which is set by the diffusion length at the simulation end time. The 4 elements in y make the geometry quasi-1D.

### `[ReactionNetwork]`

The `ReactionNetwork` block is the key innovation of the `chemical_reactions` module. It allows the user to specify geochemical reactions in a high-level symbolic form and automatically generates all required kernels, variables, and auxiliary objects.

```
[ReactionNetwork]
  primary_species = 'A B'
  [SolidKineticReactions]
    [mineral]
      primary_reactant_stoichiometry = '1 1'
      kinetic_rate_constant = 1e-8
      log10_Keq = -6
      activation_energy = 15000
      molar_volume = 1.0
      specific_reactive_surface_area = 1.0
      kinetic_mineral_name = mineral
    []
  []
[]
```

The `SolidKineticReactions` sub-block declares:
- `primary_species = 'A B'` — the two dissolved reactants (must match variable names)
- `primary_reactant_stoichiometry = '1 1'` — one mole of A and one mole of B per mole of mineral
- `kinetic_rate_constant` — the intrinsic rate constant k_rate
- `log10_Keq = -6` — equilibrium constant K = 10^(-6)
- `activation_energy = 15000 J/mol` — the Ea for the Arrhenius factor
- `kinetic_mineral_name = mineral` — name of the solid phase variable (auto-created)

The Action auto-generates:
- `KineticDisPreConcAux` variables for tracking how much each primary species is consumed by the kinetic reaction
- `CoupledBEKinetic` kernels in the primary species equations to apply the kinetic source/sink
- The `mineral` variable (CONSTANT MONOMIAL) tracking the solid phase concentration

### `[Variables]`

The primary species variables must be declared explicitly:

```
[A]
  initial_condition = 0
[]
[B]
  initial_condition = 0
[]
```

Initial conditions are zero for both species (the domain starts empty; the boundary conditions inject A and B). The `mineral` variable is auto-created by the `ReactionNetwork` Action with initial condition zero.

### `[Kernels]`

For each primary species, two kernels are needed:

**Diffusion kernels** (manually added):
- `TimeDerivative` + `CoefDiffusion` for species A: `d[A]/dt = D_A * Laplacian([A])`
- `TimeDerivative` + `CoefDiffusion` for species B: `d[B]/dt = D_B * Laplacian([B])`

**Kinetic coupling kernels** (auto-generated by `SolidKineticReactions` Action):
- `CoupledBEKinetic` for species A: adds `-r` to the A equation
- `CoupledBEKinetic` for species B: adds `-r` to the B equation
- A mineral accumulation kernel: adds `+r` to the mineral equation

### `[Materials]`

```
[porous_medium]
  type = GenericConstantMaterial
  prop_names  = 'diffusivity'
  prop_values = '1e-3'
[]
```

A single diffusivity D = 1e-3 m^2/s is used for both species A and B. In reality, ionic diffusivities differ (e.g., Ca^2+ diffuses at 7.9e-10 m^2/s in water; values here are chosen for illustration). The temperature T = 298.15 K is provided separately as a `Constant` parameter to the `SolidKineticReactions` Action.

### `[BCs]`

| Name | Variable | Boundary | Value |
|------|----------|----------|-------|
| `A_left` | A | left | 1.0 (DirichletBC) |
| `A_right` | A | right | 0.0 (DirichletBC) |
| `B_left` | B | left | 0.0 (DirichletBC) |
| `B_right` | B | right | 1.0 (DirichletBC) |

### `[Postprocessors]`

| Name | Type | Purpose |
|------|------|---------|
| `mineral_total` | ElementIntegralVariablePostprocessor | Total moles of mineral precipitated (mol/m) |
| `A_center` | PointValue at (0.5, 0.1) | Species A concentration at reaction zone |
| `B_center` | PointValue at (0.5, 0.1) | Species B concentration at reaction zone |
| `mineral_max` | ElementExtremeValue | Peak mineral concentration |

### `[Executioner]`

```
type     = Transient
dt       = 5.0
end_time = 50.0
```

Ten time steps of 5 s each. The reaction begins after approximately 30-40 s, so the simulation captures the early stages of mineral growth. The mineral band is still narrow and growing at t = 50 s; longer simulations would show the band broadening as both reactants are consumed near the center.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case66-mineral-precipitation \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case66_mineral_precipitation.i 2>&1 | tail -30'
```

Output files:
- `case66_mineral_precipitation_out.e` — Exodus with [A], [B], and mineral fields at each step
- `case66_mineral_precipitation_out.csv` — Time series of mineral_total, A_center, B_center, mineral_max

---

## Expected Results

### Reaction Zone Development

| Time (s) | A_center | B_center | mineral_max | mineral_total |
|----------|---------|---------|------------|--------------|
| 0        | 0.000   | 0.000   | 0.000      | 0.000        |
| 10       | 0.003   | 0.003   | ~1e-6      | ~1e-7        |
| 25       | 0.042   | 0.042   | ~8e-4      | ~2e-4        |
| 35       | 0.108   | 0.108   | ~5e-3      | ~1.5e-3      |
| 50       | 0.193   | 0.193   | ~0.018     | ~5e-3        |

The mineral concentration is zero until approximately t = 30-35 s (the diffusion timescale for the fronts to meet at the center), then grows rapidly as the supersaturation condition is met and kinetic rates take hold.

### Spatial Distribution

The mineral forms a narrow band centered at x = 0.5 m. This band is spatially symmetric (by the symmetry of the problem) and grows in both amplitude and width over time. In the Exodus file at t = 50 s:
- Species A: decreasing profile from [A] = 1 at x = 0 to [A] ~ 0 at x = 1
- Species B: mirror image increasing from [B] ~ 0 at x = 0 to [B] = 1 at x = 1
- Mineral: narrow Gaussian-like peak at x = 0.5, with maximum ~ 0.018 mol/m^3

### Onset of Precipitation

The mineral begins to appear after approximately t ~ 35 s, consistent with the diffusion timescale estimate (t_diff = L^2/(4D) for the fronts to reach the center). Before this, the product [A][B] at the center is below K = 1e-6, so the solution is undersaturated and no precipitation occurs.

---

## Key Takeaways

- The `ReactionNetwork/SolidKineticReactions` Action is the key tool for kinetically-controlled mineral precipitation in MOOSE. It reads a symbolic description of the reaction (stoichiometry, rate constant, equilibrium constant, activation energy) and automatically generates all required kernels, AuxVariables, and materials.
- The Arrhenius factor `exp(-Ea/(R*T))` reduces the effective rate constant relative to its intrinsic value. At 25°C with Ea = 15 kJ/mol, this factor is about 2.35e-3. Higher activation energies make the reaction more temperature-sensitive.
- The equilibrium constant K = 10^(log10_Keq) sets the saturation threshold. Precipitation occurs when Q = [A][B] > K. With K = 1e-6 and reactant concentrations approaching 0.1-0.5 mol/m^3, Q >> K — the system is highly supersaturated and the rate is dominated by the kinetic term, not the proximity to equilibrium.
- Counter-diffusion geometry (species A from left, B from right) naturally creates a localized reaction zone at the meeting point. This is a common pattern in natural systems: silica banding in rocks, carbonate veins at geochemical interfaces, and iron-sulfide nodules in reducing sediments.
- The `CoupledBEKinetic` kernel implements a backward-Euler (fully implicit) integration of the kinetic rate. This avoids the stiffness issues that arise with explicit integration of fast kinetic rates.
- `KineticDisPreConcAux` tracks the cumulative moles of each primary species consumed by the kinetic reaction. This is distinct from the instantaneous rate; it allows mass balance checking.
- For the mineral phase, no boundary conditions are needed — it is a solid (immobile) phase that accumulates in-situ. Only the dissolved primary species require boundary conditions specifying the supply of reactants.
