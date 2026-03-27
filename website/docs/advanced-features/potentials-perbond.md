---
sidebar_position: 3
---

# Perbond & Potentials

Per-bond properties define physical characteristics assigned individually to each bond in the network. Each perbond property uses the [Assignment Mode](../assignment-modes/overview) framework to define its distribution — the same system used by `strand_typology`.

---

## Perbond properties

### `perbond.kuhn`

Controls the Kuhn segment count assigned per bond. The Kuhn count determines bond extensibility in the WLC/Langevin potential — bonds with more Kuhn segments are more extensible and less stiff.

`perbond.kuhn` uses the full Assignment Mode framework. Configure it exactly as you would `strand_typology`, using `mono`, `uniform`, `poly`, or `bimodal` sub-settings.

#### `perbond.kuhn.auto`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, copies the exact same distribution used by `strand_typology`. This ensures bond stiffness is consistently paired with bond length across the network.

When `false`, the Kuhn distribution is configured independently using the standard assignment mode sub-settings under `perbond.kuhn`.

:::tip
Leave `auto = true` for most simulations. Only set to `false` if you want to decouple stiffness heterogeneity from length heterogeneity — for example, to isolate the mechanical effect of one variable at a time.
:::

```matlab
% Auto — copies strand_typology distribution exactly
net.perbond.kuhn.auto = true;

% Manual — independent bimodal kuhn distribution
net.perbond.kuhn.auto = false;
net.perbond.kuhn.method = 'bimodal';
net.perbond.kuhn.bimodal.mean_1 = 5;
net.perbond.kuhn.bimodal.mean_2 = 15;
net.perbond.kuhn.bimodal.std_1 = 1.0;
net.perbond.kuhn.bimodal.std_2 = 2.0;
```

---

### `perbond.criticalstretch` *(future)*

Controls the critical stretch assigned per bond — the extension at which a bond fails. Uses the same Assignment Mode framework as `perbond.kuhn`.

```matlab
% Example (once implemented)
net.perbond.criticalstretch.auto = false;
net.perbond.criticalstretch.method = 'gaussian';
net.perbond.criticalstretch.mean = 2.5;
net.perbond.criticalstretch.std = 0.3;
```

---

## Potential & log output

NetworkGen generates two types of force information for use in LAMMPS:

### Cohesion — `pair local/density` potential file

Activated by `ipotential = true`. Writes a bond potential lookup table used with the LAMMPS `pair_style local/density` command. This governs the **cohesive** interaction between bonded pairs.

The settings below control the resolution and range of the generated table.

### Repulsion — equilibrium separation from log file

When `ilog = true`, the log file outputs an **equilibrium separation length** for each bond type. This value should be used as the equilibrium distance parameter for the LAMMPS `pair_style bpm/spring` command, which governs the **repulsive** interaction between atoms.

:::note Two-potential workflow
A typical NetworkGen-based LAMMPS simulation uses both:
- `pair_style bpm/spring` with the equilibrium separation from the log file (repulsion)
- `pair_style local/density` with the generated potential table (cohesion)
:::

---

## Potential table settings

### `pot.k_LD`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | — |

Spring constant for the `local/density` bond potential, in LAMMPS force units.

```matlab
net.pot.k_LD = 1.0;
```

---

### `pot.N_rho`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | — |

Number of discrete points in the potential lookup table. Higher values produce a smoother potential at the cost of a larger table file.

```matlab
net.pot.N_rho = 1000;
```

---

### `pot.rho_min`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, `rho_max`) | — |

Minimum density value in the potential table.

```matlab
net.pot.rho_min = 0.01;
```

---

### `pot.rho_max`

| Type | Args | Default |
|------|------|---------|
| `double` | (`rho_min`, ∞) | — |

Maximum density value in the potential table.

```matlab
net.pot.rho_max = 0.99;
```

---

## Example

```matlab
net.flags.ipotential = true;
net.flags.ilog = true; % outputs equilibrium separation for bpm/spring

% Kuhn auto-coupled to strand_typology
net.perbond.kuhn.auto = true;

% local/density potential table settings
net.pot.k_LD = 1.0;
net.pot.N_rho = 1000;
net.pot.rho_min = 0.01;
net.pot.rho_max = 0.99;
```
