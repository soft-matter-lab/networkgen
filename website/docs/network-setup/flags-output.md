---
custom_edit_url: null
sidebar_position: 2
---

# Flags & Output

These boolean settings control what NetworkGen saves, displays, and logs during network generation.

---

### `isave`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, writes the LAMMPS data file to `write_location`. This is the primary output file used to run LAMMPS simulations.

```matlab
net.flags.isave = true;
```

---

### `iplot`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, opens a MATLAB figure showing a visualization of the generated network. Useful for quick inspection during development but should be disabled for batch generation.

:::tip Batch runs
Set `iplot = false` when generating many networks in a loop — rendering plots significantly slows down generation.
:::

```matlab
net.flags.iplot = true;
```

---

### `ilog`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, writes log files containing generation parameters and network statistics to `write_location`.

```matlab
net.flags.ilog = true;
```

---

### `savemode`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, automatically appends `smp_number` to output filenames to prevent overwriting. Recommended when generating multiple network realizations.

```matlab
net.flags.savemode = true;
net.smp_number = 1;
```

---

### `imanualseed`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, uses the value of `seed` as the random number generator seed. When `false`, a random seed is used each run.

```matlab
net.flags.imanualseed = true;
net.seed = 42;
```

---

### `idefect`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, activates defect (void) creation in the network. Requires the **Defect** settings to be configured.

:::note
See [Defects](../advanced-features/defects) for full configuration options.
:::

```matlab
net.flags.idefect = true;
net.defect.n_voids = 5;
```

---

### `ipotential`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, writes bond potential lookup tables to a file. Required when using custom non-linear bond potentials in LAMMPS.

:::note
See [Potentials & Perbond](../advanced-features/potentials-perbond) for full configuration options.
:::

```matlab
net.flags.ipotential = true;
```
