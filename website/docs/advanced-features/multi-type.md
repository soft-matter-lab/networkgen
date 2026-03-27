---
sidebar_position: 2
---

# Multi-type Networks

NetworkGen supports networks with multiple atom and bond types, allowing heterogeneous networks where different regions or connectivity classes have distinct physical properties in LAMMPS.

---

### `types.natom_type`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | `1` |

Number of distinct atom types in the network.

```matlab
net.architecture.types.natom_type = 2;
```

---

### `types.nbond_type`

| Type | Args | Default |
|------|------|---------|
| `int` | [1, ∞) | `1` |

Number of distinct bond types in the network.

```matlab
net.architecture.types.nbond_type = 2;
```

---

### `types.atype_mode`

| Type | Args | Default |
|------|------|---------|
| `string` | `'fixed'` \| `'frac'` | `'fixed'` |

Controls how atom type counts are specified.

- **fixed** — specify exact counts via `atom_count`
- **frac** — specify fractions via `atom_frac`

```matlab
net.architecture.types.atype_mode = 'frac';
```

---

### `types.btype_mode`

| Type | Args | Default |
|------|------|---------|
| `string` | `'fixed'` \| `'frac'` | `'fixed'` |

Controls how bond type counts are specified.

- **fixed** — specify exact counts via `bond_count`
- **frac** — specify fractions via `bond_frac`

```matlab
net.architecture.types.btype_mode = 'frac';
```

---

### `types.atom_count`

| Type | Args | Default |
|------|------|---------|
| `int array` | size: [1 x `natom_type`] | — |

Array of atom counts per type. Only used when `atype_mode = 'fixed'`.

```matlab
net.architecture.types.atom_count = [100, 50];
```

---

### `types.bond_count`

| Type | Args | Default |
|------|------|---------|
| `int array` | size: [1 x `nbond_type`] | — |

Array of bond counts per type. Only used when `btype_mode = 'fixed'`.

```matlab
net.architecture.types.bond_count = [200, 80];
```

---

### `types.atom_frac`

| Type | Args | Default |
|------|------|---------|
| `double array` | size: [1 x `natom_type`], sum = 1 | — |

Array of atom type fractions. Must sum to 1. Only used when `atype_mode = 'frac'`.

```matlab
net.architecture.types.atom_frac = [0.7, 0.3];
```

---

### `types.bond_frac`

| Type | Args | Default |
|------|------|---------|
| `double array` | size: [1 x `nbond_type`], sum = 1 | — |

Array of bond type fractions. Must sum to 1. Only used when `btype_mode = 'frac'`.

```matlab
net.architecture.types.bond_frac = [0.6, 0.4];
```

---

### `types.connectivity`

| Type | Args | Default |
|------|------|---------|
| `int array` | size: [N x 3] | `[]` (empty) |

An exclusion rule table that restricts which bond types can form between which atom types. Each row is `[bond_type, atom_type_1, atom_type_2]` and defines a forbidden connection.

By default this is empty, meaning there are no restrictions — any bond type can connect any atom types.

:::note Exclusion rules
Each row specifies a connection that is **not allowed**, not one that is required:
- `[1, 1, 2]` — bond type 1 **cannot** form between atom type 1 and atom type 2
- `[1, 1, 1]` — bond type 1 **cannot** bridge two atoms of type 1 (i.e. type 1 bonds are restricted to cross-links only)
:::

```matlab
% Bond type 1 cannot connect atom type 1 to atom type 2
% Bond type 1 cannot connect two atoms of type 1
net.architecture.types.connectivity = [1, 1, 2;
                           1, 1, 1];
```

---

### `types.atype_sel_method`

| Type | Args | Default |
|------|------|---------|
| `string array` | `'random'`, size: [1 x `natom_type`] | `'random'` |

Method used to select which atoms are assigned each type. Currently only `'random'` is supported.

---

### `types.btype_sel_method`

| Type | Args | Default |
|------|------|---------|
| `string array` | `'random'`, size: [1 x `nbond_type`] | `'random'` |

Method used to select which bonds are assigned each type. Currently only `'random'` is supported.

---

## Example

```matlab
net.architecture.types.natom_type = 2;
net.architecture.types.nbond_type = 2;
net.architecture.types.atype_mode = 'frac';
net.architecture.types.btype_mode = 'frac';
net.architecture.types.atom_frac = [0.7, 0.3];
net.architecture.types.bond_frac = [0.6, 0.4];
% No connectivity restrictions — any bond type can connect any atom types
net.architecture.types.connectivity = [];
```
