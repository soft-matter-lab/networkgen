---
custom_edit_url: null
sidebar_position: 3
---

# Peratom

The `peratom` settings control per-atom bonding constraints, defining the maximum number of bonds any single atom can form.

---

### `max_peratom_bond`

| Type | Args | Default |
|------|------|---------|
| `int` | [3, ∞) | — |

The maximum number of bonds allowed per atom. This sets a hard cap on node connectivity during network generation. It also defines the upper bound for `lattice_min_degree_keep` and `lattice_max_del_per_node` in the geometry settings.

:::tip Choosing a value
For 2D random networks, a value of `4` to `6` is typical. For hexagonal lattices the natural connectivity is `6`, so setting this to `6` or higher is recommended. Setting it too low can cause the generator to fail when placing bonds in dense regions.
:::

```matlab
net.peratom.Max_peratom_bond = 6;
```
