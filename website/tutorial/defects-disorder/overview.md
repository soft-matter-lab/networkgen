---
sidebar_position: 1
---

# Adding Defects and Disorder

In NetworkGen, disorder and defects are related but not identical.

## Two ways to move beyond the base network

### 1. Disorder in a lattice network

This modifies the regularity of a lattice network without switching to the explicit defect pipeline.

Typical settings are:

- `lattice_disorder_level`
- `lattice_disorder_maxfrac`
- `lattice_max_del_per_node`
- `lattice_min_degree_keep`

Example:

```matlab
net = network();

net.domain.Lx = 150;
net.domain.Ly = 150;
net.peratom.Max_peratom_bond = 6;

net.architecture.geometry = 'hex_lattice';
net.architecture.lattice_disorder_level = 0.25;
net.architecture.lattice_disorder_maxfrac = 0.30;
net.architecture.lattice_max_del_per_node = 1;
net.architecture.lattice_min_degree_keep = 5;

net.generateNetwork();
```

Use this when you want a network that is still recognizably lattice-based but no longer perfectly regular.

### 2. Explicit void defects

This activates the dedicated defect object and removes material from the network.

The key switch is:

```matlab
net.flags.idefect = true;
```

Then the main settings become:

- `defect.density_mode`
- `defect.n_voids` or `defect.void_area_frac`
- `defect.size_dist`
- `defect.radius_mean`
- `defect.center_distribution`
- `defect.prune_isolated`

Example:

```matlab
net = network();

net.domain.Lx = 150;
net.domain.Ly = 150;
net.peratom.Max_peratom_bond = 6;

net.architecture.geometry = 'random';
net.architecture.strand_typology.mode = 'mono';

net.flags.idefect = true;
net.defect.density_mode = 'count';
net.defect.n_voids = 5;
net.defect.size_dist = 'gaussian';
net.defect.radius_mean = 12;
net.defect.radius_std = 4;
net.defect.center_distribution = 'clustered';
net.defect.n_cluster_parents = 2;

net.generateNetwork();
```

## Recommended order of operations

Do not begin with defects.

Instead:

1. get the base lattice or random network generating cleanly
2. confirm the strand typology and perbond settings are what you want
3. add disorder or explicit defects
4. only then tune advanced defect parameters such as clustering, roughness, overlap, or cleanup rules

## When to use which approach

- Use lattice disorder when you want to perturb a regular scaffold.
- Use explicit defects when you want missing material, void populations, or clustered damage.
- Use both only after each one works separately.