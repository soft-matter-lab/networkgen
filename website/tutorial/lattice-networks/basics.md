---
sidebar_position: 1
---

# Lattice Networks - Basics

This is the simplest controlled tutorial path in NetworkGen.

The goal is to generate a hexagonal lattice network before adding any disorder or defects.

## Core idea

The defining switch is the geometry:

```matlab
net.architecture.geometry = 'hex_lattice';
```

From there, the most important lattice settings are:

- `lattice_spacing`
- `spacing_multiplier_mode`
- `lattice_disorder_level`
- `lattice_disorder_maxfrac`
- `lattice_max_del_per_node`
- `lattice_min_degree_keep`

For a basic lattice tutorial run, keep the disorder settings mild or leave them at their defaults.

## Minimal example

```matlab
net = network();

net.domain.Lx = 150;
net.domain.Ly = 150;
net.peratom.Max_peratom_bond = 6;

net.architecture.geometry = 'hex_lattice';
net.architecture.lattice_spacing = 6;
net.architecture.spacing_multiplier_mode = 'auto';

net.architecture.strand_typology.mode = 'mono';

net.generateNetwork();
```

## What this gives you

- a regular lattice-based network
- mono strand typology
- the cleanest starting point for testing boundary conditions, output files, and post-processing

## Recommended next changes

Once this works, vary only one setting group at a time:

1. boundary and domain size
2. lattice spacing
3. lattice disorder settings
4. strand typology

## When to move on

Move to [Random Networks - Basics](../random-networks/basics) if you no longer want a regular scaffold.

Move to [Defects & Disorder](../defects-disorder/overview) if you want to keep the lattice base but introduce disorder or voids.