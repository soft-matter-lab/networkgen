---
sidebar_position: 1
---

# Random Networks - Basics

Random networks are the default NetworkGen workflow.

You use them when you want an amorphous polymer network rather than a regular lattice scaffold.

## Default starting point

The basic random-network setup is:

```matlab
net = network();

net.domain.Lx = 150;
net.domain.Ly = 150;
net.peratom.Max_peratom_bond = 6;

net.architecture.geometry = 'random';
net.architecture.strand_typology.mode = 'mono';

net.generateNetwork();
```

This is the simplest random baseline:

- random node placement
- mono strand typology
- no defects
- no manual perbond override

## The three strand-typology variants

For random networks, the next decision is the strand typology.

### 1. Mono

Use mono when you want the narrowest and simplest target strand distribution.

```matlab
net.architecture.strand_typology.mode = 'mono';
net.architecture.strand_typology.mono.value = 20;
```

### 2. Polydisperse

Use polydisperse when you want a continuous spread of strand lengths.

```matlab
net.architecture.strand_typology.mode = 'polydisperse';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.poly.method = 'pmf';
net.architecture.strand_typology.poly.pmf_mean = 40;
net.architecture.strand_typology.poly.pmf_min = 20;
net.architecture.strand_typology.poly.pmf_max = 120;
```

### 3. Bimodal

Use bimodal when you want two strand populations in the same network.

```matlab
net.architecture.strand_typology.mode = 'bimodal';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.bimodal.method = 'gaussian';
net.architecture.strand_typology.bimodal.mean_1 = 20;
net.architecture.strand_typology.bimodal.mean_2 = 60;
net.architecture.strand_typology.bimodal.std_1 = 4;
net.architecture.strand_typology.bimodal.std_2 = 8;
```

## What to keep fixed while learning

When comparing mono, polydisperse, and bimodal random networks, keep these fixed at first:

- `domain.Lx` and `domain.Ly`
- `peratom.Max_peratom_bond`
- `architecture.rho_atom`
- output flags

That makes it much easier to see what comes from the strand assignment itself.

## Next step

Continue to [Random Networks - Advanced](./advanced) to control the relationship between `strand_typology` and `perbond.kuhn`, especially the difference between automatic and manual perbond assignment.