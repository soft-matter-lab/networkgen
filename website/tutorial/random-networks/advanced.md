---
sidebar_position: 2
---

# Random Networks - Advanced

This page covers the most important advanced concept for random networks:

how `architecture.strand_typology` and `perbond.kuhn` interact.

## The rule to remember

When `net.perbond.kuhn.auto = true`, the perbond Kuhn assignment follows the strand-typology assignment.

When `net.perbond.kuhn.auto = false`, NetworkGen reads the settings from the `net.perbond.kuhn` object instead.

That means you can either:

- keep length and Kuhn statistics coupled
- or decouple them and study them independently

## Case 1: Mono strands with automatic Kuhn assignment

This is the simplest coupled configuration.

```matlab
net = network();

net.domain.Lx = 150;
net.domain.Ly = 150;
net.peratom.Max_peratom_bond = 6;

net.architecture.geometry = 'random';
net.architecture.strand_typology.mode = 'mono';
net.architecture.strand_typology.mono.value = 20;

net.perbond.kuhn.auto = true;

net.generateNetwork();
```

Use this when you want the Kuhn assignment to stay consistent with the strand assignment with minimal manual setup.

## Case 2: Mono strands with manual mono Kuhn assignment

Here the network strands stay mono, but the perbond Kuhn distribution is explicitly read from `net.perbond.kuhn`.

```matlab
net = network();

net.domain.Lx = 150;
net.domain.Ly = 150;
net.peratom.Max_peratom_bond = 6;

net.architecture.geometry = 'random';
net.architecture.strand_typology.mode = 'mono';
net.architecture.strand_typology.mono.value = 20;

net.perbond.kuhn.auto = false;
net.perbond.kuhn.mode = 'mono';
net.perbond.kuhn.mono.value = 35;

net.generateNetwork();
```

Use this when you want a single strand length but a different fixed Kuhn assignment.

## Case 3: Mono strands with a distributed Kuhn assignment

This is the first genuinely decoupled example.

```matlab
net = network();

net.domain.Lx = 150;
net.domain.Ly = 150;
net.peratom.Max_peratom_bond = 6;

net.architecture.geometry = 'random';
net.architecture.strand_typology.mode = 'mono';
net.architecture.strand_typology.mono.value = 20;

net.perbond.kuhn.auto = false;
net.perbond.kuhn.mode = 'bimodal';
net.perbond.kuhn.bimodal.method = 'gaussian';
net.perbond.kuhn.bimodal.mean_1 = 15;
net.perbond.kuhn.bimodal.mean_2 = 35;
net.perbond.kuhn.bimodal.std_1 = 3;
net.perbond.kuhn.bimodal.std_2 = 5;

net.generateNetwork();
```

This is useful when you want to isolate the effect of perbond heterogeneity without simultaneously changing the strand-length distribution.

## Case 4: Bimodal strands with automatic Kuhn assignment

```matlab
net.architecture.strand_typology.mode = 'bimodal';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.bimodal.method = 'gaussian';
net.architecture.strand_typology.bimodal.mean_1 = 20;
net.architecture.strand_typology.bimodal.mean_2 = 60;
net.architecture.strand_typology.bimodal.std_1 = 4;
net.architecture.strand_typology.bimodal.std_2 = 8;

net.perbond.kuhn.auto = true;
```

Use this when you want the bimodal structure mirrored by the perbond Kuhn assignment.

## Case 5: Bimodal strands with manual Kuhn assignment

```matlab
net.architecture.strand_typology.mode = 'bimodal';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.bimodal.method = 'gaussian';
net.architecture.strand_typology.bimodal.mean_1 = 20;
net.architecture.strand_typology.bimodal.mean_2 = 60;

net.perbond.kuhn.auto = false;
net.perbond.kuhn.mode = 'mono';
net.perbond.kuhn.mono.value = 25;
```

Use this when you want the network topology to be bimodal while holding Kuhn assignment fixed.

## Practical workflow

For advanced studies, the safest order is:

1. get the random network generating reliably in mono mode
2. switch the strand typology to polydisperse or bimodal
3. decide whether `perbond.kuhn` should follow automatically
4. only then turn `auto` off and introduce a manual perbond distribution

That order makes debugging much easier because only one coupling changes at a time.