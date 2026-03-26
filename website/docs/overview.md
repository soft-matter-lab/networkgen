---
sidebar_position: 1
---

# Overview

NetworkGen is a MATLAB-based mesoscale polymer network generator that produces simulation-ready input files for [LAMMPS](https://www.lammps.org/). It is designed to give researchers fine-grained control over network topology, strand length distributions, and defect structures. A key objective of the tool is that generated networks are informed from the statistics of the networks themselves, either obtained from experimental data or molecular dynamics simulations. This could be for example to generate a network with strand lengths following the Flory-Schulz distribution.

NetworkGen is not restricted to just polymer networks, and can be used to generate basic spring-networks as well.

## How it works

The workflow follows three steps:

1. **Create** a network object from the `Network` class
2. **Configure** the object's settings to define the network geometry, strand typology, and output options
3. **Generate** the network by passing the object to `genNetwork()`

```matlab
net = network;
net.Lx = 10;
net.Ly = 10;
net.strand_typology = 'polydisperse';
GenerateNetwork(net);
```

:::note
Not all settings need to be modified before generation, as they may be ignored depending on the type of network. Only modify the properties needed. 
:::

## Output files

Depending on your flags, NetworkGen writes:

- **LAMMPS data file** — atom positions, bond connectivity, and type assignments
- **LAMMPS visualization file** — for use with OVITO or VMD
- **Bond table file** — pairwise potential lookup tables
- **Log files** — record of generation parameters and statistics

## Settings structure

Settings are grouped into the following categories:

| Category | Description |
|---|---|
| [Domain & Boundary](./network-setup/domain-boundary) | Simulation box size and boundary conditions |
| [Flags & Output](./network-setup/flags-output) | Control saving, plotting, and logging |
| [Geometry](./architecture/geometry) | Node placement pattern |
| [Strand Typology](./architecture/strand-typology) | Bond length distribution type |
| [Assignment Modes](./assignment-modes/overview) | Distribution shape parameters |
| [Defects](./advanced-features/defects) | Void and damage creation |
| [Multi-type Networks](./advanced-features/multi-type) | Multiple atom and bond types |
| [Potentials & Perbond](./advanced-features/potentials-perbond) | Bond stiffness and potential tables |
