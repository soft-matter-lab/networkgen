---
sidebar_position: 2
---

# Getting Started

## Requirements

- MATLAB R2020a or later
- LAMMPS (for running generated files)

## Installation

See the [Installation guide](../install) for full setup instructions.

## Your first network

The minimal configuration to generate a network requires only a domain size and a call to `genNetwork()`:

```matlab
net = Network();
net.Lx = 10;
net.Ly = 10;
net.boundary = 'periodic';
net.isave = true;
genNetwork(net);
```

This creates a 2D random network with monodisperse strands and writes a LAMMPS data file to the current directory.

## Recommended next steps

- Set [Domain & Boundary](./network-setup/domain-boundary) options to define your simulation box
- Choose a [Strand Typology](./architecture/strand-typology) to control bond length distribution
- Configure [Flags & Output](./network-setup/flags-output) to control what files are written
- Explore [Defects](./advanced-features/defects) to introduce voids and damage
