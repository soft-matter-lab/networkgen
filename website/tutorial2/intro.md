---
sidebar_position: 1
---

# Tutorial overview

These tutorials walk through complete working examples of how to use NetworkGen to build different types of polymer networks. Each tutorial covers the full MATLAB workflow from creating the network object to generating the output files.

## Before you start

Make sure NetworkGen is installed and on your MATLAB path. See the [Installation guide](/docs/install) if you haven't done this yet.

---

## Lattice networks

Lattice networks place nodes on a regular geometric pattern. They are useful for controlled studies where you want a well-defined starting topology.

| Tutorial | Description |
|---|---|
| [Perfect lattice](./lattice/perfect-lattice) | Generate a hexagonal lattice network with optional positional and topological disorder |

---

## Random networks — Basic

Random networks place nodes stochastically. These tutorials cover the three core strand typology modes.

| Tutorial | Description |
|---|---|
| [Mono](./random-basic/mono) | All strands attempt equal length — the simplest random network |
| [Polydisperse](./random-basic/poly) | Strands drawn from a continuous length distribution |
| [Bimodal](./random-basic/bimodal) | Two distinct strand populations mixed in a single network |

---

## Random networks — Advanced

More complex configurations combining multiple features.

| Tutorial | Description |
|---|---|
| [Double networks](./random-advanced/double-networks) | Interpenetrating bimodal networks with controlled mesh ratio |
| [Mixing assignment modes](./random-advanced/mixing-modes) | Decouple strand typology from perbond properties for independent control |
| [Advanced polydisperse](./random-advanced/advanced-poly) | Fine-grained PMF control, alignment, and rounding strategies |
