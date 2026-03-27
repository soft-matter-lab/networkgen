---
sidebar_position: 2
---

# Supported Networks

For the tutorial, it is most useful to think about NetworkGen in terms of three practical workflow tracks rather than one long flat list of settings.

## 1. Lattice networks

Lattice networks begin from a regular geometric scaffold. In NetworkGen this means using a hexagonal lattice and then optionally adding positional or topological disorder.

![Lattice network overview](/img/network_arch_main.svg)

Use this track when you want:

- a highly controlled starting geometry
- a clear baseline before introducing disorder
- a network whose topology is easier to reason about visually

The guided page for this track is [Lattice Networks - Basics](./lattice-networks/basics).

## 2. Random networks

Random networks begin from stochastic node placement. This is the default style of NetworkGen and is usually the best starting point for amorphous polymer networks.

![Random network overview](/img/network_gauss_main2.png)

Within the random-network workflow, the main question becomes how you want strand statistics assigned:

- mono
- polydisperse
- bimodal

The tutorial first covers the default/basic setup and then the advanced controls for strand typology and perbond assignment in [Random Networks - Basics](./random-networks/basics) and [Random Networks - Advanced](./random-networks/advanced).

## 3. Defected or disordered networks

Once the base architecture is working, you can add either geometric disorder or explicit defects.

![Defected network overview](/img/defect_network.svg)

This includes two different ideas:

- lattice disorder, where a lattice network is perturbed or partially rewired
- void defects, where material is removed from the network using the `defect` settings

The guided page for this track is [Defects & Disorder](./defects-disorder/overview).

## Which track should you start with?

- Start with lattice networks if you want the cleanest and most interpretable baseline.
- Start with random networks if you want the default NetworkGen workflow and realistic heterogeneous topologies.
- Add defects or disorder only after the base network generates reliably.

## Minimal decision tree

1. Do you want a regular scaffold or a stochastic one?
2. If stochastic, do you want mono, polydisperse, or bimodal strand statistics?
3. Once the base network works, do you want to add disorder or explicit defects?

That is the same order followed by the rest of the tutorial sidebar.