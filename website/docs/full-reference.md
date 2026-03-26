---
sidebar_position: 99
---

# Full Settings Reference

A complete flat reference of all NetworkGen settings. For detailed descriptions and examples see the individual setting pages linked in each section header.

## [Domain & Boundary](./network-setup/domain-boundary)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `b` | `double` | (0, ∞) | `1.0` |
| `dimension` | `double` | `2` \| `3` | `2` |
| `Lx` | `double` | (0, ∞) | — |
| `Ly` | `double` | (0, ∞) | — |
| `Lz` | `double` | (0, ∞) | — |
| `scale` | `double` | (0, ∞) | `1.0` |
| `boundary` | `string` | `'fixed'` \| `'periodic'` | `'fixed'` |
| `seed` | `int` | [1, ∞) | — |
| `write_location` | `string` | any valid path | `'./'` |
| `lammps_data_file` | `string` | any string | `'network'` |
| `lammps_viz_file` | `string` | any string | `'network_viz'` |
| `bond_table_file` | `string` | any string | `'bond_table'` |
| `smp_number` | `int` | [1, ∞) | `1` |

## [Flags & Output](./network-setup/flags-output)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `isave` | `boolean` | `true` \| `false` | `false` |
| `iplot` | `boolean` | `true` \| `false` | `false` |
| `ilog` | `boolean` | `true` \| `false` | `false` |
| `savemode` | `boolean` | `true` \| `false` | `false` |
| `imanualseed` | `boolean` | `true` \| `false` | `false` |
| `idefect` | `boolean` | `true` \| `false` | `false` |
| `ipotential` | `boolean` | `true` \| `false` | `false` |

## [Geometry](./architecture/geometry)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `geometry` | `string` | `'random'` \| `'hex'` | `'random'` |
| `lattice_spacing` | `double` | (0, ∞) | — |
| `spacing_multiplier_mode` | `string` | `'auto'` \| `'manual'` | `'auto'` |
| `spacing_multiplier` | `double` | [0, ∞) | — |
| `lattice_disorder_level` | `double` | [0, 1] | `0` |
| `lattice_disorder_maxfrac` | `double` | [0, 1] | — |
| `lattice_topo_disorder_flag` | `boolean` | `true` \| `false` | `false` |
| `lattice_max_del_per_node` | `int` | [0, `max_peratom_bond`-2] | `0` |
| `lattice_min_degree_keep` | `int` | [3, `max_peratom_bond`] | `3` |
| `rho_atom` | `double` | (0, ∞) | — |

## [Strand Typology](./architecture/strand-typology)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `strand_typology` | `string` | `'mono'` \| `'polydisperse'` \| `'bimodal'` | `'mono'` |

## [Assignment Modes](./assignment-modes/overview)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `assignment_mode.auto` | `boolean` | `true` \| `false` | `true` |
| `assignment_mode.mono.value` | `double` | (-∞, ∞) | — |
| `assignment_mode.uniform.min_value` | `double` | (-∞, `max_value`] | — |
| `assignment_mode.uniform.max_value` | `double` | [`min_value`, ∞) | — |
| `assignment_mode.poly.method` | `string` | `'geom'` \| `'range'` \| `'pmf'` | — |
| `assignment_mode.poly.min_value` | `double` | (0, ∞) | — |
| `assignment_mode.poly.align_to_length` | `string` | `'ascend'` \| `'none'` | `'none'` |
| `assignment_mode.poly.rounding` | `string` | `'round'` \| `'ceil'` \| `'floor'` | `'round'` |
| `assignment_mode.poly.pmf_mean` | `double` | (`pmf_min`, `pmf_max`) | — |
| `assignment_mode.poly.pmf_min` | `double` | (0, ∞) | — |
| `assignment_mode.poly.pmf_max` | `double` | (`pmf_min`, ∞) | — |
| `assignment_mode.bimodal.mean_1` | `double` | (0, `mean_2`) | — |
| `assignment_mode.bimodal.mean_2` | `double` | (`mean_1`, ∞) | — |
| `assignment_mode.bimodal.std_1` | `double` | (0, ∞) | — |
| `assignment_mode.bimodal.std_2` | `double` | (0, ∞) | — |
| `assignment_mode.bimodal.method` | `string` | `'single'` \| `'geom'` \| `'gaussian'` | — |
| `assignment_mode.bimodal.height_mode` | `string` | `'prob'` \| `'fixed'` | — |
| `assignment_mode.bimodal.height_prob` | `double` | [0, 1] | — |
| `assignment_mode.bimodal.height_count` | `int` | [1, ∞) | — |
| `assignment_mode.bimodal.long_first` | `boolean` | `true` \| `false` | `false` |
| `assignment_mode.bimodal.stdR_1` | `double` | (0, ∞) | — |
| `assignment_mode.bimodal.stdR_2` | `double` | (0, ∞) | — |
| `assignment_mode.bimodal.lam_1` | `double` | (0, 1) | — |
| `assignment_mode.bimodal.lam_2` | `double` | (0, 1) | — |
| `assignment_mode.bimodal.double_net_flag` | `boolean` | `true` \| `false` | `false` |
| `assignment_mode.bimodal.alpha` | `double` | (0, 1) | — |
| `assignment_mode.bimodal.auto_1_flag` | `boolean` | `true` \| `false` | `false` |
| `assignment_mode.bimodal.auto_2_flag` | `boolean` | `true` \| `false` | `false` |
| `assignment_mode.bimodal.bin_window_method` | `string` | `'manual'` \| `'adaptive'` | `'adaptive'` |
| `assignment_mode.bimodal.manual_deviation_type` | `string` | `'kuhn'` \| `'both'` | — |

## [Perbond & Potentials](./advanced-features/potentials-perbond)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `perbond.kuhn.auto` | `boolean` | `true` \| `false` | `true` |
| `pot.k_LD` | `double` | (0, ∞) | — |
| `pot.N_rho` | `int` | [1, ∞) | — |
| `pot.rho_min` | `double` | (0, `rho_max`) | — |
| `pot.rho_max` | `double` | (`rho_min`, ∞) | — |

## [Defects](./advanced-features/defects)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `defect.density_mode` | `string` | `'count'` \| `'area_fraction'` | `'count'` |
| `defect.N_voids` | `int` | [0, ∞) | `0` |
| `defect.void_area_frac` | `double` | [0, 1] | — |
| `defect.size_dist` | `string` | `'fixed'` \| `'gaussian'` \| `'exponential'` | `'fixed'` |
| `defect.radius_mean` | `double` | (0, ∞) | — |
| `defect.radius_std` | `double` | (0, ∞) | — |
| `defect.radius_min` | `double` | (0, `radius_max`) | — |
| `defect.radius_max` | `double` | (`radius_min`, ∞) | — |
| `defect.shape_roughness` | `double` | [0, 1] | `0` |
| `defect.shape_n_modes` | `int` | [1, ∞) | — |
| `defect.void_overlap` | `boolean` | `true` \| `false` | `false` |
| `defect.center_distribution` | `string` | `'random'` \| `'uniform'` \| `'clustered'` | `'random'` |
| `defect.n_cluster_parents` | `double` | [1, ∞) | — |
| `defect.cluster_spread` | `double` | (0, ∞) | — |
| `defect.margin_frac` | `double` | [0, 0.5) | — |
| `defect.prune_isolated` | `boolean` | `true` \| `false` | `true` |
| `defect.sparse_network` | `boolean` | `true` \| `false` | `false` |
| `defect.bridge_width` | `double` | (0, ∞) | — |
| `defect.wall_thickness` | `double` | (0, ∞) | — |
| `defect.clamp_thickness` | `double` | (0, ∞) | — |

## [Multi-type Networks](./advanced-features/multi-type)

| Setting | Type | Args | Default |
|---------|------|------|---------|
| `types.natom_type` | `int` | [1, ∞) | `1` |
| `types.nbond_type` | `int` | [1, ∞) | `1` |
| `types.atype_mode` | `string` | `'fixed'` \| `'frac'` | `'fixed'` |
| `types.btype_mode` | `string` | `'fixed'` \| `'frac'` | `'fixed'` |
| `types.atom_count` | `int array` | size [1 x `natom_type`] | — |
| `types.bond_count` | `int array` | size [1 x `nbond_type`] | — |
| `types.atom_frac` | `double array` | size [1 x `natom_type`], sum=1 | — |
| `types.bond_frac` | `double array` | size [1 x `nbond_type`], sum=1 | — |
| `types.connectivity` | `int array` | size [N x 3] | — |
