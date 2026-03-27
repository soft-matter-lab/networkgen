---
sidebar_position: 1
---

# Assignment Modes Overview

Assignment modes define statistical distributions used throughout NetworkGen. They are not specific to bond length — any property that requires a distribution uses this same framework. Currently this includes:

- **`strand_typology`** — controls the target bond length distribution
- **`perbond.kuhn`** — controls the Kuhn segment count distribution
- **`perbond.criticalstretch`** — controls the critical stretch distribution *(future)*

Each of these properties can independently be assigned any of the available distribution modes. The modes and their sub-settings are described in the pages below.

---

### `assignment_mode.auto`

| Type | Args | Default |
|------|------|---------|
| `boolean` | `'true'` \| `'false'` | `'true'` |

When `true`, NetworkGen automatically selects and configures the assignment mode to match the parent property. For example, when configuring `perbond.kuhn`, setting `auto = true` causes it to copy the exact same distribution used by `strand_typology`.

When `false`, you must manually configure the sub-settings for the chosen mode.

:::tip Recommended for most users
Leave `auto = true` unless you need an explicit decoupling between properties — for example, to study the effect of stiffness heterogeneity independently of length heterogeneity.
:::

```matlab
net.architecture.strand_typology.auto = true;
```

---

## Available modes

| Mode | Description |
|---|---|
| [Mono](./mono) | Single fixed value |
| [Uniform](./uniform) | Flat distribution between min and max |
| [Polydisperse](./polydisperse) | Continuous distribution with fine control |
| [Bimodal](./bimodal) | Mixture of two distributions |

---

:::note Top-level settings required
Whichever mode is chosen, the top-level distribution statistics (`mean`, `std`, bounds) must always be configured. The advanced sub-settings within each mode are optional and primarily relevant for `strand_typology` and `perbond.kuhn` when generating double networks. See each mode's page for details.
:::
