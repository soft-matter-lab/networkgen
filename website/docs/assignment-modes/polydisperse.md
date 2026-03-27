---
sidebar_position: 4
---

# Polydisperse

The `poly` assignment mode provides fine-grained control over a continuous bond length distribution. Used when `strand_typology = 'polydisperse'`.

---

### `assignment_mode.poly.method`

| Type | Args | Default |
|------|------|---------|
| `string` | `'geom'` \| `'range'` \| `'pmf'` | — |

The method used to generate the polydisperse distribution.

- **geom** — geometric distribution
- **range** — uniform range distribution, uses `uniform.min_value` and `uniform.max_value`
- **pmf** — probability mass function defined by mean and bounds

```matlab
net.architecture.strand_typology.poly.method = 'pmf';
```

---

### `assignment_mode.poly.min_value`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | — |

Minimum allowable bond value in the distribution.

```matlab
net.architecture.strand_typology.poly.min_value = 1;
```

---

### `assignment_mode.poly.align_to_length`

| Type | Args | Default |
|------|------|---------|
| `string` | `'ascend'` \| `'none'` | `'none'` |

When set to `'ascend'`, bonds are sorted and assigned values in ascending order of their length — longer bonds receive larger values.

```matlab
net.architecture.strand_typology.poly.align_to_length = 'ascend';
```

---

### `assignment_mode.poly.rounding`

| Type | Args | Default |
|------|------|---------|
| `string` | `'round'` \| `'ceil'` \| `'floor'` | `'round'` |

Controls how continuous distribution samples are rounded to integer bond values.

```matlab
net.architecture.strand_typology.poly.rounding = 'ceil';
```

---

### `assignment_mode.poly.range_method`

| Type | Args | Default |
|------|------|---------|
| `string` | `'rank'` \| `'linear'` | — |

Controls how values are mapped across the range when `method = 'range'`.

- **rank** — values are assigned by rank ordering
- **linear** — values are assigned by linear interpolation across the range

```matlab
net.architecture.strand_typology.poly.range_method = 'rank';
```

---

### `assignment_mode.poly.target_min`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | — |

Target minimum value for the distribution output.

```matlab
net.architecture.strand_typology.poly.target_min = 2;
```

---

### `assignment_mode.poly.target_max`

| Type | Args | Default |
|------|------|---------|
| `double` | (`target_min`, ∞) | — |

Target maximum value for the distribution output.

```matlab
net.architecture.strand_typology.poly.target_max = 15;
```

---

### PMF sub-settings

The following settings are used when `method = 'pmf'`:

### `assignment_mode.poly.pmf_min`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | — |

Lower cutoff for the PMF distribution.

---

### `assignment_mode.poly.pmf_max`

| Type | Args | Default |
|------|------|---------|
| `double` | (`pmf_min`, ∞) | — |

Upper cutoff for the PMF distribution.

---

### `assignment_mode.poly.pmf_mean`

| Type | Args | Default |
|------|------|---------|
| `double` | (`pmf_min`, `pmf_max`) | — |

Mean of the PMF distribution.

---

### `assignment_mode.poly.pmf_cut_mode`

| Type | Args | Default |
|------|------|---------|
| `string` | `'cap'` | `'cap'` |

Controls how values outside `[pmf_min, pmf_max]` are handled. Currently only `'cap'` is supported, which clamps out-of-range values to the nearest bound.

---

### `assignment_mode.poly.pmf_integerize_rule`

| Type | Args | Default |
|------|------|---------|
| `string` | `'largest_remainder'` | `'largest_remainder'` |

Method used to convert continuous PMF probabilities to integer counts while preserving the total.

---

## Example

```matlab
net.architecture.strand_typology.mode = 'polydisperse';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.poly.method = 'pmf';
net.architecture.strand_typology.poly.pmf_min = 2;
net.architecture.strand_typology.poly.pmf_max = 20;
net.architecture.strand_typology.poly.pmf_mean = 7;
net.architecture.strand_typology.poly.rounding = 'round';
net.architecture.strand_typology.poly.pmf_cut_mode = 'cap';
```
