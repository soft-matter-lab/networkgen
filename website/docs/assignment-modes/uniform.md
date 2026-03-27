---
custom_edit_url: null
sidebar_position: 3
---

# Uniform

The `uniform` assignment mode draws bond values from a flat distribution between a minimum and maximum value.

---

### `assignment_mode.uniform.min_value`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `5` |

The lower bound of the uniform distribution.

```matlab
net.architecture.strand_typology.uniform.min_value = 2;
```

---

### `assignment_mode.uniform.max_value`

| Type | Args | Default |
|------|------|---------|
| `double` | (0, ∞) | `40` |

The upper bound of the uniform distribution.

```matlab
net.architecture.strand_typology.uniform.max_value = 10;
```

---

## Example

```matlab
net.architecture.strand_typology.mode = 'polydisperse';
net.architecture.strand_typology.auto = false;
net.architecture.strand_typology.poly.method = 'range';
net.architecture.strand_typology.uniform.min_value = 3;
net.architecture.strand_typology.uniform.max_value = 12;
```
