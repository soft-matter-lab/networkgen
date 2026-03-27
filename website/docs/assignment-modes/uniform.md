---
sidebar_position: 3
---

# Uniform

The `uniform` assignment mode draws bond values from a flat distribution between a minimum and maximum value.

---

### `assignment_mode.uniform.min_value`

| Type | Args | Default |
|------|------|---------|
| `double` | (-∞, `max_value`] | — |

The lower bound of the uniform distribution.

```matlab
net.architecture.strand_typology.uniform.min_value = 2;
```

---

### `assignment_mode.uniform.max_value`

| Type | Args | Default |
|------|------|---------|
| `double` | [`min_value`, ∞) | — |

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
