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
net.assignment_mode.uniform.min_value = 2;
```

---

### `assignment_mode.uniform.max_value`

| Type | Args | Default |
|------|------|---------|
| `double` | [`min_value`, ∞) | — |

The upper bound of the uniform distribution.

```matlab
net.assignment_mode.uniform.max_value = 10;
```

---

## Example

```matlab
net.strand_typology = 'polydisperse';
net.assignment_mode.auto = false;
net.assignment_mode.poly.method = 'range';
net.assignment_mode.uniform.min_value = 3;
net.assignment_mode.uniform.max_value = 12;
```
