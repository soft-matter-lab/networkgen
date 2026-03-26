---
sidebar_position: 2
---

# Mono

The `mono` assignment mode assigns a single fixed value to all bonds. Used when `strand_typology = 'mono'`.

---

### `assignment_mode.mono.value`

| Type | Args | Default |
|------|------|---------|
| `double` | (-∞, ∞) | — |

The fixed value assigned to all bonds. For bond lengths this is typically a positive integer representing the number of Kuhn segments.

```matlab
net.strand_typology = 'mono';
net.assignment_mode.auto = false;
net.assignment_mode.mono.value = 5;
```
