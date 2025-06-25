# SymPhas PFC Coefficient File Format

This format is supported for specifying model coefficients in JSON configs using the `"coefficients_file"` key.

## Format Overview
- **CSV-like** (comma or whitespace separated)
- **Section headers**: `[default]`, `[non_coupling]`, `[coupling]`
- **Field count**: `fields: N` at the top
- **Comments**: lines starting with `#`
- **Order**: non-coupling first, then coupling, flattened for the model

## Example (3 fields)
```
# PFC Coefficient File Example
fields: 3

[default]
1.0

[non_coupling]
# alpha, beta, nu, gamma, delta for each field
-0.3, 1.0, 0.0, 0.0, 0.0
-0.3, 1.0, 0.0, 0.0, 0.0
-0.3, 1.0, 0.0, 0.0, 0.0

[coupling]
# alpha, beta, nu, gamma, delta, eps for each unique field pair (i<j)
0.1, 0.2, 0.3, 0.4, 0.5, 1.0
0.1, 0.2, 0.3, 0.4, 0.5, 1.0
0.1, 0.2, 0.3, 0.4, 0.5, 1.0
```

## Usage in JSON
```json
{
  "model": {
    "name": "PFC_CNC",
    "coefficients_file": "my_coeffs.txt"
  }
}
```

## Notes
- The parser will flatten all values in the order: all non-coupling, then all coupling.
- If a value is missing, the `[default]` value is used.
- The number of rows in `[non_coupling]` must match `fields:`.
- The number of rows in `[coupling]` must be `fields * (fields - 1) / 2`.
