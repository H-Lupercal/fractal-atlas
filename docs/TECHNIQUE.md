# Draft method — Fractal Atlas v0.1

This framework maps **numerical data** to parameterized Julia-set portraits.

## Input

- **Entity features:** one unique row ID and selected numerical columns.
- **Supplied relationships:** source, target, score edge list, optionally directional.

Feature data are normalized and mapped to pairwise similarities. These do **not** constitute biological interaction measurements.

## Hierarchy

A relationship matrix is divided into quadrants recursively, including uneven dimensions. Each rectangle is addressable by row and column bounds. At a terminal cell the tool displays exact source pair IDs and stored score.

## Rendering

For a region R, calculate mean (m) and standard deviation (s) across finite values. Use a fixed, bounded demonstration mapping:

- real(c) = -0.745 + 0.17 tanh(m)
- imaginary(c) = 0.186 + 0.17 tanh(s)

Render Julia iterations: z[n+1] = z[n]² + c, with escape radius 4. Escape iteration determines pixel intensity.

These transformations are **artistic choices**, not a claim of inherent fractality or biological significance. Alternative parameter mappings, normalization choices and color palettes may produce different images.

## Two different zoom operations

**Data zoom** selects a smaller rectangle; it progressively reaches individual original relationships.

**Mathematical zoom** enlarges the generated Julia image, exposing new mathematical detail without adding biological data.

## Reproducibility

The configuration determines normalization, feature selection, relationship similarity, and iteration count. Preserve dataset version, row ordering, input schema and configuration when comparing outputs.

## Limitations

This is an exploratory prototype. A Julia image is not invertible back to all measurements, so drill-down exposes the underlying data. Rendering cost increases with output resolution and iteration count; dense matrices cost O(N²) storage. Current simple renderer rerenders entirely on interaction.
