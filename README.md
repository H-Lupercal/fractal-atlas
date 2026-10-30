# Fractal Atlas

A dataset-agnostic Python toolkit for translating high-dimensional numerical data into **interactive Julia-set fractal landscapes**. Genetics inspired the project; biological data is **not required**.

**Research question:** How can we visualize the complexity of genetic data through a single mathematical construct?

## Features

- Import numeric CSV/TSV tables, optional Parquet, or measured/predicted pairwise relationships.
- Choose identifiers, columns, missing-value treatment, normalization, similarity metric, and computational limits in YAML or JSON.
- Build a traceable N × N relationship atlas with recursive quadrants.
- Navigate the **data hierarchy** independently of **mathematical Julia zoom**.
- Inspect original entity identifiers, measurements, and pair scores.
- Export a self-contained offline HTML explorer and provenance summary.
- Run reproducible synthetic and real-data demonstration datasets.

## Quick start

Requires Python 3.10 or newer.

~~~bash
python -m pip install -e '.[preview]'
fractal-atlas demo -o examples/generated --entities 256 --features 12
fractal-atlas build -c examples/generated/config.json -o outputs/synthetic
~~~

Open outputs/synthetic/index.html in your browser.

For your own data, create a YAML configuration:

~~~yaml
name: My data
data:
  path: example.csv
  kind: features
  id_column: id
  normalization: zscore
  missing: median
  similarity: cosine
  max_cells: 4000000
render:
  max_iter: 112
  feature_influence: 0.075
~~~

Run:

~~~bash
fractal-atlas build -c config.yml -o outputs/mine
~~~

For measured relationships set kind to relationships and configure the source, target, and value column names. Missing relationship measurements remain missing.

## Important distinctions

A similarity calculated from features is **not** an experimentally observed genetic interaction. Likewise, fractal zoom reveals increasingly fine **mathematical** details, while data zoom moves toward more specific **measurements**. The visualization is an artistic encoding, **not** evidence of inherent fractality or scientific discovery.

Exported explorers embed the data used to generate them. **Never publish confidential data inside a public explorer.**

## Documentation

See docs/TECHNIQUE.md for the mapping, docs/REAL_DATA.md for biological examples, and docs/VALIDATION.md for test limits. The software is a **v0.1 research prototype** under the MIT license.
