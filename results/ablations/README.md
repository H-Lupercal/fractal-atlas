# Input-size experiments

All images use real data and the same Julia rendering rule.
For the **feature-count** test we hold the 180 mouse cortex samples fixed and introduce more measured protein features.
For the **entity-count** test we hold six maize expression features fixed and grow nested, reproducible gene subsets.
The results test visual sensitivity, **not** improved information recovery or biology.

| Experiment | Portrait | Entities | Features | Gray entropy (bits) | Source |
|---|---|---:|---:|---:|---|
| feature-count | [Mouse cortex: 8 features, 180 observations](mouse_features_08.png) | 180 | 8 | 3.7517 | [Source](https://archive.ics.uci.edu/dataset/342/mice+protein+expression) |
| feature-count | [Mouse cortex: 24 features, 180 observations](mouse_features_24.png) | 180 | 24 | 3.7696 | [Source](https://archive.ics.uci.edu/dataset/342/mice+protein+expression) |
| feature-count | [Mouse cortex: 55 features, 180 observations](mouse_features_55.png) | 180 | 55 | 3.7947 | [Source](https://archive.ics.uci.edu/dataset/342/mice+protein+expression) |
| feature-count | [Mouse cortex: 77 features, 180 observations](mouse_features_77.png) | 180 | 77 | 3.802 | [Source](https://archive.ics.uci.edu/dataset/342/mice+protein+expression) |
| entity-count | [Maize: 64 genes, six expression features](maize_genes_064.png) | 64 | 6 | 3.7028 | [Source](https://github.com/corn2code/TWAStutorial/blob/main/merged_gene_tpms.csv) |
| entity-count | [Maize: 128 genes, six expression features](maize_genes_128.png) | 128 | 6 | 3.6997 | [Source](https://github.com/corn2code/TWAStutorial/blob/main/merged_gene_tpms.csv) |
| entity-count | [Maize: 180 genes, six expression features](maize_genes_180.png) | 180 | 6 | 3.707 | [Source](https://github.com/corn2code/TWAStutorial/blob/main/merged_gene_tpms.csv) |

## Reproduce

    python -m pip install -e '.[gallery]'
    python scripts/run_ablations.py

Only derived images and metadata are committed, never the third-party raw datasets.
Feature selection is deterministic but not a controlled biological experiment.
The order of rows and feature selection can affect results, as can normalization.
