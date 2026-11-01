# Real-data fractal gallery

Only computed images and metadata are stored here, not redistributed raw source datasets.

## Reproduce

    python -m pip install -e '.[gallery]'
    python scripts/run_gallery.py

All portraits use the same preprocessing and quadratic Julia parameters. Results are illustrative visual comparisons, not biological inference.

| Group | Dataset | Entities | Features | Mean gray | Entropy (bits) | Source |
|---|---|---:|---:|---:|---:|---|
| Plants | [Iris — setosa](iris_setosa.png) | 50 | 4 | 16.27 | 3.7263 | [Source](https://archive.ics.uci.edu/dataset/53/iris) |
| Plants | [Iris — versicolor](iris_versicolor.png) | 50 | 4 | 15.6268 | 3.7008 | [Source](https://archive.ics.uci.edu/dataset/53/iris) |
| Plants | [Iris — virginica](iris_virginica.png) | 50 | 4 | 15.7838 | 3.7065 | [Source](https://archive.ics.uci.edu/dataset/53/iris) |
| Human | [Human breast nuclei — malignant](human_malignant.png) | 180 | 30 | 17.7346 | 3.7783 | [Source](https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic) |
| Human | [Human breast nuclei — benign](human_benign.png) | 180 | 30 | 18.464 | 3.8018 | [Source](https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic) |
| Plants | [Maize RNA expression — lower](maize_lower.png) | 180 | 6 | 16.9323 | 3.7484 | [Source](https://github.com/corn2code/TWAStutorial/blob/main/merged_gene_tpms.csv) |
| Plants | [Maize RNA expression — higher](maize_higher.png) | 180 | 6 | 15.7278 | 3.707 | [Source](https://github.com/corn2code/TWAStutorial/blob/main/merged_gene_tpms.csv) |
| Animals | [Mouse cortex proteins — control](mouse_control.png) | 180 | 77 | 18.5071 | 3.802 | [Source](https://archive.ics.uci.edu/dataset/342/mice+protein+expression) |
| Animals | [Mouse cortex proteins — trisomy](mouse_trisomy.png) | 180 | 77 | 18.4652 | 3.8015 | [Source](https://archive.ics.uci.edu/dataset/342/mice+protein+expression) |
| Bacteria | [E. coli iModulon activity — all](ecoli_all.png) | 92 | 278 | 19.8507 | 3.8374 | [Source](https://github.com/SBRG/precise-db) |
| Bacteria | [E. coli iModulon activity — iron](ecoli_iron.png) | 92 | 8 | 16.3493 | 3.7314 | [Source](https://github.com/SBRG/precise-db) |
| Viruses | [SARS-CoV-2 reference — Wuhan/Hu-1/2019](sars_ref_1.png) | 185 | 16 | 19.5378 | 3.8319 | [Source](https://github.com/nextstrain/ncov/blob/master/data/references_sequences.fasta) |
| Viruses | [SARS-CoV-2 reference — 21L](sars_ref_2.png) | 186 | 16 | 19.5258 | 3.8317 | [Source](https://github.com/nextstrain/ncov/blob/master/data/references_sequences.fasta) |

## Interpretation

The image differences in pairwise_visual_differences.csv compare pixel intensities with fixed iteration count and palette. These are **not** biological comparisons.
Input types include plant flower morphology, plant RNA expression, human nucleus morphology, mouse cortex protein abundance, bacterial iModulon activity, and viral nucleotide-window composition.
Source matrices may have different numbers of features and records. Viral genomes are not cells. Feature similarity is not experimentally measured genetic interaction.
The square Julia-set previews are only single-region portraits; the generic Fractal Atlas CLI supports full data and mathematical navigation.

See [interactive gallery](../../docs/gallery/index.html), [montage](contact_sheet.png), [manifest](manifest.json), and [pairwise metrics](pairwise_visual_differences.csv).
