# Power BI dashboard notes

The Power BI dashboard itself has not been built; Power BI Desktop is a GUI application, so that step is for Mohamed to do by hand, not something Claude Code can run. This file documents which exported table maps to which recommended visual, so building it is a matter of importing five CSVs and wiring up the visuals below, not re-deriving the data.

| Table | Suggested visual | Notes |
|---|---|---|
| `results/tables/model_performance_summary.csv` | Bar chart: `model` on the axis, `test_accuracy` as the value | The headline chart: shows the 3 single-omics bars next to the integrated DIABLO bar |
| `results/tables/pca_coordinates_by_block.csv` | Scatter plot: `PC1` vs `PC2`, colored by `subtype`, with a slicer on `block` | Lets a viewer flip between mRNA, miRNA, and protein and see how cleanly each one separates subtypes alone |
| `results/tables/diablo_top_features_comp1.csv` | Table or bar chart: `feature`, `loading`, filtered/sliced by `block` | The integrated model's top drivers per omics layer |
| `results/tables/diablo_test_confusion_matrix.csv` | Matrix visual: `true` (rows) vs `predicted` (columns), value `Freq` | The integrated model's held-out test set errors |
| `results/tables/sample_metadata.csv` | Dimension table | Join key (`sample_id`) and a `subtype`/`split` slicer source for the other visuals |

All five tables share `sample_id` or `subtype` as join keys, so they can be related in Power BI's model view rather than needing a manual merge beforehand.
