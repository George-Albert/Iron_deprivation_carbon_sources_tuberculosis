# Script Inventory

This inventory summarizes the current role of the active scripts in `Codes/`.
Scripts under `Codes/trash_codes/` are treated as legacy/exploratory and are not
part of the tracked workflow by default.

## Main Workflow

| Script range | Current role |
| --- | --- |
| `001`-`003` | Read comparison, metadata creation, feature annotation cleanup, and final 32-sample subset construction. |
| `004`-`007` | DESeq2 differential expression analysis and consolidation of contrast-level statistics. |
| `008`-`014` | Figure-level summaries for DE gene counts, correlograms, heatmaps, volcano plots, and interaction bar plots. |
| `015`-`024` | Coherent response classification, non-coding gene summaries, GO enrichment, response trend summaries, and Fe/carbon-source effect summaries. |
| `025`-`029` | Carbon-source interaction summaries, regulon ORA, figure gene panels, volcano plot panels, and final GO gene-set heatmaps. |

## Technical Notes

- Active scripts currently assume execution from the project root.
- Inputs and generated outputs are local-only and intentionally ignored.
- Several scripts write generated intermediates back into `Inputs/002_Processed_data`; this is part of the current workflow but should be documented before publication.
- A lightweight syntax validation script is available at `scripts/validate_r_parse.R`.
- Package loading is still repeated across scripts; centralization can be a later refactor after the repository baseline is stable.
