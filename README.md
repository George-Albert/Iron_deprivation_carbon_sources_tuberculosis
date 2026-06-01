# Iron Deprivation and Carbon-Source Growth Arrest in *Mycobacterium tuberculosis*

This repository contains the R analysis workflow for an RNA-seq study of how
carbon-source availability reshapes the transcriptional response of
*Mycobacterium tuberculosis* H37Rv to iron deprivation during growth arrest.

The project supports a manuscript in preparation on the interaction between
iron limitation, carbon metabolism, and the transition from exponential growth
to stationary phase. The local repository is the canonical working copy for the
analysis.

## Table of Contents

- [Scientific Scope](#scientific-scope)
- [Main Questions](#main-questions)
- [Repository Layout](#repository-layout)
- [Workflow Overview](#workflow-overview)
- [Key Inputs](#key-inputs)
- [Key Outputs](#key-outputs)
- [Reproducibility](#reproducibility)
- [Current Validation Status](#current-validation-status)
- [Traceability](#traceability)
- [Public Release Notes](#public-release-notes)

## Scientific Scope

Iron restriction is a major host-imposed stress during tuberculosis infection.
At the same time, persistent infection environments can expose *M. tuberculosis*
to host-derived lipid carbon sources. This study asks whether lipid-rich carbon
conditions alter, dampen, or redirect the classical iron-starvation response.

The analysis compares cultures grown under four carbon-source regimes:

| Code | Carbon source |
| --- | --- |
| `G-D` | glycerol plus dextrose |
| `G` | glycerol |
| `G-LCFA` | glycerol plus long-chain fatty acids |
| `LCFA` | long-chain fatty acids as the sole carbon source |

Each carbon-source condition is represented under iron-replete and
iron-deprived conditions, across exponential and stationary phases. The lipid
condition includes even long-chain fatty acids such as palmitic, oleic, and
stearic acids.

## Main Questions

The workflow is organized around these analysis questions:

- How does iron deprivation affect the transcriptome in exponential and
  stationary phases?
- How does growth arrest change expression under iron-replete and
  iron-deprived conditions?
- Does the interaction between iron availability and growth arrest depend on
  carbon source?
- Do lipid carbon sources attenuate or rewire the iron-starvation response?
- Which coding genes, non-coding RNAs, GO terms, regulons, and pathway-level
  gene sets explain the observed transcriptional patterns?

## Repository Layout

```text
.
|-- Codes/                       # Ordered R analysis scripts
|-- Docs/                        # Repo notes, script inventory, data review
|-- scripts/                     # Validation helpers
|-- Inputs/                      # Local raw/processed inputs, currently ignored
|-- Outputs/                     # Generated outputs and figures, currently ignored
|-- README.md
|-- .gitignore
`-- .gitattributes
```

Important documentation files:

- `Docs/script_inventory.md`: high-level map of the active R scripts.
- `Docs/publication_data_review.md`: review of local data/output groups for
  public release planning.

## Workflow Overview

The active scripts in `Codes/` are numbered in approximate execution order.
Scripts assume they are run from the project root.

| Stage | Scripts | Purpose |
| --- | --- | --- |
| Input preparation | `001` to `003` | Read comparison, metadata creation, feature annotation cleanup, and 32-sample subset construction. |
| Differential expression | `004` to `007` | DESeq2 modeling, contrast export, normalized data, and statistical summary tables. |
| Core figures | `008` to `014` | DE gene count plots, correlograms, heatmaps, volcano plots, and interaction summaries. |
| Functional interpretation | `015` to `024` | Coherent response classes, non-coding RNA summaries, GO enrichment, route trends, Venn diagrams, and density plots. |
| Final analyses | `025` to `029` | Carbon-source interaction summaries, regulon ORA, gene expression panels, volcano panels, and GO gene-set heatmaps. |
| Helper | `adding_col_to_read.R` | Adds feature biotype metadata to the read table export. |

The current pipeline writes some reproducibility intermediates back into
`Inputs/002_Processed_data/` and manuscript-style outputs into `Outputs/`.
This is intentional for the current script-based workflow and should be kept in
mind when rerunning the analysis.

## Key Inputs

The working local analysis expects:

```text
Inputs/001_Raw_data/
Inputs/002_Processed_data/
```

Major input classes include:

- raw and mapped read-count tables;
- metadata for the 32 selected RNA-seq samples;
- feature annotation tables for coding and non-coding genes;
- DESeq2 contrast metadata;
- GO annotation matrices and enrichment objects;
- regulon annotation files used by the ORA workflow.

The current Git baseline does not track `Inputs/` because the data-release
strategy is still being reviewed for provenance, size, and manuscript timing.
The repository is intended to be public, and there is no privacy objection to
publishing the data. The remaining decision is where each data group belongs:
Git, GitHub Releases, manuscript supplements, or an external archive.

## Key Outputs

Generated outputs are written mainly under:

```text
Outputs/Data/
Outputs/001_Figures_paper/
Outputs/6_Enrichment_GO/
Outputs/ORA_regulons/
Outputs/Final_GO_gene_heatmaps/
```

Representative outputs include:

- DESeq2 contrast tables and Excel exports;
- PCA diagnostics and correlation plots;
- DE gene-count summaries;
- coherent response classifications;
- non-coding RNA summaries;
- GO enrichment tables and heatmaps;
- regulon over-representation results;
- volcano plots and final figure panels.

`Outputs/` is currently ignored by Git. See
`Docs/publication_data_review.md` before deciding which generated outputs should
be versioned directly.

## Reproducibility

### R version

The workflow has been validated locally with R 4.6.0. Any recent R installation
with `Rscript` available on `PATH` should be able to run the validation helper
and the analysis scripts once the required packages are installed.

### Package requirements

The scripts use CRAN and Bioconductor packages including:

```text
ashr, ComplexHeatmap, corrplot, cowplot, DESeq2, dendextend, dendsort,
edgeR, EnhancedVolcano, eulerr, fdrtool, ggplot2, ggdendro, ggrepel,
ggsci, ggthemes, gridExtra, gtools, igraph, limma, locfdr, openxlsx,
patchwork, plotly, qvalue, reshape2, rgl, shadowtext, svglite, tidyverse,
umap, writexl
```

Package versions are not frozen yet. A future improvement would be to add
`renv` or a small dependency preflight script.

### Validate script syntax

From the repository root:

```powershell
Rscript scripts\validate_r_parse.R
```

Expected result:

```text
34 active R scripts parse OK.
```

### Run the workflow

Run scripts from the project root in numeric order:

```powershell
Rscript Codes\001_Corr_new_reads_vs_bowtie_reads.R
Rscript Codes\002_Create_metadata.R
```

Continue through the numbered scripts as needed. Several later scripts depend
on intermediate files produced by earlier scripts, especially `004_DESeq2.R`,
`006_Load_DE_data.R`, `007_Create_stat_data.R`, `014_Figure_3C_bar_plot_enhanced_interactions.R`,
`017_GO_function.R`, `018_enrichments_heatmap_Fig_4A_B.R`, and
`026_ORA_regulons.R`.

## Current Validation Status

The codebase has been reviewed and stabilized in blocks:

- `001` to `003`: input preparation scripts run successfully.
- `004` to `007`: DESeq2/statistical workflow runs successfully.
- `008` to `014`: figure and correlation scripts run successfully.
- `015` to `019`: non-coding and GO analysis scripts run successfully.
- `020`: metabolic route trends now rebuild route gene lists from current GO
  outputs when legacy `Routes_GO` files are absent.
- `021`, `022`, `024`: downstream export, Venn, and density scripts run
  successfully.
- `025` to `029` and `adding_col_to_read.R`: final helper and figure scripts
  run successfully.

Known notes:

- `020_Clusters_genes_trends_per_route.R` detects no pentose phosphate gene set
  from the current GO output sources.
- Some scripts still print package startup messages and diagnostic tables.
- The workflow is script-based rather than a formal pipeline manager.

## Traceability

The repository is organized so that each major result can be traced to a script
range:

| Result type | Main scripts |
| --- | --- |
| Metadata and processed feature tables | `001` to `003` |
| Differential expression contrasts | `004`, `006`, `007`, `021` |
| PCA and QC plots | `005` |
| DE count and interaction summaries | `008`, `012`, `014`, `023`, `025` |
| Correlation and correlogram figures | `009`, `010`, `011` |
| Volcano plots | `013`, `028` |
| Non-coding RNA summaries | `016`, `023`, `adding_col_to_read.R` |
| GO enrichment and heatmaps | `017`, `018`, `019`, `029*` |
| Metabolic route trend summaries | `020` |
| Regulon ORA | `026` |
| Gene-level expression panels | `027` |

For Git-level traceability, the cleanup and stabilization history is preserved
as incremental commits on `main`.

## Public Release Notes

This repository is intended for public scientific sharing. The current tracked
baseline prioritizes code, documentation, validation helpers, and metadata about
publication readiness. Large input/output groups remain ignored until the final
release strategy is chosen.

Before publishing data or generated outputs, review:

```text
Docs/publication_data_review.md
```

The main decision is practical rather than privacy-related: which files should
be committed to Git, which should be attached to a release, and which should be
placed in a manuscript-associated archive.
