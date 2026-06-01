# Iron Deprivation and Carbon-Source Growth Arrest

R analysis repository for studying the transcriptional response of
*Mycobacterium tuberculosis* to iron deprivation across growth phase and
carbon-source contexts.

This local repository is intended to become the canonical working version for
the current iron-deprivation/carbon-source manuscript analysis. An older GitHub
repository named `Iron_deprivation_carbon_sources_tuberculosis` exists, but it
is not treated as the source of truth for this cleanup.

## Scientific Scope

The analysis compares transcriptional responses associated with:

- iron availability;
- growth phase or growth-arrest state;
- carbon-source conditions, including glycerol, dextrose/glycerol, LCFA, and
  mixed carbon-source contexts;
- coding and non-coding RNA-associated features;
- enriched, damped, and coherent response patterns across contrasts.

The workflow generates differential expression summaries, PCA and correlogram
diagnostics, volcano plots, GO enrichment summaries, regulon ORA, gene panels,
and manuscript-oriented figure outputs.

## Repository Layout

- `Codes/`: active R scripts, ordered by approximate execution stage.
- `Codes/trash_codes/`: legacy/exploratory scripts, ignored by git.
- `Inputs/`: local input data and processed intermediates, ignored by git.
- `Outputs/`: generated figures, tables, RDS objects, and manuscript artifacts,
  ignored by git.
- `Docs/`: lightweight repository documentation.
- `scripts/`: validation and utility scripts.

## Workflow Overview

The active workflow is script based:

1. compare read tables and build metadata;
2. clean feature annotations and subset the final sample set;
3. run DESeq2 differential expression models;
4. consolidate contrast statistics;
5. generate figure-level analyses and functional summaries;
6. produce GO/regulon enrichment heatmaps and manuscript-ready panels.

See `Docs/script_inventory.md` for the current script inventory.

## Local Data Policy

This repository tracks code and documentation only. Raw inputs, processed
intermediates, generated figures, RDS/RData objects, and root-level manuscript
tables are intentionally excluded from version control.

## Validation

Check that all active R scripts parse:

```r
source("scripts/validate_r_parse.R")
```

On this machine, `Rscript` is available at:

```text
C:/Users/JorgeAlbertoCardenas/AppData/Local/Programs/R/R-4.6.0/bin/Rscript.exe
```

## Reproducibility Notes

- Scripts currently assume they are run from the project root.
- The workflow is order-dependent because several scripts write intermediate
  tables consumed by later scripts.
- Package versions are not yet frozen; adding `renv` or a documented package
  preflight is a recommended future step.
- Path/bootstrap centralization should be handled in a later refactor once this
  initial clean repository baseline is committed.
