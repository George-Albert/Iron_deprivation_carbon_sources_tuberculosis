# Iron Deprivation and Carbon-Source Growth Arrest in *Mycobacterium tuberculosis*

This repository contains the R analysis code for an RNA-seq study of how carbon
source availability reshapes the transcriptional response of *Mycobacterium
tuberculosis* H37Rv to iron deprivation during the transition from exponential
growth to stationary phase.

The project supports a manuscript in preparation on the transcriptomic effect of
lipid carbon sources during in vitro iron deprivation. The central biological
question is whether lipid-rich environments, used here as an in vitro proxy for
conditions encountered during long-lasting infection, modify or dampen the
classical iron-starvation response of *M. tuberculosis*.

## Scientific Context

Iron restriction is a major host-imposed stress during tuberculosis infection.
At the same time, chronic and latent infection environments can be enriched in
host-derived lipids. This analysis combines those two environmental dimensions:
iron availability and carbon-source composition.

The study compares axenic cultures grown with four carbon-source profiles:

- `G_D`: glycerol plus dextrose;
- `G`: glycerol;
- `G_L`: glycerol plus long-chain fatty acids;
- `L`: long-chain fatty acids as the sole carbon source.

For each carbon-source context, the analysis compares iron-replete and
iron-deprived cultures across exponential and stationary phases. The lipid
condition uses even long-chain fatty acids, including palmitic, oleic, and
stearic acids.

## Analysis Goals

The scripts are organized to address several connected questions:

- how iron deprivation affects the transcriptome at exponential and stationary
  phases;
- how growth arrest changes expression under iron-rich and iron-depleted
  conditions;
- how the interaction between iron availability and growth arrest depends on
  carbon source;
- whether lipid availability reduces, redirects, or rewires the global
  iron-starvation response;
- which coding and non-coding RNA features, functional categories, GO terms,
  regulons, and pathway-level gene sets explain those differences.

Major outputs include differential expression tables, PCA diagnostics,
correlograms, volcano plots, functional-category summaries, coherent-response
classifications, non-coding RNA summaries, GO enrichment analyses, regulon
over-representation analyses, and manuscript-oriented figure panels.

## Experimental and Statistical Overview

The working analysis uses 32 RNA-seq samples from *M. tuberculosis* H37Rv
cultures spanning carbon source, iron availability, and growth phase. Reads were
mapped against the *M. tuberculosis* H37Rv ASM19595v2 reference, with custom
feature annotation used for coding and non-coding RNA analyses.

Differential expression is modeled in R with DESeq2. The main design considers:

1. iron deprivation effects during exponential phase;
2. iron deprivation effects during stationary phase;
3. growth-arrest effects under iron-replete conditions;
4. growth-arrest effects under iron-deprived conditions;
5. the interaction between iron availability and growth arrest.

The repository currently preserves the existing script-based workflow before
deeper refactoring. The next development steps are to review individual scripts,
centralize repeated setup code, and improve reproducibility without changing the
scientific results.

## Repository Layout

- `Codes/`: active R scripts, ordered by approximate execution stage.
- `Codes/trash_codes/`: legacy or exploratory scripts kept local and ignored by
  git.
- `Inputs/`: local raw inputs and processed intermediates.
- `Outputs/`: generated figures, tables, RDS objects, and manuscript artifacts.
- `Docs/`: repository notes, script inventory, and publication-data review.
- `scripts/`: validation and utility scripts.

See `Docs/script_inventory.md` for a current map of the active scripts.

## Data and Outputs

The current public commit tracks code and documentation. Data and generated
outputs are available locally and are intended to be reviewed for publication
readiness step by step.

Because these data are associated with a manuscript in preparation, the decision
to add files to git, attach them as GitHub release assets, or deposit them in an
external repository should be made deliberately. See
`Docs/publication_data_review.md` for the current inventory of local-only data
and output groups.

## Validation

Check that all active R scripts parse:

```r
source("scripts/validate_r_parse.R")
```

On this machine, `Rscript` is available at:

```text
C:/Users/JorgeAlbertoCardenas/AppData/Local/Programs/R/R-4.6.0/bin/Rscript.exe
```

The latest validation pass parsed 34 active R scripts successfully.

## Reproducibility Notes

- Scripts currently assume they are run from the project root.
- Several scripts write intermediate files that are consumed by later scripts.
- Package versions are not yet frozen; adding `renv` or a documented package
  preflight is a good future step.
- The manuscript draft is used locally for scientific context but is not tracked
  in this repository at this stage.
