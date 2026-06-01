# Publication Data Review

This repository currently publishes the analysis code, documentation, project metadata, and validation helpers. Data inputs, generated outputs, figures, and large result tables are currently kept out of Git until a publication decision is made. The review focus is reproducibility, file size, provenance, and manuscript timing rather than privacy.

## Current Local-Only Candidates

| Path or group | Approximate size | Notes for review |
| --- | ---: | --- |
| `Inputs/001_Raw_data/` | 24.9 MB | Raw mapping/count inputs and reference files. Review licensing, provenance, and whether raw tables should be redistributed. |
| `Inputs/002_Processed_data/` | 71.0 MB | Processed differential expression, GO, RDS, and intermediate analysis inputs. Likely useful for reproducibility, but derived from the working analysis. |
| `Outputs/Data/` | 117.2 MB | Large exported result tables, including many `.txt`/`.xlsx` files and `data_MJ.zip`. Useful for reproducibility but high volume. |
| `Outputs/001_Figures_paper/` | 7.9 MB | Generated paper/supplement figures. Consider publishing only final figures or source data needed to regenerate them. |
| `Outputs/Volcano_plots_shared/` | 6.0 MB | Shared SVG volcano panels. Good candidate if these are final reusable figures. |
| `Outputs/gene_names_by_class/` | 12.8 MB | FASTA/GFF/reference-style outputs. Review source/provenance before publishing. |
| Root `Table_1_Section_2.xlsx` | 8.6 MB | Manuscript/reporting table. Review whether it is final and safe for public release. |
| Root `GO_enrichment_*.xlsx` | small | Generated enrichment exports. Potentially useful to publish if they correspond to final reported analyses. |
| Root `Figure_*.pdf` | small | Generated figure export. Review whether it is final. |
| `Docs/pnas.1718003115.sd01.xlsx` | unknown | Third-party supplemental data candidate; review redistribution rights before publishing. |

## Current Recommendation

The first public baseline is code-only. For a reproducible public release, choose one of these follow-up strategies:

1. Publish a minimal reproducibility bundle with selected processed input tables and final result tables.
2. Publish final figures/tables only, while documenting how private/local inputs are expected to be placed.
3. Publish all non-sensitive processed outputs if licensing and file size are acceptable.
4. Keep large data outside Git and link to a release, Zenodo/OSF record, or manuscript supplement.

## Decision Needed

Before changing `.gitignore` or adding data, decide which local-only groups should be tracked in Git, attached as GitHub release assets, or deposited externally.
