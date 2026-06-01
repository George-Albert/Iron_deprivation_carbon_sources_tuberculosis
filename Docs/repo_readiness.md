# Repository Readiness Notes

## Current Git State

The local `.git` directory was incomplete and had to be reinitialized so Git
could recognize the repository.

The public GitHub repository
`George-Albert/Iron_deprivation_carbon_sources_tuberculosis` has since been
replaced with this cleaned local history and is now the active remote for the
project.

## Publishing Policy

Track:

- analysis scripts in `Codes/`;
- project README and lightweight documentation;
- validation/configuration files.

Review before tracking:

- `Inputs/`;
- `Outputs/`;
- root-level generated PDFs/XLSX tables;
- manuscript drafts and internal working documents.

Do not track by default:

- local IDE/session files;
- `Codes/trash_codes/`.

## Validation

Run from the project root:

```r
source("scripts/validate_r_parse.R")
```

At the time of this readiness pass, all active root-level scripts in `Codes/`
parse successfully.
