# Repository Readiness Notes

## Current Git State

The local `.git` directory was incomplete and had to be reinitialized so Git
could recognize the repository. This repository currently has no historical
commits and no remote configured.

The old GitHub repository `Iron_deprivation_carbon_sources_tuberculosis` should
not be treated as canonical for this local project unless explicitly reconnected
later.

## Publishing Policy

Track:

- analysis scripts in `Codes/`;
- project README and lightweight documentation;
- validation/configuration files.

Do not track by default:

- `Inputs/`;
- `Outputs/`;
- root-level generated PDFs/XLSX tables;
- local IDE/session files;
- `Codes/trash_codes/`.

## Validation

Run from the project root:

```r
source("scripts/validate_r_parse.R")
```

At the time of this readiness pass, all active root-level scripts in `Codes/`
parse successfully.
