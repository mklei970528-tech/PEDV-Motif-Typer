# Deployment Manifest

This folder is the clean deployment package for the PEDV S Polymorphism Typing Platform.

## Included Runtime Files

| Path | Purpose |
|---|---|
| `app.R` | Main Shiny application. |
| `Information_current_typing.xlsx` | Typing descriptions, reference strains, and reference links. |
| `sequence_database.csv` | Searchable curated PEDV S-gene sequence database. |
| `www/` | Static spatiotemporal and phylogenetic image assets. |
| `rsconnect/` | Existing shinyapps.io deployment metadata, if reused. |

## Included Documentation and Templates

| Path | Purpose |
|---|---|
| `README_DEPLOYMENT.md` | Deployment and update instructions. |
| `OUTPUT_FIELD_DICTIONARY.md` | Explanation of output and database fields. |
| `VALIDATION_REPORT_TEMPLATE.md` | Template for independent validation reporting. |
| `WEBSITE_METHODS_DRAFT_CN.md` | Chinese thesis-methods draft for the website section. |
| `validation_expected_labels_template.csv` | Template for expected labels in future validation datasets. |

## Files intentionally excluded

The original working directory may contain diagnostics, failed-case reports, and historical backups. These are not required for deployment and were excluded from this clean package.

## Pre-deployment checklist

- [x] `app.R` parses without error.
- [x] Local Shiny launch works from this directory.
- [x] The following pages are visible: Guide & typing logic, Sequence typing, Database search, Information, Validation & thesis use.
- [x] Static images under `www/` render in the Information and Database search pages.
- [x] `Information_current_typing.xlsx` and `sequence_database.csv` are in the same directory as `app.R`.
- [ ] Upload the full folder rather than only `app.R`.

Last local package check: 2026-05-28, using `http://127.0.0.1:3842/`.
