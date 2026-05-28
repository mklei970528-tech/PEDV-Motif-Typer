# PEDV Motif Typer

PEDV Motif Typer is a Shiny web application for identifying PEDV S-protein polymorphic motif types from nucleotide or amino-acid sequences. The app supports sequence typing, strain lookup, curated reference information, phylogenetic context, and spatiotemporal distribution panels for the current PEDV S-gene motif typing system.

## Main Features

- Unified input for FASTA sequence typing and accession or strain lookup.
- Supports nucleotide S genes, amino-acid S proteins, longer nucleotide fragments containing S, lowercase input, aligned sequences with gaps, and common reading-frame offsets.
- Reports five polymorphic-locus states:
  - N57/N62
  - 135/136 motif
  - N1192/N1194
  - G1157
  - N718/N722
- Reports four typing systems:
  - S1-type: N57/N62 + 135/136 motif
  - Genotype: N57/N62 + N1192/N1194
  - Geo-type: N57/N62 + 135/136 motif + N1192/N1194 + G1157
  - Haplotype: N57/N62 + 135/136 motif + N1192/N1194 + G1157 + N718/N722
- Provides curated group descriptions, reference strains with NCBI links, static tree panels, and dynamic spatiotemporal trend plots.

## Repository Contents

- `app.R`: Main Shiny application.
- `sequence_database.csv`: Curated searchable sequence database.
- `typing_description_editable.csv`: Editable description table used by the information panel.
- `four_tree_legend_branch_colors.csv`: Color mapping for tree and trend labels.
- `www/`: Static tree and spatiotemporal image assets.

## Running Locally

Open R in the repository folder and run:

```r
shiny::runApp(".", host = "127.0.0.1", port = 3838)
```

Then open:

```text
http://127.0.0.1:3838/
```

## Required R Packages

The app uses:

- `shiny`
- `bslib`
- `shinycssloaders`
- `Biostrings`
- `dplyr`
- `tidyr`
- `purrr`
- `DT`
- `readxl`

Install missing packages before deployment. `Biostrings` is distributed through Bioconductor.

## Sequence Handling

The backend cleans input sequences, detects nucleotide versus amino-acid input, translates nucleotide input in three reading frames, localizes motif decision windows, and reports only confidently detected loci. Incomplete sequences are allowed; missing downstream regions are reported as undetermined rather than forced into false-negative calls.

Supported input examples include:

- Complete S-gene nucleotide sequences.
- Longer genome fragments containing S.
- Amino-acid S-protein sequences.
- FASTA sequences with line breaks, lowercase letters, or alignment gaps.
- Sequences with limited ambiguity, where key motif windows remain interpretable.

Highly ambiguous sequences or fragments that do not cover the decision loci may be reported as `Undetermined`.

## Validation Summary

The current backend was validated against a curated 3042-sequence S-gene dataset. Final validation showed agreement with curated expected labels for the four typing systems and the five polymorphic-locus states.

## Deployment Notes

Deploy the complete repository folder, not only `app.R`. The application expects the database files, description files, color table, and `www/` assets to be present beside `app.R`.

Before deploying, verify:

```r
invisible(parse("app.R"))
shiny::runApp(".", host = "127.0.0.1", port = 3838)
```

## Updating the App

- Update `sequence_database.csv` to refresh searchable strain metadata and typing results.
- Update `typing_description_editable.csv` to revise group descriptions and reference strains.
- Replace files in `www/` to update static trees or distribution figures.
- Update `four_tree_legend_branch_colors.csv` to keep tree and spatiotemporal colors consistent.

## Citation

If you use this tool in research, please cite or acknowledge this repository:

- Repository: <https://github.com/mklei970528-tech/PEDV-Motif-Typer>
- Authors: Mingkai Lei and Joao Paulo Herrera da Silva, V-EPI lab

Related peer-reviewed article:

Lei, M., Li, H., Chen, X., Li, X., Yu, X., Ruan, S., Wu, H., Ghonaim, A.H., Yan, Z., Li, W. and He, Q. (2025). A novel genotyping system based on site polymorphism on spike gene reveals the evolutionary pathway of porcine epidemic diarrhea virus. *iMetaOmics*, 2, e70013. <https://doi.org/10.1002/imo2.70013>
