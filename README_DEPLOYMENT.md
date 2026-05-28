# PEDV S Polymorphism Typing Platform

This folder contains the deployable Shiny application for PEDV S-gene polymorphism typing.

## Core Files

- `app.R`: Main Shiny application.
- `Information_current_typing.xlsx`: Curated typing information, reference strains, descriptions, and reference links.
- `sequence_database.csv`: Curated sequence database used by the database-search module.
- `www/`: Static image assets used by the information and database panels.

## Main Website Modules

1. **Unified sequence typing / strain search**
   - Uses one input box for both de novo sequence typing and database lookup.
   - FASTA-like or long nucleotide/amino-acid input is analyzed as sequence data.
   - Short accession, strain, country/region, year, locus-state, or typing-label queries are treated as database searches.

2. **Guide & typing logic**
   - Shows the four typing-system layers as a compact logic chart.
   - Explains which polymorphic loci define S1-type, genotype, geo-type, and five-locus haplotype.

3. **Results & group information**
   - Uses one shared result panel for both sequence input and strain search.
   - Displays the five polymorphic-locus states plus four typing-system assignments.
   - Clicking or uniquely matching a result updates the reference strains, group description, tree panel, and spatiotemporal trend panel.

4. **Phylogenetic and spatiotemporal context**
   - Displays one context tree selected from the highest confidently available typing level.
   - Shows dynamic global or region/country-specific trends for the selected typing result.

## Sequence Analysis Behavior

The sequence module:

   - Accepts pasted or uploaded FASTA sequences.
   - Supports nucleotide S genes, amino-acid S proteins, and longer nucleotide fragments containing S.
   - Cleans input sequences, translates nucleotide input in three reading frames, localizes decision frames, and reports locus states plus four typing-system assignments.
   - For incomplete sequences, reports only confidently detected loci and avoids converting missing regions into false-negative states.

## Typing Systems

- **S1-type**: N57/N62 + 135/136 motif.
- **Genotype**: N57/N62 + N1192/N1194.
- **Geo-type**: N57/N62 + 135/136 motif + N1192/N1194 + G1157.
- **Haplotype**: N57/N62 + 135/136 motif + N1192/N1194 + G1157 + N718/N722.

## Deployment Notes

Upload the entire `D:/test/2.3_deploy` folder, not only `app.R`. The app expects `Information_current_typing.xlsx`, `sequence_database.csv`, and the `www/` image assets to be present in the same deployment directory.

Before deployment, run:

```r
invisible(parse("D:/test/2.3_deploy/app.R"))
shiny::runApp("D:/test/2.3_deploy", host = "127.0.0.1", port = 3838)
```

Then confirm that the following modules appear in the browser:

- Guide & typing logic
- Unified sequence typing / strain search
- Results and group information
- Phylogenetic context
- Spatiotemporal dynamics

## Validation

The current deployable backend was validated against the curated 3042-sequence S-gene dataset:

```powershell
Rscript D:/test/2.3_deploy_ready/validation_scripts/validate_backend.R `
  --app-dir D:/test/2.3_deploy_ready `
  --input-fasta C:/Users/LeiMingkai/OneDrive/PEDV_S_selection_portable_JN599150_20260526/inputs/S_merged_modified.five_site_main_3042.protein.fasta `
  --expected-labels D:/test/2.3_deploy_ready/validation_inputs/expected_labels_3042_from_sequence_database.csv `
  --out-dir D:/test/2.3_deploy_ready/validation_output_3042_20260528_final
```

Final validation result: all 3042 sequences matched the curated expected labels for S1-type, genotype, geo-type, haplotype, and all five polymorphic-locus states. The final report is `validation_output_3042_20260528_final/validation_result_and_analysis.md`.

## Updating Data

- To update typing descriptions, edit `Information_current_typing.xlsx`.
- To update the searchable database, replace `sequence_database.csv` with a table using the same column names.
- To update static figures, replace the matching PNG files in `www/`.
