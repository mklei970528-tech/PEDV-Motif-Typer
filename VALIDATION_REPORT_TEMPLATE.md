# PEDV S Typing Website Validation Report Template

## 1. Validation Dataset

- Dataset name:
- Source:
- Number of sequences:
- Sequence formats included:
  - Nucleotide S gene:
  - Amino-acid S protein:
  - Aligned sequence:
  - Unaligned sequence:
  - Genome-length fragment:
  - Partial sequence:
- Independent from the searchable website database: yes/no

## 2. Validation Objectives

1. Confirm that the website can read and process representative PEDV S sequences.
2. Confirm that five polymorphic loci are detected or conservatively marked as not detected.
3. Confirm that four typing systems are assigned consistently with manually curated or database-derived labels.
4. Confirm that partial sequences do not produce false negative calls for uncovered loci.
5. Confirm that accession/strain database search returns the expected background information.

## 3. Required Output Files

- `validation_input.fasta`
- `validation_expected_labels.csv`
- `validation_web_output.csv`
- `validation_locus_detection_summary.csv`
- `validation_typing_agreement_summary.csv`
- `validation_failed_cases.csv`
- `validation_result_and_analysis.md`

## 4. Locus Detection Summary

| Locus | Detected | Not detected | Other/ambiguous | Detection rate |
|---|---:|---:|---:|---:|
| N57/N62 |  |  |  |  |
| 135/136 motif |  |  |  |  |
| N1192/N1194 |  |  |  |  |
| G1157 |  |  |  |  |
| N718/N722 |  |  |  |  |

## 5. Typing Agreement Summary

| Typing system | Compared sequences | Matched | Mismatched | Agreement rate |
|---|---:|---:|---:|---:|
| S1-type |  |  |  |  |
| Genotype |  |  |  |  |
| Geo-type |  |  |  |  |
| Haplotype |  |  |  |  |

## 6. Partial Sequence Handling

Summarize whether short or incomplete sequences were typed only with confidently detected loci.

## 7. Database Search Checks

| Query type | Example query | Expected match | Observed match | Status |
|---|---|---|---|---|
| Accession |  |  |  |  |
| Strain name |  |  |  |  |
| Region |  |  |  |  |
| Typing label |  |  |  |  |

## 8. Thesis-Ready Result Draft

Write a concise paragraph describing validation dataset composition, locus detection performance, typing consistency, and limitations.

