# Output Field Dictionary

## Sequence Typing Output

| Field | Meaning |
|---|---|
| `sequence_ID` | Input sequence identifier parsed from the FASTA header. |
| `Five_locus_haplotype` | Five-locus haplotype assignment based on N57/N62, 135/136, N1192/N1194, G1157, and N718/N722. |
| `S_4_locus_geo_type` | Four-locus geography-associated type based on N57/N62, 135/136, N1192/N1194, and G1157. |
| `S_two_locus_type` | S-gene two-locus genotype based on N57/N62 and N1192/N1194. |
| `S1_two_locus_type` | S1-region type based on N57/N62 and the 135/136 motif. |
| `N57_N62_state` | N-terminal glycosylation state around N57/N62. |
| `Site_135_136_motif` | Local amino-acid motif state around positions 135/136. |
| `N1192_N1194_state` | S2/CD-region glycosylation state around N1192/N1194. |
| `G1157_state` | Amino-acid state at the reference-anchored G1157 position. |
| `N718_N722_state` | Paired glycosylation state around N718/N722. |

## Common State Labels

| Label | Meaning |
|---|---|
| `NXS` | N-linked glycosylation motif with serine as the third residue. |
| `NXT` | N-linked glycosylation motif with threonine as the third residue. |
| `No glycosylation` | The local decision frame was detected, but the required glycosylation motif was absent. |
| `Not detected` | The local decision frame was not confidently detected or not covered by the input sequence. |
| `Not available` | The typing system cannot be assigned because one or more required loci are missing. |
| `Other` | A detected motif/state falls outside the predefined typing categories. |

## Search Database Fields

| Field | Meaning |
|---|---|
| `Accession` | GenBank or curated accession identifier. |
| `Sequence ID` | Sequence name used in the curated S-gene database. |
| `Strain/Isolate` | Strain or isolate name when available. |
| `Country`, `Region`, `Location` | Standardized spatial metadata. |
| `Collection date`, `Year` | Standardized temporal metadata. |
| `Haplotype`, `Geo-type`, `Genotype`, `S1-type` | Four typing-system assignments. |
| `Haplotype pattern`, `Geo-type pattern`, `Genotype pattern`, `S1-type pattern` | Underlying polymorphic state combinations. |
| `Notes` | Processing or typing notes. |

