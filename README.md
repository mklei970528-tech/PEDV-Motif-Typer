PEDV Motif Typing Platform (Shiny)
=================================

An online Shiny application for Porcine Epidemic Diarrhea Virus (PEDV) classification
based on motif haplotypes.

Live demo (shinyapps.io):  https://raysapp.shinyapps.io/pedvmotiffinder_v2/


FEATURES
--------
- Accepts FASTA sequences as nucleotide (NT) or amino acid (AA) (auto-detected).
- Robust preprocessing to remove illegal characters (e.g., '?', unicode, punctuation).
- NT mode: translates in three frames and selects the best frame by motif alignment score.
- AA mode: aligns and extracts motifs directly.
- Motif extraction via local alignment (DECIPHER AlignSeqs).
- Detects N-linked glycosylation motifs (N-X-[S/T], X≠P) for M1–M3.
- Assigns Motif Type by matching the combined glyco/motif pattern to a curated lookup table.
- Interactive results table and an Information panel linked to references from Information.xlsx.
- Fully compatible with shinyapps.io (≤20MB upload limit set in code).


INPUT & OUTPUT
--------------
Input options:
1) Paste sequences in FASTA format, or
2) Upload a FASTA file (≤20MB)

FASTA examples:

Nucleotide
>seq1
ATGCTG...

Amino acid
>seq2
MKTIIALSYIFCLV...

Output table includes:
- Status (Success / Failed)
- Motif_Type (e.g., G2c L10)
- Typing_Combination (e.g., 12_NXT | 5_NXT | 12_NXT | S-type)
- Motif sequences (M1_AA–M4_AA)
Download: PEDV_motif_typing_results.csv


HOW IT WORKS (HIGH LEVEL)
-------------------------
1) Preprocess input sequence (strip illegal characters, remove gaps).
2) Detect sequence type:
   - NT-like if ≥90% A/T/C/G and low non-standard characters
   - otherwise AA-like
3) NT sequences are translated into three reading frames.
4) For each frame (or AA input directly), the app:
   - aligns the full protein against motif templates (M1–M4)
   - extracts the aligned motif region and computes identity (%)
5) If any motif fails to align above threshold -> Status = Failed
6) Otherwise:
   - calls N-glyco motifs for M1–M3
   - calls M4 residue type at position 4 (G/S/R/N -> *-type)
   - concatenates into a Typing_Combination
   - assigns Motif_Type using the mapping table


REQUIREMENTS
------------
R (recommended: R >= 4.2)

R packages:
- shiny, bslib, shinycssloaders, DT
- Biostrings, DECIPHER
- dplyr, tidyr, purrr
- readxl


RUN LOCALLY
-----------
1) Clone the repository:
   git clone https: https://github.com/mklei970528-tech/Motif-Typer.git

2) Install dependencies in R:
   install.packages(c("shiny","bslib","shinycssloaders","DT","dplyr","tidyr","purrr","readxl"))
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install(c("Biostrings","DECIPHER"))

3) Ensure required file exists at project root:
   - Information.xlsx
     Sheet 1: information table (must include a Motif_Type column)
     Sheet 2: reference mapping table (columns: RefID, Link)

4) Start the app:
   shiny::runApp("app.R")
   (or shiny::runApp("."))


RECOMMENDED PROJECT STRUCTURE
-----------------------------
.
├── app.R
├── Information.xlsx
├── README.txt
└── LICENSE


COMMON ISSUES
-------------
- “Information.xlsx not found”
  Make sure Information.xlsx is in the same directory as app.R (project root).

- “Invalid FASTA format”
  Input must contain at least one header line starting with '>'.

- Sequences fail typing
  If any motif (M1–M4) cannot be aligned above the identity cutoff, the sequence is labeled Failed.
  Verify sequence coverage/length for motif regions.


CITATION
--------
If you use this tool in research, please cite/acknowledge this repository:
- Repository: https://github.com/mklei970528-tech/Motif-Typer.git
- Authors: <Mingkai Lei / V-EPI lab>
(Peer reviewed article citation: Lei, M., Li, H., Chen, X., Li, X., Yu, X., Ruan, S., Wu, H., Ghonaim, A.H., Yan, Z., Li, W. and He, Q. (2025), A novel genotyping system based on site polymorphism on spike gene reveals the evolutionary pathway of porcine epidemic diarrhea virus. iMetaOmics, 2: e70013. https://doi.org/10.1002/imo2.70013)


LICENSE
-------
MIT License

Copyright (c) 2026 mklei970528-tech

CONTACT
-------
For questions or bug reports, please open an Issue on GitHub.
