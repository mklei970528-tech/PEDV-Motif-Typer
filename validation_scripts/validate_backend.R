#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readxl)
  library(DT)
  library(shiny)
  library(bslib)
  library(shinycssloaders)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[[hit + 1]]
}

app_dir <- normalizePath(get_arg("--app-dir", "."), winslash = "/", mustWork = TRUE)
input_fasta <- normalizePath(get_arg("--input-fasta"), winslash = "/", mustWork = TRUE)
expected_csv <- get_arg("--expected-labels", NA_character_)
out_dir <- normalizePath(get_arg("--out-dir", file.path(app_dir, "validation_output")),
                         winslash = "/", mustWork = FALSE)
max_n <- as.integer(get_arg("--max-n", "0"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load_app_functions <- function(app_file) {
  lines <- readLines(app_file, warn = FALSE, encoding = "UTF-8")
  marker <- grep("^# RUN APP", lines)
  if (length(marker) == 0) stop("Cannot find '# RUN APP' marker in app.R")
  code <- paste(lines[seq_len(marker[1] - 1)], collapse = "\n")
  eval(parse(text = code), envir = .GlobalEnv)
}

load_app_functions(file.path(app_dir, "app.R"))

dss <- parse_fasta_file(input_fasta)
if (is.null(dss)) stop("Invalid FASTA input")
raw_seq <- attr(dss, "raw_seq")
raw_id <- attr(dss, "raw_id")
if (!is.null(max_n) && max_n > 0 && length(raw_seq) > max_n) {
  keep <- seq_len(max_n)
  attr(dss, "raw_seq") <- raw_seq[keep]
  attr(dss, "raw_id") <- raw_id[keep]
  raw_seq <- raw_seq[keep]
  raw_id <- raw_id[keep]
}

web_output <- purrr::map2_dfr(raw_seq, raw_id, function(seq, id) {
  analyze_one_sequence(seq, id, motif_AA)
})

display_cols <- c(
  "sequence_ID",
  "Five_locus_haplotype",
  "S_4_locus_geo_type",
  "S_two_locus_type",
  "S1_two_locus_type",
  "N57_N62_state",
  "Site_135_136_motif",
  "N1192_N1194_state",
  "G1157_state",
  "N718_N722_state"
)
display_output <- web_output %>% dplyr::select(dplyr::any_of(display_cols))
write.csv(display_output, file.path(out_dir, "validation_web_output.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

locus_cols <- c(
  "N57_N62_state",
  "Site_135_136_motif",
  "N1192_N1194_state",
  "G1157_state",
  "N718_N722_state"
)

locus_summary <- purrr::map_dfr(locus_cols, function(col) {
  x <- as.character(web_output[[col]])
  total <- length(x)
  not_detected <- sum(x == NOT_DETECTED, na.rm = TRUE)
  other <- sum(x %in% c("Other", "Unassigned"), na.rm = TRUE)
  detected <- total - not_detected
  data.frame(
    Locus = col,
    Total = total,
    Detected_or_classified = detected,
    Not_detected = not_detected,
    Other_or_unassigned = other,
    Detection_rate = ifelse(total > 0, detected / total, NA_real_)
  )
})
write.csv(locus_summary, file.path(out_dir, "validation_locus_detection_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

agreement_summary <- data.frame()
failed_cases <- data.frame()

if (!is.na(expected_csv) && file.exists(expected_csv)) {
  expected <- read.csv(expected_csv, check.names = FALSE, stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8-BOM")
  joined <- display_output %>%
    left_join(expected, by = "sequence_ID")
  write.csv(joined, file.path(out_dir, "validation_joined_expected_observed.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  comparisons <- list(
    S1_type = c("S1_two_locus_type", "expected_S1_type"),
    Genotype = c("S_two_locus_type", "expected_Genotype"),
    Geo_type = c("S_4_locus_geo_type", "expected_Geo_type"),
    Haplotype = c("Five_locus_haplotype", "expected_Haplotype"),
    N57_N62_state = c("N57_N62_state", "expected_N57_N62_state"),
    Site_135_136_motif = c("Site_135_136_motif", "expected_Site_135_136_motif"),
    N1192_N1194_state = c("N1192_N1194_state", "expected_N1192_N1194_state"),
    G1157_state = c("G1157_state", "expected_G1157_state"),
    N718_N722_state = c("N718_N722_state", "expected_N718_N722_state")
  )

  agreement_summary <- purrr::imap_dfr(comparisons, function(cols, label) {
    obs <- cols[[1]]
    exp <- cols[[2]]
    if (!(exp %in% names(joined))) {
      return(data.frame(Typing_system = label, Compared = 0, Matched = 0,
                        Mismatched = 0, Agreement_rate = NA_real_))
    }
    ok <- nzchar(as.character(joined[[exp]]))
    compared <- sum(ok, na.rm = TRUE)
    matched <- sum(as.character(joined[[obs]]) == as.character(joined[[exp]]) & ok, na.rm = TRUE)
    data.frame(
      Typing_system = label,
      Compared = compared,
      Matched = matched,
      Mismatched = compared - matched,
      Agreement_rate = ifelse(compared > 0, matched / compared, NA_real_)
    )
  })

  failed_cases <- joined %>%
    mutate(
      S1_match = if ("expected_S1_type" %in% names(.)) S1_two_locus_type == expected_S1_type else NA,
      Genotype_match = if ("expected_Genotype" %in% names(.)) S_two_locus_type == expected_Genotype else NA,
      Geo_type_match = if ("expected_Geo_type" %in% names(.)) S_4_locus_geo_type == expected_Geo_type else NA,
      Haplotype_match = if ("expected_Haplotype" %in% names(.)) Five_locus_haplotype == expected_Haplotype else NA,
      N57_N62_match = if ("expected_N57_N62_state" %in% names(.)) N57_N62_state == expected_N57_N62_state else NA,
      Site_135_136_match = if ("expected_Site_135_136_motif" %in% names(.)) Site_135_136_motif == expected_Site_135_136_motif else NA,
      N1192_N1194_match = if ("expected_N1192_N1194_state" %in% names(.)) N1192_N1194_state == expected_N1192_N1194_state else NA,
      G1157_match = if ("expected_G1157_state" %in% names(.)) G1157_state == expected_G1157_state else NA,
      N718_N722_match = if ("expected_N718_N722_state" %in% names(.)) N718_N722_state == expected_N718_N722_state else NA
    ) %>%
    filter(if_any(ends_with("_match"), ~ !is.na(.) & !.))

  write.csv(agreement_summary, file.path(out_dir, "validation_typing_agreement_summary.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(failed_cases, file.path(out_dir, "validation_failed_cases.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}

summary_text <- c(
  "# Backend Validation Summary",
  "",
  paste0("- Input FASTA: `", input_fasta, "`"),
  paste0("- Sequences processed: ", nrow(web_output)),
  paste0("- Output directory: `", out_dir, "`"),
  "",
  "## Locus Detection",
  "",
  paste(capture.output(print(locus_summary, row.names = FALSE)), collapse = "\n")
)

if (nrow(agreement_summary) > 0) {
  summary_text <- c(
    summary_text,
    "",
    "## Expected-Label Agreement",
    "",
    paste(capture.output(print(agreement_summary, row.names = FALSE)), collapse = "\n"),
    "",
    paste0("- Mismatched cases: ", nrow(failed_cases))
  )
}

writeLines(summary_text, file.path(out_dir, "validation_result_and_analysis.md"), useBytes = TRUE)

cat("Validation complete\n")
cat("Output directory:", out_dir, "\n")
