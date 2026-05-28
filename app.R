###############################
## PEDV Motif Typing 鈥?FIXED SHINYAPPS VERSION
## - Illegal char filtering (fix '?')
## - NT & AA independent processing
## - Three-frame translation
## - AA local alignment
## - Fully shinyapps.io compatible
###############################


options(shiny.maxRequestSize = 20 * 1024^2)

library(shiny)
library(bslib)
library(shinycssloaders)
library(Biostrings)
library(dplyr)
library(tidyr)
library(purrr)
library(DT)
library(readxl)


###############################
# Valid AA symbols
###############################
AA_ALL <- c(
  "A","R","N","D","C","Q","E","G","H","I",
  "L","K","M","F","P","S","T","W","Y","V",
  "B","Z","J","X","U","O"
)


###############################
# FASTA parser (raw, no filtering)
###############################
parse_fasta_raw <- function(text_lines){
  
  # remove blank lines
  txt <- trimws(text_lines[nzchar(trimws(text_lines))])
  headers <- which(startsWith(txt, ">"))
  if(length(headers)==0) {
    seq <- toupper(paste0(txt, collapse=""))
    seq <- gsub("\\s+", "", seq)
    if(!nzchar(seq)) return(NULL)
    
    dss <- DNAStringSet("A")
    names(dss) <- "input_sequence"
    attr(dss,"raw_seq") <- seq
    attr(dss,"raw_id")  <- "input_sequence"
    return(dss)
  }
  
  ids <- sub("^>", "", txt[headers])
  seqs <- character(length(headers))
  
  for(i in seq_along(headers)){
    s <- headers[i] + 1
    e <- if(i==length(headers)) length(txt) else headers[i+1]-1
    seqs[i] <- toupper(paste0(txt[s:e], collapse=""))
  }
  
  dss <- DNAStringSet(rep("A", length(seqs)))
  names(dss) <- ids
  
  attr(dss,"raw_seq") <- seqs
  attr(dss,"raw_id")  <- ids
  return(dss)
}

parse_fasta_file <- function(path){
  txt <- readLines(path, warn=FALSE)
  parse_fasta_raw(txt)
}


###############################
# Enhanced preprocess: FIX '?' + illegal chars
###############################
preprocess_sequence <- function(seq_raw){
  
  # Remove whitespace, unicode, '?', digits, punctuation, and alignment gaps.
  s <- gsub("[^A-Z\\-]", "", toupper(seq_raw))
  s <- gsub("-", "", s)
  s <- gsub("U", "T", s, fixed = TRUE)
  
  if(nchar(s)==0)
    return(list(clean_seq="", seq_type="NT"))
  
  chars <- strsplit(s,"")[[1]]
  nt_all <- c("A","T","C","G","R","Y","S","W","K","M","B","D","H","V","N")
  nt <- sum(chars %in% nt_all)
  aa   <- sum(chars %in% AA_ALL)
  non  <- length(chars) - nt - aa
  aa_only <- sum(chars %in% c("E","F","I","L","P","Q","J","O","Z","X"))
  
  prop_nt <- nt / length(chars)
  prop_non  <- non  / length(chars)
  
  # NT-like. Ambiguous IUPAC bases are kept as N and translated to X later.
  if(aa_only == 0 && prop_nt >= 0.85 && prop_non <= 0.05){
    clean <- gsub("[^ACGT]", "N", s)
    return(list(clean_seq=clean, seq_type="NT"))
  }
  
  # AA-like
  clean <- gsub(paste0("[^", paste(AA_ALL, collapse=""), "]"), "A", s)
  return(list(clean_seq=clean, seq_type="AA"))
}


###############################
# Clean AA: FIX '?' and all illegal AA residues
###############################
clean_aa_for_alignment <- function(aa_string){
  aa <- strsplit(aa_string, "")[[1]]
  aa[!(aa %in% AA_ALL)] <- "A"
  paste0(aa, collapse="")
}


###############################
# Motif AA list
###############################
motif_AA <- list(
  M1A = "GYLPIGENQGVNSTWYC",
  M1B = "GYLPSMNSSSWYC",
  M135 = "KTLGPTVNDVTTGRN",
  M2  = "FTHELQNHTATEY",
  M2B = "FTHELQDTATEY",
  M3  = "GVISSLSSSTFNSTRELP",
  M4  = "TVLVPGDFVNVIAIDG"
)


###############################
# Enhanced NT 3-frame translation: fixes unicode
###############################
translate_three_frames <- function(seq_raw){
  
  seq_raw <- gsub("[^ACGTN]", "N", toupper(seq_raw))
  dna <- DNAString(seq_raw)
  len <- length(dna)
  
  rf1 <- if(len>=3) subseq(dna,1,len-len%%3) else DNAString("")
  rf2 <- if(len>=4) subseq(dna,2,len-(len-1)%%3) else DNAString("")
  rf3 <- if(len>=5) subseq(dna,3,len-(len-2)%%3) else DNAString("")
  
  list(
    AA1 = translate(rf1, if.fuzzy.codon="X"),
    AA2 = translate(rf2, if.fuzzy.codon="X"),
    AA3 = translate(rf3, if.fuzzy.codon="X")
  )
}


###############################
# Fast AA motif extraction by local decision windows
###############################
empty_motif <- function(identity = 0) {
  list(aa="Not aligned", identity=identity)
}

best_window_motif <- function(fullAA, motifAA, search_start = NULL, search_end = NULL, cut = 50) {
  fullA <- clean_aa_for_alignment(as.character(fullAA))
  motA  <- clean_aa_for_alignment(motifAA)
  n <- nchar(fullA)
  m <- nchar(motA)
  
  if(n == 0 || m == 0 || n < m)
    return(empty_motif())
  if (!is.null(search_start) && n < search_start + m - 1)
    return(empty_motif())
  
  exact <- regexpr(motA, fullA, fixed=TRUE)[[1]]
  if(exact > 0)
    return(list(aa=substr(fullA, exact, exact + m - 1), identity=100, start=exact))
  
  s <- if (is.null(search_start)) 1 else max(1, search_start)
  e <- if (is.null(search_end)) n else min(n, search_end)
  s <- min(s, max(1, n - m + 1))
  e <- max(e, s + m - 1)
  if (e - s + 1 < m) return(empty_motif())
  
  best_identity <- -1
  best_start <- NA_integer_
  mot_chars <- strsplit(motA, "")[[1]]
  
  last_start <- min(e - m + 1, n - m + 1)
  for (pos in s:last_start) {
    win <- substr(fullA, pos, pos + m - 1)
    win_chars <- strsplit(win, "")[[1]]
    identity <- sum(win_chars == mot_chars) / m * 100
    if (identity > best_identity) {
      best_identity <- identity
      best_start <- pos
    }
  }
  
  if (!is.na(best_start) && best_identity >= cut) {
    return(list(
      aa = substr(fullA, best_start, best_start + m - 1),
      identity = best_identity,
      start = best_start
    ))
  }
  
  empty_motif()
}

best_fixed_width_motif <- function(fullAA, motifAA){
  best_window_motif(fullAA, motifAA, cut = 100)
}

local_align_motif <- function(fullAA, motifAA, cut=30){
  fast <- best_window_motif(fullAA, motifAA, cut = cut)
  if(fast$identity >= cut)
    return(fast)
  
  empty_motif(fast$identity)
}



###############################
# Extract all motifs
###############################
extract_all_motifs <- function(fullAA, motif_AA){
  
  M1A <- best_window_motif(fullAA, motif_AA$M1A, cut = 50)
  M1B <- best_window_motif(fullAA, motif_AA$M1B, cut = 50)
  M1  <- if(M1A$identity >= M1B$identity) M1A else M1B
  
  list(
    M1 = M1,
    M135 = best_window_motif(fullAA, motif_AA$M135, cut = 50),
    M2 = {
      M2A <- best_window_motif(fullAA, motif_AA$M2, cut = 50)
      M2B <- best_window_motif(fullAA, motif_AA$M2B, cut = 50)
      if (M2A$identity >= M2B$identity) M2A else M2B
    },
    M3 = best_window_motif(fullAA, motif_AA$M3, cut = 50),
    M4 = best_window_motif(fullAA, motif_AA$M4, cut = 50),
    M1A_ident = M1A$identity,
    M1B_ident = M1B$identity
  )
}


###############################
# Frame scoring
###############################
select_best_frame <- function(frames, motif_AA){
  
  score_frame <- function(fullAA){
    aa <- as.character(fullAA)
    stop_hits <- gregexpr("\\*", aa, fixed=TRUE)[[1]]
    stop_count <- if (stop_hits[1] < 0) 0 else length(stop_hits)
    id <- c(
      max(best_window_motif(aa, motif_AA$M1A, cut = 0)$identity,
          best_window_motif(aa, motif_AA$M1B, cut = 0)$identity),
      best_window_motif(aa, motif_AA$M135, cut = 0)$identity,
      max(best_window_motif(aa, motif_AA$M2, cut = 0)$identity,
          best_window_motif(aa, motif_AA$M2B, cut = 0)$identity),
      best_window_motif(aa, motif_AA$M3, cut = 0)$identity,
      best_window_motif(aa, motif_AA$M4, cut = 0)$identity
    )
    exact_hits <- sum(vapply(
      c(motif_AA$M1A, motif_AA$M1B, motif_AA$M135, motif_AA$M2,
        motif_AA$M2B, motif_AA$M3, motif_AA$M4, "GPT", "GPS",
        "DFVNVIAI", "DFVDVIAI", "DFIDVIAI"),
      function(m) grepl(m, aa, fixed=TRUE),
      logical(1)
    ))
    sum(id, na.rm = TRUE) * 10 + exact_hits * 100 - stop_count * 0.1
  }
  
  sc <- c(score_frame(frames$AA1), score_frame(frames$AA2), score_frame(frames$AA3))
  best <- which.max(sc)
  best_motifs <- extract_all_motifs(frames[[best]], motif_AA)
  
  list(motifs=best_motifs, best_frame=best, frame_scores=sc)
}


###############################
# Glycan detection
###############################
find_nglyco <- function(aa){
  if(is.na(aa) || aa=="Not aligned" || nchar(aa) < 3) return(tibble())
  a <- strsplit(aa,"")[[1]]
  out <- tibble(pos=integer(), type=character())
  for(i in 1:(length(a)-2)){
    if(a[i]=="N" && a[i+1]!="P" && a[i+2] %in% c("S","T"))
      out <- add_row(out,pos=i,type=paste0("NX",a[i+2]))
  }
  out
}

glyco_pattern <- function(g){
  if (nrow(g) == 0) return("not detected")
  return(paste(paste0(g$pos,"_",g$type), collapse=";"))
}

align_to_reference_frame <- function(ref, query) {
  ref <- clean_aa_for_alignment(gsub("-", "", ref, fixed = TRUE))
  query <- clean_aa_for_alignment(gsub("-", "", query, fixed = TRUE))
  if (!nzchar(ref) || !nzchar(query)) {
    return(list(ref = "", query = ""))
  }
  
  r <- strsplit(ref, "")[[1]]
  q <- strsplit(query, "")[[1]]
  nr <- length(r)
  nq <- length(q)
  score <- matrix(0, nr + 1, nq + 1)
  trace <- matrix("", nr + 1, nq + 1)
  gap <- -1
  mismatch <- -2
  match <- 2
  
  for (i in 1:nr) {
    score[i + 1, 1] <- score[i, 1] + gap
    trace[i + 1, 1] <- "U"
  }
  for (j in 1:nq) {
    score[1, j + 1] <- score[1, j] + gap
    trace[1, j + 1] <- "L"
  }
  
  for (i in 1:nr) {
    for (j in 1:nq) {
      diag_score <- score[i, j] + ifelse(r[i] == q[j], match, mismatch)
      up_score <- score[i, j + 1] + gap
      left_score <- score[i + 1, j] + gap
      scores <- c(diag_score, up_score, left_score)
      best <- which.max(scores)
      score[i + 1, j + 1] <- scores[best]
      trace[i + 1, j + 1] <- c("D", "U", "L")[best]
    }
  }
  
  i <- nr + 1
  j <- nq + 1
  ar <- character()
  aq <- character()
  while (i > 1 || j > 1) {
    step <- trace[i, j]
    if (step == "D") {
      ar <- c(r[i - 1], ar)
      aq <- c(q[j - 1], aq)
      i <- i - 1
      j <- j - 1
    } else if (step == "U") {
      ar <- c(r[i - 1], ar)
      aq <- c("-", aq)
      i <- i - 1
    } else {
      ar <- c("-", ar)
      aq <- c(q[j - 1], aq)
      j <- j - 1
    }
  }
  
  list(ref = paste0(ar, collapse = ""), query = paste0(aq, collapse = ""))
}

query_at_ref_pos <- function(aln, ref_pos) {
  ar <- strsplit(aln$ref, "")[[1]]
  aq <- strsplit(aln$query, "")[[1]]
  ref_count <- 0
  for (i in seq_along(ar)) {
    if (ar[i] != "-") {
      ref_count <- ref_count + 1
      if (ref_count == ref_pos) {
        if (aq[i] == "-") return("")
        return(aq[i])
      }
    }
  }
  ""
}

insertions_before_ref_pos <- function(aln, ref_pos) {
  ar <- strsplit(aln$ref, "")[[1]]
  aq <- strsplit(aln$query, "")[[1]]
  ref_count <- 0
  target_col <- NA_integer_
  for (i in seq_along(ar)) {
    if (ar[i] != "-") {
      ref_count <- ref_count + 1
      if (ref_count == ref_pos) {
        target_col <- i
        break
      }
    }
  }
  if (is.na(target_col) || target_col <= 1) return("")
  ins <- character()
  i <- target_col - 1
  while (i >= 1 && ar[i] == "-") {
    if (aq[i] != "-") ins <- c(aq[i], ins)
    i <- i - 1
  }
  paste0(ins, collapse = "")
}

nglyco_at_ref_pos <- function(aln, ref_pos) {
  triplet <- paste0(
    query_at_ref_pos(aln, ref_pos),
    query_at_ref_pos(aln, ref_pos + 1),
    query_at_ref_pos(aln, ref_pos + 2)
  )
  if (nchar(triplet) < 3) return(NO_GLYCO)
  a <- strsplit(triplet, "")[[1]]
  if (a[1] == "N" && a[2] != "P" && a[3] %in% c("S", "T")) {
    return(paste0("NX", a[3]))
  }
  NO_GLYCO
}

motif_type_M4 <- function(aa){
  if (aa == "Not aligned" || nchar(aa) < 4) 
    return("not detected")
  
  r4 <- substr(aa,4,4)
  if (r4 %in% c("G","S","R","N")) 
    return(paste0(r4, "-type"))
  
  return("Other-type")
}

###############################
# Motif Type Mapping Table
###############################
motif_type_map <- tribble(
  ~pattern,                                       ~Motif_Type,
  "7_NXS | none | 8_NXT;12_NXT | G-type",         "G1a L6",
  "7_NXS | none | 12_NXT | G-type",               "G1a L7",
  "7_NXS | 5_NXT | 12_NXT | S-type",              "G1b L12",
  "7_NXS | 5_NXT | 12_NXT | G-type",              "G1b L8",
  "7_NXS | 5_NXT | 8_NXT;12_NXT | G-type",        "G1b L9",
  "7_NXS | 5_NXT | 8_NXT | G-type",               "G1b L9.1",
  "none | none | 8_NXT;12_NXT | G-type",          "G2a L1",
  "12_NXT | none | 12_NXT | G-type",              "G2b L2",
  "12_NXT | none | 8_NXT;12_NXT | G-type",        "G2b L3",
  "12_NXT | none | 8_NXT;12_NXT | S-type",        "G2b L3.1",
  "12_NXT | 5_NXT | 12_NXT | S-type",             "G2c L10",
  "10_NXT | 5_NXT | 12_NXT | S-type",             "G2c L10.1",
  "12_NXT | 5_NXT | 8_NXT;12_NXT | S-type",       "G2c L11",
  "12_NXT | 5_NXT | 12_NXT | G-type",             "G2c L4",
  "12_NXT | 6_NXT | 12_NXT | G-type",             "G2c L4.1",
  "12_NXT | 5_NXT | 8_NXT;12_NXT | G-type",       "G2c L5",
  "10_NXT | 5_NXT | 8_NXT;12_NXT | G-type",       "G2c L5.1"
)

get_motif_type <- function(tc) {
  row <- motif_type_map %>% filter(pattern == tc)
  if(nrow(row)==0) return("other")
  return(row$Motif_Type[[1]])
}


###############################
# Current S polymorphism typing system
###############################
NO_GLYCO <- "No glycosylation"
NOT_AVAILABLE <- "Undetermined"
NOT_DETECTED <- "Not detected"

has_motif <- function(motif) {
  !is.null(motif) &&
    !is.null(motif$aa) &&
    !is.na(motif$aa) &&
    motif$aa != "Not aligned"
}

extract_local_region_with_insertions <- function(fullAA, motifAA, cut = 50) {
  fullA <- clean_aa_for_alignment(as.character(fullAA))
  motA  <- clean_aa_for_alignment(motifAA)
  
  if (!nzchar(fullA) || !nzchar(motA))
    return(empty_motif())
  
  # The 135/136 motif may contain an insertion immediately before the
  # dipeptide. A regex window preserves GPS/GPT contexts when present.
  if (identical(motA, clean_aa_for_alignment(motif_AA$M135))) {
    hit <- regexpr("[AG]P[A-Z]{1,6}[VA]T[TS]G", fullA, perl = TRUE)
    if (hit[[1]] > 0) {
      aa <- substr(fullA, hit[[1]], hit[[1]] + attr(hit, "match.length") - 1)
      return(list(aa = aa, identity = 100))
    }
  }
  
  fast <- best_fixed_width_motif(fullA, motA)
  if (fast$identity >= cut) return(fast)
  empty_motif(fast$identity)
}

call_n57_n62 <- function(m1_aa) {
  if (is.na(m1_aa) || m1_aa == "Not aligned") return(NOT_DETECTED)
  gp <- glyco_pattern(find_nglyco(m1_aa))
  if (grepl("NXS", gp, fixed = TRUE)) return("NXS")
  if (grepl("NXT", gp, fixed = TRUE)) return("NXT")
  NO_GLYCO
}

call_n57_n62_from_window <- function(fullAA) {
  aa <- clean_aa_for_alignment(as.character(fullAA))
  if (nchar(aa) < 55) return(NOT_DETECTED)
  window <- substr(aa, 45, min(95, nchar(aa)))
  if (!grepl("GYLP|GEN|NSS|NST", window, perl = TRUE)) return(NOT_DETECTED)
  gp <- glyco_pattern(find_nglyco(window))
  if (grepl("NXS", gp, fixed = TRUE)) return("NXS")
  if (grepl("NXT", gp, fixed = TRUE)) return("NXT")
  NO_GLYCO
}

call_135_136 <- function(m135_aa) {
  if (is.na(m135_aa) || m135_aa == "Not aligned") return(NOT_DETECTED)
  aa <- gsub("[^A-Z]", "", toupper(m135_aa))
  if (nchar(aa) < 9) return(NOT_DETECTED)

  # In some SINDEL-derived frames, a residue immediately before 135/136
  # shifts the local string to GPTA-NDD/DND/NND-VTT. Prefer these curated
  # triplet states before falling back to the two-residue readout.
  triplet_patterns <- c(
    "GP[A-Z]{1,3}((?:NND|DND|NDD))[VA]TT",
    "GP[A-Z]{1,3}((?:NND|DND|NDD))[VA]T",
    "GPT[A-Z]?((?:NND|DND|NDD))VTT",
    "GPT[A-Z]?((?:NND|DND|NDD))VT",
    "PT[A-Z]?((?:NND|DND|NDD))VTT",
    "PT[A-Z]?((?:NND|DND|NDD))VT"
  )
  for (pattern in triplet_patterns) {
    hit <- regmatches(aa, regexec(pattern, aa, perl = TRUE))[[1]]
    if (length(hit) >= 2) return(hit[2])
  }
  
  # Align the extracted decision frame to the CV777 reference frame, then
  # read the local ref-position 7-9 block. CV777 is VND, so it is reported
  # as ND; insertion/SINDEL variants map to NND/DND/NDD or SG/GH.
  aln <- align_to_reference_frame(motif_AA$M135, aa)
  q7 <- query_at_ref_pos(aln, 7)
  q8 <- query_at_ref_pos(aln, 8)
  q9 <- query_at_ref_pos(aln, 9)
  motif3 <- paste0(q7, q8, q9)
  motif2_from_7 <- paste0(q7, q8)
  motif2 <- paste0(q8, q9)
  if (motif3 %in% c("NND", "DND", "NDD")) return(motif3)
  if (motif2_from_7 %in% c("SG", "GH")) return(motif2_from_7)
  if (motif2 %in% c("ND", "SG", "GH")) return(motif2)
  
  patterns <- c(
    "GPT[A-Z]?([A-Z]{2,3})VTT",
    "GPT[A-Z]?([A-Z]{2,3})VT",
    "GPS([A-Z]{2,3})VTT",
    "GPS([A-Z]{2,3})VT",
    "PT[A-Z]?([A-Z]{2,3})VTT",
    "PT[A-Z]?([A-Z]{2,3})VT",
    "PS([A-Z]{2,3})VTT",
    "PS([A-Z]{2,3})VT"
  )
  motif <- NOT_DETECTED
  for (pattern in patterns) {
    hit <- regmatches(aa, regexec(pattern, aa, perl = TRUE))[[1]]
    if (length(hit) >= 2) {
      motif <- hit[2]
      break
    }
  }
  if (motif %in% c("ND", "NND", "DND", "SG", "GH", "NDD")) return(motif)
  if (motif == NOT_DETECTED) return(NOT_DETECTED)
  "Other"
}

call_135_136_from_window <- function(fullAA) {
  aa <- clean_aa_for_alignment(as.character(fullAA))
  if (nchar(aa) < 130) return(NOT_DETECTED)
  window <- substr(aa, 105, min(170, nchar(aa)))
  if (!grepl("KTL|CQF|GRN|VTT|VTS|ATT", window, perl = TRUE)) return(NOT_DETECTED)
  known_hit <- regexpr("(NND|DND|NDD|ND|SG|GH)[VA]T[TS]", window, perl = TRUE)
  known <- if (known_hit[[1]] > 0) regmatches(window, known_hit)[[1]] else ""
  if (nzchar(known)) {
    motif <- sub("[VA]T[TS]$", "", known, perl = TRUE)
    if (motif %in% c("ND", "NND", "DND", "SG", "GH", "NDD")) return(motif)
  }
  if (grepl("KTL|CQF|GRN", window, perl = TRUE)) {
    local_known_hit <- regexpr("(NND|DND|NDD|ND|SG|GH)", window, perl = TRUE)
    local_known <- if (local_known_hit[[1]] > 0) regmatches(window, local_known_hit)[[1]] else ""
    if (nzchar(local_known)) return(local_known)
  }
  if (grepl("[A-Z]{2,6}[VA]T[TS]G", window, perl = TRUE)) return("Other")
  NOT_DETECTED
}

call_n1192_n1194 <- function(m2_aa) {
  if (is.na(m2_aa) || m2_aa == "Not aligned") return(NOT_DETECTED)
  gp <- glyco_pattern(find_nglyco(m2_aa))
  if (grepl("NXT", gp, fixed = TRUE)) return("NXT")
  NO_GLYCO
}

call_n718_n722 <- function(m3_aa) {
  if (is.na(m3_aa) || m3_aa == "Not aligned") return(NOT_DETECTED)
  gp <- glyco_pattern(find_nglyco(m3_aa))
  if (grepl("8_NXT", gp, fixed = TRUE) && grepl("12_NXT", gp, fixed = TRUE)) return("NXT")
  NO_GLYCO
}

call_g1157 <- function(m4_aa) {
  if (is.na(m4_aa) || m4_aa == "Not aligned" || nchar(m4_aa) < 6) return(NOT_DETECTED)
  aln <- align_to_reference_frame(motif_AA$M4, m4_aa)
  aa <- query_at_ref_pos(aln, 6)
  if (aa %in% c("G", "S")) return(aa)
  if (aa %in% c("", "X", "*")) return(NOT_DETECTED)
  "Other"
}

call_g1157_from_position <- function(fullAA) {
  aa <- clean_aa_for_alignment(as.character(fullAA))
  if (nchar(aa) < 16) return(NOT_DETECTED)
  # Fallback scans the whole input for the local G1157 context. This keeps
  # fragment inputs usable when their sequence does not start at S position 1.
  window <- aa
  
  # Typical local context: TVLVP[G/S]DFVNVIAI; common variants include
  # DFVDV, DFIDV and non-G/S residues at the focal position.
  pat1 <- "TVL[A-Z]{2}[A-Z]D[FIY][VIL][A-Z][A-Z]IAI"
  hit <- regexpr(pat1, window, perl = TRUE)[[1]]
  if (hit > 0) {
    residue <- substr(window, hit + 5, hit + 5)
    if (residue %in% c("G", "S")) return(residue)
    if (residue %in% c("X", "*", "")) return(NOT_DETECTED)
    return("Other")
  }
  
  pat2 <- "[A-Z]D[FIY][VIL][A-Z][A-Z]IAI"
  hits <- gregexpr(pat2, window, perl = TRUE)[[1]]
  if (hits[1] > 0) {
    for (hit in hits) {
      local <- substr(window, max(1, hit - 12), hit + 10)
      if (grepl("TVL|LV|VP|EP", local, perl = TRUE)) {
        residue <- substr(window, hit, hit)
        if (residue %in% c("G", "S")) return(residue)
        if (residue %in% c("X", "*", "")) return(NOT_DETECTED)
        return("Other")
      }
    }
  }
  
  NOT_DETECTED
}

lookup_type <- function(pattern, map) {
  if (any(grepl(NOT_DETECTED, pattern, fixed = TRUE)) || any(grepl("Other", pattern, fixed = TRUE)))
    return(NOT_AVAILABLE)
  hit <- unname(map[pattern])
  if (length(hit) == 0 || is.na(hit)) return("Unassigned")
  hit[[1]]
}

s1_type_map <- c(
  "NXT|NND" = "European/American non-SINDEL",
  "NXT|DND" = "Asian non-SINDEL",
  "NXS|ND" = "SINDEL",
  "NXS|SG" = "Prevalent SINDEL",
  "NXS|GH" = "European SINDEL",
  "No glycosylation|NDD" = "Korean G2"
)

s2_type_map <- c(
  "NXT|NXT" = "G2c",
  "NXT|No glycosylation" = "G2b",
  "NXS|NXT" = "G1b",
  "NXS|No glycosylation" = "G1a",
  "No glycosylation|No glycosylation" = "G2a"
)

s4_type_map <- c(
  "NXT|DND|G|NXT" = "Prevalent G2c",
  "NXT|NND|S|NXT" = "European/American G2c",
  "NXT|NND|G|NXT" = "Asian G2c",
  "NXT|NND|G|No glycosylation" = "Asian G2b",
  "NXT|DND|G|No glycosylation" = "Prevalent G2b",
  "NXS|ND|G|NXT" = "European/American G1b",
  "NXS|SG|G|NXT" = "Prevalent G1b",
  "NXS|ND|G|No glycosylation" = "G1a",
  "NXS|ND|S|NXT" = "G1b",
  "NXS|GH|G|NXT" = "European G1b",
  "No glycosylation|NDD|G|No glycosylation" = "G2a"
)

h5_type_map <- c(
  "NXT|DND|No glycosylation|G|NXT" = "G2c: H1",
  "NXT|NND|No glycosylation|S|NXT" = "G2c: H2",
  "NXT|NND|No glycosylation|G|NXT" = "G2c: H3",
  "NXT|NND|NXT|G|NXT" = "G2c: H4",
  "NXT|NND|NXT|S|NXT" = "G2c: H5",
  "NXT|NND|No glycosylation|G|No glycosylation" = "G2b: H6",
  "NXT|DND|No glycosylation|G|No glycosylation" = "G2b: H7",
  "NXS|ND|No glycosylation|G|NXT" = "G1b: H8",
  "NXT|NND|NXT|G|No glycosylation" = "G2b: H9",
  "NXT|DND|NXT|G|NXT" = "G2c: H10",
  "NXS|SG|No glycosylation|G|NXT" = "G1b: H11",
  "NXT|DND|NXT|G|No glycosylation" = "G2b: H12",
  "NXS|ND|No glycosylation|S|NXT" = "G1b: H13",
  "NXS|ND|NXT|G|No glycosylation" = "G1a: H14",
  "NXS|GH|No glycosylation|G|NXT" = "G1b: H15",
  "NXS|ND|No glycosylation|G|No glycosylation" = "G1a: H16",
  "No glycosylation|NDD|NXT|G|No glycosylation" = "G2a: H17",
  "NXS|SG|NXT|G|NXT" = "G1b: H18"
)

assign_current_typing <- function(n57, site135, n718, g1157, n1192) {
  s1_pattern <- paste(n57, site135, sep = "|")
  s2_pattern <- paste(n57, n1192, sep = "|")
  s4_pattern <- paste(n57, site135, g1157, n1192, sep = "|")
  h5_pattern <- paste(n57, site135, n718, g1157, n1192, sep = "|")
  
  s1_type <- lookup_type(s1_pattern, s1_type_map)
  s2_type <- lookup_type(s2_pattern, s2_type_map)
  s4_type <- lookup_type(s4_pattern, s4_type_map)
  h5_type <- lookup_type(h5_pattern, h5_type_map)
  
  main_type <- dplyr::case_when(
    h5_type != NOT_AVAILABLE & h5_type != "Unassigned" ~ h5_type,
    s4_type != NOT_AVAILABLE & s4_type != "Unassigned" ~ s4_type,
    s2_type != NOT_AVAILABLE & s2_type != "Unassigned" ~ s2_type,
    s1_type != NOT_AVAILABLE & s1_type != "Unassigned" ~ s1_type,
    TRUE ~ "Unassigned"
  )
  
  list(
    main_type = main_type,
    s1_pattern = s1_pattern,
    s2_pattern = s2_pattern,
    s4_pattern = s4_pattern,
    h5_pattern = h5_pattern,
    s1_type = s1_type,
    s2_type = s2_type,
    s4_type = s4_type,
    h5_type = h5_type
  )
}


###############################
# Analyze one sequence (FIXED)
###############################
analyze_one_sequence <- function(seq_raw, seq_id, motif_AA){
  
  pre <- preprocess_sequence(seq_raw)
  seq_clean <- pre$clean_seq
  seq_type  <- pre$seq_type
  
  if(seq_type=="NT"){
    frames <- translate_three_frames(seq_clean)
    best_frame <- select_best_frame(frames, motif_AA)
    motifs <- best_frame$motifs
    best_aa <- frames[[best_frame$best_frame]]
  } else {
    fullAA <- AAString(clean_aa_for_alignment(seq_clean))
    motifs <- extract_all_motifs(fullAA, motif_AA)
    best_aa <- fullAA
  }
  
  n57_state <- call_n57_n62(motifs$M1$aa)
  if (n57_state == NOT_DETECTED) {
    n57_state <- call_n57_n62_from_window(best_aa)
  }
  site135_state <- call_135_136(motifs$M135$aa)
  if (site135_state == NOT_DETECTED) {
    m135_fallback <- extract_local_region_with_insertions(best_aa, motif_AA$M135, 50)
    site135_state <- call_135_136(m135_fallback$aa)
  }
  if (site135_state == NOT_DETECTED) {
    site135_state <- call_135_136_from_window(best_aa)
  }
  n1192_state <- call_n1192_n1194(motifs$M2$aa)
  g1157_state <- call_g1157(motifs$M4$aa)
  if (g1157_state == NOT_DETECTED) {
    g1157_state <- call_g1157_from_position(best_aa)
  }
  n718_state <- call_n718_n722(motifs$M3$aa)
  
  typing <- assign_current_typing(
    n57 = n57_state,
    site135 = site135_state,
    n718 = n718_state,
    g1157 = g1157_state,
    n1192 = n1192_state
  )
  
  s1_ready <- n57_state != NOT_DETECTED && site135_state != NOT_DETECTED
  full_ready <- s1_ready && n1192_state != NOT_DETECTED && g1157_state != NOT_DETECTED && n718_state != NOT_DETECTED
  
  missing <- c()
  if (n57_state == NOT_DETECTED) missing <- c(missing, "N57/N62")
  if (site135_state == NOT_DETECTED) missing <- c(missing, "135/136")
  if (n1192_state == NOT_DETECTED) missing <- c(missing, "N1192/N1194")
  if (g1157_state == NOT_DETECTED) missing <- c(missing, "G1157")
  if (n718_state == NOT_DETECTED) missing <- c(missing, "N718/N722")
  
  status <- dplyr::case_when(
    !s1_ready ~ "Failed",
    full_ready ~ "Success",
    TRUE ~ "Partial"
  )
  
  note <- dplyr::case_when(
    status == "Success" ~ "All five polymorphic loci were detected.",
    s1_ready ~ paste0(
      "Sequence is too short or lacks downstream coverage; only available loci were typed. Missing loci: ",
      paste(missing, collapse = ", "), "."
    ),
    TRUE ~ paste0("Insufficient S1 coverage for S1 typing. Missing loci: ", paste(missing, collapse = ", "), ".")
  )
  
  tibble(
    sequence_ID = seq_id,
    Status = status,
    Input_Type = seq_type,
    Motif_Type = typing$main_type,
    N57_N62_state = n57_state,
    Site_135_136_motif = site135_state,
    N1192_N1194_state = n1192_state,
    G1157_state = g1157_state,
    N718_N722_state = n718_state,
    S1_two_locus_type = typing$s1_type,
    S1_two_locus_pattern = typing$s1_pattern,
    S_two_locus_type = typing$s2_type,
    S_two_locus_pattern = typing$s2_pattern,
    S_4_locus_geo_type = typing$s4_type,
    S_4_locus_pattern = typing$s4_pattern,
    Five_locus_haplotype = typing$h5_type,
    Five_locus_pattern = typing$h5_pattern,
    Notes = note
  )
}


###############################
# Batch processing
###############################
process_all_sequences <- function(dss){
  
  raw_seq <- attr(dss,"raw_seq")
  raw_id  <- attr(dss,"raw_id")
  out <- vector("list", length(raw_seq))
  
  withProgress(message="Processing sequences...", value=0,{
    for(i in seq_along(raw_seq)){
      incProgress(1/length(raw_seq), detail=paste("Processing",raw_id[i]))
      out[[i]] <- analyze_one_sequence(raw_seq[i],raw_id[i],motif_AA)
    }
  })
  
  bind_rows(out)
}


###############################
# UI
###############################
ui <- fluidPage(
  tags$script(HTML("
Shiny.addCustomMessageHandler('scrollInfoRow', function(message) {
  var table = document.querySelector('#info_table table');
  if (!table) return;

  var rows = table.querySelectorAll('tbody tr');
  var idx = message.row - 1;

  if (rows[idx]) {
    rows[idx].scrollIntoView({behavior: 'smooth', block: 'center'});
  }
});
")),
  
  tags$style(HTML("
  @keyframes flashRow {
    0%   { background-color: #ffeaa7; }
    50%  { background-color: #fab1a0; }
    100% { background-color: #ffeaa7; }
  }

  tr.flash-highlight td {
    animation: flashRow 1s ease-in-out;
  }

  .static-figure-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
    gap: 18px;
    margin-top: 12px;
  }

  .static-figure-card {
    border: 1px solid #d9e2ec;
    border-radius: 6px;
    padding: 12px;
    background: #ffffff;
  }

  .static-figure-card h5 {
    margin-top: 0;
    font-weight: 700;
  }

  .static-figure-card img,
  .tree-panel img {
    width: 100%;
    height: auto;
    border: 1px solid #e5e7eb;
    border-radius: 4px;
    background: #ffffff;
  }

  .guide-card {
    border: 1px solid #d9e2ec;
    border-radius: 6px;
    padding: 14px 16px;
    margin-bottom: 14px;
    background: #ffffff;
  }

  .guide-card h4,
  .guide-card h5 {
    font-weight: 700;
    margin-top: 0;
  }

  .logic-flow {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
    gap: 12px;
    margin: 12px 0 16px 0;
  }

  .logic-step {
    border: 1px solid #b8c7d9;
    border-left: 5px solid #2c3e50;
    border-radius: 6px;
    padding: 12px;
    background: #f8fafc;
    min-height: 118px;
  }

  .logic-step .step-title {
    font-weight: 700;
    color: #1f2937;
    margin-bottom: 6px;
  }

  .logic-step .step-detail {
    color: #4b5563;
    font-size: 13px;
    line-height: 1.35;
  }

  .compact-table {
    width: 100%;
    border-collapse: collapse;
    margin-top: 8px;
    font-size: 13px;
  }

  .compact-table th,
  .compact-table td {
    border: 1px solid #d9e2ec;
    padding: 7px 8px;
    vertical-align: top;
  }

  .compact-table th {
    background: #edf2f7;
    font-weight: 700;
  }

  .overview-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 12px;
    margin-bottom: 14px;
  }

  .overview-card {
    border: 1px solid #d9e2ec;
    border-radius: 6px;
    padding: 12px 14px;
    background: #f8fafc;
  }

  .overview-value {
    font-size: 26px;
    line-height: 1.1;
    font-weight: 800;
    color: #1f2937;
  }

  .overview-label {
    margin-top: 4px;
    color: #52616f;
    font-size: 13px;
    font-weight: 600;
  }
")
             ),
  
  
  theme = bs_theme(bootswatch = "flatly"),
  
  titlePanel("PEDV S Polymorphism Typing Platform"),
  
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        "paste_fasta", "Paste FASTA sequences:",
        width = "100%", height = "200px",
        placeholder = ">seq1\nATGC... or AA..."
      ),
      
      fileInput("seqfile", "Upload FASTA file (up to 20 MB)"),
      tags$div(
        class = "alert alert-info",
        style = "font-size: 13px; line-height: 1.35; padding: 10px 12px;",
        tags$strong("Incomplete sequence handling"),
        tags$p(
          style = "margin: 6px 0 0 0;",
          "Only confidently detected loci are used for typing. ",
          "If only S1 loci are covered, only S1 typing is reported; missing loci are shown as Not detected."
        )
      ),
      actionButton("run", "Start Analysis", class = "btn btn-primary"),
      hr(),
      downloadButton("download_res", "Download Results CSV")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Guide & typing logic",
          br(),
          div(
            class = "guide-card",
            h4("PEDV S polymorphism typing platform"),
            tags$p(
              "This web application identifies PEDV S-gene polymorphic motif types from pasted or uploaded FASTA sequences, ",
              "supports accession/strain lookup in the curated sequence database, and provides static reference panels for the ",
              "spatiotemporal and phylogenetic distribution of the typing systems."
            ),
            tags$p(
              "Input sequences may be nucleotide S genes, amino-acid S proteins, or longer genome fragments containing S. ",
              "The app cleans illegal characters and alignment gaps, translates nucleotide sequences in three reading frames, ",
              "locates local decision frames by reference-motif matching, and reports only confidently detected loci."
            )
          ),
          div(
            class = "guide-card",
            h4("Typing architecture"),
            div(
              class = "logic-flow",
              div(
                class = "logic-step",
                div(class = "step-title", "1. Polymorphic loci"),
                div(class = "step-detail",
                    "Five local states are detected: N57/N62, 135/136 motif, N1192/N1194, G1157, and N718/N722.")
              ),
              div(
                class = "logic-step",
                div(class = "step-title", "2. S1-type"),
                div(class = "step-detail",
                    "S1-type is assigned from N57/N62 plus the 135/136 motif. This is the fallback system for short S1-only sequences.")
              ),
              div(
                class = "logic-step",
                div(class = "step-title", "3. Genotype"),
                div(class = "step-detail",
                    "Genotype is assigned from N57/N62 plus N1192/N1194, linking the N-terminal and S2/CD glycosylation states.")
              ),
              div(
                class = "logic-step",
                div(class = "step-title", "4. Geo-type"),
                div(class = "step-detail",
                    "Geo-type combines N57/N62, 135/136, N1192/N1194, and G1157 to capture geography-associated S-gene signatures.")
              ),
              div(
                class = "logic-step",
                div(class = "step-title", "5. Haplotype"),
                div(class = "step-detail",
                    "Haplotype combines all five loci, adding N718/N722 to resolve finer motif-level variation.")
              )
            ),
            tags$table(
              class = "compact-table",
              tags$thead(
                tags$tr(
                  tags$th("Typing system"),
                  tags$th("Loci used"),
                  tags$th("Output")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("S1-type"),
                  tags$td("N57/N62 + 135/136 motif"),
                  tags$td("S1-domain classification; available for S1-only partial sequences.")
                ),
                tags$tr(
                  tags$td("Genotype"),
                  tags$td("N57/N62 + N1192/N1194"),
                  tags$td("Two-locus S-gene genotype.")
                ),
                tags$tr(
                  tags$td("Geo-type"),
                  tags$td("N57/N62 + 135/136 motif + N1192/N1194 + G1157"),
                  tags$td("Geography-associated S-gene grouping.")
                ),
                tags$tr(
                  tags$td("Haplotype"),
                  tags$td("N57/N62 + 135/136 motif + N1192/N1194 + G1157 + N718/N722"),
                  tags$td("Fine-scale five-locus haplotype.")
                )
              )
            )
          ),
          div(
            class = "guide-card",
            h4("Named typing groups"),
            tags$p(
              "The table below is loaded from the current Information workbook and summarizes the named groups used by the platform."
            ),
            DTOutput("typing_logic_table")
          ),
          div(
            class = "guide-card",
            h4("Sequence processing and locus rules"),
            tags$table(
              class = "compact-table",
              tags$thead(
                tags$tr(
                  tags$th("Step"),
                  tags$th("Rule"),
                  tags$th("Purpose")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("Input cleaning"),
                  tags$td("Whitespace, alignment gaps, illegal symbols, and unsupported characters are removed or normalized before analysis."),
                  tags$td("Allows aligned S genes, unaligned S genes, and longer genome fragments to be submitted through the same interface.")
                ),
                tags$tr(
                  tags$td("Input type detection"),
                  tags$td("The sequence is classified as nucleotide-like or amino-acid-like from its character composition."),
                  tags$td("Determines whether translation is required before motif localization.")
                ),
                tags$tr(
                  tags$td("Translation"),
                  tags$td("Nucleotide input is translated in three reading frames; the best frame is chosen from reference motif hits and stop-codon burden."),
                  tags$td("Improves compatibility with full genomes, S genes, and frame-shifted pasted fragments.")
                ),
                tags$tr(
                  tags$td("Decision-frame localization"),
                  tags$td("Local frames around N57/N62, 135/136, N718/N722, G1157, and N1192/N1194 are located by imperfect reference-motif matching."),
                  tags$td("Keeps site extraction stable when nearby insertions or deletions are present.")
                ),
                tags$tr(
                  tags$td("Conservative missing calls"),
                  tags$td("If a local decision frame is not confidently detected, the locus is reported as Not detected."),
                  tags$td("Prevents missing sequence coverage from being misinterpreted as absence of a motif.")
                )
              )
            ),
            br(),
            tags$table(
              class = "compact-table",
              tags$thead(
                tags$tr(
                  tags$th("Locus"),
                  tags$th("Reported states"),
                  tags$th("Decision logic")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("N57/N62"),
                  tags$td("NXS, NXT, No glycosylation, Not detected"),
                  tags$td("The localized N-terminal frame is scanned for N-linked glycosylation motifs. NXS is prioritized when an NXS motif is present; otherwise NXT is reported when an NXT motif is detected.")
                ),
                tags$tr(
                  tags$td("135/136 motif"),
                  tags$td("SG, ND, NND, DND, NDD, GH, Other, Not detected"),
                  tags$td("The local frame is aligned to the CV777 reference decision frame, and the motif around the 135/136 region is read from the aligned position.")
                ),
                tags$tr(
                  tags$td("N1192/N1194"),
                  tags$td("NXT, No glycosylation, Not detected"),
                  tags$td("The S2/CD decision frame is scanned using N1192/N1194 reference contexts, including glycosylated and non-glycosylated states.")
                ),
                tags$tr(
                  tags$td("G1157"),
                  tags$td("G, S, Other, Not detected"),
                  tags$td("The G1157 frame is aligned to the reference context, and the amino acid at the reference-anchored position is read directly.")
                ),
                tags$tr(
                  tags$td("N718/N722"),
                  tags$td("NXT, No glycosylation, Not detected"),
                  tags$td("The N718/N722 frame is scanned for paired glycosylation signals. If the decision frame is missing, the site remains Not detected.")
                )
              )
            )
          ),
          div(
            class = "guide-card",
            h4("Output interpretation"),
            tags$table(
              class = "compact-table",
              tags$thead(
                tags$tr(
                  tags$th("Output field"),
                  tags$th("Interpretation")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("Five_locus_haplotype"),
                  tags$td("Fine-scale five-locus assignment using all detected polymorphic loci.")
                ),
                tags$tr(
                  tags$td("S_4_locus_geo_type"),
                  tags$td("Four-locus geography-associated type based on N57/N62, 135/136, N1192/N1194, and G1157.")
                ),
                tags$tr(
                  tags$td("S_two_locus_type"),
                  tags$td("Two-locus S-gene genotype based on N57/N62 and N1192/N1194.")
                ),
                tags$tr(
                  tags$td("S1_two_locus_type"),
                  tags$td("S1-region type based on N57/N62 and the 135/136 motif.")
                ),
                tags$tr(
                  tags$td("No glycosylation"),
                  tags$td("The local decision frame was detected, but the required glycosylation motif was absent.")
                ),
                tags$tr(
                  tags$td("Not detected"),
                  tags$td("The local decision frame was missing or not confidently detected; the state is not inferred.")
                ),
                tags$tr(
                  tags$td("Undetermined"),
                  tags$td("The corresponding typing system cannot be assigned because one or more required loci are missing.")
                )
              )
            )
          ),
          div(
            class = "guide-card",
            h4("How to use"),
            tags$ol(
              tags$li("Use the Sequence typing tab to paste or upload FASTA sequences and download the typing results."),
              tags$li("Use the Database search tab to find curated strains by accession, strain name, region, year, or typing result."),
              tags$li("Use the Information tab to inspect reference descriptions, literature links, static spatiotemporal panels, and tree distributions.")
            ),
            tags$p(
              tags$strong("Partial sequences: "),
              "If only S1 loci are detected, the app reports S1-type and marks downstream systems as Undetermined. ",
              "If a decision frame is not confidently detected, that locus is reported as Not detected rather than inferred as a negative state."
            )
          )
        ),
        tabPanel(
          "Sequence typing",
          br(),
          fluidRow(
            column(
              width = 12,
              withSpinner(
                DTOutput("res_table"),
                type = 6, color = "#2c3e50"
              )
            )
          )
        ),
        tabPanel(
          "Database search",
          br(),
          fluidRow(
            column(
              width = 12,
              uiOutput("db_overview_cards"),
              textInput(
                "db_query",
                "Search by accession, sequence name, strain, country, region, or typing result:",
                value = "",
                placeholder = "e.g. AB857233, CV777, GDS01, China, G2b"
              ),
              withSpinner(
                DTOutput("db_table"),
                type = 6, color = "#2c3e50"
              ),
              tags$div(
                style = "margin-top: 10px;",
                downloadButton("download_db_filtered", "Download Search Results CSV"),
                tags$span(" "),
                downloadButton("download_db_full", "Download Full Database CSV")
              )
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 12,
              h4("Selected strain"),
              DTOutput("db_selected_table")
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 12,
              h4("Group summary"),
              DTOutput("db_group_summary_table")
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 12,
              h4("Typing distribution"),
              uiOutput("db_group_panel")
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 12,
              h4("Phylogenetic distribution"),
              uiOutput("db_tree_panel")
            )
          )
        ),
        tabPanel(
          "Information",
          br(),
          fluidRow(
            column(
              width = 12,
              h4("Information"),
              div(
                style = "height: 350px; overflow-y: auto;",
                div(
                  id = "info_section",
                  DTOutput("info_table")
                )
              )
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 12,
              h4("Static reference panels"),
              uiOutput("static_info_panel")
            )
          )
        ),
        tabPanel(
          "Validation & thesis use",
          br(),
          div(
            class = "guide-card",
            h4("Planned validation workflow"),
            tags$p(
              "This section defines the validation framework for the typing platform. ",
              "It is intended for documenting independent test-set performance before the website is used as a thesis result."
            ),
            tags$table(
              class = "compact-table",
              tags$thead(
                tags$tr(
                  tags$th("Validation item"),
                  tags$th("Evidence to collect"),
                  tags$th("Expected output")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("Input compatibility"),
                  tags$td("Test nucleotide S genes, amino-acid S proteins, aligned sequences, unaligned sequences, and genome-length fragments."),
                  tags$td("A compatibility summary table listing sequence type, status, and detected loci.")
                ),
                tags$tr(
                  tags$td("Locus detection"),
                  tags$td("Count detected, not detected, and ambiguous states for N57/N62, 135/136, N1192/N1194, G1157, and N718/N722."),
                  tags$td("A locus-level detection-rate table.")
                ),
                tags$tr(
                  tags$td("Typing consistency"),
                  tags$td("Compare web output with curated database labels or manually checked reference labels."),
                  tags$td("An agreement table for S1-type, genotype, geo-type, and haplotype.")
                ),
                tags$tr(
                  tags$td("Partial sequence behavior"),
                  tags$td("Submit sequences covering only S1, only downstream S regions, and truncated S genes."),
                  tags$td("A missing-locus report confirming that unavailable loci are not converted into false negative states.")
                ),
                tags$tr(
                  tags$td("Database search"),
                  tags$td("Query representative accessions, strain names, regions, years, and typing labels."),
                  tags$td("A search-function checklist and exported search-result examples.")
                )
              )
            )
          ),
          div(
            class = "guide-card",
            h4("Thesis result structure"),
            tags$ol(
              tags$li("Describe the website construction: Shiny framework, input parser, motif localization, typing systems, and database resources."),
              tags$li("Report validation data: test-set composition, success rate, locus detection rate, and typing consistency."),
              tags$li("Summarize website application: database search, reference information, static spatiotemporal panels, and tree distribution panels."),
              tags$li("Discuss utility and limitations: support for partial sequences, conservative missing calls, and the need to update the database with new strains.")
            )
          ),
          div(
            class = "guide-card",
            h4("Recommended validation files"),
            tags$table(
              class = "compact-table",
              tags$thead(
                tags$tr(
                  tags$th("File"),
                  tags$th("Content")
                )
              ),
              tags$tbody(
                tags$tr(
                  tags$td("validation_input.fasta"),
                  tags$td("Independent sequences used to test website typing.")
                ),
                tags$tr(
                  tags$td("validation_expected_labels.csv"),
                  tags$td("Manual or curated reference labels for expected locus states and typing assignments.")
                ),
                tags$tr(
                  tags$td("validation_web_output.csv"),
                  tags$td("Downloaded website output for the validation input.")
                ),
                tags$tr(
                  tags$td("validation_summary.csv"),
                  tags$td("Detection rates, typing consistency, and failed-case summaries.")
                ),
                tags$tr(
                  tags$td("validation_result_and_analysis.md"),
                  tags$td("Thesis-ready result and analysis draft based on validation output.")
                )
              )
            )
          )
        )
      )
    )
  )
)


###############################
# UI override: two-entry application layout
###############################
ui <- fluidPage(
  tags$script(HTML("
    Shiny.addCustomMessageHandler('scrollToResult', function(message) {
      var target = document.getElementById('result_anchor');
      if (target) {
        target.scrollIntoView({ behavior: 'smooth', block: 'start' });
      }
    });
  ")),
  tags$style(HTML("
    body { background: #f5f7fb; }
    .app-title { margin: 6px 0 0 0; font-weight: 800; color: #1f2937; font-size: 23px; }
    .app-subtitle { color: #52616f; margin-bottom: 6px; font-size: 12px; }
    .main-row { display: flex; align-items: stretch; min-height: calc(100vh - 62px); }
    .left-stack, .right-stack { display: flex; flex-direction: column; min-height: calc(100vh - 62px); }
    .right-stack > .right-panel, .right-stack > .result-box { min-height: 0; overflow: hidden; display: flex; flex-direction: column; }
    .right-stack > .right-panel { flex: 1 1 0; }
    .right-stack > .result-box { flex: 3 1 0; }
    .right-stack > .trend-box { flex: 4 1 0; }
    .left-stack > .tree-box { flex: 1 1 auto; min-height: 0; overflow: hidden; display: flex; flex-direction: column; }
    .left-panel, .right-panel, .tree-box, .result-box {
      border: 1px solid #d9e2ec;
      border-radius: 6px;
      background: #ffffff;
      padding: 7px 10px;
      margin-bottom: 7px;
    }
    .left-panel h4, .right-panel h4, .tree-box h4, .result-box h4 {
      margin-top: 0;
      font-weight: 800;
      color: #1f2937;
    }
    .mini-help {
      border: 1px solid #d9e2ec;
      border-left: 4px solid #2c3e50;
      border-radius: 6px;
      padding: 7px 10px;
      background: #f8fafc;
      color: #4b5563;
      font-size: 12px;
      line-height: 1.25;
      margin-bottom: 7px;
    }
    .two-entry-grid { display: grid; grid-template-columns: 1fr; gap: 7px; }
    .tree-box img { width: 100%; max-height: none; height: 100%; object-fit: contain; border: 1px solid #e5e7eb; border-radius: 4px; background: #ffffff; }
    .tree-box .shiny-html-output { flex: 1 1 auto; min-height: 0; display: flex; flex-direction: column; }
    .trend-plot-wrap { flex: 1 1 auto; min-height: 0; }
    .trend-plot-wrap .shiny-plot-output { height: 100% !important; min-height: 210px; }
    .compact-table { width: 100%; border-collapse: collapse; margin-top: 4px; font-size: 11px; }
    .compact-table th, .compact-table td { border: 1px solid #d9e2ec; padding: 3px 5px; vertical-align: top; }
    .compact-table th { background: #edf2f7; font-weight: 700; }
    .overview-grid { display: grid; grid-template-columns: repeat(4, minmax(80px, 1fr)); gap: 6px; margin-bottom: 4px; }
    .overview-card { border: 1px solid #d9e2ec; border-radius: 6px; padding: 5px 7px; background: #f8fafc; }
    .overview-value { font-size: 16px; line-height: 1; font-weight: 800; color: #1f2937; }
    .overview-label { margin-top: 2px; color: #52616f; font-size: 10px; font-weight: 600; }
    .shiny-input-container { margin-bottom: 6px; }
    h4 { font-size: 15px; margin-bottom: 7px; }
    h5 { font-size: 13px; font-weight: 800; margin-top: 4px; margin-bottom: 4px; }
    hr { margin-top: 6px; margin-bottom: 6px; }
    .dataTables_wrapper { font-size: 13px; }
    .result-box .form-group { margin-bottom: 4px; }
    .result-box .radio { margin-top: 0; margin-bottom: 0; }
    .result-box table.dataTable { width: 100% !important; }
    .result-box table.dataTable tbody td, .result-box table.dataTable thead th {
      padding: 7px 8px;
      font-size: 13px;
      line-height: 1.32;
      white-space: normal;
      vertical-align: top;
    }
    .result-box table.dataTable thead th { font-size: 13px; font-weight: 800; }
    .result-box .dataTables_paginate {
      padding-top: 4px;
      text-align: center;
      font-size: 12px;
      font-weight: 700;
    }
    .result-box .dataTables_paginate .paginate_button {
      padding: 2px 7px !important;
      margin: 0 2px;
      border-radius: 4px !important;
    }
    .left-panel .input-group { width: 100%; max-width: 100%; }
    .left-panel .input-group .form-control { min-width: 0; font-size: 12px; }
    .left-panel .btn { max-width: 100%; white-space: normal; }
    .button-row { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 6px; }
    .button-row .btn { width: 100%; font-size: 12px; padding-left: 4px; padding-right: 4px; }
    .logic-flow { display: grid; grid-template-columns: minmax(66px, 0.72fr) 12px minmax(70px, 0.72fr) 12px minmax(86px, 0.86fr) 12px minmax(160px, 1.7fr); gap: 3px; align-items: stretch; }
    .logic-step { border: 1px solid #d9e2ec; border-radius: 6px; background: #f8fafc; padding: 5px; min-height: 70px; }
    .logic-step-title { font-weight: 800; color: #1f2937; margin-bottom: 3px; font-size: 13px; }
    .logic-step-loci { font-size: 9.5px; color: #52616f; line-height: 1.16; overflow-wrap: anywhere; }
    .logic-arrow { display: flex; align-items: center; justify-content: center; color: #2c3e50; font-size: 16px; font-weight: 800; }
    .logic-note { margin-top: 5px; color: #52616f; font-size: 11px; }
    .result-box .selectize-input { min-height: 32px; padding: 5px 8px; white-space: nowrap; overflow: hidden; }
    .result-box .selectize-input > * { white-space: nowrap; }
    .trend-box label,
    .trend-box .control-label,
    .trend-box .selectize-input,
    .trend-box .selectize-dropdown,
    .trend-box .btn {
      font-size: 13px;
    }
    .tree-switch { display: grid; grid-template-columns: repeat(4, 1fr); gap: 5px; margin-bottom: 6px; }
    .tree-switch .btn { font-size: 12px; padding: 5px 4px; white-space: normal; }
    .tree-switch .active-tree-btn { background: #2c3e50; color: #ffffff; border-color: #2c3e50; }
  ")),
  theme = bs_theme(bootswatch = "flatly"),
  fluidRow(
    column(
      width = 12,
      h2(class = "app-title", "PEDV S Polymorphism Typing Platform"),
      div(class = "app-subtitle", "Sequence typing and strain lookup based on PEDV S-gene polymorphic motifs.")
    )
  ),
  fluidRow(
    class = "main-row",
    column(
      width = 5,
      class = "left-stack",
      div(
        class = "left-panel",
        h4("Sequence typing or strain search"),
        textAreaInput(
          "paste_fasta", "Paste FASTA sequence, or enter accession/strain",
          width = "100%", height = "115px",
          placeholder = ">seq1\nATGC... / amino-acid sequence...\n\nor: CV777, AB857233, JN599150"
        ),
        fileInput("seqfile", "Upload FASTA file (up to 20 MB)"),
        div(
          class = "button-row",
          actionButton("run", "Analyze / Search", class = "btn btn-primary"),
          downloadButton("download_res", "Typing CSV"),
          downloadButton("download_db_filtered", "Search CSV")
        )
      ),
      div(
        class = "mini-help",
        tags$strong("Quick guide. "),
        "Paste/upload a sequence for de novo typing, or enter a known accession/strain for database search. ",
        "For incomplete sequences, only confidently detected loci are reported."
      ),
      div(
        class = "tree-box",
        h4("Phylogenetic context"),
        uiOutput("active_tree_panel")
      )
    ),
    column(
      width = 7,
      class = "right-stack",
      div(
        class = "right-panel",
        h4("Guide and typing logic"),
        uiOutput("concise_guide_panel")
      ),
      div(
        id = "result_anchor",
        class = "result-box",
        h4("Results and group information"),
        withSpinner(DTOutput("unified_result_table"), type = 6, color = "#2c3e50"),
        h5("Reference strains and group description"),
        DTOutput("info_table")
      ),
      div(
        class = "result-box trend-box",
        h4("Spatiotemporal dynamics"),
        fluidRow(
          column(width = 3, uiOutput("trend_region_filter")),
          column(width = 3, uiOutput("trend_country_filter")),
          column(width = 3, uiOutput("trend_system_filter")),
          column(width = 3, uiOutput("trend_type_filter"))
        ),
        div(class = "trend-plot-wrap", plotOutput("dynamic_trend_plot", height = "100%"))
      )
    )
  )
)


###############################
# SERVER (FINAL 鈥?FILTER MODE)
###############################
server <- function(input, output, session) {
  
  ###############################
  # 鏂囨湰瑙勮寖鍖栧嚱鏁帮紙鐢ㄤ簬绋冲仴鍖归厤锛?
  ###############################
  normalize_motif <- function(x) {
    x <- as.character(x)
    x <- gsub("<[^>]+>", "", x)                 # 鍘?HTML
    x <- gsub("\u00A0", " ", x, fixed = TRUE)   # 鍘?NBSP
    x <- gsub("\\s+", " ", x)                   # 澶氱┖鏍煎帇缂?
    trimws(x)
  }

  result_display <- function(df) {
    df %>%
      dplyr::select(
        sequence_ID,
        Five_locus_haplotype,
        S_4_locus_geo_type,
        S_two_locus_type,
        S1_two_locus_type,
        N57_N62_state,
        Site_135_136_motif,
        N1192_N1194_state,
        G1157_state,
        N718_N722_state
      )
  }
  
  ###############################
  # 淇濆瓨鈥滆鐐瑰嚮鐨?Motif_Type鈥?  ###############################
  selected_motif <- reactiveVal(NULL)

  ###############################
  # Sequence database search
  ###############################
  translate_display_terms <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("Not available", "Undetermined", x, fixed = TRUE)
    x <- gsub("\u4e2d\u56fd\u5927\u9646", "Chinese Mainland", x, fixed = TRUE)
    x <- gsub("Mainland China", "Chinese Mainland", x, fixed = TRUE)
    x <- gsub("\u4e2d\u56fd\u53f0\u6e7e", "Chinese Taiwan", x, fixed = TRUE)
    x <- gsub("China Taiwan", "Chinese Taiwan", x, fixed = TRUE)
    x
  }

  translate_description_terms <- function(x) {
    x <- translate_display_terms(x)
    from <- c(
      "\u4e9a\u6d32\u591a\u4e2a\u56fd\u5bb62010\u5e74\u540e\u66b4\u53d1\u7684\u9ad8\u81f4\u75c5\u6027non-SINDEL\u6bd2\u682a",
      "\u6b27\u6d32\u6d41\u884c\u7684SINDEL\u6bd2\u682a",
      "2013\u5e74\u540e\u5728\u5317\u7f8e\u3001\u97e9\u56fd\u3001\u65e5\u672c\u3001Chinese Taiwan\u7b49\u5730\u533a\u6d41\u884c\u7684\u9ad8\u81f4\u75c5\u6027non-SINDEL\u6bd2\u682a",
      "\u97e9\u56fd2010\u5e74\u524d\u6d41\u884c\u7684\u6bd2\u682a\uff0c\u4f4d\u4e8e\u72ec\u7acb\u8fdb\u5316\u5206\u652f\uff0c\u53ef\u80fd\u662fSINDEL\u5230non-SINDEL\u6bd2\u682a\u7684\u8fc7\u6e21",
      "\u5728Chinese Mainland\u6d41\u884c\u7684SINDEL\u6bd2\u682a",
      "\u5305\u62ec\u7ecf\u5178\u6bd2\u682a\u5728\u5185\u7684SINDEL\u6bd2\u682a",
      "\u7ecf\u5178\u6bd2\u682a\u7fa4\uff1b",
      "\u4e2d\u7b49\u81f4\u75c5\u6bd2\u682a\uff0c\u65e9\u671f\u5206\u5316\u4e8e\u7ecf\u5178\u6bd2\u682a\uff1b",
      "\u97e9\u56fd2010\u5e74\u524d\u6d41\u884c\u7684\u6bd2\u682a\uff0c\u4f4d\u4e8e\u72ec\u7acb\u8fdb\u5316\u5206\u652f\uff0c\u53ef\u80fd\u662fG1\u6bd2\u682a\u5230G2\u6bd2\u682a\u7684\u8fc7\u6e21\uff1b",
      "\u9ad8\u81f4\u75c5\u6027\u6bd2\u682a\uff0c\u57282010\u540e\u4e2d\u56fd\u548c\u4e1c\u5357\u4e9a\u56fd\u5bb6\u6d41\u884c\uff1b",
      "\u9ad8\u81f4\u75c5\u6027\u6bd2\u682a\uff0c\u57282010\u5e74\u4e2d\u56fd\u66b4\u53d1\u65e9\u671f\u6bd2\u682a\u4e2d\u5360\u636e\u4e3b\u5bfc\uff1b",
      "\u5728\u6b27\u6d32\u6d41\u884c\u7684\u5206\u5316\u6bd2\u682a\uff1b",
      "2010\u5e74\u540e\u5728\u591a\u4e2a\u5730\u533a\u6d41\u884c\uff1b",
      "2013\u5e74\u540e\u5728\u5317\u7f8e\u3001\u97e9\u56fd\u3001\u65e5\u672c\u3001Chinese Taiwan\u7b49\u5730\u533a\u66b4\u53d1\u4e2d\u5360\u636e\u4e3b\u5bfc\uff1b",
      "2018\u5e74\u540e\u5728Chinese Mainland\u518d\u6b21\u51fa\u73b0\u5206\u5316\u5e76\u6d41\u884c\uff1b",
      "2016\u5e74\u540e\u5728Chinese Mainland\u6d41\u884c\uff0c\u53ef\u80fd\u8d77\u6e90\u4e8e\u91cd\u7ec4\u76f8\u5173\u7684\u8c31\u7cfb\u5206\u5316\uff1b",
      "2014\u5e74\u540e\u5728Chinese Mainland\u6d41\u884c\uff0c\u53ef\u80fd\u8d77\u6e90\u4e8e\u91cd\u7ec4\u76f8\u5173\u7684\u8c31\u7cfb\u5206\u5316",
      "N718/N722\u7cd6\u57fa\u5316\u6d88\u5931"
    )
    to <- c(
      "Highly pathogenic non-SINDEL strains that emerged in multiple Asian countries after 2010",
      "SINDEL strains circulating in Europe",
      "Highly pathogenic non-SINDEL strains circulating in North America, South Korea, Japan, Chinese Taiwan, and related regions after 2013",
      "Strains circulating in South Korea before 2010, forming an independent branch and possibly representing a transition from SINDEL to non-SINDEL strains",
      "SINDEL strains circulating in Chinese Mainland",
      "SINDEL strains including classical strains",
      "Classical strain group",
      "Moderately pathogenic strains that diverged early from classical strains",
      "Strains circulating in South Korea before 2010, forming an independent branch and possibly representing a transition from G1 to G2 strains",
      "Highly pathogenic strains circulating in China and Southeast Asian countries after 2010",
      "Highly pathogenic strains that dominated early outbreak strains in China in 2010",
      "Divergent strains circulating in Europe",
      "Strains circulating across multiple regions after 2010",
      "Dominant outbreak strains in North America, South Korea, Japan, Chinese Taiwan, and related regions after 2013",
      "A divergent lineage that re-emerged and circulated in Chinese Mainland after 2018",
      "A lineage circulating in Chinese Mainland after 2016, possibly derived from recombination-associated lineage divergence",
      "A lineage circulating in Chinese Mainland after 2014, possibly derived from recombination-associated lineage divergence",
      "N718/N722 glycosylation loss"
    )
    for (i in seq_along(from)) x <- gsub(from[i], to[i], x, fixed = TRUE)
    x
  }

  database_data <- reactive({
    validate(need(file.exists("sequence_database.csv"),
                  "sequence_database.csv not found"))
    db <- read.csv(
      "sequence_database.csv",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fileEncoding = "UTF-8-BOM"
    )
    if ("Country" %in% names(db)) {
      db$Country <- ifelse(
        grepl("Taiwan|鍙版咕", db$Country, ignore.case = TRUE),
        "China",
        db$Country
      )
    }
    if ("Region" %in% names(db)) {
      region_map <- c(
        "\u4e2d\u56fd\u5927\u9646" = "Chinese Mainland",
        "\u4e2d\u56fd\u53f0\u6e7e" = "Chinese Taiwan",
        "\u97e9\u56fd" = "South Korea",
        "\u65e5\u672c" = "Japan",
        "\u7f8e\u6d32" = "Americas",
        "\u6b27\u6d32" = "Europe",
        "\u4e1c\u5357\u4e9a" = "Southeast Asia"
      )
      hit <- match(db$Region, names(region_map))
      db$Region[!is.na(hit)] <- unname(region_map[hit[!is.na(hit)]])
    }
    char_cols <- names(db)[vapply(db, is.character, logical(1))]
    db[char_cols] <- lapply(db[char_cols], translate_display_terms)
    db
  })

  is_probable_sequence_text <- function(x) {
    x <- trimws(as.character(x))
    if (!nzchar(x)) return(FALSE)
    if (startsWith(x, ">")) return(TRUE)
    compact <- toupper(gsub("[^A-Z\\-]", "", x))
    if (nchar(compact) < 80) return(FALSE)
    nt_chars <- gregexpr("[ACGTRYSWKMBDHVN\\-]", compact, perl = TRUE)[[1]]
    nt_n <- if (nt_chars[1] < 0) 0 else length(nt_chars)
    aa_chars <- gregexpr("[ARNDCQEGHILKMFPSTWYVBZJXUO\\-]", compact, perl = TRUE)[[1]]
    aa_n <- if (aa_chars[1] < 0) 0 else length(aa_chars)
    max(nt_n, aa_n) / nchar(compact) >= 0.90
  }

  search_query_text <- reactive({
    q <- if (is.null(input$paste_fasta)) "" else trimws(input$paste_fasta)
    if (!nzchar(q) || is_probable_sequence_text(q)) return("")
    q
  })

  db_filtered <- reactive({
    db <- database_data()
    q <- search_query_text()
    if (!nzchar(q)) return(utils::head(db, 5))

    searchable <- intersect(
      c(
        "Accession", "Sequence ID", "Strain/Isolate"
      ),
      names(db)
    )
    haystack <- apply(db[, searchable, drop = FALSE], 1, paste, collapse = " | ")
    db[grepl(q, haystack, ignore.case = TRUE, fixed = TRUE), , drop = FALSE]
  })

  selected_database_row <- reactive({
    db <- db_filtered()
    if (nrow(db) == 0) return(NULL)
    idx <- input$db_table_rows_selected
    if (length(idx) == 0) {
      if (nrow(db) == 1) return(db[1, , drop = FALSE])
      return(NULL)
    }
    db[idx[1], , drop = FALSE]
  })

  value_or_na <- function(row, nm) {
    if (is.null(row) || !(nm %in% names(row))) return(NOT_AVAILABLE)
    x <- as.character(row[[nm]][1])
    if (!nzchar(x)) NOT_AVAILABLE else translate_display_terms(x)
  }

  first_valid_type <- function(row) {
    candidates <- c(
      "Five_locus_haplotype", "Haplotype",
      "S_4_locus_geo_type", "Geo-type",
      "S_two_locus_type", "Genotype",
      "S1_two_locus_type", "S1-type"
    )
    for (nm in candidates) {
      if (nm %in% names(row)) {
        x <- as.character(row[[nm]][1])
        if (nzchar(x) && !(x %in% c(NOT_AVAILABLE, NOT_DETECTED, "Unassigned"))) {
          return(x)
        }
      }
    }
    NULL
  }

  standardize_db_results <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame(Message = "No matched database records.", check.names = FALSE))
    }
    out <- data.frame(
      sequence_ID = df[["Sequence ID"]],
      Five_locus_haplotype = df[["Haplotype"]],
      S_4_locus_geo_type = df[["Geo-type"]],
      S_two_locus_type = df[["Genotype"]],
      S1_two_locus_type = df[["S1-type"]],
      N57_N62_state = df[["N57/N62"]],
      Site_135_136_motif = df[["135/136 motif"]],
      N1192_N1194_state = df[["N1192/N1194"]],
      G1157_state = df[["G1157"]],
      N718_N722_state = df[["N718/N722"]],
      check.names = FALSE
    )
    out
  }

  unified_result_data <- reactive({
    q <- search_query_text()
    if (nzchar(q)) {
      return(standardize_db_results(db_filtered()))
    }
    res <- tryCatch(results(), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) {
      return(result_display(res))
    }
    data.frame(Message = "Paste a FASTA sequence or search an accession/strain to display results.", check.names = FALSE)
  })

  selected_unified_row <- reactive({
    df <- unified_result_data()
    if (nrow(df) == 0 || "Message" %in% names(df)) return(NULL)
    idx <- input$unified_result_table_rows_selected
    if (length(idx) == 0) {
      if (nrow(df) == 1) return(df[1, , drop = FALSE])
      return(NULL)
    }
    df[idx[1], , drop = FALSE]
  })

  static_distribution_cards <- function() {
    div(
      class = "static-figure-grid",
      div(
        class = "static-figure-card",
        h5("S1-type spatiotemporal distribution"),
        tags$img(src = "spatiotemporal_S1_2site_type.png")
      ),
      div(
        class = "static-figure-card",
        h5("Genotype spatiotemporal distribution"),
        tags$img(src = "spatiotemporal_S_2site_type.png")
      ),
      div(
        class = "static-figure-card",
        h5("Geo-type spatiotemporal distribution"),
        tags$img(src = "spatiotemporal_S_4site_geotype.png")
      ),
      div(
        class = "static-figure-card",
        h5("Haplotype spatiotemporal distribution"),
        tags$img(src = "spatiotemporal_five_site_haplotype.png")
      )
    )
  }

  static_tree_panel <- function() {
    div(
      class = "tree-panel",
      tags$p(
        "The static tree summarizes the n500 PARNAS-pruned dataset colored by the four typing systems."
      ),
      tags$img(src = "tree_four_typing_combined_new_n500.png")
    )
  }
  
  ###############################
  # 1锔忊儯 Analysis results锛堜笉鍔級
  ###############################
  results <- eventReactive(input$run, {
    
    validate(
      need(input$paste_fasta != "" || !is.null(input$seqfile),
           "Please paste or upload FASTA sequences.")
    )
    
    if (input$paste_fasta != "" && !is_probable_sequence_text(input$paste_fasta)) {
      session$sendCustomMessage("scrollToResult", list())
      return(NULL)
    }

    if (input$paste_fasta != "") {
      txt <- strsplit(input$paste_fasta, "\n")[[1]]
      dss <- parse_fasta_raw(txt)
    } else {
      dss <- parse_fasta_file(input$seqfile$datapath)
    }
    
    validate(need(!is.null(dss), "Invalid FASTA format"))
    process_all_sequences(dss)
  })
  
  ###############################
  # 2锔忊儯 Results table锛堣鐐瑰嚮锛?
  ###############################
  output$res_table <- renderDT({
    if (is.null(input$run) || input$run == 0) {
      return(datatable(
        data.frame(Message = "No sequence submitted yet.", check.names = FALSE),
        rownames = FALSE,
        options = list(
          paging = FALSE,
          searching = FALSE,
          info = FALSE,
          autoWidth = FALSE,
          scrollX = TRUE,
          dom = "t"
        )
      ))
    }
    req(results())
    
    datatable(
      result_display(results()),
      selection = "single",   # 鈽?琛岀偣鍑伙紙鏈€绋冲畾锛?
      options = list(
        pageLength = 2,
        lengthChange = FALSE, # 鈽?涓嶅厑璁哥敤鎴锋敼椤靛ぇ灏?
        paging = TRUE,        # 鈽?寮€鍚垎椤?
        searching = FALSE,
        info = FALSE,
        autoWidth = FALSE,
        scrollX = TRUE,
        scrollY = "95px",
        dom = "t"
      )
    )
  })
  
  ###############################
  # 3锔忊儯 琛岀偣鍑?鈫?璁板綍 Motif_Type
  ###############################
  observeEvent(input$res_table_rows_selected, {
    
    idx <- input$res_table_rows_selected
    if (length(idx) == 0) return()
    
    motif_clicked <- normalize_motif(results()$Motif_Type[idx])
    if (motif_clicked == "") return()
    
    selected_motif(motif_clicked)
    
  }, ignoreInit = TRUE)

  observeEvent(results(), {
    res <- results()
    if (!is.null(res) && nrow(res) == 1) {
      selected_motif(normalize_motif(res$Motif_Type[1]))
    }
    session$sendCustomMessage("scrollToResult", list())
  }, ignoreInit = TRUE)

  observeEvent(input$db_table_rows_selected, {
    row <- selected_database_row()
    if (is.null(row)) return()
    selected_motif(normalize_motif(row$Haplotype[1]))
    session$sendCustomMessage("scrollToResult", list())
  }, ignoreInit = TRUE)

  observeEvent(input$unified_result_table_rows_selected, {
    row <- selected_unified_row()
    motif_clicked <- first_valid_type(row)
    if (is.null(motif_clicked)) {
      selected_motif(NULL)
      return()
    }
    selected_motif(normalize_motif(motif_clicked))
  }, ignoreInit = TRUE)

  observeEvent(unified_result_data(), {
    df <- unified_result_data()
    if ("Message" %in% names(df)) {
      selected_motif(NULL)
      return()
    }
    if (!("Message" %in% names(df)) && nrow(df) == 1) {
      motif_clicked <- first_valid_type(df[1, , drop = FALSE])
      if (!is.null(motif_clicked)) {
        selected_motif(normalize_motif(motif_clicked))
      } else {
        selected_motif(NULL)
      }
    }
  }, ignoreInit = TRUE)

  valid_type_label <- function(x) {
    x <- as.character(x)
    nzchar(x) && !(x %in% c(NOT_AVAILABLE, NOT_DETECTED, "Unassigned", "Not available", "Undetermined", "Not detected", ""))
  }

  active_context <- reactive({
    row <- selected_unified_row()
    if (!is.null(row)) {
      return(list(
        source = "Unified result",
        id = value_or_na(row, "sequence_ID"),
        haplotype = value_or_na(row, "Five_locus_haplotype"),
        geo = value_or_na(row, "S_4_locus_geo_type"),
        genotype = value_or_na(row, "S_two_locus_type"),
        s1 = value_or_na(row, "S1_two_locus_type")
      ))
    }
    NULL
  })

  active_system <- reactive({
    ctx <- active_context()
    if (is.null(ctx)) return(list(system = "Geo-type", value = NULL, tree = "Geo_tree_en_right.svg"))
    if (valid_type_label(ctx$haplotype)) return(list(system = "Haplotype", value = ctx$haplotype, tree = "tree_five_site_haplotype.png"))
    if (valid_type_label(ctx$geo)) return(list(system = "Geo-type", value = ctx$geo, tree = "Geo_tree_en_right.svg"))
    if (valid_type_label(ctx$genotype)) return(list(system = "Genotype", value = ctx$genotype, tree = "tree_S_2site_type.png"))
    if (valid_type_label(ctx$s1)) return(list(system = "S1-type", value = ctx$s1, tree = "S1_tree_en_right.svg"))
    list(system = "Geo-type", value = NULL, tree = "Geo_tree_en_right.svg")
  })

  selected_tree_view <- reactiveVal("Geo-type")
  observeEvent(input$tree_geo, selected_tree_view("Geo-type"), ignoreInit = TRUE)
  observeEvent(input$tree_genotype, selected_tree_view("Genotype"), ignoreInit = TRUE)
  observeEvent(input$tree_haplotype, selected_tree_view("Haplotype"), ignoreInit = TRUE)
  observeEvent(input$tree_s1, selected_tree_view("S1-type"), ignoreInit = TRUE)
  observeEvent(unified_result_data(), selected_tree_view("Geo-type"), ignoreInit = TRUE)

  output$active_tree_panel <- renderUI({
    ctx <- active_context()
    view <- selected_tree_view()
    value <- NULL
    tree <- "Geo_tree_en_right.svg"
    if (identical(view, "Haplotype")) {
      value <- if (!is.null(ctx) && valid_type_label(ctx$haplotype)) ctx$haplotype else NULL
      tree <- "tree_five_site_haplotype.png"
    } else if (identical(view, "Genotype")) {
      value <- if (!is.null(ctx) && valid_type_label(ctx$genotype)) ctx$genotype else NULL
      tree <- "Genotype_tree_en_right.svg"
    } else if (identical(view, "S1-type")) {
      value <- if (!is.null(ctx) && valid_type_label(ctx$s1)) ctx$s1 else NULL
      tree <- "S1_tree_en_right.svg"
    } else {
      view <- "Geo-type"
      value <- if (!is.null(ctx) && valid_type_label(ctx$geo)) ctx$geo else NULL
      tree <- "Geo_tree_en_right.svg"
    }
    label <- if (is.null(value)) paste0("Default: global ", view, " tree") else paste0(view, ": ", value)
    btn_class <- function(x) {
      paste("btn btn-default", if (identical(view, x)) "active-tree-btn" else "")
    }
    tagList(
      div(
        class = "tree-switch",
        actionButton("tree_geo", "Geo-type", class = btn_class("Geo-type")),
        actionButton("tree_genotype", "Genotype", class = btn_class("Genotype")),
        actionButton("tree_haplotype", "Haplotype", class = btn_class("Haplotype")),
        actionButton("tree_s1", "S1-type", class = btn_class("S1-type"))
      ),
      tags$p(style = "color:#52616f; font-size:13px;", label),
      tags$img(src = tree)
    )
  })

  output$unified_result_table <- renderDT({
    df <- unified_result_data()
    datatable(
      df,
      rownames = FALSE,
      selection = if ("Message" %in% names(df)) "none" else "single",
      options = list(
        pageLength = 10,
        lengthChange = FALSE,
        paging = nrow(df) > 10,
        searching = FALSE,
        info = FALSE,
        autoWidth = TRUE,
        scrollX = TRUE,
        scrollY = "155px",
        pagingType = "simple_numbers",
        dom = "tp",
        language = list(
          paginate = list(previous = "Prev", `next` = "Next")
        )
      )
    )
  })

  output$concise_guide_panel <- renderUI({
    tagList(
      div(
        class = "logic-flow",
        div(
          class = "logic-step",
          div(class = "logic-step-title", "S1-type"),
          div(class = "logic-step-loci", "N57/N62 + 135/136")
        ),
        div(class = "logic-arrow", HTML("&rarr;")),
        div(
          class = "logic-step",
          div(class = "logic-step-title", "Genotype"),
          div(class = "logic-step-loci", "N57/N62 + N1192/N1194")
        ),
        div(class = "logic-arrow", HTML("&rarr;")),
        div(
          class = "logic-step",
          div(class = "logic-step-title", "Geo-type"),
          div(class = "logic-step-loci", "N57/N62 + 135/136 + G1157 + N1192/N1194")
        ),
        div(class = "logic-arrow", HTML("&rarr;")),
        div(
          class = "logic-step",
          div(class = "logic-step-title", "Haplotype"),
          div(class = "logic-step-loci", "N57/N62 + 135/136 + G1157 + N1192/N1194 + N718/N722")
        )
      ),
      div(
        class = "logic-note",
        "Use sequence typing for FASTA input, or strain search for known accessions. Short sequences are classified only at confidently detected typing levels."
      )
    )
  })

  output$trend_region_filter <- renderUI({
    db <- database_data()
    choices <- unique(c(
      "All",
      sort(unique(db$Region[db$Region != ""]))
    ))
    selectInput("trend_region", "Region", choices = choices, selected = "All")
  })

  output$trend_country_filter <- renderUI({
    db <- database_data()
    region <- if (is.null(input$trend_region)) "All" else input$trend_region
    if (!identical(region, "All")) {
      db <- db[db$Region == region, , drop = FALSE]
    }
    countries <- sort(unique(db$Country[!is.na(db$Country) & db$Country != ""]))
    countries <- countries[!grepl("Taiwan", countries, ignore.case = TRUE)]
    choices <- c("All", countries)
    selected <- if (!is.null(input$trend_country) && input$trend_country %in% choices) input$trend_country else "All"
    selectInput("trend_country", "Country", choices = choices, selected = selected)
  })

  trend_system_col <- function(system) {
    switch(
      system,
      "S1-type" = "S1-type",
      "Genotype" = "Genotype",
      "Haplotype" = "Haplotype",
      "Geo-type" = "Geo-type",
      "Geo-type"
    )
  }

  trend_color_tree_type <- function(system) {
    switch(
      system,
      "S1-type" = "S1 Tree",
      "Genotype" = "Genotype Tree",
      "Haplotype" = "Haplotype Tree",
      "Geo-type" = "Geo Tree",
      "Geo Tree"
    )
  }

  normalize_trend_color_label <- function(x, system = NULL) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("European/American", "Euro-American", x, fixed = TRUE)
    if (!is.null(system) && identical(system, "S1-type")) {
      x[x == "Asian non-SINDEL"] <- "non-SINDEL"
      x[x == "European/American non-SINDEL"] <- "Prevalent non-SINDEL"
      x[x == "Euro-American non-SINDEL"] <- "Prevalent non-SINDEL"
    }
    x
  }

  trend_color_map <- reactive({
    path <- if (file.exists("four_tree_legend_branch_colors.csv")) {
      "four_tree_legend_branch_colors.csv"
    } else {
      file.path("www", "four_tree_legend_branch_colors.csv")
    }
    if (!file.exists(path)) return(data.frame())
    cmap <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("tree_type", "type_label", "color_hex")
    if (!all(required %in% names(cmap))) return(data.frame())
    cmap <- cmap[stats::complete.cases(cmap[, required]), required, drop = FALSE]
    cmap$type_key <- paste(cmap$tree_type, normalize_trend_color_label(cmap$type_label), sep = "||")
    cmap
  })

  trend_type_colors <- function(type_labels, system) {
    fallback <- c(
      "#2563EB", "#DC2626", "#059669", "#D97706", "#7C3AED", "#0891B2",
      "#DB2777", "#65A30D", "#111827", "#F97316", "#0F766E", "#9333EA",
      "#64748B", "#B91C1C", "#1D4ED8", "#15803D"
    )
    cols <- rep(fallback, length.out = length(type_labels))
    cmap <- trend_color_map()
    if (nrow(cmap) > 0) {
      tree_type <- trend_color_tree_type(system)
      keys <- paste(tree_type, normalize_trend_color_label(type_labels, system), sep = "||")
      hit <- match(keys, cmap$type_key)
      cols[!is.na(hit)] <- cmap$color_hex[hit[!is.na(hit)]]
    }
    cols[type_labels == "Other"] <- "#CBD5E1"
    cols
  }

  output$trend_system_filter <- renderUI({
    choices <- c("Geo-type", "Genotype", "S1-type", "Haplotype")
    active <- active_system()$system
    selected <- if (!is.null(input$trend_system) && input$trend_system %in% choices) {
      input$trend_system
    } else if (active %in% choices) {
      active
    } else {
      "Geo-type"
    }
    selectInput("trend_system", "Typing system", choices = choices, selected = selected)
  })

  output$trend_type_filter <- renderUI({
    db <- database_data()
    region <- if (is.null(input$trend_region)) "All" else input$trend_region
    country <- if (is.null(input$trend_country)) "All" else input$trend_country
    if (!identical(region, "All")) {
      db <- db[db$Region == region, , drop = FALSE]
    }
    if (!identical(country, "All")) {
      db <- db[db$Country == country, , drop = FALSE]
    }
    system <- if (is.null(input$trend_system)) active_system()$system else input$trend_system
    col <- trend_system_col(system)
    types <- sort(unique(as.character(db[[col]][vapply(as.character(db[[col]]), valid_type_label, logical(1))])))
    choices <- c("All" = "__ALL__", stats::setNames(types, types))
    active <- active_system()
    default <- if (identical(active$system, system) && !is.null(active$value) && active$value %in% types) {
      active$value
    } else {
      "__ALL__"
    }
    selected <- if (!is.null(input$trend_type) && input$trend_type %in% choices) input$trend_type else default
    selectInput("trend_type", "Specific type", choices = choices, selected = selected)
  })

  trend_period <- function(year) {
    y <- suppressWarnings(as.integer(year))
    ifelse(is.na(y), NA_character_, ifelse(y <= 2010, "<=2010", as.character(y)))
  }

  output$dynamic_trend_plot <- renderPlot({
    db <- database_data()
    db$Period <- trend_period(db$Year)
    db <- db[!is.na(db$Period), , drop = FALSE]
    region <- if (is.null(input$trend_region)) "All" else input$trend_region
    country <- if (is.null(input$trend_country)) "All" else input$trend_country
    if (!identical(region, "All")) {
      db <- db[db$Region == region, , drop = FALSE]
    }
    if (!identical(country, "All")) {
      db <- db[db$Country == country, , drop = FALSE]
    }
    validate(need(nrow(db) > 0, "No records are available for the selected region/country."))

    system <- if (is.null(input$trend_system)) active_system()$system else input$trend_system
    system_col <- trend_system_col(system)
    selected_type <- if (is.null(input$trend_type)) "__ALL__" else input$trend_type
    periods <- unique(db$Period)
    periods <- periods[order(ifelse(periods == "<=2010", 2010L, suppressWarnings(as.integer(periods))))]

    type_vec <- as.character(db[[system_col]])
    keep <- vapply(type_vec, valid_type_label, logical(1))
    db <- db[keep, , drop = FALSE]
    type_vec <- type_vec[keep]
    validate(need(nrow(db) > 0, "No valid typing records are available for this selection."))

    if (!identical(selected_type, "__ALL__") && valid_type_label(selected_type)) {
      type_vec <- ifelse(type_vec == selected_type, selected_type, "Other")
      type_levels <- c(selected_type, "Other")
      plot_title <- paste0(system_col, ": ", selected_type)
    } else {
      totals <- sort(table(type_vec), decreasing = TRUE)
      type_levels <- names(totals)
      plot_title <- paste0(system_col, " distribution")
    }
    mat <- table(factor(type_vec, levels = type_levels), factor(db$Period, levels = periods))
    mat <- as.matrix(mat)
    mat <- mat[rowSums(mat) > 0, , drop = FALSE]
    validate(need(nrow(mat) > 0, "No valid typing records are available for this selection."))
    col_total <- colSums(mat)
    pct_mat <- sweep(mat, 2, pmax(col_total, 1), "/") * 100

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(6.0, 4.4, 4.0, 1.0), mgp = c(2.2, 0.72, 0), xpd = NA, cex.axis = 1.0, cex.lab = 1.0)

    cols <- trend_type_colors(rownames(mat), system)
    bp <- barplot(
      pct_mat,
      col = cols,
      border = NA,
      names.arg = rep("", length(periods)),
      ylim = c(0, 100),
      ylab = "Percentage (%)",
      main = "",
      cex.main = 1.0,
      space = 0.35
    )
    title(main = paste0(plot_title, " (till 2025)"), line = 2.2, cex.main = 1.0)
    axis(1, at = bp, labels = periods, las = 1, cex.axis = 0.9)
    grid(nx = NA, ny = NULL, col = "#e5e7eb", lty = 1)
    cumulative <- apply(pct_mat, 2, cumsum)
    bottoms <- rbind(rep(0, ncol(pct_mat)), cumulative[-nrow(pct_mat), , drop = FALSE])
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        val <- mat[i, j]
        pct_val <- pct_mat[i, j]
        if (val > 0 && pct_val >= 10) {
          text(bp[j], bottoms[i, j] + pct_val / 2, labels = val, cex = 0.78, col = "#111827", font = 2)
        }
      }
    }
    legend(
      x = "bottom",
      inset = c(0, -0.25),
      legend = rownames(mat),
      fill = cols,
      border = NA,
      bty = "n",
      cex = 0.9,
      ncol = min(6, max(1, ceiling(nrow(mat) / 2))),
      x.intersp = 0.38,
      y.intersp = 0.85
    )
  })

  ###############################
  # Database search table and static panels
  ###############################
  output$db_table <- renderDT({
    df <- db_filtered()
    display_cols <- intersect(
      c(
        "Accession", "Sequence ID", "Strain/Isolate", "Country", "Region",
        "Collection date", "Year", "Haplotype", "Geo-type", "Genotype",
        "S1-type", "N57/N62", "135/136 motif", "N1192/N1194",
        "G1157", "N718/N722"
      ),
      names(df)
    )
    datatable(
      df[, display_cols, drop = FALSE],
      selection = "single",
      filter = "top",
      options = list(
        pageLength = 3,
        lengthChange = FALSE,
        searching = FALSE,
        info = FALSE,
        autoWidth = FALSE,
        scrollX = TRUE,
        scrollY = "95px",
        dom = "t"
      )
    )
  })

  output$db_overview_cards <- renderUI({
    db <- database_data()
    metric_card <- function(value, label) {
      div(
        class = "overview-card",
        div(class = "overview-value", value),
        div(class = "overview-label", label)
      )
    }
    div(
      class = "overview-grid",
      metric_card(format(nrow(db), big.mark = ","), "Sequences"),
      metric_card(length(unique(db$Country[db$Country != ""])), "Countries"),
      metric_card(length(unique(db$Region[db$Region != ""])), "Regions"),
      metric_card(length(unique(db$Haplotype[db$Haplotype != ""])), "Haplotypes")
    )
  })

  output$db_selected_table <- renderDT({
    row <- selected_database_row()
    validate(need(!is.null(row),
                  "Select one database row, or enter a query that returns one sequence."))
    fields <- intersect(
      c(
        "Accession", "Sequence ID", "Strain/Isolate", "Country", "Region",
        "Location", "Collection date", "Year", "Haplotype", "Geo-type",
        "Genotype", "S1-type", "N57/N62", "135/136 motif", "N1192/N1194",
        "G1157", "N718/N722", "Notes"
      ),
      names(row)
    )
    df <- data.frame(
      Field = fields,
      Value = as.character(row[1, fields, drop = TRUE]),
      check.names = FALSE
    )
    datatable(
      df,
      rownames = FALSE,
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        autoWidth = FALSE
      )
    )
  })

  output$db_group_summary_table <- renderDT({
    row <- selected_database_row()
    validate(need(!is.null(row),
                  "Select one database row, or enter a query that returns one sequence."))
    db <- database_data()
    systems <- c(
      "Haplotype" = "Haplotype",
      "Geo-type" = "Geo-type",
      "Genotype" = "Genotype",
      "S1-type" = "S1-type"
    )
    summary_df <- purrr::map_dfr(names(systems), function(label) {
      col <- systems[[label]]
      selected_type <- value_or_na(row, col)
      n <- if (selected_type == NOT_AVAILABLE) {
        0
      } else {
        sum(db[[col]] == selected_type, na.rm = TRUE)
      }
      data.frame(
        `Typing system` = label,
        `Selected group` = selected_type,
        `Database count` = n,
        `Database percent` = sprintf("%.2f%%", n / nrow(db) * 100),
        check.names = FALSE
      )
    })
    datatable(
      summary_df,
      rownames = FALSE,
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        autoWidth = FALSE,
        scrollX = TRUE,
        scrollY = "95px",
        dom = "t"
      )
    )
  })

  output$db_group_panel <- renderUI({
    row <- selected_database_row()
    if (is.null(row)) {
      return(tagList(
        tags$p("Select one strain to display its typing assignment and the corresponding static distribution panels."),
        static_distribution_cards()
      ))
    }
    tagList(
      tags$p(
        tags$strong("Selected assignment: "),
        paste0(
          "Haplotype = ", value_or_na(row, "Haplotype"),
          "; Geo-type = ", value_or_na(row, "Geo-type"),
          "; Genotype = ", value_or_na(row, "Genotype"),
          "; S1-type = ", value_or_na(row, "S1-type"), "."
        )
      ),
      static_distribution_cards()
    )
  })

  output$db_tree_panel <- renderUI({
    static_tree_panel()
  })

  output$static_info_panel <- renderUI({
    tagList(
      tags$p("Static panels are based on the curated 3042-sequence S-gene dataset and the n500 PARNAS-pruned tree dataset."),
      static_distribution_cards(),
      hr(),
      static_tree_panel()
    )
  })

  output$download_db_filtered <- downloadHandler(
    filename = function() {
      paste0("PEDV_S_database_search_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      write.csv(db_filtered(), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  output$download_db_full <- downloadHandler(
    filename = function() {
      paste0("PEDV_S_full_typing_database_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      write.csv(database_data(), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  ncbi_link_accessions <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    pattern <- "\\b[A-Z]{1,2}[0-9]{5,8}(?:\\.[0-9]+)?\\b"
    vapply(
      x,
      function(s) {
        if (!nzchar(s) || grepl("<a\\s", s, ignore.case = TRUE)) return(s)
        hits <- gregexpr(pattern, s, perl = TRUE)[[1]]
        if (hits[1] < 0) return(s)
        lens <- attr(hits, "match.length")
        starts <- as.integer(hits)
        ends <- starts + lens - 1L
        out <- character(length(starts) * 2L + 1L)
        cursor <- 1L
        k <- 1L
        for (i in seq_along(starts)) {
          if (starts[i] > cursor) {
            out[k] <- substr(s, cursor, starts[i] - 1L)
            k <- k + 1L
          }
          acc <- substr(s, starts[i], ends[i])
          out[k] <- paste0(
            "<a href=\"https://www.ncbi.nlm.nih.gov/nuccore/",
            acc,
            "\" target=\"_blank\">",
            acc,
            "</a>"
          )
          k <- k + 1L
          cursor <- ends[i] + 1L
        }
        if (cursor <= nchar(s)) out[k] <- substr(s, cursor, nchar(s))
        paste0(out[nzchar(out)], collapse = "")
      },
      character(1)
    )
  }

  read_editable_description_csv <- function(path) {
    read_with <- function(enc) {
      suppressWarnings(
        tryCatch(
          read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = enc),
          error = function(e) NULL
        )
      )
    }
    candidates <- list(
      read_with("UTF-8-BOM"),
      read_with("UTF-8"),
      read_with("GB18030"),
      suppressWarnings(tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL))
    )
    candidates <- candidates[vapply(candidates, function(x) !is.null(x) && nrow(x) > 0 && ncol(x) > 1, logical(1))]
    validate(need(length(candidates) > 0, "typing_description_editable.csv could not be read."))
    candidates[[which.max(vapply(candidates, nrow, integer(1)))]]
  }
  
  ###############################
  # 4. Load editable typing descriptions
  ###############################
  info_data <- reactive({
    csv_path <- "typing_description_editable.csv"
    if (file.exists(csv_path)) {
      df <- read_editable_description_csv(csv_path)
    } else {
      path <- if (file.exists("Information_current_typing.xlsx")) {
        "Information_current_typing.xlsx"
      } else {
        "Information.xlsx"
      }
      validate(need(file.exists(path), "typing_description_editable.csv or Information.xlsx not found"))
      src <- read_excel(path, sheet = 1) %>% as.data.frame()
      levels <- c("S1-type", "Genotype", "Geo-type", "Haplotype")
      df <- purrr::map_dfr(levels, function(level) {
        rows <- src[!is.na(src[[level]]) & nzchar(as.character(src[[level]])), , drop = FALSE]
        split_rows <- split(rows, as.character(rows[[level]]))
        purrr::imap_dfr(split_rows, function(z, type) {
          data.frame(
            `Typing system` = level,
            Type = type,
            `Polymorphic pattern` = paste(unique(as.character(z[["Polymorphic pattern"]])), collapse = " | "),
            `Selected reference strains` = paste(unique(as.character(z[["Reference strains"]])), collapse = " | "),
            Description = unique(as.character(z[["Description"]]))[1],
            `Parent type` = "",
            `Display order` = 0,
            check.names = FALSE
          )
        })
      })
    }
    if ("Reference strains" %in% names(df) && !("Selected reference strains" %in% names(df))) {
      names(df)[names(df) == "Reference strains"] <- "Selected reference strains"
    }
    for (nm in c("Typing system", "Type", "Polymorphic pattern", "Selected reference strains", "Description", "Parent type")) {
      if (!(nm %in% names(df))) df[[nm]] <- ""
      df[[nm]] <- as.character(df[[nm]])
      df[[nm]][is.na(df[[nm]])] <- ""
      df[[nm]] <- translate_display_terms(df[[nm]])
    }
    df[["Description"]] <- translate_description_terms(df[["Description"]])
    if (!("Display order" %in% names(df))) df[["Display order"]] <- seq_len(nrow(df))
    df[["Selected reference strains"]] <- ncbi_link_accessions(df[["Selected reference strains"]])
    df[order(suppressWarnings(as.numeric(df[["Display order"]])), df[["Typing system"]], df[["Type"]]), , drop = FALSE]
  })

  ###############################
  # 5. Match descriptions to the active result hierarchy
  ###############################
  info_filtered <- reactive({
    info <- info_data()
    ctx <- active_context()
    if (is.null(ctx)) {
      return(data.frame(
        Message = "Select a result row to view the matched polymorphic pattern, reference strains, and group description.",
        check.names = FALSE
      ))
    }
    if (valid_type_label(ctx$haplotype)) {
      wanted <- data.frame(
        `Typing system` = c("Genotype", "Geo-type", "Haplotype"),
        Type = c(ctx$genotype, ctx$geo, ctx$haplotype),
        Rank = c(1, 2, 3),
        OutputSystem = "Haplotype",
        OutputType = ctx$haplotype,
        check.names = FALSE
      )
    } else if (valid_type_label(ctx$geo)) {
      wanted <- data.frame(
        `Typing system` = c("Genotype", "Geo-type"),
        Type = c(ctx$genotype, ctx$geo),
        Rank = c(1, 2),
        OutputSystem = "Geo-type",
        OutputType = ctx$geo,
        check.names = FALSE
      )
    } else if (valid_type_label(ctx$genotype)) {
      wanted <- data.frame(
        `Typing system` = "Genotype",
        Type = ctx$genotype,
        Rank = 1,
        OutputSystem = "Genotype",
        OutputType = ctx$genotype,
        check.names = FALSE
      )
    } else {
      wanted <- data.frame(
        `Typing system` = "S1-type",
        Type = ctx$s1,
        Rank = 1,
        OutputSystem = "S1-type",
        OutputType = ctx$s1,
        check.names = FALSE
      )
    }
    keep <- vapply(wanted$Type, valid_type_label, logical(1))
    wanted <- wanted[keep, , drop = FALSE]
    if (nrow(wanted) == 0) {
      return(data.frame(
        Message = "No assigned typing level is available for this result.",
        check.names = FALSE
      ))
    }
    key_info <- paste(normalize_motif(info[["Typing system"]]), normalize_motif(info[["Type"]]), sep = "||")
    key_wanted <- paste(normalize_motif(wanted[["Typing system"]]), normalize_motif(wanted[["Type"]]), sep = "||")
    hit <- match(key_wanted, key_info)
    matched <- info[hit[!is.na(hit)], , drop = FALSE]
    if (nrow(matched) == 0) {
      return(data.frame(
        Message = "No editable description row matched this typing result.",
        check.names = FALSE
      ))
    }
    matched <- matched[order(wanted$Rank[!is.na(hit)]), , drop = FALSE]
    split_items <- function(x) {
      x <- gsub("\uff1b", ";", as.character(x), fixed = TRUE)
      x <- unlist(strsplit(x, "\\s*;\\s*|\\s*\\|\\s*", perl = TRUE), use.names = FALSE)
      x <- trimws(x)
      unique(x[nzchar(x)])
    }
    collapse_field <- function(x, sep = "; ") {
      vals <- unique(unlist(lapply(x, split_items), use.names = FALSE))
      vals <- vals[nzchar(vals)]
      paste(vals, collapse = sep)
    }
    collapse_description <- function(x) {
      vals <- trimws(as.character(x))
      vals <- vals[nzchar(vals)]
      vals <- unique(vals)
      paste(vals, collapse = "; ")
    }
    data.frame(
      `Typing system` = wanted$OutputSystem[1],
      Type = wanted$OutputType[1],
      `Polymorphic pattern` = collapse_field(matched[["Polymorphic pattern"]]),
      `Selected reference strains` = collapse_field(matched[["Selected reference strains"]], sep = " | "),
      Description = collapse_description(matched[["Description"]]),
      check.names = FALSE
    )
  })
  
  ###############################
  # 6锔忊儯 Render Information table
  ###############################
  output$info_table <- renderDT({
    
    df <- info_filtered()

    if ("Message" %in% names(df)) {
      return(datatable(
        df,
        rownames = FALSE,
        options = list(
          paging = FALSE,
          searching = FALSE,
          info = FALSE,
          autoWidth = FALSE,
          scrollX = TRUE,
          dom = "t"
        )
      ))
    }
    
    # 灞曠ず鐢ㄦ暟鎹細鍘绘帀鍐呴儴鍒?
    ctx <- active_context()
    s1_only <- !is.null(ctx) &&
      valid_type_label(ctx$s1) &&
      !valid_type_label(ctx$genotype) &&
      !valid_type_label(ctx$geo) &&
      !valid_type_label(ctx$haplotype)
    info_cols <- intersect(
      c("Typing system", "Type", "Polymorphic pattern", "Selected reference strains", if (!s1_only) "Description"),
      names(df)
    )
    df_display <- df[, info_cols, drop = FALSE]
    if ("Description" %in% names(df_display)) {
      df_display$Description <- gsub(
        "\\s*Current typing assignment:.*?(?=The defining polymorphic profile is|$)",
        " ",
        df_display$Description,
        perl = TRUE
      )
      df_display$Description <- gsub(
        "\\s*The defining polymorphic profile is.*$",
        "",
        df_display$Description,
        perl = TRUE
      )
      df_display$Description <- trimws(df_display$Description)
    }
    
    n <- ncol(df_display)
    
    datatable(
      df_display,
      escape = FALSE,
      rownames = FALSE,
      selection = "none",
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        autoWidth = TRUE,
        scrollY = "180px",
        dom = "t",
        columnDefs = list(
          list(targets = 0, width = "80px"),
          list(targets = 1, width = "100px"),
          list(targets = n - 2, width = "220px"),
          list(targets = n - 1, width = "220px")
        )
      )
    )
  })

  output$typing_logic_table <- renderDT({
    df <- info_data()
    display_cols <- intersect(
      c(
        "Typing system", "Type", "Polymorphic pattern", "Selected reference strains", "Description"
      ),
      names(df)
    )
    df_display <- df[, display_cols, drop = FALSE]
    datatable(
      df_display,
      escape = FALSE,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 8,
        lengthChange = FALSE,
        searching = FALSE,
        info = TRUE,
        scrollX = TRUE,
        autoWidth = FALSE
      )
    )
  })
  
  ###############################
  # 7锔忊儯 Download results
  ###############################
  output$download_res <- downloadHandler(
    filename = "PEDV_S_polymorphism_typing_results.csv",
    content = function(file) {
      write.csv(result_display(results()), file, row.names = FALSE)
    }
  )
}



###############################
# RUN APP
###############################
shinyApp(ui, server)

