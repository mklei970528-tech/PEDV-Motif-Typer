###############################
## PEDV Motif Typing – FIXED SHINYAPPS VERSION
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
library(DECIPHER)
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
  txt <- text_lines[nzchar(text_lines)]
  headers <- which(startsWith(txt, ">"))
  if(length(headers)==0) return(NULL)
  
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
  
  # Remove whitespace, unicode, '?', digits, punctuation
  s <- toupper(gsub("[^A-Z\\-]", "", seq_raw))
  s <- gsub("-", "", s)
  
  if(nchar(s)==0)
    return(list(clean_seq="", seq_type="NT"))
  
  chars <- strsplit(s,"")[[1]]
  atcg <- sum(chars %in% c("A","T","C","G"))
  aa   <- sum(chars %in% AA_ALL)
  non  <- length(chars) - atcg - aa
  
  prop_atcg <- atcg / length(chars)
  prop_non  <- non  / length(chars)
  
  # NT-like
  if(prop_atcg >= 0.90 && prop_non <= 0.05){
    clean <- gsub("[^ACGT]", "A", s)
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
  M2  = "HELQNHTATEYFV",
  M3  = "GVISSLSSSTFNSTRELP",
  M4  = "LVPGDFV"
)


###############################
# Enhanced NT 3-frame translation: fixes unicode
###############################
translate_three_frames <- function(seq_raw){
  
  seq_raw <- gsub("[^ACGT]", "A", toupper(seq_raw))
  dna <- DNAString(seq_raw)
  len <- length(dna)
  
  rf1 <- if(len>=3) subseq(dna,1,len-len%%3) else DNAString("")
  rf2 <- if(len>=4) subseq(dna,2,len-(len-1)%%3) else DNAString("")
  rf3 <- if(len>=5) subseq(dna,3,len-(len-2)%%3) else DNAString("")
  
  list(
    AA1 = translate(rf1),
    AA2 = translate(rf2),
    AA3 = translate(rf3)
  )
}


###############################
# AA local alignment
###############################
local_align_motif <- function(fullAA, motifAA, cut=30){
  
  fullA <- clean_aa_for_alignment(as.character(fullAA))
  motA  <- clean_aa_for_alignment(motifAA)
  
  aa_set <- AAStringSet(c(motif=motA, query=fullA))
  aln <- AlignSeqs(aa_set, refinements=0, iterations=0)
  
  a1 <- strsplit(as.character(aln[[1]]),"")[[1]]
  a2 <- strsplit(as.character(aln[[2]]),"")[[1]]
  
  motif_pos <- which(a1 != "-")
  extracted <- paste(a2[motif_pos][a2[motif_pos]!="-" ], collapse="")
  
  m1 <- a1[motif_pos]
  m2 <- a2[motif_pos]
  valid <- which(m2!="-")
  
  if(length(valid)==0) 
    return(list(aa="Not aligned", identity=0))
  
  ident <- sum(m1[valid]==m2[valid]) / length(valid) * 100
  ident <- round(ident,2)
  
  if(ident < cut) 
    return(list(aa="Not aligned", identity=ident))
  
  return(list(aa=extracted, identity=ident))
}



###############################
# Extract all motifs
###############################
extract_all_motifs <- function(fullAA, motif_AA){
  
  M1A <- local_align_motif(fullAA, motif_AA$M1A,50)
  M1B <- local_align_motif(fullAA, motif_AA$M1B,50)
  M1  <- if(M1A$identity >= M1B$identity) M1A else M1B
  
  list(
    M1 = M1,
    M2 = local_align_motif(fullAA, motif_AA$M2,50),
    M3 = local_align_motif(fullAA, motif_AA$M3,50),
    M4 = local_align_motif(fullAA, motif_AA$M4,50),
    M1A_ident = M1A$identity,
    M1B_ident = M1B$identity
  )
}


###############################
# Frame scoring
###############################
select_best_frame <- function(frames, motif_AA){
  
  r1 <- extract_all_motifs(frames$AA1, motif_AA)
  r2 <- extract_all_motifs(frames$AA2, motif_AA)
  r3 <- extract_all_motifs(frames$AA3, motif_AA)
  
  score <- function(r){
    ids <- c(r$M1$identity, r$M2$identity, r$M3$identity, r$M4$identity)
    good <- ids[ids>0]
    if(length(good)==0) return(0)
    mean(good)
  }
  
  sc <- c(score(r1), score(r2), score(r3))
  best <- which.max(sc)
  
  list(motifs=list(r1,r2,r3)[[best]])
}


###############################
# Glycan detection
###############################
find_nglyco <- function(aa){
  if(aa=="Not aligned") return(tibble())
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
# Analyze one sequence (FIXED)
###############################
analyze_one_sequence <- function(seq_raw, seq_id, motif_AA){
  
  pre <- preprocess_sequence(seq_raw)
  seq_clean <- pre$clean_seq
  seq_type  <- pre$seq_type
  
  if(seq_type=="NT"){
    frames <- translate_three_frames(seq_clean)
    motifs <- select_best_frame(frames, motif_AA)$motifs
  } else {
    fullAA <- AAString(clean_aa_for_alignment(seq_clean))
    motifs <- extract_all_motifs(fullAA, motif_AA)
  }
  
  failed <- c(motifs$M1$aa, motifs$M2$aa, motifs$M3$aa, motifs$M4$aa)
  
  if(any(failed=="Not aligned")){
    return(tibble(
      sequence_ID = seq_id,
      Status      = "Failed",
      Typing_Combination = "Not typed",
      M1_type = NA, M2_type = NA, M3_type = NA, M4_type = NA,
      M1_identity = motifs$M1$identity, M1_AA = motifs$M1$aa,
      M2_identity = motifs$M2$identity, M2_AA = motifs$M2$aa,
      M3_identity = motifs$M3$identity, M3_AA = motifs$M3$aa,
      M4_identity = motifs$M4$identity, M4_AA = motifs$M4$aa
    ))
  }
  
  g1 <- glyco_pattern(find_nglyco(motifs$M1$aa))
  g2 <- glyco_pattern(find_nglyco(motifs$M2$aa))
  g3 <- glyco_pattern(find_nglyco(motifs$M3$aa))
  t4 <- motif_type_M4(motifs$M4$aa)
  Type_Combination <- paste(g1, g2, g3, t4, sep=" | ")
  
  tibble(
    sequence_ID = seq_id,
    Status      = "Success",
    Motif_Type = get_motif_type(Type_Combination),
    Typing_Combination = Type_Combination,
    M1_type = g1,
    M2_type = g2,
    M3_type = g3,
    M4_type = t4,
    M1_AA = motifs$M1$aa,
    M2_AA = motifs$M2$aa,
    M3_AA = motifs$M3$aa,
    M4_AA = motifs$M4$aa
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
")
             ),
  
  
  theme = bs_theme(bootswatch = "flatly"),
  
  titlePanel("PEDV Motif Typing Platform"),
  
  sidebarLayout(
    sidebarPanel(
      textAreaInput(
        "paste_fasta", "Paste FASTA sequences:",
        width = "100%", height = "200px",
        placeholder = ">seq1\nATGC... or AA..."
      ),
      
      fileInput("seqfile", "Upload FASTA file (≤20MB)"),
      actionButton("run", "Start Analysis", class = "btn btn-primary"),
      hr(),
      downloadButton("download_res", "Download Results CSV")
    ),
    
    mainPanel(
      fluidRow(
        column(
          width = 12,
          withSpinner(
            DTOutput("res_table"),
            type = 6, color = "#2c3e50"
          )
        )
      ),
      
      hr(),
      
      ### NEW: Information table (right-bottom style) ###
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
      )
    )
  )
)



###############################
# SERVER (FINAL – FILTER MODE)
###############################
server <- function(input, output, session) {
  
  ###############################
  # 文本规范化函数（用于稳健匹配）
  ###############################
  normalize_motif <- function(x) {
    x <- as.character(x)
    x <- gsub("<[^>]+>", "", x)                 # 去 HTML
    x <- gsub("\u00A0", " ", x, fixed = TRUE)   # 去 NBSP
    x <- gsub("\\s+", " ", x)                   # 多空格压缩
    trimws(x)
  }
  
  ###############################
  # 保存“被点击的 Motif_Type”
  ###############################
  selected_motif <- reactiveVal(NULL)
  
  ###############################
  # 1️⃣ Analysis results（不动）
  ###############################
  results <- eventReactive(input$run, {
    
    validate(
      need(input$paste_fasta != "" || !is.null(input$seqfile),
           "Please paste or upload FASTA sequences.")
    )
    
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
  # 2️⃣ Results table（行点击）
  ###############################
  output$res_table <- renderDT({
    req(results())
    
    datatable(
      results(),
      selection = "single",   # ★ 行点击（最稳定）
      options = list(
        pageLength = 6,        # ★ 每页 6 行
        lengthChange = FALSE, # ★ 不允许用户改页大小
        paging = TRUE,        # ★ 开启分页
        searching = FALSE,
        info = FALSE,
        autoWidth = FALSE,
        scrollX = TRUE
      )
    )
  })
  
  ###############################
  # 3️⃣ 行点击 → 记录 Motif_Type
  ###############################
  observeEvent(input$res_table_rows_selected, {
    
    idx <- input$res_table_rows_selected
    if (length(idx) == 0) return()
    
    motif_clicked <- normalize_motif(results()$Motif_Type[idx])
    if (motif_clicked == "") return()
    
    selected_motif(motif_clicked)
    
  }, ignoreInit = TRUE)
  
  ###############################
  # 4️⃣ Load Information.xlsx
  ###############################
  info_data <- reactive({
    
    path <- "Information.xlsx"
    validate(need(file.exists(path), "Information.xlsx not found"))
    
    sheet1 <- read_excel(path, sheet = 1)
    sheet2 <- read_excel(path, sheet = 2)
    
    ref_map <- sheet2 %>%
      dplyr::mutate(
        RefID = as.character(RefID),
        Link  = as.character(Link)
      )
    
    sheet1 <- sheet1 %>%
      dplyr::mutate(
        dplyr::across(
          everything(),
          function(x) {
            
            x <- as.character(x)
            
            for (i in seq_len(nrow(ref_map))) {
              
              rid  <- ref_map$RefID[i]
              link <- ref_map$Link[i]
              
              pattern <- paste0("\\(", rid, "\\)")
              replace <- paste0(
                " (",
                "<a href=\"", link, "\" target=\"_blank\">",
                rid,
                "</a>",
                ")"
              )
              
              x <- gsub(pattern, replace, x)
            }
            
            x
          }
        )
      ) %>%
      dplyr::mutate(
        # 内部匹配用字段
        Motif_Type_plain = normalize_motif(Motif_Type)
      )
    
    sheet1
  })
  
  ###############################
  # 5️⃣ 只保留匹配行的 Information
  ###############################
  info_filtered <- reactive({
    
    info <- info_data()
    m <- selected_motif()
    
    # 尚未点击结果表：显示全部
    if (is.null(m)) {
      return(info)
    }
    
    keep <- normalize_motif(info$Motif_Type_plain) == m
    info[keep, , drop = FALSE]
  })
  
  ###############################
  # 6️⃣ Render Information table
  ###############################
  output$info_table <- renderDT({
    
    df <- info_filtered()
    
    # 展示用数据：去掉内部列
    df_display <- df %>% dplyr::select(-Motif_Type_plain)
    
    n <- ncol(df_display)
    
    datatable(
      df_display,
      escape = FALSE,
      selection = "none",
      options = list(
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        autoWidth = FALSE,
        columnDefs = list(
          list(targets = n - 1, width = "220px"),
          list(targets = n - 0, width = "260px")
        )
      )
    )
  })
  
  ###############################
  # 7️⃣ Download results
  ###############################
  output$download_res <- downloadHandler(
    filename = "PEDV_motif_typing_results.csv",
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
}



###############################
# RUN APP
###############################
shinyApp(ui, server)

