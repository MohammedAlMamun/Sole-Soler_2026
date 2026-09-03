Detect_ScriptDir <- function(default = getwd()){
  sourced_files <- unlist(lapply(sys.frames(), function(frame){
    if(is.null(frame$ofile)){
      NA_character_
    } else {
      frame$ofile
    }
  }))
  sourced_files <- sourced_files[!is.na(sourced_files)]

  if(length(sourced_files) > 0){
    dirname(normalizePath(sourced_files[length(sourced_files)], mustWork = FALSE))
  } else {
    file_args <- commandArgs(trailingOnly = FALSE)
    file_args <- file_args[startsWith(file_args, "--file=")]
    if(length(file_args) > 0){
      script_file <- sub("^--file=", "", file_args[length(file_args)])
      dirname(normalizePath(script_file, mustWork = FALSE))
    } else {
      default
    }
  }
}

WholeChromosome_Plotter_ScriptDir <- Detect_ScriptDir()

Validate_Logical <- function(x, arg_name){
  if(!is.logical(x) || length(x) != 1 || is.na(x)){
    stop(arg_name, " must be TRUE or FALSE.")
  }
}

Validate_Positive_Number <- function(x, arg_name){
  if(!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0){
    stop(arg_name, " must be a positive number.")
  }
}

Validate_Nonnegative_Number <- function(x, arg_name){
  if(!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0){
    stop(arg_name, " must be zero or a positive number.")
  }
}

Validate_Fraction <- function(x, arg_name){
  if(!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x > 1){
    stop(arg_name, " must be a number between 0 and 1.")
  }
}

Validate_Margin_Vector <- function(x, arg_name){
  if(!is.numeric(x) || length(x) != 4 || any(is.na(x)) || any(x < 0)){
    stop(arg_name, " must be a numeric vector of four non-negative values.")
  }
}

Safe_Divide <- function(numerator, denominator, replacement = 0){
  answer <- numerator / denominator
  answer[!is.finite(answer)] <- replacement
  answer
}

Safe_Log2 <- function(values, zero_replacement = 1){
  values <- suppressWarnings(as.numeric(values))
  values[!is.finite(values) | values <= 0] <- zero_replacement
  log2(values)
}

Safe_Smooth <- function(x, y, spar = NULL){
  ok <- is.finite(x) & is.finite(y)
  x_ok <- x[ok]
  y_ok <- y[ok]

  if(length(x_ok) < 4 || length(unique(x_ok)) < 4 || length(unique(y_ok)) < 2){
    return(y)
  }

  sm <- tryCatch(
    {
      if(is.null(spar)){
        stats::smooth.spline(x_ok, y_ok)
      } else {
        stats::smooth.spline(x_ok, y_ok, spar = spar)
      }
    },
    error = function(e) NULL
  )

  if(is.null(sm)){
    return(y)
  }

  stats::predict(sm, x)$y
}

Get_Threshold <- function(score_vals, Y_axis_scale = 8){
  score_vals <- score_vals[is.finite(score_vals)]
  if(length(score_vals) == 0){
    return(1)
  }

  q_vals <- stats::quantile(score_vals, probs = c(.01, .99), na.rm = TRUE)
  iqr_val <- stats::IQR(score_vals, na.rm = TRUE)
  upper <- q_vals[2] + 1.5 * iqr_val
  lower <- q_vals[1] - 1.5 * iqr_val
  trimmed <- score_vals[score_vals < upper & score_vals > lower]

  if(length(trimmed) < 2){
    thd <- max(abs(score_vals), na.rm = TRUE)
  } else {
    thd <- mean(trimmed, na.rm = TRUE) + (Y_axis_scale * stats::sd(trimmed, na.rm = TRUE))
  }

  if(!is.finite(thd) || thd <= 0){
    thd <- max(abs(score_vals), na.rm = TRUE)
  }
  if(!is.finite(thd) || thd <= 0){
    thd <- 1
  }

  thd
}

Read_Table_Fast <- function(file, header = FALSE, ...){
  if(requireNamespace("data.table", quietly = TRUE)){
    return(data.table::fread(file, header = header, data.table = FALSE, ...))
  }

  read.table(file, header = header, stringsAsFactors = FALSE,
             quote = "", comment.char = "", ...)
}

Looks_Like_Header <- function(file){
  first_line <- readLines(file, n = 1, warn = FALSE)
  if(length(first_line) == 0){
    stop("File is empty: ", file)
  }
  fields <- strsplit(trimws(first_line), "\\s+")[[1]]
  if(length(fields) < 3){
    return(FALSE)
  }

  is.na(suppressWarnings(as.numeric(fields[2]))) ||
    is.na(suppressWarnings(as.numeric(fields[3])))
}

Read_Coverage_Bed <- function(file, CoverageColumn = 5){
  file <- path.expand(file)
  if(!file.exists(file)){
    stop("Missing coverage file: ", file)
  }

  has_header <- Looks_Like_Header(file)
  cov <- Read_Table_Fast(file, header = has_header)

  if(ncol(cov) < CoverageColumn){
    stop("Coverage file must have at least ", CoverageColumn, " columns: ", file)
  }

  if(has_header){
    lower_names <- tolower(colnames(cov))
    chrom_col <- match(TRUE, lower_names %in% c("chrom", "chr", "seqnames", "seqname"))
    start_col <- match(TRUE, lower_names %in% c("chromstart", "start", "chrom_start"))
    end_col <- match(TRUE, lower_names %in% c("chromend", "end", "chrom_end"))

    if(any(is.na(c(chrom_col, start_col, end_col)))){
      chrom_col <- 1
      start_col <- 2
      end_col <- 3
    }
  } else {
    chrom_col <- 1
    start_col <- 2
    end_col <- 3
  }

  name_val <- if(ncol(cov) >= 4) as.character(cov[[4]]) else basename(file)

  out <- data.frame(
    chrom = as.character(cov[[chrom_col]]),
    chromStart = as.numeric(cov[[start_col]]),
    chromEnd = as.numeric(cov[[end_col]]),
    name = name_val,
    score = suppressWarnings(as.numeric(cov[[CoverageColumn]])),
    stringsAsFactors = FALSE
  )
  out$score[!is.finite(out$score)] <- 0
  out
}

Find_Normal_Counts_File <- function(SampleDir, SampleName = NULL, CountsFile = NULL){
  if(!is.null(CountsFile)){
    CountsFile <- path.expand(CountsFile)
    if(file.exists(CountsFile)){
      return(CountsFile)
    }
    stop("CountsFile was supplied but not found: ", CountsFile)
  }

  SampleDir <- normalizePath(path.expand(SampleDir), mustWork = TRUE)
  SampleName <- Get_SampleName(SampleDir, SampleName)

  candidates <- unique(c(
    file.path(SampleDir, paste0(SampleName, "_normal_counts.bed")),
    file.path(SampleDir, "normal_counts.bed"),
    Sys.glob(file.path(SampleDir, "*_normal_counts.bed"))
  ))

  candidates <- candidates[file.exists(candidates)]
  if(length(candidates) == 0){
    return(NA_character_)
  }

  candidates[1]
}

Read_Normal_Counts_Bed <- function(file, SampleName = NULL){
  file <- path.expand(file)
  if(!file.exists(file)){
    stop("Missing normal counts file: ", file)
  }

  has_header <- Looks_Like_Header(file)
  counts <- Read_Table_Fast(file, header = has_header)
  if(ncol(counts) < 5){
    stop("Normal counts file must have at least five columns: ", file)
  }

  if(has_header){
    lower_names <- tolower(colnames(counts))
    chrom_col <- match(TRUE, lower_names %in% c("chrom", "chr", "seqnames", "seqname"))
    start_col <- match(TRUE, lower_names %in% c("chromstart", "start", "chrom_start"))
    end_col <- match(TRUE, lower_names %in% c("chromend", "end", "chrom_end"))
    input_col <- match(TRUE, lower_names %in% c("input", "input.score",
                                                "input_score", "true.input",
                                                "norm.input", "normal.input"))
    chip_col <- match(TRUE, lower_names %in% c("chip", "chip.score",
                                               "chip_score", "ip", "ip.score",
                                               "ip_score", "chipseq"))

    if(any(is.na(c(chrom_col, start_col, end_col, input_col, chip_col)))){
      chrom_col <- 1
      start_col <- 2
      end_col <- 3
      input_col <- 4
      chip_col <- 5
    }
  } else {
    chrom_col <- 1
    start_col <- 2
    end_col <- 3
    input_col <- 4
    chip_col <- 5
  }

  if(is.null(SampleName)){
    SampleName <- tools::file_path_sans_ext(basename(file))
  }

  out <- data.frame(
    chrom = as.character(counts[[chrom_col]]),
    chromStart = as.numeric(counts[[start_col]]),
    chromEnd = as.numeric(counts[[end_col]]),
    name = paste0(SampleName, "_normal_counts"),
    ChIP = suppressWarnings(as.numeric(counts[[chip_col]])),
    Input = suppressWarnings(as.numeric(counts[[input_col]])),
    ChIP_noise = NA_real_,
    Input_noise = NA_real_,
    stringsAsFactors = FALSE
  )

  out$ChIP[!is.finite(out$ChIP)] <- 0
  out$Input[!is.finite(out$Input)] <- 0
  out
}

Resolve_Ratio_Profile_File <- function(SampleDir, SampleName = NULL,
                                       RatioFile = NULL){
  if(is.null(RatioFile)){
    return(NA_character_)
  }

  SampleDir <- normalizePath(path.expand(SampleDir), mustWork = TRUE)
  SampleName <- Get_SampleName(SampleDir, SampleName)

  if(length(RatioFile) == 1 && identical(tolower(as.character(RatioFile)), "auto")){
    candidates <- unique(c(
      file.path(SampleDir, paste0(SampleName, "_Ratio_pro.bed")),
      file.path(SampleDir, paste0(SampleName, "_Ratio.bed")),
      file.path(SampleDir, "Ratio_pro.bed"),
      file.path(SampleDir, "Ratio.bed"),
      Sys.glob(file.path(SampleDir, "*_Ratio_pro.bed")),
      Sys.glob(file.path(SampleDir, "*_Ratio.bed"))
    ))
    candidates <- candidates[file.exists(candidates)]
    if(length(candidates) == 0){
      return(NA_character_)
    }
    return(candidates[1])
  }

  RatioFile <- path.expand(as.character(RatioFile)[1])
  if(file.exists(RatioFile)){
    return(RatioFile)
  }

  stop("RatioFile was supplied but not found: ", RatioFile)
}

Read_Ratio_Profile_Bed <- function(file, SampleName = NULL){
  file <- path.expand(file)
  if(!file.exists(file)){
    stop("Missing ratio profile file: ", file)
  }

  has_header <- Looks_Like_Header(file)
  ratio <- Read_Table_Fast(file, header = has_header)
  if(ncol(ratio) < 5){
    stop("Ratio profile file must have at least five columns: ", file)
  }

  if(has_header){
    lower_names <- tolower(colnames(ratio))
    chrom_col <- match(TRUE, lower_names %in% c("chrom", "chr", "seqnames", "seqname"))
    start_col <- match(TRUE, lower_names %in% c("chromstart", "start", "chrom_start"))
    end_col <- match(TRUE, lower_names %in% c("chromend", "end", "chrom_end"))
    name_col <- match(TRUE, lower_names %in% c("name", "sample", "sample_name"))
    ip_col <- match(TRUE, lower_names %in% c("ip.score", "ip_score", "ip",
                                             "chip.score", "chip_score", "chip"))
    true_input_col <- match(TRUE, lower_names %in% c("true.input", "true_input",
                                                     "input", "input.score",
                                                     "input_score"))
    norm_input_col <- match(TRUE, lower_names %in% c("norm.input", "norm_input",
                                                     "corr.input", "corr_input"))
    log2_col <- match(TRUE, lower_names %in% c("log2.ratio", "log2_ratio",
                                               "log2", "ratio"))
    pvalue_col <- match(TRUE, lower_names %in% c("pvalue", "p.value",
                                                 "p_value", "pval"))

    if(any(is.na(c(chrom_col, start_col, end_col)))){
      chrom_col <- 1
      start_col <- 2
      end_col <- 3
    }
    if(is.na(name_col) && ncol(ratio) >= 4){
      name_col <- 4
    }
    if(is.na(ip_col) && ncol(ratio) >= 6){
      ip_col <- 6
    }
    if(is.na(true_input_col) && ncol(ratio) >= 7){
      true_input_col <- 7
    }
    if(is.na(norm_input_col) && ncol(ratio) >= 8){
      norm_input_col <- 8
    }
    if(is.na(log2_col) && ncol(ratio) >= 9){
      log2_col <- 9
    }
    if(is.na(pvalue_col) && ncol(ratio) >= 10){
      pvalue_col <- 10
    }
  } else {
    chrom_col <- 1
    start_col <- 2
    end_col <- 3
    name_col <- if(ncol(ratio) >= 4) 4 else NA_integer_
    ip_col <- if(ncol(ratio) >= 6) 6 else 5
    true_input_col <- if(ncol(ratio) >= 7) 7 else NA_integer_
    norm_input_col <- if(ncol(ratio) >= 8) 8 else NA_integer_
    log2_col <- if(ncol(ratio) >= 9) 9 else NA_integer_
    pvalue_col <- if(ncol(ratio) >= 10) 10 else NA_integer_
  }

  profile_name <- if(is.na(name_col)){
    tools::file_path_sans_ext(basename(file))
  } else {
    as.character(ratio[[name_col]])
  }

  if(is.null(SampleName)){
    names_seen <- unique(profile_name[!is.na(profile_name) & nzchar(profile_name)])
    if(length(names_seen) > 0){
      SampleName <- names_seen[1]
    } else {
      SampleName <- tools::file_path_sans_ext(basename(file))
    }
  }

  ip_score <- suppressWarnings(as.numeric(ratio[[ip_col]]))
  true_input <- if(is.na(true_input_col)){
    rep(NA_real_, nrow(ratio))
  } else {
    suppressWarnings(as.numeric(ratio[[true_input_col]]))
  }
  norm_input <- if(is.na(norm_input_col)){
    true_input
  } else {
    suppressWarnings(as.numeric(ratio[[norm_input_col]]))
  }

  ratio_log2 <- if(is.na(log2_col)){
    Safe_Log2(Safe_Divide(ip_score, norm_input))
  } else {
    suppressWarnings(as.numeric(ratio[[log2_col]]))
  }
  ratio_pvalue <- if(is.na(pvalue_col)){
    rep(NA_real_, nrow(ratio))
  } else {
    suppressWarnings(as.numeric(ratio[[pvalue_col]]))
  }

  out <- data.frame(
    chrom = as.character(ratio[[chrom_col]]),
    chromStart = as.numeric(ratio[[start_col]]),
    chromEnd = as.numeric(ratio[[end_col]]),
    name = as.character(profile_name),
    ChIP = ip_score,
    Input = norm_input,
    ChIP_noise = NA_real_,
    Input_noise = NA_real_,
    ratio_log2 = ratio_log2,
    ratio_pvalue = ratio_pvalue,
    true_input = true_input,
    norm_input = norm_input,
    stringsAsFactors = FALSE
  )

  out$ChIP[!is.finite(out$ChIP)] <- 0
  out$Input[!is.finite(out$Input)] <- 0
  out$ratio_log2[!is.finite(out$ratio_log2)] <- 0
  attr(out, "SampleName") <- SampleName
  out
}

Check_Matching_Coverage <- function(reference, query, query_name){
  same_rows <- nrow(reference) == nrow(query)
  same_coords <- same_rows &&
    identical(reference$chrom, query$chrom) &&
    identical(reference$chromStart, query$chromStart) &&
    identical(reference$chromEnd, query$chromEnd)

  if(!same_coords){
    stop("Coverage coordinates do not match for ", query_name,
         ". Check chromosome/start/end order before strand collapsing.")
  }

  invisible(TRUE)
}

Get_SampleName <- function(SampleDir, SampleName = NULL){
  if(!is.null(SampleName)){
    return(as.character(SampleName))
  }

  basename(normalizePath(path.expand(SampleDir), mustWork = FALSE))
}

Resolve_Output_Directory <- function(OutputDir = NULL, SampleDir){
  SampleDir <- normalizePath(path.expand(SampleDir), mustWork = TRUE)

  use_sample_dir <- is.null(OutputDir) ||
    length(OutputDir) == 0 ||
    (length(OutputDir) == 1 &&
       (is.na(OutputDir) ||
          OutputDir == "" ||
          tolower(OutputDir) %in% c("sampledir", "sample_dir", "sample")))

  if(use_sample_dir){
    OutputDir <- SampleDir
  } else if(length(OutputDir) != 1){
    stop("OutputDir must be NULL, 'SampleDir', or one output folder path.")
  }

  OutputDir <- normalizePath(path.expand(OutputDir), mustWork = FALSE)
  if(!dir.exists(OutputDir)){
    created <- dir.create(OutputDir, recursive = TRUE, showWarnings = FALSE)
    if(!created && !dir.exists(OutputDir)){
      stop("Could not create OutputDir: ", OutputDir)
    }
  }

  OutputDir
}

Resolve_Command <- function(command = NULL, command_name, candidates = character()){
  candidates <- unique(c(command, candidates, command_name))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  for(candidate in candidates){
    candidate <- path.expand(candidate)
    if(grepl(.Platform$file.sep, candidate, fixed = TRUE) && file.exists(candidate)){
      return(candidate)
    }

    found <- Sys.which(candidate)
    if(nzchar(found)){
      return(unname(found))
    }
  }

  stop(command_name, " was not found. Install it, add it to PATH, or pass ",
       command_name, "Command with the full executable path.", call. = FALSE)
}

Build_Coverage_File_Paths <- function(SampleDir, SampleName = NULL, CoverageDir = NULL){
  SampleDir <- normalizePath(path.expand(SampleDir), mustWork = TRUE)
  SampleName <- Get_SampleName(SampleDir, SampleName)

  if(is.null(CoverageDir)){
    CoverageDir <- file.path(SampleDir, "Coverage")
  }
  CoverageDir <- normalizePath(path.expand(CoverageDir), mustWork = FALSE)

  list(
    SampleDir = SampleDir,
    SampleName = SampleName,
    CoverageDir = CoverageDir,
    ChIPWatson = file.path(CoverageDir, paste0(SampleName, "_ChIP_watson.bed")),
    ChIPCrick = file.path(CoverageDir, paste0(SampleName, "_ChIP_crick.bed")),
    InputWatson = file.path(CoverageDir, paste0(SampleName, "_Input_watson.bed")),
    InputCrick = file.path(CoverageDir, paste0(SampleName, "_Input_crick.bed"))
  )
}

Read_StrandCollapsed_Coverage <- function(SampleDir, SampleName = NULL,
                                          CoverageDir = NULL, CoverageColumn = 5){
  paths <- Build_Coverage_File_Paths(SampleDir = SampleDir, SampleName = SampleName,
                                     CoverageDir = CoverageDir)

  chip_watson <- Read_Coverage_Bed(paths$ChIPWatson, CoverageColumn = CoverageColumn)
  chip_crick <- Read_Coverage_Bed(paths$ChIPCrick, CoverageColumn = CoverageColumn)
  input_watson <- Read_Coverage_Bed(paths$InputWatson, CoverageColumn = CoverageColumn)
  input_crick <- Read_Coverage_Bed(paths$InputCrick, CoverageColumn = CoverageColumn)

  Check_Matching_Coverage(chip_watson, chip_crick, "ChIP crick")
  Check_Matching_Coverage(chip_watson, input_watson, "Input watson")
  Check_Matching_Coverage(chip_watson, input_crick, "Input crick")

  out <- data.frame(
    chrom = chip_watson$chrom,
    chromStart = chip_watson$chromStart,
    chromEnd = chip_watson$chromEnd,
    name = paste0(paths$SampleName, "_strandcollapsed"),
    ChIP = chip_watson$score + chip_crick$score,
    Input = input_watson$score + input_crick$score,
    ChIP_noise = NA_real_,
    Input_noise = NA_real_,
    stringsAsFactors = FALSE
  )

  list(coverage = out, paths = paths)
}

Get_SacCer3_Chromosome_Table <- function(){
  fallback <- data.frame(
    chrom = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII",
              "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII",
              "chrXIV", "chrXV", "chrXVI", "chrM"),
    length = c(230218, 813184, 316620, 1531933, 576874, 270161, 1090940,
               562643, 439888, 745751, 666816, 1078177, 924431, 784333,
               1091291, 948066, 85779),
    stringsAsFactors = FALSE
  )

  if(requireNamespace("BSgenome.Scerevisiae.UCSC.sacCer3", quietly = TRUE) &&
     requireNamespace("GenomeInfoDb", quietly = TRUE)){
    suppressPackageStartupMessages(
      library("BSgenome.Scerevisiae.UCSC.sacCer3", character.only = TRUE)
    )
    chrom_names <- as.character(GenomeInfoDb::seqnames(Scerevisiae))
    chrom_lengths <- as.numeric(GenomeInfoDb::seqlengths(Scerevisiae))
    return(data.frame(chrom = chrom_names, length = chrom_lengths,
                      stringsAsFactors = FALSE))
  }

  fallback
}

Resolve_Chromosome <- function(Chromosome, ChromosomeInfo = Get_SacCer3_Chromosome_Table()){
  if(is.numeric(Chromosome)){
    k <- as.integer(Chromosome[1])
    if(k < 1 || k > nrow(ChromosomeInfo)){
      stop("Chromosome number is outside the available sacCer3 chromosomes.")
    }
    return(list(k = k, chrom = ChromosomeInfo$chrom[k],
                length = ChromosomeInfo$length[k]))
  }

  chr <- as.character(Chromosome[1])
  if(grepl("^[0-9]+$", chr)){
    return(Resolve_Chromosome(as.numeric(chr), ChromosomeInfo))
  }

  candidates <- unique(c(chr, paste0("chr", chr),
                         gsub("^chr", "", chr, ignore.case = TRUE)))
  idx <- match(candidates, ChromosomeInfo$chrom)
  idx <- idx[!is.na(idx)]
  if(length(idx) == 0){
    stop("Could not resolve chromosome '", chr,
         "'. Use a number like 3 or a sacCer3 name like chrIII.")
  }

  list(k = idx[1], chrom = ChromosomeInfo$chrom[idx[1]],
       length = ChromosomeInfo$length[idx[1]])
}

Build_WholeChromosome_Windows <- function(Chromosome, WindowSizeKb = 100,
                                          ChromosomeInfo = Get_SacCer3_Chromosome_Table()){
  chr_info <- Resolve_Chromosome(Chromosome, ChromosomeInfo)
  window_size <- WindowSizeKb * 1000
  starts <- seq(0, chr_info$length - 1, by = window_size)
  ends <- pmin(starts + window_size, chr_info$length)

  data.frame(k = chr_info$k, chrom = chr_info$chrom,
             S = starts, E = ends, stringsAsFactors = FALSE)
}

Resolve_Plot_Chromosomes <- function(Chromosomes = "all",
                                     ChromosomeInfo = Get_SacCer3_Chromosome_Table()){
  if(is.null(Chromosomes) || length(Chromosomes) == 0){
    Chromosomes <- "all"
  }

  if(length(Chromosomes) == 1){
    request <- tolower(as.character(Chromosomes))
    if(request %in% c("all", "genome", "whole_genome", "whole genome")){
      out <- ChromosomeInfo[seq_len(min(16, nrow(ChromosomeInfo))),
                            c("chrom", "length"), drop = FALSE]
      out$k <- seq_len(nrow(out))
      return(out[, c("k", "chrom", "length")])
    }
  }

  resolved <- lapply(Chromosomes, Resolve_Chromosome,
                     ChromosomeInfo = ChromosomeInfo)
  out <- data.frame(
    k = vapply(resolved, function(x) x$k, integer(1)),
    chrom = vapply(resolved, function(x) x$chrom, character(1)),
    length = vapply(resolved, function(x) x$length, numeric(1)),
    stringsAsFactors = FALSE
  )
  out[!duplicated(out$chrom), , drop = FALSE]
}

Find_Support_File <- function(File = NULL, names, AnnotationDir = NULL){
  if(!is.null(File)){
    File <- path.expand(File)
    if(file.exists(File)){
      return(File)
    }
    stop("Annotation file was supplied but not found: ", File)
  }

  candidate_dirs <- unique(c(
    AnnotationDir,
    file.path(WholeChromosome_Plotter_ScriptDir, "WholeChromosome_Plotter_0626_support"),
    file.path(WholeChromosome_Plotter_ScriptDir, "Local_Plotter_0626_support"),
    WholeChromosome_Plotter_ScriptDir,
    file.path(getwd(), "WholeChromosome_Plotter_0626_support"),
    file.path(getwd(), "Local_Plotter_0626_support"),
    getwd()
  ))
  candidate_dirs <- candidate_dirs[!is.na(candidate_dirs)]

  for(candidate_dir in candidate_dirs){
    for(nm in names){
      candidate <- file.path(path.expand(candidate_dir), nm)
      if(file.exists(candidate)){
        return(candidate)
      }
    }
  }

  NA_character_
}

Find_ngsAnalyser_Origin_File <- function(ARSFile = NULL, AnnotationDir = NULL){
  if(!is.null(ARSFile)){
    ARSFile <- path.expand(ARSFile)
    if(file.exists(ARSFile)){
      return(ARSFile)
    }
    stop("ARSFile was supplied but not found: ", ARSFile)
  }

  candidate_dirs <- unique(c(
    AnnotationDir,
    Sys.getenv("NGSANALYSER_APP_DIR"),
    "/Volumes/rbmData/ngsAnalyser1.1.4/dist/mac/ngsAnalyser.app/Contents/Resources/app",
    "/Volumes/rbmData/ngsAnalyser1.1.4",
    "/Volumes/rbmData/ngsAnalyser/Application/ngsAnalyser-darwin-x64/ngsAnalyser.app/Contents/Resources/app",
    "/Volumes/rbmData/ngsAnalyser",
    "/Applications/ngsAnalyser.app/Contents/Resources/app",
    "~/Desktop/ngsAnalyser1.1.4"
  ))
  candidate_dirs <- candidate_dirs[!is.na(candidate_dirs) & nzchar(candidate_dirs)]

  for(candidate_dir in candidate_dirs){
    candidate <- file.path(path.expand(candidate_dir), "OriginList_Full.bed")
    if(file.exists(candidate)){
      return(candidate)
    }
  }

  NA_character_
}

Normalise_Feature_Strands <- function(features){
  if(is.null(features) || nrow(features) == 0){
    return(features)
  }

  features$strand[features$strand %in% c("1", "plus", "Plus", "Watson", "W")] <- "+"
  features$strand[features$strand %in% c("-1", "minus", "Minus", "Crick", "C")] <- "-"
  features$strand[!(features$strand %in% c("+", "-", "."))] <- "."
  features
}

Read_Feature_File <- function(FeatureFile, default_type){
  empty <- data.frame(chrom = character(), chromStart = numeric(), chromEnd = numeric(),
                      name = character(), score = numeric(), strand = character(),
                      type = character(), stat = character(), stringsAsFactors = FALSE)
  if(is.na(FeatureFile) || !file.exists(FeatureFile)){
    warning("Missing ", default_type, " annotation file; this track will be skipped.")
    return(empty)
  }

  has_header <- Looks_Like_Header(FeatureFile)
  features <- Read_Table_Fast(FeatureFile, header = has_header)
  if(ncol(features) < 3){
    warning("Annotation file has fewer than 3 columns and will be skipped: ", FeatureFile)
    return(empty)
  }

  if(has_header){
    lower_names <- tolower(colnames(features))
    chrom_col <- match(TRUE, lower_names %in% c("chrom", "chr", "seqnames", "seqname"))
    start_col <- match(TRUE, lower_names %in% c("chromstart", "start", "chrom_start"))
    end_col <- match(TRUE, lower_names %in% c("chromend", "end", "chrom_end"))
    name_col <- match(TRUE, lower_names %in% c("name", "orf", "orfiname", "gene", "feature"))
    score_col <- match(TRUE, lower_names %in% c("score"))
    strand_col <- match(TRUE, lower_names %in% c("strand"))
    type_col <- match(TRUE, lower_names %in% c("type", "featuretype"))
    stat_col <- match(TRUE, lower_names %in% c("stat", "status", "timing"))

    if(any(is.na(c(chrom_col, start_col, end_col)))){
      chrom_col <- 1
      start_col <- 2
      end_col <- 3
    }
  } else {
    chrom_col <- 1
    start_col <- 2
    end_col <- 3
    name_col <- if(ncol(features) >= 4) 4 else NA_integer_
    score_col <- if(ncol(features) >= 5) 5 else NA_integer_
    strand_col <- if(ncol(features) >= 6) 6 else NA_integer_
    type_col <- if(ncol(features) >= 7) 7 else NA_integer_
    stat_col <- NA_integer_
  }

  feature_name <- if(is.na(name_col)){
    paste0(default_type, "_", seq_len(nrow(features)))
  } else {
    features[[name_col]]
  }
  feature_score <- if(is.na(score_col)) 0 else features[[score_col]]
  feature_strand <- if(is.na(strand_col)) "." else features[[strand_col]]
  feature_type <- if(is.na(type_col)) default_type else features[[type_col]]
  feature_stat <- if(is.na(stat_col)) "" else features[[stat_col]]

  out <- data.frame(
    chrom = as.character(features[[chrom_col]]),
    chromStart = as.numeric(features[[start_col]]),
    chromEnd = as.numeric(features[[end_col]]),
    name = as.character(feature_name),
    score = suppressWarnings(as.numeric(feature_score)),
    strand = as.character(feature_strand),
    type = as.character(feature_type),
    stat = as.character(feature_stat),
    stringsAsFactors = FALSE
  )
  out$score[is.na(out$score)] <- 0
  out <- out[is.finite(out$chromStart) & is.finite(out$chromEnd), ]
  Normalise_Feature_Strands(out)
}

Get_WholeChromosome_Annotations <- function(AnnotationDir = NULL, ORFFile = NULL,
                                            ARSFile = NULL, TERFile = NULL,
                                            TyFile = NULL, tRNAFile = NULL,
                                            CentromereFile = NULL){
  list(
    ORFs = Read_Feature_File(Find_Support_File(ORFFile,
                                               c("sacCer3_S288C_ORFs.bed",
                                                 "sacCer3_ORFs.bed", "ORFs.bed"),
                                               AnnotationDir), "ORF"),
    ARS = Read_Feature_File(Find_ngsAnalyser_Origin_File(ARSFile,
                                                         AnnotationDir), "ARS"),
    TER = Read_Feature_File(Find_Support_File(TERFile,
                                              c("sacCer3_TER.bed", "TER.bed",
                                                "TER_List.bed"),
                                              AnnotationDir), "TER"),
    Ty = Read_Feature_File(Find_Support_File(TyFile,
                                             c("sacCer3_TyElements.bed",
                                               "TyElements.bed",
                                               "TyElement_list.bed"),
                                             AnnotationDir), "Ty"),
    tRNA = Read_Feature_File(Find_Support_File(tRNAFile,
                                               c("sacCer3_tRNAs.bed",
                                                 "tRNAs.bed", "tRNA_list.bed"),
                                               AnnotationDir), "tRNA"),
    Centromeres = Read_Feature_File(Find_Support_File(CentromereFile,
                                                      c("sacCer3_centromeres.bed",
                                                        "centromeres.bed",
                                                        "CEN.bed"),
                                                      AnnotationDir), "centromere")
  )
}

Get_Window_Features <- function(features, chrom, S, E){
  if(is.null(features) || nrow(features) == 0){
    return(features)
  }

  features[features$chrom == chrom & features$chromEnd >= S & features$chromStart <= E, ]
}

Read_Genomewide_Peaks <- function(PeakFile){
  empty <- data.frame(chrom = character(), peakStart = numeric(),
                      peakEnd = numeric(), peakSummit = numeric(),
                      stringsAsFactors = FALSE)
  if(is.null(PeakFile) || is.na(PeakFile) || !file.exists(PeakFile)){
    return(empty)
  }

  has_header <- Looks_Like_Header(PeakFile)
  peaks <- Read_Table_Fast(PeakFile, header = has_header)
  if(nrow(peaks) == 0 || ncol(peaks) < 3){
    return(empty)
  }

  if(has_header){
    lower_names <- tolower(colnames(peaks))
    chrom_col <- match(TRUE, lower_names %in% c("chrom", "chr"))
    start_col <- match(TRUE, lower_names %in% c("peakstart", "start"))
    end_col <- match(TRUE, lower_names %in% c("peakend", "end"))
    summit_col <- match(TRUE, lower_names %in% c("peaksummit", "abs_summit",
                                                 "summit", "peak_summit"))

    if(any(is.na(c(chrom_col, start_col, end_col)))){
      chrom_col <- 1
      start_col <- 2
      end_col <- 3
    }
  } else {
    chrom_col <- 1
    start_col <- 2
    end_col <- 3
    summit_col <- if(ncol(peaks) >= 5) 5 else NA_integer_
  }

  out <- data.frame(
    chrom = as.character(peaks[[chrom_col]]),
    peakStart = as.numeric(peaks[[start_col]]),
    peakEnd = as.numeric(peaks[[end_col]]),
    peakSummit = if(is.na(summit_col)){
      (as.numeric(peaks[[start_col]]) + as.numeric(peaks[[end_col]])) / 2
    } else {
      suppressWarnings(as.numeric(peaks[[summit_col]]))
    },
    stringsAsFactors = FALSE
  )
  out$peakSummit[!is.finite(out$peakSummit)] <-
    (out$peakStart[!is.finite(out$peakSummit)] +
       out$peakEnd[!is.finite(out$peakSummit)]) / 2
  out <- out[is.finite(out$peakStart) & is.finite(out$peakEnd), ]
  out <- out[out$peakEnd > out$peakStart, ]
  out[order(out$chrom, out$peakStart, out$peakEnd), ]
}

Command_Failed <- function(status){
  if(is.null(status)){
    return(FALSE)
  }
  if(length(status) == 0){
    return(FALSE)
  }
  any(as.integer(status) != 0L)
}

Read_Command_Log <- function(file, n = 30){
  if(is.null(file) || !file.exists(file)){
    return(character())
  }
  lines <- readLines(file, warn = FALSE)
  if(length(lines) == 0){
    return(character())
  }
  utils::tail(lines, n)
}

Stop_With_Log <- function(message_text, log_file = NULL){
  log_lines <- Read_Command_Log(log_file)
  if(length(log_lines) > 0){
    stop(message_text, "\nLast command log lines:\n",
         paste(log_lines, collapse = "\n"), call. = FALSE)
  }
  stop(message_text, call. = FALSE)
}

Process_ChIP_Peaks_From_Bam <- function(SampleDir, SampleName = NULL, PeakDir = NULL,
                                        ChIPBam = NULL, InputBam = NULL,
                                        BedtoolsCommand = NULL,
                                        Macs2Command = "macs2",
                                        GenomeSize = 12157105,
                                        PeakPValue = "10e-6"){
  SampleDir <- normalizePath(path.expand(SampleDir), mustWork = TRUE)
  SampleName <- Get_SampleName(SampleDir, SampleName)

  if(is.null(PeakDir)){
    PeakDir <- file.path(SampleDir, "ChIP_Peaks")
  }
  PeakDir <- path.expand(PeakDir)
  if(!dir.exists(PeakDir)){
    dir.create(PeakDir, recursive = TRUE, showWarnings = FALSE)
  }

  if(is.null(ChIPBam)){
    ChIPBam <- file.path(SampleDir, "Bam", paste0(SampleName, "_ChIP.bam"))
  }
  if(is.null(InputBam)){
    InputBam <- file.path(SampleDir, "Bam", paste0(SampleName, "_Input.bam"))
  }
  ChIPBam <- path.expand(ChIPBam)
  InputBam <- path.expand(InputBam)

  if(!file.exists(ChIPBam)){
    stop("Cannot process ChIP peaks; missing ChIP BAM: ", ChIPBam)
  }
  if(!file.exists(InputBam)){
    stop("Cannot process ChIP peaks; missing Input BAM: ", InputBam)
  }

  app_bedtools <- "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools"
  BedtoolsCommand <- Resolve_Command(
    command = BedtoolsCommand,
    command_name = "Bedtools",
    candidates = c(app_bedtools, "bedtools")
  )
  Macs2Command <- Resolve_Command(
    command = Macs2Command,
    command_name = "Macs2",
    candidates = c("/Library/Frameworks/Python.framework/Versions/3.7/bin/macs2",
                   "/opt/homebrew/bin/macs2",
                   "/usr/local/bin/macs2",
                   "macs2")
  )

  ip_bed <- tempfile(fileext = ".bed")
  input_bed <- tempfile(fileext = ".bed")
  macs_out <- tempfile("macs2_chip_peaks_")
  dir.create(macs_out, showWarnings = FALSE)
  on.exit(unlink(c(ip_bed, input_bed, macs_out), recursive = TRUE), add = TRUE)

  chip_bedtools_log <- file.path(PeakDir, paste0(SampleName, "_bedtools_ChIP.log"))
  input_bedtools_log <- file.path(PeakDir, paste0(SampleName, "_bedtools_Input.log"))
  macs2_log <- file.path(PeakDir, paste0(SampleName, "_MACS2_callpeak.log"))

  message("Converting BAM to BED for MACS2 peak calling.")
  ip_status <- system2(BedtoolsCommand, args = c("bamtobed", "-i", ChIPBam),
                       stdout = ip_bed, stderr = chip_bedtools_log)
  input_status <- system2(BedtoolsCommand, args = c("bamtobed", "-i", InputBam),
                          stdout = input_bed, stderr = input_bedtools_log)
  if(Command_Failed(ip_status)){
    Stop_With_Log("bedtools failed for ChIP BAM.", chip_bedtools_log)
  }
  if(Command_Failed(input_status)){
    Stop_With_Log("bedtools failed for Input BAM.", input_bedtools_log)
  }
  if(!file.exists(ip_bed) || file.info(ip_bed)$size == 0){
    Stop_With_Log("bedtools produced an empty ChIP BED file.", chip_bedtools_log)
  }
  if(!file.exists(input_bed) || file.info(input_bed)$size == 0){
    Stop_With_Log("bedtools produced an empty Input BED file.", input_bedtools_log)
  }

  message("Calling ChIP peaks with MACS2.")
  macs_args <- c("callpeak", "-t", ip_bed, "-c", input_bed,
                 "-f", "BED", "-g", as.character(GenomeSize),
                 "-p", as.character(PeakPValue), "--nomodel",
                 "-n", "Peak", "--outdir", macs_out)
  macs_status <- system2(Macs2Command, args = macs_args,
                         stdout = macs2_log, stderr = macs2_log)
  if(Command_Failed(macs_status)){
    Stop_With_Log("MACS2 peak calling failed. Check MACS2 log: ", macs2_log)
  }

  macs_peak_file <- file.path(macs_out, "Peak_peaks.xls")
  if(!file.exists(macs_peak_file)){
    Stop_With_Log("MACS2 completed but did not produce Peak_peaks.xls.", macs2_log)
  }

  macs_peaks <- read.delim(macs_peak_file, comment.char = "#",
                           stringsAsFactors = FALSE, check.names = FALSE)
  if(nrow(macs_peaks) == 0){
    warning("MACS2 produced no peaks.")
    genomewide <- data.frame(chrom = character(), peakStart = numeric(),
                             peakEnd = numeric(), peakLength = numeric(),
                             peakSummit = numeric())
  } else {
    genomewide <- data.frame(
      chrom = macs_peaks$chr,
      peakStart = macs_peaks$start,
      peakEnd = macs_peaks$end,
      peakLength = macs_peaks$length,
      peakSummit = macs_peaks$abs_summit,
      stringsAsFactors = FALSE
    )
  }

  out_file <- file.path(PeakDir, paste0(SampleName, "_Genomewide_Peaks.bed"))
  write.table(genomewide, file = out_file, quote = FALSE,
              row.names = FALSE, sep = "\t")
  out_file
}

Get_Or_Create_Genomewide_Peaks <- function(SampleDir, SampleName = NULL,
                                           PeakDir = NULL, PeakFile = NULL,
                                           ProcessPeaks = FALSE,
                                           ChIPBam = NULL, InputBam = NULL,
                                           BedtoolsCommand = NULL,
                                           Macs2Command = "macs2",
                                           GenomeSize = 12157105,
                                           PeakPValue = "10e-6"){
  SampleDir <- normalizePath(path.expand(SampleDir), mustWork = TRUE)
  SampleName <- Get_SampleName(SampleDir, SampleName)

  if(is.null(PeakDir)){
    PeakDir <- file.path(SampleDir, "ChIP_Peaks")
  }
  PeakDir <- path.expand(PeakDir)

  if(!is.null(PeakFile)){
    PeakFile <- path.expand(PeakFile)
    if(file.exists(PeakFile)){
      return(list(file = PeakFile, peaks = Read_Genomewide_Peaks(PeakFile)))
    }
    stop("PeakFile was supplied but not found: ", PeakFile)
  }

  expected_files <- unique(c(
    file.path(PeakDir, paste0(SampleName, "_Genomewide_Peaks.bed")),
    file.path(PeakDir, paste0(SampleName, "_Peaks.bed")),
    file.path(SampleDir, paste0(SampleName, "_Genomewide_Peaks.bed")),
    file.path(SampleDir, paste0(SampleName, "_Peaks.bed")),
    Sys.glob(file.path(PeakDir, "*_Genomewide_Peaks.bed")),
    Sys.glob(file.path(PeakDir, "*_Peaks.bed")),
    Sys.glob(file.path(SampleDir, "*_Genomewide_Peaks.bed")),
    Sys.glob(file.path(SampleDir, "*_Peaks.bed"))
  ))
  expected_files <- expected_files[file.exists(expected_files)]
  if(length(expected_files) > 0){
    return(list(file = expected_files[1],
                peaks = Read_Genomewide_Peaks(expected_files[1])))
  }

  if(ProcessPeaks == TRUE){
    created <- Process_ChIP_Peaks_From_Bam(
      SampleDir = SampleDir,
      SampleName = SampleName,
      PeakDir = PeakDir,
      ChIPBam = ChIPBam,
      InputBam = InputBam,
      BedtoolsCommand = BedtoolsCommand,
      Macs2Command = Macs2Command,
      GenomeSize = GenomeSize,
      PeakPValue = PeakPValue
    )
    return(list(file = created, peaks = Read_Genomewide_Peaks(created)))
  }

  list(file = NA_character_,
       peaks = data.frame(chrom = character(), peakStart = numeric(),
                          peakEnd = numeric(), peakSummit = numeric(),
                          stringsAsFactors = FALSE))
}

Merge_Intervals <- function(intervals){
  if(is.null(intervals) || nrow(intervals) == 0){
    return(intervals)
  }

  intervals <- intervals[order(intervals$start, intervals$end), ]
  merged <- intervals[1, , drop = FALSE]

  for(i in seq_len(nrow(intervals))[-1]){
    last <- nrow(merged)
    if(intervals$start[i] <= merged$end[last]){
      merged$end[last] <- max(merged$end[last], intervals$end[i])
    } else {
      merged <- rbind(merged, intervals[i, , drop = FALSE])
    }
  }

  merged
}

Build_NonPeak_Gaps <- function(chrom, chrom_length, peaks, PeakPaddingBp = 0,
                               MinNoiseGapBp = 1){
  peaks_chr <- peaks[peaks$chrom == chrom, , drop = FALSE]
  if(nrow(peaks_chr) == 0){
    return(data.frame(start = 0, end = chrom_length,
                      length = chrom_length, stringsAsFactors = FALSE))
  }

  intervals <- data.frame(
    start = pmax(0, peaks_chr$peakStart - PeakPaddingBp),
    end = pmin(chrom_length, peaks_chr$peakEnd + PeakPaddingBp),
    stringsAsFactors = FALSE
  )
  intervals <- intervals[is.finite(intervals$start) & is.finite(intervals$end), ]
  intervals <- intervals[intervals$end > intervals$start, ]
  if(nrow(intervals) == 0){
    return(data.frame(start = 0, end = chrom_length,
                      length = chrom_length, stringsAsFactors = FALSE))
  }

  intervals <- Merge_Intervals(intervals)
  gaps <- data.frame(start = numeric(), end = numeric())
  cursor <- 0
  for(i in seq_len(nrow(intervals))){
    if(intervals$start[i] > cursor){
      gaps <- rbind(gaps, data.frame(start = cursor, end = intervals$start[i]))
    }
    cursor <- max(cursor, intervals$end[i])
  }
  if(cursor < chrom_length){
    gaps <- rbind(gaps, data.frame(start = cursor, end = chrom_length))
  }

  if(nrow(gaps) == 0){
    return(data.frame(start = numeric(), end = numeric(),
                      length = numeric(), stringsAsFactors = FALSE))
  }

  gaps$length <- gaps$end - gaps$start
  gaps[gaps$length >= MinNoiseGapBp, , drop = FALSE]
}

Auto_Noise_ChunkSize <- function(all_gaps, NoiseChunkMinBp = 500,
                                 NoiseChunkMaxBp = 10000,
                                 fallback = 2000){
  usable <- all_gaps$length[is.finite(all_gaps$length) &
                              all_gaps$length >= NoiseChunkMinBp]
  if(length(usable) == 0){
    return(fallback)
  }

  raw_chunk <- as.numeric(stats::quantile(usable, probs = 0.10, na.rm = TRUE))
  raw_chunk <- max(NoiseChunkMinBp, min(raw_chunk, NoiseChunkMaxBp))
  max(10, round(raw_chunk / 10) * 10)
}

Auto_Noise_Iterations <- function(NoiseChunkSizeBp){
  if(NoiseChunkSizeBp <= 500){
    6000L
  } else if(NoiseChunkSizeBp <= 1000){
    4000L
  } else if(NoiseChunkSizeBp <= 2500){
    2500L
  } else if(NoiseChunkSizeBp <= 5000){
    1500L
  } else {
    1000L
  }
}

Summarise_Values <- function(values, NoiseStatistic = "median"){
  values <- values[is.finite(values)]
  if(length(values) == 0){
    return(NA_real_)
  }

  if(NoiseStatistic == "mean"){
    mean(values)
  } else {
    stats::median(values)
  }
}

Sample_NonPeak_Noise_Chunks <- function(chrom_df, gaps, NoiseChunkSizeBp,
                                        iterations, NoiseStatistic = "median"){
  empty <- data.frame(center = numeric(), ChIP = numeric(), Input = numeric())
  if(nrow(chrom_df) == 0 || nrow(gaps) == 0 || iterations < 1){
    return(empty)
  }

  centers <- (chrom_df$chromStart + chrom_df$chromEnd) / 2
  gap_lengths <- gaps$length
  gap_prob <- gap_lengths / sum(gap_lengths)

  sampled_gap_index <- sample(seq_len(nrow(gaps)), size = iterations,
                              replace = TRUE, prob = gap_prob)
  out <- vector("list", iterations + nrow(gaps))
  out_i <- 0L

  add_chunk <- function(start, end){
    left <- findInterval(start, centers) + 1
    right <- findInterval(end, centers)
    left <- max(1, left)
    right <- min(length(centers), right)
    if(right < left){
      return(c(center = (start + end) / 2, ChIP = NA_real_, Input = NA_real_))
    }
    idx <- left:right
    c(center = (start + end) / 2,
      ChIP = Summarise_Values(chrom_df$ChIP[idx], NoiseStatistic),
      Input = Summarise_Values(chrom_df$Input[idx], NoiseStatistic))
  }

  for(i in seq_len(nrow(gaps))){
    out_i <- out_i + 1L
    out[[out_i]] <- add_chunk(gaps$start[i], gaps$end[i])
  }

  for(i in seq_len(iterations)){
    gap <- gaps[sampled_gap_index[i], ]
    chunk_size <- min(NoiseChunkSizeBp, gap$length)
    if(gap$length > chunk_size){
      start <- stats::runif(1, min = gap$start, max = gap$end - chunk_size)
    } else {
      start <- gap$start
    }
    end <- start + chunk_size
    out_i <- out_i + 1L
    out[[out_i]] <- add_chunk(start, end)
  }

  sampled <- as.data.frame(do.call(rbind, out[seq_len(out_i)]))
  sampled <- sampled[is.finite(sampled$center), , drop = FALSE]
  sampled
}

Predict_Noise_Profile <- function(bin_centers, sample_centers, sample_values,
                                  fallback_values, NoiseSmoothingSpar = 0.65,
                                  NoiseFloor = 1e-6){
  sample_ok <- is.finite(sample_centers) & is.finite(sample_values)
  sample_centers <- sample_centers[sample_ok]
  sample_values <- sample_values[sample_ok]

  positive_background <- fallback_values[
    is.finite(fallback_values) & fallback_values > NoiseFloor
  ]
  background_floor <- if(length(positive_background) > 0){
    max(
      NoiseFloor,
      as.numeric(stats::quantile(positive_background, probs = 0.01,
                                 na.rm = TRUE, names = FALSE))
    )
  } else {
    NoiseFloor
  }

  fallback <- Summarise_Values(fallback_values, "median")
  if(!is.finite(fallback) || fallback < background_floor){
    fallback <- background_floor
  }

  if(length(sample_values) < 2 || length(unique(sample_centers)) < 2){
    return(rep(fallback, length(bin_centers)))
  }

  ## Background coverage is positive and multiplicative. Smooth it on the
  ## log scale so the spline cannot create negative noise estimates.
  sample_values <- log(pmax(sample_values, background_floor))

  median_no_na <- function(x){
    stats::median(x, na.rm = TRUE)
  }

  sample_df <- data.frame(center = sample_centers, value = sample_values)
  sample_df$center <- round(sample_df$center)
  sample_df <- stats::aggregate(value ~ center, data = sample_df,
                                FUN = median_no_na)
  sample_df <- sample_df[order(sample_df$center), ]

  if(length(unique(sample_df$value)) < 2 || nrow(sample_df) < 4){
    pred <- stats::approx(sample_df$center, sample_df$value, xout = bin_centers,
                          rule = 2, ties = median_no_na)$y
  } else {
    sm <- tryCatch(
      stats::smooth.spline(sample_df$center, sample_df$value,
                           spar = NoiseSmoothingSpar),
      error = function(e) NULL
    )
    if(is.null(sm)){
      pred <- stats::approx(sample_df$center, sample_df$value, xout = bin_centers,
                            rule = 2, ties = median_no_na)$y
    } else {
      pred <- stats::predict(sm, bin_centers)$y
    }
  }

  log_fallback <- log(fallback)
  pred[!is.finite(pred)] <- log_fallback
  pred <- pmin(pmax(pred, min(sample_df$value)), max(sample_df$value))
  pred <- exp(pred)
  pred[pred < background_floor] <- background_floor
  pred
}

Estimate_WholeChromosome_Noise <- function(coverage, peaks, ChromosomeInfo,
                                           NoiseChunkSizeBp = "auto",
                                           NoiseIterations = "auto",
                                           NoiseSeed = 123,
                                           NoiseStatistic = "median",
                                           NoiseSmoothingSpar = 0.65,
                                           NoiseFloor = 1e-6,
                                           PeakPaddingBp = 0,
                                           MinNoiseGapBp = 200,
                                           NoiseChunkMinBp = 500,
                                           NoiseChunkMaxBp = 10000){
  NoiseStatistic <- match.arg(NoiseStatistic, choices = c("median", "mean"))

  all_gaps <- data.frame(chrom = character(), start = numeric(), end = numeric(),
                         length = numeric(), stringsAsFactors = FALSE)
  for(chr in unique(coverage$chrom)){
    chr_len <- ChromosomeInfo$length[match(chr, ChromosomeInfo$chrom)]
    if(!is.finite(chr_len)){
      chr_len <- max(coverage$chromEnd[coverage$chrom == chr], na.rm = TRUE)
    }
    gaps_chr <- Build_NonPeak_Gaps(chr, chr_len, peaks,
                                   PeakPaddingBp = PeakPaddingBp,
                                   MinNoiseGapBp = MinNoiseGapBp)
    if(nrow(gaps_chr) > 0){
      all_gaps <- rbind(all_gaps, data.frame(chrom = chr, gaps_chr))
    }
  }

  if(nrow(all_gaps) == 0){
    stop("No non-peak regions are available for noise estimation.")
  }

  if(identical(NoiseChunkSizeBp, "auto")){
    NoiseChunkSizeBp <- Auto_Noise_ChunkSize(
      all_gaps,
      NoiseChunkMinBp = NoiseChunkMinBp,
      NoiseChunkMaxBp = NoiseChunkMaxBp
    )
  }
  NoiseChunkSizeBp <- as.numeric(NoiseChunkSizeBp)
  Validate_Positive_Number(NoiseChunkSizeBp, "NoiseChunkSizeBp")

  if(identical(NoiseIterations, "auto")){
    NoiseIterations <- Auto_Noise_Iterations(NoiseChunkSizeBp)
  }
  NoiseIterations <- as.integer(NoiseIterations)
  if(!is.finite(NoiseIterations) || NoiseIterations < 1){
    stop("NoiseIterations must be a positive integer or 'auto'.")
  }

  set.seed(NoiseSeed)
  total_gap_bp <- sum(all_gaps$length, na.rm = TRUE)
  coverage$ChIP_noise <- NA_real_
  coverage$Input_noise <- NA_real_

  for(chr in unique(coverage$chrom)){
    chr_idx <- which(coverage$chrom == chr)
    chrom_df <- coverage[chr_idx, , drop = FALSE]
    gaps_chr <- all_gaps[all_gaps$chrom == chr, c("start", "end", "length"),
                         drop = FALSE]
    if(nrow(gaps_chr) == 0){
      next
    }

    chr_gap_bp <- sum(gaps_chr$length, na.rm = TRUE)
    iterations_chr <- max(50L, round(NoiseIterations * chr_gap_bp / total_gap_bp))
    sampled <- Sample_NonPeak_Noise_Chunks(
      chrom_df = chrom_df,
      gaps = gaps_chr,
      NoiseChunkSizeBp = NoiseChunkSizeBp,
      iterations = iterations_chr,
      NoiseStatistic = NoiseStatistic
    )

    bin_centers <- (chrom_df$chromStart + chrom_df$chromEnd) / 2
    nonpeak_idx <- rep(FALSE, nrow(chrom_df))
    for(i in seq_len(nrow(gaps_chr))){
      nonpeak_idx <- nonpeak_idx |
        (bin_centers >= gaps_chr$start[i] & bin_centers <= gaps_chr$end[i])
    }

    coverage$ChIP_noise[chr_idx] <- Predict_Noise_Profile(
      bin_centers = bin_centers,
      sample_centers = sampled$center,
      sample_values = sampled$ChIP,
      fallback_values = chrom_df$ChIP[nonpeak_idx],
      NoiseSmoothingSpar = NoiseSmoothingSpar,
      NoiseFloor = NoiseFloor
    )

    coverage$Input_noise[chr_idx] <- Predict_Noise_Profile(
      bin_centers = bin_centers,
      sample_centers = sampled$center,
      sample_values = sampled$Input,
      fallback_values = chrom_df$Input[nonpeak_idx],
      NoiseSmoothingSpar = NoiseSmoothingSpar,
      NoiseFloor = NoiseFloor
    )
  }

  list(
    coverage = coverage,
    parameters = list(
      NoiseChunkSizeBp = NoiseChunkSizeBp,
      NoiseIterations = NoiseIterations,
      NoiseSeed = NoiseSeed,
      NoiseStatistic = NoiseStatistic,
      NoiseSmoothingSpar = NoiseSmoothingSpar,
      NoiseFloor = NoiseFloor,
      PeakPaddingBp = PeakPaddingBp,
      MinNoiseGapBp = MinNoiseGapBp,
      NonPeakGapCount = nrow(all_gaps),
      NonPeakGapBp = total_gap_bp
    )
  )
}

Prepare_WholeChromosome_ChIP_Data <- function(SampleDir, SampleName = NULL,
                                              CoverageDir = NULL,
                                              CountsFile = NULL,
                                              RatioFile = NULL,
                                              OutputDir = NULL,
                                              CoverageColumn = 5,
                                              FilterNoise = FALSE,
                                              PeakDir = NULL,
                                              PeakFile = NULL,
                                              ProcessPeaks = FALSE,
                                              ChIPBam = NULL,
                                              InputBam = NULL,
                                              BedtoolsCommand = NULL,
                                              Macs2Command = "macs2",
                                              GenomeSize = 12157105,
                                              PeakPValue = "10e-6",
                                              NoiseChunkSizeBp = "auto",
                                              NoiseIterations = "auto",
                                              NoiseSeed = 123,
                                              NoiseStatistic = "median",
                                              NoiseSmoothingSpar = 0.65,
                                              NoiseFloor = 1e-6,
                                              PeakPaddingBp = 0,
                                              MinNoiseGapBp = 200,
                                              NoiseChunkMinBp = 500,
                                              NoiseChunkMaxBp = 10000,
                                              SaveCollapsedTable = TRUE){
  Validate_Logical(FilterNoise, "FilterNoise")
  Validate_Logical(ProcessPeaks, "ProcessPeaks")
  Validate_Logical(SaveCollapsedTable, "SaveCollapsedTable")

  RequestedSampleName <- SampleName
  SampleDir <- normalizePath(path.expand(SampleDir), mustWork = TRUE)
  SampleName <- Get_SampleName(SampleDir, SampleName)
  ratio_file <- Resolve_Ratio_Profile_File(
    SampleDir = SampleDir,
    SampleName = SampleName,
    RatioFile = RatioFile
  )

  if(!is.na(ratio_file)){
    message("Using ratio profile table: ", ratio_file)
    coverage <- Read_Ratio_Profile_Bed(ratio_file,
                                       SampleName = RequestedSampleName)
    if(is.null(RequestedSampleName)){
      SampleName <- attr(coverage, "SampleName")
    }
    paths <- list(
      SampleDir = SampleDir,
      SampleName = SampleName,
      CoverageDir = dirname(ratio_file),
      CountsFile = NA_character_,
      RatioFile = ratio_file,
      ChIPWatson = NA_character_,
      ChIPCrick = NA_character_,
      InputWatson = NA_character_,
      InputCrick = NA_character_
    )
  } else {
    normal_counts_file <- Find_Normal_Counts_File(
      SampleDir = SampleDir,
      SampleName = SampleName,
      CountsFile = CountsFile
    )
  }

  if(is.na(ratio_file) && !is.na(normal_counts_file)){
    message("Using normal counts table: ", normal_counts_file)
    coverage <- Read_Normal_Counts_Bed(normal_counts_file,
                                       SampleName = SampleName)
    paths <- list(
      SampleDir = SampleDir,
      SampleName = SampleName,
      CoverageDir = dirname(normal_counts_file),
      CountsFile = normal_counts_file,
      RatioFile = NA_character_,
      ChIPWatson = NA_character_,
      ChIPCrick = NA_character_,
      InputWatson = NA_character_,
      InputCrick = NA_character_
    )
  } else if(is.na(ratio_file)) {
    collapsed <- Read_StrandCollapsed_Coverage(
      SampleDir = SampleDir,
      SampleName = SampleName,
      CoverageDir = CoverageDir,
      CoverageColumn = CoverageColumn
    )
    coverage <- collapsed$coverage
    paths <- collapsed$paths
    paths$CountsFile <- NA_character_
    paths$RatioFile <- NA_character_
  }

  OutputDir <- Resolve_Output_Directory(OutputDir, paths$SampleDir)

  peak_result <- Get_Or_Create_Genomewide_Peaks(
    SampleDir = paths$SampleDir,
    SampleName = paths$SampleName,
    PeakDir = PeakDir,
    PeakFile = PeakFile,
    ProcessPeaks = ProcessPeaks,
    ChIPBam = ChIPBam,
    InputBam = InputBam,
    BedtoolsCommand = BedtoolsCommand,
    Macs2Command = Macs2Command,
    GenomeSize = GenomeSize,
    PeakPValue = PeakPValue
  )

  noise_parameters <- NULL
  if(FilterNoise == TRUE){
    if(nrow(peak_result$peaks) == 0){
      stop("FilterNoise=TRUE requires a genomewide peak list. Provide PeakFile, ",
           "place <SampleName>_Genomewide_Peaks.bed in ChIP_Peaks, or set ProcessPeaks=TRUE.")
    }

    noise_result <- Estimate_WholeChromosome_Noise(
      coverage = coverage,
      peaks = peak_result$peaks,
      ChromosomeInfo = Get_SacCer3_Chromosome_Table(),
      NoiseChunkSizeBp = NoiseChunkSizeBp,
      NoiseIterations = NoiseIterations,
      NoiseSeed = NoiseSeed,
      NoiseStatistic = NoiseStatistic,
      NoiseSmoothingSpar = NoiseSmoothingSpar,
      NoiseFloor = NoiseFloor,
      PeakPaddingBp = PeakPaddingBp,
      MinNoiseGapBp = MinNoiseGapBp,
      NoiseChunkMinBp = NoiseChunkMinBp,
      NoiseChunkMaxBp = NoiseChunkMaxBp
    )
    coverage <- noise_result$coverage
    noise_parameters <- noise_result$parameters
  } else {
    message("FilterNoise=FALSE: ChIP_noise and Input_noise columns are left as NA by design.")
  }

  table_suffix <- if(!is.na(paths$RatioFile)){
    "_ratio_profile_table.bed"
  } else if(FilterNoise == TRUE){
    "_strandcollapsed_ChIP_Input_noise_table.bed"
  } else {
    "_strandcollapsed_ChIP_Input_table.bed"
  }
  data_file <- file.path(
    OutputDir,
    paste0(paths$SampleName, table_suffix)
  )
  if(SaveCollapsedTable == TRUE){
    write.table(coverage, file = data_file, quote = FALSE,
                row.names = FALSE, sep = "\t")
  }

  list(
    coverage = coverage,
    peaks = peak_result$peaks,
    peak_file = peak_result$file,
    counts_file = paths$CountsFile,
    ratio_file = paths$RatioFile,
    data_file = if(SaveCollapsedTable) data_file else NA_character_,
    paths = paths,
    noise_parameters = noise_parameters
  )
}

Profile_Values <- function(coverage, PlotMode, Log2Profile = TRUE,
                           UseRatioLog2 = TRUE,
                           ApplyRatioPValueSign = TRUE,
                           RatioPValueCutoff = 10e-2){
  PlotMode <- match.arg(PlotMode, choices = c("ChIP", "ChIP_Input",
                                              "ChIP_Noise", "Clean"))
  Validate_Logical(Log2Profile, "Log2Profile")
  Validate_Logical(UseRatioLog2, "UseRatioLog2")
  Validate_Logical(ApplyRatioPValueSign, "ApplyRatioPValueSign")
  Validate_Nonnegative_Number(RatioPValueCutoff, "RatioPValueCutoff")

  if(PlotMode == "ChIP_Input" && UseRatioLog2 == TRUE &&
     "ratio_log2" %in% colnames(coverage)){
    values <- coverage$ratio_log2
    if(ApplyRatioPValueSign == TRUE &&
       "ratio_pvalue" %in% colnames(coverage)){
      nonsig <- is.finite(coverage$ratio_pvalue) &
        coverage$ratio_pvalue > RatioPValueCutoff
      values[nonsig] <- -abs(values[nonsig])
    }
    return(list(values = values,
                label = "log2(IP / Input)",
                baseline = 0,
                is_log2 = TRUE))
  }

  if(PlotMode == "ChIP"){
    values <- coverage$ChIP
    label <- "IP coverage"
    baseline <- NA_real_
  } else if(PlotMode == "ChIP_Input"){
    values <- Safe_Divide(coverage$ChIP, coverage$Input)
    label <- "IP / Input"
    baseline <- 1
  } else if(PlotMode == "ChIP_Noise"){
    if(all(is.na(coverage$ChIP_noise))){
      stop("PlotMode='ChIP_Noise' requires FilterNoise=TRUE.")
    }
    values <- Safe_Divide(coverage$ChIP, coverage$ChIP_noise)
    label <- "IP / Noise"
    baseline <- 1
  } else {
    if(all(is.na(coverage$ChIP_noise)) || all(is.na(coverage$Input_noise))){
      stop("PlotMode='Clean' requires FilterNoise=TRUE.")
    }
    chip_noise <- Safe_Divide(coverage$ChIP, coverage$ChIP_noise)
    input_noise <- Safe_Divide(coverage$Input, coverage$Input_noise)
    values <- Safe_Divide(chip_noise, input_noise)
    label <- "(IP / Noise) / (Input / Noise)"
    baseline <- 1
  }

  if(Log2Profile == TRUE){
    values <- Safe_Log2(values)
    label <- paste0("log2(", label, ")")
    baseline <- 0
  }

  list(values = values, label = label, baseline = baseline,
       is_log2 = Log2Profile)
}

ARS_Label_Color <- function(stat){
  stat <- tolower(as.character(stat))
  col <- rep("gray10", length(stat))
  col[stat == "early"] <- "red"
  col[stat == "late"] <- "blue"
  col
}

ARS_Label_Levels <- function(x, min_gap = 50, max_levels = 5){
  if(length(x) == 0){
    return(integer())
  }

  levels <- integer(length(x))
  last_x <- rep(-Inf, max_levels)
  for(idx in order(x)){
    available <- which((x[idx] - last_x) >= min_gap)
    if(length(available) == 0){
      level <- which.min(last_x)
    } else {
      level <- available[1]
    }
    levels[idx] <- level - 1L
    last_x[level] <- x[idx]
  }

  levels
}

Draw_ARS_Origins <- function(ARS_chr, window, steps, Ylim_sc,
                             BaselineValue = 0,
                             ARSColor = "purple",
                             ARSFill = "yellow",
                             ShowARSLabels = TRUE,
                             ARSLabelCex = 0.75,
                             ARSLabelMinGap = 350,
                             ARSLabelLevels = 5){
  if(is.null(ARS_chr) || nrow(ARS_chr) == 0){
    return(invisible(NULL))
  }

  y_range <- diff(Ylim_sc)
  if(!is.finite(y_range) || y_range <= 0){
    y_range <- 1
  }

  ars_x <- ((ARS_chr$chromStart + ARS_chr$chromEnd) / 2 - window$S) / steps
  ars_y <- if(is.finite(BaselineValue) &&
              BaselineValue >= min(Ylim_sc) &&
              BaselineValue <= max(Ylim_sc)){
    BaselineValue
  } else if(0 >= min(Ylim_sc) && 0 <= max(Ylim_sc)){
    0
  } else {
    min(Ylim_sc) + 0.08 * y_range
  }

  points(ars_x, rep(ars_y, length(ars_x)), pch = 21,
         bg = ARSFill, col = ARSColor, lwd = 1.4, cex = 1.05)

  if(ShowARSLabels == FALSE){
    return(invisible(NULL))
  }

  label_levels <- ARS_Label_Levels(ars_x, min_gap = ARSLabelMinGap,
                                   max_levels = ARSLabelLevels)
  label_y <- max(Ylim_sc) + (0.16 + 0.12 * label_levels) * y_range
  triangle_y <- max(Ylim_sc) + (0.055 + 0.12 * label_levels) * y_range
  label_col <- if("stat" %in% colnames(ARS_chr)){
    ARS_Label_Color(ARS_chr$stat)
  } else {
    rep("gray10", nrow(ARS_chr))
  }

  segments(x0 = ars_x, x1 = ars_x, y0 = max(Ylim_sc), y1 = triangle_y,
           col = grDevices::adjustcolor(ARSColor, alpha.f = 0.35),
           lwd = 1.2, xpd = NA)
  points(ars_x, triangle_y, pch = 24, bg = ARSFill, col = ARSColor,
         lwd = 1.3, cex = 0.95, xpd = NA)
  text(ars_x, label_y, labels = ARS_chr$name, xpd = NA,
       cex = ARSLabelCex, col = label_col)

  invisible(NULL)
}

Bins_Overlap_Peaks <- function(chrom, bin_start, bin_end, peaks){
  overlap <- rep(FALSE, length(bin_start))
  if(is.null(peaks) || nrow(peaks) == 0 || length(bin_start) == 0){
    return(overlap)
  }

  peaks_chr <- peaks[peaks$chrom == chrom &
                       peaks$peakEnd >= min(bin_start, na.rm = TRUE) &
                       peaks$peakStart <= max(bin_end, na.rm = TRUE),
                     , drop = FALSE]
  if(nrow(peaks_chr) == 0){
    return(overlap)
  }

  for(i in seq_len(nrow(peaks_chr))){
    overlap <- overlap |
      (bin_end >= peaks_chr$peakStart[i] & bin_start <= peaks_chr$peakEnd[i])
  }

  overlap
}

Peak_Supported_Positive_Bars <- function(chrom, bin_start, bin_end, y_plot,
                                         peaks,
                                         PeakIslandBoundaryFraction = 0.25,
                                         PeakIslandMaxFlankBp = 250){
  Validate_Fraction(PeakIslandBoundaryFraction, "PeakIslandBoundaryFraction")
  Validate_Nonnegative_Number(PeakIslandMaxFlankBp, "PeakIslandMaxFlankBp")

  positive <- is.finite(y_plot) & y_plot > 0
  qualified <- rep(FALSE, length(positive))
  if(!any(positive) || is.null(peaks) || nrow(peaks) == 0 ||
     length(bin_start) == 0){
    return(qualified)
  }

  search_start <- min(bin_start, na.rm = TRUE) - PeakIslandMaxFlankBp
  search_end <- max(bin_end, na.rm = TRUE) + PeakIslandMaxFlankBp
  peaks_chr <- peaks[peaks$chrom == chrom &
                       peaks$peakEnd >= search_start &
                       peaks$peakStart <= search_end,
                     , drop = FALSE]
  if(nrow(peaks_chr) == 0){
    return(qualified)
  }

  for(i in seq_len(nrow(peaks_chr))){
    peak_start <- peaks_chr$peakStart[i]
    peak_end <- peaks_chr$peakEnd[i]
    peak_overlap <- positive &
      bin_end >= peak_start & bin_start <= peak_end
    if(!any(peak_overlap)){
      next
    }

    peak_height <- max(y_plot[peak_overlap], na.rm = TRUE)
    if(!is.finite(peak_height) || peak_height <= 0){
      next
    }

    in_peak_flank <- bin_end >= (peak_start - PeakIslandMaxFlankBp) &
      bin_start <= (peak_end + PeakIslandMaxFlankBp)
    min_signal <- peak_height * PeakIslandBoundaryFraction
    candidate <- positive & in_peak_flank & y_plot >= min_signal
    if(!any(candidate)){
      next
    }

    candidate_runs <- rle(candidate)
    candidate_end <- cumsum(candidate_runs$lengths)
    candidate_start <- candidate_end - candidate_runs$lengths + 1L
    for(j in seq_along(candidate_runs$values)){
      if(candidate_runs$values[j] == TRUE){
        idx <- candidate_start[j]:candidate_end[j]
        if(any(peak_overlap[idx])){
          qualified[idx] <- TRUE
        }
      }
    }
  }

  qualified
}

Plot_WholeChromosome_Profile_Window <- function(ProfileTable, annotations, peaks,
                                                window, steps = 10,
                                                FullWindowBp = NULL,
                                                Ylim_sc, yAxis_reads,
                                                PlotStyle = "lines",
                                                SmoothSignal = TRUE,
                                                SmoothingSpar = NULL,
                                                ShowBaseline = TRUE,
                                                BaselineValue = NA_real_,
                                                SignalColor = "gray25",
                                                PositiveSignalColor = "red",
                                                NegativeSignalColor = "gray70",
                                                ProfileBarLwd = 0.12,
                                                ProfileLineLwd = 1.15,
                                                RedRequiresPeakOverlap = TRUE,
                                                PeakIslandBoundaryFraction = 0.35,
                                                PeakIslandMaxFlankBp = 500,
                                                PlotPeaks = TRUE,
                                                PeakColor = "firebrick2",
                                                ARSColor = "purple",
                                                ARSFill = "yellow",
                                                ShowARSLabels = TRUE,
                                                ARSLabelCex = 0.75,
                                                ARSLabelMinGap = 350,
                                                ARSLabelLevels = 5,
                                                CentromereColor = "darkgreen",
                                                CentromereLwd = 1.6,
                                                WindowTitle = "",
                                                YLabel = "",
                                                ProfileMar = c(0.55, 4.3, 1.0, 1.0)){
  Coverage <- ProfileTable[ProfileTable$chrom == window$chrom &
                             ProfileTable$chromStart >= window$S &
                             ProfileTable$chromStart < window$E, ]
  panel_width <- if(is.null(FullWindowBp)) window$E - window$S else FullWindowBp
  x_max <- panel_width / steps
  y_range <- diff(Ylim_sc)
  if(!is.finite(y_range) || y_range <= 0){
    y_range <- 1
  }
  Validate_Logical(RedRequiresPeakOverlap, "RedRequiresPeakOverlap")
  Validate_Fraction(PeakIslandBoundaryFraction, "PeakIslandBoundaryFraction")
  Validate_Nonnegative_Number(PeakIslandMaxFlankBp, "PeakIslandMaxFlankBp")

  par(mar = ProfileMar)
  if(nrow(Coverage) == 0){
    plot(NA, xlim = c(0, x_max), ylim = Ylim_sc, ylab = " ", xlab = " ",
         xaxt = "n", yaxt = "n", bty = "n", xaxs = "i")
  } else {
    x_vals <- (Coverage$chromStart - window$S) / steps
    y_vals <- Coverage$plot_signal
    if(SmoothSignal == TRUE){
      y_plot <- Safe_Smooth(x_vals, y_vals, spar = SmoothingSpar)
      x_plot <- x_vals
    } else {
      x_plot <- x_vals
      y_plot <- y_vals
    }

    plot_type <- if(PlotStyle == "lines") "l" else "h"
    plot_lwd <- if(PlotStyle == "lines") ProfileLineLwd else ProfileBarLwd
    plot_col <- if(PlotStyle == "lines"){
      SignalColor
    } else {
      positive_bar <- is.finite(y_plot) & y_plot > 0
      if(RedRequiresPeakOverlap == TRUE){
        positive_bar <- Peak_Supported_Positive_Bars(
          chrom = window$chrom,
          bin_start = Coverage$chromStart,
          bin_end = Coverage$chromEnd,
          y_plot = y_plot,
          peaks = peaks,
          PeakIslandBoundaryFraction = PeakIslandBoundaryFraction,
          PeakIslandMaxFlankBp = PeakIslandMaxFlankBp
        )
      }
      ifelse(positive_bar, PositiveSignalColor, NegativeSignalColor)
    }
    suppressWarnings(
      plot(x_plot, y_plot, type = plot_type, xlim = c(0, x_max), ylim = Ylim_sc,
           col = plot_col, ylab = " ", xlab = " ", xaxt = "n", yaxt = "n",
           lwd = plot_lwd, bty = "n", cex.lab = 1, las = 2, xaxs = "i")
    )
    segments(x0 = x_vals, x1 = (Coverage$chromEnd - window$S) / steps,
             y0 = 0, y1 = 0, lwd = 0.25,
             col = grDevices::adjustcolor("gray40", alpha.f = 0.35))
  }

  axis(side = 2, at = yAxis_reads, labels = signif(yAxis_reads, 3),
       line = 0, tick = TRUE, lwd.ticks = 1.2, las = 2, cex.axis = 0.8)
  abline(h = yAxis_reads, lwd = 0.05,
         col = grDevices::rgb(112, 128, 144, alpha = 150, maxColorValue = 255))

  if(ShowBaseline == TRUE && is.finite(BaselineValue) &&
     BaselineValue >= min(Ylim_sc) && BaselineValue <= max(Ylim_sc)){
    abline(h = BaselineValue, lwd = 0.9, lty = 3,
           col = grDevices::rgb(40, 40, 40, alpha = 220, maxColorValue = 255))
  }

  if(PlotPeaks == TRUE && !is.null(peaks) && nrow(peaks) > 0){
    peaks_chr <- peaks[peaks$chrom == window$chrom &
                         peaks$peakEnd >= window$S &
                         peaks$peakStart <= window$E, , drop = FALSE]
    if(nrow(peaks_chr) > 0){
      peak_x0 <- (pmax(peaks_chr$peakStart, window$S) - window$S) / steps
      peak_x1 <- (pmin(peaks_chr$peakEnd, window$E) - window$S) / steps
      rect(peak_x0, max(Ylim_sc) - 0.08 * y_range,
           peak_x1, max(Ylim_sc) - 0.015 * y_range,
           col = grDevices::adjustcolor(PeakColor, alpha.f = 0.35),
           border = NA)
    }
  }

  ARS_chr <- Get_Window_Features(annotations$ARS, window$chrom, window$S, window$E)
  Draw_ARS_Origins(ARS_chr = ARS_chr, window = window, steps = steps,
                   Ylim_sc = Ylim_sc, BaselineValue = BaselineValue,
                   ARSColor = ARSColor, ARSFill = ARSFill,
                   ShowARSLabels = ShowARSLabels,
                   ARSLabelCex = ARSLabelCex,
                   ARSLabelMinGap = ARSLabelMinGap,
                   ARSLabelLevels = ARSLabelLevels)

  CEN_chr <- Get_Window_Features(annotations$Centromeres, window$chrom,
                                 window$S, window$E)
  if(nrow(CEN_chr) > 0){
    cen_x <- ((CEN_chr$chromStart + CEN_chr$chromEnd) / 2 - window$S) / steps
    segments(x0 = cen_x, x1 = cen_x, y0 = min(Ylim_sc), y1 = max(Ylim_sc),
             lwd = CentromereLwd, col = CentromereColor)
    text(x = cen_x, y = max(Ylim_sc) - 0.11 * y_range, labels = "CEN",
         cex = 0.65, col = CentromereColor, pos = 4)
  }

  title(main = WindowTitle, col = "gray35", adj = 0, cex.main = 0.95, line = 0.25)
  mtext(side = 2, line = 2.75, at = mean(Ylim_sc), cex = 0.85, YLabel)
  box(col = "gray45")
}

Draw_Annotation_Blocks <- function(features, chrom, S, E, steps, y, height,
                                   col, border = col, min_width = 0.6){
  features_chr <- Get_Window_Features(features, chrom, S, E)
  if(nrow(features_chr) == 0){
    return(invisible(NULL))
  }

  x0 <- (pmax(features_chr$chromStart, S) - S) / steps
  x1 <- (pmin(features_chr$chromEnd, E) - S) / steps
  too_small <- (x1 - x0) < min_width
  x_mid <- (x0 + x1) / 2
  x0[too_small] <- x_mid[too_small] - min_width / 2
  x1[too_small] <- x_mid[too_small] + min_width / 2
  rect(xleft = x0, ybottom = y - height, xright = x1, ytop = y + height,
       col = col, border = border, lwd = 0.4)

  invisible(NULL)
}

Plot_WholeChromosome_Annotation_Window <- function(annotations, window,
                                                   steps = 10, ORFLwd = 2.5,
                                                   FullWindowBp = NULL,
                                                   ShowXAxis = TRUE,
                                                   ShowXLabel = TRUE,
                                                   FeatureTracks = c("ORF", "Ty",
                                                                     "TER", "tRNA"),
                                                   AnnotationMar = c(1.6, 4.3, 0.35, 1.0)){
  FeatureTracks <- intersect(FeatureTracks, c("ORF", "Ty", "TER", "tRNA"))
  panel_width <- if(is.null(FullWindowBp)) window$E - window$S else FullWindowBp
  x_max <- panel_width / steps

  par(mar = AnnotationMar)
  plot(NA, xlim = c(0, x_max), ylim = c(0, 1), ylab = " ", xlab = " ",
       xaxt = "n", yaxt = "n", bty = "n", xaxs = "i", yaxs = "i")
  abline(h = c(0.16, 0.30, 0.48, 0.66, 0.82), lwd = 0.05,
         col = grDevices::rgb(112, 128, 144, alpha = 150, maxColorValue = 255))

  if("ORF" %in% FeatureTracks){
    ORFs_chr <- Get_Window_Features(annotations$ORFs, window$chrom, window$S, window$E)
    if(nrow(ORFs_chr) > 0){
      ORFs_chr$plotStart <- pmax(ORFs_chr$chromStart, window$S)
      ORFs_chr$plotEnd <- pmin(ORFs_chr$chromEnd, window$E)
      ORFs_chr <- ORFs_chr[ORFs_chr$plotEnd > ORFs_chr$plotStart, ]
      for(i in seq_len(nrow(ORFs_chr))){
        x0 <- (ORFs_chr$plotStart[i] - window$S) / steps
        x1 <- (ORFs_chr$plotEnd[i] - window$S) / steps
        y <- if(ORFs_chr$strand[i] == "+") 0.82 else 0.66
        col <- if(ORFs_chr$strand[i] == "+") "brown3" else "cornflowerblue"
        if(ORFs_chr$strand[i] == "+"){
          arrows(x0 = x0, y0 = y, x1 = x1, y1 = y, length = 0.04,
                 angle = 25, code = 2, lwd = ORFLwd, col = col)
        } else {
          arrows(x0 = x1, y0 = y, x1 = x0, y1 = y, length = 0.04,
                 angle = 25, code = 2, lwd = ORFLwd, col = col)
        }
      }
    }
  }

  if("Ty" %in% FeatureTracks){
    Draw_Annotation_Blocks(annotations$Ty, window$chrom, window$S, window$E,
                           steps, y = 0.48, height = 0.055,
                           col = grDevices::adjustcolor("mediumpurple3",
                                                        alpha.f = 0.75),
                           border = "mediumpurple4")
  }
  if("TER" %in% FeatureTracks){
    Draw_Annotation_Blocks(annotations$TER, window$chrom, window$S, window$E,
                           steps, y = 0.30, height = 0.055,
                           col = grDevices::adjustcolor("orange2",
                                                        alpha.f = 0.65),
                           border = "orange4")
  }
  if("tRNA" %in% FeatureTracks){
    Draw_Annotation_Blocks(annotations$tRNA, window$chrom, window$S, window$E,
                           steps, y = 0.16, height = 0.045,
                           col = grDevices::adjustcolor("seagreen3",
                                                        alpha.f = 0.75),
                           border = "seagreen4")
  }

  axis(side = 2, at = c(0.82, 0.66, 0.48, 0.30, 0.16),
       labels = c("ORF+", "ORF-", "Ty", "TER", "tRNA"),
       line = 0, tick = TRUE, lwd.ticks = 1.0, las = 2, cex.axis = 0.75)
  if(ShowXAxis == TRUE){
    tick_by <- max(10000, round(panel_width / 5, -3))
    ticks_bp <- seq(ceiling(window$S / tick_by) * tick_by, window$E, by = tick_by)
    ticks_x <- (ticks_bp - window$S) / steps
    axis(side = 1, at = ticks_x, labels = round(ticks_bp / 1000),
         line = 0, tick = FALSE, cex.axis = 0.8)
    if(ShowXLabel == TRUE){
      title(xlab = "Chromosomal Coordinates (Kbp)", col = "gray",
            cex.lab = 0.9, line = 0.75)
    }
  }
  box(col = "gray45")
}

Plot_WholeChromosome_Spacer <- function(ShowLegend = FALSE,
                                        PeakColor = "firebrick2",
                                        ARSColor = "purple",
                                        ARSFill = "yellow",
                                        PlotPeaks = TRUE){
  par(mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", bty = "n", xaxs = "i", yaxs = "i")

  if(ShowLegend == TRUE){
    y <- 0.55
    x_ars <- if(PlotPeaks == TRUE) 0.40 else 0.46
    points(x_ars, y, pch = 21, bg = grDevices::adjustcolor(ARSFill, alpha.f = 0),
           col = ARSColor,
           lwd = 1.4, cex = 1.1)
    text(x_ars + 0.035, y, "ARS", adj = 0, cex = 0.8, col = "gray25")

    if(PlotPeaks == TRUE){
      x_peak <- 0.56
      rect(x_peak, y - 0.035, x_peak + 0.035, y + 0.035,
           col = grDevices::adjustcolor(PeakColor, alpha.f = 0.35),
           border = NA)
      text(x_peak + 0.055, y, "Peaks", adj = 0, cex = 0.8, col = "gray25")
    }
  }
}

Plot_WholeChromosome_Blank_Window <- function(ProfileMar = c(0.55, 4.3, 1.0, 1.0),
                                              AnnotationMar = c(1.6, 4.3, 0.35, 1.0),
                                              ShowLegend = FALSE,
                                              PeakColor = "firebrick2",
                                              PlotPeaks = TRUE){
  par(mar = ProfileMar)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", bty = "n", xaxs = "i", yaxs = "i")

  par(mar = AnnotationMar)
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", bty = "n", xaxs = "i", yaxs = "i")

  Plot_WholeChromosome_Spacer(
    ShowLegend = ShowLegend,
    PeakColor = PeakColor,
    ARSColor = "purple",
    ARSFill = "yellow",
    PlotPeaks = PlotPeaks
  )
}

WholeChromosome_ChIP_Plotter <- function(SampleDir,
                                         Chromosome,
                                         SampleName = NULL,
                                         CoverageDir = NULL,
                                         CountsFile = NULL,
                                         RatioFile = NULL,
                                         OutputDir = NULL,
                                         AnnotationDir = NULL,
                                         ORFFile = NULL,
                                         ARSFile = NULL,
                                         TERFile = NULL,
                                         TyFile = NULL,
                                         tRNAFile = NULL,
                                         CentromereFile = NULL,
                                         PeakDir = NULL,
                                         PeakFile = NULL,
                                         ProcessPeaks = FALSE,
                                         ChIPBam = NULL,
                                         InputBam = NULL,
                                         BedtoolsCommand = NULL,
                                         Macs2Command = "macs2",
                                         GenomeSize = 12157105,
                                         PeakPValue = "10e-6",
                                         CoverageColumn = 5,
                                         FilterNoise = FALSE,
                                         NoiseChunkSizeBp = "auto",
                                         NoiseIterations = "auto",
                                         NoiseSeed = 123,
                                         NoiseStatistic = "median",
                                         NoiseSmoothingSpar = 0.65,
                                         NoiseFloor = 1e-6,
                                         PeakPaddingBp = 0,
                                         MinNoiseGapBp = 200,
                                         NoiseChunkMinBp = 500,
                                         NoiseChunkMaxBp = 10000,
                                         PlotMode = c("ChIP", "ChIP_Input",
                                                     "ChIP_Noise", "Clean"),
                                         WindowSizeKb = 50,
                                         PanelsPerPage = 3,
                                         PlotStyle = "hist",
                                         Log2Profile = TRUE,
                                         UseRatioLog2 = TRUE,
                                         ApplyRatioPValueSign = TRUE,
                                         RatioPValueCutoff = 10e-2,
                                         Log2YMin = -1,
                                         SmoothSignal = TRUE,
                                         SmoothingSpar = NULL,
                                         ShowBaseline = TRUE,
                                         BaselineValue = NULL,
                                         FixedYlim = NULL,
                                         Y_axis_scale = 8,
                                         SignalColor = "gray25",
                                         PositiveSignalColor = "red",
                                         NegativeSignalColor = "gray70",
                                         ProfileLineLwd = 1.15,
                                         ProfileBarLwd = 0.35,
                                         RedRequiresPeakOverlap = TRUE,
                                         PeakIslandBoundaryFraction = 0.35,
                                         PeakIslandMaxFlankBp = 500,
                                         ORFLwd = 2.5,
                                         PlotPeaks = TRUE,
                                         PeakColor = "firebrick2",
                                         ShowARSLabels = TRUE,
                                         ARSLabelCex = 0.75,
                                         ARSLabelMinGap = 350,
                                         ARSLabelLevels = 5,
                                         ShowFeatureLegend = TRUE,
                                         FeatureTracks = c("ORF", "Ty",
                                                           "TER", "tRNA"),
                                         SaveCollapsedTable = TRUE,
                                         PdfWidth = 11,
                                         PdfHeight = 11,
                                         ProfileMar = c(0.55, 4.3, 1.75, 1.0),
                                         AnnotationMar = c(1.6, 4.3, 0.35, 1.0),
                                         OuterMar = c(3.2, 2.8, 2.8, 1.4),
                                         PanelGapHeight = 0.45){
  Validate_Logical(ProcessPeaks, "ProcessPeaks")
  Validate_Logical(FilterNoise, "FilterNoise")
  Validate_Logical(Log2Profile, "Log2Profile")
  Validate_Logical(UseRatioLog2, "UseRatioLog2")
  Validate_Logical(ApplyRatioPValueSign, "ApplyRatioPValueSign")
  Validate_Logical(SmoothSignal, "SmoothSignal")
  Validate_Logical(ShowBaseline, "ShowBaseline")
  Validate_Logical(PlotPeaks, "PlotPeaks")
  Validate_Logical(RedRequiresPeakOverlap, "RedRequiresPeakOverlap")
  Validate_Fraction(PeakIslandBoundaryFraction, "PeakIslandBoundaryFraction")
  Validate_Nonnegative_Number(PeakIslandMaxFlankBp, "PeakIslandMaxFlankBp")
  Validate_Logical(ShowARSLabels, "ShowARSLabels")
  Validate_Logical(ShowFeatureLegend, "ShowFeatureLegend")
  Validate_Logical(SaveCollapsedTable, "SaveCollapsedTable")
  Validate_Positive_Number(WindowSizeKb, "WindowSizeKb")
  Validate_Positive_Number(PanelsPerPage, "PanelsPerPage")
  Validate_Positive_Number(Y_axis_scale, "Y_axis_scale")
  Validate_Positive_Number(ProfileLineLwd, "ProfileLineLwd")
  Validate_Positive_Number(ProfileBarLwd, "ProfileBarLwd")
  Validate_Positive_Number(ORFLwd, "ORFLwd")
  Validate_Positive_Number(ARSLabelCex, "ARSLabelCex")
  Validate_Positive_Number(ARSLabelMinGap, "ARSLabelMinGap")
  Validate_Positive_Number(ARSLabelLevels, "ARSLabelLevels")
  Validate_Nonnegative_Number(RatioPValueCutoff, "RatioPValueCutoff")
  Validate_Positive_Number(PdfWidth, "PdfWidth")
  Validate_Positive_Number(PdfHeight, "PdfHeight")
  Validate_Margin_Vector(ProfileMar, "ProfileMar")
  Validate_Margin_Vector(AnnotationMar, "AnnotationMar")
  Validate_Margin_Vector(OuterMar, "OuterMar")
  Validate_Nonnegative_Number(PanelGapHeight, "PanelGapHeight")
  if(!is.numeric(Log2YMin) || length(Log2YMin) != 1 || !is.finite(Log2YMin)){
    stop("Log2YMin must be one finite number.")
  }

  PlotMode <- match.arg(PlotMode)
  PlotStyle <- match.arg(PlotStyle, choices = c("hist", "bars", "lines"))
  if(PlotStyle == "hist"){
    PlotStyle <- "bars"
  }
  PanelsPerPage <- as.integer(PanelsPerPage)
  ARSLabelLevels <- as.integer(ARSLabelLevels)

  if(PlotMode %in% c("ChIP_Noise", "Clean") && FilterNoise == FALSE){
    stop("PlotMode='", PlotMode, "' requires FilterNoise=TRUE.")
  }

  prepared <- Prepare_WholeChromosome_ChIP_Data(
    SampleDir = SampleDir,
    SampleName = SampleName,
    CoverageDir = CoverageDir,
    CountsFile = CountsFile,
    RatioFile = RatioFile,
    OutputDir = OutputDir,
    CoverageColumn = CoverageColumn,
    FilterNoise = FilterNoise,
    PeakDir = PeakDir,
    PeakFile = PeakFile,
    ProcessPeaks = ProcessPeaks,
    ChIPBam = ChIPBam,
    InputBam = InputBam,
    BedtoolsCommand = BedtoolsCommand,
    Macs2Command = Macs2Command,
    GenomeSize = GenomeSize,
    PeakPValue = PeakPValue,
    NoiseChunkSizeBp = NoiseChunkSizeBp,
    NoiseIterations = NoiseIterations,
    NoiseSeed = NoiseSeed,
    NoiseStatistic = NoiseStatistic,
    NoiseSmoothingSpar = NoiseSmoothingSpar,
    NoiseFloor = NoiseFloor,
    PeakPaddingBp = PeakPaddingBp,
    MinNoiseGapBp = MinNoiseGapBp,
    NoiseChunkMinBp = NoiseChunkMinBp,
    NoiseChunkMaxBp = NoiseChunkMaxBp,
    SaveCollapsedTable = SaveCollapsedTable
  )

  SampleName <- prepared$paths$SampleName
  SampleDir <- prepared$paths$SampleDir
  OutputDir <- Resolve_Output_Directory(OutputDir, SampleDir)

  annotations <- Get_WholeChromosome_Annotations(
    AnnotationDir = AnnotationDir,
    ORFFile = ORFFile,
    ARSFile = ARSFile,
    TERFile = TERFile,
    TyFile = TyFile,
    tRNAFile = tRNAFile,
    CentromereFile = CentromereFile
  )

  ChromosomeInfo <- Get_SacCer3_Chromosome_Table()
  plot_chromosomes <- Resolve_Plot_Chromosomes(Chromosome, ChromosomeInfo)
  full_window_bp <- WindowSizeKb * 1000

  profile <- Profile_Values(prepared$coverage, PlotMode,
                            Log2Profile = Log2Profile,
                            UseRatioLog2 = UseRatioLog2,
                            ApplyRatioPValueSign = ApplyRatioPValueSign,
                            RatioPValueCutoff = RatioPValueCutoff)
  profile_table <- prepared$coverage
  profile_table$plot_signal <- profile$values

  chr_signal <- profile_table$plot_signal[
    profile_table$chrom %in% plot_chromosomes$chrom
  ]
  if(is.null(FixedYlim)){
    if(isTRUE(profile$is_log2)){
      y_max <- Get_Threshold(chr_signal, Y_axis_scale = Y_axis_scale)
      y_max <- max(1, y_max, na.rm = TRUE)
      if(ShowBaseline == TRUE && is.finite(profile$baseline)){
        y_max <- max(y_max, profile$baseline + 1)
      }
      Ylim_sc <- c(Log2YMin, ceiling(y_max))
    } else {
      y_max <- Get_Threshold(chr_signal, Y_axis_scale = Y_axis_scale)
      if(ShowBaseline == TRUE && is.finite(profile$baseline)){
        y_max <- max(y_max, profile$baseline * 1.05)
      }
      Ylim_sc <- c(0, y_max)
    }
  } else {
    Ylim_sc <- FixedYlim
  }

  if(is.null(BaselineValue)){
    BaselineValue <- profile$baseline
  }

  yAxis_reads <- pretty(Ylim_sc, n = 5)
  yAxis_reads <- yAxis_reads[yAxis_reads >= min(Ylim_sc) &
                               yAxis_reads <= max(Ylim_sc)]

  pdf_suffix <- switch(
    PlotMode,
    ChIP = "chip",
    ChIP_Input = "chip_input",
    ChIP_Noise = "chip_noise",
    Clean = "clean"
  )
  if(isTRUE(profile$is_log2)){
    pdf_suffix <- paste0("log2_", pdf_suffix)
  }
  chrom_suffix <- if(nrow(plot_chromosomes) == 1){
    plot_chromosomes$chrom[1]
  } else if(nrow(plot_chromosomes) == min(16, nrow(ChromosomeInfo)) &&
            identical(plot_chromosomes$chrom,
                      ChromosomeInfo$chrom[seq_len(nrow(plot_chromosomes))])){
    "whole_genome"
  } else {
    paste(plot_chromosomes$chrom, collapse = "_")
  }
  out_file <- file.path(
    OutputDir,
    paste0(SampleName, "_", chrom_suffix, "_", pdf_suffix, ".pdf")
  )

  grDevices::pdf(file = out_file, width = PdfWidth, height = PdfHeight)
  on.exit(grDevices::dev.off(), add = TRUE)

  for(chr_i in seq_len(nrow(plot_chromosomes))){
    chr_info <- plot_chromosomes[chr_i, , drop = FALSE]
    windows <- Build_WholeChromosome_Windows(
      chr_info$chrom,
      WindowSizeKb = WindowSizeKb,
      ChromosomeInfo = ChromosomeInfo
    )

    n_pages <- ceiling(nrow(windows) / PanelsPerPage)
    for(page in seq_len(n_pages)){
      page_start <- (page - 1) * PanelsPerPage + 1
      page_end <- min(page * PanelsPerPage, nrow(windows))
      page_windows <- windows[page_start:page_end, , drop = FALSE]
      page_windows_n <- nrow(page_windows)
      page_slots <- PanelsPerPage

      layout(matrix(seq_len(page_slots * 3), ncol = 1),
             heights = rep(c(3.2, 1.35, max(PanelGapHeight, 0.001)),
                           page_slots))
      par(oma = OuterMar)

      for(i in seq_len(page_slots)){
        if(i > page_windows_n){
          Plot_WholeChromosome_Blank_Window(
            ProfileMar = ProfileMar,
            AnnotationMar = AnnotationMar,
            ShowLegend = ShowFeatureLegend && i == page_slots,
            PeakColor = PeakColor,
            PlotPeaks = PlotPeaks
          )
          next
        }

        w <- page_windows[i, ]
        title_text <- paste0(
          SampleName, " ", w$chrom, ":",
          format(w$S, scientific = FALSE, trim = TRUE), "-",
          format(w$E, scientific = FALSE, trim = TRUE)
        )

        Plot_WholeChromosome_Profile_Window(
          ProfileTable = profile_table,
          annotations = annotations,
          peaks = prepared$peaks,
          window = w,
          steps = 10,
          FullWindowBp = full_window_bp,
          Ylim_sc = Ylim_sc,
          yAxis_reads = yAxis_reads,
          PlotStyle = PlotStyle,
          SmoothSignal = SmoothSignal,
          SmoothingSpar = SmoothingSpar,
          ShowBaseline = ShowBaseline,
          BaselineValue = BaselineValue,
          SignalColor = SignalColor,
          PositiveSignalColor = PositiveSignalColor,
          NegativeSignalColor = NegativeSignalColor,
          ProfileBarLwd = ProfileBarLwd,
          ProfileLineLwd = ProfileLineLwd,
          RedRequiresPeakOverlap = RedRequiresPeakOverlap,
          PeakIslandBoundaryFraction = PeakIslandBoundaryFraction,
          PeakIslandMaxFlankBp = PeakIslandMaxFlankBp,
          PlotPeaks = PlotPeaks,
          PeakColor = PeakColor,
          ShowARSLabels = ShowARSLabels,
          ARSLabelCex = ARSLabelCex,
          ARSLabelMinGap = ARSLabelMinGap,
          ARSLabelLevels = ARSLabelLevels,
          WindowTitle = title_text,
          YLabel = profile$label,
          ProfileMar = ProfileMar
        )
        Plot_WholeChromosome_Annotation_Window(
          annotations = annotations,
          window = w,
          steps = 10,
          FullWindowBp = full_window_bp,
          ORFLwd = ORFLwd,
          ShowXAxis = TRUE,
          ShowXLabel = i == page_windows_n,
          FeatureTracks = FeatureTracks,
          AnnotationMar = AnnotationMar
        )
        Plot_WholeChromosome_Spacer(
          ShowLegend = ShowFeatureLegend && i == page_slots,
          PeakColor = PeakColor,
          ARSColor = "purple",
          ARSFill = "yellow",
          PlotPeaks = PlotPeaks
        )
      }

      title_prefix <- if(nrow(plot_chromosomes) == 1){
        "Whole chromosome "
      } else {
        "Whole genome "
      }
      mtext(paste0(title_prefix, profile$label, " | ", chr_info$chrom,
                   " | page ", page, "/", n_pages),
            outer = TRUE, side = 3, line = 0.5, cex = 0.95, col = "gray35")
    }
  }

  invisible(list(
    pdf = out_file,
    data_file = prepared$data_file,
    peak_file = prepared$peak_file,
    counts_file = prepared$counts_file,
    ratio_file = prepared$ratio_file,
    chromosomes = plot_chromosomes$chrom,
    noise_parameters = prepared$noise_parameters
  ))
}

WholeGenome_ChIP_Plotter <- WholeChromosome_ChIP_Plotter

WholeGenome_Plotter <- function(SampleDir,
                                Chromosomes = "all",
                                PlotMode = c("ChIP_Input", "ChIP",
                                             "ChIP_Noise", "Clean"),
                                SampleName = NULL,
                                OutputDir = "SampleDir",
                                CountsFile = NULL,
                                RatioFile = NULL,
                                PeakFile = NULL,
                                WindowSizeKb = 50,
                                PlotStyle = "hist",
                                ProcessPeaks = "auto",
                                FilterNoise = "auto",
                                ...){
  if(missing(SampleDir) ||
     !nzchar(as.character(SampleDir)) ||
     identical(as.character(SampleDir), "EDIT_SAMPLE_DIRECTORY_HERE")){
    stop("Edit SampleDir in the run script to the sample folder path.")
  }

  PlotMode <- match.arg(PlotMode)

  if(identical(FilterNoise, "auto")){
    FilterNoise <- PlotMode %in% c("ChIP_Noise", "Clean")
  } else {
    Validate_Logical(FilterNoise, "FilterNoise")
  }

  if(identical(ProcessPeaks, "auto")){
    ProcessPeaks <- FilterNoise
  } else {
    Validate_Logical(ProcessPeaks, "ProcessPeaks")
  }

  result <- WholeGenome_ChIP_Plotter(
    SampleDir = SampleDir,
    SampleName = SampleName,
    OutputDir = OutputDir,
    CountsFile = CountsFile,
    RatioFile = RatioFile,
    PeakFile = PeakFile,
    Chromosome = Chromosomes,
    PlotMode = PlotMode,
    WindowSizeKb = WindowSizeKb,
    PlotStyle = PlotStyle,
    FilterNoise = FilterNoise,
    ProcessPeaks = ProcessPeaks,
    ...
  )

  message("PDF: ", result$pdf)
  message("RATIO: ", result$ratio_file)
  message("COUNTS: ", result$counts_file)
  message("DATA: ", result$data_file)
  message("PEAKS: ", result$peak_file)

  invisible(result)
}

WholeChromosome_Plotter <- function(SampleDir,
                                    Chromosome = "chrI",
                                    PlotMode = c("ChIP_Input", "ChIP",
                                                 "ChIP_Noise", "Clean"),
                                    ...){
  WholeGenome_Plotter(
    SampleDir = SampleDir,
    Chromosomes = Chromosome,
    PlotMode = PlotMode,
    ...
  )
}
