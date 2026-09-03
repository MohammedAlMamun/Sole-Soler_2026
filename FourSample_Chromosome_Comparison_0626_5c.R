source("/Volumes/rbmData/Smc5_project/whole_chr_map/ChIP_Chr_Ann_Map_Script/WholeChromosome_Plotter_0626_5c.R")

FourSample_Expand_Arg <- function(x, arg_name, n = 4){
  if(is.null(x)){
    return(vector("list", n))
  }
  if(is.list(x)){
    if(length(x) == 1){
      return(rep(x, n))
    }
    if(length(x) != n){
      stop(arg_name, " must be NULL, length 1, or length ", n, ".")
    }
    return(x)
  }
  if(length(x) == 1){
    return(as.list(rep(x, n)))
  }
  if(length(x) != n){
    stop(arg_name, " must be NULL, length 1, or length ", n, ".")
  }
  as.list(x)
}

FourSample_Null_If_NA <- function(x){
  if(is.null(x) || length(x) == 0){
    return(NULL)
  }
  if(length(x) == 1 && is.na(x)){
    return(NULL)
  }
  x
}

FourSample_Safe_File_Name <- function(x){
  x <- as.character(x)
  x[!nzchar(x) | is.na(x)] <- "sample"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

FourSample_Output_Dir <- function(OutputDir, SampleDirs){
  first_sample_dir <- normalizePath(path.expand(SampleDirs[1]), mustWork = TRUE)
  common_parent <- dirname(first_sample_dir)

  use_first_sample <- is.null(OutputDir) ||
    length(OutputDir) == 0 ||
    (length(OutputDir) == 1 &&
       (is.na(OutputDir) ||
          OutputDir == "" ||
          tolower(OutputDir) %in% c("sampledir", "sample_dir",
                                    "sample", "firstsampledir",
                                    "first_sample_dir")))
  use_common_parent <- length(OutputDir) == 1 &&
    !is.na(OutputDir) &&
    tolower(OutputDir) %in% c("commonparent", "common_parent",
                              "parent", "whole_chr_map")

  if(use_first_sample){
    OutputDir <- first_sample_dir
  } else if(use_common_parent){
    OutputDir <- common_parent
  } else if(length(OutputDir) != 1){
    stop("OutputDir must be NULL, 'SampleDir', 'CommonParent', or one folder path.")
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

FourSample_PlotMode_Suffix <- function(PlotMode, IsLog2 = TRUE){
  suffix <- switch(
    PlotMode,
    ChIP = "chip",
    ChIP_Input = "chip_input",
    ChIP_Noise = "chip_noise",
    Clean = "clean"
  )
  if(IsLog2 == TRUE){
    suffix <- paste0("log2_", suffix)
  }
  suffix
}

FourSample_Blank_Gap <- function(){
  par(mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", bty = "n", xaxs = "i", yaxs = "i")
}

FourSample_Prepare_One <- function(i, SampleDirs, SampleNames, CoverageDirs,
                                   CountsFiles, RatioFiles, PeakFiles,
                                   OutputDir, PlotMode, Log2Profile,
                                   UseRatioLog2, ApplyRatioPValueSign,
                                   RatioPValueCutoff, CoverageColumn,
                                   FilterNoise, ProcessPeaks, ChIPBams,
                                   InputBams, BedtoolsCommand, Macs2Command,
                                   GenomeSize, PeakPValue, NoiseChunkSizeBp,
                                   NoiseIterations, NoiseSeed, NoiseStatistic,
                                   NoiseSmoothingSpar, NoiseFloor,
                                   PeakPaddingBp, MinNoiseGapBp,
                                   NoiseChunkMinBp, NoiseChunkMaxBp,
                                   SaveCollapsedTable){
  prepared <- Prepare_WholeChromosome_ChIP_Data(
    SampleDir = SampleDirs[[i]],
    SampleName = FourSample_Null_If_NA(SampleNames[[i]]),
    CoverageDir = FourSample_Null_If_NA(CoverageDirs[[i]]),
    CountsFile = FourSample_Null_If_NA(CountsFiles[[i]]),
    RatioFile = FourSample_Null_If_NA(RatioFiles[[i]]),
    OutputDir = OutputDir,
    CoverageColumn = CoverageColumn,
    FilterNoise = FilterNoise,
    PeakFile = FourSample_Null_If_NA(PeakFiles[[i]]),
    ProcessPeaks = ProcessPeaks,
    ChIPBam = FourSample_Null_If_NA(ChIPBams[[i]]),
    InputBam = FourSample_Null_If_NA(InputBams[[i]]),
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

  profile <- Profile_Values(
    prepared$coverage,
    PlotMode,
    Log2Profile = Log2Profile,
    UseRatioLog2 = UseRatioLog2,
    ApplyRatioPValueSign = ApplyRatioPValueSign,
    RatioPValueCutoff = RatioPValueCutoff
  )
  profile_table <- prepared$coverage
  profile_table$plot_signal <- profile$values

  list(
    label = prepared$paths$SampleName,
    prepared = prepared,
    profile = profile,
    profile_table = profile_table
  )
}

FourSample_Make_Ylim <- function(signal_values, profile, Log2YMin,
                                 Y_axis_scale, ShowBaseline, FixedYlim){
  if(!is.null(FixedYlim)){
    return(FixedYlim)
  }

  if(isTRUE(profile$is_log2)){
    y_max <- Get_Threshold(signal_values, Y_axis_scale = Y_axis_scale)
    y_max <- max(1, y_max, na.rm = TRUE)
    if(ShowBaseline == TRUE && is.finite(profile$baseline)){
      y_max <- max(y_max, profile$baseline + 1)
    }
    return(c(Log2YMin, ceiling(y_max)))
  }

  y_max <- Get_Threshold(signal_values, Y_axis_scale = Y_axis_scale)
  if(ShowBaseline == TRUE && is.finite(profile$baseline)){
    y_max <- max(y_max, profile$baseline * 1.05)
  }
  c(0, y_max)
}

FourSample_Chromosome_Comparison_Plotter <- function(
    SampleDirs,
    Chromosome = "chrI",
    StartCoordinate = NULL,
    EndCoordinate = NULL,
    PlotMode = c("ChIP_Input", "ChIP", "ChIP_Noise", "Clean"),
    SampleNames = NULL,
    OutputDir = "CommonParent",
    ComparisonName = NULL,
    CoverageDirs = NULL,
    CountsFiles = NULL,
    RatioFiles = "auto",
    PeakFiles = NULL,
    AnnotationDir = NULL,
    ORFFile = NULL,
    ARSFile = NULL,
    TERFile = NULL,
    TyFile = NULL,
    tRNAFile = NULL,
    CentromereFile = NULL,
    ProcessPeaks = "auto",
    FilterNoise = "auto",
    ChIPBams = NULL,
    InputBams = NULL,
    BedtoolsCommand = NULL,
    Macs2Command = "macs2",
    GenomeSize = 12157105,
    PeakPValue = "10e-6",
    CoverageColumn = 5,
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
    WindowSizeKb = 50,
    PlotStyle = "hist",
    Log2Profile = TRUE,
    UseRatioLog2 = TRUE,
    ApplyRatioPValueSign = TRUE,
    RatioPValueCutoff = 10e-2,
    SharedYScale = TRUE,
    Log2YMin = -1,
    FixedYlim = NULL,
    Y_axis_scale = 8,
    SmoothSignal = TRUE,
    SmoothingSpar = NULL,
    ShowBaseline = TRUE,
    BaselineValue = NULL,
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
    ARSLabelCex = 0.62,
    ARSLabelMinGap = 350,
    ARSLabelLevels = 4,
    FeatureTracks = c("ORF", "Ty", "TER", "tRNA"),
    SaveCollapsedTable = FALSE,
    PdfWidth = 11,
    PdfHeight = 11,
    TopProfileMar = c(0.30, 4.3, 2.15, 1.0),
    ProfileMar = c(0.30, 4.3, 1.15, 1.0),
    AnnotationMar = c(1.65, 4.3, 0.35, 1.0),
    OuterMar = c(1.4, 1.2, 1.7, 0.8),
    PanelGapHeight = 0.30){
  if(missing(SampleDirs) || length(SampleDirs) != 4){
    stop("SampleDirs must contain exactly four sample folders.")
  }
  if(length(Chromosome) != 1){
    stop("Chromosome must be one chromosome for this comparison plotter.")
  }
  coordinate_mode <- !is.null(StartCoordinate) || !is.null(EndCoordinate)
  if(coordinate_mode){
    if(is.null(StartCoordinate) || is.null(EndCoordinate)){
      stop("StartCoordinate and EndCoordinate must either both be NULL or both be supplied.")
    }
    if(!is.numeric(StartCoordinate) || length(StartCoordinate) != 1 ||
       !is.finite(StartCoordinate) || StartCoordinate < 0){
      stop("StartCoordinate must be one finite, non-negative number.")
    }
    if(!is.numeric(EndCoordinate) || length(EndCoordinate) != 1 ||
       !is.finite(EndCoordinate) || EndCoordinate <= StartCoordinate){
      stop("EndCoordinate must be one finite number greater than StartCoordinate.")
    }
    if(StartCoordinate != floor(StartCoordinate) ||
       EndCoordinate != floor(EndCoordinate)){
      stop("StartCoordinate and EndCoordinate must be whole-base coordinates.")
    }
  }

  PlotMode <- match.arg(PlotMode)
  PlotStyle <- match.arg(PlotStyle, choices = c("hist", "bars", "lines"))
  if(PlotStyle == "hist"){
    PlotStyle <- "bars"
  }

  Validate_Logical(Log2Profile, "Log2Profile")
  Validate_Logical(UseRatioLog2, "UseRatioLog2")
  Validate_Logical(ApplyRatioPValueSign, "ApplyRatioPValueSign")
  Validate_Logical(SharedYScale, "SharedYScale")
  Validate_Logical(SmoothSignal, "SmoothSignal")
  Validate_Logical(ShowBaseline, "ShowBaseline")
  Validate_Logical(PlotPeaks, "PlotPeaks")
  Validate_Logical(RedRequiresPeakOverlap, "RedRequiresPeakOverlap")
  Validate_Fraction(PeakIslandBoundaryFraction, "PeakIslandBoundaryFraction")
  Validate_Nonnegative_Number(PeakIslandMaxFlankBp, "PeakIslandMaxFlankBp")
  Validate_Logical(ShowARSLabels, "ShowARSLabels")
  Validate_Logical(SaveCollapsedTable, "SaveCollapsedTable")
  Validate_Positive_Number(WindowSizeKb, "WindowSizeKb")
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
  Validate_Margin_Vector(TopProfileMar, "TopProfileMar")
  Validate_Margin_Vector(ProfileMar, "ProfileMar")
  Validate_Margin_Vector(AnnotationMar, "AnnotationMar")
  Validate_Margin_Vector(OuterMar, "OuterMar")
  Validate_Nonnegative_Number(PanelGapHeight, "PanelGapHeight")
  if(!is.numeric(Log2YMin) || length(Log2YMin) != 1 || !is.finite(Log2YMin)){
    stop("Log2YMin must be one finite number.")
  }
  if(!is.null(FixedYlim) &&
     (!is.numeric(FixedYlim) || length(FixedYlim) != 2 ||
        any(!is.finite(FixedYlim)) || FixedYlim[2] <= FixedYlim[1])){
    stop("FixedYlim must be NULL or a numeric vector c(ymin, ymax).")
  }

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
  if(PlotMode %in% c("ChIP_Noise", "Clean") && FilterNoise == FALSE){
    stop("PlotMode='", PlotMode, "' requires FilterNoise=TRUE.")
  }

  SampleDirs <- as.list(normalizePath(path.expand(as.character(SampleDirs)),
                                      mustWork = TRUE))
  SampleNames <- FourSample_Expand_Arg(SampleNames, "SampleNames")
  CoverageDirs <- FourSample_Expand_Arg(CoverageDirs, "CoverageDirs")
  CountsFiles <- FourSample_Expand_Arg(CountsFiles, "CountsFiles")
  RatioFiles <- FourSample_Expand_Arg(RatioFiles, "RatioFiles")
  PeakFiles <- FourSample_Expand_Arg(PeakFiles, "PeakFiles")
  ChIPBams <- FourSample_Expand_Arg(ChIPBams, "ChIPBams")
  InputBams <- FourSample_Expand_Arg(InputBams, "InputBams")

  OutputDir <- FourSample_Output_Dir(OutputDir, unlist(SampleDirs))

  sample_results <- vector("list", 4)
  for(i in seq_len(4)){
    sample_results[[i]] <- FourSample_Prepare_One(
      i = i,
      SampleDirs = SampleDirs,
      SampleNames = SampleNames,
      CoverageDirs = CoverageDirs,
      CountsFiles = CountsFiles,
      RatioFiles = RatioFiles,
      PeakFiles = PeakFiles,
      OutputDir = OutputDir,
      PlotMode = PlotMode,
      Log2Profile = Log2Profile,
      UseRatioLog2 = UseRatioLog2,
      ApplyRatioPValueSign = ApplyRatioPValueSign,
      RatioPValueCutoff = RatioPValueCutoff,
      CoverageColumn = CoverageColumn,
      FilterNoise = FilterNoise,
      ProcessPeaks = ProcessPeaks,
      ChIPBams = ChIPBams,
      InputBams = InputBams,
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
  }

  sample_labels <- vapply(sample_results, function(x) x$label, character(1))
  if(is.null(ComparisonName)){
    ComparisonName <- paste(FourSample_Safe_File_Name(sample_labels), collapse = "_")
  }

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
  chr_info <- Resolve_Plot_Chromosomes(Chromosome, ChromosomeInfo)
  if(nrow(chr_info) != 1){
    stop("Chromosome must resolve to exactly one chromosome.")
  }
  if(coordinate_mode){
    if(EndCoordinate > chr_info$length){
      stop("EndCoordinate (", EndCoordinate, ") exceeds the length of ",
           chr_info$chrom, " (", chr_info$length, " bp).")
    }
    windows <- data.frame(
      k = chr_info$k,
      chrom = chr_info$chrom,
      S = StartCoordinate,
      E = EndCoordinate,
      stringsAsFactors = FALSE
    )
    full_window_bp <- EndCoordinate - StartCoordinate
  } else {
    windows <- Build_WholeChromosome_Windows(
      chr_info$chrom,
      WindowSizeKb = WindowSizeKb,
      ChromosomeInfo = ChromosomeInfo
    )
    full_window_bp <- WindowSizeKb * 1000
  }

  signal_in_plot_range <- function(profile_table){
    profile_table$chrom == chr_info$chrom &
      profile_table$chromStart >= windows$S[1] &
      profile_table$chromStart < windows$E[nrow(windows)]
  }

  profile_labels <- unique(vapply(sample_results, function(x) x$profile$label,
                                  character(1)))
  profile_label <- paste(profile_labels, collapse = " / ")
  is_log2 <- all(vapply(sample_results, function(x) isTRUE(x$profile$is_log2),
                        logical(1)))
  if(is.null(BaselineValue)){
    BaselineValue <- sample_results[[1]]$profile$baseline
  }

  if(SharedYScale == TRUE){
    all_signal <- unlist(lapply(sample_results, function(x){
      x$profile_table$plot_signal[signal_in_plot_range(x$profile_table)]
    }), use.names = FALSE)
    shared_ylim <- FourSample_Make_Ylim(
      all_signal,
      sample_results[[1]]$profile,
      Log2YMin = Log2YMin,
      Y_axis_scale = Y_axis_scale,
      ShowBaseline = ShowBaseline,
      FixedYlim = FixedYlim
    )
    ylims <- rep(list(shared_ylim), 4)
  } else {
    ylims <- lapply(sample_results, function(x){
      signal <- x$profile_table$plot_signal[signal_in_plot_range(x$profile_table)]
      FourSample_Make_Ylim(
        signal,
        x$profile,
        Log2YMin = Log2YMin,
        Y_axis_scale = Y_axis_scale,
        ShowBaseline = ShowBaseline,
        FixedYlim = FixedYlim
      )
    })
  }

  y_axes <- lapply(ylims, function(ylim_now){
    axis_vals <- pretty(ylim_now, n = 4)
    axis_vals[axis_vals >= min(ylim_now) & axis_vals <= max(ylim_now)]
  })

  pdf_suffix <- FourSample_PlotMode_Suffix(PlotMode, IsLog2 = is_log2)
  region_suffix <- if(coordinate_mode){
    paste0(
      format(StartCoordinate, scientific = FALSE, trim = TRUE), "-",
      format(EndCoordinate, scientific = FALSE, trim = TRUE), "bp"
    )
  } else {
    paste0(WindowSizeKb, "kb")
  }
  out_file <- file.path(
    OutputDir,
    paste0(FourSample_Safe_File_Name(ComparisonName), "_",
           chr_info$chrom, "_", region_suffix, "_", pdf_suffix, ".pdf")
  )

  grDevices::pdf(file = out_file, width = PdfWidth, height = PdfHeight)
  on.exit(grDevices::dev.off(), add = TRUE)

  for(page in seq_len(nrow(windows))){
    w <- windows[page, ]
    layout(matrix(seq_len(9), ncol = 1),
           heights = c(2.35, max(PanelGapHeight, 0.001),
                       2.05, max(PanelGapHeight, 0.001),
                       2.05, max(PanelGapHeight, 0.001),
                       2.05, max(PanelGapHeight, 0.001),
                       1.35))
    par(oma = OuterMar)

    for(i in seq_len(4)){
      panel_title <- paste0(
        sample_labels[i], " ", w$chrom, ":",
        format(w$S, scientific = FALSE, trim = TRUE), "-",
        format(w$E, scientific = FALSE, trim = TRUE)
      )

      Plot_WholeChromosome_Profile_Window(
        ProfileTable = sample_results[[i]]$profile_table,
        annotations = annotations,
        peaks = sample_results[[i]]$prepared$peaks,
        window = w,
        steps = 10,
        FullWindowBp = full_window_bp,
        Ylim_sc = ylims[[i]],
        yAxis_reads = y_axes[[i]],
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
        ShowARSLabels = ShowARSLabels && i == 1,
        ARSLabelCex = ARSLabelCex,
        ARSLabelMinGap = ARSLabelMinGap,
        ARSLabelLevels = ARSLabelLevels,
        WindowTitle = panel_title,
        YLabel = if(i == 2) profile_label else "",
        ProfileMar = if(i == 1) TopProfileMar else ProfileMar
      )
      FourSample_Blank_Gap()
    }

    Plot_WholeChromosome_Annotation_Window(
      annotations = annotations,
      window = w,
      steps = 10,
      FullWindowBp = full_window_bp,
      ORFLwd = ORFLwd,
      ShowXAxis = TRUE,
      ShowXLabel = TRUE,
      FeatureTracks = FeatureTracks,
      AnnotationMar = AnnotationMar
    )

    title_text <- paste0(
      ComparisonName, " | ", chr_info$chrom, ":",
      format(round(w$S / 1000, 1), scientific = FALSE, trim = TRUE), "-",
      format(round(w$E / 1000, 1), scientific = FALSE, trim = TRUE),
      " kbp",
      if(!coordinate_mode) paste0(" | page ", page, "/", nrow(windows)) else ""
    )
    mtext(title_text, outer = TRUE, side = 3, line = 0.35,
          cex = 0.92, col = "gray35")
  }

  message("PDF: ", out_file)
  message("RATIO: ",
          paste(vapply(sample_results, function(x) x$prepared$ratio_file,
                       character(1)), collapse = " | "))
  message("COUNTS: ",
          paste(vapply(sample_results, function(x) x$prepared$counts_file,
                       character(1)), collapse = " | "))
  message("PEAKS: ",
          paste(vapply(sample_results, function(x) x$prepared$peak_file,
                       character(1)), collapse = " | "))

  invisible(list(
    pdf = out_file,
    samples = sample_labels,
    chromosome = chr_info$chrom,
    start_coordinate = if(coordinate_mode) StartCoordinate else NULL,
    end_coordinate = if(coordinate_mode) EndCoordinate else NULL,
    ratio_files = vapply(sample_results, function(x) x$prepared$ratio_file,
                         character(1)),
    counts_files = vapply(sample_results, function(x) x$prepared$counts_file,
                          character(1)),
    peak_files = vapply(sample_results, function(x) x$prepared$peak_file,
                        character(1))
  ))
}

FourSample_Plotter <- FourSample_Chromosome_Comparison_Plotter

## Edit only this block for routine four-sample comparison plotting.
## Set Run_FourSample_Comparison <- FALSE before sourcing if you only want functions.
if(!exists("Run_FourSample_Comparison", inherits = FALSE)){
  Run_FourSample_Comparison <- TRUE
}

SampleDirs <- c(
  "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-WT-July23",
  "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-DA-July23",
  "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-EQ-July23",
  "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-noTag-July23"
)
SampleNames <- NULL
Chromosome <- "chrI"
## Leave both as NULL for the current full-chromosome, multi-page plot.
## Set both (0-based start, end-exclusive) for one exact region on one PDF page.
StartCoordinate <- NULL
EndCoordinate <- NULL
PlotMode <- "ChIP_Input"
RatioFiles <- "auto"
PeakFiles <- NULL
WindowSizeKb <- 50
OutputDir <- "CommonParent"
ComparisonName <- "Smc5_July23_4sample"
SaveDataTable <- FALSE

if(isTRUE(Run_FourSample_Comparison)){
  FourSample_Chromosome_Comparison_Plotter(
    SampleDirs = SampleDirs,
    SampleNames = SampleNames,
    Chromosome = Chromosome,
    StartCoordinate = StartCoordinate,
    EndCoordinate = EndCoordinate,
    PlotMode = PlotMode,
    RatioFiles = RatioFiles,
    PeakFiles = PeakFiles,
    WindowSizeKb = WindowSizeKb,
    OutputDir = OutputDir,
    ComparisonName = ComparisonName,
    SaveCollapsedTable = SaveDataTable
  )
}
