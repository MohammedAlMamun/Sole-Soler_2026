Run_FourSample_Comparison <- FALSE
source("/Volumes/rbmData/Smc5_project/whole_chr_map/ChIP_Chr_Ann_Map_Script/FourSample_Chromosome_Comparison_0626_5c.R")

## Edit only this block for routine four-sample chromosome comparisons.
## RatioFiles = "auto" uses *_Ratio_pro.bed / *_Ratio.bed from each SampleDir.
## Leave PeakFiles NULL to auto-use *_Peaks.bed from each SampleDir.
## PlotMode choices: "ChIP_Input", "ChIP", "ChIP_Noise", "Clean".
## OutputDir = "CommonParent" writes into the shared parent folder.

SampleDirs <- c(
  "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-WT-July23",
  "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-DA-July23",
  "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-EQ-July23",
  "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-noTag-July23"
)

SampleNames <- NULL
Chromosome <- "chrVI"
StartCoordinate <- NULL
EndCoordinate <- NULL
PlotMode <- "Clean"
RatioFiles <- "auto"
PeakFiles <- NULL
WindowSizeKb <- 50
OutputDir <- "CommonParent"
ComparisonName <- "Smc5_July23_4sample"
SaveDataTable <- FALSE

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
