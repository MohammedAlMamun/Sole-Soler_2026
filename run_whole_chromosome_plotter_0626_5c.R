

source("/Volumes/rbmData/Smc5_project/whole_chr_map/ChIP_Chr_Ann_Map_Script/WholeChromosome_Plotter_0626_5c.R")

## Keep WholeChromosome_Plotter_0626_support beside the main script.
## Edit only this block for routine plotting.
## RatioFile = "auto" uses <SampleName>_Ratio_pro.bed / *_Ratio_pro.bed
## from SampleDir. This exactly reuses old log2.ratio profiles when present.
## Set RatioFile <- NULL to plot from *_normal_counts.bed or strand coverage files.
## Leave PeakFile NULL to auto-use *_Peaks.bed in SampleDir.
## Use "all" for chrI-chrXVI. For testing one chromosome, use "chrI", "chrII", etc.
## PlotMode choices: "ChIP", "ChIP_Input", "ChIP_Noise", "Clean".

SampleDir <- "/Volumes/rbmData/Smc5_project/whole_chr_map/Smc5-noTag-July23"
SampleName <- NULL
RatioFile <- "auto"
PeakFile <- NULL
Chromosomes <- "all"
PlotMode <- "Clean"
WindowSizeKb <- 50
PanelsPerPage <- 3
OutputDir <- "SampleDir"
SaveDataTable <- FALSE

WholeGenome_Plotter(
  SampleDir = SampleDir,
  SampleName = SampleName,
  RatioFile = RatioFile,
  PeakFile = PeakFile,
  Chromosomes = Chromosomes,
  PlotMode = PlotMode,
  WindowSizeKb = WindowSizeKb,
  PanelsPerPage = PanelsPerPage,
  OutputDir = OutputDir,
  SaveCollapsedTable = SaveDataTable
)
