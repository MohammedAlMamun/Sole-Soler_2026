

source("/Volumes/rbmData/Smc5_project/Smc5_0826.R")

## ChIP-seq primary analysis
## Alignment = "generic", "malign" or "mrdna"
## generic is our usual unique proper pairs, malign is multiple alignments genomewide and mrdna is multiple alignment at rDNA
## Directory e.g. "/Volumes/rbmData/Smc5_project/August_26" will create the result folder in that directory, not in desktop like before
  
# 
# ChIPseq_Primary_Analysis( Input_R1 = "/Volumes/rbmData_2/17_March_2025/fastq/fastq/61_S61_R1_001.fastq.gz",
#                           Input_R2 = "/Volumes/rbmData_2/17_March_2025/fastq/fastq/61_S61_R2_001.fastq.gz",
#                           ChIP_R1 = "/Volumes/rbmData_2/17_March_2025/fastq/fastq/63_S63_R1_001.fastq.gz",
#                           ChIP_R2 = "/Volumes/rbmData_2/17_March_2025/fastq/fastq/63_S63_R2_001.fastq.gz",
#                           ExpTitle = "Rpb3PK-Sen1deg-EtOH-ChIP",
#                           Alignment = "generic",
#                           Directory = "/Volumes/rbmData_2/CD_HO_26" )
# 
# 
# BrDUseq_Primary_Analysis( Input_R1 = "/Volumes/rbmData_2/17_March_2025/fastq/fastq/65_S65_R1_001.fastq.gz",
#                           Input_R2 = "/Volumes/rbmData_2/17_March_2025/fastq/fastq/65_S65_R2_001.fastq.gz",
#                           BrDU_R1 = "/Volumes/rbmData_2/17_March_2025/fastq/fastq/66_S66_R1_001.fastq.gz",
#                           BrDU_R2 = "/Volumes/rbmData_2/17_March_2025/fastq/fastq/66_S66_R2_001.fastq.gz",
#                           ExpTitle = "Rpb3PK-Sen1deg-IAA-BrDU",
#                           Alignment = "malign",
#                           Directory = "/Volumes/rbmData_2/CD_HO_26" )


## ChIP-seq average profiles at early and late origins
## Alignment = "generic" or "malign"
## StrandMode = "collapsed" or "separated"

  ChIPseq_Analysis_early_late(ExpTitle = "Smc5-EQ-60HU-2025-ChIP", 
                              Window = 12000,
                              Alignment = "malign", 
                              StrandMode = "collapsed",
                              Directory = "/Volumes/rbmData_2/SMC5_2026/August_26/2025" )


## BrDU-seq average profiles at early and late origins
## Alignment = "generic" or "malign"
## StrandMode = "collapsed" or "separated"

  BrDUseq_Analysis_early_late(ExpTitle = "Smc5-EQ-60HU-2025-BrDU", 
                              Window = 12000,
                              Alignment = "malign", 
                              StrandMode = "collapsed",
                              Directory = "/Volumes/rbmData_2/SMC5_2026/August_26/2025" )
  

## dual four panel final figures 
  
  ChIP_BrDU_Analysis_early_late(Directory = "/Volumes/rbmData_2/SMC5_2026/August_26/2025",
                                Genotype = "WT",
                                Window = 50000,
                                Alignment = "malign", 
                                yLim = "per-plot",
                                Metric = "ratio.ipin.noise")
  
  ChIP_BrDU_Analysis_early_late(Directory = "/Volumes/rbmData_2/SMC5_2026/August_26/2025",
                                Genotype = "EQ",
                                Window = 50000,
                                Alignment = "malign", 
                                yLim = "per-plot",
                                Metric = "ratio.ipin.noise")
  
  ChIP_BrDU_Analysis_early_late_old_approach(Directory = "/Volumes/rbmData_2/SMC5_2026/August_26/2025",
                                             Genotype = "WT", Window = 50000,
                                             Alignment = "generic", yLim = "per-plot")
  
  ChIP_BrDU_Analysis_early_late_old_approach(Directory = "/Volumes/rbmData_2/SMC5_2026/August_26/2025",
                                             Genotype = "EQ", Window = 50000,
                                             Alignment = "generic", yLim = "per-plot")
  
## Pairwise ChIP enrichment by distance from origins
## origin_file = "E_Ori", "L_Ori" 
## Alignment = "generic" or "malign"
## Signal = "chip_input" or "clean"

  origin_distance_pair_analysis(sample_1_folder = "/Volumes/rbmData_2/SMC5_2026/August_26/2026/Smc5-WT-75HU-2026-ChIP",
                                sample_2_folder = "/Volumes/rbmData_2/SMC5_2026/August_26/2026/Smc5-EQ-75HU-2026-ChIP",
                                sample_1_label = "WT", sample_2_label = "EQ",
                                origin_file = "E_Ori", Alignment = "malign", Signal = "clean",
                                distance_breaks_kb = c(0, 0.75, 1.5, 3, 4.5, 6) )
  
  origin_distance_pair_analysis_old_approach( sample_1_folder = "/Volumes/rbmData_2/SMC5_2026/August_26/2026/Smc5-WT-75HU-2026-ChIP",
                                              sample_2_folder = "/Volumes/rbmData_2/SMC5_2026/August_26/2026/Smc5-EQ-75HU-2026-ChIP",
                                              sample_1_label = "WT", sample_2_label = "EQ",
                                              origin_file = "E_Ori", Alignment = "generic",
                                              yLim = "fixed", distance_breaks_kb = c(0, 0.75, 1.5, 3, 4.5, 6) )
  
  