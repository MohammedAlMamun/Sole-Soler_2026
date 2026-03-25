
### ------------ Run the alignmemnt with multiple mapping reads with Rsubread ----------
### ----- mask and unmask the reference index appropriately for rDNA or whole genome ---

SubRead_Multiple_Aligner_rDNA <- function(Input_R1, 
                                          Input_R2, 
                                          ChIP_R1, 
                                          ChIP_R2,
                                          binSize = 300, 
                                          stepSize = 10, 
                                          Sample_Name, 
                                          AlignmentMode = "multimapping"){
  
  
  library(Rsubread)
  library(GenomicAlignments)
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  library(rasterpdf)
  library(plotrix)
  library(ORFik)
  library(IRanges)
  library(MACSr)
  library(readxl)
  
  #create a folder by the sample name in the desktop
  suppressWarnings(dir.create(paste0("~/Desktop/", Sample_Name)))
  
  #create a folder named BAM to save the bam files
  suppressWarnings(dir.create(paste0("~/Desktop/", Sample_Name, "/", "Bam")))
  
  runAlignment <- {
    
    ###read prebuilt indexed reference genome
    
    ##for rDNA as reference
    indexed_reference = "/Volumes/rbmData/Smc5_project/rDNA_2nts1"
    
    ##for whole genome as reference
    #indexed_reference = "/Users/mohammed/Desktop/ngsAnalyser1.1.4/bowtie2-2.4.4-macos-x86_64/indexes/subread_ref"
    
    inputBam <- paste0("~/Desktop/", Sample_Name, "/", "Bam", "/", Sample_Name, "_", "Input.bam")
    chipBam <- paste0("~/Desktop/", Sample_Name, "/", "Bam", "/", Sample_Name, "_", "ChIP.bam")
    
    ##Input 
    if(!file.exists(inputBam)){
      align(index = indexed_reference,
            readfile1 = Input_R1,
            readfile2 = Input_R2,
            output_format = "BAM",
            output_file = inputBam,
            color2base = F,
            type = "dna",
            unique = FALSE,
            nBestLocations = 16,
            nthreads = 8,
            sortReadsByCoordinates = TRUE)
    }
    
    ##ChIP 
    if(!file.exists(chipBam)){
      align(index = indexed_reference,
            readfile1 = ChIP_R1,
            readfile2 = ChIP_R2,
            output_format = "BAM",
            output_file = chipBam,
            color2base = F,
            type = "dna",
            unique = FALSE,
            nBestLocations = 16,
            nthreads = 8,
            sortReadsByCoordinates = TRUE)
    }
    
  }
  
  rm(list=ls())
  gc()
  
}

### ----- Sample_Name should be adjusted when running for rDNA and the for whole genome again ---
### ----- we need both the alignments for rDNA an whole genome for downstream analysis ---

#WT
SubRead_Multiple_Aligner_rDNA(  Input_R1 = "/Volumes/rbmData/00_October_2021/1_S1_R1_001.fastq.gz", 
                                Input_R2 = "/Volumes/rbmData/00_October_2021/1_S1_R2_001.fastq.gz", 
                                ChIP_R1 = "/Volumes/rbmData/00_October_2021/2_S2_R1_001.fastq.gz", 
                                ChIP_R2 = "/Volumes/rbmData/00_October_2021/2_S2_R2_001.fastq.gz",
                                Sample_Name = "Smc5-WT-rDNA"    )

#EQ
SubRead_Multiple_Aligner_rDNA(  Input_R1 = "/Volumes/rbmData/00_October_2021/5_S5_R1_001.fastq.gz", 
                                Input_R2 = "/Volumes/rbmData/00_October_2021/5_S5_R2_001.fastq.gz", 
                                ChIP_R1 = "/Volumes/rbmData/00_October_2021/6_S6_R1_001.fastq.gz", 
                                ChIP_R2 = "/Volumes/rbmData/00_October_2021/6_S6_R2_001.fastq.gz",
                                Sample_Name = "Smc5-EQ-rDNA" )

#DA
SubRead_Multiple_Aligner_rDNA(  Input_R1 = "/Volumes/rbmData/00_October_2021/3_S3_R1_001.fastq.gz", 
                                Input_R2 = "/Volumes/rbmData/00_October_2021/3_S3_R2_001.fastq.gz", 
                                ChIP_R1 = "/Volumes/rbmData/00_October_2021/4_S4_R1_001.fastq.gz", 
                                ChIP_R2 = "/Volumes/rbmData/00_October_2021/4_S4_R2_001.fastq.gz",
                                Sample_Name = "Smc5-DA-rDNA"  )

#noTag
SubRead_Multiple_Aligner_rDNA(  Input_R1 = "/Volumes/rbmData/12_Aug_21/Raw/83_1.fastq.gz", 
                                Input_R2 = "/Volumes/rbmData/12_Aug_21/Raw/83_2.fastq.gz", 
                                ChIP_R1 = "/Volumes/rbmData/12_Aug_21/Raw/84_1.fastq.gz", 
                                ChIP_R2 = "/Volumes/rbmData/12_Aug_21/Raw/84_2.fastq.gz",
                                Sample_Name = "Smc5-noTag-rDNA"  )



### ------------ distribute reads mapping to multiple locations with equal propability & calculate ratio ----------
### ----- mask and unmask inside appropriately for rDNA or whole genome ---

EqualDistRatio <- function(binSize = 300, stepSize = 10, Name, AlignmentMode = "multimapping"){
  
  library(Rsubread)
  library(GenomicAlignments)
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  library(ORFik)
  library(rasterpdf)
  library(plotrix)
  
  #function to count coverage
  BinChIPseq = function( reads, bins ){
    mcols(bins)$score = countOverlapsW( bins, reads, weight = "weight" ) 
    return( bins ) 
  }
  
  bamIP <- paste0("/Users/mohammed/Desktop/", Name, "/", "Bam", "/", Name, "_ChIP.bam")
  bamInput <- paste0("/Users/mohammed/Desktop/", Name, "/", "Bam", "/", Name, "_Input.bam")
  
  # binning
  
  #for rDNA reference only
  gen <- Seqinfo(seqnames="Scer_2xrDNA_unit",
                      seqlengths=18274,
                      isCircular=TRUE,
                      genome="Scer_2xrDNA_unit")
  Windows = tileGenome(gen, tilewidth=stepSize)
  Windows = unlist(Windows)
  Windows = suppressWarnings(trim(IRanges::resize(Windows, width=binSize)))
  
  #for whole genome
  # Windows = tileGenome(seqlengths(Scerevisiae)[-17], tilewidth=stepSize)
  # Windows = unlist(Windows)
  # Windows = suppressWarnings(trim(IRanges::resize(Windows, width=binSize)))
  
  # count coverage values
  
  Counts_ChIP <- 0
  Counts_Input <- 0
  
  if(AlignmentMode == "multimapping"){
    
    for(i in 1:16){
      
      paramS <- ScanBamParam(flag = scanBamFlag(isProperPair = TRUE,
                                                isUnmappedQuery = FALSE,
                                                hasUnmappedMate = FALSE), 
                             what = c("qname", "mapq", "isize", "seq"), 
                             tag = c("NH", "NM"), 
                             tagFilter = list(NH = i, NM = c(0:2)))
      
      #input
      Fragments <- suppressWarnings(granges(readGAlignmentPairs(bamInput, param = paramS, strandMode = 1)))
      
      if(length(seqnames(Fragments)) == 0) next
      
      if(length(seqnames(Fragments)) > 0){
        start(Fragments) <- round((start(Fragments) + end(Fragments))/2)
        end(Fragments) <- start(Fragments)
        mcols(Fragments)$weight <- 1/i
      } else {
        next
      }
      
      Counts_Input <- Counts_Input + countOverlapsW(Windows, Fragments, "weight")
      
      #chip
      Fragments <- suppressWarnings(granges(readGAlignmentPairs(bamIP, param = paramS, strandMode = 1)))
      
      if(length(seqnames(Fragments)) == 0) next
      
      if(length(seqnames(Fragments)) > 0){
        start(Fragments) <- round((start(Fragments) + end(Fragments))/2)
        end(Fragments) <- start(Fragments)
        mcols(Fragments)$weight <- 1/i
      } else {
        next
      }
      
      Counts_ChIP <- Counts_ChIP + countOverlapsW(Windows, Fragments, "weight")
      # 
      
    }
    
    ##
    mcols(Windows)$input <- Counts_Input/(binSize/stepSize)
    mcols(Windows)$chip <- Counts_ChIP/(binSize/stepSize)
    
  }
  
  #funtion to calculate ration
  CalculateRatio <- function(Coverages){
    
    Coverages <- as.data.frame(Coverages)
    
    IP_Sum <- sum(as.numeric(Coverages[,"chip"]))
    In_Sum <- sum(as.numeric(Coverages[,"input"]))
    
    corrFactor <- IP_Sum/In_Sum
    
    Ratio <- round(Coverages[,"chip"]/Coverages[,"input"]/corrFactor, 4)
    Ratio[!is.finite(Ratio)] <- 0
    #Ratio <- log2((abs(Ratio))^(sign(Ratio)))
    
    Coverages$corr.input <- round(Coverages[,"input"]*corrFactor)
    Coverages$ratio <- signif(Ratio, 3)
    Coverages$pvalue <- ppois( q = Coverages[,"chip"] - 1, 
                               lambda=Coverages[,"corr.input"], 
                               lower.tail=FALSE, log=FALSE      )
    
    Coverages[,4] <- Name
    
    colnames(Coverages) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'ip.score', 
                             'true.input', 'norm.input', 'ratio', 'pvalue')
    
    return(Coverages)
    
  }
  
  RatioFile <- CalculateRatio(Windows)
  
  #for rDNA 2nts1
  write.table(RatioFile,
              paste0("~/Desktop/", Name, "/", Name, "_Ratio_primary", ".bed"),
              quote=FALSE, row.names=FALSE, sep="\t")
  
  #for whole genome
  # write.table(RatioFile, 
  #             paste0("~/Desktop/", Name, "/", Name, "_Ratio_primary_wG_new", ".bed"),
  #             quote=FALSE, row.names=FALSE, sep="\t")
  
  ###
  
}

EqualDistRatio(Name = "Smc5-WT-rDNA")
EqualDistRatio(Name = "Smc5-noTag-rDNA")
EqualDistRatio(Name = "Smc5-EQ-rDNA")
EqualDistRatio(Name = "Smc5-DA-rDNA")

# run the same function for the whole genome alignments too!


### ------------ calculate averages at genomic features and plotting ----------

# load required packages

packages <- c("shiny", "basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
              "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw", "rtracklayer", "GenomicRanges",
              "BSgenome.Scerevisiae.UCSC.sacCer3", "readxl")
suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))

# parameters 
binSize = 300
stepSize = 10
Colors = c("lightgray", "cornflowerblue", "forestgreen", "brown1", "yellow", "lightgreen")
Fragment <- 300


#rF_2nts1
NumOfSamples <- 4

Sample_1 <- "Smc5-noTag-rDNA"
Sample_2 <- "Smc5-WT-rDNA"
Sample_3 <- "Smc5-DA-rDNA"
Sample_4 <- "Smc5-EQ-rDNA"

#RDN1
rDNA_rF <- data.frame(chrom = "Scer_2xrDNA_unit", chromStart = 1, chromEnd = 18274, name = 'rDNA', mid = 9137)

### Read ratio
# read the resultant ratio files from either within same result folder containing both rDNA and wG (whole-genome) file 
# or from appropriate result folders if they are distinct for rDNA and wG

for(i in 1:NumOfSamples){
  
  assign(paste0("S", i, "_", "Ratio_rF"), read.table(paste0("~/Desktop/", 
                                                            paste0(get(paste0("Sample", "_", i)), "/", get(paste0("Sample", "_", i)), "_", "Ratio_primary.bed")), header = T))
  assign(paste0("S", i, "_", "Ratio_wC"), read.table(paste0("~/Desktop/", 
                                                            paste0(get(paste0("Sample", "_", i)), "/", get(paste0("Sample", "_", i)), "_", "Ratio_primary_wG_new.bed")), header = T))
  
}

### ----- following function calculates per feature average and 
### ----- also simulates the background by sampling across the genome 
### ----- simulated background is then used to normalised the signal against noise

Calc_parallel <- function(Sample){
  
  IP_Combo_rF <- get(paste0(Sample, "_Ratio_rF"))
  IP_Combo_wC <- get(paste0(Sample, "_Ratio_wC"))
  
  Average_Enrichment <- function(
    IP_Combo_1,
    IP_Combo_2,
    PeakList_1,
    Window = 9137,
    NumberOfRepetitions = 100
  ) {
    
    library(matrixStats)
    
    V <- "ratio"
    
    ### ================================
    ### PART 1 — GENOME-WIDE SIMULATION
    ### ================================
    
    RegionSize <- round(2 * Window)                  # ~18274 bp (2x rDNA)
    n_bins <- round((2 * Window) / stepSize)
    
    ChromNames <- seqnames(Scerevisiae)[-17]
    ChromLengths <- seqlengths(Scerevisiae)[-17]
    
    Averages <- matrix(NA, nrow = n_bins, ncol = NumberOfRepetitions)
    
    pb <- txtProgressBar(min = 0, max = NumberOfRepetitions, style = 3)
    
    for (a in seq_len(NumberOfRepetitions)) {
      
      ChunkList <- list()
      
      for (j in seq_along(ChromNames)) {
        
        chr <- ChromNames[j]
        chr_len <- ChromLengths[j]
        
        if (chr_len <= RegionSize) next
        
        start <- sample.int(chr_len - RegionSize + 1, 1)
        end <- start + RegionSize
        
        IP_chr <- IP_Combo_2[IP_Combo_2$chrom == chr, ]
        
        IP_Z <- IP_chr[
          IP_chr$chromStart >= start &
            IP_chr$chromStart <= end,
          V
        ]
        
        length(IP_Z) <- n_bins
        ChunkList[[chr]] <- IP_Z
      }
      
      Mat <- do.call(cbind, ChunkList)
      Averages[, a] <- rowMedians(Mat, na.rm = TRUE)
      
      setTxtProgressBar(pb, a)
    }
    close(pb)
    
    RanSamp <- rowMeans(Averages, na.rm = TRUE)
    RanSamp <- cbind.data.frame(RanSamp)
    rownames(RanSamp) <- paste0("bin", seq_len(nrow(RanSamp)))
    
    ### ==================================
    ### PART 2 — TRUE rDNA ENRICHMENT
    ### ==================================
    
    PeakList_1$AvBstart <- PeakList_1$mid - Window
    PeakList_1$AvBend   <- PeakList_1$mid + Window
    
    # force rDNA peaks onto rDNA pseudo-chromosome
    PeakList_1$chrom <- "Scer_2xrDNA_unit"
    
    chrS <- c(paste0("chr", as.roman(1:16)), "Scer_2xrDNA_unit")
    
    IP_R <- IP_Combo_1
    IP_W <- NULL
    
    for (i in seq_along(chrS)) {
      
      OriList <- PeakList_1[PeakList_1$chrom == chrS[i], ]
      IP_ROs  <- IP_R[IP_R$chrom == chrS[i], ]
      
      if (nrow(OriList) == 0) next
      
      IP_C <- NULL
      
      for (y in seq_len(nrow(OriList))) {
        
        IP_Z <- IP_ROs[
          IP_ROs$chromStart >= OriList$AvBstart[y] &
            IP_ROs$chromStart <= OriList$AvBend[y],
          V
        ]
        
        if (length(IP_Z) >= n_bins) {
          Ratios <- IP_Z[seq_len(n_bins)]
        } else {
          Ratios <- c(IP_Z, rep(0, n_bins - length(IP_Z)))
        }
        
        IP_C <- cbind(IP_C, matrix(Ratios, ncol = 1))
      }
      
      IP_W <- cbind(IP_W, IP_C)
    }
    
    IP_W <- as.data.frame(IP_W)
    
    MedianSignal <- rowMedians(as.matrix(IP_W), na.rm = TRUE)
    Average <- cbind.data.frame(median = MedianSignal)
    rownames(Average) <- paste0("bin", seq_len(nrow(Average)))
    
    ### =========================
    ### PART 3 — NORMALISATION
    ### =========================
    
    Results <- signif(Average / (RanSamp * 200), 5)
    Results$median[!is.finite(Results$median)] <- 0
    
    return(Results)
  }
  
  IP_AvE_rDNA <- Average_Enrichment(IP_Combo_1 = IP_Combo_rF, IP_Combo_2 = IP_Combo_wC, PeakList_1 = rDNA_rF)
  
  Res <- list(IP_AvE_rDNA)
  
  names(Res) <- c('IP_AvE_rDNA')
  
  return(Res)
  
}

### as rDNA is only one genomic feature - we reiterate the above function 10, 100, or 1000 times 

for(k in 1:100){
  print(paste0('k = ', k))
  for(i in 1:NumOfSamples){
    assign(paste0("S", k, i), Calc_parallel(Sample = paste0("S", i)))
  }
}

for(i in 1:NumOfSamples){
  ETTT <- NULL
  for(k in 1:100){
    median <- get('IP_AvE_rDNA', eval(as.symbol(paste0("S", k, i))))[,1]
    ETTT <- cbind(ETTT, median)
  }
  assign(paste0("S", i), list(IP_AvE_rDNA = cbind.data.frame(median = rowMeans(ETTT))))
}


#plotting

Name <- paste0(unlist(strsplit(Sample_1, "-"))[1])

par(oma=c(0,0,0,0))

PlotMat <- {matrix(c(
  0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
  0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
  0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
  0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,
  0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,
  0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,
  0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,
  0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,
  0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,
  0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
  3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), 
  
  19,20,byrow=TRUE)}

layout(PlotMat, c(1,1), c(1,1), TRUE)

PlotEnrichments <- function(DataType, Window = 9137){
  
  SampName <- c("noTag", "WT", "DA", "EQ")
  
  DataRange <- c()
  
  for(i in 1:NumOfSamples){
    
    assign(paste0("AvE", "_", i), get(DataType, eval(as.symbol(paste0("S", i)))))
    
    Med <- get(paste0("AvE", "_", i))$median
    Med[!is.finite(Med)] <- 0
    
    DataRange <- c(DataRange, Med)
    
    Y <- max(round(abs(range(DataRange))+0.5))
    
    Y <- 0.1
    
  }
  
  S <- "rDNA 2xNTS1 18.2 kb | ChIP/Input over background"
  
  X <- (Window/stepSize)*2
  
  plot(NULL,
       ylim = c(0, Y),
       main = S,
       ylab = "Normalised density", cex.main=1.25, cex.lab=1.25, cex.axis=1.25, xlab = "Distance from center (Kbp)", 
       xaxt = "n", bty = 'n', xaxs = 'i', yaxs = 'i', xlim=c(0, X), las=1)
  legend("topleft", legend = SampName, lwd = 2, bty = "n", 
         col = Colors[1:NumOfSamples], horiz = T, cex = 1.25, pt.cex = 1.25)
  
  
  # rect((10000-8755)/10, -0.25, (10000+8755)/10, 0.25, col = adjustcolor('red', 0.3), border = "transparent", xpd = T)
  
  for(i in 1:NumOfSamples){
    
    Counts <- get(paste0("AvE", "_", i))$median
    Counts[!is.finite(Counts)] <- 0
    
    Counts <- as.numeric(smooth.spline(1:length(Counts), Counts)$y)
    # smooth.spline(1:length(Counts), Counts, spar = 0.5)$y
    points(Counts, col = Colors[i], type = 'l', lwd = 2)
    
  }
  
  axisLabels <- seq(0,
                    18274,
                    length.out = 9)
  axisLabels[c(2,4,6,8)] <- NA
  At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
  axis(1, at=At, labels = signif(axisLabels/1000, 2))
  
  # segments((10000-8755)/10, -0.5, (10000-8755)/10, 0.5, lwd = 2, xpd = T)
  # segments((10000+8755)/10, -0.5, (10000+8755)/10, 0.5, lwd = 2, xpd = T)
  # segments((10000+380)/10, -0.5, (10000+380)/10, 0.5, lwd = 2, xpd = T)
  # 
  # segments((10000-498)/10, -0.5, (10000-498)/10, 0.5, lwd = 2, xpd = T)
  # #segments((10000-498)/10, 0.5, (10000-1300)/10, 1, lwd = 1, lty = 2, xpd = T)
  # 
  # segments((10000-1190)/10, -0.5, (10000-1190)/10, 0.5, lwd = 2, xpd = T)
  # segments((10000-2442)/10, -0.5, (10000-2442)/10, 0.5, lwd = 2, xpd = T)
  # #segments((10000-2442)/10, 0.5, (10000-3000)/10, 1, lwd = 1, lty = 2, xpd = T)
  # 
  # draw.circle((10000-1190)/10, 0, radius = 20, col = 'blue')
  # mtext('rARS', side = 1, line = 1, at = (10000-1190)/10, cex = 0.75)
  # 
  #segments((10000+380)/10, 0.5, (10000+1200)/10, 1, lwd = 1, lty = 2, xpd = T)
  
  #text((10000-8700)/10, 4, "left rDNA\njunction", srt = 45)
  #text((10000+8800)/10, 4, "right rDNA\njunction", srt = 45)
  #text((10000-1200)/10, 4, "5S rRNA", srt = 45)
  #text((10000-3000)/10, 4, "35S rRNA", srt = 45)
  #text((10000+1900)/10, 4, "start of last\nrDNA repeat", srt = 45)
  
}

pdf(paste0("~/Desktop/", Sample_1, "/", Name, "_rDNA_2nts1_chunk_ratio.pdf"), width = 11, height = 10)

PlotEnrichments(DataType = "IP_AvE_rDNA")

dev.off()

