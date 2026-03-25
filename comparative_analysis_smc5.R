# Comparative analysis

NumOfSamples <- 4

Sample_1 <- "Smc5-noTag-ChIP"
Sample_2 <- "Smc5-WT-ChIP"
Sample_3 <- "Smc5-DA-ChIP"
Sample_4 <- "Smc5-EQ-ChIP"


# load required packages

packages <- c("shiny", "basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
              "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw", "rtracklayer", "GenomicRanges",
              "BSgenome.Scerevisiae.UCSC.sacCer3", "readxl")
suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))

# read genomic coords for specific elements of interest

#ARS
E_Ori <- read.table("E_Rep.bed", header = F, quote = "/t")[ ,1:4]
L_Ori <- read.table("L_Rep.bed", header = F, quote = "/t")[ ,1:4]
All_Ori <- read.table("OriginList_Full.bed", header = TRUE, quote = "/t")[ ,1:4]

colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); E_Ori$mid <- round((E_Ori$chromStart + E_Ori$chromEnd)/2)
colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); L_Ori$mid <- round((L_Ori$chromStart + L_Ori$chromEnd)/2)
colnames(All_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); All_Ori$mid <- round((All_Ori$chromStart + All_Ori$chromEnd)/2)

All_Ori$name <- "ARS"; L_Ori$name <- "ARS"; E_Ori$name <- "ARS"

#tRNA
tRNA_list <- as.data.frame(read_excel("tRNA_list.xlsx", col_names=F))[,c(9, 10, 11, 12)]
colnames(tRNA_list) <- c("chrom", "chromStart", "chromEnd", "strand"); tRNA_list$mid <- round((tRNA_list$chromStart + tRNA_list$chromEnd)/2)
Chromosome <- paste0('chr', as.roman(1:16))
for(i in 1:length(Chromosome)){
  tRNA_list$chrom[tRNA_list$chrom==i] <- Chromosome[i]
}
tRNA_list <- tRNA_list[tRNA_list$chrom != 17, ]
colnames(tRNA_list)[4] <- "name"
tRNA_list$name <- "tRNA"

#TyElement
TyElement_list <- as.data.frame(read_excel("TyElement_list.xlsx", col_names=F))[,c(9, 10, 11, 12)]
colnames(TyElement_list) <- c("chrom", "chromStart", "chromEnd", "strand"); TyElement_list$mid <- round((TyElement_list$chromStart + TyElement_list$chromEnd)/2)
Chromosome <- paste0('chr', as.roman(1:16))
for(i in 1:length(Chromosome)){
  TyElement_list$chrom[TyElement_list$chrom==i] <- Chromosome[i]
}
TyElement_list <- TyElement_list[TyElement_list$chrom != 17, ]
colnames(TyElement_list)[4] <- "name"
TyElement_list$name <- "TyElement"

#Covergent
CONVERGENT_list <- as.data.frame(read_excel("CONVERGENT_list.xlsx", col_names=T))
colnames(CONVERGENT_list) <- c("chrom", "chromStart", "chromEnd")
Chromosome <- paste0('chr', as.roman(1:16))
for(i in 1:length(Chromosome)){
  CONVERGENT_list$chrom[CONVERGENT_list$chrom==i] <- Chromosome[i]
}
CONVERGENT_list <- CONVERGENT_list[CONVERGENT_list$chrom != 17, ]
CONVERGENT_list$name <- "CONVERGENT"
CONVERGENT_list$mid <- round((CONVERGENT_list$chromStart + CONVERGENT_list$chromEnd)/2)

#Divergent
DIVERGENT_list <- as.data.frame(read_excel("DIVERGENT_list.xlsx", col_names=T))
colnames(DIVERGENT_list) <- c("chrom", "chromStart", "chromEnd")
Chromosome <- paste0('chr', as.roman(1:16))
for(i in 1:length(Chromosome)){
  DIVERGENT_list$chrom[DIVERGENT_list$chrom==i] <- Chromosome[i]
}
DIVERGENT_list <- DIVERGENT_list[DIVERGENT_list$chrom != 17, ]
DIVERGENT_list$name <- "DIVERGENT"
DIVERGENT_list$mid <- round((DIVERGENT_list$chromStart + DIVERGENT_list$chromEnd)/2)

#Ctrans
CTrans_list <- as.data.frame(read_excel("CTrans_list.xlsx", col_names=T))
colnames(CTrans_list) <- c("chrom", "chromStart", "chromEnd")
Chromosome <- paste0('chr', as.roman(1:16))
for(i in 1:length(Chromosome)){
  CTrans_list$chrom[CTrans_list$chrom==i] <- Chromosome[i]
}
CTrans_list <- CTrans_list[CTrans_list$chrom != 17, ]
CTrans_list$name <- "CTrans"
CTrans_list$mid <- round((CTrans_list$chromStart + CTrans_list$chromEnd)/2)

#WTrans
WTrans_list <- as.data.frame(read_excel("WTrans_list.xlsx", col_names=T))
colnames(WTrans_list) <- c("chrom", "chromStart", "chromEnd")
Chromosome <- paste0('chr', as.roman(1:16))
for(i in 1:length(Chromosome)){
  WTrans_list$chrom[WTrans_list$chrom==i] <- Chromosome[i]
}
WTrans_list <- WTrans_list[WTrans_list$chrom != 17, ]
WTrans_list$name <- "WTrans"
WTrans_list$mid <- round((WTrans_list$chromStart + WTrans_list$chromEnd)/2)


#CENTROMERES

library(Repliscope)
data("sacCer3")
Cen <- sacCer3$cen
Cen$mid <- Cen$chromEnd

# parameters 

Normalise2Input = "NO"
Window = 2000
stepSize = 10
Colors = c("lightgray", "cornflowerblue", "forestgreen", "brown1", "yellow", "lightgreen")



# Read strand-wise ratio files from the previously generated data

for(i in 1:NumOfSamples){
  
  assign(paste0("S", i, "_", "IP_watson"), read.table(paste0(get(paste0("Sample", "_", i)), "/", "Ratios", "/", get(paste0("Sample", "_", i)), "_", "ChIP_watson.bed"), header = T))
  assign(paste0("S", i, "_", "IP_crick"), read.table(paste0(get(paste0("Sample", "_", i)), "/", "Ratios", "/", get(paste0("Sample", "_", i)), "_", "ChIP_crick.bed"), header = T))
  
}

# read peak files

for(i in 1:NumOfSamples){
  assign(paste0("S", i, "_", "Peaks"), read.table(paste0(get(paste0("Sample", "_", i)), "/", "Peaks", "/", get(paste0("Sample", "_", i)), "_", "Primary_ALL_Peaks.bed"), header = T))
}


# consensus peaks
peak_files <- c()
for(i in 1:NumOfSamples){
  peak_files <- c(peak_files, paste0(get(paste0("Sample", "_", i)), "/", "Peaks", "/", get(paste0("Sample", "_", i)), "_", "Primary_ALL_Peaks.bed"))
}

peak_granges <- lapply(peak_files[-1], import)

peak_grangeslist <- GRangesList(peak_granges)

peak_coverage <- coverage(peak_grangeslist)

covered_ranges <- slice(peak_coverage, lower=2, rangesOnly=T)

covered_granges <- GRanges(covered_ranges)

covered_granges <- reduce(covered_granges, min.gapwidth=51)

# export(covered_granges, "consensus.bed")

MetaPeaks <- as.data.frame(covered_granges)

chrOrder <- paste0('chr', as.roman(1:16))
MetaPeaks[,1] <- factor(MetaPeaks[,1], levels=chrOrder)
MetaPeaks <- MetaPeaks[order(MetaPeaks[,1]),]
MetaPeaks <- MetaPeaks[,c(1,2,3,5)]


colnames(MetaPeaks) <- c("chrom", "chromStart", "chromEnd", "name"); MetaPeaks$mid <- round((MetaPeaks$chromStart + MetaPeaks$chromEnd)/2)
MetaPeaks$name <- "MetaPeaks"

# check fragment size

bamFragments <- c()
for(i in 1:NumOfSamples){
  bamFragments <- c(bamFragments,
                    mean(getPESizes(paste0(get(paste0("Sample", "_", i)), "/", "Bam", "/", get(paste0("Sample", "_", i)), "_", "ChIP.bam"),
                                    param=readParam(pe="both"))$sizes))
}
Fragment <- round_any(mean(bamFragments), 100)

# Calculations

Calc <- function(Sample){
  
  IP_watson <- get(paste0(Sample, "_IP_watson"))
  IP_crick <- get(paste0(Sample, "_IP_crick"))
  
  IP_Combo <- data.frame(chrom = IP_watson$chrom,
                         chromStart = IP_watson$chromStart,
                         chromEnd = IP_watson$chromEnd,
                         name = unlist(strsplit(IP_watson$name[1], "_"))[1],
                         ip.score = IP_watson$ip.score + IP_crick$ip.score)
  
  
  Average_Enrichment <- function(IP_Combo, PeakList){
    
    Window = Window
    V <- 5
    
    #Random Sampling
    Total_NumberOfelements <- length(PeakList$chrom)
    
    GetDistances <- function(Positions){
      return(Positions[2: length(Positions)] - Positions[1:(length(Positions)-1)])
    }
    
    rowStat <- function(DF){
      
      Quantiles <- apply(DF, 1, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
      Ses <- apply(DF, 1, std.error)
      Means <- apply(DF, 1, mean)
      Sds <- apply(DF, 1, sd)
      
      return(list(q.25 = Quantiles["25%",],
                  Median = Quantiles["50%",],
                  q.75 = Quantiles["75%",],
                  Std.err = Ses,
                  Mean = Means,
                  Sd = Sds))
    }
    
    sampleByMinDistance <- function(Vect, Size, MinDistance){
      
      Reps <- 0
      
      Samp <- TRUE
      while(Samp){
        
        Reps <- Reps + 1
        # if(Reps %% 1000 == 1)
        #   print(paste("Working on", Reps))
        
        x1 <- sort(round(runif(Size, min=0+MinDistance, max=Vect-MinDistance)))
        xdist <- GetDistances(x1)
        
        while(min(xdist)<MinDistance){
          # print('working' )
          x2 <- which.min(xdist) + sample(x = c(0,1), size = 1)
          x1 <- x1[-x2]
          x1 <- sort(c(x1, round(runif(1, min=0+MinDistance, max=Vect-MinDistance))))
          xdist <- GetDistances(x1)
        }
        
        #x1 <- c(0, x1, Vect)
        
        return(x1)
      }
    }
    
    ChromLengths <- seqlengths(Scerevisiae)[-17]
    
    Chrom_NumberOfelements <- c()
    for (i in 1:length(unique(PeakList$chrom))){
      Chrom_NumberOfelements <- c(Chrom_NumberOfelements, length(PeakList[PeakList$chrom==unique(PeakList$chrom)[i], ]$chrom))
    }
    
    ChrEleDF <- data.frame(chrom = unique(PeakList$chrom), numElements = Chrom_NumberOfelements)
    
    MinDists <- c()
    for (i in 1:length(unique(PeakList$chrom))){
      if(Chrom_NumberOfelements[i]==0) next 
      MinDists <- c(MinDists, min(GetDistances(PeakList[PeakList$chrom==unique(PeakList$chrom)[i], ]$mid)))
    }
    
    
    NumberOfRepetitions <- 2
    
    Averages <- as.data.frame(matrix(0, ncol = NumberOfRepetitions, nrow = (Window*2)/stepSize))
    
    for(a in 1: NumberOfRepetitions){
      
      
      Sampled.Elements <- list()
      
      for(j in 1:length(ChromLengths)){
        
        print(paste("SAMPLE", "-", unique(IP_Combo$name), ",", "Element", "-", unique(PeakList$name), ",", "Running simulation", a, "for chromosome", j))
        
        Rep <- TRUE
        while(Rep){
          
          SampElements <- sampleByMinDistance(ChromLengths[j], Chrom_NumberOfelements[j], MinDists[j])
          UniVec <- length(unique(SampElements))
          
          if(length(SampElements)==UniVec){
            Rep <- FALSE
          } else {
            print("Repeating the sampling")
          }
        }
        Sampled.Elements[[j]] <- SampElements
      }
      names(Sampled.Elements) <- seqnames(Scerevisiae)[-17]
      
      ElementsList <- NULL
      for(k in 1:length(names(Sampled.Elements))){
        List <- data.frame(chrom=names(Sampled.Elements)[k], elements = Sampled.Elements[[k]])
        ElementsList <- rbind.data.frame(ElementsList, List)
      }
      
      ElementsList$AvBstart <- ElementsList$elements - Window
      ElementsList$AvBend <- ElementsList$elements + Window
      
      
      IP_R <- IP_Combo
      IP_W <- NULL
      for(i in 1:length(names(Sampled.Elements))){
        
        i <- i
        
        Elements <- ElementsList[ElementsList$chrom == names(Sampled.Elements)[i], ]
        IP_ROs <- IP_R[IP_R$chrom == names(Sampled.Elements)[i], ]
        
        if(length(Elements$chrom)==0) next
        
        IP_C <- NULL
        for(y in 1:length(Elements$chrom)){
          
          y <- y
          
          IP_Z <- IP_ROs[IP_ROs$chromStart>=Elements$AvBstart[y]-150 & IP_ROs$chromStart<=Elements$AvBend[y]-150, ]
          
          if(length(IP_Z[,V])==round(2*Window/stepSize)){
            Ratios <- IP_Z[,V]
          }
          
          if(length(IP_Z[,V]) < round(2*Window/stepSize)){
            if(length(1:(length(IP_Z[,V])/2)) < round(2*Window/stepSize)/2){
              Ratios <- c(rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)), IP_Z[,V] )
            }
            if(length((length(IP_Z[,V])/2+1):length(IP_Z[,V])) < round(2*Window/stepSize)/2){
              Ratios <- c(IP_Z[,V], rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)) )
            }
          }
          
          if(length(IP_Z[,V]) > round(2*Window/stepSize)){
            Ratios <- IP_Z[,V][-c((round(2*Window/stepSize)+1):length(IP_Z[,V]))]
          }
          
          Rat <- matrix(Ratios, ncol = 1)
          
          IP_C <- cbind(IP_C, Rat)
        }
        IP_W <- cbind(IP_W, IP_C)
      }
      IP_W <- as.data.frame(IP_W)
      colnames(IP_W) <- c(1:length(ElementsList$chrom))
      
      Median <- cbind.data.frame(median = rowStat(IP_W)$Median)
      
      Averages <- cbind.data.frame(Averages, Median)
      
      Averages <- Averages[,(dim(Averages)[2]-NumberOfRepetitions+1):dim(Averages)[2]]
      
    }
    
    RanSamp <- rowMeans(Averages)
    RanSamp <- cbind.data.frame(RanSamp)
    rownames(RanSamp) <- paste0("bin", 1:length(RanSamp[,1]))
    
    ###
    
    PeakList$AvBstart <- PeakList$mid - Window
    PeakList$AvBend <- PeakList$mid + Window
    
    chrS <- paste0("chr", as.roman(1:16))
    
    IP_R <- IP_Combo
    IP_W <- NULL
    for(i in 1:length(chrS)){
      
      i <- i
      
      OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
      IP_ROs <- IP_R[IP_R$chrom == chrS[i], ]
      
      if(length(OriList_ROs$chrom)==0) next 
      
      IP_C <- NULL
      for(y in 1:length(OriList_ROs$chrom)){
        
        y <- y
        
        IP_Z <- IP_ROs[IP_ROs$chromStart>=OriList_ROs$AvBstart[y]-150 & IP_ROs$chromStart<=OriList_ROs$AvBend[y]-150, ]
        
        if(length(IP_Z[,V])==round(2*Window/stepSize)){
          Ratios <- IP_Z[,V]
        } 
        
        if(length(IP_Z[,V]) < round(2*Window/stepSize)){
          if(length(1:(length(IP_Z[,V])/2)) < round(2*Window/stepSize)/2){
            Ratios <- c(rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)), IP_Z[,V] )
          }
          if(length((length(IP_Z[,V])/2+1):length(IP_Z[,V])) < round(2*Window/stepSize)/2){
            Ratios <- c(IP_Z[,V], rep(0, (round(2*Window/stepSize))-length(IP_Z$chromStart)) )
          }
        } 
        
        if(length(IP_Z[,V]) > round(2*Window/stepSize)){
          Ratios <- IP_Z[,V][-c((round(2*Window/stepSize)+1):length(IP_Z[,V]))]
        }
        
        Rat <- matrix(Ratios, ncol = 1)
        
        IP_C <- cbind(IP_C, Rat)
      }
      IP_W <- cbind(IP_W, IP_C)
    }
    IP_W <- as.data.frame(IP_W)
    colnames(IP_W) <- c(1:length(PeakList$chrom))
    
    rowStat <- function(DF){
      
      Quantiles <- apply(DF, 1, quantile, probs = c(0.25, 0.50, 0.75), na.rm = T)
      Ses <- apply(DF, 1, std.error)
      Means <- apply(DF, 1, mean)
      Sds <- apply(DF, 1, sd)
      
      return(list(q.25 = Quantiles["25%",],
                  Median = Quantiles["50%",],
                  q.75 = Quantiles["75%",],
                  Std.err = Ses,
                  Mean = Means,
                  Sd = Sds))
    }
    
    Average <- cbind.data.frame(median = rowStat(IP_W)$Median)
    
    rownames(Average) <- paste0("bin", 1:length(Average[,1]))
    
    ###
    
    Results <- signif(Average/RanSamp, 3)
    
    return(Results)
  }
  
  IP_AvE_Peaks <- Average_Enrichment(IP_Combo=IP_Combo, PeakList = MetaPeaks)
  IP_AvE_ARS <- Average_Enrichment(IP_Combo=IP_Combo, PeakList = All_Ori)
  IP_AvE_tRNA <- Average_Enrichment(IP_Combo=IP_Combo, PeakList = tRNA_list)
  IP_AvE_TyElement <- Average_Enrichment(IP_Combo=IP_Combo, PeakList = TyElement_list)
  
  IP_AvE_Convergent <- Average_Enrichment(IP_Combo=IP_Combo, PeakList = CONVERGENT_list)
  IP_AvE_Divergent <- Average_Enrichment(IP_Combo=IP_Combo, PeakList = DIVERGENT_list)
  IP_AvE_Ctrans <- Average_Enrichment(IP_Combo=IP_Combo, PeakList = CTrans_list)
  IP_AvE_Wtrans <- Average_Enrichment(IP_Combo=IP_Combo, PeakList = WTrans_list)
  
  Res <- list(IP_AvE_Peaks, IP_AvE_ARS, IP_AvE_tRNA, IP_AvE_TyElement,
              IP_AvE_Convergent, IP_AvE_Divergent, IP_AvE_Ctrans, IP_AvE_Wtrans)
  
  names(Res) <- c('IP_AvE_Peaks', 'IP_AvE_ARS', 'IP_AvE_tRNA', 'IP_AvE_TyElement', 
                  "IP_AvE_Convergent", "IP_AvE_Divergent", "IP_AvE_Ctrans", "IP_AvE_Wtrans")
  
  return(Res)
  
}

for(i in 1:NumOfSamples){
  assign(paste0("S", i), Calc(Sample = paste0("S", i)))
}


# Plotting

PlotEnrichments <- function(DataType){
  
  DataRange <- c()
  
  for(i in 1:NumOfSamples){
    
    assign(paste0("AvE", "_", i), get(DataType, eval(as.symbol(paste0("S", i)))))
    
    DataRange <- c(DataRange, get(paste0("AvE", "_", i))$median)
    
    Y <- max(round(abs(range(DataRange))+0.5))
    
  }
  
  S <- unlist(strsplit(DataType, "_"))[3]
  
  X <- (Window/stepSize)*2
  
  plot(NULL,
       ylim = c(0, Y),
       main = paste0(S, "-", "Enrichments"),
       ylab = "Normalised density", cex.main=1.25, cex.lab=1.25, cex.axis=1.25, xlab = "Distance from center (Kbp)", 
       xaxt = "n", bty = 'n', xaxs = 'i', yaxs = 'i', xlim=c(0, X), las=1)
  
 
  for(i in 1:NumOfSamples){
    
    Counts <- round(as.numeric(smooth.spline(1:length(get(paste0("AvE", "_", i))$median), get(paste0("AvE", "_", i))$median)$y), 2)
    
    points(Counts, col = Colors[i], type = 'l', lwd = 2)
    
  }
  
  abline(h=0,lwd=0.4); abline(v=(length(Counts))/2,lwd=0.4)
  axisLabels <- seq(-Window,
                    +Window,
                    length.out = 9)
  axisLabels[c(2,4,6,8)] <- NA
  At <- (Window/stepSize)*seq(0,2,0.25); At[1] <- 1
  axis(1, at=At, labels = signif(axisLabels/1000, 2))
  
}

Name <- paste0(unlist(strsplit(Sample_1, "-"))[1])

pdf(paste0(Name, "_shifted.pdf"), width = 11, height = 10)

par(oma=c(0,0,0,0))

PlotMat <- {matrix(c(
                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                          0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
                          0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
                          0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,
                          0,2,2,2,2,2,2,2,2,0,0,3,3,3,3,3,3,3,3,0,
                          0,2,2,2,2,2,2,2,2,0,0,3,3,3,3,3,3,3,3,0,
                          0,2,2,2,2,2,2,2,2,0,0,3,3,3,3,3,3,3,3,0,
                          0,2,2,2,2,2,2,2,2,0,0,3,3,3,3,3,3,3,3,0,
                          0,2,2,2,2,2,2,2,2,0,0,3,3,3,3,3,3,3,3,0,
                          0,2,2,2,2,2,2,2,2,0,0,3,3,3,3,3,3,3,3,0,
                          0,2,2,2,2,2,2,2,2,0,0,3,3,3,3,3,3,3,3,0,
                          0,4,4,4,4,4,4,4,4,0,0,5,5,5,5,5,5,5,5,0,
                          0,4,4,4,4,4,4,4,4,0,0,5,5,5,5,5,5,5,5,0,
                          0,4,4,4,4,4,4,4,4,0,0,5,5,5,5,5,5,5,5,0,
                          0,4,4,4,4,4,4,4,4,0,0,5,5,5,5,5,5,5,5,0,
                          0,4,4,4,4,4,4,4,4,0,0,5,5,5,5,5,5,5,5,0,
                          0,4,4,4,4,4,4,4,4,0,0,5,5,5,5,5,5,5,5,0,
                          0,4,4,4,4,4,4,4,4,0,0,5,5,5,5,5,5,5,5,0,
                          
                          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                   
                   19,20,byrow=TRUE)}

layout(PlotMat, c(1,1), c(1,1), TRUE)

SampName <- c()
for(i in 1:NumOfSamples){
  SampName <- c(SampName, (get(paste0("Sample_", i))))
}

plot(NULL, xaxt = "n", yaxt = "n", xlab = " ", ylab = " ", bty = "n", xlim=c(0,2), ylim=c(0,2))
     legend("center", legend = SampName, lwd = 2, bty = "n", 
            col = Colors[1:NumOfSamples], horiz = T, cex = 1.25, pt.cex = 1.25)

PlotEnrichments(DataType = "IP_AvE_Peaks")
PlotEnrichments(DataType = "IP_AvE_ARS")
PlotEnrichments(DataType = "IP_AvE_tRNA")
PlotEnrichments(DataType = "IP_AvE_TyElement")

plot(NULL, xaxt = "n", yaxt = "n", xlab = " ", ylab = " ", bty = "n", xlim=c(0,2), ylim=c(0,2))
legend("center", legend = SampName, lwd = 2, bty = "n", 
       col = Colors[1:NumOfSamples], horiz = T, cex = 1.25, pt.cex = 1.25)

PlotEnrichments(DataType = "IP_AvE_Convergent")
PlotEnrichments(DataType = "IP_AvE_Divergent")
PlotEnrichments(DataType = "IP_AvE_Ctrans")
PlotEnrichments(DataType = "IP_AvE_Wtrans")


dev.off()


###############################
###############################

Sample_1 <- "Smc5-noTag-ChIP"
Sample_2 <- "Smc5-WT-ChIP"
Sample_3 <- "Smc5-DA-ChIP"
Sample_4 <- "Smc5-EQ-ChIP"


Y_Lim <- 10
stepSize <- 10

All_Ori <- read.table("OriginList_Full.bed", header = TRUE, quote = "/t")

ChIP_Plot_function <- function(CoverageFile, peakFile, DataType, SamplType){
  
  CoverageFile = CoverageFile
  peakFile = peakFile
  
  scoreVals <- c(CoverageFile$wat.score, CoverageFile$cri.score)
  
  # Q <- quantile(scoreVals, probs=c(.01, .99), na.rm = FALSE)
  # I <- IQR(scoreVals)
  # up  <-  Q[2]+1.5*I # Upper Range
  # low <- Q[1]-1.5*I # Lower Range
  # Scores <- scoreVals[which(scoreVals < up & scoreVals > low)]
  # thd <- round(mean(Scores)+(12*sd(Scores)))
  
  thd <- Y_Lim
  
  Ylim_sc <- c(-thd, thd)
  
  yAxis_reads <- c(-0.3*thd, -0.6*thd, -0.9*thd, 0*thd, 0.3*thd, 0.6*thd, 0.9*thd)
  
  Coverage_chr <- CoverageFile[CoverageFile$chrom == seqnames(Scerevisiae)[k], ]
  All_Ori_chr <- All_Ori[All_Ori$chrom == seqnames(Scerevisiae)[k], ]
  Peaks2Plot_chr <- peakFile[peakFile$chrom == seqnames(Scerevisiae)[k], ]
  
  Coverage <-  Coverage_chr[Coverage_chr$chromStart>=S & Coverage_chr$chromStart<=E, ]
  Ori_chr <- All_Ori_chr[All_Ori_chr$chromStart>=S & All_Ori_chr$chromStart<=E, ]
  PeakReg <- Peaks2Plot_chr[Peaks2Plot_chr$peakStart>=S & Peaks2Plot_chr$peakStart<=E, ]
  
  CovWat <- Coverage$wat.score
  CovCri <- Coverage$cri.score
  
  ###
  steps <- 10
  
  if(stepSize == steps){
    CovWat <- CovWat
    CovCri <- CovCri
  }
  if(stepSize > steps){
    s <- stepSize/steps
    CovWat <- as.vector(sapply(CovWat, function (x) rep(x,s)))
    CovCri <- as.vector(sapply(CovCri, function (x) rep(x,s)))
  }
  if(stepSize < steps){
    s <- round(steps/stepSize)
    CovWat <- round(as.vector(tapply(CovWat, gl(length(CovWat)/s, s), mean)))
    CovCri <- round(as.vector(tapply(CovCri, gl(length(CovCri)/s, s), mean)))
  }
  
  ###
  
  WCrat <- log2(CovWat/CovCri)
  WCrat[!is.finite(WCrat)] <- 0
  
  Cov.wat <- round(as.numeric(smooth.spline(1:length(CovWat), CovWat)$y), 3)
  Cov.cri <- round(as.numeric(smooth.spline(1:length(CovCri), CovCri)$y), 3)
  
  WC.rat <- round(as.numeric(smooth.spline(1:length(WCrat), WCrat)$y), 3)
  
  if(length(Cov.wat) < (Length_per_Row/steps - 1)){
    Cov.wat <- c(Cov.wat, rep(NA, (Length_per_Row/steps - 1)-length(Cov.wat)) )
  }
  
  if(length(Cov.cri) < (Length_per_Row/steps - 1)){
    Cov.cri <- c(Cov.cri, rep(NA, (Length_per_Row/steps - 1)-length(Cov.cri)) )
  }
  
  if(length(WC.rat) < (Length_per_Row/steps - 1)){
    WC.rat <- c(WC.rat, rep(NA, (Length_per_Row/steps - 1)-length(WC.rat)) )
  }
  
  #plot
  par(mar = c(0,0,0,0))
  suppressWarnings(
    plot(Cov.wat, type='h', ylim=Ylim_sc, col =  rgb(100,0,0,alpha=180, maxColorValue=255), ylab=' ',
         xlab=' ', xaxt='n', yaxt='n', lwd=0.07, bty = 'n',  cex.lab=1, las = 2, xaxs='i')
  )
  lines(Cov.cri*(-1), type='h', col =  rgb(0,100,0,alpha=180, maxColorValue=255), lwd=0.07)
  
  segments(x0 = c(Coverage$chromStart - S)/steps, x1 = c(Coverage$chromEnd - S)/steps, y0 = 0, y1= 0, lwd = 0.5)
  
  #plot left y axis
  axis(side = 2, at = yAxis_reads, labels = round(yAxis_reads), line = 0, tick = TRUE, lwd.ticks = 1.5, las = 2, cex.axis = 0.8)
  abline(h=yAxis_reads, lwd=0.05, col =  rgb(112,128,144,alpha=225, maxColorValue=255))
  
  # #draw the peaks
  # if(DataType=="strandedRatio"){
  #   if(length(PeakReg$peakStart) > 0){
  #     for(i in 1:length(PeakReg$peakStart)){
  #       segments(x0 = c(PeakReg$peakStart[i] - S)/steps, x1 = c(PeakReg$peakEnd[i] - S)/steps,
  #                y0 = par('usr')[4]-(thd*0.1), y1= par('usr')[4]-(thd*0.1), lwd = 4, col = "red", xpd = TRUE)
  #       draw.circle(x = (PeakReg$peakSummit[i] - S)/steps, y = par('usr')[4]-(thd*0.1), radius = 50, border = "yellow", lwd=2, col="blue")
  #     }
  # 
  #     # if((PeakReg$peakSummit[length(PeakReg$peakStart)] - S)/steps<9000){
  #     #   mtext(side=3, line=-0.70, at=9900, adj=1, cex=0.5, "macs2-peaks", col = 'red')
  #     # }
  #   }
  # }
  
  #draw Replication origins
  
  if(length(Ori_chr$chromStart) > 0){
    for(i in 1:length(Ori_chr$chromStart)){
      draw.circle(x = ((Ori_chr$chromStart[i]+Ori_chr$chromEnd[i])/2 - S)/steps, y = 0, radius = 60, border = "purple", lwd=2, col="yellow")
    }
  }
  
  #put chromosome name
  
  if(SamplType==Sample_1 & DataType=="strandedRatio"){
    title(main = paste("Chromosome", gsub("[[:punct:]]*chr[[:punct:]]*", "", seqnames(Scerevisiae)[k])), col="gray", adj = 0, cex.main=1.5, line = 0, outer = TRUE)
  }
  
  #put the datatype
  if(DataType=="strandedRatio"){
    mtext(side=3, line=-1.40, at=100, adj=0, cex=0.85, SamplType)
  }
  # if(DataType=="strandedRatio"){
  #   if(length(PeakReg$peakStart) == 0){
  #     mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
  #   }
  #   if(length(PeakReg$peakStart) > 0){
  #     if((PeakReg$peakSummit[1] - S)/steps<1500){
  #       mtext(side=1, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
  #     } else {
  #       mtext(side=3, line=-1.25, at=100, adj=0, cex=0.85, "normalised-Ratio")
  #     }
  #   }
  # }
  
  #put origin names
  
  S1 <- Starts[seq(1, length(Starts), 1)]
  #S2 <- Starts[seq(2, length(Starts), 2)]
  #S3 <- Starts[seq(3, length(Starts), 3)]
  
  if(length(Ori_chr$chromStart) > 0){
    
    if(DataType=="strandedRatio" & SamplType==Sample_1){
      
      Draw_name <- function(Ss, x){
        
        colors_ars <- rep("gray1", length(Ori_chr$name))
        colors_ars[which(Ori_chr$stat=='early')] <- "red"
        colors_ars[which(Ori_chr$stat=='late')] <- "blue"
        
        Ori_ticks <- round((c(Ori_chr$chromStart+Ori_chr$chromEnd)/2 - S)/steps)
        
        GetDistances <- function(x){
          x <- sort(x)
          Distances <- c()
          for(i in 1:length(x)-1){
            Distances <- c(Distances, x[i + 1] - x[i])
          }
          return(Distances)
        }
        
        Ori_Dists <- GetDistances(Ori_ticks)
        
        DistIndex <- which(Ori_Dists<500)
        
        CloseOriIndices <- unique(sort(c(DistIndex, DistIndex+1)))
        
        for(i in 1:length(Ss)){
          if(S == Ss[i]){
            
            if(length(Ori_ticks)>0){
              
              if(length(CloseOriIndices)==0){
                Ori_ticks_distal <- Ori_ticks
                Ori_name_distal <- Ori_chr$name
                color_distal <- colors_ars
              } else {
                Ori_ticks_distal <- Ori_ticks[-CloseOriIndices]
                Ori_name_distal <- Ori_chr$name[-CloseOriIndices]
                color_distal <- colors_ars[-CloseOriIndices]
              }
              
              if(Ori_ticks_distal[1]>300){
                axis(side = 3, at = Ori_ticks_distal, labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                axis(side = 3, at = Ori_ticks_distal, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal)), at = Ori_ticks_distal, side = 3, line=0, col = "purple", cex = 1)
                mtext(rep(intToUtf8(9650), length(Ori_ticks_distal)), at = Ori_ticks_distal, side = 3, line=0.075, col = "yellow", cex = 0.65)
                
                text(x = Ori_ticks_distal, y = grconvertY(x, from = "ndc"), labels = Ori_name_distal, xpd = NA, srt = 0, col = color_distal, cex = 0.9 )
              } else {
                if(length(Ori_ticks_distal)>1){
                  axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=-0.07, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  axis(side = 3, at = Ori_ticks_distal[-1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0, col = "purple", cex = 1)
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[-1])), at = Ori_ticks_distal[-1], side = 3, line=0.075, col = "yellow", cex = 0.65)
                  
                  text(x = Ori_ticks_distal[-1], y = grconvertY(x, from = "ndc"), labels = Ori_name_distal[-1], xpd = NA, srt = 0, col = color_distal[-1], cex = 0.9 )
                  ###
                  axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
                  
                  text(x = Ori_ticks_distal[1], y = grconvertY(x+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
                } else {
                  axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=-(0.07+0.21), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  axis(side = 3, at = Ori_ticks_distal[1], labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                  
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0+2.1, col = "purple", cex = 1)
                  mtext(rep(intToUtf8(9650), length(Ori_ticks_distal[1])), at = Ori_ticks_distal[1], side = 3, line=0.075+2.1, col = "yellow", cex = 0.65)
                  
                  text(x = Ori_ticks_distal[1], y = grconvertY(x+0.025, from = "ndc"), labels = Ori_name_distal[1], xpd = NA, srt = 0, col = color_distal[1], cex = 0.9 )
                }
                
              }
              
            }
            
            if(length(CloseOriIndices)>1){
              Consec_Oris <- split(CloseOriIndices, cumsum(c(1, diff(CloseOriIndices) != 1)))
              
              if(length(Consec_Oris)>=1){
                for(j in 1:length(Consec_Oris)){
                  for(s in 0:(length(Consec_Oris[[j]])-1)){
                    
                    Pos <- Ori_ticks[Consec_Oris[[j]]][s+1]
                    Nam <- Ori_chr$name[Consec_Oris[[j]]][s+1]
                    Col <- colors_ars[Consec_Oris[[j]]][s+1]
                    
                    tckL <- seq(0.05, 0.75, length.out = 5)
                    
                    axis(side = 3, at = Pos, labels = F, line = 0, tick = TRUE, tck=-(tckL[s+1]), lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    axis(side = 3, at = Pos, labels = F, line = 0, tick = TRUE, tck=0.03, lwd.ticks = 2, col.ticks = rgb(128,0,128,alpha=100, maxColorValue=255))
                    
                    lineL <- seq(0, 6.20, length.out = 5)
                    
                    mtext(rep(intToUtf8(9650), length(Pos)), at = Pos, side = 3, line=(lineL[s+1])+0, col = "purple", cex = 0.90)
                    mtext(rep(intToUtf8(9650), length(Pos)), at = Pos, side = 3, line=(lineL[s+1])+0.075, col = "yellow", cex = 0.585)
                    
                    texL <- seq(0, 0.0725, length.out = 5)
                    
                    text(x = Pos, y = grconvertY(texL[s+1]+x, from = "ndc"), labels = Nam, xpd = NA, srt = 0, col = Col, cex = 0.80 )
                  }
                }
              }
            }
            
          }
        }
        
      }
      
      Draw_name(S1, 0.82)
      
      #Draw_name(S2, 0.47)
      
      #Draw_name(S3, 0.31)
    }
  }
  
  
  
  
  #plot x axis
  if(E>seqlengths(Scerevisiae)[k]){
    axis_ticks <-  seq(0, (seqlengths(Scerevisiae)[k] - S)/steps, ((seqlengths(Scerevisiae)[k] - S)/(steps*(seqlengths(Scerevisiae)[k] - S)/10000)))
  } else {
    axis_ticks <-  seq(0, (E - S)/steps, ((E - S)/(steps*(Length_per_Row/10000))))
  }
  axis_ticks[1] <- 1
  axis(side = 1, at = axis_ticks, labels = rep(NA, 11), line = 0, tick = TRUE, tck=0.03, lwd.ticks = 1.5)
  axis(side = 1, at = axis_ticks, labels = rep(NA, 11), line = 0, tick = TRUE, tck=-0.03, lwd.ticks = 1.5)
  title(ylab='readDensity', col="gray", cex.lab=1.25, line = 0, outer = T)
  
  if(SamplType==Sample_4){
    axis(side = 1, at = axis_ticks, labels = round((S/1000)+axis_ticks*(steps/1000)), line = 0, tick = F)
    title(xlab="Chromosomal Coordinates (Kbp)", col="gray", cex.lab=1.25, line = 0, outer = T)
  }
  
  # #plot the watson to crick ratio in a new plot
  # par(new=TRUE)
  # plot(WC.rat, ylim=Ylim_wc, col =  rgb(0,0,205,alpha=90, maxColorValue=255), ann=FALSE, axes=FALSE, type='l', lty = 1, lwd = 1.5, xaxs='i')
  # mtext("log2(watson/crick)", side=4, line=0, outer = T, cex = 0.75)
  #
  # #plot right y axis
  # axis(side = 4, at = yAxis_wc, labels = round(yAxis_wc), line = 0, tick = TRUE, lwd.ticks = 1.5, las = 2, cex.axis = 0.8)
  #
  #
  box("figure", col="forestgreen")
}

S1_Ratio <- cbind.data.frame(chrom=S1_IP_watson$chrom, chromStart=S1_IP_watson$chromStart, chromEnd=S1_IP_watson$chromEnd, wat.score=S1_IP_watson$ratio, cri.score=S1_IP_crick$ratio)
S2_Ratio <- cbind.data.frame(chrom=S2_IP_watson$chrom, chromStart=S2_IP_watson$chromStart, chromEnd=S2_IP_watson$chromEnd, wat.score=S2_IP_watson$ratio, cri.score=S2_IP_crick$ratio)
S3_Ratio <- cbind.data.frame(chrom=S3_IP_watson$chrom, chromStart=S3_IP_watson$chromStart, chromEnd=S3_IP_watson$chromEnd, wat.score=S3_IP_watson$ratio, cri.score=S3_IP_crick$ratio)
S4_Ratio <- cbind.data.frame(chrom=S4_IP_watson$chrom, chromStart=S4_IP_watson$chromStart, chromEnd=S4_IP_watson$chromEnd, wat.score=S4_IP_watson$ratio, cri.score=S4_IP_crick$ratio)

Name <- paste0(unlist(strsplit(Sample_1, "-"))[1])

raster_pdf(file = paste0(Name, "_global_profiles.pdf"), width = 9, height = 11, units = "in", res = 400)
par(oma=c(2,2,2,2))
for(k in 1:16){

  #define plot layout for global profiles
  LayOut.Dims <- function(x, y){
    xx <- c()
    y <- y
    for(i in 1:length(x)){
      xx <- c(xx, rep(x[i], y))
    }
    return(xx)
  }
  Plot.Nums <- c(1:4)
  Cols <- 14

  mat <- matrix(LayOut.Dims(Plot.Nums, Cols),
                length(Plot.Nums), Cols, byrow=TRUE)

  vec <- c(1,2,7,8)
  new_mat <- matrix(0,nrow=length(Plot.Nums)+4,ncol=Cols)
  new_mat[-vec,] <- mat

  fin_mat <- new_mat[rep(1:nrow(new_mat), times = c(2,2,3,3,3,3,2,2)), ]

  layout(fin_mat, c(1,1), c(1,1), TRUE)
  #layout.show(n)


  k <- k

  #define plot intervals by chromosomes
  Length_per_Row <- 100000
  Plotting_Rows <- c(seq(0, seqlengths(Scerevisiae)[k]+Length_per_Row, Length_per_Row))
  Starts <- Plotting_Rows[-length(Plotting_Rows)]
  Ends <- Plotting_Rows[-1]

  #plot global profiles
  for(i in 1:length(Plotting_Rows[-length(Plotting_Rows)])){

    S <- Starts[i]
    E <- Ends[i]

    ChIP_Plot_function(CoverageFile = S1_Ratio, peakFile = S1_Peaks, DataType="strandedRatio", SamplType = Sample_1)
    ChIP_Plot_function(CoverageFile = S2_Ratio, peakFile = S2_Peaks, DataType="strandedRatio", SamplType = Sample_2)
    ChIP_Plot_function(CoverageFile = S3_Ratio, peakFile = S3_Peaks, DataType="strandedRatio", SamplType = Sample_3)
    ChIP_Plot_function(CoverageFile = S4_Ratio, peakFile = S4_Peaks, DataType="strandedRatio", SamplType = Sample_4)
    
  }
}
dev.off()











