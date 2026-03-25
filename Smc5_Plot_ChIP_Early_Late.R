
### ------------ This function works on previously available results from ngsAnalyser-1.1.4 ---

Plot_ChIP_Early_Late <- function(Samplink, AvgWindow, WCwindow){
  
  packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
                "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw",
                "BSgenome.Scerevisiae.UCSC.sacCer3", "shiny")
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  E_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/E_Rep.bed", header = TRUE, quote = "/t")
  L_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/L_Rep.bed", header = TRUE, quote = "/t")
  
  image <- load.image("/Applications/ngsAnalyser.app/Contents/Resources/app/ForkPic_1.jpeg")
  
  colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); E_Ori$mid <- round((E_Ori$chromStart + E_Ori$chromEnd)/2)
  colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "name"); L_Ori$mid <- round((L_Ori$chromStart + L_Ori$chromEnd)/2)
  
  binSize <- 300
  stepSize <- 10
  
  #Samplink <- "/Volumes/rbmData_2.0/Smc5_WT_vs_EQ_time_course/analysis_folders/Smc5-Flag-30UNT-ChIP"
  #Samplink_B <- "/Volumes/rbmData_2.0/Smc5_WT_vs_EQ_time_course/analysis_folders/Smc5-Flag-30UNTBrDU"
  
  #peakFiles <- Sys.glob(paste0("~/Desktop/", Name, "/", "Peaks", "/", "*", ".bed"))
  watsonFiles <-  Sys.glob(paste0(Samplink, "/", "Ratios", "/", "*", "watson.bed"))
  crickFiles <-  Sys.glob(paste0(Samplink, "/", "Ratios", "/", "*", "crick.bed"))
  
  Average_Enrichment <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = TRUE){
    
    #Window <- as.numeric(read.table(peakFiles[6], row.names = 1)["AveragingWindow", ])
    
    AvgWindow <- AvgWindow
    
    #BrDU_Input_Peaks <- read.table(peakFiles[5], header = T)
    
    watsonRatio <- read.table(watsonRatio, header = T)
    crickRatio <- read.table(crickRatio, header = T)
    #PeakList <- read.table(PeakList, header = T)
    
    
    if(Normalise2Input == TRUE){
      V <- 7
    } else {
      V <- 5
    }
    
    ##
    #
    #
    # RemList <- anti_join(BrDU_Input_Peaks, PeakList, by = c("chrom", "oriName", "oriCenter"))
    # NamList <- anti_join(BrDU_Input_Peaks, RemList, by = c("chrom", "oriName", "oriCenter"))
    #
    # NamList <- NamList[!duplicated(NamList$oriName), ]
    #
    # PeakList$BrDUSummit <- NamList$peakSummit
    ###
    
    PeakList$AvBstart <- PeakList$mid - AvgWindow
    PeakList$AvBend <- PeakList$mid + AvgWindow
    
    chrS <- paste0("chr", as.roman(1:16))
    
    #watson
    IP_R <- watsonRatio
    IP_W <- NULL
    for(i in 1:length(chrS)){
      
      OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
      IP_ROs <- IP_R[IP_R$chrom == chrS[i], ]
      
      if(length(OriList_ROs$chrom)==0) next
      
      IP_C <- NULL
      for(y in 1:length(OriList_ROs$chrom)){
        
        IP_Z <- IP_ROs[IP_ROs$chromStart>=OriList_ROs$AvBstart[y] & IP_ROs$chromStart<=OriList_ROs$AvBend[y], ]
        
        if(length(IP_Z[,V])==round(2*AvgWindow/stepSize)){
          Ratios <- IP_Z[,V]
        }
        
        if(length(IP_Z[,V]) < round(2*AvgWindow/stepSize)){
          if(length(1:(length(IP_Z[,V])/2)) < round(2*AvgWindow/stepSize)/2){
            Ratios <- c(rep(0, (round(2*AvgWindow/stepSize))-length(IP_Z$chromStart)), IP_Z[,V] )
          }
          if(length((length(IP_Z[,V])/2+1):length(IP_Z[,V])) < round(2*AvgWindow/stepSize)/2){
            Ratios <- c(IP_Z[,V], rep(0, (round(2*AvgWindow/stepSize))-length(IP_Z$chromStart)) )
          }
        }
        
        if(length(IP_Z[,V]) > round(2*AvgWindow/stepSize)){
          Ratios <- IP_Z[,V][-c((round(2*AvgWindow/stepSize)+1):length(IP_Z[,V]))]
        }
        
        Rat <- matrix(Ratios, ncol = 1)
        
        IP_C <- cbind(IP_C, Rat)
      }
      IP_W <- cbind(IP_W, IP_C)
    }
    IP_W <- as.data.frame(IP_W)
    colnames(IP_W) <- c(1:length(PeakList$chrom))
    
    #crick
    IP_R <- crickRatio
    IP_Cr <- NULL
    for(i in 1:length(chrS)){
      OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
      IP_ROs <- IP_R[IP_R$chrom == chrS[i], ]
      
      if(length(OriList_ROs$chrom)==0) next
      
      IP_C <- NULL
      for(y in 1:length(OriList_ROs$chrom)){
        
        IP_Z <- IP_ROs[IP_ROs$chromStart>=OriList_ROs$AvBstart[y] & IP_ROs$chromStart<=OriList_ROs$AvBend[y], ]
        
        if(length(IP_Z[,V])==round(2*AvgWindow/stepSize)){
          Ratios <- IP_Z[,V]
        }
        
        if(length(IP_Z[,V]) < round(2*AvgWindow/stepSize)){
          if(length(1:(length(IP_Z[,V])/2)) < round(2*AvgWindow/stepSize)/2){
            Ratios <- c(rep(0, (round(2*AvgWindow/stepSize))-length(IP_Z$chromStart)), IP_Z[,V] )
          }
          if(length((length(IP_Z[,V])/2+1):length(IP_Z[,V])) < round(2*AvgWindow/stepSize)/2){
            Ratios <- c(IP_Z[,V], rep(0, (round(2*AvgWindow/stepSize))-length(IP_Z$chromStart)) )
          }
        }
        
        if(length(IP_Z[,V]) > round(2*AvgWindow/stepSize)){
          Ratios <- IP_Z[,V][-c((round(2*AvgWindow/stepSize)+1):length(IP_Z[,V]))]
        }
        
        Rat <- matrix(Ratios, ncol = 1)
        
        IP_C <- cbind(IP_C, Rat)
      }
      IP_Cr <- cbind(IP_Cr, IP_C)
    }
    IP_Cr <- as.data.frame(IP_Cr)
    colnames(IP_Cr) <- c(1:length(PeakList$chrom))
    
    ###
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
    
    watsonAverage <- cbind.data.frame(watson.q25 = rowStat(IP_W)$q.25, watson.median = rowStat(IP_W)$Median, watson.q75 = rowStat(IP_W)$q.75, watson.se = rowStat(IP_W)$Std.err, watson.mean = rowStat(IP_W)$Mean, watson.sd = rowStat(IP_W)$Sd)
    crickAverage <- cbind.data.frame(crick.q25 = rowStat(IP_Cr)$q.25, crick.median = rowStat(IP_Cr)$Median, crick.q75 = rowStat(IP_Cr)$q.75, crick.se = rowStat(IP_Cr)$Std.err, crick.mean = rowStat(IP_Cr)$Mean, crick.sd = rowStat(IP_Cr)$Sd)
    
    TwoStrands <- cbind.data.frame(watsonAverage, crickAverage)
    
    rownames(TwoStrands) <- paste0("bin", 1:length(TwoStrands[,1]))
    
    return(TwoStrands)
  }
  WatsonOverCrick_Average <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = TRUE){
    
    #Window <- as.numeric(read.table(peakFiles[6], row.names = 1)["AveragingWindow", ])
    
    WCwindow <- WCwindow
    
    # BrDU_Input_Peaks <- read.table(peakFiles[5], header = T)
    
    
    Watson <- read.table(watsonRatio, header = T)
    Crick <- read.table(crickRatio, header = T)
    #PeakList <- read.table(PeakList, header = T)
    
    
    if(Normalise2Input == TRUE){
      V <- 7
    } else {
      V <- 5
    }
    
    ##
    
    #
    # RemList <- anti_join(BrDU_Input_Peaks, PeakList, by = c("chrom", "oriName", "oriCenter"))
    # NamList <- anti_join(BrDU_Input_Peaks, RemList, by = c("chrom", "oriName", "oriCenter"))
    #
    # NamList <- NamList[!duplicated(NamList$oriName), ]
    #
    # PeakList$BrDUSummit <- NamList$peakSummit
    ###
    
    PeakList$AvBstart <- PeakList$mid - WCwindow
    PeakList$AvBend <- PeakList$mid + WCwindow
    
    chrS <- paste0("chr", as.roman(1:16))
    
    IP_T <- NULL
    for(i in 1:length(chrS)){
      
      i <- i
      
      OriList_ROs <- PeakList[PeakList$chrom == chrS[i], ]
      Crick_ROs <- Crick[Crick$chrom == chrS[i], ]
      Watson_ROs <- Watson[Watson$chrom == chrS[i], ]
      
      if(length(OriList_ROs$chrom)==0) next
      
      IP_C <- NULL
      
      for(y in 1:length(OriList_ROs$chrom)){
        
        y <- y
        
        Crick_S <- Crick_ROs[Crick_ROs$chromStart>=OriList_ROs$AvBstart[y] & Crick_ROs$chromStart<=OriList_ROs$AvBend[y], ]
        Watson_S <- Watson_ROs[Watson_ROs$chromStart>=OriList_ROs$AvBstart[y] & Watson_ROs$chromStart<=OriList_ROs$AvBend[y], ]
        
        IP_Z <- log2(Watson_S[,V] / Crick_S[,V]); IP_Z[!is.finite(IP_Z)] <- 0
        
        if(length(IP_Z)==round(2*WCwindow/stepSize)){
          Ratios <- IP_Z
        }
        
        if(length(IP_Z) < round(2*WCwindow/stepSize)){
          if(length(1:(length(IP_Z)/2)) < round(2*WCwindow/stepSize)/2){
            Ratios <- c(rep(0, (round(2*WCwindow/stepSize))-length(Crick_S$chromStart)), IP_Z )
          }
          if(length((length(IP_Z)/2+1):length(IP_Z)) < round(2*WCwindow/stepSize)/2){
            Ratios <- c(IP_Z, rep(0, (round(2*WCwindow/stepSize))-length(Crick_S$chromStart)) )
          }
        }
        
        if(length(IP_Z) > round(2*WCwindow/stepSize)){
          Ratios <- IP_Z[-c((round(2*WCwindow/stepSize)+1):length(IP_Z))]
        }
        
        Rat <- matrix(Ratios, ncol = 1)
        
        IP_C <- cbind(IP_C, Rat)
        
      }
      
      IP_T <- cbind(IP_T, IP_C)
      IP_T <- as.data.frame(IP_T)
    }
    colnames(IP_T) <- c(1:length(PeakList$chrom))
    
    ##extract bin stats
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
    
    watsonOvercrickAvg <- cbind.data.frame(q25 = rowStat(IP_T)$q.25, median = rowStat(IP_T)$Median, q75 = rowStat(IP_T)$q.75, se = rowStat(IP_T)$Std.err, mean = rowStat(IP_T)$Mean, sd = rowStat(IP_T)$Sd)
    
    rownames(watsonOvercrickAvg) <- paste0("bin", 1:length(watsonOvercrickAvg[,1]))
    
    return(watsonOvercrickAvg)
    
  }
  Bias_at_individual_peaks <- function(watsonRatio, crickRatio, PeakList, Normalise2Input = F){
    
    watsonRatio <- read.table(watsonRatio, header = T)
    crickRatio <- read.table(crickRatio, header = T)
    #PeakList <- read.table(PeakList, header = T)
    
    
    chrS <- paste0("chr", as.roman(1:16))
    
    IP_T <- NULL
    
    for(i in 1:length(chrS)){
      
      i <- i
      
      IP_ROs_pos <- watsonRatio[watsonRatio$chrom == chrS[i], ]
      IP_ROs_neg <- crickRatio[crickRatio$chrom == chrS[i], ]
      Chr_Peaks <- PeakList[PeakList$chrom == chrS[i], ]
      
      if(length(Chr_Peaks$chrom)==0) next
      
      ###
      
      IP_C <- NULL
      
      for(y in 1:length(Chr_Peaks$chrom)){
        
        y <- y
        
        
        IP_Z_pos_left <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$chromStart[y] & IP_ROs_pos$chromStart<=Chr_Peaks$mid[y], ]
        IP_Z_pos_right <- IP_ROs_pos[IP_ROs_pos$chromStart>=Chr_Peaks$mid[y] & IP_ROs_pos$chromStart<=Chr_Peaks$chromEnd[y], ]
        
        IP_Z_neg_left <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$chromStart[y] & IP_ROs_neg$chromStart<=Chr_Peaks$mid[y], ]
        IP_Z_neg_right <- IP_ROs_neg[IP_ROs_neg$chromStart>=Chr_Peaks$mid[y] & IP_ROs_neg$chromStart<=Chr_Peaks$chromEnd[y], ]
        
        
        
        IP_Z_CL <- sum(IP_Z_pos_left[,5], na.rm = TRUE)
        IP_Z_CR <- sum(IP_Z_pos_right[,5], na.rm = TRUE)
        IP_Z_WL <- sum(IP_Z_neg_left[,5], na.rm = TRUE)
        IP_Z_WR <- sum(IP_Z_neg_right[,5], na.rm = TRUE)
        
        In_Z_CL <- sum(IP_Z_pos_left[,6], na.rm = TRUE)
        In_Z_CR <- sum(IP_Z_pos_right[,6], na.rm = TRUE)
        In_Z_WL <- sum(IP_Z_neg_left[,6], na.rm = TRUE)
        In_Z_WR <- sum(IP_Z_neg_right[,6], na.rm = TRUE)
        
        IP_Z_CL_WR <- IP_Z_CL + IP_Z_WR; IP_Z_CL_WR[is.na(IP_Z_CL_WR)]<-0
        IP_Z_WL_CR <- IP_Z_WL + IP_Z_CR; IP_Z_WL_CR[is.na(IP_Z_WL_CR)]<-0
        
        IP_Z_Pvalue <- binom.test(round(c(IP_Z_CL_WR, IP_Z_WL_CR)))$p.value
        
        In_Z_CL_WR <- In_Z_CL + In_Z_WR; In_Z_CL_WR[is.na(In_Z_CL_WR)]<-0
        In_Z_WL_CR <- In_Z_WL + In_Z_CR; In_Z_WL_CR[is.na(In_Z_WL_CR)]<-0
        
        In_Z_LaLe <- In_Z_CL_WR/In_Z_WL_CR; In_Z_LaLe[!is.finite(In_Z_LaLe)] <- 0
        IP_Z_LaLe <- IP_Z_CL_WR/IP_Z_WL_CR; IP_Z_LaLe[!is.finite(IP_Z_LaLe)] <- 0
        
        if(Normalise2Input==FALSE){
          IP_Z_LagLead <- log2(IP_Z_LaLe); IP_Z_LagLead[!is.finite(IP_Z_LagLead)] <- 0
        } else {
          IP_Z_LagLead <- log2(IP_Z_LaLe/In_Z_LaLe); IP_Z_LagLead[!is.finite(IP_Z_LagLead)] <- 0
        }
        
        IP_Z <- cbind.data.frame(round(IP_Z_CL_WR), round(IP_Z_WL_CR), IP_Z_Pvalue, IP_Z_LagLead)
        IP_C <- rbind.data.frame(IP_C, IP_Z)
      }
      IP_C <- cbind.data.frame(Chr_Peaks$chrom, IP_C)
      colnames(IP_C) <- c('chrom', "Lagg.sum", "Lead.sum", "p_value", "Bias")
      
      status <- rep('null', length(Chr_Peaks$chrom))
      signif <- rep('null', length(Chr_Peaks$chrom))
      status[which(IP_C$Bias>=0)] <- 'lagging_bias'; status[which(IP_C$Bias<0)] <- 'leading_bias'
      signif[which(IP_C$p_value<=10e-6)] <- 'significant'; signif[which(IP_C$p_value>10e-6)] <- 'not_signif'
      
      IP_C$description <- paste0(signif, "_", status)
      
      IP_T <- rbind(IP_T, IP_C)
      IP_T <- as.data.frame(IP_T)
    }
    return(IP_T)
  }
  
  ##
  # LeftLim <- as.numeric(read.table(peakFiles[7], row.names = 1)["Left", ])
  # RightLim <- as.numeric(read.table(peakFiles[7], row.names = 1)["Right", ])
  # LeadingAverage <-  read.table(peakFiles[7], row.names = 1)["LeadSynthesis", ]
  # LaggingAverage <-  read.table(peakFiles[7], row.names = 1)["LaggSynthesis", ]
  # binSize <- as.numeric(read.table(peakFiles[7], row.names = 1)["bin", ])
  # stepSize <- as.numeric(read.table(peakFiles[7], row.names = 1)["slide", ])
  # AveragingWindow <- as.numeric(read.table(peakFiles[7], row.names = 1)["AveragingWindow", ])
  
  SavePlot <- paste0(dirname(dirname(Samplink)), "/", "Early_Late_Origins_Plots")
  Name <- basename(Samplink)
  
  PlotEnrichments <- function(DataFile, PlotHeader){
    
    Watson <- round(as.numeric(smooth.spline(1:length(DataFile$watson.median), DataFile$watson.median)$y), 2)
    Crick <- round(as.numeric(smooth.spline(1:length(DataFile$crick.median), DataFile$crick.median)$y), 2)
    
    Wat25 <- round(as.numeric(smooth.spline(1:length(DataFile$watson.q25), DataFile$watson.q25)$y), 2)
    Wat75 <- round(as.numeric(smooth.spline(1:length(DataFile$watson.q75), DataFile$watson.q75)$y), 2)
    
    Cri25 <- round(as.numeric(smooth.spline(1:length(DataFile$crick.q25), DataFile$crick.q25)$y), 2)
    Cri75 <- round(as.numeric(smooth.spline(1:length(DataFile$crick.q75), DataFile$crick.q75)$y), 2)
    
    Y <- max(round(abs(range(c(Wat75, Cri75*(-1))))+0.5))
    
    plot(Watson,
         ylim = c(-Y, Y),
         main = PlotHeader,
         ylab = "Average Enrichment", cex.main=0.8, xlab = "Distance from Ori Center (Kbp)",
         xaxt = "n", col = 'brown3', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs = 'i', yaxs = 'i')
    
    lines(Wat25, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
    lines(Wat75, lwd=0.5, col = adjustcolor("red", alpha.f = 0.2))
    
    polygon(x = c(1:length(Watson), rev(1:length(Watson))),
            y = c(Wat25, rev(Wat75)),
            col = adjustcolor("red", alpha.f = 0.2), border = NA)
    
    lines((Crick)*(-1), lwd=2, col = 'cornflowerblue', type = 'l')
    lines((Cri25)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
    lines((Cri75)*(-1), lwd=0.5, col = adjustcolor("cornflowerblue", alpha.f = 0.2))
    
    polygon(x = c(1:length(Crick*(-1)), rev(1:length(Crick*(-1)))),
            y = c(Cri25*(-1), rev(Cri75*(-1))),
            col = adjustcolor("cornflowerblue", alpha.f = 0.2), border = NA)
    
    # rect(0, -Y, (length(Watson))/2 - (LeftLim/stepSize), Y, col = adjustcolor("grey", alpha.f = 0.65), border = NA)
    # rect((length(Watson))/2 + (RightLim/stepSize), -Y, length(Watson), Y,
    #      col = adjustcolor("grey", alpha.f = 0.65), border = NA)
    
    text((AvgWindow)*2-500, Y-1, labels = "Watson", cex = 0.9, col = 'brown3')
    text((AvgWindow)*2-500, -Y+1, labels = "Crick", cex = 0.9, col = 'cornflowerblue')
    
    abline(h=0,lwd=0.4); abline(v=(length(Watson))/2,lwd=0.4)
    axisLabels <- seq(-AvgWindow,
                      +AvgWindow,
                      length.out = 9)
    axisLabels[c(2,4,6,8)] <- NA
    At <- (AvgWindow/stepSize)*seq(0,2,0.25); At[1] <- 1
    axis(1, at=At, labels = signif(axisLabels/1000, 2))
    
  }
  PlotAverages <- function(DataFile, PlotHeader){
    
    
    Med <- round(as.numeric(smooth.spline(1:length(DataFile$median), DataFile$median)$y), 2)
    q25 <- round(as.numeric(smooth.spline(1:length(DataFile$q25), DataFile$q25)$y), 2)
    q75 <- round(as.numeric(smooth.spline(1:length(DataFile$q75), DataFile$q75)$y), 2)
    
    Y <- max(round(abs(range(q75))+0.5))
    
    # Y <- 1
    
    plot(Med,
         ylim = c(-Y, +Y),
         main = PlotHeader,
         ylab = "log2 watson/crick", cex.main=0.8, xlab = "Distance from Ori Center (Kbp)",
         xaxt = "n", col = 'blue', type = 'l', lwd = 2, bty = 'n', las = 2, xaxs='i', yaxs='i')
    
    polygon(x = c(1:length(Med), rev(1:length(Med))),
            y = c(q25, rev(q75)),
            col = adjustcolor("blue", alpha.f = 0.2), border = NA)
    
    abline(h=0, lwd=0.4); abline(v=(length(Med))/2,lwd=0.4)
    axisLabels <- seq(-WCwindow,
                      +WCwindow,
                      length.out = 9)
    axisLabels[c(2,4,6,8)] <- NA
    At <- (WCwindow/stepSize)*seq(0,2,0.25); At[1] <- 1
    axis(1, at=At, labels = signif(axisLabels/1000, 2))
    
  }
  PlotIdBias <- function(DataFile, PlotHeader){
    
    DataFile$Bias[!is.finite(DataFile$Bias)] <- 0
    DataFile$Bias[is.na(DataFile$Bias)] <- 0
    
    boxplot(DataFile$Bias, ylim = c(-2,2), las = 2, ylab = 'log2 lagging/leading',
            cex.main=0.8, xlab = " ", bty = 'n', main=PlotHeader, col=adjustcolor("grey", alpha.f = 0.25), lwd=0.5)
    
    if(length(DataFile$Bias[which(DataFile$p_value > 10e-6)])>0){
      spreadPoints(values=DataFile$Bias[which(DataFile$p_value > 10e-6)], position=1.0, pointCex=0.65, col="blue", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
    }
    if(length(DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)])>0){
      spreadPoints(values=DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)], position=1.0, pointCex=0.65, col="orange", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
    }
    if(length(DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)])>0){
      spreadPoints(values=DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)], position=1.0, pointCex=0.65, col="green", pch=19, alpha=0.5, plotOutliers=T, fitToBoxWidth=TRUE, xpd=FALSE, widthCex=1)
    }
    
    
    Lagg <- length(DataFile$Bias[which(DataFile$Bias >=0 & DataFile$p_value <= 10e-6)])
    Lead <- length(DataFile$Bias[which(DataFile$Bias < 0 & DataFile$p_value <= 10e-6)])
    Indt <- length(DataFile$Bias[which(DataFile$p_value > 10e-6)])
    
    legend("bottomright", legend = c(paste0("lagg", " (", Lagg, ")"),
                                     paste0("lead", " (", Lead, ")"),
                                     paste0("inde", " (", Indt, ")")),
           col = c("orange", "green", "blue"), pch = 19, pt.cex=0.65, bty = "n", cex = 0.65)
    
    
    InName <- PlotHeader
    
    #calculate p value
    #Decision Tree
    biasP <- binom.test(round(c(Lagg+Lead, Indt)))$p.value
    
    if(InName=="ChIP_Input" || InName=="ChIP"){
      
      if(biasP <= 10e-3 & (Lagg+Lead) > Indt){
        
        biasQ <- binom.test(round(c(Lagg, Lead)))$p.value
        
        if(biasQ <= 10e-6){
          if(Lagg > Lead){
            conc <- paste0("Strong lagging bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
          }
          if(Lead > Lagg){
            conc <- paste0("Strong leading bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
          }
        }
        
        if(biasQ > 10e-6 & biasQ <= 10e-4){
          if(Lagg > Lead){
            conc <- paste0("Weak lagging bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
          }
          if(Lead > Lagg){
            conc <- paste0("Weak leading bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
          }
        }
        
        if(biasQ > 10e-4  & biasQ <= 10e-2){
          if(Lagg > Lead){
            conc <- paste0("Very weak lagging bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
          }
          if(Lead > Lagg){
            conc <- paste0("Very weak leading bias detected", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
          }
        }
        
        if(biasQ > 10e-2){
          conc <- paste0("Factor may bind one or the other strand", "\n","p-value (Lagg~Lead): ", signif(biasQ, 3))
        }
      }
      else
      {
        conc <- paste0("No decisive bias pattern detected", "\n","too many origins show indeterminable bias")
      }
    }
    
    mtext(conc, side = 1, line = 2, cex = 0.65)
  }
  
  pdf(paste0(SavePlot, "/", Name, "_early_late.pdf"), width = 9, height = 12)
  
  par(oma=c(0,0,0,0))
  
  PlotMat <- {matrix(c(     0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,
                            
                            0,0,0,0,2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,
                            0,0,0,0,2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,
                            0,0,0,0,2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,
                            0,0,0,0,2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,
                            0,0,0,0,2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,
                            0,0,0,0,2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,
                            
                            0,0,0,0,4,4,4,4,4,4,5,5,5,5,5,5,0,0,0,0,
                            0,0,0,0,4,4,4,4,4,4,5,5,5,5,5,5,0,0,0,0,
                            0,0,0,0,4,4,4,4,4,4,5,5,5,5,5,5,0,0,0,0,
                            0,0,0,0,4,4,4,4,4,4,5,5,5,5,5,5,0,0,0,0,
                            0,0,0,0,4,4,4,4,4,4,5,5,5,5,5,5,0,0,0,0,
                            0,0,0,0,4,4,4,4,4,4,5,5,5,5,5,5,0,0,0,0,
                            
                            0,0,0,0,6,6,6,6,6,6,7,7,7,7,7,7,0,0,0,0,
                            0,0,0,0,6,6,6,6,6,6,7,7,7,7,7,7,0,0,0,0,
                            0,0,0,0,6,6,6,6,6,6,7,7,7,7,7,7,0,0,0,0,
                            0,0,0,0,6,6,6,6,6,6,7,7,7,7,7,7,0,0,0,0,
                            0,0,0,0,6,6,6,6,6,6,7,7,7,7,7,7,0,0,0,0,
                            0,0,0,0,6,6,6,6,6,6,7,7,7,7,7,7,0,0,0,0),
                     
                     23,20,byrow=TRUE)}
  
  #early origins
  
  ChIP_AvE_peaks <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = F)
  ChIP_Input_AvE_peaks <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = T)
  
  ChIP_WoC_peaks <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = F)
  ChIP_Input_WoC_peaks <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = T)
  
  ChIP_IdB_peaks <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = F)
  ChIP_Input_IdB_peaks <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = E_Ori, Normalise2Input = T)
  
  layout(PlotMat, c(1,1), c(1,1), TRUE)
  plot(image, axes = FALSE, main = paste0("Early_ROs", " ", "(", dim(E_Ori)[1], ")"))
  PlotEnrichments(ChIP_AvE_peaks, "ChIP"); PlotEnrichments(ChIP_Input_AvE_peaks, "ChIP_Input")
  PlotAverages(ChIP_WoC_peaks, "ChIP"); PlotAverages(ChIP_Input_WoC_peaks, "ChIP_Input")
  PlotIdBias(ChIP_IdB_peaks, "ChIP"); PlotIdBias(ChIP_Input_IdB_peaks, "ChIP_Input")
  
  #late origins
  
  ChIP_AvE_peaks <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = F)
  ChIP_Input_AvE_peaks <- Average_Enrichment(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = T)
  
  ChIP_WoC_peaks <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = F)
  ChIP_Input_WoC_peaks <- WatsonOverCrick_Average(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = T)
  
  ChIP_IdB_peaks <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = F)
  ChIP_Input_IdB_peaks <- Bias_at_individual_peaks(watsonRatio=watsonFiles[1], crickRatio=crickFiles[1], PeakList = L_Ori, Normalise2Input = T)
  
  layout(PlotMat, c(1,1), c(1,1), TRUE)
  plot(image, axes = FALSE, main = paste0("Late_ROs", " ", "(", dim(L_Ori)[1], ")"))
  PlotEnrichments(ChIP_AvE_peaks, "ChIP"); PlotEnrichments(ChIP_Input_AvE_peaks, "ChIP_Input")
  PlotAverages(ChIP_WoC_peaks, "ChIP"); PlotAverages(ChIP_Input_WoC_peaks, "ChIP_Input")
  PlotIdBias(ChIP_IdB_peaks, "ChIP"); PlotIdBias(ChIP_Input_IdB_peaks, "ChIP_Input")
  
  
  dev.off()
  
}

### example run

Plot_ChIP_Early_Late(
  Samplink = "/Volumes/rbmData_2.0/Smc5_WT_vs_EQ_time_course/analysis_folders/Smc5-Flag-30UNT-ChIP",
  AvgWindow = 3000,
  WCwindow = 3000)