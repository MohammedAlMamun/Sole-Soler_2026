

ChIPseq_Primary_Analysis <- function(  Input_R1 = "/full/path/to/file_R1.fastq.gz", 
                                       Input_R2 = "/full/path/to/file_R2.fastq.gz", 
                                       
                                       ChIP_R1 = "/full/path/to/file_R1.fastq.gz", 
                                       ChIP_R2 = "/full/path/to/file_R2.fastq.gz",
                                       
                                       Alignment = "generic",  # "generic" (bowtie2) / "malign" (Rsubread multi-aligner) / "mrdna" (Rsubread multi-aligner at rDNA locus) 
                                       
                                       ExpTitle = "Smc5-trial-ChIP", 
                                       Directory = "None", 
                                       slidingWindow = "YES" ) { 
  
  ## load packages
  
  packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
                "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw", "shiny",
                "BSgenome.Scerevisiae.UCSC.sacCer3", "Rsubread", "GenomicAlignments",
                "IRanges", "readxl", "data.table", "ORFik")
  
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  All_Ori_Link <- "/Applications/ngsAnalyser.app/Contents/Resources/app/OriginList_Full.bed"
  E_Ori_Link <- "/Applications/ngsAnalyser.app/Contents/Resources/app/E_Rep.bed"
  L_Ori_Link <- "/Applications/ngsAnalyser.app/Contents/Resources/app/L_Rep.bed"
  
  # Sequencing Alignment & Binned Coverage Calculation
  
  #
  useDef <- function(a,d) ifelse(isTruthy(a), a,d)
  
  ExpTitle = useDef(ExpTitle, "None")
  
  
  if(ExpTitle == "None"){
    Pro_1 <- unlist(strsplit(basename(Input_R1), split='_', fixed=TRUE))[[1]]
  } else {
    Pro_1 <- ExpTitle
  }
  
  #
  message(paste0("✅ Experiment: ", Pro_1))
  #
  
  #
  Directory = useDef(Directory, "None")
  
  if(Directory == "None"){
    dir <- "~/Desktop/"
  } else {
    dir <- paste0(Directory, "/")
  }
  #
  
  suppressWarnings(dir.create(paste0(dir, Pro_1)))  
  
  ## Quality check of fastqs'
  # #
  # message("⏳ Running QC ...")
  # #
  # if(!file.exists(paste0(dir, Pro_1, "/", Pro_1, "_", "QR", ".html"))){
  #   
  #   fls = c(Input_R1, Input_R2, ChIP_R1, ChIP_R2)
  #   
  #   names(fls) = sub(".fastq", "", basename(fls))
  #   
  #   qas = lapply(seq_along(fls),
  #                function(i, fls) qa(readFastq(fls[i]), names(fls)[i]),
  #                fls)
  #   qa = do.call(rbind, qas)
  #   rpt = report(qa, dest = paste0(dir, Pro_1, "/", Pro_1, "_", "QR", ".html"))
  #   
  # }
  # 
  ## Run alignment
  #
  message("✅ Reference yeast genome : S288C")
  message("⏳ Running alignments...")
  #
  if(Alignment == "generic"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Bam"))){
      
      RunAlignment_bowtie2 <- function(File_R1, File_R2, SampName){
        
        message(paste0("➤ Running alignment for ", Pro_1, "_", SampName))
        
        tempdir(check = TRUE)
        
        Sam <- tempfile(fileext = ".sam")
        Bam <- tempfile(fileext = ".bam")
        nmCollate <- tempfile(fileext = ".bam")
        fixMat <- tempfile(fileext = ".bam")
        SrtBam <- tempfile(fileext = ".bam")
        
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Bam")))
        
        AlnLog <- paste0(dir, Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".log")
        SFBam <- paste0(dir, Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome for the alignment of sequenced data
        ref_index <- "~/Desktop/chip_support/genome_index_ref/S288C_Ref"
        
        #following commands will run the alignemnt, check quality, sort, filter and index the resultant bam file 
        
        system(sprintf("(/Applications/ngsAnalyser.app/Contents/Resources/app/bowtie2-2.4.4-macos-x86_64/bowtie2 -p 8  --no-discordant --fr -x %s -1 %s -2 %s -S %s) 2> %s", 
                       ref_index, File_R1, File_R2, Sam, AlnLog))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -bS -@ 15 -q 30 -f 2 %s > %s", Sam, Bam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools collate -@ 15 -o %s %s", nmCollate, Bam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools fixmate -@ 15 -m %s %s", nmCollate, fixMat))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools sort -l 9 -@ 15 -m 1024M  -O bam -o %s %s", SrtBam, fixMat))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools markdup -@ 15 %s %s", SrtBam, SFBam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools index -@ 15 %s", SFBam))
        
        unlink(c(Sam, Bam, nmCollate, fixMat, SrtBam), recursive = T, force = T)
        
      }
      
      RunAlignment_bowtie2(Input_R1, Input_R2, "Input")
      RunAlignment_bowtie2(ChIP_R1, ChIP_R2, "ChIP")
      
    }
    
  }
  
  if(Alignment == "malign"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Bam_ma"))){
      
      RunAlignment_subread_malign <- function(File_R1, File_R2, SampName){
        
        message(paste0("➤ Running Rsubread-based multiple-alignment for ", Pro_1, "_", SampName))
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Bam_ma")))
        
        resultBam <- paste0(dir, Pro_1, "/", "Bam_ma", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome
        
        ref_index <- "~/Desktop/chip_support/genome_index_ref/subread_ref"
        
        #run alignment
        
        align(index = ref_index,
              readfile1 = File_R1,
              readfile2 = File_R2,
              output_format = "BAM",
              output_file = resultBam,
              color2base = F,
              type = "dna",
              unique = FALSE,
              nBestLocations = 16,
              nthreads = 8,
              sortReadsByCoordinates = TRUE)
        
      }
      
      RunAlignment_subread_malign(Input_R1, Input_R2, "Input")
      RunAlignment_subread_malign(ChIP_R1, ChIP_R2, "ChIP")
      
    }
    
  }
  
  if(Alignment == "mrdna"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Bam_ma_rdna"))){
      
      RunAlignment_subread_malign_rDNA <- function(File_R1, File_R2, SampName){
        
        message(paste0("➤ Running Rsubread-based multiple-alignment at rDNA for ", Pro_1, "_", SampName))
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Bam_ma_rdna")))
        
        resultBam <- paste0(dir, Pro_1, "/", "Bam_ma_rdna", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome
        
        ref_index <- "~/Desktop/chip_support/genome_index_ref/rDNA_2nts1"
        
        #run alignment
        
        align(index = ref_index,
              readfile1 = File_R1,
              readfile2 = File_R2,
              output_format = "BAM",
              output_file = resultBam,
              color2base = F,
              type = "dna",
              unique = FALSE,
              nBestLocations = 16,
              nthreads = 8,
              sortReadsByCoordinates = TRUE)
        
      }
      
      RunAlignment_subread_malign_rDNA(Input_R1, Input_R2, "Input_rDNA")
      RunAlignment_subread_malign_rDNA(ChIP_R1, ChIP_R2, "ChIP_rDNA")
      
      RunAlignment_subread_malign <- function(File_R1, File_R2, SampName){
        
        message(paste0("➤ Running Rsubread-based multiple-alignment for ", Pro_1, "_", SampName))
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Bam_ma_rdna")))
        
        resultBam <- paste0(dir, Pro_1, "/", "Bam_ma_rdna", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome
        
        ref_index <- "~/Desktop/chip_support/genome_index_ref/subread_ref"
        
        #run alignment
        
        align(index = ref_index,
              readfile1 = File_R1,
              readfile2 = File_R2,
              output_format = "BAM",
              output_file = resultBam,
              color2base = F,
              type = "dna",
              unique = FALSE,
              nBestLocations = 16,
              nthreads = 8,
              sortReadsByCoordinates = TRUE)
        
      }
      
      RunAlignment_subread_malign(Input_R1, Input_R2, "Input_ma")
      RunAlignment_subread_malign(ChIP_R1, ChIP_R2, "ChIP_ma")
      
    }
    
  }
  
  #
  ## Calculate genome-wide binned coverage
  #
  message("⏳ Calculating read coverage ...")
  #
  if(Alignment == "generic"){ BamFolder = "Bam" }
  if(Alignment == "malign"){ BamFolder = "Bam_ma" }
  if(Alignment == "mrdna"){ BamFolder = "Bam_ma_rdna" }
  #
  
  if(Alignment == "generic"){
    
    bamFiles <- c(paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_ChIP.bam"),
                  paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input.bam"))
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Coverage"))){
      
      #coverage
      BamCoverage <- function(bamFile, binSize = 300, stepSize = 10, slidingWindow = "YES", byReads_5p = T){
        
        SampleName = tools::file_path_sans_ext(basename(bamFile))
        Pro_2 = substring(SampleName, nchar(Pro_1)+2)
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Coverage"))) 
        
        tempdir(check = TRUE)
        
        GenomFile <- tempfile(fileext = ".txt")
        binFile <- tempfile(fileext = ".bed")
        
        command_1 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools idxstats %s | awk 'BEGIN {OFS=\"\\t\"} {if ($2>0) print ($1,$2)}' >  %s"
        system(sprintf(command_1, bamFile, GenomFile))
        
        if(slidingWindow=="YES"){
          command_2 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools makewindows -g %s -w %s -s %s | sort -k1,1V -k2,2n > %s"
          system(sprintf(command_2, GenomFile, binSize, stepSize, binFile))
        } else {
          command_2 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools makewindows -g %s -w %s | sort -k1,1V -k2,2n > %s"
          system(sprintf(command_2, GenomFile, binSize, binFile))
        }
        
        pncFiles_watson <- tempfile(fileext = ".bed")
        pncFiles_crick <- tempfile(fileext = ".bed")
        
        if(byReads_5p == TRUE){
          
          message(paste0("➤ Calculating coverage with 5' ends of first-mate reads for", " ", 
                         tools::file_path_sans_ext(basename(bamFile)) ))
          
          # calculate coverage at watson strand by 5' end of the first mate reads
          command_3 <- paste0(
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s |",
            "grep -v XS:i: |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -@ 8 -b - |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools genomecov -ibam stdin -strand + -d -5 |",
            "awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' | sort -k1,1V -k2,2n > %s"
          )
          
          # calculate coverage at crick strand by 5' end of the first mate reads
          command_4 <- paste0(
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s |",
            "grep -v XS:i: |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -@ 8 -b - |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools genomecov -ibam stdin -strand - -d -5 |",
            "awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' | sort -k1,1V -k2,2n > %s"
          )
        }
        
        #
        # Watson (+)
        system(sprintf(command_3, binFile, bamFile,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"),
                       pncFiles_watson))
        
        # Crick (-)
        system(sprintf(command_4, binFile, bamFile,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"),
                       pncFiles_crick))
        
        #
        # sum the counts per bin for watson and crick separately and store in finFiles
        
        finFiles_watson <- paste0(dir, Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "watson.bed")
        finFiles_crick <- paste0(dir, Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "crick.bed")
        
        command_5 <- paste(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools map",
          "-a %s -b %s -null 0 -o sum |",
          "awk 'BEGIN {OFS=\"\\t\"} {if ($4>=0) print $1,$2,$3,\"%s\",$4}' > %s"
        )
        
        # Watson (+)
        system(sprintf(command_5, binFile, pncFiles_watson,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"),
                       finFiles_watson))
        # Crick (-)
        system(sprintf(command_5, binFile, pncFiles_crick,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"),
                       finFiles_crick))
        
        #
        
        unlink(c(GenomFile, binFile, pncFiles_watson, pncFiles_crick), recursive = T, force = T)
        
      }
      
      for(i in 1:length(bamFiles)){
        BamCoverage(bamFile = bamFiles[i])
      }
      
    }
    
  }
  
  if(Alignment == "malign"){
    
    bamFiles <- c(paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_ChIP.bam"),
                  paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input.bam"))
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Coverage_ma"))){
      
      BamCoverage <- function(bamFile, binSize=300, stepSize=10, byReads_5p=TRUE){
        
        SampleName = tools::file_path_sans_ext(basename(bamFile))
        Pro_2 = substring(SampleName, nchar(Pro_1)+2)
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Coverage_ma")))
        
        # whole-genome sliding windows
        Windows = tileGenome(seqlengths(Scerevisiae)[-17], tilewidth=stepSize)
        Windows = unlist(Windows)
        Windows = suppressWarnings(trim(IRanges::resize(Windows, width=binSize)))
        
        Counts_watson = numeric(length(Windows))
        Counts_crick = numeric(length(Windows))
        
        if(byReads_5p == TRUE){
          
          message(paste0("➤ Calculating multiple-alignment coverage using 5' ends of first-mate reads for ",
                         tools::file_path_sans_ext(basename(bamFile))))
          
        } else {
          
          message(paste0("➤ Calculating multiple-alignment coverage using fragment midpoints for ",
                         tools::file_path_sans_ext(basename(bamFile))))
        }
        
        for(i in 1:16){
          
          paramS = ScanBamParam(flag=scanBamFlag(isProperPair=TRUE, isUnmappedQuery=FALSE,
                                                 hasUnmappedMate=FALSE),
                                what=c("qname", "mapq", "isize", "seq"),
                                tag=c("NH", "NM"),
                                tagFilter=list(NH=i, NM=c(0:2)))
          
          AlignmentPairs = suppressWarnings(readGAlignmentPairs(bamFile, param=paramS, strandMode=1))
          
          if(length(AlignmentPairs) == 0) next
          
          if(byReads_5p == TRUE){
            
            # use the strand-specific 5' end of the first mate
            CoveragePoints = granges(GenomicAlignments::first(AlignmentPairs))
            
            fivePrimePosition = ifelse(as.character(strand(CoveragePoints)) == "+",
                                       start(CoveragePoints), end(CoveragePoints))
            
            start(CoveragePoints) = fivePrimePosition
            end(CoveragePoints) = fivePrimePosition
            
          } else {
            
            # use the midpoint of the complete paired-end fragment
            CoveragePoints = granges(AlignmentPairs)
            
            midpoint = round((start(CoveragePoints) + end(CoveragePoints))/2)
            
            start(CoveragePoints) = midpoint
            end(CoveragePoints) = midpoint
          }
          
          # fractional contribution of multiple alignments
          mcols(CoveragePoints)$weight = 1/i
          
          CoveragePoints_watson = CoveragePoints[as.character(strand(CoveragePoints)) == "+"]
          CoveragePoints_crick = CoveragePoints[as.character(strand(CoveragePoints)) == "-"]
          
          if(length(CoveragePoints_watson) > 0){
            
            Counts_watson = Counts_watson + countOverlapsW(Windows, CoveragePoints_watson,
                                                           weight="weight")
          }
          
          if(length(CoveragePoints_crick) > 0){
            
            Counts_crick = Counts_crick + countOverlapsW(Windows, CoveragePoints_crick,
                                                         weight="weight")
          }
        }
        
        finFile_watson = paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_", Pro_2, "_watson.bed")
        finFile_crick = paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_", Pro_2, "_crick.bed")
        
        bed_watson = data.frame(seqnames=as.character(seqnames(Windows)),
                                start=start(Windows)-1L,
                                end=end(Windows),
                                name=paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"),
                                score=Counts_watson)
        
        bed_crick = data.frame(seqnames=as.character(seqnames(Windows)),
                               start=start(Windows)-1L,
                               end=end(Windows),
                               name=paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"),
                               score=Counts_crick)
        
        write.table(bed_watson, file=finFile_watson, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)
        
        write.table(bed_crick, file=finFile_crick, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)
      }
      
      for(i in 1:length(bamFiles)){
        BamCoverage(bamFile = bamFiles[i])
      }
      
    }
    
  }
  
  if(Alignment == "mrdna"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Coverage_ma_rdna"))){
      
      BamCoverage <- function(bamFile, binSize=300, stepSize=10, byReads_5p=TRUE, rDNA=FALSE){
        
        SampleName = tools::file_path_sans_ext(basename(bamFile))
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Coverage_ma_rdna")))
        
        if(rDNA == FALSE){
          
          # whole-genome sliding windows
          Windows = tileGenome(seqlengths(Scerevisiae)[-17], tilewidth=stepSize)
          Windows = unlist(Windows)
          Windows = suppressWarnings(trim(IRanges::resize(Windows, width=binSize)))
        }
        
        if(rDNA == TRUE){
          
          # rDNA-reference sliding windows
          gen = Seqinfo(seqnames="Scer_2xrDNA_unit", seqlengths=18274,
                        isCircular=TRUE, genome="Scer_2xrDNA_unit")
          
          Windows = tileGenome(gen, tilewidth=stepSize)
          Windows = unlist(Windows)
          Windows = suppressWarnings(trim(IRanges::resize(Windows, width=binSize)))
        }
        
        Counts_watson = numeric(length(Windows))
        Counts_crick = numeric(length(Windows))
        
        if(byReads_5p == TRUE){
          
          message(paste0("➤ Calculating multiple-alignment coverage using 5' ends of first-mate reads for ",
                         SampleName))
          
        } else {
          
          message(paste0("➤ Calculating multiple-alignment coverage using fragment midpoints for ",
                         SampleName))
        }
        
        for(i in 1:16){
          
          paramS = ScanBamParam(flag=scanBamFlag(isProperPair=TRUE, isUnmappedQuery=FALSE,
                                                 hasUnmappedMate=FALSE),
                                what=c("qname", "mapq", "isize", "seq"),
                                tag=c("NH", "NM"),
                                tagFilter=list(NH=i, NM=c(0:2)))
          
          AlignmentPairs = suppressWarnings(readGAlignmentPairs(bamFile, param=paramS, strandMode=1))
          
          if(length(AlignmentPairs) == 0) next
          
          if(byReads_5p == TRUE){
            
            # use the strand-specific 5' end of the first mate
            CoveragePoints = granges(GenomicAlignments::first(AlignmentPairs))
            
            fivePrimePosition = ifelse(as.character(strand(CoveragePoints)) == "+",
                                       start(CoveragePoints), end(CoveragePoints))
            
            start(CoveragePoints) = fivePrimePosition
            end(CoveragePoints) = fivePrimePosition
            
          } else {
            
            # use the midpoint of the complete paired-end fragment
            CoveragePoints = granges(AlignmentPairs)
            
            midpoint = round((start(CoveragePoints) + end(CoveragePoints))/2)
            
            start(CoveragePoints) = midpoint
            end(CoveragePoints) = midpoint
          }
          
          # fractional contribution of multiple alignments
          mcols(CoveragePoints)$weight = 1/i
          
          CoveragePoints_watson = CoveragePoints[as.character(strand(CoveragePoints)) == "+"]
          CoveragePoints_crick = CoveragePoints[as.character(strand(CoveragePoints)) == "-"]
          
          if(length(CoveragePoints_watson) > 0){
            
            Counts_watson = Counts_watson + countOverlapsW(Windows, CoveragePoints_watson,
                                                           weight="weight")
          }
          
          if(length(CoveragePoints_crick) > 0){
            
            Counts_crick = Counts_crick + countOverlapsW(Windows, CoveragePoints_crick,
                                                         weight="weight")
          }
        }
        
        finFile_watson = paste0(dir, Pro_1, "/Coverage_ma_rdna/", SampleName, "_watson.bed")
        finFile_crick = paste0(dir, Pro_1, "/Coverage_ma_rdna/", SampleName, "_crick.bed")
        
        bed_watson = data.frame(seqnames=as.character(seqnames(Windows)),
                                start=start(Windows)-1L,
                                end=end(Windows),
                                name=paste0(SampleName, "_watson"),
                                score=Counts_watson)
        
        bed_crick = data.frame(seqnames=as.character(seqnames(Windows)),
                               start=start(Windows)-1L,
                               end=end(Windows),
                               name=paste0(SampleName, "_crick"),
                               score=Counts_crick)
        
        write.table(bed_watson, file=finFile_watson, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)
        
        write.table(bed_crick, file=finFile_crick, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)
      }
      
      BamCoverage(bamFile=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_ChIP_rDNA.bam"), rDNA=TRUE)
      BamCoverage(bamFile=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input_rDNA.bam"), rDNA=TRUE)
      BamCoverage(bamFile=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_ChIP_ma.bam"), rDNA=FALSE)
      BamCoverage(bamFile=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input_ma.bam"), rDNA=FALSE)
      
    }
    
  }
  
  
  ## Define and process peaks
  #
  message("⏳ Processing ChIP peaks ...")
  #
  
  if(Alignment == "generic"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Peaks"))){
      
      ChIP_Peaks_Processing <- function(IPBam, InBam){
        
        message(paste0("➤ Processing ChIP peaks for ", Pro_1))
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Peaks")))
        
        tempdir(check=TRUE)
        
        IP = tempfile(fileext=".bed")
        Input = tempfile(fileext=".bed")
        peakFile = tempfile(fileext=".bed")
        outDir = tempfile(pattern="MACS2_")
        
        suppressWarnings(dir.create(outDir))
        
        # convert ChIP and Input BAM files to BED
        command_1 = "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s > %s"
        
        system(sprintf(command_1, IPBam, IP))
        system(sprintf(command_1, InBam, Input))
        
        # call genome-wide ChIP peaks
        command_2 = "macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n Peak --outdir %s 2> /dev/null"
        
        system(sprintf(command_2, IP, Input, outDir))
        
        allPeaks = read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#")
        
        write.table(allPeaks, file=peakFile, quote=FALSE, row.names=FALSE,
                    sep="\t", col.names=FALSE)
        
        ColHeads = "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\""
        
        # peaks associated with all known origins
        command_3 = paste0(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools intersect ",
          "-wa -wb -a %s -b %s | ",
          "awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'"
        )
        
        Peaks_at_all_Origins = read.table(
          pipe(sprintf(command_3, peakFile, All_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with early origins
        Peaks_at_early_Origins = read.table(
          pipe(sprintf(command_3, peakFile, E_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with late origins
        Peaks_at_late_Origins = read.table(
          pipe(sprintf(command_3, peakFile, L_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # standardize the MACS2 peak columns
        allPeaks = dplyr::rename(allPeaks,
                                 chrom=chr,
                                 peakStart=start,
                                 peakEnd=end,
                                 peakLength=length,
                                 peakSummit=abs_summit)
        
        # peaks not associated with known origins
        Peaks_at_no_Origins = dplyr::anti_join(
          allPeaks,
          Peaks_at_all_Origins,
          by=c("chrom", "peakStart", "peakEnd")
        )
        
        # remove duplicate peak entries
        Peaks_at_all_Origins = Peaks_at_all_Origins[ !duplicated(Peaks_at_all_Origins$peakSummit), ]
        
        Peaks_at_early_Origins = Peaks_at_early_Origins[ !duplicated(Peaks_at_early_Origins$peakSummit), ]
        
        Peaks_at_late_Origins = Peaks_at_late_Origins[ !duplicated(Peaks_at_late_Origins$peakSummit), ]
        
        Peaks_at_no_Origins = Peaks_at_no_Origins[ !duplicated(Peaks_at_no_Origins$peakSummit), ]
        
        # print peak summary
        message("There are", " ",
                nrow(allPeaks), " genome-wide ChIP peaks", "\n",
                nrow(Peaks_at_all_Origins), " were at known origins", "\n",
                nrow(Peaks_at_early_Origins), " were at early origins", "\n",
                nrow(Peaks_at_late_Origins), " were at late origins and", "\n",
                nrow(Peaks_at_no_Origins), " were at non-origin positions.")
        
        # save all five peak classes
        write.table(
          allPeaks[,1:5],
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_Genomewide_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_all_Origins,
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_Origin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_early_Origins,
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_EarlyOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_late_Origins,
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_LateOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_no_Origins[,1:5],
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_NonOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        unlink(c(IP, Input, peakFile), recursive=TRUE, force=TRUE)
        unlink(outDir, recursive=TRUE, force=TRUE)
      }
      
      ChIP_Peaks_Processing(IPBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_ChIP.bam"),
                            InBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input.bam"))
      
    }
    
  }
  
  if(Alignment == "malign"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Peaks_ma"))){
      
      ChIP_Peaks_Processing <- function(IPBam, InBam){
        
        message(paste0("➤ Processing ChIP peaks for ", Pro_1))
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Peaks_ma")))
        
        tempdir(check=TRUE)
        
        IP = tempfile(fileext=".bed")
        Input = tempfile(fileext=".bed")
        peakFile = tempfile(fileext=".bed")
        outDir = tempfile(pattern="MACS2_")
        
        suppressWarnings(dir.create(outDir))
        
        # convert ChIP and Input BAM files to BED
        command_1 = "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s > %s"
        
        system(sprintf(command_1, IPBam, IP))
        system(sprintf(command_1, InBam, Input))
        
        # call genome-wide ChIP peaks
        command_2 = "macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n Peak --outdir %s 2> /dev/null"
        
        system(sprintf(command_2, IP, Input, outDir))
        
        allPeaks = read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#")
        
        write.table(allPeaks, file=peakFile, quote=FALSE, row.names=FALSE,
                    sep="\t", col.names=FALSE)
        
        ColHeads = "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\""
        
        # peaks associated with all known origins
        command_3 = paste0(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools intersect ",
          "-wa -wb -a %s -b %s | ",
          "awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'"
        )
        
        Peaks_at_all_Origins = read.table(
          pipe(sprintf(command_3, peakFile, All_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with early origins
        Peaks_at_early_Origins = read.table(
          pipe(sprintf(command_3, peakFile, E_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with late origins
        Peaks_at_late_Origins = read.table(
          pipe(sprintf(command_3, peakFile, L_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # standardize the MACS2 peak columns
        allPeaks = dplyr::rename(allPeaks,
                                 chrom=chr,
                                 peakStart=start,
                                 peakEnd=end,
                                 peakLength=length,
                                 peakSummit=abs_summit)
        
        # peaks not associated with known origins
        Peaks_at_no_Origins = dplyr::anti_join(
          allPeaks,
          Peaks_at_all_Origins,
          by=c("chrom", "peakStart", "peakEnd")
        )
        
        # remove duplicate peak entries
        Peaks_at_all_Origins = Peaks_at_all_Origins[ !duplicated(Peaks_at_all_Origins$peakSummit), ]
        
        Peaks_at_early_Origins = Peaks_at_early_Origins[ !duplicated(Peaks_at_early_Origins$peakSummit), ]
        
        Peaks_at_late_Origins = Peaks_at_late_Origins[ !duplicated(Peaks_at_late_Origins$peakSummit), ]
        
        Peaks_at_no_Origins = Peaks_at_no_Origins[ !duplicated(Peaks_at_no_Origins$peakSummit), ]
        
        # print peak summary
        message("There are", " ",
                nrow(allPeaks), " genome-wide ChIP peaks", "\n",
                nrow(Peaks_at_all_Origins), " were at known origins", "\n",
                nrow(Peaks_at_early_Origins), " were at early origins", "\n",
                nrow(Peaks_at_late_Origins), " were at late origins and", "\n",
                nrow(Peaks_at_no_Origins), " were at non-origin positions.")
        
        # save all five peak classes
        write.table(
          allPeaks[,1:5],
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_Genomewide_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_all_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_Origin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_early_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_EarlyOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_late_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_LateOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_no_Origins[,1:5],
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_NonOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        unlink(c(IP, Input, peakFile), recursive=TRUE, force=TRUE)
        unlink(outDir, recursive=TRUE, force=TRUE)
      }
      
      ChIP_Peaks_Processing(IPBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_ChIP.bam"),
                            InBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input.bam"))
      
    }
    
  }
  
  if(Alignment == "mrdna"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Peaks_ma"))){
      
      ChIP_Peaks_Processing <- function(IPBam, InBam){
        
        message(paste0("➤ Processing ChIP peaks for ", Pro_1))
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Peaks_ma")))
        
        tempdir(check=TRUE)
        
        IP = tempfile(fileext=".bed")
        Input = tempfile(fileext=".bed")
        peakFile = tempfile(fileext=".bed")
        outDir = tempfile(pattern="MACS2_")
        
        suppressWarnings(dir.create(outDir))
        
        # convert ChIP and Input BAM files to BED
        command_1 = "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s > %s"
        
        system(sprintf(command_1, IPBam, IP))
        system(sprintf(command_1, InBam, Input))
        
        # call genome-wide ChIP peaks
        command_2 = "macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n Peak --outdir %s 2> /dev/null"
        
        system(sprintf(command_2, IP, Input, outDir))
        
        allPeaks = read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#")
        
        write.table(allPeaks, file=peakFile, quote=FALSE, row.names=FALSE,
                    sep="\t", col.names=FALSE)
        
        ColHeads = "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\""
        
        # peaks associated with all known origins
        command_3 = paste0(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools intersect ",
          "-wa -wb -a %s -b %s | ",
          "awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'"
        )
        
        Peaks_at_all_Origins = read.table(
          pipe(sprintf(command_3, peakFile, All_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with early origins
        Peaks_at_early_Origins = read.table(
          pipe(sprintf(command_3, peakFile, E_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with late origins
        Peaks_at_late_Origins = read.table(
          pipe(sprintf(command_3, peakFile, L_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # standardize the MACS2 peak columns
        allPeaks = dplyr::rename(allPeaks,
                                 chrom=chr,
                                 peakStart=start,
                                 peakEnd=end,
                                 peakLength=length,
                                 peakSummit=abs_summit)
        
        # peaks not associated with known origins
        Peaks_at_no_Origins = dplyr::anti_join(
          allPeaks,
          Peaks_at_all_Origins,
          by=c("chrom", "peakStart", "peakEnd")
        )
        
        # remove duplicate peak entries
        Peaks_at_all_Origins = Peaks_at_all_Origins[ !duplicated(Peaks_at_all_Origins$peakSummit), ]
        
        Peaks_at_early_Origins = Peaks_at_early_Origins[ !duplicated(Peaks_at_early_Origins$peakSummit), ]
        
        Peaks_at_late_Origins = Peaks_at_late_Origins[ !duplicated(Peaks_at_late_Origins$peakSummit), ]
        
        Peaks_at_no_Origins = Peaks_at_no_Origins[ !duplicated(Peaks_at_no_Origins$peakSummit), ]
        
        # print peak summary
        message("There are", " ",
                nrow(allPeaks), " genome-wide ChIP peaks", "\n",
                nrow(Peaks_at_all_Origins), " were at known origins", "\n",
                nrow(Peaks_at_early_Origins), " were at early origins", "\n",
                nrow(Peaks_at_late_Origins), " were at late origins and", "\n",
                nrow(Peaks_at_no_Origins), " were at non-origin positions.")
        
        # save all five peak classes
        write.table(
          allPeaks[,1:5],
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_Genomewide_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_all_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_Origin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_early_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_EarlyOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_late_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_LateOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_no_Origins[,1:5],
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_NonOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        unlink(c(IP, Input, peakFile), recursive=TRUE, force=TRUE)
        unlink(outDir, recursive=TRUE, force=TRUE)
      }
      
      ChIP_Peaks_Processing(IPBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_ChIP_ma.bam"),
                            InBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input_ma.bam"))
      
    }
    
  }
  
  
  ## Calculate Ratio
  #
  message("⏳ Calculating enrichment ratios...")
  #
  
  if(Alignment == "generic"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Ratios"))){
      
      CalculateRatio <- function(IP_coverage, Input_coverage, RatioName,
                                 NoiseChunkSizeBp=2000, NoiseIterations=2500,
                                 NoiseSeed=123, NoiseSmoothingSpar=0.65,
                                 NoiseFloor=1e-6){
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Ratios")))
        
        # read and optionally collapse coverage files
        ReadCoverage <- function(coverageFiles){
          
          coverage = read.table(coverageFiles[1], header=FALSE)
          
          if(length(coverageFiles) > 1){
            
            for(j in 2:length(coverageFiles)){
              
              additionalCoverage = read.table(coverageFiles[j], header=FALSE)
              
              coordinatesMatch = identical(as.character(coverage[,1]),
                                           as.character(additionalCoverage[,1])) &&
                identical(as.numeric(coverage[,2]), as.numeric(additionalCoverage[,2])) &&
                identical(as.numeric(coverage[,3]), as.numeric(additionalCoverage[,3]))
              
              if(coordinatesMatch == FALSE){
                stop("Coverage coordinates do not match for strand collapsing.")
              }
              
              coverage[,5] = as.numeric(coverage[,5]) + as.numeric(additionalCoverage[,5])
            }
          }
          
          coverage
        }
        
        IP.df = ReadCoverage(IP_coverage)
        In.df = ReadCoverage(Input_coverage)
        
        coordinatesMatch = identical(as.character(IP.df[,1]), as.character(In.df[,1])) &&
          identical(as.numeric(IP.df[,2]), as.numeric(In.df[,2])) &&
          identical(as.numeric(IP.df[,3]), as.numeric(In.df[,3]))
        
        if(coordinatesMatch == FALSE){
          stop("ChIP and Input coverage coordinates do not match.")
        }
        
        IP.df[,5] = as.numeric(IP.df[,5])
        In.df[,5] = as.numeric(In.df[,5])
        
        # read genome-wide ChIP peaks
        PeakFile = paste0(dir, Pro_1, "/", "Peaks", "/",
                          Pro_1, "_Genomewide_Peaks.bed")
        
        if(!file.exists(PeakFile)){
          stop("Genome-wide peak file is missing: ", PeakFile)
        }
        
        Peaks = read.table(PeakFile, header=TRUE)
        
        Peaks$peakStart = as.numeric(Peaks$peakStart)
        Peaks$peakEnd = as.numeric(Peaks$peakEnd)
        
        # merge overlapping peak intervals
        MergeIntervals <- function(intervals){
          
          if(nrow(intervals) == 0){
            return(intervals)
          }
          
          intervals = intervals[order(intervals$start, intervals$end), ]
          merged = intervals[1, , drop=FALSE]
          
          if(nrow(intervals) > 1){
            
            for(j in 2:nrow(intervals)){
              
              last = nrow(merged)
              
              if(intervals$start[j] <= merged$end[last]){
                
                merged$end[last] = max(merged$end[last], intervals$end[j])
                
              } else {
                
                merged = rbind(merged, intervals[j, , drop=FALSE])
              }
            }
          }
          
          merged
        }
        
        # construct non-peak regions
        BuildNonPeakGaps <- function(chrom, chromLength){
          
          Peaks_chr = Peaks[Peaks$chrom == chrom, , drop=FALSE]
          
          if(nrow(Peaks_chr) == 0){
            
            return(data.frame(start=0, end=chromLength, length=chromLength))
          }
          
          intervals = data.frame(start=pmax(0, Peaks_chr$peakStart),
                                 end=pmin(chromLength, Peaks_chr$peakEnd))
          
          intervals = intervals[
            is.finite(intervals$start) &
              is.finite(intervals$end) &
              intervals$end > intervals$start, ]
          
          if(nrow(intervals) == 0){
            
            return(data.frame(start=0, end=chromLength, length=chromLength))
          }
          
          intervals = MergeIntervals(intervals)
          
          gaps = data.frame(start=numeric(), end=numeric())
          cursor = 0
          
          for(j in 1:nrow(intervals)){
            
            if(intervals$start[j] > cursor){
              
              gaps = rbind(gaps, data.frame(start=cursor,
                                            end=intervals$start[j]))
            }
            
            cursor = max(cursor, intervals$end[j])
          }
          
          if(cursor < chromLength){
            
            gaps = rbind(gaps, data.frame(start=cursor, end=chromLength))
          }
          
          if(nrow(gaps) == 0){
            
            return(data.frame(start=numeric(), end=numeric(),
                              length=numeric()))
          }
          
          gaps$length = gaps$end - gaps$start
          gaps = gaps[gaps$length >= 200, , drop=FALSE]
          
          gaps
        }
        
        SummariseValues <- function(values){
          
          values = values[is.finite(values)]
          
          if(length(values) == 0){
            return(NA_real_)
          }
          
          median(values)
        }
        
        # predict chromosome-position-dependent noise
        PredictNoise <- function(binCenters, sampleCenters, sampleValues,
                                 fallbackValues){
          
          sampleOK = is.finite(sampleCenters) & is.finite(sampleValues)
          
          sampleCenters = sampleCenters[sampleOK]
          sampleValues = sampleValues[sampleOK]
          
          positiveBackground = fallbackValues[
            is.finite(fallbackValues) & fallbackValues > NoiseFloor
          ]
          
          if(length(positiveBackground) > 0){
            
            backgroundFloor = max(
              NoiseFloor,
              as.numeric(quantile(positiveBackground, probs=0.01,
                                  na.rm=TRUE, names=FALSE))
            )
            
          } else {
            
            backgroundFloor = NoiseFloor
          }
          
          fallback = median(fallbackValues, na.rm=TRUE)
          
          if(!is.finite(fallback) || fallback < backgroundFloor){
            fallback = backgroundFloor
          }
          
          if(length(sampleValues) < 2 ||
             length(unique(sampleCenters)) < 2){
            
            return(rep(fallback, length(binCenters)))
          }
          
          # background smoothing is performed on the log scale
          sampleValues = log(pmax(sampleValues, backgroundFloor))
          
          sample.df = data.frame(center=round(sampleCenters),
                                 value=sampleValues)
          
          sample.df = aggregate(value~center, data=sample.df,
                                FUN=function(x) median(x, na.rm=TRUE))
          
          sample.df = sample.df[order(sample.df$center), ]
          
          if(nrow(sample.df) < 4 ||
             length(unique(sample.df$value)) < 2){
            
            prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                rule=2, ties="ordered")$y
            
          } else {
            
            splineObject = tryCatch(
              smooth.spline(sample.df$center, sample.df$value,
                            spar=NoiseSmoothingSpar),
              error=function(e) NULL
            )
            
            if(is.null(splineObject)){
              
              prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                  rule=2, ties="ordered")$y
              
            } else {
              
              prediction = predict(splineObject, binCenters)$y
            }
          }
          
          prediction[!is.finite(prediction)] = log(fallback)
          
          prediction = pmin(
            pmax(prediction, min(sample.df$value)),
            max(sample.df$value)
          )
          
          prediction = exp(prediction)
          prediction[prediction < backgroundFloor] = backgroundFloor
          
          prediction
        }
        
        # estimate ChIP and Input noise separately
        EstimateNoise <- function(){
          
          set.seed(NoiseSeed)
          
          ChIP_noise = numeric(nrow(IP.df))
          Input_noise = numeric(nrow(In.df))
          
          chroms = unique(as.character(IP.df[,1]))
          
          for(chr in chroms){
            
            chrIndex = which(IP.df[,1] == chr)
            
            chrStart = as.numeric(IP.df[chrIndex,2])
            chrEnd = as.numeric(IP.df[chrIndex,3])
            chrCenters = (chrStart + chrEnd)/2
            chrLength = max(chrEnd, na.rm=TRUE)
            
            gaps = BuildNonPeakGaps(chr, chrLength)
            
            if(nrow(gaps) == 0){
              
              ChIP_noise[chrIndex] = max(NoiseFloor,
                                         median(IP.df[chrIndex,5], na.rm=TRUE))
              
              Input_noise[chrIndex] = max(NoiseFloor,
                                          median(In.df[chrIndex,5], na.rm=TRUE))
              
              next
            }
            
            gapProbability = gaps$length/sum(gaps$length)
            
            sampledGapIndex = sample(1:nrow(gaps), size=NoiseIterations,
                                     replace=TRUE, prob=gapProbability)
            
            sampleCenters = numeric(NoiseIterations)
            sampledChIP = numeric(NoiseIterations)
            sampledInput = numeric(NoiseIterations)
            
            for(j in 1:NoiseIterations){
              
              gap = gaps[sampledGapIndex[j], ]
              chunkSize = min(NoiseChunkSizeBp, gap$length)
              
              if(gap$length > chunkSize){
                
                chunkStart = runif(1, min=gap$start,
                                   max=gap$end-chunkSize)
                
              } else {
                
                chunkStart = gap$start
              }
              
              chunkEnd = chunkStart + chunkSize
              
              chunkIndex = which(chrCenters >= chunkStart &
                                   chrCenters <= chunkEnd)
              
              sampleCenters[j] = (chunkStart + chunkEnd)/2
              
              sampledChIP[j] = SummariseValues(
                IP.df[chrIndex[chunkIndex],5]
              )
              
              sampledInput[j] = SummariseValues(
                In.df[chrIndex[chunkIndex],5]
              )
            }
            
            nonPeakIndex = rep(FALSE, length(chrIndex))
            
            for(j in 1:nrow(gaps)){
              
              nonPeakIndex = nonPeakIndex |
                (chrCenters >= gaps$start[j] &
                   chrCenters <= gaps$end[j])
            }
            
            ChIP_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledChIP,
              fallbackValues=IP.df[chrIndex[nonPeakIndex],5]
            )
            
            Input_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledInput,
              fallbackValues=In.df[chrIndex[nonPeakIndex],5]
            )
          }
          
          list(ChIP_noise=ChIP_noise, Input_noise=Input_noise)
        }
        
        noise = EstimateNoise()
        
        # library-size normalization
        IP_Sum = sum(IP.df[,5], na.rm=TRUE)
        In_Sum = sum(In.df[,5], na.rm=TRUE)
        
        if(!is.finite(In_Sum) || In_Sum <= 0){
          stop("Input coverage sum is zero.")
        }
        
        corrFactor = IP_Sum/In_Sum
        
        In.score.norm = In.df[,5]*corrFactor
        In.noise.norm = noise$Input_noise*corrFactor
        
        # calculate the three enrichment ratios
        Ratio.ipin = IP.df[,5]/In.score.norm
        
        Ratio.ipnoise = IP.df[,5]/noise$ChIP_noise
        
        Ratio.ipin.noise = Ratio.ipnoise/
          (In.score.norm/In.noise.norm)
        
        Ratio.ipin[!is.finite(Ratio.ipin)] = 0
        Ratio.ipnoise[!is.finite(Ratio.ipnoise)] = 0
        Ratio.ipin.noise[!is.finite(Ratio.ipin.noise)] = 0
        
        # Poisson enrichment probability
        pvalue = ppois(q=IP.df[,5]-1, lambda=In.score.norm,
                       lower.tail=FALSE, log=FALSE)
        
        ratio.df = data.frame(
          chrom=IP.df[,1],
          chromStart=IP.df[,2],
          chromEnd=IP.df[,3],
          name=RatioName,
          ip.score=IP.df[,5],
          in.score=round(In.score.norm, 4),
          ip.noise=round(noise$ChIP_noise, 4),
          in.noise=round(In.noise.norm, 4),
          ratio.ipin=round(Ratio.ipin, 4),
          ratio.ipnoise=round(Ratio.ipnoise, 4),
          ratio.ipin.noise=round(Ratio.ipin.noise, 4),
          pvalue=pvalue
        )
        
        write.table(
          ratio.df,
          file=paste0(dir, Pro_1, "/", "Ratios", "/", RatioName, ".bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
      }
      
      # Watson-strand ratio
      CalculateRatio(IP_coverage=paste0(dir, Pro_1, "/Coverage/", Pro_1, "_ChIP_watson.bed"),
                     Input_coverage=paste0(dir, Pro_1, "/Coverage/", Pro_1, "_Input_watson.bed"),
                     RatioName=paste0(Pro_1, "_ChIP_watson"))
      
      # Crick-strand ratio
      CalculateRatio(IP_coverage=paste0(dir, Pro_1, "/Coverage/", Pro_1, "_ChIP_crick.bed"),
                     Input_coverage=paste0(dir, Pro_1, "/Coverage/", Pro_1, "_Input_crick.bed"),
                     RatioName=paste0(Pro_1, "_ChIP_crick"))
      
      # strand-collapsed ratio
      CalculateRatio(IP_coverage=c(paste0(dir, Pro_1, "/Coverage/", Pro_1, "_ChIP_watson.bed"),
                                   paste0(dir, Pro_1, "/Coverage/", Pro_1, "_ChIP_crick.bed")),
                     Input_coverage=c(paste0(dir, Pro_1, "/Coverage/", Pro_1, "_Input_watson.bed"),
                                      paste0(dir, Pro_1, "/Coverage/", Pro_1, "_Input_crick.bed")),
                     RatioName=paste0(Pro_1, "_ChIP_collapsed"))
    }
  }
  
  if(Alignment == "malign"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Ratios_ma"))){
      
      CalculateRatio <- function(IP_coverage, Input_coverage, RatioName,
                                 NoiseChunkSizeBp=2000, NoiseIterations=2500,
                                 NoiseSeed=123, NoiseSmoothingSpar=0.65,
                                 NoiseFloor=1e-6){
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Ratios_ma")))
        
        # read and optionally collapse coverage files
        ReadCoverage <- function(coverageFiles){
          
          coverage = read.table(coverageFiles[1], header=FALSE)
          
          if(length(coverageFiles) > 1){
            
            for(j in 2:length(coverageFiles)){
              
              additionalCoverage = read.table(coverageFiles[j], header=FALSE)
              
              coordinatesMatch = identical(as.character(coverage[,1]),
                                           as.character(additionalCoverage[,1])) &&
                identical(as.numeric(coverage[,2]), as.numeric(additionalCoverage[,2])) &&
                identical(as.numeric(coverage[,3]), as.numeric(additionalCoverage[,3]))
              
              if(coordinatesMatch == FALSE){
                stop("Coverage coordinates do not match for strand collapsing.")
              }
              
              coverage[,5] = as.numeric(coverage[,5]) + as.numeric(additionalCoverage[,5])
            }
          }
          
          coverage
        }
        
        IP.df = ReadCoverage(IP_coverage)
        In.df = ReadCoverage(Input_coverage)
        
        coordinatesMatch = identical(as.character(IP.df[,1]), as.character(In.df[,1])) &&
          identical(as.numeric(IP.df[,2]), as.numeric(In.df[,2])) &&
          identical(as.numeric(IP.df[,3]), as.numeric(In.df[,3]))
        
        if(coordinatesMatch == FALSE){
          stop("ChIP and Input coverage coordinates do not match.")
        }
        
        IP.df[,5] = as.numeric(IP.df[,5])
        In.df[,5] = as.numeric(In.df[,5])
        
        # read genome-wide ChIP peaks
        PeakFile = paste0(dir, Pro_1, "/", "Peaks_ma", "/", Pro_1, "_Genomewide_Peaks.bed")
        
        if(!file.exists(PeakFile)){
          stop("Genome-wide peak file is missing: ", PeakFile)
        }
        
        Peaks = read.table(PeakFile, header=TRUE)
        
        Peaks$peakStart = as.numeric(Peaks$peakStart)
        Peaks$peakEnd = as.numeric(Peaks$peakEnd)
        
        # merge overlapping peak intervals
        MergeIntervals <- function(intervals){
          
          if(nrow(intervals) == 0){
            return(intervals)
          }
          
          intervals = intervals[order(intervals$start, intervals$end), ]
          merged = intervals[1, , drop=FALSE]
          
          if(nrow(intervals) > 1){
            
            for(j in 2:nrow(intervals)){
              
              last = nrow(merged)
              
              if(intervals$start[j] <= merged$end[last]){
                
                merged$end[last] = max(merged$end[last], intervals$end[j])
                
              } else {
                
                merged = rbind(merged, intervals[j, , drop=FALSE])
              }
            }
          }
          
          merged
        }
        
        # construct non-peak regions
        BuildNonPeakGaps <- function(chrom, chromLength){
          
          Peaks_chr = Peaks[Peaks$chrom == chrom, , drop=FALSE]
          
          if(nrow(Peaks_chr) == 0){
            
            return(data.frame(start=0, end=chromLength, length=chromLength))
          }
          
          intervals = data.frame(start=pmax(0, Peaks_chr$peakStart),
                                 end=pmin(chromLength, Peaks_chr$peakEnd))
          
          intervals = intervals[
            is.finite(intervals$start) &
              is.finite(intervals$end) &
              intervals$end > intervals$start, ]
          
          if(nrow(intervals) == 0){
            
            return(data.frame(start=0, end=chromLength, length=chromLength))
          }
          
          intervals = MergeIntervals(intervals)
          
          gaps = data.frame(start=numeric(), end=numeric())
          cursor = 0
          
          for(j in 1:nrow(intervals)){
            
            if(intervals$start[j] > cursor){
              
              gaps = rbind(gaps, data.frame(start=cursor,
                                            end=intervals$start[j]))
            }
            
            cursor = max(cursor, intervals$end[j])
          }
          
          if(cursor < chromLength){
            
            gaps = rbind(gaps, data.frame(start=cursor, end=chromLength))
          }
          
          if(nrow(gaps) == 0){
            
            return(data.frame(start=numeric(), end=numeric(),
                              length=numeric()))
          }
          
          gaps$length = gaps$end - gaps$start
          gaps = gaps[gaps$length >= 200, , drop=FALSE]
          
          gaps
        }
        
        SummariseValues <- function(values){
          
          values = values[is.finite(values)]
          
          if(length(values) == 0){
            return(NA_real_)
          }
          
          median(values)
        }
        
        # predict chromosome-position-dependent noise
        PredictNoise <- function(binCenters, sampleCenters, sampleValues,
                                 fallbackValues){
          
          sampleOK = is.finite(sampleCenters) & is.finite(sampleValues)
          
          sampleCenters = sampleCenters[sampleOK]
          sampleValues = sampleValues[sampleOK]
          
          positiveBackground = fallbackValues[
            is.finite(fallbackValues) & fallbackValues > NoiseFloor
          ]
          
          if(length(positiveBackground) > 0){
            
            backgroundFloor = max(
              NoiseFloor,
              as.numeric(quantile(positiveBackground, probs=0.01,
                                  na.rm=TRUE, names=FALSE))
            )
            
          } else {
            
            backgroundFloor = NoiseFloor
          }
          
          fallback = median(fallbackValues, na.rm=TRUE)
          
          if(!is.finite(fallback) || fallback < backgroundFloor){
            fallback = backgroundFloor
          }
          
          if(length(sampleValues) < 2 ||
             length(unique(sampleCenters)) < 2){
            
            return(rep(fallback, length(binCenters)))
          }
          
          # background smoothing is performed on the log scale
          sampleValues = log(pmax(sampleValues, backgroundFloor))
          
          sample.df = data.frame(center=round(sampleCenters),
                                 value=sampleValues)
          
          sample.df = aggregate(value~center, data=sample.df,
                                FUN=function(x) median(x, na.rm=TRUE))
          
          sample.df = sample.df[order(sample.df$center), ]
          
          if(nrow(sample.df) < 4 ||
             length(unique(sample.df$value)) < 2){
            
            prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                rule=2, ties="ordered")$y
            
          } else {
            
            splineObject = tryCatch(
              smooth.spline(sample.df$center, sample.df$value,
                            spar=NoiseSmoothingSpar),
              error=function(e) NULL
            )
            
            if(is.null(splineObject)){
              
              prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                  rule=2, ties="ordered")$y
              
            } else {
              
              prediction = predict(splineObject, binCenters)$y
            }
          }
          
          prediction[!is.finite(prediction)] = log(fallback)
          
          prediction = pmin(
            pmax(prediction, min(sample.df$value)),
            max(sample.df$value)
          )
          
          prediction = exp(prediction)
          prediction[prediction < backgroundFloor] = backgroundFloor
          
          prediction
        }
        
        # estimate ChIP and Input noise separately
        EstimateNoise <- function(){
          
          set.seed(NoiseSeed)
          
          ChIP_noise = numeric(nrow(IP.df))
          Input_noise = numeric(nrow(In.df))
          
          chroms = unique(as.character(IP.df[,1]))
          
          for(chr in chroms){
            
            chrIndex = which(IP.df[,1] == chr)
            
            chrStart = as.numeric(IP.df[chrIndex,2])
            chrEnd = as.numeric(IP.df[chrIndex,3])
            chrCenters = (chrStart + chrEnd)/2
            chrLength = max(chrEnd, na.rm=TRUE)
            
            gaps = BuildNonPeakGaps(chr, chrLength)
            
            if(nrow(gaps) == 0){
              
              ChIP_noise[chrIndex] = max(NoiseFloor,
                                         median(IP.df[chrIndex,5], na.rm=TRUE))
              
              Input_noise[chrIndex] = max(NoiseFloor,
                                          median(In.df[chrIndex,5], na.rm=TRUE))
              
              next
            }
            
            gapProbability = gaps$length/sum(gaps$length)
            
            sampledGapIndex = sample(1:nrow(gaps), size=NoiseIterations,
                                     replace=TRUE, prob=gapProbability)
            
            sampleCenters = numeric(NoiseIterations)
            sampledChIP = numeric(NoiseIterations)
            sampledInput = numeric(NoiseIterations)
            
            for(j in 1:NoiseIterations){
              
              gap = gaps[sampledGapIndex[j], ]
              chunkSize = min(NoiseChunkSizeBp, gap$length)
              
              if(gap$length > chunkSize){
                
                chunkStart = runif(1, min=gap$start,
                                   max=gap$end-chunkSize)
                
              } else {
                
                chunkStart = gap$start
              }
              
              chunkEnd = chunkStart + chunkSize
              
              chunkIndex = which(chrCenters >= chunkStart &
                                   chrCenters <= chunkEnd)
              
              sampleCenters[j] = (chunkStart + chunkEnd)/2
              
              sampledChIP[j] = SummariseValues(
                IP.df[chrIndex[chunkIndex],5]
              )
              
              sampledInput[j] = SummariseValues(
                In.df[chrIndex[chunkIndex],5]
              )
            }
            
            nonPeakIndex = rep(FALSE, length(chrIndex))
            
            for(j in 1:nrow(gaps)){
              
              nonPeakIndex = nonPeakIndex |
                (chrCenters >= gaps$start[j] &
                   chrCenters <= gaps$end[j])
            }
            
            ChIP_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledChIP,
              fallbackValues=IP.df[chrIndex[nonPeakIndex],5]
            )
            
            Input_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledInput,
              fallbackValues=In.df[chrIndex[nonPeakIndex],5]
            )
          }
          
          list(ChIP_noise=ChIP_noise, Input_noise=Input_noise)
        }
        
        noise = EstimateNoise()
        
        # library-size normalization
        IP_Sum = sum(IP.df[,5], na.rm=TRUE)
        In_Sum = sum(In.df[,5], na.rm=TRUE)
        
        if(!is.finite(In_Sum) || In_Sum <= 0){
          stop("Input coverage sum is zero.")
        }
        
        corrFactor = IP_Sum/In_Sum
        
        In.score.norm = In.df[,5]*corrFactor
        In.noise.norm = noise$Input_noise*corrFactor
        
        # calculate the three enrichment ratios
        Ratio.ipin = IP.df[,5]/In.score.norm
        
        Ratio.ipnoise = IP.df[,5]/noise$ChIP_noise
        
        Ratio.ipin.noise = Ratio.ipnoise/(In.score.norm/In.noise.norm)
        
        Ratio.ipin[!is.finite(Ratio.ipin)] = 0
        Ratio.ipnoise[!is.finite(Ratio.ipnoise)] = 0
        Ratio.ipin.noise[!is.finite(Ratio.ipin.noise)] = 0
        
        
        ratio.df = data.frame(
          chrom=IP.df[,1],
          chromStart=IP.df[,2],
          chromEnd=IP.df[,3],
          name=RatioName,
          ip.score=IP.df[,5],
          in.score=round(In.score.norm, 4),
          ip.noise=round(noise$ChIP_noise, 4),
          in.noise=round(In.noise.norm, 4),
          ratio.ipin=round(Ratio.ipin, 4),
          ratio.ipnoise=round(Ratio.ipnoise, 4),
          ratio.ipin.noise=round(Ratio.ipin.noise, 4)
        )
        
        write.table(
          ratio.df,
          file=paste0(dir, Pro_1, "/", "Ratios_ma", "/", RatioName, ".bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
      }
      
      # Watson-strand ratio
      CalculateRatio(IP_coverage=paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_ChIP_watson.bed"),
                     Input_coverage=paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_Input_watson.bed"),
                     RatioName=paste0(Pro_1, "_ChIP_watson"))
      
      # Crick-strand ratio
      CalculateRatio(IP_coverage=paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_ChIP_crick.bed"),
                     Input_coverage=paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_Input_crick.bed"),
                     RatioName=paste0(Pro_1, "_ChIP_crick"))
      
      # strand-collapsed ratio
      CalculateRatio(IP_coverage=c(paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_ChIP_watson.bed"),
                                   paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_ChIP_crick.bed")),
                     Input_coverage=c(paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_Input_watson.bed"),
                                      paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_Input_crick.bed")),
                     RatioName=paste0(Pro_1, "_ChIP_collapsed"))
    }
    
  }
  
  if(Alignment == "mrdna"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Ratios_ma_rdna"))){
      
      CalculateRatio <- function(IP_coverage, Input_coverage,
                                 Noise_IP_coverage, Noise_Input_coverage,
                                 RatioName, NoiseChunkSizeBp=2000,
                                 NoiseIterations=2500, NoiseSeed=123,
                                 NoiseSmoothingSpar=0.65, NoiseFloor=1e-6){
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Ratios_ma_rdna")))
        
        # read and optionally collapse coverage files
        ReadCoverage <- function(coverageFiles){
          
          coverage = read.table(coverageFiles[1], header=FALSE)
          
          if(length(coverageFiles) > 1){
            
            for(j in 2:length(coverageFiles)){
              
              additionalCoverage = read.table(coverageFiles[j], header=FALSE)
              
              coordinatesMatch = identical(as.character(coverage[,1]),
                                           as.character(additionalCoverage[,1])) &&
                identical(as.numeric(coverage[,2]), as.numeric(additionalCoverage[,2])) &&
                identical(as.numeric(coverage[,3]), as.numeric(additionalCoverage[,3]))
              
              if(coordinatesMatch == FALSE){
                stop("Coverage coordinates do not match for strand collapsing.")
              }
              
              coverage[,5] = as.numeric(coverage[,5]) +
                as.numeric(additionalCoverage[,5])
            }
          }
          
          coverage
        }
        
        # rDNA-reference coverage used for final ratio calculation
        IP.df = ReadCoverage(IP_coverage)
        In.df = ReadCoverage(Input_coverage)
        
        coordinatesMatch = identical(as.character(IP.df[,1]),
                                     as.character(In.df[,1])) &&
          identical(as.numeric(IP.df[,2]), as.numeric(In.df[,2])) &&
          identical(as.numeric(IP.df[,3]), as.numeric(In.df[,3]))
        
        if(coordinatesMatch == FALSE){
          stop("rDNA ChIP and Input coverage coordinates do not match.")
        }
        
        IP.df[,5] = as.numeric(IP.df[,5])
        In.df[,5] = as.numeric(In.df[,5])
        
        # whole-genome coverage used for background estimation
        Noise.IP.df = ReadCoverage(Noise_IP_coverage)
        Noise.In.df = ReadCoverage(Noise_Input_coverage)
        
        noiseCoordinatesMatch = identical(as.character(Noise.IP.df[,1]),
                                          as.character(Noise.In.df[,1])) &&
          identical(as.numeric(Noise.IP.df[,2]),
                    as.numeric(Noise.In.df[,2])) &&
          identical(as.numeric(Noise.IP.df[,3]),
                    as.numeric(Noise.In.df[,3]))
        
        if(noiseCoordinatesMatch == FALSE){
          stop("Whole-genome ChIP and Input coverage coordinates do not match.")
        }
        
        Noise.IP.df[,5] = as.numeric(Noise.IP.df[,5])
        Noise.In.df[,5] = as.numeric(Noise.In.df[,5])
        
        # read peaks obtained from whole-genome multiple alignments
        PeakFile = paste0(dir, Pro_1, "/", "Peaks_ma", "/",
                          Pro_1, "_Genomewide_Peaks.bed")
        
        if(!file.exists(PeakFile)){
          stop("Genome-wide peak file is missing: ", PeakFile)
        }
        
        Peaks = read.table(PeakFile, header=TRUE)
        
        Peaks$peakStart = as.numeric(Peaks$peakStart)
        Peaks$peakEnd = as.numeric(Peaks$peakEnd)
        
        # merge overlapping peak intervals
        MergeIntervals <- function(intervals){
          
          if(nrow(intervals) == 0){
            return(intervals)
          }
          
          intervals = intervals[order(intervals$start, intervals$end), ]
          merged = intervals[1, , drop=FALSE]
          
          if(nrow(intervals) > 1){
            
            for(j in 2:nrow(intervals)){
              
              last = nrow(merged)
              
              if(intervals$start[j] <= merged$end[last]){
                
                merged$end[last] = max(merged$end[last],
                                       intervals$end[j])
                
              } else {
                
                merged = rbind(merged, intervals[j, , drop=FALSE])
              }
            }
          }
          
          merged
        }
        
        # construct whole-genome non-peak regions
        BuildNonPeakGaps <- function(chrom, chromLength){
          
          Peaks_chr = Peaks[Peaks$chrom == chrom, , drop=FALSE]
          
          if(nrow(Peaks_chr) == 0){
            
            return(data.frame(start=0, end=chromLength,
                              length=chromLength))
          }
          
          intervals = data.frame(start=pmax(0, Peaks_chr$peakStart),
                                 end=pmin(chromLength, Peaks_chr$peakEnd))
          
          intervals = intervals[
            is.finite(intervals$start) &
              is.finite(intervals$end) &
              intervals$end > intervals$start, ]
          
          if(nrow(intervals) == 0){
            
            return(data.frame(start=0, end=chromLength,
                              length=chromLength))
          }
          
          intervals = MergeIntervals(intervals)
          
          gaps = data.frame(start=numeric(), end=numeric())
          cursor = 0
          
          for(j in 1:nrow(intervals)){
            
            if(intervals$start[j] > cursor){
              
              gaps = rbind(
                gaps,
                data.frame(start=cursor, end=intervals$start[j])
              )
            }
            
            cursor = max(cursor, intervals$end[j])
          }
          
          if(cursor < chromLength){
            
            gaps = rbind(
              gaps,
              data.frame(start=cursor, end=chromLength)
            )
          }
          
          if(nrow(gaps) == 0){
            
            return(data.frame(start=numeric(), end=numeric(),
                              length=numeric()))
          }
          
          gaps$length = gaps$end - gaps$start
          gaps = gaps[gaps$length >= 200, , drop=FALSE]
          
          gaps
        }
        
        SummariseValues <- function(values){
          
          values = values[is.finite(values)]
          
          if(length(values) == 0){
            return(NA_real_)
          }
          
          median(values)
        }
        
        # predict chromosome-position-dependent background
        PredictNoise <- function(binCenters, sampleCenters, sampleValues,
                                 fallbackValues){
          
          sampleOK = is.finite(sampleCenters) &
            is.finite(sampleValues)
          
          sampleCenters = sampleCenters[sampleOK]
          sampleValues = sampleValues[sampleOK]
          
          positiveBackground = fallbackValues[
            is.finite(fallbackValues) &
              fallbackValues > NoiseFloor
          ]
          
          if(length(positiveBackground) > 0){
            
            backgroundFloor = max(
              NoiseFloor,
              as.numeric(
                quantile(
                  positiveBackground,
                  probs=0.01,
                  na.rm=TRUE,
                  names=FALSE
                )
              )
            )
            
          } else {
            
            backgroundFloor = NoiseFloor
          }
          
          fallback = median(fallbackValues, na.rm=TRUE)
          
          if(!is.finite(fallback) || fallback < backgroundFloor){
            fallback = backgroundFloor
          }
          
          if(length(sampleValues) < 2 ||
             length(unique(sampleCenters)) < 2){
            
            return(rep(fallback, length(binCenters)))
          }
          
          # background smoothing is performed on the log scale
          sampleValues = log(pmax(sampleValues, backgroundFloor))
          
          sample.df = data.frame(center=round(sampleCenters),
                                 value=sampleValues)
          
          sample.df = aggregate(
            value~center,
            data=sample.df,
            FUN=function(x) median(x, na.rm=TRUE)
          )
          
          sample.df = sample.df[order(sample.df$center), ]
          
          if(nrow(sample.df) < 4 ||
             length(unique(sample.df$value)) < 2){
            
            prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                rule=2, ties="ordered")$y
            
          } else {
            
            splineObject = tryCatch(
              smooth.spline(
                sample.df$center,
                sample.df$value,
                spar=NoiseSmoothingSpar
              ),
              error=function(e) NULL
            )
            
            if(is.null(splineObject)){
              
              prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                  rule=2, ties="ordered")$y
              
            } else {
              
              prediction = predict(splineObject, binCenters)$y
            }
          }
          
          prediction[!is.finite(prediction)] = log(fallback)
          
          prediction = pmin(
            pmax(prediction, min(sample.df$value)),
            max(sample.df$value)
          )
          
          prediction = exp(prediction)
          prediction[prediction < backgroundFloor] = backgroundFloor
          
          prediction
        }
        
        # estimate background from whole-genome non-peak coverage
        EstimateNoise <- function(){
          
          set.seed(NoiseSeed)
          
          ChIP_noise = numeric(nrow(Noise.IP.df))
          Input_noise = numeric(nrow(Noise.In.df))
          
          chroms = unique(as.character(Noise.IP.df[,1]))
          
          for(chr in chroms){
            
            chrIndex = which(Noise.IP.df[,1] == chr)
            
            chrStart = as.numeric(Noise.IP.df[chrIndex,2])
            chrEnd = as.numeric(Noise.IP.df[chrIndex,3])
            chrCenters = (chrStart + chrEnd)/2
            chrLength = max(chrEnd, na.rm=TRUE)
            
            gaps = BuildNonPeakGaps(chr, chrLength)
            
            if(nrow(gaps) == 0){
              
              ChIP_noise[chrIndex] = max(
                NoiseFloor,
                median(Noise.IP.df[chrIndex,5], na.rm=TRUE)
              )
              
              Input_noise[chrIndex] = max(
                NoiseFloor,
                median(Noise.In.df[chrIndex,5], na.rm=TRUE)
              )
              
              next
            }
            
            gapProbability = gaps$length/sum(gaps$length)
            
            sampledGapIndex = sample(
              1:nrow(gaps),
              size=NoiseIterations,
              replace=TRUE,
              prob=gapProbability
            )
            
            sampleCenters = numeric(NoiseIterations)
            sampledChIP = numeric(NoiseIterations)
            sampledInput = numeric(NoiseIterations)
            
            for(j in 1:NoiseIterations){
              
              gap = gaps[sampledGapIndex[j], ]
              chunkSize = min(NoiseChunkSizeBp, gap$length)
              
              if(gap$length > chunkSize){
                
                chunkStart = runif(
                  1,
                  min=gap$start,
                  max=gap$end-chunkSize
                )
                
              } else {
                
                chunkStart = gap$start
              }
              
              chunkEnd = chunkStart + chunkSize
              
              chunkIndex = which(
                chrCenters >= chunkStart &
                  chrCenters <= chunkEnd
              )
              
              sampleCenters[j] = (chunkStart + chunkEnd)/2
              
              sampledChIP[j] = SummariseValues(
                Noise.IP.df[chrIndex[chunkIndex],5]
              )
              
              sampledInput[j] = SummariseValues(
                Noise.In.df[chrIndex[chunkIndex],5]
              )
            }
            
            nonPeakIndex = rep(FALSE, length(chrIndex))
            
            for(j in 1:nrow(gaps)){
              
              nonPeakIndex = nonPeakIndex |
                (chrCenters >= gaps$start[j] &
                   chrCenters <= gaps$end[j])
            }
            
            ChIP_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledChIP,
              fallbackValues=Noise.IP.df[
                chrIndex[nonPeakIndex],5
              ]
            )
            
            Input_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledInput,
              fallbackValues=Noise.In.df[
                chrIndex[nonPeakIndex],5
              ]
            )
          }
          
          list(ChIP_noise=ChIP_noise,
               Input_noise=Input_noise)
        }
        
        noise = EstimateNoise()
        
        # whole-genome library-size normalization
        IP_Sum = sum(Noise.IP.df[,5], na.rm=TRUE)
        In_Sum = sum(Noise.In.df[,5], na.rm=TRUE)
        
        if(!is.finite(In_Sum) || In_Sum <= 0){
          stop("Whole-genome Input coverage sum is zero.")
        }
        
        corrFactor = IP_Sum/In_Sum
        
        # convert whole-genome background profiles into transferable levels
        positive.IP.noise = noise$ChIP_noise[
          is.finite(noise$ChIP_noise) & noise$ChIP_noise > NoiseFloor
        ]
        
        positive.In.noise = noise$Input_noise[
          is.finite(noise$Input_noise) & noise$Input_noise > NoiseFloor
        ]
        
        IP.noise = median(positive.IP.noise, na.rm=TRUE)
        In.noise = median(positive.In.noise, na.rm=TRUE)*corrFactor
        
        if(!is.finite(IP.noise) || IP.noise < NoiseFloor){
          IP.noise = NoiseFloor
        }
        
        if(!is.finite(In.noise) || In.noise < NoiseFloor){
          In.noise = NoiseFloor
        }
        
        # normalize rDNA Input using the whole-genome correction factor
        In.score.norm = In.df[,5]*corrFactor
        
        # calculate rDNA enrichment ratios
        Ratio.ipin = IP.df[,5]/In.score.norm
        
        Ratio.ipnoise = IP.df[,5]/IP.noise
        
        Ratio.ipin.noise = Ratio.ipnoise/
          (In.score.norm/In.noise)
        
        Ratio.ipin[!is.finite(Ratio.ipin)] = 0
        Ratio.ipnoise[!is.finite(Ratio.ipnoise)] = 0
        Ratio.ipin.noise[!is.finite(Ratio.ipin.noise)] = 0
        
        # output contains rDNA coordinates only
        ratio.df = data.frame(
          chrom=IP.df[,1],
          chromStart=IP.df[,2],
          chromEnd=IP.df[,3],
          name=RatioName,
          ip.score=IP.df[,5],
          in.score=round(In.score.norm, 4),
          ip.noise=round(rep(IP.noise, nrow(IP.df)), 4),
          in.noise=round(rep(In.noise, nrow(IP.df)), 4),
          ratio.ipin=round(Ratio.ipin, 4),
          ratio.ipnoise=round(Ratio.ipnoise, 4),
          ratio.ipin.noise=round(Ratio.ipin.noise, 4)
        )
        
        write.table(
          ratio.df,
          file=paste0(dir, Pro_1, "/", "Ratios_ma_rdna", "/",
                      RatioName, ".bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
      }
      
      # Watson-strand rDNA ratio
      CalculateRatio(
        IP_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_ChIP_rDNA_watson.bed"
        ),
        Input_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_Input_rDNA_watson.bed"
        ),
        Noise_IP_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_ChIP_ma_watson.bed"
        ),
        Noise_Input_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_Input_ma_watson.bed"
        ),
        RatioName=paste0(Pro_1, "_ChIP_rDNA_watson")
      )
      
      # Crick-strand rDNA ratio
      CalculateRatio(
        IP_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_ChIP_rDNA_crick.bed"
        ),
        Input_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_Input_rDNA_crick.bed"
        ),
        Noise_IP_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_ChIP_ma_crick.bed"
        ),
        Noise_Input_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_Input_ma_crick.bed"
        ),
        RatioName=paste0(Pro_1, "_ChIP_rDNA_crick")
      )
      
      # strand-collapsed rDNA ratio
      CalculateRatio(
        IP_coverage=c(
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_ChIP_rDNA_watson.bed"
          ),
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_ChIP_rDNA_crick.bed"
          )
        ),
        Input_coverage=c(
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_Input_rDNA_watson.bed"
          ),
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_Input_rDNA_crick.bed"
          )
        ),
        Noise_IP_coverage=c(
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_ChIP_ma_watson.bed"
          ),
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_ChIP_ma_crick.bed"
          )
        ),
        Noise_Input_coverage=c(
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_Input_ma_watson.bed"
          ),
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_Input_ma_crick.bed"
          )
        ),
        RatioName=paste0(Pro_1, "_ChIP_rDNA_collapsed")
      )
    }
    
  }
  
  ##
  
  rm(list=ls())
  gc()
  
  #
  message("✅ Alignment & Primary Analysis complete!")
  #
  
}


BrDUseq_Primary_Analysis <- function(  Input_R1 = "/full/path/to/file_R1.fastq.gz", 
                                       Input_R2 = "/full/path/to/file_R2.fastq.gz", 
                                       
                                       BrDU_R1 = "/full/path/to/file_R1.fastq.gz", 
                                       BrDU_R2 = "/full/path/to/file_R2.fastq.gz",
                                       
                                       Alignment = "generic",  # "generic" (bowtie2) / "malign" (Rsubread multi-aligner) / "mrdna" (Rsubread multi-aligner at rDNA locus) 
                                       
                                       ExpTitle = "Smc5-trial-BrDU", 
                                       Directory = "None", 
                                       slidingWindow = "YES" ) { 
  
  ## load packages
  
  packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
                "VennDiagram", "grid", "gridBase", "gridExtra", "ShortRead", "csaw", "shiny",
                "BSgenome.Scerevisiae.UCSC.sacCer3", "Rsubread", "GenomicAlignments",
                "IRanges", "readxl", "data.table", "ORFik")
  
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  All_Ori_Link <- "/Applications/ngsAnalyser.app/Contents/Resources/app/OriginList_Full.bed"
  E_Ori_Link <- "/Applications/ngsAnalyser.app/Contents/Resources/app/E_Rep.bed"
  L_Ori_Link <- "/Applications/ngsAnalyser.app/Contents/Resources/app/L_Rep.bed"
  
  # Sequencing Alignment & Binned Coverage Calculation
  
  #
  useDef <- function(a,d) ifelse(isTruthy(a), a,d)
  
  ExpTitle = useDef(ExpTitle, "None")
  
  
  if(ExpTitle == "None"){
    Pro_1 <- unlist(strsplit(basename(Input_R1), split='_', fixed=TRUE))[[1]]
  } else {
    Pro_1 <- ExpTitle
  }
  
  #
  message(paste0("✅ Experiment: ", Pro_1))
  #
  
  #
  Directory = useDef(Directory, "None")
  
  if(Directory == "None"){
    dir <- "~/Desktop/"
  } else {
    dir <- paste0(Directory, "/")
  }
  #
  
  suppressWarnings(dir.create(paste0(dir, Pro_1)))  
  
  ## Quality check of fastqs'
  # #
  # message("⏳ Running QC ...")
  # #
  # if(!file.exists(paste0(dir, Pro_1, "/", Pro_1, "_", "QR", ".html"))){
  #   
  #   fls = c(Input_R1, Input_R2, BrDU_R1, BrDU_R2)
  #   
  #   names(fls) = sub(".fastq", "", basename(fls))
  #   
  #   qas = lapply(seq_along(fls),
  #                function(i, fls) qa(readFastq(fls[i]), names(fls)[i]),
  #                fls)
  #   qa = do.call(rbind, qas)
  #   rpt = report(qa, dest = paste0(dir, Pro_1, "/", Pro_1, "_", "QR", ".html"))
  #   
  # }
  # 
  ## Run alignment
  #
  message("✅ Reference yeast genome : S288C")
  message("⏳ Running alignments...")
  #
  if(Alignment == "generic"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Bam"))){
      
      RunAlignment_bowtie2 <- function(File_R1, File_R2, SampName){
        
        message(paste0("➤ Running alignment for ", Pro_1, "_", SampName))
        
        tempdir(check = TRUE)
        
        Sam <- tempfile(fileext = ".sam")
        Bam <- tempfile(fileext = ".bam")
        nmCollate <- tempfile(fileext = ".bam")
        fixMat <- tempfile(fileext = ".bam")
        SrtBam <- tempfile(fileext = ".bam")
        
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Bam")))
        
        AlnLog <- paste0(dir, Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".log")
        SFBam <- paste0(dir, Pro_1, "/", "Bam", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome for the alignment of sequenced data
        ref_index <- "~/Desktop/chip_support/genome_index_ref/S288C_Ref"
        
        #following commands will run the alignemnt, check quality, sort, filter and index the resultant bam file 
        
        system(sprintf("(/Applications/ngsAnalyser.app/Contents/Resources/app/bowtie2-2.4.4-macos-x86_64/bowtie2 -p 8  --no-discordant --fr -x %s -1 %s -2 %s -S %s) 2> %s", 
                       ref_index, File_R1, File_R2, Sam, AlnLog))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -bS -@ 15 -q 30 -f 2 %s > %s", Sam, Bam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools collate -@ 15 -o %s %s", nmCollate, Bam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools fixmate -@ 15 -m %s %s", nmCollate, fixMat))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools sort -l 9 -@ 15 -m 1024M  -O bam -o %s %s", SrtBam, fixMat))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools markdup -@ 15 %s %s", SrtBam, SFBam))
        
        system(sprintf("/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools index -@ 15 %s", SFBam))
        
        unlink(c(Sam, Bam, nmCollate, fixMat, SrtBam), recursive = T, force = T)
        
      }
      
      RunAlignment_bowtie2(Input_R1, Input_R2, "Input")
      RunAlignment_bowtie2(BrDU_R1, BrDU_R2, "BrDU")
      
    }
    
  }
  
  if(Alignment == "malign"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Bam_ma"))){
      
      RunAlignment_subread_malign <- function(File_R1, File_R2, SampName){
        
        message(paste0("➤ Running Rsubread-based multiple-alignment for ", Pro_1, "_", SampName))
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Bam_ma")))
        
        resultBam <- paste0(dir, Pro_1, "/", "Bam_ma", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome
        
        ref_index <- "~/Desktop/chip_support/genome_index_ref/subread_ref"
        
        #run alignment
        
        align(index = ref_index,
              readfile1 = File_R1,
              readfile2 = File_R2,
              output_format = "BAM",
              output_file = resultBam,
              color2base = F,
              type = "dna",
              unique = FALSE,
              nBestLocations = 16,
              nthreads = 8,
              sortReadsByCoordinates = TRUE)
        
      }
      
      RunAlignment_subread_malign(Input_R1, Input_R2, "Input")
      RunAlignment_subread_malign(BrDU_R1, BrDU_R2, "BrDU")
      
    }
    
  }
  
  if(Alignment == "mrdna"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Bam_ma_rdna"))){
      
      RunAlignment_subread_malign_rDNA <- function(File_R1, File_R2, SampName){
        
        message(paste0("➤ Running Rsubread-based multiple-alignment at rDNA for ", Pro_1, "_", SampName))
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Bam_ma_rdna")))
        
        resultBam <- paste0(dir, Pro_1, "/", "Bam_ma_rdna", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome
        
        ref_index <- "~/Desktop/chip_support/genome_index_ref/rDNA_2nts1"
        
        #run alignment
        
        align(index = ref_index,
              readfile1 = File_R1,
              readfile2 = File_R2,
              output_format = "BAM",
              output_file = resultBam,
              color2base = F,
              type = "dna",
              unique = FALSE,
              nBestLocations = 16,
              nthreads = 8,
              sortReadsByCoordinates = TRUE)
        
      }
      
      RunAlignment_subread_malign_rDNA(Input_R1, Input_R2, "Input_rDNA")
      RunAlignment_subread_malign_rDNA(BrDU_R1, BrDU_R2, "BrDU_rDNA")
      
      RunAlignment_subread_malign <- function(File_R1, File_R2, SampName){
        
        message(paste0("➤ Running Rsubread-based multiple-alignment for ", Pro_1, "_", SampName))
        
        Pro_1 <- Pro_1
        Pro_2 <- SampName
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Bam_ma_rdna")))
        
        resultBam <- paste0(dir, Pro_1, "/", "Bam_ma_rdna", "/", Pro_1, "_", Pro_2, ".bam")
        
        #read the indexed reference genome
        
        ref_index <- "~/Desktop/chip_support/genome_index_ref/subread_ref"
        
        #run alignment
        
        align(index = ref_index,
              readfile1 = File_R1,
              readfile2 = File_R2,
              output_format = "BAM",
              output_file = resultBam,
              color2base = F,
              type = "dna",
              unique = FALSE,
              nBestLocations = 16,
              nthreads = 8,
              sortReadsByCoordinates = TRUE)
        
      }
      
      RunAlignment_subread_malign(Input_R1, Input_R2, "Input_ma")
      RunAlignment_subread_malign(BrDU_R1, BrDU_R2, "BrDU_ma")
      
    }
    
  }
  
  #
  ## Calculate genome-wide binned coverage
  #
  message("⏳ Calculating read coverage ...")
  #
  if(Alignment == "generic"){ BamFolder = "Bam" }
  if(Alignment == "malign"){ BamFolder = "Bam_ma" }
  if(Alignment == "mrdna"){ BamFolder = "Bam_ma_rdna" }
  #
  
  if(Alignment == "generic"){
    
    bamFiles <- c(paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_BrDU.bam"),
                  paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input.bam"))
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Coverage"))){
      
      #coverage
      BamCoverage <- function(bamFile, binSize = 300, stepSize = 10, slidingWindow = "YES", byReads_5p = T){
        
        SampleName = tools::file_path_sans_ext(basename(bamFile))
        Pro_2 = substring(SampleName, nchar(Pro_1)+2)
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Coverage"))) 
        
        tempdir(check = TRUE)
        
        GenomFile <- tempfile(fileext = ".txt")
        binFile <- tempfile(fileext = ".bed")
        
        command_1 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools idxstats %s | awk 'BEGIN {OFS=\"\\t\"} {if ($2>0) print ($1,$2)}' >  %s"
        system(sprintf(command_1, bamFile, GenomFile))
        
        if(slidingWindow=="YES"){
          command_2 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools makewindows -g %s -w %s -s %s | sort -k1,1V -k2,2n > %s"
          system(sprintf(command_2, GenomFile, binSize, stepSize, binFile))
        } else {
          command_2 <- "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools makewindows -g %s -w %s | sort -k1,1V -k2,2n > %s"
          system(sprintf(command_2, GenomFile, binSize, binFile))
        }
        
        pncFiles_watson <- tempfile(fileext = ".bed")
        pncFiles_crick <- tempfile(fileext = ".bed")
        
        if(byReads_5p == TRUE){
          
          message(paste0("➤ Calculating coverage with 5' ends of first-mate reads for", " ", 
                         tools::file_path_sans_ext(basename(bamFile)) ))
          
          # calculate coverage at watson strand by 5' end of the first mate reads
          command_3 <- paste0(
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s |",
            "grep -v XS:i: |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -@ 8 -b - |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools genomecov -ibam stdin -strand + -d -5 |",
            "awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' | sort -k1,1V -k2,2n > %s"
          )
          
          # calculate coverage at crick strand by 5' end of the first mate reads
          command_4 <- paste0(
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -h -@ 8 -q 30 -F 3840 -f 64 -L %s %s |",
            "grep -v XS:i: |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/samtools-1.13/samtools view -@ 8 -b - |",
            "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools genomecov -ibam stdin -strand - -d -5 |",
            "awk 'BEGIN {OFS=\"\\t\"} {if ($3>0) print $1,$2,$2,\"%s\",$3}' | sort -k1,1V -k2,2n > %s"
          )
        }
        
        #
        # Watson (+)
        system(sprintf(command_3, binFile, bamFile,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"),
                       pncFiles_watson))
        
        # Crick (-)
        system(sprintf(command_4, binFile, bamFile,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"),
                       pncFiles_crick))
        
        #
        # sum the counts per bin for watson and crick separately and store in finFiles
        
        finFiles_watson <- paste0(dir, Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "watson.bed")
        finFiles_crick <- paste0(dir, Pro_1, "/", "Coverage", "/", Pro_1, "_", Pro_2, "_", "crick.bed")
        
        command_5 <- paste(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools map",
          "-a %s -b %s -null 0 -o sum |",
          "awk 'BEGIN {OFS=\"\\t\"} {if ($4>=0) print $1,$2,$3,\"%s\",$4}' > %s"
        )
        
        # Watson (+)
        system(sprintf(command_5, binFile, pncFiles_watson,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"),
                       finFiles_watson))
        # Crick (-)
        system(sprintf(command_5, binFile, pncFiles_crick,
                       paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"),
                       finFiles_crick))
        
        #
        
        unlink(c(GenomFile, binFile, pncFiles_watson, pncFiles_crick), recursive = T, force = T)
        
      }
      
      for(i in 1:length(bamFiles)){
        BamCoverage(bamFile = bamFiles[i])
      }
      
    }
    
  }
  
  if(Alignment == "malign"){
    
    bamFiles <- c(paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_BrDU.bam"),
                  paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input.bam"))
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Coverage_ma"))){
      
      BamCoverage <- function(bamFile, binSize=300, stepSize=10, byReads_5p=TRUE){
        
        SampleName = tools::file_path_sans_ext(basename(bamFile))
        Pro_2 = substring(SampleName, nchar(Pro_1)+2)
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Coverage_ma")))
        
        # whole-genome sliding windows
        Windows = tileGenome(seqlengths(Scerevisiae)[-17], tilewidth=stepSize)
        Windows = unlist(Windows)
        Windows = suppressWarnings(trim(IRanges::resize(Windows, width=binSize)))
        
        Counts_watson = numeric(length(Windows))
        Counts_crick = numeric(length(Windows))
        
        if(byReads_5p == TRUE){
          
          message(paste0("➤ Calculating multiple-alignment coverage using 5' ends of first-mate reads for ",
                         tools::file_path_sans_ext(basename(bamFile))))
          
        } else {
          
          message(paste0("➤ Calculating multiple-alignment coverage using fragment midpoints for ",
                         tools::file_path_sans_ext(basename(bamFile))))
        }
        
        for(i in 1:16){
          
          paramS = ScanBamParam(flag=scanBamFlag(isProperPair=TRUE, isUnmappedQuery=FALSE,
                                                 hasUnmappedMate=FALSE),
                                what=c("qname", "mapq", "isize", "seq"),
                                tag=c("NH", "NM"),
                                tagFilter=list(NH=i, NM=c(0:2)))
          
          AlignmentPairs = suppressWarnings(readGAlignmentPairs(bamFile, param=paramS, strandMode=1))
          
          if(length(AlignmentPairs) == 0) next
          
          if(byReads_5p == TRUE){
            
            # use the strand-specific 5' end of the first mate
            CoveragePoints = granges(GenomicAlignments::first(AlignmentPairs))
            
            fivePrimePosition = ifelse(as.character(strand(CoveragePoints)) == "+",
                                       start(CoveragePoints), end(CoveragePoints))
            
            start(CoveragePoints) = fivePrimePosition
            end(CoveragePoints) = fivePrimePosition
            
          } else {
            
            # use the midpoint of the complete paired-end fragment
            CoveragePoints = granges(AlignmentPairs)
            
            midpoint = round((start(CoveragePoints) + end(CoveragePoints))/2)
            
            start(CoveragePoints) = midpoint
            end(CoveragePoints) = midpoint
          }
          
          # fractional contribution of multiple alignments
          mcols(CoveragePoints)$weight = 1/i
          
          CoveragePoints_watson = CoveragePoints[as.character(strand(CoveragePoints)) == "+"]
          CoveragePoints_crick = CoveragePoints[as.character(strand(CoveragePoints)) == "-"]
          
          if(length(CoveragePoints_watson) > 0){
            
            Counts_watson = Counts_watson + countOverlapsW(Windows, CoveragePoints_watson,
                                                           weight="weight")
          }
          
          if(length(CoveragePoints_crick) > 0){
            
            Counts_crick = Counts_crick + countOverlapsW(Windows, CoveragePoints_crick,
                                                         weight="weight")
          }
        }
        
        finFile_watson = paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_", Pro_2, "_watson.bed")
        finFile_crick = paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_", Pro_2, "_crick.bed")
        
        bed_watson = data.frame(seqnames=as.character(seqnames(Windows)),
                                start=start(Windows)-1L,
                                end=end(Windows),
                                name=paste0(tools::file_path_sans_ext(basename(bamFile)), "_watson"),
                                score=Counts_watson)
        
        bed_crick = data.frame(seqnames=as.character(seqnames(Windows)),
                               start=start(Windows)-1L,
                               end=end(Windows),
                               name=paste0(tools::file_path_sans_ext(basename(bamFile)), "_crick"),
                               score=Counts_crick)
        
        write.table(bed_watson, file=finFile_watson, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)
        
        write.table(bed_crick, file=finFile_crick, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)
      }
      
      for(i in 1:length(bamFiles)){
        BamCoverage(bamFile = bamFiles[i])
      }
      
    }
    
  }
  
  if(Alignment == "mrdna"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Coverage_ma_rdna"))){
      
      BamCoverage <- function(bamFile, binSize=300, stepSize=10, byReads_5p=TRUE, rDNA=FALSE){
        
        SampleName = tools::file_path_sans_ext(basename(bamFile))
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Coverage_ma_rdna")))
        
        if(rDNA == FALSE){
          
          # whole-genome sliding windows
          Windows = tileGenome(seqlengths(Scerevisiae)[-17], tilewidth=stepSize)
          Windows = unlist(Windows)
          Windows = suppressWarnings(trim(IRanges::resize(Windows, width=binSize)))
        }
        
        if(rDNA == TRUE){
          
          # rDNA-reference sliding windows
          gen = Seqinfo(seqnames="Scer_2xrDNA_unit", seqlengths=18274,
                        isCircular=TRUE, genome="Scer_2xrDNA_unit")
          
          Windows = tileGenome(gen, tilewidth=stepSize)
          Windows = unlist(Windows)
          Windows = suppressWarnings(trim(IRanges::resize(Windows, width=binSize)))
        }
        
        Counts_watson = numeric(length(Windows))
        Counts_crick = numeric(length(Windows))
        
        if(byReads_5p == TRUE){
          
          message(paste0("➤ Calculating multiple-alignment coverage using 5' ends of first-mate reads for ",
                         SampleName))
          
        } else {
          
          message(paste0("➤ Calculating multiple-alignment coverage using fragment midpoints for ",
                         SampleName))
        }
        
        for(i in 1:16){
          
          paramS = ScanBamParam(flag=scanBamFlag(isProperPair=TRUE, isUnmappedQuery=FALSE,
                                                 hasUnmappedMate=FALSE),
                                what=c("qname", "mapq", "isize", "seq"),
                                tag=c("NH", "NM"),
                                tagFilter=list(NH=i, NM=c(0:2)))
          
          AlignmentPairs = suppressWarnings(readGAlignmentPairs(bamFile, param=paramS, strandMode=1))
          
          if(length(AlignmentPairs) == 0) next
          
          if(byReads_5p == TRUE){
            
            # use the strand-specific 5' end of the first mate
            CoveragePoints = granges(GenomicAlignments::first(AlignmentPairs))
            
            fivePrimePosition = ifelse(as.character(strand(CoveragePoints)) == "+",
                                       start(CoveragePoints), end(CoveragePoints))
            
            start(CoveragePoints) = fivePrimePosition
            end(CoveragePoints) = fivePrimePosition
            
          } else {
            
            # use the midpoint of the complete paired-end fragment
            CoveragePoints = granges(AlignmentPairs)
            
            midpoint = round((start(CoveragePoints) + end(CoveragePoints))/2)
            
            start(CoveragePoints) = midpoint
            end(CoveragePoints) = midpoint
          }
          
          # fractional contribution of multiple alignments
          mcols(CoveragePoints)$weight = 1/i
          
          CoveragePoints_watson = CoveragePoints[as.character(strand(CoveragePoints)) == "+"]
          CoveragePoints_crick = CoveragePoints[as.character(strand(CoveragePoints)) == "-"]
          
          if(length(CoveragePoints_watson) > 0){
            
            Counts_watson = Counts_watson + countOverlapsW(Windows, CoveragePoints_watson,
                                                           weight="weight")
          }
          
          if(length(CoveragePoints_crick) > 0){
            
            Counts_crick = Counts_crick + countOverlapsW(Windows, CoveragePoints_crick,
                                                         weight="weight")
          }
        }
        
        finFile_watson = paste0(dir, Pro_1, "/Coverage_ma_rdna/", SampleName, "_watson.bed")
        finFile_crick = paste0(dir, Pro_1, "/Coverage_ma_rdna/", SampleName, "_crick.bed")
        
        bed_watson = data.frame(seqnames=as.character(seqnames(Windows)),
                                start=start(Windows)-1L,
                                end=end(Windows),
                                name=paste0(SampleName, "_watson"),
                                score=Counts_watson)
        
        bed_crick = data.frame(seqnames=as.character(seqnames(Windows)),
                               start=start(Windows)-1L,
                               end=end(Windows),
                               name=paste0(SampleName, "_crick"),
                               score=Counts_crick)
        
        write.table(bed_watson, file=finFile_watson, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)
        
        write.table(bed_crick, file=finFile_crick, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)
      }
      
      BamCoverage(bamFile=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_BrDU_rDNA.bam"), rDNA=TRUE)
      BamCoverage(bamFile=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input_rDNA.bam"), rDNA=TRUE)
      BamCoverage(bamFile=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_BrDU_ma.bam"), rDNA=FALSE)
      BamCoverage(bamFile=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input_ma.bam"), rDNA=FALSE)
      
    }
    
  }
  
  
  ## Define and process peaks
  #
  message("⏳ Processing BrDU peaks ...")
  #
  
  if(Alignment == "generic"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Peaks"))){
      
      BrDU_Peaks_Processing <- function(IPBam, InBam){
        
        message(paste0("➤ Processing BrDU peaks for ", Pro_1))
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Peaks")))
        
        tempdir(check=TRUE)
        
        IP = tempfile(fileext=".bed")
        Input = tempfile(fileext=".bed")
        peakFile = tempfile(fileext=".bed")
        outDir = tempfile(pattern="MACS2_")
        
        suppressWarnings(dir.create(outDir))
        
        # convert BrDU and Input BAM files to BED
        command_1 = "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s > %s"
        
        system(sprintf(command_1, IPBam, IP))
        system(sprintf(command_1, InBam, Input))
        
        # call genome-wide BrDU peaks
        command_2 = "macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n Peak --outdir %s 2> /dev/null"
        
        system(sprintf(command_2, IP, Input, outDir))
        
        allPeaks = read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#")
        
        write.table(allPeaks, file=peakFile, quote=FALSE, row.names=FALSE,
                    sep="\t", col.names=FALSE)
        
        ColHeads = "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\""
        
        # peaks associated with all known origins
        command_3 = paste0(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools intersect ",
          "-wa -wb -a %s -b %s | ",
          "awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'"
        )
        
        Peaks_at_all_Origins = read.table(
          pipe(sprintf(command_3, peakFile, All_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with early origins
        Peaks_at_early_Origins = read.table(
          pipe(sprintf(command_3, peakFile, E_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with late origins
        Peaks_at_late_Origins = read.table(
          pipe(sprintf(command_3, peakFile, L_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # standardize the MACS2 peak columns
        allPeaks = dplyr::rename(allPeaks,
                                 chrom=chr,
                                 peakStart=start,
                                 peakEnd=end,
                                 peakLength=length,
                                 peakSummit=abs_summit)
        
        # peaks not associated with known origins
        Peaks_at_no_Origins = dplyr::anti_join(
          allPeaks,
          Peaks_at_all_Origins,
          by=c("chrom", "peakStart", "peakEnd")
        )
        
        # remove duplicate peak entries
        Peaks_at_all_Origins = Peaks_at_all_Origins[ !duplicated(Peaks_at_all_Origins$peakSummit), ]
        
        Peaks_at_early_Origins = Peaks_at_early_Origins[ !duplicated(Peaks_at_early_Origins$peakSummit), ]
        
        Peaks_at_late_Origins = Peaks_at_late_Origins[ !duplicated(Peaks_at_late_Origins$peakSummit), ]
        
        Peaks_at_no_Origins = Peaks_at_no_Origins[ !duplicated(Peaks_at_no_Origins$peakSummit), ]
        
        # print peak summary
        message("There are", " ",
                nrow(allPeaks), " genome-wide BrDU peaks", "\n",
                nrow(Peaks_at_all_Origins), " were at known origins", "\n",
                nrow(Peaks_at_early_Origins), " were at early origins", "\n",
                nrow(Peaks_at_late_Origins), " were at late origins and", "\n",
                nrow(Peaks_at_no_Origins), " were at non-origin positions.")
        
        # save all five peak classes
        write.table(
          allPeaks[,1:5],
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_Genomewide_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_all_Origins,
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_Origin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_early_Origins,
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_EarlyOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_late_Origins,
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_LateOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_no_Origins[,1:5],
          file=paste0(dir, Pro_1, "/Peaks/", Pro_1, "_NonOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        unlink(c(IP, Input, peakFile), recursive=TRUE, force=TRUE)
        unlink(outDir, recursive=TRUE, force=TRUE)
      }
      
      BrDU_Peaks_Processing(IPBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_BrDU.bam"),
                            InBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input.bam"))
      
    }
    
  }
  
  if(Alignment == "malign"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Peaks_ma"))){
      
      BrDU_Peaks_Processing <- function(IPBam, InBam){
        
        message(paste0("➤ Processing BrDU peaks for ", Pro_1))
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Peaks_ma")))
        
        tempdir(check=TRUE)
        
        IP = tempfile(fileext=".bed")
        Input = tempfile(fileext=".bed")
        peakFile = tempfile(fileext=".bed")
        outDir = tempfile(pattern="MACS2_")
        
        suppressWarnings(dir.create(outDir))
        
        # convert BrDU and Input BAM files to BED
        command_1 = "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s > %s"
        
        system(sprintf(command_1, IPBam, IP))
        system(sprintf(command_1, InBam, Input))
        
        # call genome-wide BrDU peaks
        command_2 = "macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n Peak --outdir %s 2> /dev/null"
        
        system(sprintf(command_2, IP, Input, outDir))
        
        allPeaks = read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#")
        
        write.table(allPeaks, file=peakFile, quote=FALSE, row.names=FALSE,
                    sep="\t", col.names=FALSE)
        
        ColHeads = "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\""
        
        # peaks associated with all known origins
        command_3 = paste0(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools intersect ",
          "-wa -wb -a %s -b %s | ",
          "awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'"
        )
        
        Peaks_at_all_Origins = read.table(
          pipe(sprintf(command_3, peakFile, All_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with early origins
        Peaks_at_early_Origins = read.table(
          pipe(sprintf(command_3, peakFile, E_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with late origins
        Peaks_at_late_Origins = read.table(
          pipe(sprintf(command_3, peakFile, L_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # standardize the MACS2 peak columns
        allPeaks = dplyr::rename(allPeaks,
                                 chrom=chr,
                                 peakStart=start,
                                 peakEnd=end,
                                 peakLength=length,
                                 peakSummit=abs_summit)
        
        # peaks not associated with known origins
        Peaks_at_no_Origins = dplyr::anti_join(
          allPeaks,
          Peaks_at_all_Origins,
          by=c("chrom", "peakStart", "peakEnd")
        )
        
        # remove duplicate peak entries
        Peaks_at_all_Origins = Peaks_at_all_Origins[ !duplicated(Peaks_at_all_Origins$peakSummit), ]
        
        Peaks_at_early_Origins = Peaks_at_early_Origins[ !duplicated(Peaks_at_early_Origins$peakSummit), ]
        
        Peaks_at_late_Origins = Peaks_at_late_Origins[ !duplicated(Peaks_at_late_Origins$peakSummit), ]
        
        Peaks_at_no_Origins = Peaks_at_no_Origins[ !duplicated(Peaks_at_no_Origins$peakSummit), ]
        
        # print peak summary
        message("There are", " ",
                nrow(allPeaks), " genome-wide BrDU peaks", "\n",
                nrow(Peaks_at_all_Origins), " were at known origins", "\n",
                nrow(Peaks_at_early_Origins), " were at early origins", "\n",
                nrow(Peaks_at_late_Origins), " were at late origins and", "\n",
                nrow(Peaks_at_no_Origins), " were at non-origin positions.")
        
        # save all five peak classes
        write.table(
          allPeaks[,1:5],
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_Genomewide_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_all_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_Origin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_early_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_EarlyOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_late_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_LateOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_no_Origins[,1:5],
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_NonOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        unlink(c(IP, Input, peakFile), recursive=TRUE, force=TRUE)
        unlink(outDir, recursive=TRUE, force=TRUE)
      }
      
      BrDU_Peaks_Processing(IPBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_BrDU.bam"),
                            InBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input.bam"))
      
    }
    
  }
  
  if(Alignment == "mrdna"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Peaks_ma"))){
      
      BrDU_Peaks_Processing <- function(IPBam, InBam){
        
        message(paste0("➤ Processing BrDU peaks for ", Pro_1))
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Peaks_ma")))
        
        tempdir(check=TRUE)
        
        IP = tempfile(fileext=".bed")
        Input = tempfile(fileext=".bed")
        peakFile = tempfile(fileext=".bed")
        outDir = tempfile(pattern="MACS2_")
        
        suppressWarnings(dir.create(outDir))
        
        # convert BrDU and Input BAM files to BED
        command_1 = "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools bamtobed -i %s > %s"
        
        system(sprintf(command_1, IPBam, IP))
        system(sprintf(command_1, InBam, Input))
        
        # call genome-wide BrDU peaks
        command_2 = "macs2 callpeak -t %s -c %s -f BED -g 12157105 -p 10e-6 --nomodel -n Peak --outdir %s 2> /dev/null"
        
        system(sprintf(command_2, IP, Input, outDir))
        
        allPeaks = read.delim2(paste0(outDir, "/Peak_peaks.xls"), comment.char="#")
        
        write.table(allPeaks, file=peakFile, quote=FALSE, row.names=FALSE,
                    sep="\t", col.names=FALSE)
        
        ColHeads = "\"chrom\\tpeakStart\\tpeakEnd\\tpeakLength\\tpeakSummit\\toriName\\toriStart\\toriEnd\""
        
        # peaks associated with all known origins
        command_3 = paste0(
          "/Applications/ngsAnalyser.app/Contents/Resources/app/bedtools2/bin/bedtools intersect ",
          "-wa -wb -a %s -b %s | ",
          "awk 'BEGIN {print %s} {OFS=\"\\t\"} {print $1,$2,$3,$4,$5,$14,$12,$13}'"
        )
        
        Peaks_at_all_Origins = read.table(
          pipe(sprintf(command_3, peakFile, All_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with early origins
        Peaks_at_early_Origins = read.table(
          pipe(sprintf(command_3, peakFile, E_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # peaks associated with late origins
        Peaks_at_late_Origins = read.table(
          pipe(sprintf(command_3, peakFile, L_Ori_Link, ColHeads)),
          header=TRUE
        )
        
        # standardize the MACS2 peak columns
        allPeaks = dplyr::rename(allPeaks,
                                 chrom=chr,
                                 peakStart=start,
                                 peakEnd=end,
                                 peakLength=length,
                                 peakSummit=abs_summit)
        
        # peaks not associated with known origins
        Peaks_at_no_Origins = dplyr::anti_join(
          allPeaks,
          Peaks_at_all_Origins,
          by=c("chrom", "peakStart", "peakEnd")
        )
        
        # remove duplicate peak entries
        Peaks_at_all_Origins = Peaks_at_all_Origins[ !duplicated(Peaks_at_all_Origins$peakSummit), ]
        
        Peaks_at_early_Origins = Peaks_at_early_Origins[ !duplicated(Peaks_at_early_Origins$peakSummit), ]
        
        Peaks_at_late_Origins = Peaks_at_late_Origins[ !duplicated(Peaks_at_late_Origins$peakSummit), ]
        
        Peaks_at_no_Origins = Peaks_at_no_Origins[ !duplicated(Peaks_at_no_Origins$peakSummit), ]
        
        # print peak summary
        message("There are", " ",
                nrow(allPeaks), " genome-wide BrDU peaks", "\n",
                nrow(Peaks_at_all_Origins), " were at known origins", "\n",
                nrow(Peaks_at_early_Origins), " were at early origins", "\n",
                nrow(Peaks_at_late_Origins), " were at late origins and", "\n",
                nrow(Peaks_at_no_Origins), " were at non-origin positions.")
        
        # save all five peak classes
        write.table(
          allPeaks[,1:5],
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_Genomewide_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_all_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_Origin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_early_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_EarlyOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_late_Origins,
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_LateOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        write.table(
          Peaks_at_no_Origins[,1:5],
          file=paste0(dir, Pro_1, "/Peaks_ma/", Pro_1, "_NonOrigin_Peaks.bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
        
        unlink(c(IP, Input, peakFile), recursive=TRUE, force=TRUE)
        unlink(outDir, recursive=TRUE, force=TRUE)
      }
      
      BrDU_Peaks_Processing(IPBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_BrDU_ma.bam"),
                            InBam=paste0(dir, Pro_1, "/", BamFolder, "/", Pro_1, "_Input_ma.bam"))
      
    }
    
  }
  
  
  ## Calculate Ratio
  #
  message("⏳ Calculating enrichment ratios...")
  #
  
  if(Alignment == "generic"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Ratios"))){
      
      CalculateRatio <- function(IP_coverage, Input_coverage, RatioName,
                                 NoiseChunkSizeBp=2000, NoiseIterations=2500,
                                 NoiseSeed=123, NoiseSmoothingSpar=0.65,
                                 NoiseFloor=1e-6){
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Ratios")))
        
        # read and optionally collapse coverage files
        ReadCoverage <- function(coverageFiles){
          
          coverage = read.table(coverageFiles[1], header=FALSE)
          
          if(length(coverageFiles) > 1){
            
            for(j in 2:length(coverageFiles)){
              
              additionalCoverage = read.table(coverageFiles[j], header=FALSE)
              
              coordinatesMatch = identical(as.character(coverage[,1]),
                                           as.character(additionalCoverage[,1])) &&
                identical(as.numeric(coverage[,2]), as.numeric(additionalCoverage[,2])) &&
                identical(as.numeric(coverage[,3]), as.numeric(additionalCoverage[,3]))
              
              if(coordinatesMatch == FALSE){
                stop("Coverage coordinates do not match for strand collapsing.")
              }
              
              coverage[,5] = as.numeric(coverage[,5]) + as.numeric(additionalCoverage[,5])
            }
          }
          
          coverage
        }
        
        IP.df = ReadCoverage(IP_coverage)
        In.df = ReadCoverage(Input_coverage)
        
        coordinatesMatch = identical(as.character(IP.df[,1]), as.character(In.df[,1])) &&
          identical(as.numeric(IP.df[,2]), as.numeric(In.df[,2])) &&
          identical(as.numeric(IP.df[,3]), as.numeric(In.df[,3]))
        
        if(coordinatesMatch == FALSE){
          stop("BrDU and Input coverage coordinates do not match.")
        }
        
        IP.df[,5] = as.numeric(IP.df[,5])
        In.df[,5] = as.numeric(In.df[,5])
        
        # read genome-wide BrDU peaks
        PeakFile = paste0(dir, Pro_1, "/", "Peaks", "/",
                          Pro_1, "_Genomewide_Peaks.bed")
        
        if(!file.exists(PeakFile)){
          stop("Genome-wide peak file is missing: ", PeakFile)
        }
        
        Peaks = read.table(PeakFile, header=TRUE)
        
        Peaks$peakStart = as.numeric(Peaks$peakStart)
        Peaks$peakEnd = as.numeric(Peaks$peakEnd)
        
        # merge overlapping peak intervals
        MergeIntervals <- function(intervals){
          
          if(nrow(intervals) == 0){
            return(intervals)
          }
          
          intervals = intervals[order(intervals$start, intervals$end), ]
          merged = intervals[1, , drop=FALSE]
          
          if(nrow(intervals) > 1){
            
            for(j in 2:nrow(intervals)){
              
              last = nrow(merged)
              
              if(intervals$start[j] <= merged$end[last]){
                
                merged$end[last] = max(merged$end[last], intervals$end[j])
                
              } else {
                
                merged = rbind(merged, intervals[j, , drop=FALSE])
              }
            }
          }
          
          merged
        }
        
        # construct non-peak regions
        BuildNonPeakGaps <- function(chrom, chromLength){
          
          Peaks_chr = Peaks[Peaks$chrom == chrom, , drop=FALSE]
          
          if(nrow(Peaks_chr) == 0){
            
            return(data.frame(start=0, end=chromLength, length=chromLength))
          }
          
          intervals = data.frame(start=pmax(0, Peaks_chr$peakStart),
                                 end=pmin(chromLength, Peaks_chr$peakEnd))
          
          intervals = intervals[
            is.finite(intervals$start) &
              is.finite(intervals$end) &
              intervals$end > intervals$start, ]
          
          if(nrow(intervals) == 0){
            
            return(data.frame(start=0, end=chromLength, length=chromLength))
          }
          
          intervals = MergeIntervals(intervals)
          
          gaps = data.frame(start=numeric(), end=numeric())
          cursor = 0
          
          for(j in 1:nrow(intervals)){
            
            if(intervals$start[j] > cursor){
              
              gaps = rbind(gaps, data.frame(start=cursor,
                                            end=intervals$start[j]))
            }
            
            cursor = max(cursor, intervals$end[j])
          }
          
          if(cursor < chromLength){
            
            gaps = rbind(gaps, data.frame(start=cursor, end=chromLength))
          }
          
          if(nrow(gaps) == 0){
            
            return(data.frame(start=numeric(), end=numeric(),
                              length=numeric()))
          }
          
          gaps$length = gaps$end - gaps$start
          gaps = gaps[gaps$length >= 200, , drop=FALSE]
          
          gaps
        }
        
        SummariseValues <- function(values){
          
          values = values[is.finite(values)]
          
          if(length(values) == 0){
            return(NA_real_)
          }
          
          median(values)
        }
        
        # predict chromosome-position-dependent noise
        PredictNoise <- function(binCenters, sampleCenters, sampleValues,
                                 fallbackValues){
          
          sampleOK = is.finite(sampleCenters) & is.finite(sampleValues)
          
          sampleCenters = sampleCenters[sampleOK]
          sampleValues = sampleValues[sampleOK]
          
          positiveBackground = fallbackValues[
            is.finite(fallbackValues) & fallbackValues > NoiseFloor
          ]
          
          if(length(positiveBackground) > 0){
            
            backgroundFloor = max(
              NoiseFloor,
              as.numeric(quantile(positiveBackground, probs=0.01,
                                  na.rm=TRUE, names=FALSE))
            )
            
          } else {
            
            backgroundFloor = NoiseFloor
          }
          
          fallback = median(fallbackValues, na.rm=TRUE)
          
          if(!is.finite(fallback) || fallback < backgroundFloor){
            fallback = backgroundFloor
          }
          
          if(length(sampleValues) < 2 ||
             length(unique(sampleCenters)) < 2){
            
            return(rep(fallback, length(binCenters)))
          }
          
          # background smoothing is performed on the log scale
          sampleValues = log(pmax(sampleValues, backgroundFloor))
          
          sample.df = data.frame(center=round(sampleCenters),
                                 value=sampleValues)
          
          sample.df = aggregate(value~center, data=sample.df,
                                FUN=function(x) median(x, na.rm=TRUE))
          
          sample.df = sample.df[order(sample.df$center), ]
          
          if(nrow(sample.df) < 4 ||
             length(unique(sample.df$value)) < 2){
            
            prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                rule=2, ties="ordered")$y
            
          } else {
            
            splineObject = tryCatch(
              smooth.spline(sample.df$center, sample.df$value,
                            spar=NoiseSmoothingSpar),
              error=function(e) NULL
            )
            
            if(is.null(splineObject)){
              
              prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                  rule=2, ties="ordered")$y
              
            } else {
              
              prediction = predict(splineObject, binCenters)$y
            }
          }
          
          prediction[!is.finite(prediction)] = log(fallback)
          
          prediction = pmin(
            pmax(prediction, min(sample.df$value)),
            max(sample.df$value)
          )
          
          prediction = exp(prediction)
          prediction[prediction < backgroundFloor] = backgroundFloor
          
          prediction
        }
        
        # estimate BrDU and Input noise separately
        EstimateNoise <- function(){
          
          set.seed(NoiseSeed)
          
          BrDU_noise = numeric(nrow(IP.df))
          Input_noise = numeric(nrow(In.df))
          
          chroms = unique(as.character(IP.df[,1]))
          
          for(chr in chroms){
            
            chrIndex = which(IP.df[,1] == chr)
            
            chrStart = as.numeric(IP.df[chrIndex,2])
            chrEnd = as.numeric(IP.df[chrIndex,3])
            chrCenters = (chrStart + chrEnd)/2
            chrLength = max(chrEnd, na.rm=TRUE)
            
            gaps = BuildNonPeakGaps(chr, chrLength)
            
            if(nrow(gaps) == 0){
              
              BrDU_noise[chrIndex] = max(NoiseFloor,
                                         median(IP.df[chrIndex,5], na.rm=TRUE))
              
              Input_noise[chrIndex] = max(NoiseFloor,
                                          median(In.df[chrIndex,5], na.rm=TRUE))
              
              next
            }
            
            gapProbability = gaps$length/sum(gaps$length)
            
            sampledGapIndex = sample(1:nrow(gaps), size=NoiseIterations,
                                     replace=TRUE, prob=gapProbability)
            
            sampleCenters = numeric(NoiseIterations)
            sampledBrDU = numeric(NoiseIterations)
            sampledInput = numeric(NoiseIterations)
            
            for(j in 1:NoiseIterations){
              
              gap = gaps[sampledGapIndex[j], ]
              chunkSize = min(NoiseChunkSizeBp, gap$length)
              
              if(gap$length > chunkSize){
                
                chunkStart = runif(1, min=gap$start,
                                   max=gap$end-chunkSize)
                
              } else {
                
                chunkStart = gap$start
              }
              
              chunkEnd = chunkStart + chunkSize
              
              chunkIndex = which(chrCenters >= chunkStart &
                                   chrCenters <= chunkEnd)
              
              sampleCenters[j] = (chunkStart + chunkEnd)/2
              
              sampledBrDU[j] = SummariseValues(
                IP.df[chrIndex[chunkIndex],5]
              )
              
              sampledInput[j] = SummariseValues(
                In.df[chrIndex[chunkIndex],5]
              )
            }
            
            nonPeakIndex = rep(FALSE, length(chrIndex))
            
            for(j in 1:nrow(gaps)){
              
              nonPeakIndex = nonPeakIndex |
                (chrCenters >= gaps$start[j] &
                   chrCenters <= gaps$end[j])
            }
            
            BrDU_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledBrDU,
              fallbackValues=IP.df[chrIndex[nonPeakIndex],5]
            )
            
            Input_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledInput,
              fallbackValues=In.df[chrIndex[nonPeakIndex],5]
            )
          }
          
          list(BrDU_noise=BrDU_noise, Input_noise=Input_noise)
        }
        
        noise = EstimateNoise()
        
        # library-size normalization
        IP_Sum = sum(IP.df[,5], na.rm=TRUE)
        In_Sum = sum(In.df[,5], na.rm=TRUE)
        
        if(!is.finite(In_Sum) || In_Sum <= 0){
          stop("Input coverage sum is zero.")
        }
        
        corrFactor = IP_Sum/In_Sum
        
        In.score.norm = In.df[,5]*corrFactor
        In.noise.norm = noise$Input_noise*corrFactor
        
        # calculate the three enrichment ratios
        Ratio.ipin = IP.df[,5]/In.score.norm
        
        Ratio.ipnoise = IP.df[,5]/noise$BrDU_noise
        
        Ratio.ipin.noise = Ratio.ipnoise/
          (In.score.norm/In.noise.norm)
        
        Ratio.ipin[!is.finite(Ratio.ipin)] = 0
        Ratio.ipnoise[!is.finite(Ratio.ipnoise)] = 0
        Ratio.ipin.noise[!is.finite(Ratio.ipin.noise)] = 0
        
        # Poisson enrichment probability
        pvalue = ppois(q=IP.df[,5]-1, lambda=In.score.norm,
                       lower.tail=FALSE, log=FALSE)
        
        ratio.df = data.frame(
          chrom=IP.df[,1],
          chromStart=IP.df[,2],
          chromEnd=IP.df[,3],
          name=RatioName,
          ip.score=IP.df[,5],
          in.score=round(In.score.norm, 4),
          ip.noise=round(noise$BrDU_noise, 4),
          in.noise=round(In.noise.norm, 4),
          ratio.ipin=round(Ratio.ipin, 4),
          ratio.ipnoise=round(Ratio.ipnoise, 4),
          ratio.ipin.noise=round(Ratio.ipin.noise, 4),
          pvalue=pvalue
        )
        
        write.table(
          ratio.df,
          file=paste0(dir, Pro_1, "/", "Ratios", "/", RatioName, ".bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
      }
      
      # Watson-strand ratio
      CalculateRatio(IP_coverage=paste0(dir, Pro_1, "/Coverage/", Pro_1, "_BrDU_watson.bed"),
                     Input_coverage=paste0(dir, Pro_1, "/Coverage/", Pro_1, "_Input_watson.bed"),
                     RatioName=paste0(Pro_1, "_BrDU_watson"))
      
      # Crick-strand ratio
      CalculateRatio(IP_coverage=paste0(dir, Pro_1, "/Coverage/", Pro_1, "_BrDU_crick.bed"),
                     Input_coverage=paste0(dir, Pro_1, "/Coverage/", Pro_1, "_Input_crick.bed"),
                     RatioName=paste0(Pro_1, "_BrDU_crick"))
      
      # strand-collapsed ratio
      CalculateRatio(IP_coverage=c(paste0(dir, Pro_1, "/Coverage/", Pro_1, "_BrDU_watson.bed"),
                                   paste0(dir, Pro_1, "/Coverage/", Pro_1, "_BrDU_crick.bed")),
                     Input_coverage=c(paste0(dir, Pro_1, "/Coverage/", Pro_1, "_Input_watson.bed"),
                                      paste0(dir, Pro_1, "/Coverage/", Pro_1, "_Input_crick.bed")),
                     RatioName=paste0(Pro_1, "_BrDU_collapsed"))
    }
  }
  
  if(Alignment == "malign"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Ratios_ma"))){
      
      CalculateRatio <- function(IP_coverage, Input_coverage, RatioName,
                                 NoiseChunkSizeBp=2000, NoiseIterations=2500,
                                 NoiseSeed=123, NoiseSmoothingSpar=0.65,
                                 NoiseFloor=1e-6){
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Ratios_ma")))
        
        # read and optionally collapse coverage files
        ReadCoverage <- function(coverageFiles){
          
          coverage = read.table(coverageFiles[1], header=FALSE)
          
          if(length(coverageFiles) > 1){
            
            for(j in 2:length(coverageFiles)){
              
              additionalCoverage = read.table(coverageFiles[j], header=FALSE)
              
              coordinatesMatch = identical(as.character(coverage[,1]),
                                           as.character(additionalCoverage[,1])) &&
                identical(as.numeric(coverage[,2]), as.numeric(additionalCoverage[,2])) &&
                identical(as.numeric(coverage[,3]), as.numeric(additionalCoverage[,3]))
              
              if(coordinatesMatch == FALSE){
                stop("Coverage coordinates do not match for strand collapsing.")
              }
              
              coverage[,5] = as.numeric(coverage[,5]) + as.numeric(additionalCoverage[,5])
            }
          }
          
          coverage
        }
        
        IP.df = ReadCoverage(IP_coverage)
        In.df = ReadCoverage(Input_coverage)
        
        coordinatesMatch = identical(as.character(IP.df[,1]), as.character(In.df[,1])) &&
          identical(as.numeric(IP.df[,2]), as.numeric(In.df[,2])) &&
          identical(as.numeric(IP.df[,3]), as.numeric(In.df[,3]))
        
        if(coordinatesMatch == FALSE){
          stop("BrDU and Input coverage coordinates do not match.")
        }
        
        IP.df[,5] = as.numeric(IP.df[,5])
        In.df[,5] = as.numeric(In.df[,5])
        
        # read genome-wide BrDU peaks
        PeakFile = paste0(dir, Pro_1, "/", "Peaks_ma", "/", Pro_1, "_Genomewide_Peaks.bed")
        
        if(!file.exists(PeakFile)){
          stop("Genome-wide peak file is missing: ", PeakFile)
        }
        
        Peaks = read.table(PeakFile, header=TRUE)
        
        Peaks$peakStart = as.numeric(Peaks$peakStart)
        Peaks$peakEnd = as.numeric(Peaks$peakEnd)
        
        # merge overlapping peak intervals
        MergeIntervals <- function(intervals){
          
          if(nrow(intervals) == 0){
            return(intervals)
          }
          
          intervals = intervals[order(intervals$start, intervals$end), ]
          merged = intervals[1, , drop=FALSE]
          
          if(nrow(intervals) > 1){
            
            for(j in 2:nrow(intervals)){
              
              last = nrow(merged)
              
              if(intervals$start[j] <= merged$end[last]){
                
                merged$end[last] = max(merged$end[last], intervals$end[j])
                
              } else {
                
                merged = rbind(merged, intervals[j, , drop=FALSE])
              }
            }
          }
          
          merged
        }
        
        # construct non-peak regions
        BuildNonPeakGaps <- function(chrom, chromLength){
          
          Peaks_chr = Peaks[Peaks$chrom == chrom, , drop=FALSE]
          
          if(nrow(Peaks_chr) == 0){
            
            return(data.frame(start=0, end=chromLength, length=chromLength))
          }
          
          intervals = data.frame(start=pmax(0, Peaks_chr$peakStart),
                                 end=pmin(chromLength, Peaks_chr$peakEnd))
          
          intervals = intervals[
            is.finite(intervals$start) &
              is.finite(intervals$end) &
              intervals$end > intervals$start, ]
          
          if(nrow(intervals) == 0){
            
            return(data.frame(start=0, end=chromLength, length=chromLength))
          }
          
          intervals = MergeIntervals(intervals)
          
          gaps = data.frame(start=numeric(), end=numeric())
          cursor = 0
          
          for(j in 1:nrow(intervals)){
            
            if(intervals$start[j] > cursor){
              
              gaps = rbind(gaps, data.frame(start=cursor,
                                            end=intervals$start[j]))
            }
            
            cursor = max(cursor, intervals$end[j])
          }
          
          if(cursor < chromLength){
            
            gaps = rbind(gaps, data.frame(start=cursor, end=chromLength))
          }
          
          if(nrow(gaps) == 0){
            
            return(data.frame(start=numeric(), end=numeric(),
                              length=numeric()))
          }
          
          gaps$length = gaps$end - gaps$start
          gaps = gaps[gaps$length >= 200, , drop=FALSE]
          
          gaps
        }
        
        SummariseValues <- function(values){
          
          values = values[is.finite(values)]
          
          if(length(values) == 0){
            return(NA_real_)
          }
          
          median(values)
        }
        
        # predict chromosome-position-dependent noise
        PredictNoise <- function(binCenters, sampleCenters, sampleValues,
                                 fallbackValues){
          
          sampleOK = is.finite(sampleCenters) & is.finite(sampleValues)
          
          sampleCenters = sampleCenters[sampleOK]
          sampleValues = sampleValues[sampleOK]
          
          positiveBackground = fallbackValues[
            is.finite(fallbackValues) & fallbackValues > NoiseFloor
          ]
          
          if(length(positiveBackground) > 0){
            
            backgroundFloor = max(
              NoiseFloor,
              as.numeric(quantile(positiveBackground, probs=0.01,
                                  na.rm=TRUE, names=FALSE))
            )
            
          } else {
            
            backgroundFloor = NoiseFloor
          }
          
          fallback = median(fallbackValues, na.rm=TRUE)
          
          if(!is.finite(fallback) || fallback < backgroundFloor){
            fallback = backgroundFloor
          }
          
          if(length(sampleValues) < 2 ||
             length(unique(sampleCenters)) < 2){
            
            return(rep(fallback, length(binCenters)))
          }
          
          # background smoothing is performed on the log scale
          sampleValues = log(pmax(sampleValues, backgroundFloor))
          
          sample.df = data.frame(center=round(sampleCenters),
                                 value=sampleValues)
          
          sample.df = aggregate(value~center, data=sample.df,
                                FUN=function(x) median(x, na.rm=TRUE))
          
          sample.df = sample.df[order(sample.df$center), ]
          
          if(nrow(sample.df) < 4 ||
             length(unique(sample.df$value)) < 2){
            
            prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                rule=2, ties="ordered")$y
            
          } else {
            
            splineObject = tryCatch(
              smooth.spline(sample.df$center, sample.df$value,
                            spar=NoiseSmoothingSpar),
              error=function(e) NULL
            )
            
            if(is.null(splineObject)){
              
              prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                  rule=2, ties="ordered")$y
              
            } else {
              
              prediction = predict(splineObject, binCenters)$y
            }
          }
          
          prediction[!is.finite(prediction)] = log(fallback)
          
          prediction = pmin(
            pmax(prediction, min(sample.df$value)),
            max(sample.df$value)
          )
          
          prediction = exp(prediction)
          prediction[prediction < backgroundFloor] = backgroundFloor
          
          prediction
        }
        
        # estimate BrDU and Input noise separately
        EstimateNoise <- function(){
          
          set.seed(NoiseSeed)
          
          BrDU_noise = numeric(nrow(IP.df))
          Input_noise = numeric(nrow(In.df))
          
          chroms = unique(as.character(IP.df[,1]))
          
          for(chr in chroms){
            
            chrIndex = which(IP.df[,1] == chr)
            
            chrStart = as.numeric(IP.df[chrIndex,2])
            chrEnd = as.numeric(IP.df[chrIndex,3])
            chrCenters = (chrStart + chrEnd)/2
            chrLength = max(chrEnd, na.rm=TRUE)
            
            gaps = BuildNonPeakGaps(chr, chrLength)
            
            if(nrow(gaps) == 0){
              
              BrDU_noise[chrIndex] = max(NoiseFloor,
                                         median(IP.df[chrIndex,5], na.rm=TRUE))
              
              Input_noise[chrIndex] = max(NoiseFloor,
                                          median(In.df[chrIndex,5], na.rm=TRUE))
              
              next
            }
            
            gapProbability = gaps$length/sum(gaps$length)
            
            sampledGapIndex = sample(1:nrow(gaps), size=NoiseIterations,
                                     replace=TRUE, prob=gapProbability)
            
            sampleCenters = numeric(NoiseIterations)
            sampledBrDU = numeric(NoiseIterations)
            sampledInput = numeric(NoiseIterations)
            
            for(j in 1:NoiseIterations){
              
              gap = gaps[sampledGapIndex[j], ]
              chunkSize = min(NoiseChunkSizeBp, gap$length)
              
              if(gap$length > chunkSize){
                
                chunkStart = runif(1, min=gap$start,
                                   max=gap$end-chunkSize)
                
              } else {
                
                chunkStart = gap$start
              }
              
              chunkEnd = chunkStart + chunkSize
              
              chunkIndex = which(chrCenters >= chunkStart &
                                   chrCenters <= chunkEnd)
              
              sampleCenters[j] = (chunkStart + chunkEnd)/2
              
              sampledBrDU[j] = SummariseValues(
                IP.df[chrIndex[chunkIndex],5]
              )
              
              sampledInput[j] = SummariseValues(
                In.df[chrIndex[chunkIndex],5]
              )
            }
            
            nonPeakIndex = rep(FALSE, length(chrIndex))
            
            for(j in 1:nrow(gaps)){
              
              nonPeakIndex = nonPeakIndex |
                (chrCenters >= gaps$start[j] &
                   chrCenters <= gaps$end[j])
            }
            
            BrDU_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledBrDU,
              fallbackValues=IP.df[chrIndex[nonPeakIndex],5]
            )
            
            Input_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledInput,
              fallbackValues=In.df[chrIndex[nonPeakIndex],5]
            )
          }
          
          list(BrDU_noise=BrDU_noise, Input_noise=Input_noise)
        }
        
        noise = EstimateNoise()
        
        # library-size normalization
        IP_Sum = sum(IP.df[,5], na.rm=TRUE)
        In_Sum = sum(In.df[,5], na.rm=TRUE)
        
        if(!is.finite(In_Sum) || In_Sum <= 0){
          stop("Input coverage sum is zero.")
        }
        
        corrFactor = IP_Sum/In_Sum
        
        In.score.norm = In.df[,5]*corrFactor
        In.noise.norm = noise$Input_noise*corrFactor
        
        # calculate the three enrichment ratios
        Ratio.ipin = IP.df[,5]/In.score.norm
        
        Ratio.ipnoise = IP.df[,5]/noise$BrDU_noise
        
        Ratio.ipin.noise = Ratio.ipnoise/(In.score.norm/In.noise.norm)
        
        Ratio.ipin[!is.finite(Ratio.ipin)] = 0
        Ratio.ipnoise[!is.finite(Ratio.ipnoise)] = 0
        Ratio.ipin.noise[!is.finite(Ratio.ipin.noise)] = 0
        
        
        ratio.df = data.frame(
          chrom=IP.df[,1],
          chromStart=IP.df[,2],
          chromEnd=IP.df[,3],
          name=RatioName,
          ip.score=IP.df[,5],
          in.score=round(In.score.norm, 4),
          ip.noise=round(noise$BrDU_noise, 4),
          in.noise=round(In.noise.norm, 4),
          ratio.ipin=round(Ratio.ipin, 4),
          ratio.ipnoise=round(Ratio.ipnoise, 4),
          ratio.ipin.noise=round(Ratio.ipin.noise, 4)
        )
        
        write.table(
          ratio.df,
          file=paste0(dir, Pro_1, "/", "Ratios_ma", "/", RatioName, ".bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
      }
      
      # Watson-strand ratio
      CalculateRatio(IP_coverage=paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_BrDU_watson.bed"),
                     Input_coverage=paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_Input_watson.bed"),
                     RatioName=paste0(Pro_1, "_BrDU_watson"))
      
      # Crick-strand ratio
      CalculateRatio(IP_coverage=paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_BrDU_crick.bed"),
                     Input_coverage=paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_Input_crick.bed"),
                     RatioName=paste0(Pro_1, "_BrDU_crick"))
      
      # strand-collapsed ratio
      CalculateRatio(IP_coverage=c(paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_BrDU_watson.bed"),
                                   paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_BrDU_crick.bed")),
                     Input_coverage=c(paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_Input_watson.bed"),
                                      paste0(dir, Pro_1, "/Coverage_ma/", Pro_1, "_Input_crick.bed")),
                     RatioName=paste0(Pro_1, "_BrDU_collapsed"))
    }
    
  }
  
  if(Alignment == "mrdna"){
    
    if(!dir.exists(paste0(dir, Pro_1, "/", "Ratios_ma_rdna"))){
      
      CalculateRatio <- function(IP_coverage, Input_coverage,
                                 Noise_IP_coverage, Noise_Input_coverage,
                                 RatioName, NoiseChunkSizeBp=2000,
                                 NoiseIterations=2500, NoiseSeed=123,
                                 NoiseSmoothingSpar=0.65, NoiseFloor=1e-6){
        
        suppressWarnings(dir.create(paste0(dir, Pro_1, "/", "Ratios_ma_rdna")))
        
        # read and optionally collapse coverage files
        ReadCoverage <- function(coverageFiles){
          
          coverage = read.table(coverageFiles[1], header=FALSE)
          
          if(length(coverageFiles) > 1){
            
            for(j in 2:length(coverageFiles)){
              
              additionalCoverage = read.table(coverageFiles[j], header=FALSE)
              
              coordinatesMatch = identical(as.character(coverage[,1]),
                                           as.character(additionalCoverage[,1])) &&
                identical(as.numeric(coverage[,2]), as.numeric(additionalCoverage[,2])) &&
                identical(as.numeric(coverage[,3]), as.numeric(additionalCoverage[,3]))
              
              if(coordinatesMatch == FALSE){
                stop("Coverage coordinates do not match for strand collapsing.")
              }
              
              coverage[,5] = as.numeric(coverage[,5]) +
                as.numeric(additionalCoverage[,5])
            }
          }
          
          coverage
        }
        
        # rDNA-reference coverage used for final ratio calculation
        IP.df = ReadCoverage(IP_coverage)
        In.df = ReadCoverage(Input_coverage)
        
        coordinatesMatch = identical(as.character(IP.df[,1]),
                                     as.character(In.df[,1])) &&
          identical(as.numeric(IP.df[,2]), as.numeric(In.df[,2])) &&
          identical(as.numeric(IP.df[,3]), as.numeric(In.df[,3]))
        
        if(coordinatesMatch == FALSE){
          stop("rDNA BrDU and Input coverage coordinates do not match.")
        }
        
        IP.df[,5] = as.numeric(IP.df[,5])
        In.df[,5] = as.numeric(In.df[,5])
        
        # whole-genome coverage used for background estimation
        Noise.IP.df = ReadCoverage(Noise_IP_coverage)
        Noise.In.df = ReadCoverage(Noise_Input_coverage)
        
        noiseCoordinatesMatch = identical(as.character(Noise.IP.df[,1]),
                                          as.character(Noise.In.df[,1])) &&
          identical(as.numeric(Noise.IP.df[,2]),
                    as.numeric(Noise.In.df[,2])) &&
          identical(as.numeric(Noise.IP.df[,3]),
                    as.numeric(Noise.In.df[,3]))
        
        if(noiseCoordinatesMatch == FALSE){
          stop("Whole-genome BrDU and Input coverage coordinates do not match.")
        }
        
        Noise.IP.df[,5] = as.numeric(Noise.IP.df[,5])
        Noise.In.df[,5] = as.numeric(Noise.In.df[,5])
        
        # read peaks obtained from whole-genome multiple alignments
        PeakFile = paste0(dir, Pro_1, "/", "Peaks_ma", "/",
                          Pro_1, "_Genomewide_Peaks.bed")
        
        if(!file.exists(PeakFile)){
          stop("Genome-wide peak file is missing: ", PeakFile)
        }
        
        Peaks = read.table(PeakFile, header=TRUE)
        
        Peaks$peakStart = as.numeric(Peaks$peakStart)
        Peaks$peakEnd = as.numeric(Peaks$peakEnd)
        
        # merge overlapping peak intervals
        MergeIntervals <- function(intervals){
          
          if(nrow(intervals) == 0){
            return(intervals)
          }
          
          intervals = intervals[order(intervals$start, intervals$end), ]
          merged = intervals[1, , drop=FALSE]
          
          if(nrow(intervals) > 1){
            
            for(j in 2:nrow(intervals)){
              
              last = nrow(merged)
              
              if(intervals$start[j] <= merged$end[last]){
                
                merged$end[last] = max(merged$end[last],
                                       intervals$end[j])
                
              } else {
                
                merged = rbind(merged, intervals[j, , drop=FALSE])
              }
            }
          }
          
          merged
        }
        
        # construct whole-genome non-peak regions
        BuildNonPeakGaps <- function(chrom, chromLength){
          
          Peaks_chr = Peaks[Peaks$chrom == chrom, , drop=FALSE]
          
          if(nrow(Peaks_chr) == 0){
            
            return(data.frame(start=0, end=chromLength,
                              length=chromLength))
          }
          
          intervals = data.frame(start=pmax(0, Peaks_chr$peakStart),
                                 end=pmin(chromLength, Peaks_chr$peakEnd))
          
          intervals = intervals[
            is.finite(intervals$start) &
              is.finite(intervals$end) &
              intervals$end > intervals$start, ]
          
          if(nrow(intervals) == 0){
            
            return(data.frame(start=0, end=chromLength,
                              length=chromLength))
          }
          
          intervals = MergeIntervals(intervals)
          
          gaps = data.frame(start=numeric(), end=numeric())
          cursor = 0
          
          for(j in 1:nrow(intervals)){
            
            if(intervals$start[j] > cursor){
              
              gaps = rbind(
                gaps,
                data.frame(start=cursor, end=intervals$start[j])
              )
            }
            
            cursor = max(cursor, intervals$end[j])
          }
          
          if(cursor < chromLength){
            
            gaps = rbind(
              gaps,
              data.frame(start=cursor, end=chromLength)
            )
          }
          
          if(nrow(gaps) == 0){
            
            return(data.frame(start=numeric(), end=numeric(),
                              length=numeric()))
          }
          
          gaps$length = gaps$end - gaps$start
          gaps = gaps[gaps$length >= 200, , drop=FALSE]
          
          gaps
        }
        
        SummariseValues <- function(values){
          
          values = values[is.finite(values)]
          
          if(length(values) == 0){
            return(NA_real_)
          }
          
          median(values)
        }
        
        # predict chromosome-position-dependent background
        PredictNoise <- function(binCenters, sampleCenters, sampleValues,
                                 fallbackValues){
          
          sampleOK = is.finite(sampleCenters) &
            is.finite(sampleValues)
          
          sampleCenters = sampleCenters[sampleOK]
          sampleValues = sampleValues[sampleOK]
          
          positiveBackground = fallbackValues[
            is.finite(fallbackValues) &
              fallbackValues > NoiseFloor
          ]
          
          if(length(positiveBackground) > 0){
            
            backgroundFloor = max(
              NoiseFloor,
              as.numeric(
                quantile(
                  positiveBackground,
                  probs=0.01,
                  na.rm=TRUE,
                  names=FALSE
                )
              )
            )
            
          } else {
            
            backgroundFloor = NoiseFloor
          }
          
          fallback = median(fallbackValues, na.rm=TRUE)
          
          if(!is.finite(fallback) || fallback < backgroundFloor){
            fallback = backgroundFloor
          }
          
          if(length(sampleValues) < 2 ||
             length(unique(sampleCenters)) < 2){
            
            return(rep(fallback, length(binCenters)))
          }
          
          # background smoothing is performed on the log scale
          sampleValues = log(pmax(sampleValues, backgroundFloor))
          
          sample.df = data.frame(center=round(sampleCenters),
                                 value=sampleValues)
          
          sample.df = aggregate(
            value~center,
            data=sample.df,
            FUN=function(x) median(x, na.rm=TRUE)
          )
          
          sample.df = sample.df[order(sample.df$center), ]
          
          if(nrow(sample.df) < 4 ||
             length(unique(sample.df$value)) < 2){
            
            prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                rule=2, ties="ordered")$y
            
          } else {
            
            splineObject = tryCatch(
              smooth.spline(
                sample.df$center,
                sample.df$value,
                spar=NoiseSmoothingSpar
              ),
              error=function(e) NULL
            )
            
            if(is.null(splineObject)){
              
              prediction = approx(sample.df$center, sample.df$value, xout=binCenters,
                                  rule=2, ties="ordered")$y
              
            } else {
              
              prediction = predict(splineObject, binCenters)$y
            }
          }
          
          prediction[!is.finite(prediction)] = log(fallback)
          
          prediction = pmin(
            pmax(prediction, min(sample.df$value)),
            max(sample.df$value)
          )
          
          prediction = exp(prediction)
          prediction[prediction < backgroundFloor] = backgroundFloor
          
          prediction
        }
        
        # estimate background from whole-genome non-peak coverage
        EstimateNoise <- function(){
          
          set.seed(NoiseSeed)
          
          BrDU_noise = numeric(nrow(Noise.IP.df))
          Input_noise = numeric(nrow(Noise.In.df))
          
          chroms = unique(as.character(Noise.IP.df[,1]))
          
          for(chr in chroms){
            
            chrIndex = which(Noise.IP.df[,1] == chr)
            
            chrStart = as.numeric(Noise.IP.df[chrIndex,2])
            chrEnd = as.numeric(Noise.IP.df[chrIndex,3])
            chrCenters = (chrStart + chrEnd)/2
            chrLength = max(chrEnd, na.rm=TRUE)
            
            gaps = BuildNonPeakGaps(chr, chrLength)
            
            if(nrow(gaps) == 0){
              
              BrDU_noise[chrIndex] = max(
                NoiseFloor,
                median(Noise.IP.df[chrIndex,5], na.rm=TRUE)
              )
              
              Input_noise[chrIndex] = max(
                NoiseFloor,
                median(Noise.In.df[chrIndex,5], na.rm=TRUE)
              )
              
              next
            }
            
            gapProbability = gaps$length/sum(gaps$length)
            
            sampledGapIndex = sample(
              1:nrow(gaps),
              size=NoiseIterations,
              replace=TRUE,
              prob=gapProbability
            )
            
            sampleCenters = numeric(NoiseIterations)
            sampledBrDU = numeric(NoiseIterations)
            sampledInput = numeric(NoiseIterations)
            
            for(j in 1:NoiseIterations){
              
              gap = gaps[sampledGapIndex[j], ]
              chunkSize = min(NoiseChunkSizeBp, gap$length)
              
              if(gap$length > chunkSize){
                
                chunkStart = runif(
                  1,
                  min=gap$start,
                  max=gap$end-chunkSize
                )
                
              } else {
                
                chunkStart = gap$start
              }
              
              chunkEnd = chunkStart + chunkSize
              
              chunkIndex = which(
                chrCenters >= chunkStart &
                  chrCenters <= chunkEnd
              )
              
              sampleCenters[j] = (chunkStart + chunkEnd)/2
              
              sampledBrDU[j] = SummariseValues(
                Noise.IP.df[chrIndex[chunkIndex],5]
              )
              
              sampledInput[j] = SummariseValues(
                Noise.In.df[chrIndex[chunkIndex],5]
              )
            }
            
            nonPeakIndex = rep(FALSE, length(chrIndex))
            
            for(j in 1:nrow(gaps)){
              
              nonPeakIndex = nonPeakIndex |
                (chrCenters >= gaps$start[j] &
                   chrCenters <= gaps$end[j])
            }
            
            BrDU_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledBrDU,
              fallbackValues=Noise.IP.df[
                chrIndex[nonPeakIndex],5
              ]
            )
            
            Input_noise[chrIndex] = PredictNoise(
              binCenters=chrCenters,
              sampleCenters=sampleCenters,
              sampleValues=sampledInput,
              fallbackValues=Noise.In.df[
                chrIndex[nonPeakIndex],5
              ]
            )
          }
          
          list(BrDU_noise=BrDU_noise,
               Input_noise=Input_noise)
        }
        
        noise = EstimateNoise()
        
        # whole-genome library-size normalization
        IP_Sum = sum(Noise.IP.df[,5], na.rm=TRUE)
        In_Sum = sum(Noise.In.df[,5], na.rm=TRUE)
        
        if(!is.finite(In_Sum) || In_Sum <= 0){
          stop("Whole-genome Input coverage sum is zero.")
        }
        
        corrFactor = IP_Sum/In_Sum
        
        # convert whole-genome background profiles into transferable levels
        positive.IP.noise = noise$BrDU_noise[
          is.finite(noise$BrDU_noise) & noise$BrDU_noise > NoiseFloor
        ]
        
        positive.In.noise = noise$Input_noise[
          is.finite(noise$Input_noise) & noise$Input_noise > NoiseFloor
        ]
        
        IP.noise = median(positive.IP.noise, na.rm=TRUE)
        In.noise = median(positive.In.noise, na.rm=TRUE)*corrFactor
        
        if(!is.finite(IP.noise) || IP.noise < NoiseFloor){
          IP.noise = NoiseFloor
        }
        
        if(!is.finite(In.noise) || In.noise < NoiseFloor){
          In.noise = NoiseFloor
        }
        
        # normalize rDNA Input using the whole-genome correction factor
        In.score.norm = In.df[,5]*corrFactor
        
        # calculate rDNA enrichment ratios
        Ratio.ipin = IP.df[,5]/In.score.norm
        
        Ratio.ipnoise = IP.df[,5]/IP.noise
        
        Ratio.ipin.noise = Ratio.ipnoise/
          (In.score.norm/In.noise)
        
        Ratio.ipin[!is.finite(Ratio.ipin)] = 0
        Ratio.ipnoise[!is.finite(Ratio.ipnoise)] = 0
        Ratio.ipin.noise[!is.finite(Ratio.ipin.noise)] = 0
        
        # output contains rDNA coordinates only
        ratio.df = data.frame(
          chrom=IP.df[,1],
          chromStart=IP.df[,2],
          chromEnd=IP.df[,3],
          name=RatioName,
          ip.score=IP.df[,5],
          in.score=round(In.score.norm, 4),
          ip.noise=round(rep(IP.noise, nrow(IP.df)), 4),
          in.noise=round(rep(In.noise, nrow(IP.df)), 4),
          ratio.ipin=round(Ratio.ipin, 4),
          ratio.ipnoise=round(Ratio.ipnoise, 4),
          ratio.ipin.noise=round(Ratio.ipin.noise, 4)
        )
        
        write.table(
          ratio.df,
          file=paste0(dir, Pro_1, "/", "Ratios_ma_rdna", "/",
                      RatioName, ".bed"),
          quote=FALSE, row.names=FALSE, sep="\t"
        )
      }
      
      # Watson-strand rDNA ratio
      CalculateRatio(
        IP_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_BrDU_rDNA_watson.bed"
        ),
        Input_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_Input_rDNA_watson.bed"
        ),
        Noise_IP_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_BrDU_ma_watson.bed"
        ),
        Noise_Input_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_Input_ma_watson.bed"
        ),
        RatioName=paste0(Pro_1, "_BrDU_rDNA_watson")
      )
      
      # Crick-strand rDNA ratio
      CalculateRatio(
        IP_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_BrDU_rDNA_crick.bed"
        ),
        Input_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_Input_rDNA_crick.bed"
        ),
        Noise_IP_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_BrDU_ma_crick.bed"
        ),
        Noise_Input_coverage=paste0(
          dir, Pro_1, "/Coverage_ma_rdna/",
          Pro_1, "_Input_ma_crick.bed"
        ),
        RatioName=paste0(Pro_1, "_BrDU_rDNA_crick")
      )
      
      # strand-collapsed rDNA ratio
      CalculateRatio(
        IP_coverage=c(
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_BrDU_rDNA_watson.bed"
          ),
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_BrDU_rDNA_crick.bed"
          )
        ),
        Input_coverage=c(
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_Input_rDNA_watson.bed"
          ),
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_Input_rDNA_crick.bed"
          )
        ),
        Noise_IP_coverage=c(
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_BrDU_ma_watson.bed"
          ),
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_BrDU_ma_crick.bed"
          )
        ),
        Noise_Input_coverage=c(
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_Input_ma_watson.bed"
          ),
          paste0(
            dir, Pro_1, "/Coverage_ma_rdna/",
            Pro_1, "_Input_ma_crick.bed"
          )
        ),
        RatioName=paste0(Pro_1, "_BrDU_rDNA_collapsed")
      )
    }
    
  }
  
  ##
  
  rm(list=ls())
  gc()
  
  #
  message("✅ Alignment & Primary Analysis complete!")
  #
  
}

##

ChIPseq_Analysis_early_late <- function(ExpTitle = "None",
                                        Window = "None",
                                        Alignment = "generic",
                                        StrandMode = "collapsed",
                                        Directory = "None"){
  
  ## load packages
  
  packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
                "VennDiagram", "grid", "gridBase", "gridExtra", "BSgenome.Scerevisiae.UCSC.sacCer3")
  
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  ## load and process basic files
  
  E_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/E_Rep.bed", header = FALSE, quote = "/t")[ ,1:4]
  L_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/L_Rep.bed", header = FALSE, quote = "/t")[ ,1:4]
  
  colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName")
  colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName")
  
  E_Ori$oriCenter <- round((E_Ori$chromStart + E_Ori$chromEnd)/2)
  L_Ori$oriCenter <- round((L_Ori$chromStart + L_Ori$chromEnd)/2)
  
  E_Ori$peakSummit <- E_Ori$oriCenter
  L_Ori$peakSummit <- L_Ori$oriCenter
  
  Alignment <- match.arg(Alignment, c("generic", "malign"))
  StrandMode <- match.arg(StrandMode, c("collapsed", "separated"))
  
  useDef <- function(a,d) if(length(a) == 0 || is.null(a) || is.na(a) || a == "") d else a
  
  ExpTitle = useDef(ExpTitle, "None")
  Directory = useDef(Directory, "None")
  
  if(ExpTitle == "None"){
    stop("ExpTitle must be provided.")
  }
  
  if(dir.exists(ExpTitle)){
    SamplePath = ExpTitle
    Pro_1 = basename(ExpTitle)
  } else {
    Pro_1 = ExpTitle
    
    if(Directory == "None"){
      SamplePath = paste0("~/Desktop/", Pro_1)
    } else {
      SamplePath = paste0(Directory, "/", Pro_1)
    }
  }
  
  if(!dir.exists(SamplePath)){
    stop("Experiment folder is missing: ", SamplePath)
  }
  
  RatioFolder = if(Alignment == "generic") "Ratios" else "Ratios_ma"
  RatioPath = paste0(SamplePath, "/", RatioFolder)
  
  if(!dir.exists(RatioPath)){
    stop("Ratio folder is missing: ", RatioPath)
  }
  
  if(Window == "None"){
    InfoCandidates = c(
      paste0(SamplePath, "/", Pro_1, "_Info.txt"),
      paste0(SamplePath, "/Info/", Pro_1, "_Info.txt"),
      paste0(SamplePath, "/Coverage/", Pro_1, "_Info.txt")
    )
    
    InfoFile = InfoCandidates[file.exists(InfoCandidates)][1]
    
    if(is.na(InfoFile)){
      stop("Window could not be resolved. Provide Window explicitly.")
    }
    
    Window = as.numeric(read.table(InfoFile, header=FALSE)[5,2])
  } else {
    Window = as.numeric(Window)
  }
  
  if(!is.finite(Window) || Window <= 0){
    stop("Window must be a positive number.")
  }
  
  metrics = c("ip.score", "ratio.ipin", "ratio.ipnoise", "ratio.ipin.noise")
  metricTitles = c("ChIP", "ChIP / Input", "ChIP filtered", "ChIP / Input filtered")
  metricYlabels = c("Coverage", "Enrichment over input", "Enrichment over noise", "Clean enrichment")
  
  ReadRatio <- function(strand){
    
    RatioFile = paste0(RatioPath, "/", Pro_1, "_ChIP_", strand, ".bed")
    
    if(!file.exists(RatioFile)){
      stop("Ratio file is missing: ", RatioFile)
    }
    
    Ratio = read.table(RatioFile, header=TRUE, check.names=FALSE)
    
    MissingColumns = metrics[!metrics %in% colnames(Ratio)]
    
    if(length(MissingColumns) > 0){
      stop("Required ratio columns are missing from ", basename(RatioFile), ": ",
           paste(MissingColumns, collapse=", "))
    }
    
    Ratio$chromStart = as.numeric(Ratio$chromStart)
    Ratio$chromEnd = as.numeric(Ratio$chromEnd)
    
    for(V in metrics){
      Ratio[,V] = as.numeric(Ratio[,V])
    }
    
    Ratio
  }
  
  if(StrandMode == "collapsed"){
    RatioData = list(collapsed=ReadRatio("collapsed"))
  } else {
    RatioData = list(watson=ReadRatio("watson"), crick=ReadRatio("crick"))
  }
  
  RatioByChrom = lapply(
    RatioData,
    function(Ratio) split(Ratio, as.character(Ratio$chrom))
  )
  
  infer_step <- function(RatioData){
    
    diffs = unlist(
      lapply(RatioData, function(Ratio){
        StartsByChrom = split(Ratio$chromStart, as.character(Ratio$chrom))
        unlist(
          lapply(StartsByChrom, function(starts){
            starts = sort(unique(starts[is.finite(starts)]))
            diff(starts)
          }),
          use.names=FALSE
        )
      }),
      use.names=FALSE
    )
    
    diffs = diffs[is.finite(diffs) & diffs > 0]
    
    if(length(diffs) == 0){
      stop("Could not infer coverage step size.")
    }
    
    DiffTable = table(diffs)
    as.numeric(names(DiffTable)[which.max(DiffTable)])
  }
  
  step = infer_step(RatioData)
  expected_n = max(1L, round(2*Window/step))
  
  extract_centered_values <- function(Ratio, center, V){
    
    regionStart = center - Window
    regionEnd = center + Window
    
    if(is.null(Ratio)){
      return(rep(0, expected_n))
    }
    
    Signal = Ratio[
      Ratio$chromStart >= regionStart & Ratio$chromStart <= regionEnd, ]
    
    Signal = Signal[order(Signal$chromStart), ]
    
    if(nrow(Signal) == 0){
      return(rep(0, expected_n))
    }
    
    Values = Signal[,V]
    Values[!is.finite(Values)] = 0
    
    if(length(Values) < expected_n){
      missing_n = expected_n - length(Values)
      left_missing = round((min(Signal$chromStart, na.rm=TRUE) - regionStart)/step)
      left_missing = min(max(left_missing, 0L), missing_n)
      right_missing = missing_n - left_missing
      Values = c(rep(0, left_missing), Values, rep(0, right_missing))
    }
    
    if(length(Values) > expected_n){
      Values = Values[seq_len(expected_n)]
    }
    
    Values
  }
  
  CalculateProfiles <- function(OriginList, originClass){
    
    if(nrow(OriginList) == 0){
      return(NULL)
    }
    
    message(paste0("Calculating ", StrandMode, " ChIP profiles over ",
                   tolower(originClass), " set"))
    
    Results = list()
    
    for(strand in names(RatioData)){
      
      Results[[strand]] = list()
      
      for(V in metrics){
        
        ProfileMatrix = vapply(
          seq_len(nrow(OriginList)),
          function(i){
            extract_centered_values(
              Ratio=RatioByChrom[[strand]][[as.character(OriginList$chrom[i])]],
              center=OriginList$peakSummit[i], V=V
            )
          },
          numeric(expected_n)
        )
        
        if(is.null(dim(ProfileMatrix))){
          ProfileMatrix = matrix(ProfileMatrix, ncol=1)
        }
        
        colnames(ProfileMatrix) = OriginList$oriName
        
        Results[[strand]][[V]] = list(
          median=apply(ProfileMatrix, 1, median, na.rm=TRUE),
          heatmap=ProfileMatrix
        )
      }
    }
    
    Results
  }
  
  SampleRes = list(
    EarlyOrigin=CalculateProfiles(E_Ori, "early origins"),
    LateOrigin=CalculateProfiles(L_Ori, "late origins")
  )
  
  SmoothProfile <- function(Profile){
    
    Profile[!is.finite(Profile)] = 0
    
    if(length(Profile) < 4 || length(unique(Profile)) < 2){
      return(Profile)
    }
    
    smooth.spline(seq_along(Profile), Profile, spar=0.5)$y
  }
  
  PlotProfile <- function(Result, V, PlotHeader, Y_label){
    
    if(is.null(Result)){
      plot.new()
      text(0.5, 0.5, "No origins in cohort", cex=1.2, font=2)
      return(invisible(NULL))
    }
    
    if(StrandMode == "collapsed"){
      
      Collapsed = SmoothProfile(Result$collapsed[[V]]$median)
      Y_min = min(Collapsed, na.rm=TRUE)
      Y_max = max(Collapsed, na.rm=TRUE)
      Y_pad = 0.08*(Y_max - Y_min)
      if(!is.finite(Y_pad) || Y_pad == 0) Y_pad = 0.5
      
      plot(Collapsed, ylim=c(Y_min-Y_pad, Y_max+Y_pad), main=PlotHeader,
           ylab=Y_label, cex.main=1, xlab="Distance from OriCenter (Kbp)",
           xaxt="n", col="darkorchid4", type="l", lwd=2,
           bty="n", las=2, xaxs="i", yaxs="i")
      
    } else {
      
      Watson = SmoothProfile(Result$watson[[V]]$median)
      Crick = -SmoothProfile(Result$crick[[V]]$median)
      Y_max = max(c(abs(Watson), abs(Crick)), na.rm=TRUE)
      if(!is.finite(Y_max) || Y_max == 0) Y_max = 0.5
      
      plot(Watson, ylim=c(-Y_max, Y_max), main=PlotHeader,
           ylab=Y_label, cex.main=1, xlab="Distance from OriCenter (Kbp)",
           xaxt="n", col="brown3", type="l", lwd=2,
           bty="n", las=2, xaxs="i", yaxs="i")
      lines(Crick, col="cornflowerblue", lwd=2)
      abline(h=0, col="grey50", lty=2)
      legend("topright", legend=c("Watson", "Crick"),
             col=c("brown3", "cornflowerblue"), lwd=2, bty="n", cex=0.8)
    }
    
    axisLabels = seq(-Window, +Window, length.out=9)
    axisLabels[c(2,4,6,8)] = NA
    At = (Window/step)*seq(0, 2, 0.25)
    At[1] = 1
    axis(1, at=At, labels=signif(axisLabels/1000, 2))
  }
  
  PlotPairwise <- function(originClass_1, originClass_2, V,
                           PlotHeader, Y_label){
    
    Result_1 = SampleRes[[originClass_1]]
    Result_2 = SampleRes[[originClass_2]]
    
    if(is.null(Result_1) || is.null(Result_2)){
      plot.new()
      text(0.5, 0.5, "Missing origin cohort", cex=1.2, font=2)
      return(invisible(NULL))
    }
    
    if(StrandMode == "collapsed"){
      
      P1 = SmoothProfile(Result_1$collapsed[[V]]$median)
      P2 = SmoothProfile(Result_2$collapsed[[V]]$median)
      Y_min = min(c(P1, P2), na.rm=TRUE)
      Y_max = max(c(P1, P2), na.rm=TRUE)
      Y_pad = 0.08*(Y_max - Y_min)
      if(!is.finite(Y_pad) || Y_pad == 0) Y_pad = 0.5
      
      plot(P1, ylim=c(Y_min-Y_pad, Y_max+Y_pad), main=PlotHeader,
           ylab=Y_label, cex.main=1, xlab="Distance from OriCenter (Kbp)",
           xaxt="n", col=adjustcolor("darkorchid4", alpha.f=0.9),
           type="l", lwd=2, bty="n", las=2, xaxs="i", yaxs="i")
      lines(P2, col=adjustcolor("darkorange3", alpha.f=0.9), lwd=2)
      legend("topright", legend=c("E", "L"),
             col=c("darkorchid4", "darkorange3"), lwd=2, bty="n", cex=0.8)
      
    } else {
      
      W1 = SmoothProfile(Result_1$watson[[V]]$median)
      C1 = -SmoothProfile(Result_1$crick[[V]]$median)
      W2 = SmoothProfile(Result_2$watson[[V]]$median)
      C2 = -SmoothProfile(Result_2$crick[[V]]$median)
      Y_max = max(abs(c(W1, C1, W2, C2)), na.rm=TRUE)
      if(!is.finite(Y_max) || Y_max == 0) Y_max = 0.5
      
      plot(W1, ylim=c(-Y_max, Y_max), main=PlotHeader,
           ylab=Y_label, cex.main=1, xlab="Distance from OriCenter (Kbp)",
           xaxt="n", col="brown3", type="l", lwd=2,
           bty="n", las=2, xaxs="i", yaxs="i")
      lines(C1, col="cornflowerblue", lwd=2)
      lines(W2, col="brown3", lwd=2, lty=2)
      lines(C2, col="cornflowerblue", lwd=2, lty=2)
      abline(h=0, col="grey50", lty=3)
      legend("topright", legend=c("E Watson", "E Crick", "L Watson", "L Crick"),
             col=c("brown3", "cornflowerblue", "brown3", "cornflowerblue"),
             lty=c(1,1,2,2), lwd=2, bty="n", cex=0.7)
    }
    
    axisLabels = seq(-Window, +Window, length.out=9)
    axisLabels[c(2,4,6,8)] = NA
    At = (Window/step)*seq(0, 2, 0.25)
    At[1] = 1
    axis(1, at=At, labels=signif(axisLabels/1000, 2))
  }
  
  OutputFile = paste0(SamplePath, "/", Pro_1, "_ChIP_early_late_",
                      Alignment, "_", StrandMode, ".pdf")
  
  pdf(OutputFile, width=12, height=10)
  
  ## first page
  
  OriginLists = list(EarlyOrigin=E_Ori, LateOrigin=L_Ori)
  
  pretty_origin_class <- function(x){
    switch(x, EarlyOrigin="Early origins", LateOrigin="Late origins", x)
  }
  
  plot_title_panel <- function(){
    
    plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
         xlab="", ylab="", bty="n")
    
    text(0.5, 0.62, labels=paste0("Experiment: ", Pro_1),
         cex=2.2, font=2, family="serif")
    
    text(0.5, 0.30,
         labels="Early and late origin ChIP signals centered on origin midpoint",
         cex=1.45, font=3, family="serif")
    
    text(0.5, 0.08, labels=date(), cex=1.4, font=3, family="serif")
  }
  
  plot_origin_title <- function(originClass){
    
    plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
         xlab="", ylab="", bty="n")
    
    text(0.5, 0.35, labels=pretty_origin_class(originClass),
         cex=1.35, font=2, family="serif")
  }
  
  plot_origin_venn <- function(OriginList, originClass){
    
    if(is.null(OriginList) || nrow(OriginList) == 0){
      plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
           xlab="", ylab="", bty="n")
      text(0.5, 0.5, labels=paste0("No origins\n", pretty_origin_class(originClass)),
           cex=1.4, font=2, family="serif")
      return(invisible(NULL))
    }
    
    invisible(futile.logger::flog.threshold(
      futile.logger::ERROR, name="VennDiagramLogger"
    ))
    
    VennObject = venn.diagram(
      list(Origins=OriginList$oriCenter), filename=NULL, disable.logging=TRUE,
      category.names="", cat.cex=0, lwd=2, lty="blank",
      fill=adjustcolor("cornflowerblue", alpha.f=0.45),
      cex=1.6, fontface="italic", fontfamily="serif"
    )
    
    plot.new()
    vps = baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    grid.draw(VennObject)
    popViewport(3)
  }
  
  par(oma=c(2,1,1,1))
  par(mar=c(0,0,0,0))
  
  PlotMat = matrix(
    c(
      0,0,0,0,0,0,0,0,
      1,1,1,1,1,1,1,1,
      0,2,2,2,3,3,3,0,
      0,4,4,4,5,5,5,0,
      0,0,0,0,0,0,0,0
    ),
    nrow=5, ncol=8, byrow=TRUE
  )
  
  layout(PlotMat, widths=rep(1,8),
         heights=c(0.2,0.7,0.3,1.15,0.25), respect=FALSE)
  
  plot_title_panel()
  plot_origin_title("EarlyOrigin")
  plot_origin_title("LateOrigin")
  
  par(mar=c(1,1,2,1), pty="s")
  plot_origin_venn(OriginLists$EarlyOrigin, "EarlyOrigin")
  plot_origin_venn(OriginLists$LateOrigin, "LateOrigin")
  
  mtext("Page 1", side=1, outer=TRUE, line=0, font=3, cex=1)
  
  ## second page
  
  layout(matrix(1:8, nrow=2, ncol=4, byrow=TRUE))
  par(oma=c(3,1,3,1))
  par(mar=c(4,3,4,3)+0.1)
  
  for(originClass in c("EarlyOrigin", "LateOrigin")){
    originLabel = if(originClass == "EarlyOrigin") "Early" else "Late"
    
    for(i in seq_along(metrics)){
      PlotProfile(SampleRes[[originClass]], metrics[i],
                  paste0(originLabel, ": ", metricTitles[i]), metricYlabels[i])
    }
  }
  
  mtext("Read enrichments at early and late origins", side=3, line=0,
        outer=TRUE, font=2, cex=2)
  mtext("Page 2", side=1, line=0, outer=TRUE, font=3, cex=1)
  
  ## third page
  
  layout(matrix(1:12, nrow=3, ncol=4, byrow=TRUE))
  par(oma=c(3,1,3,1))
  par(mar=c(4,4,4,2)+0.1)
  
  for(i in seq_along(metrics)){
    PlotPairwise("EarlyOrigin", "LateOrigin", metrics[i],
                 metricTitles[i], metricYlabels[i])
  }
  
  for(i in 1:8){
    plot.new()
  }
  
  mtext("Pairwise comparative enrichment at early and late origins",
        side=3, line=0, outer=TRUE, font=2, cex=2)
  mtext("Page 3", side=1, line=0, outer=TRUE, font=3, cex=1)
  
  dev.off()
  
  message(paste0("Plot saved as ", basename(OutputFile)))
}

BrDUseq_Analysis_early_late <- function(ExpTitle = "None",
                                        Window = "None",
                                        Alignment = "generic",
                                        StrandMode = "collapsed",
                                        Directory = "None"){
  
  ## load packages
  
  packages <- c("basicPlotteR", "plyr", "tidyverse", "dplyr", "plotrix", "rasterpdf", "imager",
                "VennDiagram", "grid", "gridBase", "gridExtra", "BSgenome.Scerevisiae.UCSC.sacCer3")
  
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  ## load and process basic files
  
  E_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/E_Rep.bed", header = FALSE, quote = "/t")[ ,1:4]
  L_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/L_Rep.bed", header = FALSE, quote = "/t")[ ,1:4]
  
  colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName")
  colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName")
  
  E_Ori$oriCenter <- round((E_Ori$chromStart + E_Ori$chromEnd)/2)
  L_Ori$oriCenter <- round((L_Ori$chromStart + L_Ori$chromEnd)/2)
  
  E_Ori$peakSummit <- E_Ori$oriCenter
  L_Ori$peakSummit <- L_Ori$oriCenter
  
  Alignment <- match.arg(Alignment, c("generic", "malign"))
  StrandMode <- match.arg(StrandMode, c("collapsed", "separated"))
  
  useDef <- function(a,d) if(length(a) == 0 || is.null(a) || is.na(a) || a == "") d else a
  
  ExpTitle = useDef(ExpTitle, "None")
  Directory = useDef(Directory, "None")
  
  if(ExpTitle == "None"){
    stop("ExpTitle must be provided.")
  }
  
  if(dir.exists(ExpTitle)){
    SamplePath = ExpTitle
    Pro_1 = basename(ExpTitle)
  } else {
    Pro_1 = ExpTitle
    
    if(Directory == "None"){
      SamplePath = paste0("~/Desktop/", Pro_1)
    } else {
      SamplePath = paste0(Directory, "/", Pro_1)
    }
  }
  
  if(!dir.exists(SamplePath)){
    stop("Experiment folder is missing: ", SamplePath)
  }
  
  RatioFolder = if(Alignment == "generic") "Ratios" else "Ratios_ma"
  RatioPath = paste0(SamplePath, "/", RatioFolder)
  
  if(!dir.exists(RatioPath)){
    stop("Ratio folder is missing: ", RatioPath)
  }
  
  if(Window == "None"){
    InfoCandidates = c(
      paste0(SamplePath, "/", Pro_1, "_Info.txt"),
      paste0(SamplePath, "/Info/", Pro_1, "_Info.txt"),
      paste0(SamplePath, "/Coverage/", Pro_1, "_Info.txt")
    )
    
    InfoFile = InfoCandidates[file.exists(InfoCandidates)][1]
    
    if(is.na(InfoFile)){
      stop("Window could not be resolved. Provide Window explicitly.")
    }
    
    Window = as.numeric(read.table(InfoFile, header=FALSE)[5,2])
  } else {
    Window = as.numeric(Window)
  }
  
  if(!is.finite(Window) || Window <= 0){
    stop("Window must be a positive number.")
  }
  
  metrics = c("ip.score", "ratio.ipin", "ratio.ipnoise", "ratio.ipin.noise")
  metricTitles = c("BrDU", "BrDU / Input", "BrDU filtered", "BrDU / Input filtered")
  metricYlabels = c("Coverage", "Enrichment over input", "Enrichment over noise", "Clean enrichment")
  
  ReadRatio <- function(strand){
    
    RatioFile = paste0(RatioPath, "/", Pro_1, "_BrDU_", strand, ".bed")
    
    if(!file.exists(RatioFile)){
      stop("Ratio file is missing: ", RatioFile)
    }
    
    Ratio = read.table(RatioFile, header=TRUE, check.names=FALSE)
    
    MissingColumns = metrics[!metrics %in% colnames(Ratio)]
    
    if(length(MissingColumns) > 0){
      stop("Required ratio columns are missing from ", basename(RatioFile), ": ",
           paste(MissingColumns, collapse=", "))
    }
    
    Ratio$chromStart = as.numeric(Ratio$chromStart)
    Ratio$chromEnd = as.numeric(Ratio$chromEnd)
    
    for(V in metrics){
      Ratio[,V] = as.numeric(Ratio[,V])
    }
    
    Ratio
  }
  
  if(StrandMode == "collapsed"){
    RatioData = list(collapsed=ReadRatio("collapsed"))
  } else {
    RatioData = list(watson=ReadRatio("watson"), crick=ReadRatio("crick"))
  }
  
  RatioByChrom = lapply(
    RatioData,
    function(Ratio) split(Ratio, as.character(Ratio$chrom))
  )
  
  infer_step <- function(RatioData){
    
    diffs = unlist(
      lapply(RatioData, function(Ratio){
        StartsByChrom = split(Ratio$chromStart, as.character(Ratio$chrom))
        unlist(
          lapply(StartsByChrom, function(starts){
            starts = sort(unique(starts[is.finite(starts)]))
            diff(starts)
          }),
          use.names=FALSE
        )
      }),
      use.names=FALSE
    )
    
    diffs = diffs[is.finite(diffs) & diffs > 0]
    
    if(length(diffs) == 0){
      stop("Could not infer coverage step size.")
    }
    
    DiffTable = table(diffs)
    as.numeric(names(DiffTable)[which.max(DiffTable)])
  }
  
  step = infer_step(RatioData)
  expected_n = max(1L, round(2*Window/step))
  
  extract_centered_values <- function(Ratio, center, V){
    
    regionStart = center - Window
    regionEnd = center + Window
    
    if(is.null(Ratio)){
      return(rep(0, expected_n))
    }
    
    Signal = Ratio[
      Ratio$chromStart >= regionStart & Ratio$chromStart <= regionEnd, ]
    
    Signal = Signal[order(Signal$chromStart), ]
    
    if(nrow(Signal) == 0){
      return(rep(0, expected_n))
    }
    
    Values = Signal[,V]
    Values[!is.finite(Values)] = 0
    
    if(length(Values) < expected_n){
      missing_n = expected_n - length(Values)
      left_missing = round((min(Signal$chromStart, na.rm=TRUE) - regionStart)/step)
      left_missing = min(max(left_missing, 0L), missing_n)
      right_missing = missing_n - left_missing
      Values = c(rep(0, left_missing), Values, rep(0, right_missing))
    }
    
    if(length(Values) > expected_n){
      Values = Values[seq_len(expected_n)]
    }
    
    Values
  }
  
  CalculateProfiles <- function(OriginList, originClass){
    
    if(nrow(OriginList) == 0){
      return(NULL)
    }
    
    message(paste0("Calculating ", StrandMode, " BrDU profiles over ",
                   tolower(originClass), " set"))
    
    Results = list()
    
    for(strand in names(RatioData)){
      
      Results[[strand]] = list()
      
      for(V in metrics){
        
        ProfileMatrix = vapply(
          seq_len(nrow(OriginList)),
          function(i){
            extract_centered_values(
              Ratio=RatioByChrom[[strand]][[as.character(OriginList$chrom[i])]],
              center=OriginList$peakSummit[i], V=V
            )
          },
          numeric(expected_n)
        )
        
        if(is.null(dim(ProfileMatrix))){
          ProfileMatrix = matrix(ProfileMatrix, ncol=1)
        }
        
        colnames(ProfileMatrix) = OriginList$oriName
        
        Results[[strand]][[V]] = list(
          median=apply(ProfileMatrix, 1, median, na.rm=TRUE),
          heatmap=ProfileMatrix
        )
      }
    }
    
    Results
  }
  
  SampleRes = list(
    EarlyOrigin=CalculateProfiles(E_Ori, "early origins"),
    LateOrigin=CalculateProfiles(L_Ori, "late origins")
  )
  
  SmoothProfile <- function(Profile){
    
    Profile[!is.finite(Profile)] = 0
    
    if(length(Profile) < 4 || length(unique(Profile)) < 2){
      return(Profile)
    }
    
    smooth.spline(seq_along(Profile), Profile, spar=0.5)$y
  }
  
  PlotProfile <- function(Result, V, PlotHeader, Y_label){
    
    if(is.null(Result)){
      plot.new()
      text(0.5, 0.5, "No origins in cohort", cex=1.2, font=2)
      return(invisible(NULL))
    }
    
    if(StrandMode == "collapsed"){
      
      Collapsed = SmoothProfile(Result$collapsed[[V]]$median)
      Y_min = min(Collapsed, na.rm=TRUE)
      Y_max = max(Collapsed, na.rm=TRUE)
      Y_pad = 0.08*(Y_max - Y_min)
      if(!is.finite(Y_pad) || Y_pad == 0) Y_pad = 0.5
      
      plot(Collapsed, ylim=c(Y_min-Y_pad, Y_max+Y_pad), main=PlotHeader,
           ylab=Y_label, cex.main=1, xlab="Distance from OriCenter (Kbp)",
           xaxt="n", col="darkorchid4", type="l", lwd=2,
           bty="n", las=2, xaxs="i", yaxs="i")
      
    } else {
      
      Watson = SmoothProfile(Result$watson[[V]]$median)
      Crick = -SmoothProfile(Result$crick[[V]]$median)
      Y_max = max(c(abs(Watson), abs(Crick)), na.rm=TRUE)
      if(!is.finite(Y_max) || Y_max == 0) Y_max = 0.5
      
      plot(Watson, ylim=c(-Y_max, Y_max), main=PlotHeader,
           ylab=Y_label, cex.main=1, xlab="Distance from OriCenter (Kbp)",
           xaxt="n", col="brown3", type="l", lwd=2,
           bty="n", las=2, xaxs="i", yaxs="i")
      lines(Crick, col="cornflowerblue", lwd=2)
      abline(h=0, col="grey50", lty=2)
      legend("topright", legend=c("Watson", "Crick"),
             col=c("brown3", "cornflowerblue"), lwd=2, bty="n", cex=0.8)
    }
    
    axisLabels = seq(-Window, +Window, length.out=9)
    axisLabels[c(2,4,6,8)] = NA
    At = (Window/step)*seq(0, 2, 0.25)
    At[1] = 1
    axis(1, at=At, labels=signif(axisLabels/1000, 2))
  }
  
  PlotPairwise <- function(originClass_1, originClass_2, V,
                           PlotHeader, Y_label){
    
    Result_1 = SampleRes[[originClass_1]]
    Result_2 = SampleRes[[originClass_2]]
    
    if(is.null(Result_1) || is.null(Result_2)){
      plot.new()
      text(0.5, 0.5, "Missing origin cohort", cex=1.2, font=2)
      return(invisible(NULL))
    }
    
    if(StrandMode == "collapsed"){
      
      P1 = SmoothProfile(Result_1$collapsed[[V]]$median)
      P2 = SmoothProfile(Result_2$collapsed[[V]]$median)
      Y_min = min(c(P1, P2), na.rm=TRUE)
      Y_max = max(c(P1, P2), na.rm=TRUE)
      Y_pad = 0.08*(Y_max - Y_min)
      if(!is.finite(Y_pad) || Y_pad == 0) Y_pad = 0.5
      
      plot(P1, ylim=c(Y_min-Y_pad, Y_max+Y_pad), main=PlotHeader,
           ylab=Y_label, cex.main=1, xlab="Distance from OriCenter (Kbp)",
           xaxt="n", col=adjustcolor("darkorchid4", alpha.f=0.9),
           type="l", lwd=2, bty="n", las=2, xaxs="i", yaxs="i")
      lines(P2, col=adjustcolor("darkorange3", alpha.f=0.9), lwd=2)
      legend("topright", legend=c("E", "L"),
             col=c("darkorchid4", "darkorange3"), lwd=2, bty="n", cex=0.8)
      
    } else {
      
      W1 = SmoothProfile(Result_1$watson[[V]]$median)
      C1 = -SmoothProfile(Result_1$crick[[V]]$median)
      W2 = SmoothProfile(Result_2$watson[[V]]$median)
      C2 = -SmoothProfile(Result_2$crick[[V]]$median)
      Y_max = max(abs(c(W1, C1, W2, C2)), na.rm=TRUE)
      if(!is.finite(Y_max) || Y_max == 0) Y_max = 0.5
      
      plot(W1, ylim=c(-Y_max, Y_max), main=PlotHeader,
           ylab=Y_label, cex.main=1, xlab="Distance from OriCenter (Kbp)",
           xaxt="n", col="brown3", type="l", lwd=2,
           bty="n", las=2, xaxs="i", yaxs="i")
      lines(C1, col="cornflowerblue", lwd=2)
      lines(W2, col="brown3", lwd=2, lty=2)
      lines(C2, col="cornflowerblue", lwd=2, lty=2)
      abline(h=0, col="grey50", lty=3)
      legend("topright", legend=c("E Watson", "E Crick", "L Watson", "L Crick"),
             col=c("brown3", "cornflowerblue", "brown3", "cornflowerblue"),
             lty=c(1,1,2,2), lwd=2, bty="n", cex=0.7)
    }
    
    axisLabels = seq(-Window, +Window, length.out=9)
    axisLabels[c(2,4,6,8)] = NA
    At = (Window/step)*seq(0, 2, 0.25)
    At[1] = 1
    axis(1, at=At, labels=signif(axisLabels/1000, 2))
  }
  
  OutputFile = paste0(SamplePath, "/", Pro_1, "_BrDU_early_late_",
                      Alignment, "_", StrandMode, ".pdf")
  
  pdf(OutputFile, width=12, height=10)
  
  ## first page
  
  OriginLists = list(EarlyOrigin=E_Ori, LateOrigin=L_Ori)
  
  pretty_origin_class <- function(x){
    switch(x, EarlyOrigin="Early origins", LateOrigin="Late origins", x)
  }
  
  plot_title_panel <- function(){
    
    plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
         xlab="", ylab="", bty="n")
    
    text(0.5, 0.62, labels=paste0("Experiment: ", Pro_1),
         cex=2.2, font=2, family="serif")
    
    text(0.5, 0.30,
         labels="Early and late origin BrDU signals centered on origin midpoint",
         cex=1.45, font=3, family="serif")
    
    text(0.5, 0.08, labels=date(), cex=1.4, font=3, family="serif")
  }
  
  plot_origin_title <- function(originClass){
    
    plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
         xlab="", ylab="", bty="n")
    
    text(0.5, 0.35, labels=pretty_origin_class(originClass),
         cex=1.35, font=2, family="serif")
  }
  
  plot_origin_venn <- function(OriginList, originClass){
    
    if(is.null(OriginList) || nrow(OriginList) == 0){
      plot(NULL, xlim=c(0,1), ylim=c(0,1), axes=FALSE,
           xlab="", ylab="", bty="n")
      text(0.5, 0.5, labels=paste0("No origins\n", pretty_origin_class(originClass)),
           cex=1.4, font=2, family="serif")
      return(invisible(NULL))
    }
    
    invisible(futile.logger::flog.threshold(
      futile.logger::ERROR, name="VennDiagramLogger"
    ))
    
    VennObject = venn.diagram(
      list(Origins=OriginList$oriCenter), filename=NULL, disable.logging=TRUE,
      category.names="", cat.cex=0, lwd=2, lty="blank",
      fill=adjustcolor("cornflowerblue", alpha.f=0.45),
      cex=1.6, fontface="italic", fontfamily="serif"
    )
    
    plot.new()
    vps = baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    grid.draw(VennObject)
    popViewport(3)
  }
  
  par(oma=c(2,1,1,1))
  par(mar=c(0,0,0,0))
  
  PlotMat = matrix(
    c(
      0,0,0,0,0,0,0,0,
      1,1,1,1,1,1,1,1,
      0,2,2,2,3,3,3,0,
      0,4,4,4,5,5,5,0,
      0,0,0,0,0,0,0,0
    ),
    nrow=5, ncol=8, byrow=TRUE
  )
  
  layout(PlotMat, widths=rep(1,8),
         heights=c(0.2,0.7,0.3,1.15,0.25), respect=FALSE)
  
  plot_title_panel()
  plot_origin_title("EarlyOrigin")
  plot_origin_title("LateOrigin")
  
  par(mar=c(1,1,2,1), pty="s")
  plot_origin_venn(OriginLists$EarlyOrigin, "EarlyOrigin")
  plot_origin_venn(OriginLists$LateOrigin, "LateOrigin")
  
  mtext("Page 1", side=1, outer=TRUE, line=0, font=3, cex=1)
  
  ## second page
  
  layout(matrix(1:8, nrow=2, ncol=4, byrow=TRUE))
  par(oma=c(3,1,3,1))
  par(mar=c(4,3,4,3)+0.1)
  
  for(originClass in c("EarlyOrigin", "LateOrigin")){
    originLabel = if(originClass == "EarlyOrigin") "Early" else "Late"
    
    for(i in seq_along(metrics)){
      PlotProfile(SampleRes[[originClass]], metrics[i],
                  paste0(originLabel, ": ", metricTitles[i]), metricYlabels[i])
    }
  }
  
  mtext("Read enrichments at early and late origins", side=3, line=0,
        outer=TRUE, font=2, cex=2)
  mtext("Page 2", side=1, line=0, outer=TRUE, font=3, cex=1)
  
  ## third page
  
  layout(matrix(1:12, nrow=3, ncol=4, byrow=TRUE))
  par(oma=c(3,1,3,1))
  par(mar=c(4,4,4,2)+0.1)
  
  for(i in seq_along(metrics)){
    PlotPairwise("EarlyOrigin", "LateOrigin", metrics[i],
                 metricTitles[i], metricYlabels[i])
  }
  
  for(i in 1:8){
    plot.new()
  }
  
  mtext("Pairwise comparative enrichment at early and late origins",
        side=3, line=0, outer=TRUE, font=2, cex=2)
  mtext("Page 3", side=1, line=0, outer=TRUE, font=3, cex=1)
  
  dev.off()
  
  message(paste0("Plot saved as ", basename(OutputFile)))
}


ChIP_BrDU_Analysis_early_late <- function(Directory = "None",
                                           Genotype = "WT",
                                           Window = "None",
                                           Alignment = "generic",
                                           yLim = "fixed",
                                           Metric = "ratio.ipin.noise"){
  
  ## load packages
  
  packages <- c("BSgenome.Scerevisiae.UCSC.sacCer3")
  
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  ## load and process basic files
  
  E_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/E_Rep.bed",
                      header=FALSE, sep="\t", quote="", comment.char="")[ ,1:4]
  L_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/L_Rep.bed",
                      header=FALSE, sep="\t", quote="", comment.char="")[ ,1:4]
  
  colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName")
  colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName")
  
  E_Ori$oriCenter <- round((E_Ori$chromStart + E_Ori$chromEnd)/2)
  L_Ori$oriCenter <- round((L_Ori$chromStart + L_Ori$chromEnd)/2)
  
  E_Ori$peakSummit <- E_Ori$oriCenter
  L_Ori$peakSummit <- L_Ori$oriCenter
  
  Alignment <- match.arg(Alignment, c("generic", "malign"))
  Genotype <- match.arg(Genotype, c("WT", "EQ"))
  yLim <- match.arg(yLim, c("fixed", "per-plot"))
  Metric <- match.arg(Metric, c("ip.score", "ratio.ipin", "ratio.ipnoise", "ratio.ipin.noise"))
  
  useDef <- function(a,d) if(length(a) == 0 || is.null(a) || is.na(a) || a == "") d else a
  
  Directory = useDef(Directory, "None")
  
  if(Directory == "None"){
    stop("Directory must be provided.")
  }
  
  if(!dir.exists(Directory)){
    stop("Analysis directory is missing: ", Directory)
  }
  
  Directory = normalizePath(Directory, winslash="/", mustWork=TRUE)
  RatioFolder = if(Alignment == "generic") "Ratios" else "Ratios_ma"
  
  SampleFolders = list.dirs(Directory, full.names=TRUE, recursive=FALSE)
  SampleNames = basename(SampleFolders)
  GenotypePattern = paste0("-", Genotype, "-")
  
  ChIPFolders = SampleFolders[
    grepl(GenotypePattern, SampleNames, fixed=TRUE) & grepl("-ChIP$", SampleNames)
  ]
  
  BrDUFolders = SampleFolders[
    grepl(GenotypePattern, SampleNames, fixed=TRUE) & grepl("-BrDU$", SampleNames)
  ]
  
  ChIPPrefixes = sub("-ChIP$", "", basename(ChIPFolders))
  BrDUPrefixes = sub("-BrDU$", "", basename(BrDUFolders))
  PairPrefixes = intersect(ChIPPrefixes, BrDUPrefixes)
  
  if(length(PairPrefixes) != 4){
    stop("Exactly four ChIP-BrDU sister pairs were expected for ", Genotype,
         " in ", Directory, "; found ", length(PairPrefixes), ".")
  }
  
  ExtractSampleInformation <- function(PairPrefix){
    Tokens = strsplit(PairPrefix, "-", fixed=TRUE)[[1]]
    GenotypePosition = which(Tokens == Genotype)[1]
    
    if(is.na(GenotypePosition) || GenotypePosition == length(Tokens)){
      stop("Could not determine the condition from sample name: ", PairPrefix)
    }
    
    Condition = Tokens[GenotypePosition+1]
    Year = if(GenotypePosition+2 <= length(Tokens)) Tokens[GenotypePosition+2] else basename(Directory)
    NumericTime = suppressWarnings(as.numeric(gsub("[^0-9]", "", Condition)))
    
    if(grepl("^G1$", Condition, ignore.case=TRUE)){
      NumericTime = -1
    }
    
    if(!is.finite(NumericTime)){
      NumericTime = Inf
    }
    
    data.frame(PairPrefix=PairPrefix, Condition=Condition, Year=Year,
               NumericTime=NumericTime, stringsAsFactors=FALSE)
  }
  
  PairInfo = do.call(rbind, lapply(PairPrefixes, ExtractSampleInformation))
  PairInfo = PairInfo[order(PairInfo$NumericTime, PairInfo$Condition), ]
  rownames(PairInfo) = NULL
  
  PairInfo$ChIPFolder = file.path(Directory, paste0(PairInfo$PairPrefix, "-ChIP"))
  PairInfo$BrDUFolder = file.path(Directory, paste0(PairInfo$PairPrefix, "-BrDU"))
  
  if(any(!dir.exists(PairInfo$ChIPFolder)) || any(!dir.exists(PairInfo$BrDUFolder))){
    stop("One or more ChIP-BrDU sister folders are missing.")
  }
  
  if(Window == "None"){
    FirstSample = basename(PairInfo$ChIPFolder[1])
    InfoCandidates = c(
      paste0(PairInfo$ChIPFolder[1], "/", FirstSample, "_Info.txt"),
      paste0(PairInfo$ChIPFolder[1], "/Info/", FirstSample, "_Info.txt"),
      paste0(PairInfo$ChIPFolder[1], "/Coverage/", FirstSample, "_Info.txt")
    )
    
    InfoFile = InfoCandidates[file.exists(InfoCandidates)][1]
    
    if(is.na(InfoFile)){
      stop("Window could not be resolved. Provide Window explicitly.")
    }
    
    Window = as.numeric(read.table(InfoFile, header=FALSE)[5,2])
  } else {
    Window = as.numeric(Window)
  }
  
  if(!is.finite(Window) || Window <= 0){
    stop("Window must be a positive number.")
  }
  
  MetricLabel = switch(
    Metric,
    "ip.score" = "IP score",
    "ratio.ipin" = "IP/Input enrichment",
    "ratio.ipnoise" = "IP/noise enrichment",
    "ratio.ipin.noise" = "clean enrichment"
  )
  OriginLists = list(EarlyOrigin=E_Ori, LateOrigin=L_Ori)
  
  ReadRatio <- function(RatioFile){
    
    if(!file.exists(RatioFile)){
      stop("Ratio file is missing: ", RatioFile)
    }
    
    Header = strsplit(readLines(RatioFile, n=1), "\t", fixed=TRUE)[[1]]
    RequiredColumns = c("chrom", "chromStart", "chromEnd", Metric)
    MissingColumns = RequiredColumns[!RequiredColumns %in% Header]
    
    if(length(MissingColumns) > 0){
      stop("Required ratio columns are missing from ", basename(RatioFile), ": ",
           paste(MissingColumns, collapse=", "))
    }
    
    ColumnClasses = ifelse(Header %in% RequiredColumns, NA, "NULL")
    Ratio = read.table(RatioFile, header=TRUE, check.names=FALSE,
                       colClasses=ColumnClasses)
    
    Ratio$chromStart = as.numeric(Ratio$chromStart)
    Ratio$chromEnd = as.numeric(Ratio$chromEnd)
    Ratio[,Metric] = as.numeric(Ratio[,Metric])
    Ratio[,Metric][!is.finite(Ratio[,Metric])] = 0
    
    Ratio
  }
  
  infer_step <- function(Ratio){
    
    StartsByChrom = split(Ratio$chromStart, as.character(Ratio$chrom))
    Diffs = unlist(
      lapply(StartsByChrom, function(Starts){
        Starts = sort(unique(Starts[is.finite(Starts)]))
        diff(Starts)
      }),
      use.names=FALSE
    )
    
    Diffs = Diffs[is.finite(Diffs) & Diffs > 0]
    
    if(length(Diffs) == 0){
      stop("Could not infer coverage step size.")
    }
    
    DiffTable = table(Diffs)
    as.numeric(names(DiffTable)[which.max(DiffTable)])
  }
  
  extract_centered_values <- function(Ratio, center, step, expected_n){
    
    regionStart = center - Window
    regionEnd = center + Window
    
    if(is.null(Ratio)){
      return(rep(0, expected_n))
    }
    
    Signal = Ratio[
      Ratio$chromStart >= regionStart & Ratio$chromStart <= regionEnd, ]
    
    Signal = Signal[order(Signal$chromStart), ]
    
    if(nrow(Signal) == 0){
      return(rep(0, expected_n))
    }
    
    Values = Signal[,Metric]
    Values[!is.finite(Values)] = 0
    
    if(length(Values) < expected_n){
      missing_n = expected_n - length(Values)
      left_missing = round((min(Signal$chromStart, na.rm=TRUE) - regionStart)/step)
      left_missing = min(max(left_missing, 0L), missing_n)
      right_missing = missing_n - left_missing
      Values = c(rep(0, left_missing), Values, rep(0, right_missing))
    }
    
    if(length(Values) > expected_n){
      Values = Values[seq_len(expected_n)]
    }
    
    Values
  }
  
  CalculateAssayProfiles <- function(RatioFile, Assay){
    
    message(paste0("Calculating ", Assay, " profiles from ", basename(RatioFile)))
    
    Ratio = ReadRatio(RatioFile)
    step = infer_step(Ratio)
    expected_n = max(1L, round(2*Window/step))
    RatioByChrom = split(Ratio, as.character(Ratio$chrom))
    
    Profiles = lapply(OriginLists, function(OriginList){
      
      ProfileMatrix = vapply(
        seq_len(nrow(OriginList)),
        function(i){
          extract_centered_values(
            Ratio=RatioByChrom[[as.character(OriginList$chrom[i])]],
            center=OriginList$peakSummit[i], step=step,
            expected_n=expected_n
          )
        },
        numeric(expected_n)
      )
      
      if(is.null(dim(ProfileMatrix))){
        ProfileMatrix = matrix(ProfileMatrix, ncol=1)
      }
      
      apply(ProfileMatrix, 1, median, na.rm=TRUE)
    })
    
    list(Profiles=Profiles, step=step)
  }
  
  PairResults = lapply(seq_len(nrow(PairInfo)), function(i){
    
    ChIPName = basename(PairInfo$ChIPFolder[i])
    BrDUName = basename(PairInfo$BrDUFolder[i])
    
    ChIPRatioFile = paste0(PairInfo$ChIPFolder[i], "/", RatioFolder, "/",
                           ChIPName, "_ChIP_collapsed.bed")
    BrDURatioFile = paste0(PairInfo$BrDUFolder[i], "/", RatioFolder, "/",
                           BrDUName, "_BrDU_collapsed.bed")
    
    message(paste0("Processing sister pair ", PairInfo$PairPrefix[i]))
    
    ChIPResult = CalculateAssayProfiles(ChIPRatioFile, "ChIP")
    BrDUResult = CalculateAssayProfiles(BrDURatioFile, "BrDU")
    
    if(!isTRUE(all.equal(ChIPResult$step, BrDUResult$step))){
      stop("ChIP and BrDU step sizes differ for ", PairInfo$PairPrefix[i], ".")
    }
    
    gc(verbose=FALSE)
    
    list(Condition=PairInfo$Condition[i], Year=PairInfo$Year[i],
         ChIP=ChIPResult$Profiles, BrDU=BrDUResult$Profiles,
         step=ChIPResult$step)
  })
  
  Steps = vapply(PairResults, function(x) x$step, numeric(1))
  
  if(any(abs(Steps - Steps[1]) > .Machine$double.eps^0.5)){
    stop("The four sister pairs do not use the same coverage step size.")
  }
  
  SmoothProfile <- function(Profile){
    
    Profile[!is.finite(Profile)] = 0
    
    if(length(Profile) < 4 || length(unique(Profile)) < 2){
      return(Profile)
    }
    
    smooth.spline(seq_along(Profile), Profile, spar=0.5)$y
  }
  
  for(i in seq_along(PairResults)){
    PairResults[[i]]$ChIP = lapply(PairResults[[i]]$ChIP, SmoothProfile)
    PairResults[[i]]$BrDU = lapply(PairResults[[i]]$BrDU, SmoothProfile)
  }
  
  ProfileYlim <- function(Profiles){
    
    Values = unlist(Profiles)
    Values = Values[is.finite(Values)]
    
    if(length(Values) == 0){
      return(c(-1, 1))
    }
    
    Y_min = min(Values)
    Y_max = max(Values)
    Y_pad = 0.08*(Y_max - Y_min)
    
    if(!is.finite(Y_pad) || Y_pad == 0){
      Y_pad = max(0.5, abs(Y_max)*0.08)
    }
    
    c(Y_min-Y_pad, Y_max+Y_pad)
  }
  
  AlignYlimsAtOne <- function(ChIPLimits, BrDULimits){
    
    ChIPLimits[1] = min(ChIPLimits[1], 1)
    ChIPLimits[2] = max(ChIPLimits[2], 1)
    BrDULimits[1] = min(BrDULimits[1], 1)
    BrDULimits[2] = max(BrDULimits[2], 1)
    
    ChIPFraction = (1-ChIPLimits[1])/diff(ChIPLimits)
    BrDUFraction = (1-BrDULimits[1])/diff(BrDULimits)
    AlignedFraction = max(ChIPFraction, BrDUFraction)
    
    if(AlignedFraction > 0 && AlignedFraction < 1){
      ChIPLimits[1] = (1-AlignedFraction*ChIPLimits[2])/(1-AlignedFraction)
      BrDULimits[1] = (1-AlignedFraction*BrDULimits[2])/(1-AlignedFraction)
    }
    
    list(ChIP=ChIPLimits, BrDU=BrDULimits)
  }
  
  PrepareYlims <- function(ChIPProfiles, BrDUProfiles){
    
    ChIPLimits = ProfileYlim(ChIPProfiles)
    BrDULimits = ProfileYlim(BrDUProfiles)
    AlignedLimits = AlignYlimsAtOne(ChIPLimits, BrDULimits)
    
    BrDULimits = AlignedLimits$BrDU
    BrDULimits[1] = min(BrDULimits[1], 1)
    BrDULimits[2] = max(BrDULimits[2], 1)
    BrDUAxisTicks = pretty(BrDULimits)
    BrDUAxisTicks = BrDUAxisTicks[
      BrDUAxisTicks >= 1 & BrDUAxisTicks <= BrDULimits[2]
    ]
    BrDUAxisTicks = sort(unique(c(1, BrDUAxisTicks)))
    BrDUAxisLabels = format(BrDUAxisTicks, trim=TRUE, nsmall=1)
    
    list(ChIP=AlignedLimits$ChIP, BrDU=BrDULimits,
         BrDUAxisTicks=BrDUAxisTicks, BrDUAxisLabels=BrDUAxisLabels)
  }
  
  FixedYlims = PrepareYlims(
    unlist(lapply(PairResults, function(x) x$ChIP), recursive=FALSE),
    unlist(lapply(PairResults, function(x) x$BrDU), recursive=FALSE)
  )

  PairYlims = lapply(PairResults, function(x){
    PrepareYlims(x$ChIP, x$BrDU)
  })
  
  EarlyColor = adjustcolor("darkorchid4", alpha.f=0.9)
  LateColor = adjustcolor("darkorange3", alpha.f=0.9)
  
  PlotPanel <- function(Result, OriginClass, PanelYlims){
    
    OriginLabel = if(OriginClass == "EarlyOrigin") "Early" else "Late"
    ProfileColor = if(OriginClass == "EarlyOrigin") EarlyColor else LateColor
    ChIPProfile = Result$ChIP[[OriginClass]]
    BrDUProfile = Result$BrDU[[OriginClass]]
    
    ChIPYlim = PanelYlims$ChIP
    BrDUYlim = PanelYlims$BrDU
    BrDUAxisTicks = PanelYlims$BrDUAxisTicks
    BrDUAxisLabels = PanelYlims$BrDUAxisLabels
    
    plot(ChIPProfile, ylim=ChIPYlim,
         main=paste0(Result$Condition, ": ", OriginLabel),
         ylab=paste0("ChIP ", MetricLabel), xlab="Distance from OriCenter (Kbp)",
         xaxt="n", col=ProfileColor, type="l", lty=1, lwd=2,
         bty="l", las=1, xaxs="i", yaxs="i",
         cex.axis=0.85, cex.lab=0.82, cex.main=0.9)
    
    axisLabels = seq(-Window, +Window, length.out=9)
    axisLabels[c(2,4,6,8)] = NA
    At = seq(1, length(ChIPProfile), length.out=9)
    axis(1, at=At, labels=signif(axisLabels/1000, 2), col="black",
         col.axis="black", col.ticks="black", cex.axis=0.85)
    
    par(new=TRUE)
    
    plot(BrDUProfile, ylim=BrDUYlim, axes=FALSE, xlab="", ylab="",
         col=ProfileColor, type="l", lty=3, lwd=2,
         bty="n", xaxs="i", yaxs="i")
    
    axis(4, at=BrDUAxisTicks, labels=BrDUAxisLabels, las=1,
         col="transparent", col.axis="black", col.ticks="black",
         lty=1, lwd=0, lwd.ticks=1.3, tck=-0.018, cex.axis=0.82)
    AxisLimits = par("usr")
    segments(AxisLimits[2], AxisLimits[3], AxisLimits[2], AxisLimits[4],
             col="black", lty=3, lwd=1.6, xpd=NA)
    mtext(paste0("BrDU ", MetricLabel), side=4, line=2.6,
          col="black", cex=0.82)
    
    legend("topright", legend=c("ChIP", "BrDU"),
           col=c(ProfileColor, ProfileColor), lty=c(1,3),
           lwd=2, bty="n", cex=0.7)
  }
  
  YearLabel = unique(PairInfo$Year)
  if(length(YearLabel) != 1) YearLabel = basename(Directory)
  MetricSuffix = if(Metric == "ratio.ipin.noise") "" else {
    paste0("_", gsub("\\.", "_", Metric))
  }
  
  OutputFile = paste0(Directory, "/Smc5-", Genotype, "-", YearLabel,
                      "_ChIP_BrDU_early_late_", Alignment,
                      if(yLim == "per-plot") "_per_plot_ylim" else "",
                      MetricSuffix, ".pdf")
  
  pdf(OutputFile, width=16, height=10.5)
  
  layout(matrix(1:12, nrow=3, ncol=4, byrow=TRUE))
  par(oma=c(2.8,1.8,3.2,1.8))
  par(mar=c(4,4,2.8,4)+0.1, mgp=c(2.3,0.7,0), tcl=-0.25)
  
  for(i in seq_along(PairResults)){
    PanelYlims = if(yLim == "fixed") FixedYlims else PairYlims[[i]]
    PlotPanel(PairResults[[i]], "EarlyOrigin", PanelYlims)
  }
  for(i in seq_along(PairResults)){
    PanelYlims = if(yLim == "fixed") FixedYlims else PairYlims[[i]]
    PlotPanel(PairResults[[i]], "LateOrigin", PanelYlims)
  }
  for(i in 1:4){
    plot.new()
  }
  
  mtext(paste0("Smc5 ", Genotype, " ", YearLabel,
               ": ChIP and BrDU enrichment at early and late origins"),
        side=3, line=1, outer=TRUE, font=2, cex=1.4)
  
  dev.off()
  
  message(paste0("Plot saved as ", OutputFile))
  
  invisible(list(results=PairResults, pair_info=PairInfo,
                 output_pdf=OutputFile, Window=Window, step=Steps[1],
                 yLim=yLim, Metric=Metric))
}


ChIP_BrDU_Analysis_early_late_old_approach <- function(Directory = "None",
                                                        Genotype = "WT",
                                                        Window = "None",
                                                        Alignment = "generic",
                                                        yLim = "fixed"){
  
  ## load packages
  
  packages <- c("BSgenome.Scerevisiae.UCSC.sacCer3")
  
  suppressWarnings(suppressPackageStartupMessages(lapply(packages, require, character.only = TRUE)))
  
  ## load and process basic files
  
  E_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/E_Rep.bed",
                      header=FALSE, sep="\t", quote="", comment.char="")[ ,1:4]
  L_Ori <- read.table("/Applications/ngsAnalyser.app/Contents/Resources/app/L_Rep.bed",
                      header=FALSE, sep="\t", quote="", comment.char="")[ ,1:4]
  
  colnames(E_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName")
  colnames(L_Ori) <- c("chrom", "chromStart", "chromEnd", "oriName")
  
  E_Ori$oriCenter <- round((E_Ori$chromStart + E_Ori$chromEnd)/2)
  L_Ori$oriCenter <- round((L_Ori$chromStart + L_Ori$chromEnd)/2)
  
  E_Ori$peakSummit <- E_Ori$oriCenter
  L_Ori$peakSummit <- L_Ori$oriCenter
  
  Alignment <- match.arg(Alignment, c("generic", "malign"))
  Genotype <- match.arg(Genotype, c("WT", "EQ"))
  yLim <- match.arg(yLim, c("fixed", "per-plot"))
  
  useDef <- function(a,d) if(length(a) == 0 || is.null(a) || is.na(a) || a == "") d else a
  
  Directory = useDef(Directory, "None")
  
  if(Directory == "None"){
    stop("Directory must be provided.")
  }
  
  if(!dir.exists(Directory)){
    stop("Analysis directory is missing: ", Directory)
  }
  
  Directory = normalizePath(Directory, winslash="/", mustWork=TRUE)
  RatioFolder = if(Alignment == "generic") "Ratios" else "Ratios_ma"
  
  SampleFolders = list.dirs(Directory, full.names=TRUE, recursive=FALSE)
  SampleNames = basename(SampleFolders)
  GenotypePattern = paste0("-", Genotype, "-")
  
  ChIPFolders = SampleFolders[
    grepl(GenotypePattern, SampleNames, fixed=TRUE) & grepl("-ChIP$", SampleNames)
  ]
  
  BrDUFolders = SampleFolders[
    grepl(GenotypePattern, SampleNames, fixed=TRUE) & grepl("-BrDU$", SampleNames)
  ]
  
  ChIPPrefixes = sub("-ChIP$", "", basename(ChIPFolders))
  BrDUPrefixes = sub("-BrDU$", "", basename(BrDUFolders))
  PairPrefixes = intersect(ChIPPrefixes, BrDUPrefixes)
  
  if(length(PairPrefixes) != 4){
    stop("Exactly four ChIP-BrDU sister pairs were expected for ", Genotype,
         " in ", Directory, "; found ", length(PairPrefixes), ".")
  }
  
  ExtractSampleInformation <- function(PairPrefix){
    Tokens = strsplit(PairPrefix, "-", fixed=TRUE)[[1]]
    GenotypePosition = which(Tokens == Genotype)[1]
    
    if(is.na(GenotypePosition) || GenotypePosition == length(Tokens)){
      stop("Could not determine the condition from sample name: ", PairPrefix)
    }
    
    Condition = Tokens[GenotypePosition+1]
    Year = if(GenotypePosition+2 <= length(Tokens)) Tokens[GenotypePosition+2] else basename(Directory)
    NumericTime = suppressWarnings(as.numeric(gsub("[^0-9]", "", Condition)))
    
    if(grepl("^G1$", Condition, ignore.case=TRUE)){
      NumericTime = -1
    }
    
    if(!is.finite(NumericTime)){
      NumericTime = Inf
    }
    
    data.frame(PairPrefix=PairPrefix, Condition=Condition, Year=Year,
               NumericTime=NumericTime, stringsAsFactors=FALSE)
  }
  
  PairInfo = do.call(rbind, lapply(PairPrefixes, ExtractSampleInformation))
  PairInfo = PairInfo[order(PairInfo$NumericTime, PairInfo$Condition), ]
  rownames(PairInfo) = NULL
  
  PairInfo$ChIPFolder = file.path(Directory, paste0(PairInfo$PairPrefix, "-ChIP"))
  PairInfo$BrDUFolder = file.path(Directory, paste0(PairInfo$PairPrefix, "-BrDU"))
  
  if(any(!dir.exists(PairInfo$ChIPFolder)) || any(!dir.exists(PairInfo$BrDUFolder))){
    stop("One or more ChIP-BrDU sister folders are missing.")
  }
  
  if(Window == "None"){
    FirstSample = basename(PairInfo$ChIPFolder[1])
    InfoCandidates = c(
      paste0(PairInfo$ChIPFolder[1], "/", FirstSample, "_Info.txt"),
      paste0(PairInfo$ChIPFolder[1], "/Info/", FirstSample, "_Info.txt"),
      paste0(PairInfo$ChIPFolder[1], "/Coverage/", FirstSample, "_Info.txt")
    )
    
    InfoFile = InfoCandidates[file.exists(InfoCandidates)][1]
    
    if(is.na(InfoFile)){
      stop("Window could not be resolved. Provide Window explicitly.")
    }
    
    Window = as.numeric(read.table(InfoFile, header=FALSE)[5,2])
  } else {
    Window = as.numeric(Window)
  }
  
  if(!is.finite(Window) || Window <= 0){
    stop("Window must be a positive number.")
  }
  
  OriginLists = list(EarlyOrigin=E_Ori, LateOrigin=L_Ori)
  
  ReadRatio <- function(RatioFile){
    
    if(!file.exists(RatioFile)){
      stop("Ratio file is missing: ", RatioFile)
    }
    
    Header = strsplit(readLines(RatioFile, n=1), "\t", fixed=TRUE)[[1]]
    RequiredColumns = c("chrom", "chromStart", "chromEnd", "ip.score", "in.score")
    MissingColumns = RequiredColumns[!RequiredColumns %in% Header]
    
    if(length(MissingColumns) > 0){
      stop("Required ratio columns are missing from ", basename(RatioFile), ": ",
           paste(MissingColumns, collapse=", "))
    }
    
    ColumnClasses = ifelse(Header %in% RequiredColumns, NA, "NULL")
    Ratio = read.table(RatioFile, header=TRUE, check.names=FALSE,
                       colClasses=ColumnClasses)
    
    Ratio$chromStart = as.numeric(Ratio$chromStart)
    Ratio$chromEnd = as.numeric(Ratio$chromEnd)
    Ratio$ip.score = as.numeric(Ratio$ip.score)
    Ratio$in.score = as.numeric(Ratio$in.score)
    Ratio$ip.score[!is.finite(Ratio$ip.score)] = 0
    Ratio$in.score[!is.finite(Ratio$in.score)] = 0
    
    Ratio
  }
  
  infer_step <- function(Ratio){
    
    StartsByChrom = split(Ratio$chromStart, as.character(Ratio$chrom))
    Diffs = unlist(
      lapply(StartsByChrom, function(Starts){
        Starts = sort(unique(Starts[is.finite(Starts)]))
        diff(Starts)
      }),
      use.names=FALSE
    )
    
    Diffs = Diffs[is.finite(Diffs) & Diffs > 0]
    
    if(length(Diffs) == 0){
      stop("Could not infer coverage step size.")
    }
    
    DiffTable = table(Diffs)
    as.numeric(names(DiffTable)[which.max(DiffTable)])
  }
  
  extract_centered_values <- function(Ratio, center, step, expected_n, ValueColumn){
    
    regionStart = center - Window
    regionEnd = center + Window
    
    if(is.null(Ratio)){
      return(rep(0, expected_n))
    }
    
    Signal = Ratio[
      Ratio$chromStart >= regionStart & Ratio$chromStart <= regionEnd, ]
    
    Signal = Signal[order(Signal$chromStart), ]
    
    if(nrow(Signal) == 0){
      return(rep(0, expected_n))
    }
    
    Values = Signal[,ValueColumn]
    Values[!is.finite(Values)] = 0
    
    if(length(Values) < expected_n){
      missing_n = expected_n - length(Values)
      left_missing = round((min(Signal$chromStart, na.rm=TRUE) - regionStart)/step)
      left_missing = min(max(left_missing, 0L), missing_n)
      right_missing = missing_n - left_missing
      Values = c(rep(0, left_missing), Values, rep(0, right_missing))
    }
    
    if(length(Values) > expected_n){
      Values = Values[seq_len(expected_n)]
    }
    
    Values
  }
  
  CalculateMedianProfile <- function(RatioByChrom, OriginList, step,
                                     expected_n, ValueColumn){
    
    ProfileMatrix = vapply(
      seq_len(nrow(OriginList)),
      function(i){
        extract_centered_values(
          Ratio=RatioByChrom[[as.character(OriginList$chrom[i])]],
          center=OriginList$peakSummit[i], step=step,
          expected_n=expected_n, ValueColumn=ValueColumn
        )
      },
      numeric(expected_n)
    )
    
    if(is.null(dim(ProfileMatrix))){
      ProfileMatrix = matrix(ProfileMatrix, ncol=1)
    }
    
    apply(ProfileMatrix, 1, median, na.rm=TRUE)
  }
  
  CalculateAssayProfiles <- function(WatsonFile, CrickFile, Assay){
    
    message(paste0("Calculating old-approach ", Assay,
                   " profiles from strand-separated IP and Input scores"))
    
    WatsonRatio = ReadRatio(WatsonFile)
    CrickRatio = ReadRatio(CrickFile)
    WatsonStep = infer_step(WatsonRatio)
    CrickStep = infer_step(CrickRatio)
    
    if(!isTRUE(all.equal(WatsonStep, CrickStep))){
      stop("Watson and Crick step sizes differ for ", Assay, ".")
    }
    
    step = WatsonStep
    expected_n = max(1L, round(2*Window/step))
    WatsonByChrom = split(WatsonRatio, as.character(WatsonRatio$chrom))
    CrickByChrom = split(CrickRatio, as.character(CrickRatio$chrom))
    
    Profiles = lapply(OriginLists, function(OriginList){
      
      WatsonIP = CalculateMedianProfile(
        WatsonByChrom, OriginList, step, expected_n, "ip.score"
      )
      WatsonInput = CalculateMedianProfile(
        WatsonByChrom, OriginList, step, expected_n, "in.score"
      )
      CrickIP = CalculateMedianProfile(
        CrickByChrom, OriginList, step, expected_n, "ip.score"
      )
      CrickInput = CalculateMedianProfile(
        CrickByChrom, OriginList, step, expected_n, "in.score"
      )

      WatsonInNorm = WatsonIP/WatsonInput
      CrickInNorm = CrickIP/CrickInput
      WatsonInNorm[!is.finite(WatsonInNorm)] = 0
      CrickInNorm[!is.finite(CrickInNorm)] = 0
      
      Watson = smooth.spline(
        1:length(WatsonInNorm), WatsonInNorm, spar=0.5
      )$y
      Crick = smooth.spline(
        1:length(CrickInNorm), CrickInNorm, spar=0.5
      )$y
      
      Watson + Crick
    })
    
    list(Profiles=Profiles, step=step)
  }
  
  PairResults = lapply(seq_len(nrow(PairInfo)), function(i){
    
    ChIPName = basename(PairInfo$ChIPFolder[i])
    BrDUName = basename(PairInfo$BrDUFolder[i])
    
    ChIPWatsonFile = paste0(PairInfo$ChIPFolder[i], "/", RatioFolder, "/",
                            ChIPName, "_ChIP_watson.bed")
    ChIPCrickFile = paste0(PairInfo$ChIPFolder[i], "/", RatioFolder, "/",
                           ChIPName, "_ChIP_crick.bed")
    BrDUWatsonFile = paste0(PairInfo$BrDUFolder[i], "/", RatioFolder, "/",
                            BrDUName, "_BrDU_watson.bed")
    BrDUCrickFile = paste0(PairInfo$BrDUFolder[i], "/", RatioFolder, "/",
                           BrDUName, "_BrDU_crick.bed")
    
    message(paste0("Processing sister pair ", PairInfo$PairPrefix[i]))
    
    ChIPResult = CalculateAssayProfiles(
      ChIPWatsonFile, ChIPCrickFile, "ChIP"
    )
    BrDUResult = CalculateAssayProfiles(
      BrDUWatsonFile, BrDUCrickFile, "BrDU"
    )
    
    if(!isTRUE(all.equal(ChIPResult$step, BrDUResult$step))){
      stop("ChIP and BrDU step sizes differ for ", PairInfo$PairPrefix[i], ".")
    }
    
    gc(verbose=FALSE)
    
    list(Condition=PairInfo$Condition[i], Year=PairInfo$Year[i],
         ChIP=ChIPResult$Profiles, BrDU=BrDUResult$Profiles,
         step=ChIPResult$step)
  })
  
  Steps = vapply(PairResults, function(x) x$step, numeric(1))
  
  if(any(abs(Steps - Steps[1]) > .Machine$double.eps^0.5)){
    stop("The four sister pairs do not use the same coverage step size.")
  }
  
  ProfileYlim <- function(Profiles){
    
    Values = unlist(Profiles)
    Values = Values[is.finite(Values)]
    
    if(length(Values) == 0){
      return(c(-1, 1))
    }
    
    Y_min = min(Values)
    Y_max = max(Values)
    Y_pad = 0.08*(Y_max - Y_min)
    
    if(!is.finite(Y_pad) || Y_pad == 0){
      Y_pad = max(0.5, abs(Y_max)*0.08)
    }
    
    c(Y_min-Y_pad, Y_max+Y_pad)
  }
  
  AlignYlimsAtOne <- function(ChIPLimits, BrDULimits){
    
    ChIPLimits[1] = min(ChIPLimits[1], 1)
    ChIPLimits[2] = max(ChIPLimits[2], 1)
    BrDULimits[1] = min(BrDULimits[1], 1)
    BrDULimits[2] = max(BrDULimits[2], 1)
    
    ChIPFraction = (1-ChIPLimits[1])/diff(ChIPLimits)
    BrDUFraction = (1-BrDULimits[1])/diff(BrDULimits)
    AlignedFraction = max(ChIPFraction, BrDUFraction)
    
    if(AlignedFraction > 0 && AlignedFraction < 1){
      ChIPLimits[1] = (1-AlignedFraction*ChIPLimits[2])/(1-AlignedFraction)
      BrDULimits[1] = (1-AlignedFraction*BrDULimits[2])/(1-AlignedFraction)
    }
    
    list(ChIP=ChIPLimits, BrDU=BrDULimits)
  }
  
  PrepareYlims <- function(ChIPProfiles, BrDUProfiles){
    
    ChIPLimits = ProfileYlim(ChIPProfiles)
    BrDULimits = ProfileYlim(BrDUProfiles)
    AlignedLimits = AlignYlimsAtOne(ChIPLimits, BrDULimits)
    
    BrDULimits = AlignedLimits$BrDU
    BrDULimits[1] = min(BrDULimits[1], 1)
    BrDULimits[2] = max(BrDULimits[2], 1)
    BrDUAxisTicks = pretty(BrDULimits)
    BrDUAxisTicks = BrDUAxisTicks[
      BrDUAxisTicks >= 1 & BrDUAxisTicks <= BrDULimits[2]
    ]
    BrDUAxisTicks = sort(unique(c(1, BrDUAxisTicks)))
    BrDUAxisLabels = format(BrDUAxisTicks, trim=TRUE, nsmall=1)
    
    list(ChIP=AlignedLimits$ChIP, BrDU=BrDULimits,
         BrDUAxisTicks=BrDUAxisTicks, BrDUAxisLabels=BrDUAxisLabels)
  }
  
  FixedYlims = PrepareYlims(
    unlist(lapply(PairResults, function(x) x$ChIP), recursive=FALSE),
    unlist(lapply(PairResults, function(x) x$BrDU), recursive=FALSE)
  )

  PairYlims = lapply(PairResults, function(x){
    PrepareYlims(x$ChIP, x$BrDU)
  })
  
  EarlyColor = adjustcolor("darkorchid4", alpha.f=0.9)
  LateColor = adjustcolor("darkorange3", alpha.f=0.9)
  
  PlotPanel <- function(Result, OriginClass, PanelYlims){
    
    OriginLabel = if(OriginClass == "EarlyOrigin") "Early" else "Late"
    ProfileColor = if(OriginClass == "EarlyOrigin") EarlyColor else LateColor
    ChIPProfile = Result$ChIP[[OriginClass]]
    BrDUProfile = Result$BrDU[[OriginClass]]
    
    ChIPYlim = PanelYlims$ChIP
    BrDUYlim = PanelYlims$BrDU
    BrDUAxisTicks = PanelYlims$BrDUAxisTicks
    BrDUAxisLabels = PanelYlims$BrDUAxisLabels
    
    plot(ChIPProfile, ylim=ChIPYlim,
         main=paste0(Result$Condition, ": ", OriginLabel),
         ylab="ChIP / Input enrichment", xlab="Distance from OriCenter (Kbp)",
         xaxt="n", col=ProfileColor, type="l", lty=1, lwd=2,
         bty="l", las=1, xaxs="i", yaxs="i",
         cex.axis=0.85, cex.lab=0.82, cex.main=0.9)
    
    axisLabels = seq(-Window, +Window, length.out=9)
    axisLabels[c(2,4,6,8)] = NA
    At = seq(1, length(ChIPProfile), length.out=9)
    axis(1, at=At, labels=signif(axisLabels/1000, 2), col="black",
         col.axis="black", col.ticks="black", cex.axis=0.85)
    
    par(new=TRUE)
    
    plot(BrDUProfile, ylim=BrDUYlim, axes=FALSE, xlab="", ylab="",
         col=ProfileColor, type="l", lty=3, lwd=2,
         bty="n", xaxs="i", yaxs="i")
    
    axis(4, at=BrDUAxisTicks, labels=BrDUAxisLabels, las=1,
         col="transparent", col.axis="black", col.ticks="black",
         lty=1, lwd=0, lwd.ticks=1.3, tck=-0.018, cex.axis=0.82)
    AxisLimits = par("usr")
    segments(AxisLimits[2], AxisLimits[3], AxisLimits[2], AxisLimits[4],
             col="black", lty=3, lwd=1.6, xpd=NA)
    mtext("BrDU / Input enrichment", side=4, line=2.6,
          col="black", cex=0.82)
    
    legend("topright", legend=c("ChIP", "BrDU"),
           col=c(ProfileColor, ProfileColor), lty=c(1,3),
           lwd=2, bty="n", cex=0.7)
  }
  
  YearLabel = unique(PairInfo$Year)
  if(length(YearLabel) != 1) YearLabel = basename(Directory)
  
  OutputFile = paste0(Directory, "/Smc5-", Genotype, "-", YearLabel,
                      "_ChIP_BrDU_early_late_", Alignment,
                      if(yLim == "per-plot") "_per_plot_ylim" else "",
                      "_old_approach.pdf")
  
  pdf(OutputFile, width=16, height=10.5)
  
  layout(matrix(1:12, nrow=3, ncol=4, byrow=TRUE))
  par(oma=c(2.8,1.8,3.2,1.8))
  par(mar=c(4,4,2.8,4)+0.1, mgp=c(2.3,0.7,0), tcl=-0.25)
  
  for(i in seq_along(PairResults)){
    PanelYlims = if(yLim == "fixed") FixedYlims else PairYlims[[i]]
    PlotPanel(PairResults[[i]], "EarlyOrigin", PanelYlims)
  }
  for(i in seq_along(PairResults)){
    PanelYlims = if(yLim == "fixed") FixedYlims else PairYlims[[i]]
    PlotPanel(PairResults[[i]], "LateOrigin", PanelYlims)
  }
  for(i in 1:4){
    plot.new()
  }
  
  mtext(paste0("Smc5 ", Genotype, " ", YearLabel,
               ": ChIP and BrDU enrichment at early and late origins"),
        side=3, line=1, outer=TRUE, font=2, cex=1.4)
  
  dev.off()
  
  message(paste0("Plot saved as ", OutputFile))
  
  invisible(list(results=PairResults, pair_info=PairInfo,
                 output_pdf=OutputFile, Window=Window, step=Steps[1],
                 yLim=yLim, method="old strand-wise ratio-of-medians approach"))
}


origin_distance_pair_analysis <- function(
    sample_1_folder,
    sample_2_folder,
    origin_file = NULL,
    Alignment = "generic",
    Signal = "chip_input",
    yLim = "fixed",
    sample_1_label = basename(normalizePath(sample_1_folder, mustWork = FALSE)),
    sample_2_label = basename(normalizePath(sample_2_folder, mustWork = FALSE)),
    distance_breaks_kb = c(0, 0.5, 1.5, 2.5, 4, 8),
    window_kb = 12,
    output_pdf = NULL,
    plot_title = NULL,
    profile_colours = c("#2166AC", "#B2182B"),
    box_colours = c("#92C5DE", "#F4A582"),
    point_alpha = 0.32,
    point_cex = 0.85,
    pdf_width = 15,
    pdf_height = 8,
    expected_origin_count = NULL,
    return_data = TRUE) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  
  Alignment <- match.arg(Alignment, c("generic", "malign"))
  Signal <- match.arg(Signal, c("chip_input", "clean"))
  yLim <- match.arg(yLim, c("fixed", "per-plot"))
  
  RatioFolder <- if (Alignment == "generic") "Ratios" else "Ratios_ma"
  RatioColumn <- if (Signal == "chip_input") "ratio.ipin" else "ratio.ipin.noise"
  ChIPYAxisLabel <- if (Signal == "chip_input") {
    "ChIP/Input enrichment"
  } else {
    "Clean ChIP enrichment"
  }
  BrDUYAxisLabel <- if (Signal == "chip_input") {
    "BrDU/Input enrichment"
  } else {
    "Clean BrDU enrichment"
  }
  
  if (length(distance_breaks_kb) < 2L ||
      any(!is.finite(distance_breaks_kb)) ||
      any(diff(distance_breaks_kb) <= 0) ||
      distance_breaks_kb[1] != 0) {
    stop("distance_breaks_kb must be strictly increasing and begin at zero.")
  }
  
  if (!is.finite(window_kb) || window_kb <= 0 ||
      window_kb < max(distance_breaks_kb)) {
    stop("window_kb must be positive and cover all distance intervals.")
  }
  
  sample_1_folder <- normalizePath(sample_1_folder, mustWork = TRUE)
  sample_2_folder <- normalizePath(sample_2_folder, mustWork = TRUE)
  
  resolve_pair_folders <- function(sample_folder) {
    sample_name <- basename(sample_folder)
    
    if (!grepl("-(ChIP|BrDU)$", sample_name)) {
      stop("Sample folder must end in -ChIP or -BrDU: ", sample_folder)
    }
    
    sample_prefix <- sub("-(ChIP|BrDU)$", "", sample_name)
    sample_parent <- dirname(sample_folder)
    ChIP_folder <- file.path(sample_parent, paste0(sample_prefix, "-ChIP"))
    BrDU_folder <- file.path(sample_parent, paste0(sample_prefix, "-BrDU"))
    
    if (!dir.exists(ChIP_folder)) {
      stop("Sister ChIP folder is missing: ", ChIP_folder)
    }
    if (!dir.exists(BrDU_folder)) {
      stop("Sister BrDU folder is missing: ", BrDU_folder)
    }
    
    list(
      prefix = sample_prefix,
      ChIP_folder = normalizePath(ChIP_folder, mustWork = TRUE),
      BrDU_folder = normalizePath(BrDU_folder, mustWork = TRUE)
    )
  }
  
  pair_1 <- resolve_pair_folders(sample_1_folder)
  pair_2 <- resolve_pair_folders(sample_2_folder)
  
  resolve_origin_file <- function(origin_file, parent_1, parent_2) {
    shorthand <- NULL
    if (!is.null(origin_file) && length(origin_file) == 1L) {
      shorthand_value <- sub("\\.bed$", "", basename(origin_file), ignore.case = TRUE)
      if (shorthand_value %in% c("E_Ori", "L_Ori")) {
        shorthand <- shorthand_value
      }
    }
    
    if (!is.null(origin_file) && is.null(shorthand)) {
      return(normalizePath(origin_file, mustWork = TRUE))
    }
    
    origin_basename <- if (is.null(shorthand)) "E_Ori" else shorthand
    app_basename <- if (origin_basename == "E_Ori") "E_Rep" else "L_Rep"
    
    sample_candidate_1 <- file.path(
      parent_1, "Peaks", paste0(origin_basename, ".bed")
    )
    sample_candidate_2 <- file.path(
      parent_2, "Peaks", paste0(origin_basename, ".bed")
    )
    app_candidate <- file.path(
      "/Applications/ngsAnalyser.app/Contents/Resources/app",
      paste0(app_basename, ".bed")
    )
    
    if (file.exists(sample_candidate_1)) {
      if (file.exists(sample_candidate_2)) {
        hash_1 <- unname(tools::md5sum(sample_candidate_1))
        hash_2 <- unname(tools::md5sum(sample_candidate_2))
        if (!identical(hash_1, hash_2)) {
          warning(
            "The two sample-local Peaks/", origin_basename, ".bed files differ. ",
            "Using the file from sample_1_folder. Supply origin_file ",
            "explicitly if another origin set is intended."
          )
        }
      }
      return(normalizePath(sample_candidate_1, mustWork = TRUE))
    }
    if (file.exists(sample_candidate_2)) {
      return(normalizePath(sample_candidate_2, mustWork = TRUE))
    }
    if (file.exists(app_candidate)) {
      return(normalizePath(app_candidate, mustWork = TRUE))
    }
    
    stop(
      "No origin BED was found automatically. Supply origin_file explicitly, ",
      "or place ", origin_basename, ".bed in <sample parent>/Peaks/."
    )
  }
  
  origin_file <- resolve_origin_file(
    origin_file,
    pair_1$ChIP_folder,
    pair_2$ChIP_folder
  )
  
  OriginName <- basename(origin_file)
  OriginClass <- if (grepl("E_(Ori|Rep)", OriginName, ignore.case = TRUE)) {
    "early"
  } else if (grepl("L_(Ori|Rep)", OriginName, ignore.case = TRUE)) {
    "late"
  } else {
    "selected"
  }
  
  if (is.null(plot_title) || !nzchar(plot_title)) {
    SampleParts <- strsplit(c(pair_1$prefix, pair_2$prefix), "-", fixed = TRUE)
    
    if (all(lengths(SampleParts) >= 4L)) {
      Target <- SampleParts[[1]][1]
      Genotypes <- vapply(SampleParts, function(x) x[2], character(1))
      Conditions <- unique(vapply(SampleParts, function(x) x[3], character(1)))
      Years <- unique(vapply(SampleParts, function(x) x[4], character(1)))
      
      plot_title <- paste0(
        Target, " ", paste(Conditions, collapse = "/"),
        " ChIP and BrDU enrichment at ", OriginClass, " origins: ",
        paste(Genotypes, collapse = " vs "),
        " (", paste(Years, collapse = "/"), ")"
      )
    } else {
      plot_title <- paste0(
        "ChIP and BrDU enrichment at ", OriginClass, " origins: ",
        sample_1_label, " vs ", sample_2_label
      )
    }
  }
  
  if (is.null(output_pdf)) {
    common_parent <- dirname(pair_1$ChIP_folder)
    
    while (
      common_parent != dirname(common_parent) &&
        !startsWith(
          paste0(pair_2$ChIP_folder, .Platform$file.sep),
          paste0(common_parent, .Platform$file.sep)
        )
    ) {
      common_parent <- dirname(common_parent)
    }
    
    safe_1 <- gsub("[^A-Za-z0-9._-]+", "_", sample_1_label)
    safe_2 <- gsub("[^A-Za-z0-9._-]+", "_", sample_2_label)
    safe_origin <- tools::file_path_sans_ext(basename(origin_file))
    safe_origin <- gsub("[^A-Za-z0-9._-]+", "_", safe_origin)
    
    output_pdf <- file.path(
      common_parent,
      paste0(
        safe_1, "_vs_", safe_2, "_", safe_origin, "_",
        Alignment, "_", Signal,
        if (yLim == "per-plot") "_per_plot_ylim" else "",
        "_ChIP_BrDU_origin_distance_analysis.pdf"
      )
    )
  }
  output_pdf <- path.expand(output_pdf)
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  
  read_origins <- function(path) {
    x <- data.table::fread(
      path, header = "auto", select = 1:4,
      data.table = FALSE, showProgress = FALSE
    )
    if (ncol(x) != 4L) {
      stop("Origin BED must contain at least four columns: ", path)
    }
    names(x) <- c("chrom", "chromStart", "chromEnd", "oriName")
    x$chromStart <- as.numeric(x$chromStart)
    x$chromEnd <- as.numeric(x$chromEnd)
    x$oriCenter <- round((x$chromStart + x$chromEnd)/2)
    
    if (any(!is.finite(x$chromStart)) || any(!is.finite(x$chromEnd))) {
      stop("Origin file contains non-numeric coordinates.")
    }
    x$origin_index <- seq_len(nrow(x))
    fallback_names <- paste0("origin_", x$origin_index)
    usable_names <- as.character(x$oriName)
    missing_name <- is.na(usable_names) | trimws(usable_names) == ""
    usable_names[missing_name] <- fallback_names[missing_name]
    x$origin_id <- make.unique(usable_names)
    x
  }
  
  locate_ratio_file <- function(parent, Assay) {
    ratios_dir <- file.path(parent, RatioFolder)
    if (!dir.exists(ratios_dir)) {
      stop(RatioFolder, " directory not found: ", ratios_dir)
    }
    
    sample_base <- basename(parent)
    preferred <- file.path(
      ratios_dir,
      paste0(sample_base, "_", Assay, "_collapsed.bed")
    )
    if (file.exists(preferred)) return(preferred)
    
    candidates <- list.files(
      ratios_dir,
      pattern = paste0("_", Assay, "_collapsed\\.bed$"),
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(candidates) != 1L) {
      stop(
        "Expected exactly one collapsed ", Assay, " BED file in ",
        ratios_dir, "; found ", length(candidates), "."
      )
    }
    candidates
  }
  
  read_ratio_table <- function(path) {
    x <- data.table::fread(
      path, header = TRUE, data.table = TRUE,
      showProgress = FALSE
    )
    required <- c("chrom", "chromStart", "chromEnd", RatioColumn)
    missing_columns <- setdiff(required, names(x))
    if (length(missing_columns) > 0L) {
      stop(
        "Missing required column(s) in ", path, ": ",
        paste(missing_columns, collapse = ", ")
      )
    }
    x <- x[, ..required]
    data.table::setnames(x, RatioColumn, "ratio")
    if (any(!is.finite(x$chromStart)) || any(!is.finite(x$ratio))) {
      stop("Non-finite coordinate or ratio in: ", path)
    }
    data.table::setkey(x, chrom, chromStart)
    x
  }
  
  load_pair <- function(pair) {
    paths <- list(
      ChIP_collapsed = locate_ratio_file(pair$ChIP_folder, "ChIP"),
      BrDU_collapsed = locate_ratio_file(pair$BrDU_folder, "BrDU")
    )
    message("Loading ", pair$prefix, " collapsed ChIP and BrDU ",
            RatioColumn, " tables...")
    tables <- lapply(paths, read_ratio_table)
    
    list(paths = paths, tables = tables)
  }
  
  infer_step <- function(...) {
    tables <- list(...)
    steps <- unlist(lapply(tables, function(z) {
      starts <- z$chromStart[z$chrom == z$chrom[1]]
      d <- diff(starts)
      min(d[is.finite(d) & d > 0])
    }))
    if (length(steps) == 0L || any(!is.finite(steps)) ||
        length(unique(steps)) != 1L) {
      stop("Could not infer one common coverage step.")
    }
    as.numeric(steps[1])
  }
  
  extract_centered_matrix <- function(coverage, origins, window_bp, step) {
    expected_n <- max(1L, round(2*window_bp/step))
    ans <- matrix(
      0, nrow = expected_n, ncol = nrow(origins),
      dimnames = list(NULL, origins$origin_id)
    )
    
    for (i in seq_len(nrow(origins))) {
      center <- origins$oriCenter[i]
      region_start <- center - window_bp
      region_end <- center + window_bp
      signal <- coverage[
        chrom == origins$chrom[i] &
          chromStart >= region_start &
          chromStart <= region_end
      ]
      data.table::setorder(signal, chromStart)
      values <- signal$ratio
      values[!is.finite(values)] <- 0
      
      if (length(values) < expected_n) {
        missing_n <- expected_n - length(values)
        if (nrow(signal) == 0L) {
          values <- numeric(expected_n)
        } else {
          left_missing <- round((min(signal$chromStart) - region_start)/step)
          left_missing <- min(max(left_missing, 0L), missing_n)
          right_missing <- missing_n - left_missing
          values <- c(rep(0, left_missing), values, rep(0, right_missing))
        }
      }
      if (length(values) > expected_n) {
        values <- values[seq_len(expected_n)]
      }
      ans[, i] <- values
    }
    ans
  }
  
  make_assay_result <- function(coverage, origins, window_bp, step) {
    per_origin_signal <- extract_centered_matrix(
      coverage, origins, window_bp, step
    )
    collapsed_profile <- apply(
      per_origin_signal, 1, stats::median, na.rm = TRUE
    )
    smooth_collapsed <- stats::smooth.spline(
      seq_along(collapsed_profile), collapsed_profile, spar = 0.5
    )$y
    
    list(profile = smooth_collapsed, per_origin_signal = per_origin_signal)
  }
  
  interval_summary <- function(signal_matrix, distance_bp, breaks_bp) {
    interval_labels <- paste0(
      format(breaks_bp[-length(breaks_bp)]/1000, trim = TRUE),
      "-",
      format(breaks_bp[-1]/1000, trim = TRUE),
      " kb"
    )
    out <- matrix(
      NA_real_, nrow = ncol(signal_matrix), ncol = length(interval_labels),
      dimnames = list(colnames(signal_matrix), interval_labels)
    )
    
    absolute_distance <- abs(distance_bp)
    for (j in seq_along(interval_labels)) {
      include <- absolute_distance >= breaks_bp[j] &
        if (j == length(interval_labels)) {
          absolute_distance <= breaks_bp[j + 1L]
        } else {
          absolute_distance < breaks_bp[j + 1L]
        }
      out[, j] <- colMeans(signal_matrix[include, , drop = FALSE], na.rm = TRUE)
    }
    out
  }
  
  origins <- read_origins(origin_file)
  if (!is.null(expected_origin_count) && nrow(origins) != expected_origin_count) {
    stop(
      "Expected ", expected_origin_count, " origins, but origin file contains ",
      nrow(origins), "."
    )
  }
  
  sample_1 <- load_pair(pair_1)
  sample_2 <- load_pair(pair_2)
  step <- infer_step(
    sample_1$tables$ChIP_collapsed,
    sample_1$tables$BrDU_collapsed,
    sample_2$tables$ChIP_collapsed,
    sample_2$tables$BrDU_collapsed
  )
  
  window_bp <- window_kb*1000
  breaks_bp <- distance_breaks_kb*1000
  expected_n <- max(1L, round(2*window_bp/step))
  distance_bp <- seq(from = -window_bp, by = step, length.out = expected_n)
  x_kb <- distance_bp/1000
  profile_limit_kb <- max(distance_breaks_kb)
  
  result_1 <- list(
    ChIP = make_assay_result(
      sample_1$tables$ChIP_collapsed, origins, window_bp, step
    ),
    BrDU = make_assay_result(
      sample_1$tables$BrDU_collapsed, origins, window_bp, step
    )
  )
  result_2 <- list(
    ChIP = make_assay_result(
      sample_2$tables$ChIP_collapsed, origins, window_bp, step
    ),
    BrDU = make_assay_result(
      sample_2$tables$BrDU_collapsed, origins, window_bp, step
    )
  )
  
  boxes_1 <- interval_summary(
    result_1$ChIP$per_origin_signal, distance_bp, breaks_bp
  )
  boxes_2 <- interval_summary(
    result_2$ChIP$per_origin_signal, distance_bp, breaks_bp
  )
  BrDU_boxes_1 <- interval_summary(
    result_1$BrDU$per_origin_signal, distance_bp, breaks_bp
  )
  BrDU_boxes_2 <- interval_summary(
    result_2$BrDU$per_origin_signal, distance_bp, breaks_bp
  )
  
  ProfileYlim <- function(Profile) {
    Values <- Profile[is.finite(Profile)]
    if (length(Values) == 0L) return(c(-1, 1))
    YRange <- range(Values)
    YPad <- diff(YRange)*0.08
    if (!is.finite(YPad) || YPad == 0) {
      YPad <- max(0.5, abs(YRange[2])*0.08)
    }
    YRange + c(-YPad, YPad)
  }
  
  AlignYlimsAtOne <- function(ChIPLimits, BrDULimits) {
    ChIPLimits[1] <- min(ChIPLimits[1], 1)
    ChIPLimits[2] <- max(ChIPLimits[2], 1)
    BrDULimits[1] <- min(BrDULimits[1], 1)
    BrDULimits[2] <- max(BrDULimits[2], 1)
    
    ChIPFraction <- (1-ChIPLimits[1])/diff(ChIPLimits)
    BrDUFraction <- (1-BrDULimits[1])/diff(BrDULimits)
    AlignedFraction <- max(ChIPFraction, BrDUFraction)
    
    if (AlignedFraction > 0 && AlignedFraction < 1) {
      ChIPLimits[1] <- (1-AlignedFraction*ChIPLimits[2])/(1-AlignedFraction)
      BrDULimits[1] <- (1-AlignedFraction*BrDULimits[2])/(1-AlignedFraction)
    }
    
    list(ChIP = ChIPLimits, BrDU = BrDULimits)
  }
  
  FixedYlims <- AlignYlimsAtOne(
    ProfileYlim(c(result_1$ChIP$profile, result_2$ChIP$profile)),
    ProfileYlim(c(result_1$BrDU$profile, result_2$BrDU$profile))
  )
  
  PlotOverlay <- function(Result, SampleLabel, ProfileColour) {
    Limits <- if (yLim == "fixed") {
      FixedYlims
    } else {
      AlignYlimsAtOne(
        ProfileYlim(Result$ChIP$profile),
        ProfileYlim(Result$BrDU$profile)
      )
    }
    BrDUAxisTicks <- pretty(Limits$BrDU)
    BrDUAxisTicks <- BrDUAxisTicks[
      BrDUAxisTicks >= 1 & BrDUAxisTicks <= Limits$BrDU[2]
    ]
    BrDUAxisTicks <- sort(unique(c(1, BrDUAxisTicks)))
    BrDUAxisLabels <- format(BrDUAxisTicks, trim = TRUE, nsmall = 1)
    
    graphics::plot(
      x_kb, Result$ChIP$profile,
      type = "l", lwd = 2.2, lty = 1, col = ProfileColour,
      xlim = c(-profile_limit_kb, profile_limit_kb), ylim = Limits$ChIP,
      xlab = "Distance from origin centre (kb)", ylab = ChIPYAxisLabel,
      main = SampleLabel, xaxt = "n", bty = "l", las = 1,
      xaxs = "i", yaxs = "i", cex.axis = 1.05,
      cex.lab = 1.08, cex.main = 1.12
    )
    ProfileAxisTicks <- pretty(c(-profile_limit_kb, profile_limit_kb), n = 5)
    ProfileAxisTicks <- ProfileAxisTicks[
      ProfileAxisTicks >= -profile_limit_kb & ProfileAxisTicks <= profile_limit_kb
    ]
    ProfileAxisTicks <- sort(unique(c(
      -profile_limit_kb, ProfileAxisTicks, 0, profile_limit_kb
    )))
    graphics::axis(1, at = ProfileAxisTicks, labels = ProfileAxisTicks)
    graphics::abline(v = 0, col = "grey75", lty = 2)
    
    graphics::par(new = TRUE)
    graphics::plot(
      x_kb, Result$BrDU$profile,
      type = "l", lwd = 2.2, lty = 3, col = ProfileColour,
      xlim = c(-profile_limit_kb, profile_limit_kb), ylim = Limits$BrDU,
      axes = FALSE, xlab = "", ylab = "", bty = "n",
      xaxs = "i", yaxs = "i"
    )
    graphics::axis(
      4, at = BrDUAxisTicks, labels = BrDUAxisLabels, las = 1,
      col = "transparent", col.axis = "black", col.ticks = "black",
      lty = 1, lwd = 0, lwd.ticks = 1.3, tck = -0.018,
      cex.axis = 1.05
    )
    AxisLimits <- graphics::par("usr")
    graphics::segments(
      AxisLimits[2], AxisLimits[3], AxisLimits[2], AxisLimits[4],
      col = "black", lty = 3, lwd = 1.6, xpd = NA
    )
    graphics::mtext(
      BrDUYAxisLabel, side = 4, line = 2.8,
      col = "black", cex = 1.08
    )
    graphics::legend(
      "topright", legend = c("ChIP", "BrDU"),
      col = c(ProfileColour, ProfileColour), lty = c(1, 3),
      lwd = 2.2, bty = "n", cex = 0.85
    )
  }
  
  interval_labels <- colnames(boxes_1)
  NumberOfIntervals <- length(interval_labels)
  plot_box_data <- c(
    lapply(seq_len(NumberOfIntervals), function(i) boxes_1[, i]),
    lapply(seq_len(NumberOfIntervals), function(i) boxes_2[, i])
  )
  sample_1_positions <- seq_len(NumberOfIntervals)
  sample_2_positions <- seq_len(NumberOfIntervals) + NumberOfIntervals + 2
  box_positions <- c(sample_1_positions, sample_2_positions)
  group_midpoints <- c(mean(sample_1_positions), mean(sample_2_positions))
  box_stats <- lapply(plot_box_data, grDevices::boxplot.stats)
  box_ylim <- range(
    unlist(lapply(box_stats, function(z) z$stats)), finite = TRUE
  )
  box_pad <- diff(box_ylim)*0.08
  if (!is.finite(box_pad) || box_pad == 0) box_pad <- 0.5
  
  grDevices::pdf(
    output_pdf, width = pdf_width, height = pdf_height,
    onefile = TRUE, useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  
  graphics::layout(
    matrix(c(0, 1, 0, 2, 0,
             0, 3, 3, 3, 0), nrow = 2, ncol = 5, byrow = TRUE),
    widths = c(1.2, 2.4, 0.8, 2.4, 1.2),
    heights = c(1, 1.3)
  )
  graphics::par(oma = c(0, 0, 3.5, 0))
  
  graphics::par(mar = c(4.3, 4.7, 3.2, 4.7))
  PlotOverlay(result_1, sample_1_label, profile_colours[1])
  
  graphics::par(mar = c(4.3, 4.7, 3.2, 4.7))
  PlotOverlay(result_2, sample_2_label, profile_colours[2])
  
  graphics::par(mar = c(6.3, 5, 3.2, 1.2))
  graphics::boxplot(
    plot_box_data, at = box_positions,
    col = c(rep(box_colours[1], NumberOfIntervals),
            rep(box_colours[2], NumberOfIntervals)),
    border = c(rep(profile_colours[1], NumberOfIntervals),
               rep(profile_colours[2], NumberOfIntervals)),
    outline = FALSE, xaxt = "n",
    ylim = box_ylim + c(-box_pad, box_pad),
    ylab = ChIPYAxisLabel,
    main = "Per-origin ChIP signal by absolute distance interval",
    bty = "n", las = 1, cex.axis = 1.12,
    cex.lab = 1.18, cex.main = 1.12
  )
  
  set.seed(123)
  for (i in seq_along(plot_box_data)) {
    visible <- plot_box_data[[i]] >= box_stats[[i]]$stats[1] &
      plot_box_data[[i]] <= box_stats[[i]]$stats[5]
    visible_values <- plot_box_data[[i]][visible]
    SampleIndex <- if (i <= NumberOfIntervals) 1L else 2L
    graphics::points(
      jitter(rep(box_positions[i], length(visible_values)), amount = 0.16),
      visible_values, pch = 16, cex = point_cex,
      col = grDevices::adjustcolor(
        profile_colours[SampleIndex], alpha.f = point_alpha
      )
    )
  }
  
  graphics::axis(
    1, at = box_positions,
    labels = rep(interval_labels, 2), tick = FALSE,
    line = 0.15, cex.axis = 0.82, las = 1
  )
  graphics::axis(
    1, at = group_midpoints,
    labels = c(sample_1_label, sample_2_label),
    tick = FALSE, line = 2, cex.axis = 1.08
  )
  graphics::mtext(
    "Absolute distance from origin centre",
    side = 1, line = 4.1, cex = 1.08
  )
  graphics::mtext(
    plot_title, side = 3, line = 1.3,
    outer = TRUE, font = 2, cex = 1.25
  )
  
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  
  MakeLongTable <- function(Assay, Values_1, Values_2) {
    do.call(
      rbind,
      lapply(seq_along(interval_labels), function(i) {
        rbind(
          data.frame(
            origin_index = origins$origin_index,
            oriName = origins$oriName, chrom = origins$chrom,
            oriCenter = origins$oriCenter, sample = sample_1_label,
            assay = Assay, interval = interval_labels[i],
            value = Values_1[, i], stringsAsFactors = FALSE
          ),
          data.frame(
            origin_index = origins$origin_index,
            oriName = origins$oriName, chrom = origins$chrom,
            oriCenter = origins$oriCenter, sample = sample_2_label,
            assay = Assay, interval = interval_labels[i],
            value = Values_2[, i], stringsAsFactors = FALSE
          )
        )
      })
    )
  }
  
  long_table <- rbind(
    MakeLongTable("ChIP", boxes_1, boxes_2),
    MakeLongTable("BrDU", BrDU_boxes_1, BrDU_boxes_2)
  )
  
  message("Saved ChIP-BrDU comparison PDF: ", output_pdf)
  
  legend_note <- paste(
    "Average profiles show collapsed ChIP as a solid line against the left",
    "axis and sister BrDU as a dotted line against the right axis. Value 1",
    "is aligned between the two axes within each sample panel. The profile",
    "x-axis ends at the largest distance_breaks_kb value. Box plots contain",
    "per-origin ChIP enrichment and are grouped as all sample 1 intervals",
    "followed by all sample 2 intervals; no pairwise p-values are displayed."
  )
  
  result <- list(
    output_pdf = output_pdf,
    alignment = Alignment,
    signal = Signal,
    yLim = yLim,
    signal_source = paste0(RatioFolder, " collapsed tables: ", RatioColumn),
    coverage_step_bp = step,
    origins = origins,
    distance_bp = distance_bp,
    profile_limit_kb = profile_limit_kb,
    interval_breaks_kb = distance_breaks_kb,
    legend_note = legend_note,
    sample_1 = list(
      label = sample_1_label,
      ChIP = list(
        profile = result_1$ChIP$profile,
        interval_values = boxes_1,
        file = sample_1$paths$ChIP_collapsed
      ),
      BrDU = list(
        profile = result_1$BrDU$profile,
        interval_values = BrDU_boxes_1,
        file = sample_1$paths$BrDU_collapsed
      )
    ),
    sample_2 = list(
      label = sample_2_label,
      ChIP = list(
        profile = result_2$ChIP$profile,
        interval_values = boxes_2,
        file = sample_2$paths$ChIP_collapsed
      ),
      BrDU = list(
        profile = result_2$BrDU$profile,
        interval_values = BrDU_boxes_2,
        file = sample_2$paths$BrDU_collapsed
      )
    ),
    interval_data = long_table
  )
  
  if (isTRUE(return_data)) invisible(result) else invisible(output_pdf)
}


origin_distance_pair_analysis_old_approach <- function(
    sample_1_folder,
    sample_2_folder,
    origin_file = NULL,
    Alignment = "generic",
    yLim = "fixed",
    sample_1_label = basename(normalizePath(sample_1_folder, mustWork = FALSE)),
    sample_2_label = basename(normalizePath(sample_2_folder, mustWork = FALSE)),
    distance_breaks_kb = c(0, 0.5, 1.5, 2.5, 4, 8),
    window_kb = 12,
    output_pdf = NULL,
    plot_title = NULL,
    profile_colours = c("#2166AC", "#B2182B"),
    box_colours = c("#92C5DE", "#F4A582"),
    point_alpha = 0.32,
    point_cex = 0.85,
    pdf_width = 15,
    pdf_height = 8,
    expected_origin_count = NULL,
    return_data = TRUE) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  
  Alignment <- match.arg(Alignment, c("generic", "malign"))
  yLim <- match.arg(yLim, c("fixed", "per-plot"))
  
  RatioFolder <- if (Alignment == "generic") "Ratios" else "Ratios_ma"
  ChIPYAxisLabel <- "ChIP/Input enrichment"
  BrDUYAxisLabel <- "BrDU/Input enrichment"
  
  if (length(distance_breaks_kb) < 2L ||
      any(!is.finite(distance_breaks_kb)) ||
      any(diff(distance_breaks_kb) <= 0) ||
      distance_breaks_kb[1] != 0) {
    stop("distance_breaks_kb must be strictly increasing and begin at zero.")
  }
  
  if (!is.finite(window_kb) || window_kb <= 0 ||
      window_kb < max(distance_breaks_kb)) {
    stop("window_kb must be positive and cover all distance intervals.")
  }
  
  sample_1_folder <- normalizePath(sample_1_folder, mustWork = TRUE)
  sample_2_folder <- normalizePath(sample_2_folder, mustWork = TRUE)
  
  resolve_pair_folders <- function(sample_folder) {
    sample_name <- basename(sample_folder)
    
    if (!grepl("-(ChIP|BrDU)$", sample_name)) {
      stop("Sample folder must end in -ChIP or -BrDU: ", sample_folder)
    }
    
    sample_prefix <- sub("-(ChIP|BrDU)$", "", sample_name)
    sample_parent <- dirname(sample_folder)
    ChIP_folder <- file.path(sample_parent, paste0(sample_prefix, "-ChIP"))
    BrDU_folder <- file.path(sample_parent, paste0(sample_prefix, "-BrDU"))
    
    if (!dir.exists(ChIP_folder)) {
      stop("Sister ChIP folder is missing: ", ChIP_folder)
    }
    if (!dir.exists(BrDU_folder)) {
      stop("Sister BrDU folder is missing: ", BrDU_folder)
    }
    
    list(
      prefix = sample_prefix,
      ChIP_folder = normalizePath(ChIP_folder, mustWork = TRUE),
      BrDU_folder = normalizePath(BrDU_folder, mustWork = TRUE)
    )
  }
  
  pair_1 <- resolve_pair_folders(sample_1_folder)
  pair_2 <- resolve_pair_folders(sample_2_folder)
  
  resolve_origin_file <- function(origin_file, parent_1, parent_2) {
    shorthand <- NULL
    if (!is.null(origin_file) && length(origin_file) == 1L) {
      shorthand_value <- sub("\\.bed$", "", basename(origin_file), ignore.case = TRUE)
      if (shorthand_value %in% c("E_Ori", "L_Ori")) {
        shorthand <- shorthand_value
      }
    }
    
    if (!is.null(origin_file) && is.null(shorthand)) {
      return(normalizePath(origin_file, mustWork = TRUE))
    }
    
    origin_basename <- if (is.null(shorthand)) "E_Ori" else shorthand
    app_basename <- if (origin_basename == "E_Ori") "E_Rep" else "L_Rep"
    
    sample_candidate_1 <- file.path(
      parent_1, "Peaks", paste0(origin_basename, ".bed")
    )
    sample_candidate_2 <- file.path(
      parent_2, "Peaks", paste0(origin_basename, ".bed")
    )
    app_candidate <- file.path(
      "/Applications/ngsAnalyser.app/Contents/Resources/app",
      paste0(app_basename, ".bed")
    )
    
    if (file.exists(sample_candidate_1)) {
      if (file.exists(sample_candidate_2)) {
        hash_1 <- unname(tools::md5sum(sample_candidate_1))
        hash_2 <- unname(tools::md5sum(sample_candidate_2))
        if (!identical(hash_1, hash_2)) {
          warning(
            "The two sample-local Peaks/", origin_basename, ".bed files differ. ",
            "Using the file from sample_1_folder. Supply origin_file ",
            "explicitly if another origin set is intended."
          )
        }
      }
      return(normalizePath(sample_candidate_1, mustWork = TRUE))
    }
    if (file.exists(sample_candidate_2)) {
      return(normalizePath(sample_candidate_2, mustWork = TRUE))
    }
    if (file.exists(app_candidate)) {
      return(normalizePath(app_candidate, mustWork = TRUE))
    }
    
    stop(
      "No origin BED was found automatically. Supply origin_file explicitly, ",
      "or place ", origin_basename, ".bed in <sample parent>/Peaks/."
    )
  }
  
  origin_file <- resolve_origin_file(
    origin_file,
    pair_1$ChIP_folder,
    pair_2$ChIP_folder
  )
  
  OriginName <- basename(origin_file)
  OriginClass <- if (grepl("E_(Ori|Rep)", OriginName, ignore.case = TRUE)) {
    "early"
  } else if (grepl("L_(Ori|Rep)", OriginName, ignore.case = TRUE)) {
    "late"
  } else {
    "selected"
  }
  
  if (is.null(plot_title) || !nzchar(plot_title)) {
    SampleParts <- strsplit(c(pair_1$prefix, pair_2$prefix), "-", fixed = TRUE)
    
    if (all(lengths(SampleParts) >= 4L)) {
      Target <- SampleParts[[1]][1]
      Genotypes <- vapply(SampleParts, function(x) x[2], character(1))
      Conditions <- unique(vapply(SampleParts, function(x) x[3], character(1)))
      Years <- unique(vapply(SampleParts, function(x) x[4], character(1)))
      
      plot_title <- paste0(
        Target, " ", paste(Conditions, collapse = "/"),
        " ChIP and BrDU enrichment at ", OriginClass, " origins: ",
        paste(Genotypes, collapse = " vs "),
        " (", paste(Years, collapse = "/"), ")"
      )
    } else {
      plot_title <- paste0(
        "ChIP and BrDU enrichment at ", OriginClass, " origins: ",
        sample_1_label, " vs ", sample_2_label
      )
    }
  }
  
  if (is.null(output_pdf)) {
    common_parent <- dirname(pair_1$ChIP_folder)
    
    while (
      common_parent != dirname(common_parent) &&
        !startsWith(
          paste0(pair_2$ChIP_folder, .Platform$file.sep),
          paste0(common_parent, .Platform$file.sep)
        )
    ) {
      common_parent <- dirname(common_parent)
    }
    
    safe_1 <- gsub("[^A-Za-z0-9._-]+", "_", sample_1_label)
    safe_2 <- gsub("[^A-Za-z0-9._-]+", "_", sample_2_label)
    safe_origin <- tools::file_path_sans_ext(basename(origin_file))
    safe_origin <- gsub("[^A-Za-z0-9._-]+", "_", safe_origin)
    
    output_pdf <- file.path(
      common_parent,
      paste0(
        safe_1, "_vs_", safe_2, "_", safe_origin, "_",
        Alignment,
        if (yLim == "per-plot") "_per_plot_ylim" else "",
        "_ChIP_BrDU_origin_distance_analysis_old_approach.pdf"
      )
    )
  }
  output_pdf <- path.expand(output_pdf)
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  
  read_origins <- function(path) {
    x <- data.table::fread(
      path, header = "auto", select = 1:4,
      data.table = FALSE, showProgress = FALSE
    )
    if (ncol(x) != 4L) {
      stop("Origin BED must contain at least four columns: ", path)
    }
    names(x) <- c("chrom", "chromStart", "chromEnd", "oriName")
    x$chromStart <- as.numeric(x$chromStart)
    x$chromEnd <- as.numeric(x$chromEnd)
    x$oriCenter <- round((x$chromStart + x$chromEnd)/2)
    
    if (any(!is.finite(x$chromStart)) || any(!is.finite(x$chromEnd))) {
      stop("Origin file contains non-numeric coordinates.")
    }
    x$origin_index <- seq_len(nrow(x))
    fallback_names <- paste0("origin_", x$origin_index)
    usable_names <- as.character(x$oriName)
    missing_name <- is.na(usable_names) | trimws(usable_names) == ""
    usable_names[missing_name] <- fallback_names[missing_name]
    x$origin_id <- make.unique(usable_names)
    x
  }
  
  locate_ratio_file <- function(parent, Assay, Strand) {
    ratios_dir <- file.path(parent, RatioFolder)
    if (!dir.exists(ratios_dir)) {
      stop(RatioFolder, " directory not found: ", ratios_dir)
    }
    
    sample_base <- basename(parent)
    preferred <- file.path(
      ratios_dir,
      paste0(sample_base, "_", Assay, "_", Strand, ".bed")
    )
    if (file.exists(preferred)) return(preferred)
    
    candidates <- list.files(
      ratios_dir,
      pattern = paste0("_", Assay, "_", Strand, "\\.bed$"),
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(candidates) != 1L) {
      stop(
        "Expected exactly one ", Assay, " ", Strand, " BED file in ",
        ratios_dir, "; found ", length(candidates), "."
      )
    }
    candidates
  }
  
  read_ratio_table <- function(path) {
    x <- data.table::fread(
      path, header = TRUE, data.table = TRUE,
      showProgress = FALSE
    )
    required <- c(
      "chrom", "chromStart", "chromEnd", "ip.score", "in.score"
    )
    missing_columns <- setdiff(required, names(x))
    if (length(missing_columns) > 0L) {
      stop(
        "Missing required column(s) in ", path, ": ",
        paste(missing_columns, collapse = ", ")
      )
    }
    x <- x[, ..required]
    if (any(!is.finite(x$chromStart)) ||
        any(!is.finite(x$ip.score)) ||
        any(!is.finite(x$in.score))) {
      stop("Non-finite coordinate, IP score or Input score in: ", path)
    }
    data.table::setkey(x, chrom, chromStart)
    x
  }
  
  load_pair <- function(pair) {
    paths <- list(
      ChIP_watson = locate_ratio_file(pair$ChIP_folder, "ChIP", "watson"),
      ChIP_crick = locate_ratio_file(pair$ChIP_folder, "ChIP", "crick"),
      BrDU_watson = locate_ratio_file(pair$BrDU_folder, "BrDU", "watson"),
      BrDU_crick = locate_ratio_file(pair$BrDU_folder, "BrDU", "crick")
    )
    message(
      "Loading ", pair$prefix,
      " strand-separated ChIP and BrDU ip.score/in.score tables..."
    )
    tables <- lapply(paths, read_ratio_table)
    
    list(paths = paths, tables = tables)
  }
  
  infer_step <- function(...) {
    tables <- list(...)
    steps <- unlist(lapply(tables, function(z) {
      starts <- z$chromStart[z$chrom == z$chrom[1]]
      d <- diff(starts)
      min(d[is.finite(d) & d > 0])
    }))
    if (length(steps) == 0L || any(!is.finite(steps)) ||
        length(unique(steps)) != 1L) {
      stop("Could not infer one common coverage step.")
    }
    as.numeric(steps[1])
  }
  
  extract_centered_matrix <- function(coverage, origins, window_bp, step,
                                      ValueColumn) {
    expected_n <- max(1L, round(2*window_bp/step))
    ans <- matrix(
      0, nrow = expected_n, ncol = nrow(origins),
      dimnames = list(NULL, origins$origin_id)
    )
    
    for (i in seq_len(nrow(origins))) {
      center <- origins$oriCenter[i]
      region_start <- center - window_bp
      region_end <- center + window_bp
      signal <- coverage[
        chrom == origins$chrom[i] &
          chromStart >= region_start &
          chromStart <= region_end
      ]
      data.table::setorder(signal, chromStart)
      values <- signal[[ValueColumn]]
      values[!is.finite(values)] <- 0
      
      if (length(values) < expected_n) {
        missing_n <- expected_n - length(values)
        if (nrow(signal) == 0L) {
          values <- numeric(expected_n)
        } else {
          left_missing <- round((min(signal$chromStart) - region_start)/step)
          left_missing <- min(max(left_missing, 0L), missing_n)
          right_missing <- missing_n - left_missing
          values <- c(rep(0, left_missing), values, rep(0, right_missing))
        }
      }
      if (length(values) > expected_n) {
        values <- values[seq_len(expected_n)]
      }
      ans[, i] <- values
    }
    ans
  }
  
  make_assay_result <- function(WatsonCoverage, CrickCoverage,
                                 origins, window_bp, step) {
    WatsonIP <- extract_centered_matrix(
      WatsonCoverage, origins, window_bp, step, "ip.score"
    )
    WatsonInput <- extract_centered_matrix(
      WatsonCoverage, origins, window_bp, step, "in.score"
    )
    CrickIP <- extract_centered_matrix(
      CrickCoverage, origins, window_bp, step, "ip.score"
    )
    CrickInput <- extract_centered_matrix(
      CrickCoverage, origins, window_bp, step, "in.score"
    )
    
    WatsonInNorm <- apply(WatsonIP, 1, stats::median, na.rm = TRUE)/
      apply(WatsonInput, 1, stats::median, na.rm = TRUE)
    CrickInNorm <- apply(CrickIP, 1, stats::median, na.rm = TRUE)/
      apply(CrickInput, 1, stats::median, na.rm = TRUE)
    
    WatsonInNorm[!is.finite(WatsonInNorm)] <- 0
    CrickInNorm[!is.finite(CrickInNorm)] <- 0
    
    Watson <- stats::smooth.spline(
      seq_along(WatsonInNorm), WatsonInNorm, spar = 0.5
    )$y
    Crick <- stats::smooth.spline(
      seq_along(CrickInNorm), CrickInNorm, spar = 0.5
    )$y
    
    list(
      profile = Watson + Crick,
      WatsonIP = WatsonIP,
      WatsonInput = WatsonInput,
      CrickIP = CrickIP,
      CrickInput = CrickInput
    )
  }
  
  interval_summary <- function(AssayResult, distance_bp, breaks_bp) {
    interval_labels <- paste0(
      format(breaks_bp[-length(breaks_bp)]/1000, trim = TRUE),
      "-",
      format(breaks_bp[-1]/1000, trim = TRUE),
      " kb"
    )
    out <- matrix(
      NA_real_, nrow = ncol(AssayResult$WatsonIP),
      ncol = length(interval_labels),
      dimnames = list(colnames(AssayResult$WatsonIP), interval_labels)
    )
    
    absolute_distance <- abs(distance_bp)
    for (j in seq_along(interval_labels)) {
      include <- absolute_distance >= breaks_bp[j] &
        if (j == length(interval_labels)) {
          absolute_distance <= breaks_bp[j + 1L]
        } else {
          absolute_distance < breaks_bp[j + 1L]
        }
      
      WatsonIP <- colMeans(
        AssayResult$WatsonIP[include, , drop = FALSE], na.rm = TRUE
      )
      WatsonInput <- colMeans(
        AssayResult$WatsonInput[include, , drop = FALSE], na.rm = TRUE
      )
      CrickIP <- colMeans(
        AssayResult$CrickIP[include, , drop = FALSE], na.rm = TRUE
      )
      CrickInput <- colMeans(
        AssayResult$CrickInput[include, , drop = FALSE], na.rm = TRUE
      )
      
      IntervalValues <- WatsonIP/WatsonInput + CrickIP/CrickInput
      IntervalValues[!is.finite(IntervalValues)] <- 0
      out[, j] <- IntervalValues
    }
    out
  }
  
  origins <- read_origins(origin_file)
  if (!is.null(expected_origin_count) && nrow(origins) != expected_origin_count) {
    stop(
      "Expected ", expected_origin_count, " origins, but origin file contains ",
      nrow(origins), "."
    )
  }
  
  sample_1 <- load_pair(pair_1)
  sample_2 <- load_pair(pair_2)
  step <- infer_step(
    sample_1$tables$ChIP_watson,
    sample_1$tables$ChIP_crick,
    sample_1$tables$BrDU_watson,
    sample_1$tables$BrDU_crick,
    sample_2$tables$ChIP_watson,
    sample_2$tables$ChIP_crick,
    sample_2$tables$BrDU_watson,
    sample_2$tables$BrDU_crick
  )
  
  window_bp <- window_kb*1000
  breaks_bp <- distance_breaks_kb*1000
  expected_n <- max(1L, round(2*window_bp/step))
  distance_bp <- seq(from = -window_bp, by = step, length.out = expected_n)
  x_kb <- distance_bp/1000
  profile_limit_kb <- max(distance_breaks_kb)
  
  result_1 <- list(
    ChIP = make_assay_result(
      sample_1$tables$ChIP_watson, sample_1$tables$ChIP_crick,
      origins, window_bp, step
    ),
    BrDU = make_assay_result(
      sample_1$tables$BrDU_watson, sample_1$tables$BrDU_crick,
      origins, window_bp, step
    )
  )
  result_2 <- list(
    ChIP = make_assay_result(
      sample_2$tables$ChIP_watson, sample_2$tables$ChIP_crick,
      origins, window_bp, step
    ),
    BrDU = make_assay_result(
      sample_2$tables$BrDU_watson, sample_2$tables$BrDU_crick,
      origins, window_bp, step
    )
  )
  
  boxes_1 <- interval_summary(
    result_1$ChIP, distance_bp, breaks_bp
  )
  boxes_2 <- interval_summary(
    result_2$ChIP, distance_bp, breaks_bp
  )
  BrDU_boxes_1 <- interval_summary(
    result_1$BrDU, distance_bp, breaks_bp
  )
  BrDU_boxes_2 <- interval_summary(
    result_2$BrDU, distance_bp, breaks_bp
  )
  
  ProfileYlim <- function(Profile) {
    Values <- Profile[is.finite(Profile)]
    if (length(Values) == 0L) return(c(-1, 1))
    YRange <- range(Values)
    YPad <- diff(YRange)*0.08
    if (!is.finite(YPad) || YPad == 0) {
      YPad <- max(0.5, abs(YRange[2])*0.08)
    }
    YRange + c(-YPad, YPad)
  }
  
  AlignYlimsAtOne <- function(ChIPLimits, BrDULimits) {
    ChIPLimits[1] <- min(ChIPLimits[1], 1)
    ChIPLimits[2] <- max(ChIPLimits[2], 1)
    BrDULimits[1] <- min(BrDULimits[1], 1)
    BrDULimits[2] <- max(BrDULimits[2], 1)
    
    ChIPFraction <- (1-ChIPLimits[1])/diff(ChIPLimits)
    BrDUFraction <- (1-BrDULimits[1])/diff(BrDULimits)
    AlignedFraction <- max(ChIPFraction, BrDUFraction)
    
    if (AlignedFraction > 0 && AlignedFraction < 1) {
      ChIPLimits[1] <- (1-AlignedFraction*ChIPLimits[2])/(1-AlignedFraction)
      BrDULimits[1] <- (1-AlignedFraction*BrDULimits[2])/(1-AlignedFraction)
    }
    
    list(ChIP = ChIPLimits, BrDU = BrDULimits)
  }
  
  FixedYlims <- AlignYlimsAtOne(
    ProfileYlim(c(result_1$ChIP$profile, result_2$ChIP$profile)),
    ProfileYlim(c(result_1$BrDU$profile, result_2$BrDU$profile))
  )
  
  PlotOverlay <- function(Result, SampleLabel, ProfileColour) {
    Limits <- if (yLim == "fixed") {
      FixedYlims
    } else {
      AlignYlimsAtOne(
        ProfileYlim(Result$ChIP$profile),
        ProfileYlim(Result$BrDU$profile)
      )
    }
    BrDUAxisTicks <- pretty(Limits$BrDU)
    BrDUAxisTicks <- BrDUAxisTicks[
      BrDUAxisTicks >= 1 & BrDUAxisTicks <= Limits$BrDU[2]
    ]
    BrDUAxisTicks <- sort(unique(c(1, BrDUAxisTicks)))
    BrDUAxisLabels <- format(BrDUAxisTicks, trim = TRUE, nsmall = 1)
    
    graphics::plot(
      x_kb, Result$ChIP$profile,
      type = "l", lwd = 2.2, lty = 1, col = ProfileColour,
      xlim = c(-profile_limit_kb, profile_limit_kb), ylim = Limits$ChIP,
      xlab = "Distance from origin centre (kb)", ylab = ChIPYAxisLabel,
      main = SampleLabel, xaxt = "n", bty = "l", las = 1,
      xaxs = "i", yaxs = "i", cex.axis = 1.05,
      cex.lab = 1.08, cex.main = 1.12
    )
    ProfileAxisTicks <- pretty(c(-profile_limit_kb, profile_limit_kb), n = 5)
    ProfileAxisTicks <- ProfileAxisTicks[
      ProfileAxisTicks >= -profile_limit_kb & ProfileAxisTicks <= profile_limit_kb
    ]
    ProfileAxisTicks <- sort(unique(c(
      -profile_limit_kb, ProfileAxisTicks, 0, profile_limit_kb
    )))
    graphics::axis(1, at = ProfileAxisTicks, labels = ProfileAxisTicks)
    graphics::abline(v = 0, col = "grey75", lty = 2)
    
    graphics::par(new = TRUE)
    graphics::plot(
      x_kb, Result$BrDU$profile,
      type = "l", lwd = 2.2, lty = 3, col = ProfileColour,
      xlim = c(-profile_limit_kb, profile_limit_kb), ylim = Limits$BrDU,
      axes = FALSE, xlab = "", ylab = "", bty = "n",
      xaxs = "i", yaxs = "i"
    )
    graphics::axis(
      4, at = BrDUAxisTicks, labels = BrDUAxisLabels, las = 1,
      col = "transparent", col.axis = "black", col.ticks = "black",
      lty = 1, lwd = 0, lwd.ticks = 1.3, tck = -0.018,
      cex.axis = 1.05
    )
    AxisLimits <- graphics::par("usr")
    graphics::segments(
      AxisLimits[2], AxisLimits[3], AxisLimits[2], AxisLimits[4],
      col = "black", lty = 3, lwd = 1.6, xpd = NA
    )
    graphics::mtext(
      BrDUYAxisLabel, side = 4, line = 2.8,
      col = "black", cex = 1.08
    )
    graphics::legend(
      "topright", legend = c("ChIP", "BrDU"),
      col = c(ProfileColour, ProfileColour), lty = c(1, 3),
      lwd = 2.2, bty = "n", cex = 0.85
    )
  }
  
  interval_labels <- colnames(boxes_1)
  NumberOfIntervals <- length(interval_labels)
  plot_box_data <- c(
    lapply(seq_len(NumberOfIntervals), function(i) boxes_1[, i]),
    lapply(seq_len(NumberOfIntervals), function(i) boxes_2[, i])
  )
  sample_1_positions <- seq_len(NumberOfIntervals)
  sample_2_positions <- seq_len(NumberOfIntervals) + NumberOfIntervals + 2
  box_positions <- c(sample_1_positions, sample_2_positions)
  group_midpoints <- c(mean(sample_1_positions), mean(sample_2_positions))
  box_stats <- lapply(plot_box_data, grDevices::boxplot.stats)
  box_ylim <- range(
    unlist(lapply(box_stats, function(z) z$stats)), finite = TRUE
  )
  box_pad <- diff(box_ylim)*0.08
  if (!is.finite(box_pad) || box_pad == 0) box_pad <- 0.5
  
  grDevices::pdf(
    output_pdf, width = pdf_width, height = pdf_height,
    onefile = TRUE, useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  
  graphics::layout(
    matrix(c(0, 1, 0, 2, 0,
             0, 3, 3, 3, 0), nrow = 2, ncol = 5, byrow = TRUE),
    widths = c(1.2, 2.4, 0.8, 2.4, 1.2),
    heights = c(1, 1.3)
  )
  graphics::par(oma = c(0, 0, 3.5, 0))
  
  graphics::par(mar = c(4.3, 4.7, 3.2, 4.7))
  PlotOverlay(result_1, sample_1_label, profile_colours[1])
  
  graphics::par(mar = c(4.3, 4.7, 3.2, 4.7))
  PlotOverlay(result_2, sample_2_label, profile_colours[2])
  
  graphics::par(mar = c(6.3, 5, 3.2, 1.2))
  graphics::boxplot(
    plot_box_data, at = box_positions,
    col = c(rep(box_colours[1], NumberOfIntervals),
            rep(box_colours[2], NumberOfIntervals)),
    border = c(rep(profile_colours[1], NumberOfIntervals),
               rep(profile_colours[2], NumberOfIntervals)),
    outline = FALSE, xaxt = "n",
    ylim = box_ylim + c(-box_pad, box_pad),
    ylab = ChIPYAxisLabel,
    main = "Per-origin ChIP signal by absolute distance interval - old approach",
    bty = "n", las = 1, cex.axis = 1.12,
    cex.lab = 1.18, cex.main = 1.12
  )
  
  set.seed(123)
  for (i in seq_along(plot_box_data)) {
    visible <- plot_box_data[[i]] >= box_stats[[i]]$stats[1] &
      plot_box_data[[i]] <= box_stats[[i]]$stats[5]
    visible_values <- plot_box_data[[i]][visible]
    SampleIndex <- if (i <= NumberOfIntervals) 1L else 2L
    graphics::points(
      jitter(rep(box_positions[i], length(visible_values)), amount = 0.16),
      visible_values, pch = 16, cex = point_cex,
      col = grDevices::adjustcolor(
        profile_colours[SampleIndex], alpha.f = point_alpha
      )
    )
  }
  
  graphics::axis(
    1, at = box_positions,
    labels = rep(interval_labels, 2), tick = FALSE,
    line = 0.15, cex.axis = 0.82, las = 1
  )
  graphics::axis(
    1, at = group_midpoints,
    labels = c(sample_1_label, sample_2_label),
    tick = FALSE, line = 2, cex.axis = 1.08
  )
  graphics::mtext(
    "Absolute distance from origin centre",
    side = 1, line = 4.1, cex = 1.08
  )
  graphics::mtext(
    plot_title, side = 3, line = 1.3,
    outer = TRUE, font = 2, cex = 1.25
  )
  
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  
  MakeLongTable <- function(Assay, Values_1, Values_2) {
    do.call(
      rbind,
      lapply(seq_along(interval_labels), function(i) {
        rbind(
          data.frame(
            origin_index = origins$origin_index,
            oriName = origins$oriName, chrom = origins$chrom,
            oriCenter = origins$oriCenter, sample = sample_1_label,
            assay = Assay, interval = interval_labels[i],
            value = Values_1[, i], stringsAsFactors = FALSE
          ),
          data.frame(
            origin_index = origins$origin_index,
            oriName = origins$oriName, chrom = origins$chrom,
            oriCenter = origins$oriCenter, sample = sample_2_label,
            assay = Assay, interval = interval_labels[i],
            value = Values_2[, i], stringsAsFactors = FALSE
          )
        )
      })
    )
  }
  
  long_table <- rbind(
    MakeLongTable("ChIP", boxes_1, boxes_2),
    MakeLongTable("BrDU", BrDU_boxes_1, BrDU_boxes_2)
  )
  
  message("Saved old-approach ChIP-BrDU comparison PDF: ", output_pdf)
  
  legend_note <- paste(
    "Profiles use the legacy strand-wise ratio-of-medians method: median IP",
    "is divided by median Input independently for Watson and Crick, each",
    "strand is smoothed and the two enrichments are added. Interval boxes",
    "apply the analogous aggregate-before-division calculation within each",
    "origin and distance interval. ChIP is solid on the left axis and BrDU",
    "is dotted on the right axis; value 1 is aligned within each panel."
  )
  
  result <- list(
    output_pdf = output_pdf,
    alignment = Alignment,
    method = "old strand-wise ratio-of-medians approach",
    yLim = yLim,
    signal_source = paste0(
      RatioFolder, " Watson/Crick tables: ip.score and in.score"
    ),
    coverage_step_bp = step,
    origins = origins,
    distance_bp = distance_bp,
    profile_limit_kb = profile_limit_kb,
    interval_breaks_kb = distance_breaks_kb,
    legend_note = legend_note,
    sample_1 = list(
      label = sample_1_label,
      ChIP = list(
        profile = result_1$ChIP$profile,
        interval_values = boxes_1,
        files = sample_1$paths[c("ChIP_watson", "ChIP_crick")]
      ),
      BrDU = list(
        profile = result_1$BrDU$profile,
        interval_values = BrDU_boxes_1,
        files = sample_1$paths[c("BrDU_watson", "BrDU_crick")]
      )
    ),
    sample_2 = list(
      label = sample_2_label,
      ChIP = list(
        profile = result_2$ChIP$profile,
        interval_values = boxes_2,
        files = sample_2$paths[c("ChIP_watson", "ChIP_crick")]
      ),
      BrDU = list(
        profile = result_2$BrDU$profile,
        interval_values = BrDU_boxes_2,
        files = sample_2$paths[c("BrDU_watson", "BrDU_crick")]
      )
    ),
    interval_data = long_table
  )
  
  if (isTRUE(return_data)) invisible(result) else invisible(output_pdf)
}


