
# Rscript extract_Yuanlong_FINAL_DOMAINS_corrected.R

# => CORRECTED VERSION:
# 1) TADs on left and right of a gap might have same TAD rank
# so now impose a threshold:
# if left and right are separated by more than "gap_threshold" bins => they receive a different CD_rank
# 2) discard TADs if they contain the mid position of the centromere

hicds = "ENCSR346DCU_LNCaP_40kb"
chromo = "chr2"

script_name <- "extract_Yuanlong_FINAL_DOMAINS_corrected.R"
startTime <- Sys.time()

options(scipen=100)

SSHFS <- FALSE
prefixDir <- ifelse(SSHFS, "/media/electron", "")
allData_file <- file.path(prefixDir, "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/bin_comp_CELL_LINEs_list.Rdata")
stopifnot(file.exists(allData_file))

nCpu <- ifelse(SSHFS, 2, 40)
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu)


suppressPackageStartupMessages(library(rCGH, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
centroStart <- setNames(hg19$centromerStart, paste0("chr", hg19$chrom))
stopifnot(is.numeric(centroStart))
centroEnd <- setNames(hg19$centromerEnd, paste0("chr", hg19$chrom))
stopifnot(is.numeric(centroEnd))


gap_threshold <- 50  # IMPORTANT SETTING !!!
binSize <- 40000
binSizeKb <- 40

cat("!!! IMPORTANT HARD-CODED SETTINGS:\n")
cat(paste0("... gap_threshold\t=\t", gap_threshold, "\n"))
cat(paste0("... binSize\t=\t", binSize, "\n"))
cat(paste0("... binSizeKb\t=\t", binSizeKb, "\n"))

allData <- eval(parse(text = load(allData_file)))

cancerDS <- names(allData)[grep("^mega_", names(allData))]

ds = cancerDS[1]
# cancerDS=cancerDS[1]
foo <- foreach(ds = cancerDS) %dopar% {
  
  cat(paste0("... start DS: ", ds, "\n"))
  
  folderName <- gsub("^mega_", "", ds) 
  outFolder <- file.path(paste0(folderName, "_",binSizeKb , "kb"), "FINAL_DOMAINS_CORRECTED")
  dir.create(outFolder, recursive = TRUE)
  
  allChr_tad_DT <- allData[[ds]]
  stopifnot( colnames(allChr_tad_DT) == c("chr", "bin_index", "CD_rank"))
  allChr_tad_DT$chr <- paste0("chr", as.character(allChr_tad_DT$chr))
  
  all_chrs <- unique(allChr_tad_DT$chr)
  
  chromo="chr1"
  for(chromo in all_chrs) {
    
    cat(paste0("...... start chromo: ", chromo, "\n"))
    
    chr_DT <- allChr_tad_DT[allChr_tad_DT$chr == chromo, ,drop=FALSE]
    ctr_start <- as.numeric(centroStart[paste0(chromo)])
    ctr_end <- as.numeric(centroEnd[paste0(chromo)])
    ctr_midPos <- (ctr_start+ctr_end)/2
    
    ##################################################### START CORRECTED PART (1)
    # CORRECTED -> CHANGED HERE BECAUSE THEY MIGHT HAVE BIG GAPS !!!
    chr_DT$CD_rank_old <- chr_DT$CD_rank
    
    i_rank <- 1
    CD_rank_new <- c(i_rank, rep(NA, nrow(chr_DT)-1))
    stopifnot(length(CD_rank_new) == nrow(chr_DT))
    
    
    
    for(k in 2:nrow(chr_DT)) {
      curr_bin <- as.numeric(chr_DT$bin_index[k])
      prev_bin <- as.numeric(chr_DT$bin_index[k-1])
      stopifnot(!is.na(curr_bin))
      stopifnot(!is.na(prev_bin))
      bin_gap <- curr_bin - prev_bin
      stopifnot(bin_gap > 0)
      curr_rank <- chr_DT$CD_rank_old[k]
      prev_rank <- chr_DT$CD_rank_old[k-1]
      if(curr_rank != prev_rank | bin_gap > gap_threshold ){
        i_rank <- i_rank + 1
      } 
      CD_rank_new[k] <- i_rank
    }
    stopifnot(range(diff(CD_rank_new)) %in% c(0,1))
    stopifnot(!is.na(CD_rank_new))
    stopifnot(length(CD_rank_new) == nrow(chr_DT))
    chr_DT$CD_rank <- CD_rank_new
    ##################################################### END CORRECTED PART (1)
    
    
    chr_DT$CD_rank <- factor(chr_DT$CD_rank, levels = unique(chr_DT$CD_rank))
    chr_DT$domainNbr <- as.numeric(chr_DT$CD_rank)
    stopifnot( diff(chr_DT$domainNbr) == 0 | diff(chr_DT$domainNbr) == 1 )
    
    all_bins_tadPos_DT <- chr_DT
    all_bins_tadPos_DT$CD_rank <- NULL
    head( all_bins_tadPos_DT)
    
    stopifnot( diff(as.numeric(all_bins_tadPos_DT$bin_index)) > 0 )
    
    tadPos_DT <- do.call(rbind, by(all_bins_tadPos_DT, all_bins_tadPos_DT$domainNbr, function(x) {
      regionIdx <- unique(x$domainNbr)
      stopifnot(length(regionIdx) == 1)
      data.frame(
        chr = chromo,
        region = paste0(chromo, "_TAD", regionIdx),
        startIdx = min(as.numeric(as.character(x$bin_index))),
        endIdx = max(as.numeric(as.character(x$bin_index))),
        stringsAsFactors = FALSE
      )
    }))
    stopifnot(nrow(tadPos_DT) == length(unique(all_bins_tadPos_DT$domainNbr)))
    tadPos_DT$startIdx <- as.numeric(as.character(tadPos_DT$startIdx))
    stopifnot(!is.na(tadPos_DT$startIdx))
    tadPos_DT$endIdx <- as.numeric(as.character(tadPos_DT$endIdx))
    stopifnot(!is.na(tadPos_DT$endIdx))
    stopifnot(tadPos_DT$endIdx >= tadPos_DT$startIdx)
    # pos_start = (bin_index - 1)*40E3 + 1; pos_end = bin_index*40E3
    # e.g. bin1 -> 1 - 40'000
    tadPos_DT$start <- (tadPos_DT$startIdx - 1)*binSize + 1 
    tadPos_DT$end <- (tadPos_DT$endIdx)*binSize 
    outDT <- tadPos_DT[,c("chr", "start", "end"), drop=FALSE]
    stopifnot(is.numeric(outDT$start))
    stopifnot(is.numeric(outDT$end))
    stopifnot(outDT$chr == chromo)
    
    ##################################################### START CORRECTED PART (2)
    # remove the TAD if it contains the midPos
    nrow0 <- nrow(tadPos_DT)
    ctr_bins <- which(tadPos_DT$start <= ctr_midPos & tadPos_DT$end >= ctr_midPos)
    if(length(ctr_bins) > 0) {
      stopifnot(length(ctr_bins) == 1)
      tadPos_DT <- tadPos_DT[-ctr_bins,]
      stopifnot(nrow(tadPos_DT) == nrow0-1)
    } else {
      stopifnot(nrow(tadPos_DT) == nrow0)
    }
    ##################################################### END CORRECTED PART (2)
    outDT <- outDT[order(outDT$start, outDT$end),,drop=FALSE]
    stopifnot(outDT$end > outDT$start)
    
    stopifnot(outDT$end[-nrow(outDT)] < outDT$start[-1])
    
    # output file name should like GSE105194_cerebellum_40kb/FINAL_DOMAINS/GSE105194_cerebellum_chr3_KR_40kb_final_domains.txt
    outFile <- file.path(outFolder, paste0(folderName, "_", chromo, "_", "YL", "_", binSizeKb, "kb_final_domains.txt"))
    write.table(outDT, file = outFile, col.names = FALSE, row.names = FALSE, append = FALSE, sep="\t", quote = FALSE)
    cat(paste0("... written: ", outFile, "\n"))

  } # end for-iterating over chromo    
} # end for-iterating over ds

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))


