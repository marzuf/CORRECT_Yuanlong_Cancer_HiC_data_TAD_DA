# Rscript empPval_avg_corrFC_rank.R 

# Rscript empPval_avg_corrFC_rank.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS # to run 1 dataset at time

withDiago <- FALSE

script_name <- "empPval_meanTADFC_rank.R"

cat(paste0("> START ", script_name, "\n"))

startTime <- Sys.time()

setDir <- ""

outFold <- "EMPPVAL_AVG_CORRFC_RANK"
dir.create(outFold)

require(foreach)
require(doMC)
registerDoMC(40)

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script6_name <- "6_runPermutationsMeanLogFC"
script7_name <- "7_runPermutationsMeanTADCorr"

# rank 1 -> high abs FC

tiesMeth <- "min"

###################
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_permutFC_files <- list.files(pipOutFolder, recursive = TRUE, pattern="meanLogFC_permDT.Rdata", full.names = TRUE)
stopifnot(length(all_permutFC_files) > 0)

fc_permut_file = all_permutFC_files[1]

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 2) {
  all_permutFC_files <- all_permutFC_files[grepl(args[1], all_permutFC_files) & grepl(args[2], all_permutFC_files)]
  stopifnot(length(all_permutFC_files) == 1)
}

for(fc_permut_file in all_permutFC_files) {

  hicds <- basename(dirname(dirname(dirname(fc_permut_file))))
  exprds <- basename(dirname(dirname(fc_permut_file)))
  
  cat("... ", hicds, " - ", exprds, "\n")
  
  corr_permut_file <- file.path(pipOutFolder, hicds, exprds, script7_name, "meanCorr_permDT.Rdata")
  stopifnot(file.exists(corr_permut_file))
  
  #################
  #### FC
  #################
  # observed data - FC
  obs_meanFC_TAD_file <- file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
  stopifnot(file.exists(obs_meanFC_TAD_file))
  meanFC_obs <- eval(parse(text = load(obs_meanFC_TAD_file)))
  # FC -> consider up- and down- equally -> take the abs value
  meanFC_obs <- abs(meanFC_obs)
  
  # take -x to have decreasing ranking
  meanFC_obs_rank <- rank(-meanFC_obs, ties = tiesMeth)
  all_regions <- names(meanFC_obs_rank)
  
  cat("...... load FC permDT", "\n")  
  meanFC_permDT <- eval(parse(text = load(fc_permut_file)))
  stopifnot( setequal(all_regions, rownames(meanFC_permDT) ))

  # FC -> consider up- and down- equally -> take the abs value
  meanFC_permDT <- abs(meanFC_permDT)
    
  # take -x to have decreasing ranking
  rank_meanFC_permDT <- apply(meanFC_permDT, 2 , function(x) rank(-x, ties=tiesMeth))
  max1_idx <- which(rownames(meanFC_permDT) == names(which(rank_meanFC_permDT[,1] == 1)))
  stopifnot(meanFC_permDT[-max1_idx,1] <  meanFC_permDT[max1_idx,1])
  
  stopifnot( setequal(all_regions, rownames(rank_meanFC_permDT) ))
  
  
  #################
  #### mean corr
  #################
  # observed data
  obs_meanCorr_TAD_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
  stopifnot(file.exists(obs_meanCorr_TAD_file))
  meanCorr_obs <- eval(parse(text = load(obs_meanCorr_TAD_file)))
  
  # take -x to have decreasing ranking
  meanCorr_obs_rank <- rank(-meanCorr_obs, ties = tiesMeth)
  all_regionsCorr <- names(meanCorr_obs_rank)
  
  stopifnot(setequal(all_regions, all_regionsCorr))
  
  cat("...... load meanCorr permDT", "\n")  
  meanCorr_permDT <- eval(parse(text = load(corr_permut_file)))
  stopifnot( setequal(all_regionsCorr, rownames(meanCorr_permDT) ))
  
  # take -x to have decreasing ranking
  rank_meanCorr_permDT <- apply(meanCorr_permDT, 2 , function(x) rank(-x, ties=tiesMeth))
  max1_idx <- which(rownames(meanCorr_permDT) == names(which(rank_meanCorr_permDT[,1] == 1)))
  stopifnot(meanCorr_permDT[-max1_idx,1] <  meanCorr_permDT[max1_idx,1])
  
  stopifnot( setequal(all_regionsCorr, rownames(rank_meanCorr_permDT) ))
  
  rank_avgCorrFC_permDT <- (rank_meanCorr_permDT[all_regions,] + rank_meanFC_permDT[all_regions,])/2
  
  stopifnot(setequal(all_regions ,rownames(rank_avgCorrFC_permDT)))
  
  # compute the empirical p-value  
  
  emp_pval_avg_corrFC_rank <- unlist(foreach(reg = all_regions, .combine='c') %dopar% {
    ######
    # FC
    ######
    # select the random values for this region
    all_shuff_meanFC_rank <- rank_meanFC_permDT[paste0(reg),]
    stopifnot(length(all_shuff_meanFC_rank) == ncol(rank_meanFC_permDT))
    # select the observed value for this region
    obs_FC_rank <- meanFC_obs_rank[reg]
    stopifnot(length(obs_FC_rank) == 1 )
    
    ######
    # meanCorr
    ######
    # select the random values for this region
    all_shuff_meanCorr_rank <- rank_meanCorr_permDT[paste0(reg),]
    stopifnot(length(all_shuff_meanCorr_rank) == ncol(rank_meanCorr_permDT))
    # select the observed value for this region
    obs_meanCorr_rank <- meanCorr_obs_rank[reg]
    stopifnot(length(obs_meanCorr_rank) == 1 )
    # the number of times the permut rank is smaller than the observed one
    # => the pvalue will be high if the permutation often a smaller rank
    
    
    obs_avgRank <- mean(obs_FC_rank, obs_meanCorr_rank)
    
    all_shuff_avg_rank <- rank_avgCorrFC_permDT[paste0(reg),]
    
    
    
    # the number of times the permut rank is smaller than the observed one
    # => the pvalue will be high if the permutation often a smaller rank
    emp_pval <- sum(all_shuff_avg_rank <= obs_avgRank) 
    # emp_pval/length(shuff_intraTADcorr_region)
    ### ADD THE +1 => THIS IS FOR THE OBSERVED VALUE -> AVOID THE 0 VALUES !!!
    (emp_pval+1)/(length(all_shuff_avg_rank)+1)
  })
  names(emp_pval_avg_corrFC_rank) <- all_regions
  
  stopifnot(all(emp_pval_avg_corrFC_rank > 0 & emp_pval_avg_corrFC_rank <= 1 ))
  
  
  outFile <- file.path(outFold, hicds, exprds, "emp_pval_avg_corrFC_rank.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(emp_pval_avg_corrFC_rank, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
}

################################****************************************************************************************
################################****************************************************************************************

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))








