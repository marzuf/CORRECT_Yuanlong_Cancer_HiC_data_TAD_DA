

# Rscript distRanking_with_permut_foreach.R 

# Rscript distRanking_with_permut_foreach.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS # to run 1 dataset at time


script_name <- "distRanking_with_permut.R"

cat(paste0("> START ", script_name, "\n"))

startTime <- Sys.time()

setDir <- ""

buildData <- FALSE

outFold <- "DISTRANKING_WITH_PERMUT_ALL"
dir.create(outFold)

require(foreach)
require(doMC)
registerDoMC(80)

source("canberra_stability_fct.R")

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script6_name <- "6_runPermutationsMeanLogFC"
script7_name <- "7_runPermutationsMeanTADCorr"


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

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

outNames1 <- basename(dirname(dirname(dirname(all_permutFC_files))))
outNames2 <- basename(dirname(dirname(all_permutFC_files)))
outNames <- file.path(outNames1, outNames2)


if(buildData){
  all_ds_canberraDist <- foreach(fc_permut_file = all_permutFC_files) %do% {
    
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
    
    stopifnot(setequal(all_regions ,rownames(rank_avgCorrFC_permDT)))
    stopifnot(setequal(all_regions ,rownames(rank_meanFC_permDT)))
    stopifnot(setequal(all_regions ,rownames(rank_meanCorr_permDT)))
    
    stopifnot(ncol(rank_meanCorr_permDT) == ncol(rank_meanFC_permDT) )
    stopifnot(ncol(rank_meanCorr_permDT) == ncol(rank_avgCorrFC_permDT) )
    
    
    nPerm <- ncol(rank_meanFC_permDT)
    
    obs_avgRank <- (meanCorr_obs_rank + meanFC_obs_rank)/2
    
    curr_ds_all_dist <- foreach(i_col = 1:nPerm) %dopar% {
      
      cat("...  ", hicds," - ", exprds, " : i_col = ", i_col, "/", nPerm, "\n")
      
      
      perm_meanCorr <- rank_meanCorr_permDT[,i_col]
      
      stopifnot(length(perm_meanCorr) == length(all_regions))
      stopifnot(length(perm_meanCorr) == length(meanCorr_obs))
      
      all_meanCorr_ranks_mat <- matrix(c(perm_meanCorr, meanCorr_obs_rank),
                                       byrow = TRUE,
                                       nrow = 2)
      stopifnot(ncol(all_meanCorr_ranks_mat) == length(all_regions))
      stopifnot(nrow(all_meanCorr_ranks_mat) == 2)
      
      stopifnot(all_meanCorr_ranks_mat[1,] %in% c(1:length(all_regions)))
      stopifnot(all_meanCorr_ranks_mat[2,] %in% c(1:length(all_regions)))
      
      # RANKING SHOULD BE 0-BASED
      all_meanCorr_ranks_mat <- all_meanCorr_ranks_mat - 1
      # EACH ROW = ranked list
      # COLUMNS = features
      obs_perm_dist_meanCorr <- canberra_stability(all_meanCorr_ranks_mat)
      
      perm_meanFC <- rank_meanFC_permDT[,i_col]
      
      stopifnot(length(perm_meanFC) == length(all_regions))
      stopifnot(length(perm_meanFC) == length(meanCorr_obs))
      
      all_meanFC_ranks_mat <- matrix(c(perm_meanFC, meanFC_obs_rank),
                                     byrow = TRUE,
                                     nrow = 2)
      stopifnot(ncol(all_meanFC_ranks_mat) == length(all_regions))
      stopifnot(nrow(all_meanFC_ranks_mat) == 2)
      
      stopifnot(all_meanFC_ranks_mat[1,] %in% c(1:length(all_regions)))
      stopifnot(all_meanFC_ranks_mat[2,] %in% c(1:length(all_regions)))
      
      # RANKING SHOULD BE 0-BASED
      all_meanFC_ranks_mat <- all_meanFC_ranks_mat - 1
      # EACH ROW = ranked list
      # COLUMNS = features
      obs_perm_dist_meanFC <- canberra_stability(all_meanFC_ranks_mat)
      
      
      perm_avgRank <- rank_avgCorrFC_permDT[,i_col]
      
      stopifnot(length(perm_avgRank) == length(all_regions))
      stopifnot(length(perm_avgRank) == length(obs_avgRank))
      
      all_avgRank_ranks_mat <- matrix(c(perm_avgRank, obs_avgRank),
                                      byrow = TRUE,
                                      nrow = 2)
      stopifnot(ncol(all_avgRank_ranks_mat) == length(all_regions))
      stopifnot(nrow(all_avgRank_ranks_mat) == 2)
      
      # RANKING SHOULD BE 0-BASED
      all_avgRank_ranks_mat <- all_avgRank_ranks_mat - 1
      # EACH ROW = ranked list
      # COLUMNS = features
      obs_perm_dist_avgRank<- canberra_stability(all_avgRank_ranks_mat, checkZero = FALSE)
      
      list(
        obs_perm_dist_meanCorr = obs_perm_dist_meanCorr,
        obs_perm_dist_meanFC = obs_perm_dist_meanFC,
        obs_perm_dist_avgRank = obs_perm_dist_avgRank
      )
    }# end-iterating over all permut for this ds
    curr_ds_all_dist
  }
  names(all_ds_canberraDist) <- outNames
  
  outFile <- file.path(outFold, "all_ds_canberraDist.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_ds_canberraDist, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else{
  outFile <- file.path(outFold, "all_ds_canberraDist.Rdata")
  all_ds_canberraDist <- eval(parse(text = load(outFile)))
}


all_obs_perm_dist_meanCorr <- unlist(lapply(all_ds_canberraDist, function(ds) {
  lapply(ds, function(x) x[["obs_perm_dist_meanCorr"]]) }))
all_obs_perm_dist_meanFC <- unlist(lapply(all_ds_canberraDist, function(ds) {
  lapply(ds, function(x) x[["obs_perm_dist_meanFC"]]) }))
all_obs_perm_dist_avgRank <- unlist(lapply(all_ds_canberraDist, function(ds) {
  lapply(ds, function(x) x[["obs_perm_dist_avgRank"]]) }))

outFile <- file.path(outFold, paste0("corr_fc_avgRank_canberraDist_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(
    meanCorr = all_obs_perm_dist_meanCorr,
    meanFC = all_obs_perm_dist_meanFC,
    avgRank = all_obs_perm_dist_avgRank
  ),
  legPos = "topleft",
  plotTit = paste0("Canberra dist. obs. <-> permut. rankings")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################****************************************************************************************
################################****************************************************************************************

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





