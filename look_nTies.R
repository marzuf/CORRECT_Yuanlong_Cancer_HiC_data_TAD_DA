### LOOK HOW MANY TIES IN THE 
# - FC
# - MEANCORR
# - AVG RANK
### in the observed and the permutation data

# Rscript look_nTies.R 

# Rscript look_nTies.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS # to run 1 dataset at time


script_name <- "look_nTies.R"

cat(paste0("> START ", script_name, "\n"))

startTime <- Sys.time()

setDir <- ""

outFold <- "LOOK_NTIES"
dir.create(outFold)

plotType <- "png"
myHeight <- 7
myWidth <- 10

require(foreach)
require(doMC)
require(ggplot2)
require(hrbrthemes)

registerDoMC(40)

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script6_name <- "6_runPermutationsMeanLogFC"
script7_name <- "7_runPermutationsMeanTADCorr"

# rank 1 -> high abs FC

tiesMeth <- "min"

buildData <- FALSE

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

outNames <- paste0(outNames1, "_", outNames2)

if(buildData) {
  all_ds_ties <- foreach(fc_permut_file = all_permutFC_files) %dopar% {
    
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
    
    ####################################
    ###### NUMBER OF TIES
    ####################################
    
    ### -> observed meanFC
    all_table_meanFC <- setNames(as.numeric(table(meanFC_obs)), names(table(meanFC_obs)))
    obs_meanFC_ties <- 1 - sum(all_table_meanFC==1)/length(meanFC_obs)
    
    
    ### -> observed meanCorr
    all_table_meanCorr <- setNames(as.numeric(table(meanCorr_obs)), names(table(meanCorr_obs)))
    obs_meanCorr_ties <- 1 - sum(all_table_meanCorr==1)/length(meanCorr_obs)
    
    stopifnot(setequal(names(meanFC_obs), names(meanCorr_obs)))
    
    ### -> observed meanRank
    obs_avgRank <- (meanCorr_obs_rank + meanFC_obs_rank)/2
    all_table_avgRank <- setNames(as.numeric(table(obs_avgRank)), names(table(obs_avgRank)))
    obs_avgRank_ties <- 1 - sum(all_table_avgRank==1)/length(all_table_avgRank)
    
    ### -> permut meanFC
    permut_meanFC_ties <- apply(meanFC_permDT, 2, function(perm_meanFC) {
      stopifnot(length(perm_meanFC) == length(all_regions))
      perm_all_table_meanFC <- setNames(as.numeric(table(perm_meanFC)), names(table(perm_meanFC)))
      1 - sum(perm_all_table_meanFC==1)/length(perm_meanFC)
    })
    
    ### -> permut meanFC
    permut_meanCorr_ties <- apply(meanCorr_permDT, 2, function(perm_meanCorr) {
      stopifnot(length(perm_meanCorr) == length(all_regions))
      perm_all_table_meanCorr <- setNames(as.numeric(table(perm_meanCorr)), names(table(perm_meanCorr)))
      1 - sum(perm_all_table_meanCorr==1)/length(perm_meanCorr)
    })
    
    ### -> permut avg rank
    permut_avgRank_ties <- apply(rank_avgCorrFC_permDT, 2, function(perm_avgRank) {
      stopifnot(length(perm_avgRank) == length(all_regions))
      perm_all_table_meanRank <- setNames(as.numeric(table(perm_avgRank)), names(table(perm_avgRank)))
      1 - sum(perm_all_table_meanRank==1)/length(perm_avgRank)
    })
    list(
      obs_meanFC_ties = obs_meanFC_ties,
      obs_meanCorr_ties = obs_meanCorr_ties,
      obs_avgRank_ties = obs_avgRank_ties,
      permut_meanFC_ties = permut_meanFC_ties,
      permut_meanCorr_ties = permut_meanCorr_ties,
      permut_avgRank_ties = permut_avgRank_ties
    )
  }
  names(all_ds_ties) <- outNames
  
  outFile <- file.path(outFold, "all_ds_ties.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_ds_ties, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFold, "all_ds_ties.Rdata")
  stopifnot(file.exists(outFile))
  all_ds_ties <- eval(parse(text = load(outFile)))
}


### PLOT
all_obs_meanFC_ties <- unlist(lapply(all_ds_ties, function(x) x[["obs_meanFC_ties"]]))
all_obs_meanCorr_ties <- unlist(lapply(all_ds_ties, function(x) x[["obs_meanCorr_ties"]]))
all_obs_avgRank_ties <- unlist(lapply(all_ds_ties, function(x) x[["obs_avgRank_ties"]]))
all_permut_meanFC_ties <- unlist(lapply(all_ds_ties, function(x) x[["permut_meanFC_ties"]]))
all_permut_meanCorr_ties <- unlist(lapply(all_ds_ties, function(x) x[["permut_meanCorr_ties"]]))
all_permut_avgRank_ties <- unlist(lapply(all_ds_ties, function(x) x[["permut_avgRank_ties"]]))
  

plot_DT <- data.frame(
  dataset = c(
    names(all_obs_meanFC_ties),
    names(all_obs_meanCorr_ties),
    names(all_obs_avgRank_ties),
    names(all_permut_meanFC_ties),
    names(all_permut_meanCorr_ties),
    names(all_permut_avgRank_ties)
    ), 
  plot_var =  c(
    rep("obs_meanFC", length(all_obs_meanFC_ties)),
    rep("obs_meanCorr", length(all_obs_meanCorr_ties)),
    rep("obs_avgRank", length(all_obs_avgRank_ties)),
    rep("permut_meanFC", length(all_permut_meanFC_ties)),
    rep("permut_meanCorr", length(all_permut_meanCorr_ties)),
    rep("permut_avgRank", length(all_permut_avgRank_ties))
    ), 
  value = c(
    as.numeric(all_obs_meanFC_ties),
    as.numeric(all_obs_meanCorr_ties),
    as.numeric(all_obs_avgRank_ties),
    as.numeric(all_permut_meanFC_ties),
    as.numeric(all_permut_meanCorr_ties),
    as.numeric(all_permut_avgRank_ties)
    ),
  stringsAsFactors = FALSE
)


outFile <- file.path(outFold, paste0("all_ratioTies_boxplot.", plotType))
p <- ggplot(plot_DT, aes(x = plot_var, y= value))+ #, fill = featureDataset, col=featureDataset))  +
  geom_boxplot()+
  scale_fill_ipsum() +
  scale_color_ipsum() +
  theme_ipsum_rc(grid="XY", axis_title_just="c", axis_title_size=14) + 
  labs(x="", 
       y=paste0("Ratio of ties"),
       title="Ties in ranking",
       # caption="Brought to you by the letter 'g'",
       # subtitle=paste0("nDS = ", length(unique(plot_DT$dataset))),
       colour = "",
       fill = ""
  ) 
ggsave(p, filename = outFile, height = myHeight, width = myWidth)
cat(paste0("... written: ", outFile, "\n"))

# THEMES FROM hrbrthemes package
# theme_ft_rc(grid="XY", axis_title_just="c") # dark
# theme_ipsum(grid="XY", axis_title_just="c")


################################****************************************************************************************
################################****************************************************************************************

cat("*** DONE ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





