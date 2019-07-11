# Rscript clean_empPval_rank_cmp.R

script_name <- "clean_empPval_rank_cmp.R"

startTime <- Sys.time()

cat("> START clean_empPval_rank_cmp.R \n")

SSHFS <- FALSE

setDir <- ifelse(SSHFS, "/media/electron", "")

pipFold <- file.path("PIPELINE", "OUTPUT_FOLDER")

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("subtype_cols.R")
source("look_tads_fct.R")
source("plot_lolliTAD_funct.R")

nTop <- 5

outFold <- "CLEAN_EMPPVAL_RANK_CMP"
dir.create(outFold, recursive = TRUE)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.4

outWidthGG <- 20
outHeightGG <- min(c(7 * nTop/2, 49))


myCexAxis <- 1.2
myCexLab <- 1.2

windowSizeBp <- 500*10^3
options(scipen=100)

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)

all_pvalComb_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = TRUE)
stopifnot(length(all_pvalComb_files) > 0)

all_pvalFC_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanLogFC.Rdata", full.names = TRUE)
stopifnot(length(all_pvalFC_files) > 0)

all_pvalCorr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanCorr.Rdata", full.names = TRUE)
stopifnot(length(all_pvalCorr_files) > 0)

all_ratio_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_ratioDown.Rdata", full.names = TRUE)
stopifnot(length(all_pvalCorr_files) > 0)

all_pvalComb_rank_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined_rank.Rdata", full.names = TRUE)
stopifnot(length(all_pvalComb_rank_files) > 0)

all_pvalFC_rank_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanLogFC_rank.Rdata", full.names = TRUE)
stopifnot(length(all_pvalFC_rank_files) > 0)

all_pvalCorr_rank_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanCorr_rank.Rdata", full.names = TRUE)
stopifnot(length(all_pvalCorr_rank_files) > 0)


rD_empPval_DT <- foreach(corr_file = all_ratio_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    ratioDown = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}


mC_empPval_DT <- foreach(corr_file = all_pvalCorr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    meanCorr_empPval = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}

FC_empPval_DT <- foreach(corr_file = all_pvalFC_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    meanFC_empPval = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}

comb_empPval_DT <- foreach(corr_file = all_pvalComb_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    comb_empPval = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}


mCrank_empPval_DT <- foreach(corr_file = all_pvalCorr_rank_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    meanCorrRank_empPval = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}

FCrank_empPval_DT <- foreach(corr_file = all_pvalFC_rank_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    meanFCrank_empPval = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}

combRank_empPval_DT <- foreach(corr_file = all_pvalCorr_rank_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    combRank_empPval = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}

stopifnot(nrow(rD_empPval_DT) == nrow(mC_empPval_DT) )
stopifnot(nrow(rD_empPval_DT) == nrow(FC_empPval_DT) )
stopifnot(nrow(rD_empPval_DT) == nrow(comb_empPval_DT) )
stopifnot(nrow(rD_empPval_DT) == nrow(mCrank_empPval_DT) )
stopifnot(nrow(rD_empPval_DT) == nrow(FCrank_empPval_DT) )
stopifnot(nrow(rD_empPval_DT) == nrow(combRank_empPval_DT) )

all_DT <- merge(combRank_empPval_DT, 
                merge(FCrank_empPval_DT, 
                      merge(mCrank_empPval_DT, 
                            merge(comb_empPval_DT, 
                                  merge(FC_empPval_DT, 
                                        merge(rD_empPval_DT, mC_empPval_DT, by =c("hicds", "exprds", "region")),
                                          by=c("hicds", "exprds", "region")), 
                                          by=c("hicds", "exprds", "region")), 
                                          by=c("hicds", "exprds", "region")), 
                                          by=c("hicds", "exprds", "region")), 
                                          by=c("hicds", "exprds", "region"))


stopifnot(nrow(all_DT) == nrow(rD_empPval_DT) )
stopifnot(nrow(all_DT) == nrow(mC_empPval_DT) )
stopifnot(nrow(all_DT) == nrow(FC_empPval_DT) )
stopifnot(nrow(all_DT) == nrow(comb_empPval_DT) )
stopifnot(nrow(all_DT) == nrow(mCrank_empPval_DT) )
stopifnot(nrow(all_DT) == nrow(FCrank_empPval_DT) )
stopifnot(nrow(all_DT) == nrow(combRank_empPval_DT) )

all_DT$dataset <- paste0(all_DT$hicds, "_", all_DT$exprds)


outFile <- file.path(outFold, "all_DT.Rdata")
save(all_DT, file=outFile)
cat("... written: ", outFile, "\n")

all_DT$meanCorrRank_empPval_log10 <- log10(all_DT$meanCorrRank_empPval)
all_DT$meanFCrank_empPval_log10 <- log10(all_DT$meanFCrank_empPval)
all_DT$combRank_empPval_log10 <- log10(all_DT$combRank_empPval)

all_DT$meanCorr_empPval_log10 <- log10(all_DT$meanCorr_empPval)
all_DT$meanFC_empPval_log10 <- log10(all_DT$meanFC_empPval)
all_DT$comb_empPval_log10 <- log10(all_DT$comb_empPval)

############################################################################
############################################################################

xvar="meanCorr_empPval_log10"
yvar="meanCorrRank_empPval_log10"

myx <- all_DT[,xvar]
myy <- all_DT[,yvar]
outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  y = myy,
  x = myx,
  cex.lab=myCexLab,
  cex.axis=myCexAxis,
  ylab=paste0(yvar),
  xlab = paste0(xvar),
  cex = 0.5,
  main = paste0(yvar, " vs. ", xvar)
)
addCorr(
  x = myx,
  y=myy,
  legPos="topleft",
  bty="n"
)
foo <- dev.off()
cat("... written: ", outFile, "\n")


############################################################################
############################################################################

xvar="meanFC_empPval_log10"
yvar="meanFCrank_empPval_log10"

myx <- all_DT[,xvar]
myy <- all_DT[,yvar]
outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  y = myy,
  x = myx,
  cex.lab=myCexLab,
  cex.axis=myCexAxis,
  ylab=paste0(yvar),
  xlab = paste0(xvar),
  cex = 0.5,
  main = paste0(yvar, " vs. ", xvar)
)
addCorr(
  x = myx,
  y=myy,
  legPos="topleft",
  bty="n"
)
foo <- dev.off()
cat("... written: ", outFile, "\n")



############################################################################
############################################################################

xvar="meanCorrRank_empPval_log10"
yvar="meanFCrank_empPval_log10"

myx <- all_DT[,xvar]
myy <- all_DT[,yvar]
outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  y = myy,
  x = myx,
  cex.lab=myCexLab,
  cex.axis=myCexAxis,
  ylab=paste0(yvar),
  xlab = paste0(xvar),
  cex = 0.5,
  main = paste0(yvar, " vs. ", xvar)
)
addCorr(
  x = myx,
  y=myy,
  legPos="topleft",
  bty="n"
)
foo <- dev.off()
cat("... written: ", outFile, "\n")


############################################################################
############################################################################

xvar="meanCorr_empPval_log10"
yvar="meanFC_empPval_log10"

myx <- all_DT[,xvar]
myy <- all_DT[,yvar]
outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  y = myy,
  x = myx,
  cex.lab=myCexLab,
  cex.axis=myCexAxis,
  ylab=paste0(yvar),
  xlab = paste0(xvar),
  cex = 0.5,
  main = paste0(yvar, " vs. ", xvar)
)
addCorr(
  x = myx,
  y=myy,
  legPos="topleft",
  bty="n"
)
foo <- dev.off()
cat("... written: ", outFile, "\n")


############################################################################
############################################################################

xvar="comb_empPval_log10"
yvar="combRank_empPval_log10"

myx <- all_DT[,xvar]
myy <- all_DT[,yvar]
outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  y = myy,
  x = myx,
  cex.lab=myCexLab,
  cex.axis=myCexAxis,
  ylab=paste0(yvar),
  xlab = paste0(xvar),
  cex = 0.5,
  main = paste0(yvar, " vs. ", xvar)
)
addCorr(
  x = myx,
  y=myy,
  legPos="topleft",
  bty="n"
)
foo <- dev.off()
cat("... written: ", outFile, "\n")



############################################################################
############################################################################

all_DT$comb_empPval_adj <- p.adjust(all_DT$comb_empPval, method="BH")
all_DT$combRank_empPval_adj <- p.adjust(all_DT$combRank_empPval, method="BH")

all_DT$comb_empPval_adj_log10 <- -log10(all_DT$comb_empPval_adj)
all_DT$combRank_empPval_adj_log10 <- -log10(all_DT$combRank_empPval_adj) 


xvar="comb_empPval_adj_log10"
yvar="combRank_empPval_adj_log10"

myx <- all_DT[,xvar]
myy <- all_DT[,yvar]
outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  y = myy,
  x = myx,
  cex.lab=myCexLab,
  cex.axis=myCexAxis,
  ylab=paste0(yvar),
  xlab = paste0(xvar),
  cex = 0.5,
  main = paste0(yvar, " vs. ", xvar)
)
addCorr(
  x = myx,
  y=myy,
  legPos="topleft",
  bty="n"
)
foo <- dev.off()
cat("... written: ", outFile, "\n")

############################################################################
############################################################################

signifThresh <- 0.05

nSignif_DT <- do.call(rbind, by(all_DT, all_DT$dataset, function(subDT) {
  nSignifComb <- sum(subDT$comb_empPval <= signifThresh)
  nSignifCombRank <- sum(subDT$combRank_empPval <= signifThresh)
  data.frame(
    #dataset=unique(subDT$dataset),
    nSignifComb = nSignifComb,
    nSignifCombRank = nSignifCombRank,
    stringsAsFactors= FALSE
  )
}))

outFile <- file.path(outFold, paste0("nSignifTADs_adjPvalComb", ".", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(nSignif_DT, ylab="# signif. TADs adj. empPval comb.", main="# signif. TADs")
foo <- dev.off()
cat("... written: ", outFile, "\n")


############################################################################
############################################################################
# load("CLEAN_EMPPVAL_RANK_CMP/all_DT.Rdata")


all_ds <- unique(all_DT$dataset)

rank_var = "combRank_empPval_adj"

all_ranks_var <- c("comb_empPval_adj", "combRank_empPval_adj")



# all_ds=all_ds[1]
ds=all_ds[1]

for(rank_var in all_ranks_var) {
  
  all_DT <- all_DT[order(all_DT$dataset, all_DT[, rank_var]),]
  
  
  for(ds in all_ds) {
    
    sub_DT <- all_DT[all_DT$dataset == ds,]
    i_hicds <- unique(sub_DT$hicds)
    stopifnot(length(i_hicds) == 1)
    i_exprds <- unique(sub_DT$exprds)
    stopifnot(length(i_exprds) == 1)
    
    i=1
    # topTable_DT <- foreach(i = 1:nTop, .combine='rbind') %dopar% {
    plotList <- list()
    topTable_DT <- foreach(i = 1:nTop, .combine='rbind') %do% {
      
      i_tad <- basename(sub_DT$region[i])
      
      
      plotList[[i]] <- plot_lolliTAD_ds(exprds = i_exprds,
                                        hicds = i_hicds,
                                        all_TADs = i_tad,
                                        orderByLolli = "startPos")
      
      tmpDT <- get_tad_symbols_DT(hicds = i_hicds,
                                  exprds = i_exprds,
                                  region = i_tad,
                                  mainPipFolder = pipFold
      )
      otherCols <- colnames(tmpDT)
      tmpDT$iTop <- i
      tmpDT[,c("iTop", otherCols)]
    } # end-foreach iterating over the "n" top TADs
    # subset_vars (subsetVar) -> subset_levels (var) -> rankingVars (rank_var)
    outFile <- file.path(outFold, paste0(ds, "_", rank_var, "_nTop", nTop, ".txt" ))
    cat(outFile, "\n")
    write.table(topTable_DT, col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE, quote=FALSE, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    mytit <- paste0(ds," - ", rank_var, " - ", nTop)
    all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nTop == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
    
    outFile <- file.path(outFold, paste0(ds,  "_", rank_var, "_nTop", nTop, ".", plotType ))
    
    ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
    cat("... written: ", outFile, "\n")
  }
}


# ######################################################################################
# ######################################################################################
# ######################################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))





