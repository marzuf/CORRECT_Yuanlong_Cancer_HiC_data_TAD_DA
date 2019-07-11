# Rscript nGenes_all_values.R

script_name <- "nGenes_all_values.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

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

outFold <- "NGENES_ALL_VALUES"
dir.create(outFold, recursive = TRUE)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

windowSizeBp <- 500*10^3
options(scipen=100)

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

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
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  tad_mC <- eval(parse(text = load(curr_file)))
  notadj <- as.numeric(tad_mC)
  adj <- p.adjust(notadj, method="BH")
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    comb_empPval = notadj,
    comb_empPval_adj = adj,
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
  notadj <- as.numeric(tad_mC)
  adj <- p.adjust(notadj, method="BH")
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    combRank_empPval = notadj,
    combRank_empPval_adj = adj,
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


pvals_all_DT <- all_DT

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanLogFC.Rdata", full.names = TRUE)
stopifnot(length(all_fc_files) > 0)

all_corr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanCorr.Rdata", full.names = TRUE)
stopifnot(length(all_corr_files) > 0)

tieMet="min"

mC_DT <- foreach(corr_file = all_corr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    meanCorr = as.numeric(tad_mC),
    meanCorr_rank = rank(-as.numeric(tad_mC), ties = tieMet),
    stringsAsFactors = FALSE
  )
}

fc_DT <- foreach(corr_file = all_fc_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  exprds <- basename(dirname(dirname(corr_file)))
  hicds <- basename(dirname(dirname(dirname(corr_file))))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mC),
    meanFC = as.numeric(tad_mC),
    meanFC_rank = rank(-abs(as.numeric(tad_mC)), ties = tieMet),
    stringsAsFactors = FALSE
  )
}

stopifnot(nrow(fc_DT) == nrow(mC_DT) )


all_DT <- merge(fc_DT, mC_DT,   by=c("hicds", "exprds", "region"))

stopifnot(nrow(all_DT) == nrow(mC_DT) )
stopifnot(nrow(all_DT) == nrow(fc_DT) )

values_all_DT <- all_DT

stopifnot(nrow(values_all_DT) == nrow(pvals_all_DT) )

all_DT <- merge(values_all_DT, pvals_all_DT, by =c("hicds", "exprds", "region"), all=TRUE)

stopifnot(nrow(all_DT) == nrow(values_all_DT) )
stopifnot(nrow(all_DT) == nrow(pvals_all_DT) )


all_DT$dataset <- paste0(all_DT$exprds, "_", all_DT$hicds)

outFile <- file.path(outFold, "all_DT.Rdata")
save(all_DT, file=outFile)
cat("... written: ", outFile, "\n")



all_DT$dspath <- file.path(all_DT$hicds, all_DT$exprds)

script0_name <- "0_prepGeneData"

pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER")

toRetrieve <- unique(all_DT$dspath)
curr_ds <- toRetrieve[1]
nGenes_DT <- foreach(curr_ds = toRetrieve, .combine='rbind') %dopar% {
  
  exprds <- basename(curr_ds)
  hicds <- dirname(curr_ds)
  
  geneListFile <- file.path(pipDir, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneListFile))
  geneList <- eval(parse(text = load(geneListFile)))
  
  
  g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(g2tFile))
  gene2tadDT <- read.delim(g2tFile, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
  gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
  stopifnot(geneList %in% gene2tadDT$entrezID)
  gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
  
  tadSize  <- setNames(as.numeric(table(gene2tadDT$region)), names((table(gene2tadDT$region))))
  
  data.frame(
    hicds=hicds,
    exprds=exprds,
    region=names(tadSize),
    nGenes=as.numeric(tadSize),
    stringsAsFactors = FALSE
  )  
  
}



all_DT_withSize <- merge(all_DT, nGenes_DT, by =c("hicds", "exprds", "region"), all.x = TRUE, all.y=FALSE)
stopifnot(nrow(all_DT_withSize) == nrow(all_DT))
stopifnot(!is.na(all_DT_withSize))

outFile <- file.path(outFold, "all_DT_withSize.Rdata")
save(all_DT_withSize, file=outFile)
cat("... written: ", outFile, "\n")


all_DT_withSize$comb_empPval_adj_log10 <- -log10(all_DT_withSize$comb_empPval_adj)
all_DT_withSize$combRank_empPval_adj_log10 <- -log10(all_DT_withSize$combRank_empPval_adj)

all_DT_withSize$meanRank <- 0.5*(all_DT_withSize$meanFC_rank + all_DT_withSize$meanCorr_rank)

############################################################################
############################################################################

all_plot_vars <- c("comb_empPval_adj_log10", "combRank_empPval_adj_log10", "meanRank")

xvar="nGenes"

nDS <- length(unique(all_DT_withSize$dataset))

for(plot_var in all_plot_vars) {
    
  yvar=paste0(plot_var)
  
  myx <- all_DT_withSize[,xvar]
  myy <- all_DT_withSize[,yvar]
  outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, "_allDS.", plotType ))
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
  mtext(side=3, text=paste0("nDS = ", nDS), font=3)
  foo <- dev.off()
  cat("... written: ", outFile, "\n")
  
    
    
  }













# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))










