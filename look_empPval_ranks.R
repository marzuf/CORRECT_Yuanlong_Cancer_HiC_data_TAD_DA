# Rscript look_empPval_ranks.R

script_name <- "look_empPval_ranks.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 90))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")

outFold <- "LOOK_EMPPVAL_RANKS"
dir.create(outFold)

plotType <- "png"
myWidth <- ifelse(plotType == "png", 400, 7)
myHeight <- myWidth

all_rankFC_files <- list.files("EMPPVAL_MEANTADFC_RANK", pattern="_meanLogFC_", full.names = TRUE, recursive = TRUE)
stopifnot(length(all_rankFC_files) > 0)

all_rankCorr_files <- list.files("EMPPVAL_MEANTADCORR_RANK", pattern="_meanCorr_", full.names = TRUE, recursive = TRUE)
stopifnot(length(all_rankCorr_files) > 0)

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_pvalComb_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = TRUE)
stopifnot(length(all_pvalComb_files) > 0)

all_pvalFC_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanLogFC.Rdata", full.names = TRUE)
stopifnot(length(all_pvalFC_files) > 0)

all_pvalCorr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanCorr.Rdata", full.names = TRUE)
stopifnot(length(all_pvalCorr_files) > 0)


### BUILD THE LOGFC TABLE
fc_file = all_pvalFC_files[1]
fc_DT <- foreach(fc_file = all_pvalFC_files, .combine = 'rbind') %dopar% {
  tad_fc <- eval(parse(text = load(fc_file)))
  hicds <- basename(dirname(dirname(dirname(fc_file))))
  exprds <- basename(dirname(dirname(fc_file)))
  data.frame(
    hicds = hicds,
    exprds = exprds, 
    region = names(tad_fc),
    meanFC_empPval = as.numeric(tad_fc),
    stringsAsFactors = FALSE
  )
}

### BUILD THE MEANCORR TABLE
meanCorr_file = all_pvalCorr_files[1]
meanCorr_DT <- foreach(meanCorr_file = all_pvalCorr_files, .combine = 'rbind') %dopar% {
  tad_meanCorr <- eval(parse(text = load(meanCorr_file)))
  hicds <- basename(dirname(dirname(dirname(meanCorr_file))))
  exprds <- basename(dirname(dirname(meanCorr_file)))
  data.frame(
    hicds = hicds,
    exprds = exprds, 
    region = names(tad_meanCorr),
    meanCorr_empPval = as.numeric(tad_meanCorr),
    stringsAsFactors = FALSE
  )
}

### BUILD THE PVAL COMBINED TABLE
pvalcomb_file = all_pvalComb_files[1]
pvalComb_DT <- foreach(pvalcomb_file = all_pvalComb_files, .combine = 'rbind') %dopar% {
  tad_pvalComb <- eval(parse(text = load(pvalcomb_file)))
  hicds <- basename(dirname(dirname(dirname(pvalcomb_file))))
  exprds <- basename(dirname(dirname(pvalcomb_file)))
  adj_pvalComb <- p.adjust(tad_pvalComb, method="BH")
  stopifnot(setequal(names(tad_pvalComb), names(adj_pvalComb)))
  data.frame(
    hicds = hicds,
    exprds = exprds, 
    region = names(tad_pvalComb),
    empPval_comb = as.numeric(tad_pvalComb),
    empPval_comb_adj = as.numeric(adj_pvalComb[names(tad_pvalComb)]),
    stringsAsFactors = FALSE
  )
}

### BUILD THE LOGFC RANK EMPPVAL TABLE
fc_file = all_rankFC_files[1]
fcRankEmpPval_DT <- foreach(fc_file = all_rankFC_files, .combine = 'rbind') %dopar% {
  stopifnot(file.exists(fc_file))
  tad_fcRank <- eval(parse(text = load(fc_file)))
  hicds <- basename(dirname(dirname(fc_file)))
  exprds <- basename(dirname(fc_file))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_fcRank),
    meanFC_rank_empPval = as.numeric(tad_fcRank),
    stringsAsFactors = FALSE
  )
}
stopifnot(fcRankEmpPval_DT$meanFC_rank_empPval > 0)
stopifnot(fcRankEmpPval_DT$meanFC_rank_empPval <= 1)


### BUILD THE MEANCORR RANK EMPPVAL TABLE
mc_file = all_rankCorr_files[1]
mcRankEmpPval_DT <- foreach(mc_file = all_rankCorr_files, .combine = 'rbind') %dopar% {
  stopifnot(file.exists(mc_file))
  tad_mcRank <- eval(parse(text = load(mc_file)))
  hicds <- basename(dirname(dirname(mc_file)))
  exprds <- basename(dirname(mc_file))
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = names(tad_mcRank),
    meanCorr_rank_empPval = as.numeric(tad_mcRank),
    stringsAsFactors = FALSE
  )
}
stopifnot(mcRankEmpPval_DT$meanCorr_rank_empPval > 0)
stopifnot(mcRankEmpPval_DT$meanCorr_rank_empPval <= 1)



stopifnot(nrow(pvalComb_DT) == nrow(mcRankEmpPval_DT) )
stopifnot(nrow(fc_DT) == nrow(mcRankEmpPval_DT) )
stopifnot(nrow(meanCorr_DT) == nrow(mcRankEmpPval_DT) )


stopifnot(nrow(fcRankEmpPval_DT) == nrow(mcRankEmpPval_DT) )
rankEmpPval_DT <- merge(fcRankEmpPval_DT, mcRankEmpPval_DT, by =c("hicds", "exprds", "region"), all = TRUE)

outFile <- file.path(outFold, "rankEmpPval_DT.Rdata")
save(rankEmpPval_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))



stopifnot(!is.na(rankEmpPval_DT))

pvalCombRank <- apply(rankEmpPval_DT, 1, function(x) {
  stouffer(c(as.numeric(as.character(x["meanCorr_rank_empPval"])), as.numeric(as.character(x["meanFC_rank_empPval"]))), two.tails=TRUE)
})

rankEmpPval_DT$rank_empPval_comb <- pvalCombRank
rankEmpPval_DT$rank_empPval_comb_log10 <- -log10(rankEmpPval_DT$rank_empPval_comb)

rankEmpPval_DT$rank_empPval_comb_adj <- p.adjust(rankEmpPval_DT$rank_empPval_comb, method="BH")
rankEmpPval_DT$rank_empPval_comb_adj_log10 <- -log10(rankEmpPval_DT$rank_empPval_comb_adj)

rankEmpPval_DT$meanCorr_rank_empPval_log10 <- -log10(rankEmpPval_DT$meanCorr_rank_empPval)
rankEmpPval_DT$meanFC_rank_empPval_log10 <- -log10(rankEmpPval_DT$meanFC_rank_empPval)

nDS <- length(unique(paste0(rankEmpPval_DT$hicds, "_",rankEmpPval_DT$exprds )))


empPval_DT <- merge(pvalComb_DT, fc_DT, by =c("hicds", "exprds", "region"), all = TRUE)
empPval_DT <- merge(empPval_DT, meanCorr_DT, by =c("hicds", "exprds", "region"), all = TRUE)

outFile <- file.path(outFold, "empPval_DT.Rdata")
save(empPval_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


stopifnot(nrow(empPval_DT) == nrow(pvalComb_DT) )
stopifnot(nrow(empPval_DT) == nrow(fc_DT) )
stopifnot(nrow(empPval_DT) == nrow(meanCorr_DT) )

stopifnot(!is.na(empPval_DT))

stopifnot(nrow(empPval_DT) == nrow(rankEmpPval_DT) )

all_empPval_DT <- merge(empPval_DT, rankEmpPval_DT, by =c("hicds", "exprds", "region"), all = TRUE)




all_empPval_DT$empPval_comb_adj_log10 <- -log10(all_empPval_DT$empPval_comb_adj)


all_empPval_DT$meanCorr_empPval_log10 <- -log10(all_empPval_DT$meanCorr_empPval)
all_empPval_DT$meanFC_empPval_log10 <- -log10(all_empPval_DT$meanFC_empPval)

outFile <- file.path(outFold, "all_empPval_DT.Rdata")
save(all_empPval_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))



##########################################################################################################
##########################################################################################################
##########################################################################################################

x_val <- "meanCorr_rank_empPval_log10"
y_val <- "meanFC_rank_empPval_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##########################################################################################################
##########################################################################################################
##########################################################################################################

x_val <- "meanCorr_empPval_log10"
y_val <- "meanFC_empPval_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


##########################################################################################################
##########################################################################################################
##########################################################################################################


x_val <- "meanFC_empPval_log10"
y_val <- "meanFC_rank_empPval_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##########################################################################################################
##########################################################################################################
##########################################################################################################

x_val <- "meanCorr_empPval_log10"
y_val <- "meanCorr_rank_empPval_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


##########################################################################################################
##########################################################################################################                       EEEEEEEEERRRRRRRRRRRRRRRROOOOOOOOOOORRRRRRRRRRr
##########################################################################################################


x_val <- "empPval_comb_adj_log10"
y_val <- "rank_empPval_comb_adj_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


##########################################################################################################
##########################################################################################################
##########################################################################################################

x_val <- "meanCorr_rank_empPval_log10"
y_val <- "rank_empPval_comb_adj_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



##########################################################################################################
##########################################################################################################
##########################################################################################################

x_val <- "meanFC_rank_empPval_log10"
y_val <- "rank_empPval_comb_adj_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##########################################################################################################
##########################################################################################################
##########################################################################################################

x_val <- "meanCorr_empPval_log10"
y_val <- "empPval_comb_adj_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



##########################################################################################################
##########################################################################################################
##########################################################################################################

x_val <- "meanFC_empPval_log10"
y_val <- "empPval_comb_adj_log10"
myx <- all_empPval_DT[,x_val]
myy <- all_empPval_DT[,y_val]
outFile <- file.path(outFold, paste0(y_val, "_vs_",x_val, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  xlab = paste0(x_val),
  ylab = paste0(y_val),
  main = paste0(y_val, " vs. ", x_val)
)
addCorr(x=myx, y=myy, legPos="topleft")
mtext(side = 3, text = paste0("nDS = ", nDS), font = 3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


# ######################################################################################
# ######################################################################################
# ######################################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))






