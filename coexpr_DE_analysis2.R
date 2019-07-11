# Rscript coexpr_DE_analysis2.R

script_name <- "coexpr_DE_analysis2.R"

startTime <- Sys.time()

cat("> START coexpr_DE_analysis2.R \n")

SSHFS <- FALSE

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("subtype_cols.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

windowSizeBp <- 500*10^3 
options(scipen=100)


outFolder <- "COEXPR_DE_ANALYSIS2"
dir.create(outFolder, recursive=TRUE)


corrMet <- "pearson"
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)


all_meanCorr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanCorr_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_meanCorr_files) > 0)

all_meanCorrCond1_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanCorr_TAD_cond1.Rdata", full.names = FALSE)
stopifnot(length(all_meanCorrCond1_files) > 0)


all_meanCorrCond2_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanCorr_TAD_cond2.Rdata", full.names = FALSE)
stopifnot(length(all_meanCorrCond2_files) > 0)

all_pvalcomb_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = FALSE)
stopifnot(length(all_pvalcomb_files) > 0)


all_ratioDown_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_ratioDown.Rdata", full.names = FALSE)
stopifnot(length(all_ratioDown_files) > 0)

stopifnot(length(all_ratioDown_files) == length(all_fc_files))
stopifnot(length(all_fc_files) == length(all_pvalcomb_files) )
  
# sort the TADs by decreasing withinCoexpr
# plot level of coexpr within and between on the same plot


### BUILD THE meanCorr TABLE
rd_file = all_meanCorrCond1_files[1]
mCcond1_DT <- foreach(rd_file = all_meanCorrCond1_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,rd_file)
  stopifnot(file.exists(curr_file))
  tad_rd <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(rd_file))
  data.frame(
    dataset = dataset,
    region = names(tad_rd),
    withinCoexpr_cond1 = as.numeric(tad_rd),
    stringsAsFactors = FALSE
  )
}


### BUILD THE meanCorr TABLE
rd_file = all_meanCorrCond2_files[1]
mCcond2_DT <- foreach(rd_file = all_meanCorrCond2_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,rd_file)
  stopifnot(file.exists(curr_file))
  tad_rd <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(rd_file))
  data.frame(
    dataset = dataset,
    region = names(tad_rd),
    withinCoexpr_cond2 = as.numeric(tad_rd),
    stringsAsFactors = FALSE
  )
}

### BUILD THE meanCorr TABLE
rd_file = all_meanCorr_files[1]
mC_DT <- foreach(rd_file = all_meanCorr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,rd_file)
  stopifnot(file.exists(curr_file))
  tad_rd <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(rd_file))
  data.frame(
    dataset = dataset,
    region = names(tad_rd),
    withinCoexpr = as.numeric(tad_rd),
    stringsAsFactors = FALSE
  )
}

tad_coexpr_DT <- merge(mC_DT, mCcond1_DT, by=c("dataset", "region"))
stopifnot(nrow(tad_coexpr_DT) == nrow(mC_DT))
stopifnot(nrow(mC_DT) == nrow(mCcond1_DT))
stopifnot(!is.na(tad_coexpr_DT))

tad_coexpr_DT <- merge(tad_coexpr_DT, mCcond2_DT, by=c("dataset", "region"))
stopifnot(nrow(tad_coexpr_DT) == nrow(mCcond2_DT))
stopifnot(!is.na(tad_coexpr_DT))



### BUILD THE ratio down TABLE
rd_file = all_ratioDown_files[1]
rD_DT <- foreach(rd_file = all_ratioDown_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,rd_file)
  stopifnot(file.exists(curr_file))
  tad_rd <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(rd_file))
  data.frame(
    dataset = dataset,
    region = names(tad_rd),
    ratioDown = as.numeric(tad_rd),
    stringsAsFactors = FALSE
  )
}

#### BUILD THE CPTMT SCORE TABLE
#score_file = all_domainScore_files[1]
#score_DT <- foreach(score_file = all_domainScore_files, .combine = 'rbind') %dopar% {
#  curr_file <- file.path(score_file)
#  stopifnot(file.exists(curr_file))
#  curr_DT <- read.delim(curr_file, header=F, 
#                        col.names = c("chromo", "start", "end", "region", "score"))
#  curr_DT$dataset <- dirname(dirname(curr_file))
#  curr_DT
#}

### BUILD THE LOGFC TABLE
fc_file = all_fc_files[1]
fc_DT <- foreach(fc_file = all_fc_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, fc_file)
  stopifnot(file.exists(curr_file))
  tad_fc <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(fc_file))
  data.frame(
    dataset = dataset,
    region = names(tad_fc),
    meanFC = as.numeric(tad_fc),
    stringsAsFactors = FALSE
  )
}


### BUILD THE LOGFC TABLE
pvalcomb_file = all_pvalcomb_files[1]
pvalcomb_DT <- foreach(pvalcomb_file = all_pvalcomb_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, pvalcomb_file)
  stopifnot(file.exists(curr_file))
  tad_pvalcomb <- eval(parse(text = load(curr_file)))
  # adj pval comb
  adj_tad_pvalcomb <- p.adjust(tad_pvalcomb, method="BH")
  dataset <- dirname(dirname(pvalcomb_file))
  data.frame(
    dataset = dataset,
    region = names(adj_tad_pvalcomb),
    adjPvalComb = as.numeric(adj_tad_pvalcomb),
    stringsAsFactors = FALSE
  )
}



############################################################################## MERGE THE TABLES

fc_adjpval_DT <- merge(pvalcomb_DT, fc_DT, by=c("dataset", "region"))
stopifnot(nrow(fc_adjpval_DT) == nrow(fc_DT))
stopifnot(nrow(fc_DT) == nrow(pvalcomb_DT))
stopifnot(!is.na(fc_adjpval_DT))

fc_rd_DT <- merge(rD_DT, fc_adjpval_DT, by=c("dataset", "region"))
stopifnot(nrow(fc_rd_DT) == nrow(rD_DT))
stopifnot(nrow(fc_DT) == nrow(rD_DT))
stopifnot(!is.na(fc_rd_DT))

stopifnot(fc_rd_DT$dataset %in% tad_coexpr_DT$dataset)
stopifnot(tad_coexpr_DT$dataset %in% fc_rd_DT$dataset)

tad_coexpr_fc_DT <- merge(tad_coexpr_DT, fc_rd_DT, by=c("dataset", "region"))
tad_coexpr_fc_DT <- tad_coexpr_fc_DT[order(tad_coexpr_fc_DT$withinCoexpr, decreasing = TRUE),]
tad_coexpr_fc_DT$TADrank <- 1:nrow(tad_coexpr_fc_DT)


tad_coexpr_fc_DT$withinDiffCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 - tad_coexpr_fc_DT$withinCoexpr_cond2)
tad_coexpr_fc_DT$withinRatioCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 / tad_coexpr_fc_DT$withinCoexpr_cond2)
tad_coexpr_fc_DT$withinChangeratioCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond2 - tad_coexpr_fc_DT$withinCoexpr_cond1)/tad_coexpr_fc_DT$withinCoexpr_cond1







# select only that with + coexpr in both condition -> I can take logFC



tad_coexpr_fc_DT$withinCond2WithinCond1LogFC <- log10(tad_coexpr_fc_DT$withinCoexpr_cond2/tad_coexpr_fc_DT$withinCoexpr_cond1)

tad_coexpr_fc_DT$adjPvalComb_log10 <- -log10(tad_coexpr_fc_DT$adjPvalComb)


# for each, plot i) densplot; ii) plot with color for subtypes
tad_coexpr_fc_DT$cmps <- basename(tad_coexpr_fc_DT$dataset)

colDT <- data.frame(
  cmps = names(all_cmps),
  cmpType = all_cmps,
stringsAsFactors = FALSE
)
colDT$cmpCol <- all_cols[colDT$cmpType]
stopifnot(!is.na(colDT))

stopifnot(tad_coexpr_fc_DT$cmps %in% colDT$cmps)

tad_coexpr_fc_DT <- merge(tad_coexpr_fc_DT, colDT, by = "cmps", all.x = TRUE, all.y = FALSE)

stopifnot(!is.na(tad_coexpr_fc_DT$cmpCol))


outFile <- file.path(outFolder, "tad_coexpr_fc_DT.Rdata")
save(tad_coexpr_fc_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


myplot_densplot <- function(xvar, yvar, addCurve=FALSE, dt = tad_coexpr_fc_DT, outPrefix="") {
  
  stopifnot(xvar %in% colnames(dt))
  myx <- dt[,xvar]
  stopifnot(yvar %in% colnames(dt))
  myy <- dt[,yvar]
  mycols <- dt$cmpCol
  
  outFile <- file.path(outFolder, paste0(outPrefix, "", yvar, "_vs_", xvar, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(y=myy,
           x=myx,
           pch=16,
           cex=0.7,
           cex.axis = myCexAxis,
           cex.lab = myCexLab,
           xlab=xvar,
           ylab=yvar,
           main = paste0(yvar, " vs. ", xvar)
  )
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  if(addCurve) {
    curve(1*x, lty=2, col="grey", add = TRUE)
  }
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

myplot_colplot <- function(xvar, yvar, mycols, addCurve = FALSE, dt = tad_coexpr_fc_DT, outPrefix="") {
  
  stopifnot(xvar %in% colnames(dt))
  myx <- dt[,xvar]
  stopifnot(yvar %in% colnames(dt))
  myy <- dt[,yvar]
  
  outFile <- file.path(outFolder, paste0(outPrefix, "", yvar, "_vs_", xvar, "_colplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(y=myy,
           x=myx,
           pch=16,
           cex=0.7,
           cex.axis = myCexAxis,
           cex.lab = myCexLab,
           xlab=xvar,
           ylab=yvar,
       col=mycols,
           main = paste0(yvar, " vs. ", xvar)
  )
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  if(addCurve) {
    curve(1*x, lty=2, col="grey",add=TRUE)
  }
  addSubtypeLeg(bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

mycols <- tad_coexpr_fc_DT$cmpCol

##########################
### adj. pval. combined and withinCoexpr diff cond1 cond2
##########################

yvar <- "adjPvalComb_log10"
xvar <- "withinDiffCond1Cond2"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)



##########################
### detect "disruption": those with change in coexpr => coexpr cond2 vs coexpr cond1
##########################

yvar <- "withinCoexpr_cond2"
xvar <- "withinCoexpr_cond1"

myplot_densplot(xvar,yvar, addCurve = TRUE)
myplot_colplot(xvar,yvar,mycols, addCurve = TRUE)

mycols_sub <- tad_coexpr_fc_DT$cmpCol[tad_coexpr_fc_DT$cmpType == "subtypes"]
mycols_tumor <- tad_coexpr_fc_DT$cmpCol[tad_coexpr_fc_DT$cmpType == "norm_vs_tumor"]
mycols_mut <- tad_coexpr_fc_DT$cmpCol[tad_coexpr_fc_DT$cmpType == "wt_vs_mut"]

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "subtypes",] , outPrefix = "subtypes_")
myplot_colplot(xvar,yvar,mycols_sub, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "subtypes",] , outPrefix = "subtypes_")

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "norm_vs_tumor",] , outPrefix = "norm_vs_tumor_")
myplot_colplot(xvar,yvar,mycols_tumor, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "norm_vs_tumor",] , outPrefix = "norm_vs_tumor_")

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "wt_vs_mut",] , outPrefix = "wt_vs_mut_")
myplot_colplot(xvar,yvar,mycols_mut, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "wt_vs_mut",] , outPrefix = "wt_vs_mut_")


##########################
### detect those that function as regulatory unit and DE: those with high FC and high diff. between-within => betw-within vs. FC expr
##########################

yvar <- "withinCoexpr"
xvar <- "meanFC"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)

yvar <- "withinCoexpr"
xvar <- "ratioDown"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)







# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





