# Rscript coexpr_DE_queryTAD_byDS.R

script_name <- "coexpr_DE_queryTAD_byDS.R"

startTime <- Sys.time()

cat("> START coexpr_DE_queryTAD_byDS.R \n")

SSHFS <- FALSE

setDir <- ifelse(SSHFS, "/media/electron", "")

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("subtype_cols.R")

source("look_tads_fct.R")

source("plot_lolliTAD_funct.R")

nTop <- 5

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

windowSizeBp <- 500*10^3
options(scipen=100)


outFolder <- paste0("COEXPR_DE_QUERYTAD_BYDS_", nTop)
dir.create(outFolder, recursive=TRUE)


corrMet <- "pearson"
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)



all_ratioDown_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_ratioDown.Rdata", full.names = FALSE)
stopifnot(length(all_ratioDown_files) > 0)

all_meanCorr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanCorr_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_meanCorr_files) > 0)


script0_name <- "0_prepGeneData"
  

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)


### BUILD THE meanCorr TABLE
corr_file = all_meanCorr_files[1]
mC_DT <- foreach(corr_file = all_meanCorr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,corr_file)
  stopifnot(file.exists(curr_file))
  tad_mC <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(corr_file))
  data.frame(
    dataset = dataset,
    region = names(tad_mC),
    meanCorr = as.numeric(tad_mC),
    stringsAsFactors = FALSE
  )
}

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

tad_coexpr_DT <- mC_DT

fc_rd_DT <- merge(rD_DT, fc_DT, by=c("dataset", "region"))
stopifnot(nrow(fc_rd_DT) == nrow(rD_DT))
stopifnot(nrow(fc_DT) == nrow(rD_DT))
stopifnot(!is.na(fc_rd_DT))

stopifnot(fc_rd_DT$dataset %in% tad_coexpr_DT$dataset)
stopifnot(tad_coexpr_DT$dataset %in% fc_rd_DT$dataset)

tad_coexpr_fc_DT <- merge(tad_coexpr_DT, fc_rd_DT, by=c("dataset", "region"))
tad_coexpr_fc_DT <- tad_coexpr_fc_DT[order(tad_coexpr_fc_DT$meanCorr, decreasing = TRUE),]
tad_coexpr_fc_DT$TADrank <- 1:nrow(tad_coexpr_fc_DT)




tad_coexpr_fc_DT <- do.call(rbind,  by(tad_coexpr_fc_DT, tad_coexpr_fc_DT$dataset, function(subDT) {
  subDT$meanCorr_dsRank <- rank(- subDT$meanCorr, ties = "min")
  subDT$meanFC_dsRank <- rank(- subDT$meanFC, ties = "min")
  subDT$avg_dsRank <- 0.5*(subDT$meanCorr_dsRank + subDT$meanFC_dsRank)
  subDT
}))



tad_coexpr_fc_DT$meanFC_rank <- rank( - abs(tad_coexpr_fc_DT$meanFC) , ties ="min")
tad_coexpr_fc_DT$meanCorr_rank <- rank( - (tad_coexpr_fc_DT$meanCorr) , ties ="min")
tad_coexpr_fc_DT$avg_rank <- 0.5 * (tad_coexpr_fc_DT$meanFC_rank + tad_coexpr_fc_DT$meanCorr_rank)

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
# load(outFile)

myx = tad_coexpr_fc_DT$ratioDown
myy = tad_coexpr_fc_DT$meanFC_dsRank

outFile <- file.path(outFolder, paste0("meanFC_dsRank", "_vs_", "ratioDown_densplot.", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
    y = myy,
    x = myx,
    cex.lab=myCexLab,
    cex.axis=myCexAxis,
    ylab="meanFC_dsRank",
    xlab = "ratioDown",
    cex = 0.5,
    main = paste0("meanFC_dsRank vs ratioDown")
)
foo <- dev.off()
cat("... written: ", outFile, "\n")

myx = tad_coexpr_fc_DT$ratioDown
myy <- tad_coexpr_fc_DT$meanCorr_dsRank
outFile <- file.path(outFolder, paste0("meanCorr_dsRank", "_vs_", "ratioDown_densplot.", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
    y = myy,
    x = myx,
    cex.lab=myCexLab,
    cex.axis=myCexAxis,
    ylab="meanCorr_dsRank",
    xlab = "ratioDown",
    cex = 0.5,
    main = paste0("meanCorr_dsRank vs ratioDown")
)
foo <- dev.off()
cat("... written: ", outFile, "\n")


# stop("--ok\n")

#################################################################
################################################################# HIGH FC and HIGH WITHIN COEXPR - all datasets
#################################################################

# iterate
# subset_vars (subsetVar) -> subset_levels (var) -> rankingVars (rank_var)

subset_vars <- c("cmpType", "all")#, "dataset_plot")
# cmpType = subtypes/norm vs tumor/wt vs mut
rankingVars <- c( "meanFC_rank", "meanCorr_rank", "avg_rank",
                  "meanFC_dsRank", "meanCorr_dsRank", "avg_dsRank")


#### > START ITERATING AND OUTPUTING

tad_coexpr_fc_DT$dataset_plot <- gsub("/", "_", tad_coexpr_fc_DT$dataset) # otherwise issue to save file


subsetVar = subset_vars[1]
for(subsetVar in subset_vars) {
  if(subsetVar == "all") {
    subset_levels <- ""
  } else {
    subset_levels <- unique(tad_coexpr_fc_DT[, subsetVar])  
  }
  var = subset_levels[1]
  for(var in subset_levels) {
    if(subsetVar == "all") {
      stopifnot(var == "")
      curr_DT <- tad_coexpr_fc_DT
    } else {
        curr_DT <- tad_coexpr_fc_DT[tad_coexpr_fc_DT[, paste0(subsetVar)] == var,]
    }
    stopifnot(nrow(curr_DT) > 0)
    rank_var = rankingVars[1]
    for(rank_var in rankingVars){
      stopifnot(rank_var %in% colnames(curr_DT))
      ranked_curr_DT <- curr_DT[order(curr_DT[,rank_var]),]
      
      
      outFile <- file.path(outFolder, paste0(subsetVar, "_", var, "_", rank_var, "_nTop", nTop, "_densplot.", plotType ))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      densplot(
        y = curr_DT$meanCorr,
        x = curr_DT$meanFC,
        cex.lab=myCexLab,
        cex.axis=myCexAxis,
        ylab="meanCorr",
        xlab = "meanFC",
        cex = 0.5,
        main = paste0(subsetVar, " - ", var, " - ", rank_var)
      )
      points(x = ranked_curr_DT$meanFC[1:nTop],
             y = ranked_curr_DT$meanCorr[1:nTop],
            pch = 4,
            col="red",
            cex=1.5
             )
      text(x = ranked_curr_DT$meanFC[1:nTop],
             y = ranked_curr_DT$meanCorr[1:nTop],
           labels = paste0(ranked_curr_DT$dataset[1:nTop],"\n", ranked_curr_DT$region[1:nTop]),
             pch = 4,
             col="red",
             cex=0.7
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      i=1
      # topTable_DT <- foreach(i = 1:nTop, .combine='rbind') %dopar% {
      plotList <- list()
        topTable_DT <- foreach(i = 1:nTop, .combine='rbind') %do% {
        i_hicds <- dirname(ranked_curr_DT$dataset[i])
        i_exprds <- basename(ranked_curr_DT$dataset[i])
        i_tad <- basename(ranked_curr_DT$region[i])
        
        
        plotList[[i]] <- plot_lolliTAD_ds(exprds = i_exprds,
                         hicds = i_hicds,
                         all_TADs = i_tad,
                         orderByLolli = "startPos")
        
        
        tmpDT <- get_tad_symbols_DT(hicds = i_hicds,
                           exprds = i_exprds,
                           region = i_tad,
                           mainPipFolder = pipOutFolder
                           )
        otherCols <- colnames(tmpDT)
        tmpDT$iTop <- i
        tmpDT[,c("iTop", otherCols)]
      } # end-foreach iterating over the "n" top TADs
      # subset_vars (subsetVar) -> subset_levels (var) -> rankingVars (rank_var)
      outFile <- file.path(outFolder, paste0(subsetVar, "_", var, "_", rank_var, "_nTop", nTop, ".txt" ))
      cat(outFile, "\n")
      write.table(topTable_DT, col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE, quote=FALSE, file = outFile)
      cat(paste0("... written: ", outFile, "\n"))
      
      mytit <- paste0(subsetVar, " - ", var, " - ", rank_var, " - ", nTop)
      all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nTop == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      

      outFile <- file.path(outFolder, paste0(subsetVar, "_", var, "_", rank_var, "_nTop", nTop, ".", plotType ))
      
      outWidthGG <- 20
      outHeightGG <- min(c(7 * nTop/2, 49))
      
      ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      cat("... written: ", outFile, "\n")
      
    } # end-for iterating over rankingVar # => rank_var in rankingVars (e.g. "meanCorr_rank")
  } # end-for iterating over levels of data subset #   for(var in subset_levels) (e.g. "subtypes")
} # end-for iterating over subset data # for(subsetVar in subset_vars) (e.g. "cmpType")

#################################################################
################################################################# HIGH FC and HIGH WITHIN COEXPR - all datasets separately
#################################################################

# iterate
# subset_vars (subsetVar) -> subset_levels (var) -> rankingVars (rank_var)

subset_vars <- c("cmpType", "all")#, "dataset_plot")
# cmpType = subtypes/norm vs tumor/wt vs mut
rankingVars <- c( "meanFC_rank", "meanCorr_rank", "avg_rank",
                  "meanFC_dsRank", "meanCorr_dsRank", "avg_dsRank")


#### > START ITERATING AND OUTPUTING

tad_coexpr_fc_DT$dataset_plot <- gsub("/", "_", tad_coexpr_fc_DT$dataset) # otherwise issue to save file

all_ds <-  unique(tad_coexpr_fc_DT$dataset_plot)



  for(ds in all_ds) {
    
    curr_DT <- tad_coexpr_fc_DT[tad_coexpr_fc_DT$dataset_plot == ds,]
 
    
    stopifnot(nrow(curr_DT) > 0)
    rank_var = rankingVars[1]
    for(rank_var in rankingVars){
      stopifnot(rank_var %in% colnames(curr_DT))
      ranked_curr_DT <- curr_DT[order(curr_DT[,rank_var]),]
      
      
      outFile <- file.path(outFolder, paste0(ds, "_", rank_var, "_nTop", nTop, "_densplot.", plotType ))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      densplot(
        y = curr_DT$meanCorr,
        x = curr_DT$meanFC,
        cex.lab=myCexLab,
        cex.axis=myCexAxis,
        ylab="meanCorr",
        xlab = "meanFC",
        cex = 0.5,
        main = paste0(ds, " - ", rank_var)
      )
      points(x = ranked_curr_DT$meanFC[1:nTop],
             y = ranked_curr_DT$meanCorr[1:nTop],
             pch = 4,
             col="red",
             cex=1.5
      )
      text(x = ranked_curr_DT$meanFC[1:nTop],
           y = ranked_curr_DT$meanCorr[1:nTop],
           labels = paste0(ranked_curr_DT$dataset[1:nTop],"\n", ranked_curr_DT$region[1:nTop]),
           pch = 4,
           col="red",
           cex=0.7
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      i=1
      # topTable_DT <- foreach(i = 1:nTop, .combine='rbind') %dopar% {
      plotList <- list()
      topTable_DT <- foreach(i = 1:nTop, .combine='rbind') %do% {
        i_hicds <- dirname(ranked_curr_DT$dataset[i])
        i_exprds <- basename(ranked_curr_DT$dataset[i])
        i_tad <- basename(ranked_curr_DT$region[i])
        
        
        plotList[[i]] <- plot_lolliTAD_ds(exprds = i_exprds,
                                          hicds = i_hicds,
                                          all_TADs = i_tad,
                                          orderByLolli = "startPos")
        
        
        tmpDT <- get_tad_symbols_DT(hicds = i_hicds,
                                    exprds = i_exprds,
                                    region = i_tad,
                                    mainPipFolder = pipOutFolder
        )
        otherCols <- colnames(tmpDT)
        tmpDT$iTop <- i
        tmpDT[,c("iTop", otherCols)]
      } # end-foreach iterating over the "n" top TADs
      # subset_vars (subsetVar) -> subset_levels (var) -> rankingVars (rank_var)
      outFile <- file.path(outFolder, paste0(ds, "_", rank_var, "_nTop", nTop, ".txt" ))
      cat(outFile, "\n")
      write.table(topTable_DT, col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE, quote=FALSE, file = outFile)
      cat(paste0("... written: ", outFile, "\n"))
      
      mytit <- paste0(ds," - ", rank_var, " - ", nTop)
      all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nTop == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      
      
      outFile <- file.path(outFolder, paste0(ds,  "_", rank_var, "_nTop", nTop, ".", plotType ))
      
      outWidthGG <- 20
      outHeightGG <- min(c(7 * nTop/2, 49))
      
      ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      cat("... written: ", outFile, "\n")
      
    } # end-for iterating over variables # => rank_var in rankingVars (e.g. "meanCorr_rank")
} # end iterating over datasets






# ######################################################################################
# ######################################################################################
# ######################################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))





