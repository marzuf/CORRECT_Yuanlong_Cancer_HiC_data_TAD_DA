
# Rscript coexpr_between_within_oneDS.R ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR

script_name <- "coexpr_between_within_oneDS.R"

startTime <- Sys.time()

cat("> START coexpr_between_within_oneDS.R \n")

SSHFS <- FALSE

buildData <- TRUE

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 90))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

args <- commandArgs(trailingOnly = TRUE)


outFolder <- "COEXPR_BETWEEN_WITHIN_ONEDS"
dir.create(outFolder, recursive=TRUE)

corrMet <- "pearson"
script0_name <- "0_prepGeneData"
mypattern <- "coexprDT.Rdata"

coexprFiles  <- list.files("CREATE_COEXPR_SORTNODUP", recursive = TRUE, pattern = mypattern, full.names = TRUE)
stopifnot(length(coexprFiles) > 0)

coexprByCondFolder <- "CREATE_COEXPR_BYCOND_SORTNODUP"
stopifnot(dir.exists(coexprByCondFolder))

coexprDataFile = coexprFiles[1]
hicds = "K562_40kb"
exprds = "TCGAlaml_wt_mutFLT3"
coexprDataFile = file.path("CREATE_COEXPR_SORTNODUP", 
                            hicds, 
                            paste0(exprds),
                            corrMet,
                            "coexprDT.Rdata")



#### LOAD INFORMATION ABOUT GENES AROUND TAD (SAME NUMBER)
aroundTADfile <- file.path("CREATE_SAMPLE_AROUND_TADS", "all_ds_sample_around_TADs.Rdata")
stopifnot(file.exists(aroundTADfile))
all_ds_sample_around_TADs <- eval(parse(text = load(aroundTADfile)))


#### LOAD INFORMATION ABOUT GENES AROUND TAD (WINDOW)
windowSizeBp <- 500*10^3
options(scipen=100)
aroundKbTADfile <- file.path("CREATE_SAMPLE_AROUNDKB_TADS", paste0(windowSizeBp), "all_ds_sample_aroundKb_TADs.Rdata")
stopifnot(file.exists(aroundKbTADfile))
all_ds_sample_aroundKb_TADs <- eval(parse(text = load(aroundKbTADfile)))


coexprFiles <-  coexprFiles[grepl(args[1], coexprFiles) & grepl(args[2], coexprFiles)]

stopifnot(length(coexprFiles) == 1)

coexprFilesOutNames <- coexprFiles

coexprDataFile = coexprFiles[1]

if(buildData) {

  allData_within_between_coexpr <- foreach(coexprDataFile = coexprFiles) %do% {
    
    cat("... start ", coexprDataFile, " \n")
    
    stopifnot(file.exists(coexprDataFile))
    
    hicds <- basename(dirname(dirname(dirname(coexprDataFile))))
    stopifnot(dir.exists(hicds))
    exprds <- basename(dirname(dirname(coexprDataFile)))
    
    curr_coexprByCondFolder <- file.path(coexprByCondFolder, hicds, exprds, corrMet)
    stopifnot(dir.exists(curr_coexprByCondFolder))

    coexprDT_cond1File <- file.path(curr_coexprByCondFolder , "coexprDT_cond1.Rdata")
    stopifnot(file.exists(coexprDT_cond1File))
    coexprDT_cond2File <- file.path(curr_coexprByCondFolder , "coexprDT_cond2.Rdata")
    stopifnot(file.exists(coexprDT_cond2File))
    
    
    cat("... start ", hicds, " - ", exprds, "\n")
    
    cat("... load coexprDT ...\n")
    coexprDT <- eval(parse(text = load(coexprDataFile)))
    coexprDT$gene1 <- as.character(coexprDT$gene1)
    coexprDT$gene2 <- as.character(coexprDT$gene2)
    
    
    cat("... load coexprDT_cond1 ...\n")
    coexprDT_cond1 <- eval(parse(text = load(paste0(coexprDT_cond1File))))
    coexprDT_cond1$gene1 <- as.character(coexprDT_cond1$gene1)
    coexprDT_cond1$gene2 <- as.character(coexprDT_cond1$gene2)
    
    
    cat("... load coexprDT_cond2 ...\n")
    coexprDT_cond2 <- eval(parse(text = load(paste0(coexprDT_cond2File))))
    coexprDT_cond2$gene1 <- as.character(coexprDT_cond2$gene1)
    coexprDT_cond2$gene2 <- as.character(coexprDT_cond2$gene2)
    

    g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, 
                         col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    dsPipOutDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
    stopifnot(dir.exists(dsPipOutDir))
    stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
    geneListFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) # 
    
    regionListFile <- file.path(dsPipOutDir, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(regionListFile))
    pipeline_regionList <- eval(parse(text = load(regionListFile))) # 
    
    stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList & g2t_DT$region %in% pipeline_regionList, ]
    stopifnot(g2t_DT$entrezID %in% coexprDT$gene1 | g2t_DT$entrezID %in% coexprDT$gene2 )
    stopifnot(pipeline_regionList %in% g2t_DT$region)
    
    reg = pipeline_regionList[1]
    #pipeline_regionList <- pipeline_regionList[1:10]
    within_between_coexpr_data <- foreach(reg = pipeline_regionList) %dopar% {
      
      cat(paste0("...... start ", hicds, " - ", exprds, " - ", reg, " \n"))
      tad_genes <- g2t_DT$entrezID[g2t_DT$region == reg]
      
      ####################################################################
      #### ALL CONDITIONS
      ####################################################################
      tad_coexprDT <- coexprDT[coexprDT$gene1 %in% tad_genes | 
                                 coexprDT$gene2 %in% tad_genes,]
    
      stopifnot(nrow(tad_coexprDT) > 0)  
      stopifnot(!is.na(tad_coexprDT$coexpr))
      stopifnot(is.numeric(tad_coexprDT$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT$gene1 | tad_genes %in% tad_coexprDT$gene2 )
      
      ### WITHIN
      within_coexprDT <- tad_coexprDT[tad_coexprDT$gene1 %in% tad_genes & tad_coexprDT$gene2 %in% tad_genes,]
      
      ### BETWEEN: ALL
      # betweenAll: 1 of the 2 is not in TAD
      betweenAll_coexprDT <- tad_coexprDT[ ! (tad_coexprDT$gene1 %in% tad_genes & tad_coexprDT$gene2 %in% tad_genes),]
      stopifnot(nrow(within_coexprDT) + nrow(betweenAll_coexprDT) == nrow(tad_coexprDT))
      
      stopifnot(betweenAll_coexprDT$gene1 %in% tad_genes | 
                  betweenAll_coexprDT$gene2 %in% tad_genes )
      
      stopifnot( ! betweenAll_coexprDT$gene1 %in% tad_genes | 
                  ! betweenAll_coexprDT$gene2 %in% tad_genes )
      
      ### BETWEEN: GENES AROUND - NUMBER
      betweenGenes_windowNbr <- all_ds_sample_around_TADs[[file.path(hicds, exprds)]][[paste0(reg)]]
      stopifnot(length(betweenGenes_windowNbr) > 0)
      
      betweenNbr_coexprDT <- tad_coexprDT[ tad_coexprDT$gene1 %in% betweenGenes_windowNbr | 
                                             tad_coexprDT$gene2 %in% betweenGenes_windowNbr,] 
      stopifnot( betweenNbr_coexprDT$gene1 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT$gene1 %in% tad_genes )
      
      stopifnot( betweenNbr_coexprDT$gene2 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenNbr_coexprDT$gene1 %in% tad_genes & 
                      betweenNbr_coexprDT$gene2 %in% tad_genes) ) 
                      
      stopifnot(betweenGenes_windowNbr %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowNbr %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowNbr %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowNbr %in% within_coexprDT$gene2)
      
      ### BETWEEN: GENES AROUND - KB
      betweenGenes_windowKb <- all_ds_sample_aroundKb_TADs[[file.path(hicds, exprds)]][[paste0(reg)]]
      # stopifnot(length(betweenGenes_windowKb) > 0) => not always TRUE
      
      betweenKb_coexprDT <- tad_coexprDT[ tad_coexprDT$gene1 %in% betweenGenes_windowKb | 
                                             tad_coexprDT$gene2 %in% betweenGenes_windowKb,] 
      stopifnot( betweenKb_coexprDT$gene1 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT$gene1 %in% tad_genes )
      
      stopifnot( betweenKb_coexprDT$gene2 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenKb_coexprDT$gene1 %in% tad_genes & 
                      betweenKb_coexprDT$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowKb %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowKb %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowKb %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowKb %in% within_coexprDT$gene2)
      
      ####################################################################
      #### COND 1
      ####################################################################
      
      tad_coexprDT_cond1 <- coexprDT_cond1[coexprDT_cond1$gene1 %in% tad_genes | 
                                 coexprDT_cond1$gene2 %in% tad_genes,]
      
      stopifnot(nrow(tad_coexprDT_cond1) > 0)  
      stopifnot(!is.na(tad_coexprDT_cond1$coexpr))
      stopifnot(is.numeric(tad_coexprDT_cond1$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT_cond1$gene1 | tad_genes %in% tad_coexprDT_cond1$gene2 )
      
      ### WITHIN
      
      within_coexprDT_cond1 <- tad_coexprDT_cond1[tad_coexprDT_cond1$gene1 %in% tad_genes & tad_coexprDT_cond1$gene2 %in% tad_genes,]
      
      
      
      ### BETWEEN: ALL
      
      
      betweenAll_coexprDT_cond1 <- tad_coexprDT_cond1[ ! (tad_coexprDT_cond1$gene1 %in% tad_genes & tad_coexprDT_cond1$gene2 %in% tad_genes),]
      
      stopifnot(nrow(within_coexprDT_cond1) + nrow(betweenAll_coexprDT_cond1) == nrow(tad_coexprDT_cond1))
      
      
      stopifnot(betweenAll_coexprDT_cond1$gene1 %in% tad_genes | 
                  betweenAll_coexprDT_cond1$gene2 %in% tad_genes )
      
      stopifnot( ! betweenAll_coexprDT_cond1$gene1 %in% tad_genes | 
                   ! betweenAll_coexprDT_cond1$gene2 %in% tad_genes )
      
      
      
      
      ### BETWEEN: GENES AROUND - NBR
    
      
      
      betweenNbr_coexprDT_cond1 <- tad_coexprDT_cond1[ tad_coexprDT_cond1$gene1 %in% betweenGenes_windowNbr | 
                                             tad_coexprDT_cond1$gene2 %in% betweenGenes_windowNbr,] 
      stopifnot( betweenNbr_coexprDT_cond1$gene1 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT_cond1$gene1 %in% tad_genes )
      
      stopifnot( betweenNbr_coexprDT_cond1$gene2 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT_cond1$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenNbr_coexprDT_cond1$gene1 %in% tad_genes & 
                      betweenNbr_coexprDT_cond1$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowNbr %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowNbr %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowNbr %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowNbr %in% within_coexprDT$gene2)
      
      
      
      
      
      ### BETWEEN: GENES AROUND - KB
      
      betweenKb_coexprDT_cond1 <- tad_coexprDT_cond1[ tad_coexprDT_cond1$gene1 %in% betweenGenes_windowKb | 
                                            tad_coexprDT_cond1$gene2 %in% betweenGenes_windowKb,] 
      stopifnot( betweenKb_coexprDT_cond1$gene1 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT_cond1$gene1 %in% tad_genes )
      
      stopifnot( betweenKb_coexprDT_cond1$gene2 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT_cond1$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenKb_coexprDT_cond1$gene1 %in% tad_genes & 
                      betweenKb_coexprDT_cond1$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowKb %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowKb %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowKb %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowKb %in% within_coexprDT$gene2)
      
      
      
      
      ####################################################################
      #### COND 2
      ####################################################################

      tad_coexprDT_cond2 <- coexprDT_cond2[coexprDT_cond2$gene1 %in% tad_genes | 
                                             coexprDT_cond2$gene2 %in% tad_genes,]
      
      stopifnot(nrow(tad_coexprDT_cond2) > 0)  
      stopifnot(!is.na(tad_coexprDT_cond2$coexpr))
      stopifnot(is.numeric(tad_coexprDT_cond2$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT_cond2$gene1 | tad_genes %in% tad_coexprDT_cond2$gene2 )
      
      ### WITHIN
      within_coexprDT_cond2 <- tad_coexprDT_cond2[tad_coexprDT_cond2$gene1 %in% tad_genes & tad_coexprDT_cond2$gene2 %in% tad_genes,]
      
      ### BETWEEN: ALL
      
      betweenAll_coexprDT_cond2 <- tad_coexprDT_cond2[ ! (tad_coexprDT_cond2$gene1 %in% tad_genes & tad_coexprDT_cond2$gene2 %in% tad_genes),]
      
      stopifnot(nrow(within_coexprDT_cond2) + nrow(betweenAll_coexprDT_cond2) == nrow(tad_coexprDT_cond2))
      
      stopifnot(betweenAll_coexprDT_cond2$gene1 %in% tad_genes | 
                  betweenAll_coexprDT_cond2$gene2 %in% tad_genes )
      
      stopifnot( ! betweenAll_coexprDT_cond2$gene1 %in% tad_genes | 
                   ! betweenAll_coexprDT_cond2$gene2 %in% tad_genes )
      
      
      ### BETWEEN: GENES AROUND - NBR
      
      
      betweenNbr_coexprDT_cond2 <- tad_coexprDT_cond2[ tad_coexprDT_cond2$gene1 %in% betweenGenes_windowNbr | 
                                             tad_coexprDT_cond2$gene2 %in% betweenGenes_windowNbr,] 
      stopifnot( betweenNbr_coexprDT_cond2$gene1 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT_cond2$gene1 %in% tad_genes )
      
      stopifnot( betweenNbr_coexprDT_cond2$gene2 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT_cond2$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenNbr_coexprDT_cond2$gene1 %in% tad_genes & 
                      betweenNbr_coexprDT_cond2$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowNbr %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowNbr %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowNbr %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowNbr %in% within_coexprDT$gene2)
      
      
      
      
      ### BETWEEN: GENES AROUND - KB
      
      
      betweenKb_coexprDT_cond2 <- tad_coexprDT_cond2[ tad_coexprDT_cond2$gene1 %in% betweenGenes_windowKb | 
                                            tad_coexprDT_cond2$gene2 %in% betweenGenes_windowKb,] 
      stopifnot( betweenKb_coexprDT_cond2$gene1 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT_cond2$gene1 %in% tad_genes )
      
      stopifnot( betweenKb_coexprDT_cond2$gene2 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT_cond2$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenKb_coexprDT_cond2$gene1 %in% tad_genes & 
                      betweenKb_coexprDT_cond2$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowKb %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowKb %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowKb %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowKb %in% within_coexprDT$gene2)
      
      
      
      
      
      
      #####################################################################################################
      
      list(withinCoexpr = mean(within_coexprDT$coexpr),
           withinCoexpr_cond1 = mean(within_coexprDT_cond1$coexpr),
           withinCoexpr_cond2 = mean(within_coexprDT_cond2$coexpr),
           betweenAllCoexpr = mean(betweenAll_coexprDT$coexpr),
           betweenAllCoexpr_cond1 = mean(betweenAll_coexprDT_cond1$coexpr),
           betweenAllCoexpr_cond2 = mean(betweenAll_coexprDT_cond2$coexpr),
           betweenKbCoexpr = mean(betweenKb_coexprDT$coexpr),
           betweenKbCoexpr_cond1 = mean(betweenKb_coexprDT_cond1$coexpr),
           betweenKbCoexpr_cond2 = mean(betweenKb_coexprDT_cond2$coexpr),
           betweenNbrCoexpr = mean(betweenNbr_coexprDT$coexpr),
           betweenNbrCoexpr_cond1 = mean(betweenNbr_coexprDT_cond1$coexpr),
           betweenNbrCoexpr_cond2 = mean(betweenNbr_coexprDT_cond2$coexpr)
           )
    } # end-foreach iterating over TADs
    names(within_between_coexpr_data) <- pipeline_regionList
    within_between_coexpr_data
  } # end-foreach iterating over dataset
  names(allData_within_between_coexpr) <- coexprFilesOutNames
  
  outFile <- file.path(outFolder, "allData_within_between_coexpr.Rdata")
  save(allData_within_between_coexpr, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else { # end if-buildData
  outFile <- file.path(outFolder, "allData_within_between_coexpr.Rdata")
  allData_within_between_coexpr <- eval(parse(text = load(outFile)))
}


stop("--ok\n")
load("COEXPR_BETWEEN_WITHIN_SUB/allData_within_between_coexpr.Rdata")

# sort the TADs by decreasing withinCoexpr
# plot level of coexpr within and between on the same plot

tad_coexpr_DT <- data.frame(
  withinCoexpr = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                          function(sublist) lapply(sublist, function(x) x[["withinCoexpr"]])))),
  betweenAllCoexpr = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                              function(sublist) lapply(sublist, function(x) x[["betweenAllCoexpr"]])))),
  betweenKbCoexpr = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                             function(sublist) lapply(sublist, function(x) x[["betweenKbCoexpr"]])))),
  betweenNbrCoexpr = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                              function(sublist) lapply(sublist, function(x) x[["betweenNbrCoexpr"]])))),

  withinCoexpr_cond1 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                function(sublist) lapply(sublist, function(x) x[["withinCoexpr_cond1"]])))),
  betweenAllCoexpr_cond1 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                    function(sublist) lapply(sublist, function(x) x[["betweenAllCoexpr_cond1"]])))),
  betweenKbCoexpr_cond1 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                   function(sublist) lapply(sublist, function(x) x[["betweenKbCoexpr_cond1"]])))),
  betweenNbrCoexpr_cond1 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                    function(sublist) lapply(sublist, function(x) x[["betweenNbrCoexpr_cond1"]])))),

  
  withinCoexpr_cond2 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                          function(sublist) lapply(sublist, function(x) x[["withinCoexpr_cond2"]])))),
  betweenAllCoexpr_cond2 = as.numeric(unlist(lapply(allData_within_between_coexpr,
                                              function(sublist) lapply(sublist, function(x) x[["betweenAllCoexpr_cond2"]])))),
  betweenKbCoexpr_cond2 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                             function(sublist) lapply(sublist, function(x) x[["betweenKbCoexpr_cond2"]])))),
  betweenNbrCoexpr_cond2 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                              function(sublist) lapply(sublist, function(x) x[["betweenNbrCoexpr_cond2"]])))),
  
  
  stringsAsFactors = FALSE
)
tad_coexpr_DT <- tad_coexpr_DT[order(tad_coexpr_DT$withinCoexpr, decreasing = TRUE),]
tad_coexpr_DT$TADrank <- 1:nrow(tad_coexpr_DT)


suffix=""
betwType="All"
for(suffix in c("", "_cond1", "_cond2")) {
  
  var1 <- paste0("withinCoexpr", suffix)
  stopifnot(var1 %in% colnames(tad_coexpr_DT))
  
  
  for(betwType in c("All", "Nbr", "Kb")) {
    
    var2 <- paste0("between", betwType, "Coexpr", suffix)
  
    stopifnot(var2 %in% colnames(tad_coexpr_DT))
      
    myy1 <- tad_coexpr_DT[, var1] 
    myy2 <- tad_coexpr_DT[, var2] 
    
    outFile <- file.path(outFolder, paste0("within_vs_between", betwType, suffix, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(y=myy1,
         x=tad_coexpr_DT$TADrank,
         ylim=c(range(myy1, myy2, na.rm = TRUE)),
         pch=16,
         cex=0.7,
         cex.axis = myCexAxis,
         cex.lab = myCexLab,
         col="black",
         xlab="TAD",
         ylab=paste0("mean pairwise ", toTitleCase(corrMet), " corr."),
         main = paste0("Between vs. within expression correlation")
    )
    points(y=myy2,
           x=tad_coexpr_DT$TADrank,
           pch=16,
           cex=0.7,
           col="red"
           )
    
    
    if(suffix != "") {
      mtext(side=3, text = paste0("Between: ", betwType), font = 3)  
    } else {
      mtext(side=3, text = paste0("Between: ", betwType, " - ", gusb("_", suffix)), font = 3)
    }
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  
}

plotvar="withinCoexpr"

for(plotvar in c("withinCoexpr", "betweenAllCoexpr", "betweenKbCoexpr", "betweenNbrCoexpr")) {
  
  var1 <- paste0(plotvar, "_cond1")
  stopifnot(var1 %in% colnames(tad_coexpr_DT))
  
  var2 <- paste0(plotvar, "_cond2")
  stopifnot(var2 %in% colnames(tad_coexpr_DT))
  
  myy1 <- tad_coexpr_DT[, var1] 
  myy2 <- tad_coexpr_DT[, var2] 
  
  outFile <- file.path(outFolder, paste0("cond1_vs_cond2_", plotvar, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(y=myy1,
       x=myy2,
       ylim=c(range(myy1, myy2, na.rm = TRUE)),
       pch=16,
       cex=0.7,
       cex.axis = myCexAxis,
       cex.lab = myCexLab,
       col="black",
       xlab="TAD",
       ylab=paste0("mean pairwise ", toTitleCase(corrMet), " corr."),
       main = paste0("Expr. corr.: ", plotvar, " - cond1 vs. cond2")
  )
  addCorr(x=myy1,y=myy2,legPos="topleft", bty='n')
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  

}


  


# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





