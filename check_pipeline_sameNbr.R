startTime <- Sys.time()

check5 <- TRUE
check7 <- TRUE
check10 <- TRUE
check11 <- TRUE
check19 <- TRUE

refFolder <- file.path("..", "2_Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(refFolder))

pipFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE","OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

script5_name <- "5sameNbr_runPermutationsCorr"
script7_name <- "7sameNbr_runPermutationsMeanTADCorr"
script10_name <- "10sameNbr_runEmpPvalMeanTADCorr"
script11_name <- "11sameNbr_runEmpPvalCombined"
script19_name <- "19sameNbr_SAM_emp_measurement"

ref_script5_folder <- file.path(refFolder, "CREATE_SAMPLE_AROUND_TADS_SAMENBR_TWOSIDED")
stopifnot(dir.exists(ref_script5_folder))

ref_script7_folder <- file.path(refFolder, "COEXPR_AROUND_TADS")
stopifnot(dir.exists(ref_script7_folder))

ref_script7_file <- file.path(ref_script7_folder, "sameNbr", "all_ds_around_TADs_corr.Rdata")
stopifnot(file.exists(ref_script7_file))

ref_script10_folder <- file.path(refFolder, "SAMPLE_MEANCORR_EMPPVALS")
stopifnot(dir.exists(ref_script10_folder))

ref_script11_folder <- file.path(refFolder, "CREATE_EMPPVALCOMB", "sameNbr")
stopifnot(dir.exists(ref_script11_folder))

ref_script19_folder <- file.path(refFolder, "SAM_EMP_MEASUREMENT_MEANCORR", "sameNbr")
stopifnot(dir.exists(ref_script19_folder))

hicds="GSE105381_HepG2_40kb"
exprds="TCGAlihc_norm_lihc"
########################################################################
### CHECK SCRIPT 5 [sampling]
########################################################################
if(check5){
  for(hicds in all_hicds){
    for(exprds in all_exprds[[paste0(hicds)]]) {
      cat("... CHECK SCRIPT 5: ", hicds, " - ", exprds, "\n")
      # script5_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/5sameNbr_runPermutationsCorr/sample_around_TADs_sameNbr.Rdata"
      script5_file <- file.path(pipFolder, hicds, exprds, script5_name, "sample_around_TADs_sameNbr.Rdata")
      stopifnot(file.exists(script5_file))
      pipData <- eval(parse(text = load(script5_file)))
      # ref_file <- "/media/electron/mnt/etemp/marie/2_Yuanlong_Cancer_HiC_data_TAD_DA/CREATE_SAMPLE_AROUND_TADS_SAMENBR_TWOSIDED/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/sample_around_TADs_sameNbr.Rdata"
      ref_file <- file.path(ref_script5_folder, hicds, exprds, "sample_around_TADs_sameNbr.Rdata")
      stopifnot(file.exists(ref_file))  
      refData <- eval(parse(text = load(ref_file)))
      stopifnot(all.equal(refData, pipData))
    }
  }
}
########################################################################
### CHECK SCRIPT 7 [sampling correlation]
########################################################################
if(check7){
  # ref_script7_file <- "../2_Yuanlong_Cancer_HiC_data_TAD_DA/COEXPR_AROUND_TADS/sameNbr/all_ds_around_TADs_corr.Rdata"
  refDataAll <- eval(parse(text = load(ref_script7_file)))
  for(hicds in all_hicds){
    for(exprds in all_exprds[[paste0(hicds)]]) {
      cat("... CHECK SCRIPT 7: ", hicds, " - ", exprds, "\n")
      # script7_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/7sameNbr_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata"
      script7_file <- file.path(pipFolder, hicds, exprds, script7_name, "meanCorr_sample_around_TADs_sameNbr.Rdata")
      stopifnot(file.exists(script7_file))
      pipData <- eval(parse(text = load(script7_file)))   
      refData <- refDataAll[[file.path(hicds, exprds)]]
      stopifnot(all.equal(refData, pipData))
    }
  }
}
########################################################################
### CHECK SCRIPT 10 [emp. pval. corr]
########################################################################
if(check10){
  for(hicds in all_hicds){
    for(exprds in all_exprds[[paste0(hicds)]]) {
      cat("... CHECK SCRIPT 10: ", hicds, " - ", exprds, "\n")
      # script10_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/10sameNbr_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata"
      script10_file <- file.path(pipFolder, hicds, exprds, script10_name, "emp_pval_meanCorr.Rdata")
      stopifnot(file.exists(script10_file))
      pipData <- eval(parse(text = load(script10_file)))   
      # SAMPLE_MEANCORR_EMPPVALS/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/all_cors_empPval_dt.Rdata
      ref_file <- file.path(ref_script10_folder, hicds, exprds, "all_cors_empPval_dt.Rdata")
      refDT <- eval(parse(text = load(ref_file)))
      refData <- setNames(refDT[, "empPval-sameNbr-meanCorr - allDS"], rownames(refDT))
      stopifnot(all.equal(refData, pipData))
    }
  }
}
########################################################################
### CHECK SCRIPT 11 [emp. pval. comb]
########################################################################
if(check11){
  for(hicds in all_hicds){
    for(exprds in all_exprds[[paste0(hicds)]]) {
      cat("... CHECK SCRIPT 11: ", hicds, " - ", exprds, "\n")
      # script11_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/11sameNbr_runEmpPvalCombined/emp_pval_combined.Rdata"
      script11_file <- file.path(pipFolder, hicds, exprds, script11_name, "emp_pval_combined.Rdata")
      stopifnot(file.exists(script11_file))
      pipData <- eval(parse(text = load(script11_file)))   
      # CREATE_EMPPVALCOMB/sameNbr/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/meanCorr_meanLogFC_notAdjCombEmpPval.Rdata
      ref_file <- file.path(ref_script11_folder, hicds, exprds, "meanCorr_meanLogFC_notAdjCombEmpPval.Rdata")
      refData <- eval(parse(text = load(ref_file)))
      stopifnot(all.equal(refData, pipData))
    }
  }
}
########################################################################
### CHECK SCRIPT 19 [SAM FDR]
########################################################################
if(check19){
  # /mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/K562_40kb/TCGAlaml_wt_mutFLT3/
  for(hicds in all_hicds){
    for(exprds in all_exprds[[paste0(hicds)]]) {
      cat("... CHECK SCRIPT 19: ", hicds, " - ", exprds, "\n")
      # script19_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/19sameNbr_SAM_emp_measurement/meanCorr_empFDR.Rdata"
      script19_file <- file.path(pipFolder, hicds, exprds, script19_name, "meanCorr_empFDR.Rdata")
      stopifnot(file.exists(script19_file))
      pipData <- eval(parse(text = load(script19_file)))   
      pipData1 <- pipData[["nbrSignif"]]
      pipData2 <- pipData[["empFDR"]]
      # SAM_EMP_MEASUREMENT_MEANCORR/sameNbr/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/all_empFDR.Rdata
      ref_file <- file.path(ref_script19_folder, hicds, exprds, "all_empFDR.Rdata")
      refData <- eval(parse(text = load(ref_file)))
      refDataMeanCorr <- refData[["sample_meanCorr_allDS"]]
      refData1 <- refDataMeanCorr[["nbrSignif"]]
      refData2 <- refDataMeanCorr[["empFDR"]]
      stopifnot(all.equal(refData1, pipData1))
      stopifnot(all.equal(refData2, pipData2))
    }
  }
}
########################################################################
########################################################################
########################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
