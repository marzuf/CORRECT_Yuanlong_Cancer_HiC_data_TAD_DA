
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

########################################################################
### CHECK SCRIPT 5 [sampling]
########################################################################

hicds="GSE105381_HepG2_40kb"
exprds="TCGAlihc_norm_lihc"

for(hicds in all_hicds){
  for(exprds in all_exprds[[paste0(hicds)]]) {
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

########################################################################
### CHECK SCRIPT 7 [sampling correlation]
########################################################################


# ref_script7_file <- "../2_Yuanlong_Cancer_HiC_data_TAD_DA/COEXPR_AROUND_TADS/sameNbr/all_ds_around_TADs_corr.Rdata"
refDataAll <- eval(parse(text = load(ref_script7_file)))

for(hicds in all_hicds){
  for(exprds in all_exprds[[paste0(hicds)]]) {
    # script7_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/7sameNbr_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata"
    script7_file <- file.path(pipFolder, hicds, exprds, script7_name, "meanCorr_sample_around_TADs_sameNbr.Rdata")
    stopifnot(file.exists(script7_file))
    pipData <- eval(parse(text = load(script7_file)))   
    
    refData <- refDataAll[[file.path(hicds, exprds)]]
    stopifnot(all.equal(refData, pipData))
    
    
  }
}



########################################################################
### CHECK SCRIPT 10 [emp. pval. corr]
########################################################################


