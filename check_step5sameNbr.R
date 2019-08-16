mainFolder <- file.path("..", "Yuanlong_Cancer_HiC_data_TAD_DA")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

script5_name <- "5sameNbr_runPermutationsCorr"

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

all_data_list <- foreach(hicds = all_hicds) %dopar% {
  
exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {

  pipfile <- file.path(pipFolder, hicds, exprds, script5_name, "sample_around_TADs_sameNbr.Rdata")
  
  pipdata <- get(pipfile)
  } 
}

# not equal
# 5sameNbr and 7sameNbr -> compute permutations and correlation for each dataset separately
# it is only in step 10sameNbr that it loads the files produced at step 7sameNbr
