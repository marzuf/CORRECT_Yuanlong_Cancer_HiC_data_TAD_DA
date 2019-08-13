
# Rscript extract_Yuanlong_FINAL_DOMAINS_corrected.R

# => CORRECTED VERSION:
# 1) TADs on left and right of a gap might have same TAD rank
# so now impose a threshold:
# if left and right are separated by more than "gap_threshold" bins => they receive a different CD_rank
# 2) discard TADs if they contain the mid position of the centromere

hicds = "ENCSR346DCU_LNCaP_40kb"
chromo = "chr2"

script_name <- "extract_Yuanlong_FINAL_DOMAINS_corrected.R"
startTime <- Sys.time()

options(scipen=100)

SSHFS <- FALSE
prefixDir <- ifelse(SSHFS, "/media/electron", "")
allData_file <- file.path(prefixDir, "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/bin_comp_CELL_LINEs_list.Rdata")
stopifnot(file.exists(allData_file))

infoFile <- file.path(prefixDir, "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/CELL_LINE_info.R")
source(infoFile)
cl_long_short_names <- setNames(CELL_LINEs,short_names)
can_short_names <- names(is_CAN[is_CAN==1])
toKeepNames <- cl_long_short_names[can_short_names]