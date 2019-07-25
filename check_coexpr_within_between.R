dt1_file <- "COEXPR_BETWEEN_WITHIN_ALL/allData_within_between_coexpr.Rdata"
dt1 <- eval(parse(text = load(dt1_file)))


# sort the TADs by decreasing withinCoexpr
# plot level of coexpr within and between on the same plot

  dataset = as.character(unlist(lapply(1:length(dt1), function(i) {
    ds_name <- names(dt1)[i]
    ds_name <- gsub("^CREATE_COEXPR_SORTNODUP/", "", ds_name)
    ds_name <- gsub("/pearson/coexprDT.Rdata$", "", ds_name)
    rep(ds_name, length(dt1[[i]]))
  })))
  
  region = as.character(unlist(lapply(1:length(dt1), function(i) {
    names(dt1[[i]])
  })))
  
  withinCoexpr = as.numeric(unlist(lapply(dt1, 
                                          function(sublist) lapply(sublist, function(x) x[["withinCoexpr"]]))))
  
  length(dataset)
  length(region)
  length(withinCoexpr)
  
  
  dt1_1 <- dt1[7]
  
  dataset_1 = as.character(unlist(lapply(1:length(dt1_1), function(i) {
    ds_name <- names(dt1_1)[i]
    ds_name <- gsub("^CREATE_COEXPR_SORTNODUP/", "", ds_name)
    ds_name <- gsub("/pearson/coexprDT.Rdata$", "", ds_name)
    rep(ds_name, length(dt1_1[[i]]))
  })))
  
  region_1 = as.character(unlist(lapply(1:length(dt1_1), function(i) {
    names(dt1_1[[i]])
  })))
  
  withinCoexpr_1 = as.numeric(unlist(lapply(dt1_1, 
                                          function(sublist) lapply(sublist, function(x) x[["withinCoexpr"]]))))
  
  length(dataset_1)
  length(region_1)
  length(withinCoexpr_1)
  
  
  
  # CREATE_COEXPR_SORTNODUP/ENCSR444WCZ_A549_40kb/TCGAluad_mutKRAS_mutEGFR/pearson/coexprDT.Rdata`$chr1_TAD71 => NULL
  
  
  
  
  dt2_file <- "COEXPR_BETWEEN_WITHIN_ALL_2/allData_within_between_coexpr.Rdata"
  dt2 <- eval(parse(text = load(dt2_file)))
  
  
  # sort the TADs by decreasing withinCoexpr
  # plot level of coexpr within and between on the same plot
  
  dataset = as.character(unlist(lapply(1:length(dt2), function(i) {
    ds_name <- names(dt2)[i]
    ds_name <- gsub("^CREATE_COEXPR_SORTNODUP/", "", ds_name)
    ds_name <- gsub("/pearson/coexprDT.Rdata$", "", ds_name)
    rep(ds_name, length(dt2[[i]]))
  })))
  
  region = as.character(unlist(lapply(1:length(dt2), function(i) {
    names(dt2[[i]])
  })))
  
  withinCoexpr = as.numeric(unlist(lapply(dt2, 
                                          function(sublist) lapply(sublist, function(x) x[["withinCoexpr"]]))))
  
  length(dataset)
  length(region)
  length(withinCoexpr)
  
  
  dt2_1 <- dt2[7]
  
  dataset_1 = as.character(unlist(lapply(1:length(dt2_1), function(i) {
    ds_name <- names(dt2_1)[i]
    ds_name <- gsub("^CREATE_COEXPR_SORTNODUP/", "", ds_name)
    ds_name <- gsub("/pearson/coexprDT.Rdata$", "", ds_name)
    rep(ds_name, length(dt2_1[[i]]))
  })))
  
  region_1 = as.character(unlist(lapply(1:length(dt2_1), function(i) {
    names(dt2_1[[i]])
  })))
  
  withinCoexpr_1 = as.numeric(unlist(lapply(dt2_1, 
                                            function(sublist) lapply(sublist, function(x) x[["withinCoexpr"]]))))
  
  length(dataset_1)
  length(region_1)
  length(withinCoexpr_1)
  
  
  
  