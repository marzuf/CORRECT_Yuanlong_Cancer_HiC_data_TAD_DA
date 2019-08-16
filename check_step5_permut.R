

g_v2_file <- "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/geneAggregExpression.Rdata"
g_v2 <- get(load(g_v2_file))

g1_v0_file <- "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/geneAggregExpression2.Rdata"
g1_v0 <- get(load(g1_v0_file))
g1_v0$shuffRegion <- NULL

all.equal(g1_v0, g_v2)

g2_v0_file <- "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/geneAggregExpression3.Rdata"
g2_v0 <- get(load(g2_v0_file))
g2_v0$shuffRegion <- NULL
all.equal(g2_v0, g_v2)

g3_v0_file <- "/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom/geneAggregExpression4.Rdata"
g3_v0 <- get(load(g3_v0_file))
g3_v0$shuffRegion <- NULL
all.equal(g3_v0, g_v2)

region_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/0_prepGeneData/pipeline_regionList.Rdata" 
regionList <- get(load(region_file))

v2_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/5testv2_runPermutationsMedian/permutationsDT.Rdata"
v2_permut <- get(load(v2_file))

v0_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/5testv0_runPermutationsMedian/permutationsDT.Rdata"
v0_permut <- get(load(v0_file))

head(v0_permut)
head(v2_permut)

all(v0_permut[,1] == v2_permut[,1])

apply(v0_permut, 2, function(x) setequal(x, regionList))
table0 <- apply(v0_permut, 2, table)
table0 <- table0[regionList,]
apply(v2_permut, 2, function(x) setequal(x, regionList))
table2 <- apply(v0_permut, 2, table)
table2 <- table2[regionList,]

all.equal(table0, table2)


x = v2_permut[,2]

apply(v0_permut, 2, function(i) all(i == x))


source("samp1.R")
myfunc(1)
myfunc(5)
source("samp2.R")
myfunc(1)
myfunc(5)


# => might not be comparable because of multiple cores !!
# even when running the same script twice with same seed does not yield same results !!!

perm_file <- "/media/electron/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/510000_runPermutationsMedian/permutationsDT.Rdata"
permDT <- get(load(perm_file))

uniqTadAssign <- apply(permDT,1, unique)