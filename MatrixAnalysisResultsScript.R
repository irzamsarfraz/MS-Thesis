assay(sce, "counts") <- as.matrix(assay(sce, "counts"))

#################################
#     Use ExperimentSubset      #
#################################

#Load PBMC4K dataset and create ExperimentSubset object:
# sce <- TENxPBMCData(dataset = "pbmc4k")
es <- ExperimentSubset(sce)

#Compute perCellQCMetrics on counts matrix:
perCellQCMetrics <- perCellQCMetrics(assay(es, "counts"))
colData(es) <- cbind(colData(es), perCellQCMetrics)

#Filter cells with low column sum and create a new subset called ‘filteredCells’:
## rank cells in order of highest sum
es <- es[, order(-es@colData$sum)]
#filteredCellsIndices <- which(colData(es)$sum > 1500)
## keep top cells only 40% in a new subset
es <- createSubset(es, "filteredCells", cols = 1:round(ncol(es) * 0.40), parentAssay = "counts")
#es <- createSubset(es, "filteredCells", cols = filteredCellsIndices, parentAssay = "counts")


#Normalize ‘filteredCells’ subset using scater library and store it back:
# assay(es, "filteredCells", subsetAssayName = "filteredCellsNormalized") <- normalizeCounts(assay(es, "filteredCells"))
es <- ExperimentSubset::setSubsetAssay(es, "filteredCells", normalizeCounts(assay(es, "filteredCells")), subsetAssayName = "filteredCellsNormalized")

#Find highly variable genes from the normalized assay in the previous step using scran library against the ‘filteredCells’ subset only:
topHVG1000 <- getTopHVGs(modelGeneVar(assay(es, "filteredCellsNormalized")), n = 1000)
es <- createSubset(es, "hvg1000", rows = topHVG1000, parentAssay = "filteredCellsNormalized")

#Run ‘PCA’ on the highly variable genes computed in the last step using scater library against the ‘filteredCells’ subset only:
reducedDim(es, type = "PCA", subsetName = "hvg1000") <- calculatePCA(assay(es, "hvg1000"))

#Find clusters from reducedDims of highly variable genes and store the results back into colData of subset
clusterPC = kmeans(reducedDim(es, "PCA", subsetName = "hvg1000"), 5)$cluster
colData(es, subsetName = "hvg1000") <- cbind(colData(es, subsetName = "hvg1000"), clusterPC)

#Show the current condition of the ExperimentSubset object:
subsetSummary(es)

# #Single ES object size
# print(object_size(es), unit = "Mb")


#####################################
# Use SingleCellExperiment (pbmc4k) #
#####################################

#Load PBMC4K dataset
# sce1 <- TENxPBMCData(dataset = "pbmc4k")
sce1 <- sce

#Compute perCellQCMetrics on counts matrix:
perCellQCMetrics <- perCellQCMetrics(assay(sce1, "counts"))
colData(sce1) <- cbind(colData(sce1), perCellQCMetrics)

#Filter cells with low column sum and create a new object:
filteredCellsIndices <- which(colData(sce1)$sum > 1500)
# sce1 <- sce1[, order(-sce1@colData$sum)]
sce2 <- sce1[,filteredCellsIndices]
# sce2 <- sce1[,1:round(ncol(es) * 0.40)]

#Normalize sce2 using scater library and store it back:
assay(sce2, "filteredCellsNormalized") <- normalizeCounts(assay(sce2, "counts"))

#Find highly variable genes from the normalized assay in the previous step using scran library and make a new object:
topHVG1000 <- getTopHVGs(modelGeneVar(assay(sce2, "filteredCellsNormalized")), n = 1000)
sce3 <- sce2[topHVG1000,]

#Run ‘PCA’ on the highly variable genes computed in the last step using scater library:
reducedDim(sce3, type = "PCA") <- calculatePCA(assay(sce3, "filteredCellsNormalized"))

#Find clusters from reducedDims of highly variable genes and store the results back into colData
clusterPC = kmeans(reducedDim(sce3, "PCA"), 5)$cluster
colData(sce3) <- cbind(colData(sce3), clusterPC)

#Three SCE object sizes
totalSize <- object_size(sce1) + object_size(sce2) + object_size(sce3) 


print("Dimensions")
print(dim(sce))

print("Original SCE Total Size:")
print(totalSize, unit = "Mb")

print("ES Total Size:")
#Single ES object size
print(object_size(es), unit = "Mb")

