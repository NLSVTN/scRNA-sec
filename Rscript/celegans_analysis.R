# For libraries installation
source("https://bioconductor.org/biocLite.R")

name_of_file <- "GSM2599701_Gene.count.matrix.celegans.cell.Rdata"
# Uncomment when unzipping '.gz' file the first time 
# path_to_zipped_file <- paste0("/home/nikita/Desktop/CElegans_raw_data/",
#                              name_of_file,
#                              ".gz", sep="")
# gunzip(path_to_zipped_file)
path_to_file <- paste0("/home/nikita/Desktop/CElegans_raw_data/",
name_of_file, sep="")
load(path_to_file, verbose=TRUE)

# Removing genes with 0 - counts in all of the columns
num_removed = 0
for(i in 1:nrow(UMI_count)) {
  idx = i - num_removed
  if(sum(UMI_count[idx,]) < 1) {
    UMI_count <- UMI_count[-idx,]
    num_removed = num_removed + 1
    print(i)
  }
}

# Removing low-abundant genes
colSums(assay(sce), na.rm = FALSE,
        dims = 1)
UMI_count <- UMI_count[rowSums(UMI_count)]

# fread is in 'data.table' package, so first need to 'library(data.table)' before 
# calling 'fread()'
library(data.table)
gene_symbols_table <- as.data.frame(fread(file=paste0("/home/nikita/Desktop/",
                                                      "CElegans_raw_data/gene_symbols_matching.csv"), 
                                                      header=TRUE))
# genes_with_symbols <- gene_symbols_table[1:3, c(1, 6)]

# filtering by briefDescription column where it contains 'mitoch' word
mt_genes <- dplyr::filter(gene_symbols_table, grepl("mitoch", briefDescription))
rownames(mt_genes) <- mt_genes$input

# Getting 'is.mito' (finding rows in sce that 
# correspond to mitochrondrial genes)
is.mito <- grepl(paste(rownames(mt_genes), collapse="|"), rownames(UMI_count))

# if we want just TRUE/FALSE values of the rows that have/don't have 'mitoch'
# in 'briefDescription' column of 'gene_symbols_table'
# is.mito <- grepl("mitoch", gene_symbols_table[, c("briefDescription")])
ncols = dim(UMI_count)[2]

biocLite("SingleCellExperiment")
biocLite("scran")
library(SingleCellExperiment)
library(scran)
list_of_sce <- list()
# Looping though the UMI_count 'split_factor' columns at a time
# split_factor = 500
# for(i in seq(1,ncols, split_factor)) {
#   num_loop = floor(i / split_factor) + 1
#   idx = ncols
#   if (i + split_factor < ncols) {
#     idx = i + split_factor
#   } 
   sce <- SingleCellExperiment(list(counts=UMI_count))
#   # Normalization of the dataset containing 
#   # heterogenous cell data (different cell types)
   clusters <- quickCluster(sce)
#   sce <- computeSumFactors(sce, cluster=clusters)
#   list_of_sce[[num_loop]] <- sce
# }
# cbind(list_of_sce)

# Quality Control
biocLite("scater")
library(scater)
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
par(mfrow=c(1,3), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)

# Assessing the percentage of filtered out, bad-quality cells
sce <- sce[,!(libsize.drop | feature.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
           Remaining=ncol(sce))

# Let's run PCA
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata") + fontsize
