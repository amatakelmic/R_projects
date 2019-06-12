#updating to latest version of Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()




#installing mzR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("mzR")




#installing faahko
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("faahKO")


#selecting package repositories
setRepositories(addURLs =
                  c(CRANxtras = "http://www.stats.ox.ac.uk/pub/RWin"))

#installing packages
install.packages("RColorBrewer")
install.packages("pander", repos="http://cran.us.r-project.org")
install.packages("magrittr")
install.packages("pheatmap", repos="http://cran.us.r-project.org")

#calling up libraries
library(xcms)
library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(pheatmap)

## Get the full path to the CDF files
cdfs <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE,
            recursive = TRUE)
## Create a phenodata data.frame
pd <- data.frame(sample_name = sub(basename(cdfs), pattern = ".CDF",
                                   replacement = "", fixed = TRUE),
                 sample_group = c(rep("KO", 6), rep("WT", 6)),
                 stringsAsFactors = FALSE)

#loading  the raw data as OnDiskMSnExp

raw_data <- readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd),
                       mode = "onDisk")

#data inspection
head(rtime(raw_data))
mzs <- mz(raw_data)

## Split the list by file
mzs_by_file <- split(mzs, f = fromFile(raw_data))
length(mzs_by_file)


## Get the base peak chromatograms. This reads data from the files.
bpis <- chromatogram(raw_data, aggregationFun = "max")
## Define colors for the two groups
group_colors <- paste0(brewer.pal(3, "Set1")[1:2], "60")
names(group_colors) <- c("KO", "WT")

## Plot all chromatograms.
plot(bpis, col = group_colors[raw_data$sample_group])



#extract of the first sample (RT and Intensity)
bpi_1 <- bpis[1, 1]
head(rtime(bpi_1))

## F01.S0001 F01.S0002 F01.S0003 F01.S0004 F01.S0005 F01.S0006 
##  2501.378  2502.943  2504.508  2506.073  2507.638  2509.203

head(intensity(bpi_1))



## Get the total ion current by file
tc <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tc, col = group_colors[raw_data$sample_group],
        ylab = "intensity", main = "Total ion current")



## Bin the BPC
bpis_bin <- bin(bpis, binSize = 2)

## Calculate correlation on the log2 transformed base peak intensities
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data$sample_name

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = raw_data$sample_group)
rownames(ann) <- raw_data$sample_name

## Perform the cluster analysis
pheatmap(cormat, annotation = ann,
         annotation_color = list(group = group_colors))





















