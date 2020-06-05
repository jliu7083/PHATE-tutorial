install.packages("phateR")
library(phateR)
data(tree.data)
plot(prcomp(tree.data$data)$x, col=tree.data$branches)

##########----------quick start----------------------------#######################
##  run PHATE on the data. We’ll just go ahead and try with the default parameters.
# runs phate
tree.phate <- phate(tree.data$data)
summary(tree.phate)
# plot embedding
palette(rainbow(10))
plot(tree.phate, col = tree.data$branches)

## runs phate with different parameters
tree.phate <- phate(tree.data$data, gamma=0, t=120, init=tree.phate)
# plot embedding
palette(rainbow(10))
plot(tree.phate, col = tree.data$branches)

## ploting in ggplot2
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 3.5.3
p<-ggplot(tree.phate, aes(x=PHATE1, y=PHATE2, color=tree.data$branches)) +
  geom_point()

ggsave('testing.png')
#########  plotly ############
library(plotly)
ggplotly(p)


#########--------------tutorial--------------------------######################
library(phateR)
library(ggplot2)
library(readr)
library(viridis)
install.packages("Rmagic")
library(Rmagic)
## load data
bmmsc <- read_csv("https://github.com/KrishnaswamyLab/PHATE/raw/master/data/BMMC_myeloid.csv.gz")
bmmsc[1:10,1:10]
bmmsc <- bmmsc[,2:ncol(bmmsc)]
bmmsc[1:5,1:10]
## filtering data
# keep genes expressed in at least 10 cells
keep_cols <- colSums(bmmsc > 0) > 10
bmmsc <- bmmsc[,keep_cols]
# look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(bmmsc)), bins=50) +
  geom_vline(xintercept = 1000, color='red')
# keep cells with at least 1000 UMIs
keep_rows <- rowSums(bmmsc) > 1000
bmmsc <- bmmsc[keep_rows,]

## Normalizing data
## We should library size normalize and transform the data prior to PHATE. 
## Many people use a log transform, which requires adding a “pseudocount” to avoid log(0). 
## We square root instead, which has a similar form but doesn’t suffer from instabilities at zero.
bmmsc <- library.size.normalize(bmmsc)
bmmsc <- sqrt(bmmsc)
bmmsc[1:5,1:10] # check
## Running PCA,  examine the raw data with PCA.
bmmsc_PCA <- as.data.frame(prcomp(bmmsc)$x)
# plot PCA, and color the plot by Mpo---->a myeloid marker. 
ggplot(bmmsc_PCA) +
  geom_point(aes(PC1, PC2, color=bmmsc$Mpo)) +
  labs(color="Mpo") +
  scale_color_viridis(option="B")


## Running PHATE with default
bmmsc_PHATE <- phate(bmmsc)
#plot
ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=bmmsc$Mpo)) +
  labs(color="Mpo") +
  scale_color_viridis(option="B")
ggsave('BMMSC_data_R_phate.png', width=5, height=5)

## Rerunning PHATE with new parameters
bmmsc_PHATE <- phate(bmmsc, knn=4, decay=100, t=10, init=bmmsc_PHATE)
ggplot(bmmsc_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=bmmsc$Mpo)) +
  labs(color="Mpo") +
  scale_color_viridis(option="C")
ggsave('BMMSC_data_R_phate.png', width=5, height=5)


##Flow data
counts<-read.table('./Desktop/counts.txt',sep="\t",header = TRUE)
rownames(counts) <-counts$Feature
counts$Feature <- NULL
counts[1:10,1:10]
counts_transpose <- as.data.frame(t(as.matrix(counts)))
counts_transpose[1:10,1:10]
## normalize
count <- library.size.normalize(counts_transpose)
count <- sqrt(count)
#count[1:10,1:10]
count_PHATE <- phate(counts_transpose,t=10)
#plot
ggplot(count_PHATE) +
  geom_point(aes(PHATE1, PHATE2, color=count$LTBR)) +
  labs(color="LTBR") +
  scale_color_viridis(option="B")
ggsave('BMMSC_data_R_phate.png', width=5, height=5)
