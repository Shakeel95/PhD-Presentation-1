#-----------------
# Load data
# EDA, plots, ect.
#-----------------

source("scripts/build_datasets.R")

# check if data has been collected, load / collect
if (file.exists("data/snp_covid.csv")){
        
        snp.raw.covid <- read.csv("data/snp_covid.csv", row.names = 1)
        
} else {
        
        # collect 
        date.from <- as.Date("2020/01/01")
        date.to <- date.from + 100
        snp.raw.covid <- SnP500.data(date.from,date.to)
        
        # check missing values
        if (anyNA(snp.raw.covid)){
                drop <- which(apply(snp.raw.covid,2,anyNA) == T)
                snp.raw.covid <- snp.raw.covid[,-drop]
        }
        
        # write out
        write.csv(snp.raw.covid,"data/snp_covid.csv")
}



# plot raw price
matplot(snp.raw.covid,type = "l", xaxt = "n",
        main = "SnP500 (raw prices)",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.raw.covid), rownames(snp.raw.covid))


# plot log prices
matplot(log(snp.raw.covid),type = "l", xaxt = "n",
        main = "SnP500 (log-prices)",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.raw.covid), rownames(snp.raw.covid))


# plot log-returns
snp.LR.covid <- apply(log(snp.raw.covid),2,diff)
matplot(snp.LR.covid, type = "l", xaxt = "n",
        main = "SnP500 (log-returns)",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.LR.covid), rownames(snp.LR.covid))



#-----------------------------------------
# Instantaneous correlation of log-returns
# 
# (model = joint multivaraite normality)
#-----------------------------------------


source("scripts/cluster_processes.R")


# make cor-sim matrix, cluster
sim.cor <- as.dist((1-cor(snp.LR.covid))/2)
hc.cor <- hclust(sim.cor)
plot(hc.cor, labels = FALSE)


# cut at a sensible height
abline(h = 0.3, col = "red", lty = 2)
clusters.cor <- ticker.cluster.info(hc.cor, h = 0.3, SnP500.info = SnP500.info)


# visually similar time series 
cluster.names <- names(subset(clusters.cor$cluster.info, clusters == 15)$Symbol)
matplot(subset(snp.raw.covid, select = cluster.names),type = "l", xaxt = "n",
        main = "Visually similar clusters",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.LR.covid), rownames(snp.LR.covid))


# random sample of time series 
set.seed(1)
matplot(snp.raw.covid[,sample(1:ncol(snp.raw.covid),5)], type = "l", xaxt = "n",
        main = "Random collection of 5 time series",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.LR.covid), rownames(snp.LR.covid))


# visually dissimilar
cluster.names <- names(subset(clusters.cor$cluster.info, clusters == 21)$Symbol)
matplot(subset(snp.raw.covid, select = cluster.names),type = "l", xaxt = "n",
        main = "Visually dissimilar clusters",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.LR.covid), rownames(snp.LR.covid))



#-----------------------------------------------------------------
# L2 distance between params in AR representation of GARCH process
#
# (model = weakly stationary time series)
#-----------------------------------------------------------------


source("ar_expansion_dist.R")


# make similarity matrix, plot dendrogram
sim.arL2 <- ar.coef.similarity(snp.LR.covid)
hc.arL2 <- hclust(100*sim.arL2)
plot(hc.arL2,  labels = FALSE)


# cut at a sensible height
abline(h = 10, col = "red", lty = 2)
clusters.arL2 <- ticker.cluster.info(hc.arL2, h = 10, SnP500.info = SnP500.info)


# find visually dissimilar time series
# plenty of nonesense results to choose from! e.g. 
# * 20 * 12 * 2
cluster.names <- names(subset(clusters.arL2$cluster.info, clusters == 12)$Symbol)
matplot(subset(snp.raw.covid, select = cluster.names),type = "l", xaxt = "n",
        main = "Large cluster of visually dissimilar time series",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.raw.covid), rownames(snp.raw.covid))

#---------------------------------
# ONE MORE MODEL BASED CLUSTER ...
#
# (model = ...)
#---------------------------------



#------------------------------------------------
# Clusters with Frechet distance + NOT skeletons
#
# (model = piecewise linear process)
#------------------------------------------------


source("frechet_dist_skeleton.R")


# get clusters 
sim.frechet <- skeleton.similarity(snp.raw.covid, set.contrast = "pcwsConstMean", parallel.comp = TRUE)
hc.frechet <- hclust(sim.frechet)
plot(hc.frechet, main = "Frechet distance clusters", labels = FALSE)


# cut at sensible height 
abline(h = 0.3, col = "red", lty = 2)
clusters.frechet <- ticker.cluster.info(hc.frechet, h = 0.3, SnP500.info = SnP500.info)


# inspect clusters 
clusters.frechet$cluster.counts


# plot some clusters
cluster.names <- names(subset(clusters.frechet$cluster.info, clusters == 16)$Symbol)
matplot(subset(snp.raw.covid, select = cluster.names),type = "l", xaxt = "n",
        main = "Large cluster of visually similar time series",
        xlab = "Date",
        ylab = "Price"
)
axis(1, at = 1:nrow(snp.raw.covid), rownames(snp.raw.covid))

cluster.names <- names(subset(clusters.frechet$cluster.info, clusters == 14)$Symbol)
matplot(subset(log(snp.raw.covid), select = cluster.names),type = "l", xaxt = "n",
        main = "Small cluster of visually simialr time series",
        xlab = "Date",
        ylab = "Log-price"
)
axis(1, at = 1:nrow(snp.raw.covid), rownames(snp.raw.covid))

