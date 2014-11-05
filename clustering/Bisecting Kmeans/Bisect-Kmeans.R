library(utils)          # for read/write csv files
library(flexclust)      # for kcca

# This function conputers the distance of each point to the centroid and assigns
# a cluster to each point w.r.t. to the closest centroid
# It also computes the cluster and clustering SSE
computeCluster <- function(graphDat,centroid,ncl){
  
  tot.withiness <- 0
  cluster <- rep(NA,nrow(graphDat))
  clusterSSE <- rep(0,ncl)
  for(i in 1:nrow(graphDat)){ # Loop runs for all points
    min = Inf
    for(j in 1:ncl){ # Loop runs for each centroids
      dist = sqrt((graphDat$length[i] - centroid[j,1])^2+(graphDat$width[i] - centroid[j,2])^2)
      if(dist < min){
        min = dist
        cluster[i] <- j  # assigns a cluster to each point
        clusterSSE[j] <- clusterSSE[j] + min  # Cluster SSE
        tot.withiness <- tot.withiness + min  # Clustering SSE
      }
    }
  }
  l = list(cluster = cluster, tot.withiness = tot.withiness, withiness = clusterSSE)
  l
}

# This function mimics the K-means algorithm 
mykmeans <- function(graphData,k, iter = 10){
  
  # Sampling initial centroids
  X <- sample(graphData$length, k)
  Y <- sample(graphData$width, k)
  centroids <- cbind(X,Y)
  
  # Computes initial clusters
  lst <- computeCluster(graphData,centroids,k)
  
  # Iterate and compute new centroids untill the no. of iterations
  while(iter > 0) {
    for(i in 1:k){
      centroids[i,1] <- mean(graphData$length[which(lst$cluster == i)])
      centroids[i,2] <- mean(graphData$width[which(lst$cluster == i)])
    }
    lst <- computeCluster(graphData,centroids,k)
    iter <- iter - 1
  }
  # List containing the clusters centroids, cluster sse and clustering sse
  l <- list(clusters = lst$cluster,centroids = centroids,tot.withiness = lst$tot.withiness, withiness = lst$withiness)
  l
}

bisectKmeans <- function(clusterData,k,ntrials = 3){
  
  # Variables initialization
  noOfClus <- 1
  subCluster <- clusterData
  splitClusterIndex <- 0
  finalkc <- list()
  clustIndex <- rep(0,k)
  
  # Performs bisecting kmeans until k clusters are obtained
  while(noOfClus != k) {
    minSSE <- Inf
    
    # Computes Kmeans for a selected cluster for n trials and
    # Selects the pair of clusters which has the minimum total 
    #clustering sse
    for(i in 1:ntrials){
      kcl <- mykmeans(subCluster,2)
      if(kcl$tot.withiness < minSSE){ 
        minSSE = kcl$tot.withiness
        splitkc <- kcl
      }
    }
    noOfClus <- noOfClus + 1
    
    # If there are only 2 clusters return the pair of clusters
    # obtained in the first split
    if(noOfClus == 2){
      finalkc <- splitkc
      clustIndex <- c(1,2)
    }# else maintain the final cluster informartion by storing all the 
    # cluster split information
    else{ 
      # Storing and Updating cluster SSE after the split
      finalkc$withiness[splitClusterIndex] = splitkc$withiness[1]
      finalkc$withiness[noOfClus] = splitkc$withiness[2]
      
      # Storing and Updating centroids after the split
      finalkc$centroids[splitClusterIndex,] = splitkc$centroids[1,]
      finalkc$centroids = rbind(finalkc$centroids,splitkc$centroids[2,])
      
      # Storing and Updating clustering SSE after the split
      finalkc$tot.withiness = sum(finalkc$withiness)
      
      # Storing and Updating cluster INDEX after the split
      clustIndex[splitClusterIndex] = splitClusterIndex
      clustIndex[noOfClus] = noOfClus
      
      # Storing and Updating cluster labels after the split
      splitkc$clusters[splitkc$clusters == 2] = noOfClus
      splitkc$clusters[splitkc$clusters == 1] = splitClusterIndex
      finalkc$clusters[finalkc$clusters == splitClusterIndex] = splitkc$clusters
    }
    
    # Select the cluster which has the highest Cluster SSE for splitting
    if(noOfClus != k){
      splitClusterIndex <- clustIndex[which.max(finalkc$withiness)]
      subCluster <- clusterData[which(finalkc$clusters == splitClusterIndex),]
    }  
  }
  finalkc
}

# Implementation

# Read data from the file
myData <- read.csv('d-c4hw2.csv')

# Subset appropriate columns
clusterData <- myData[,1:2]

# Plot the raw data as it is before clustering
plot(clusterData,xlab = "length", ylab = "width",main = "Raw Data")


# Applying K-means with K = 4 and plotting data after kmeans clustering
kc <- mykmeans(clusterData,4)
plot(clusterData,col = kc$clusters,xlab = "length", ylab = "width", main = "K-means Clustering (k = 4) Output Plot")

# Applying K-means with K = 8 and plotting data after kmeans clustering
kc <- mykmeans(clusterData,8)
plot(clusterData,col = kc$clusters,xlab = "length", ylab = "width", main = "K-means Clustering (k = 8) Output Plot")


# Applying bisecting k-means for k = 2 to 8 and applying 
# Elbow Method for selecting appropriate k
wss <- (nrow(clusterData)-1)*sum(apply(clusterData,2,var))
for (i in 2:8){
  wss[i] <- sum(bisectKmeans(clusterData,i,3)$withiness)
} 
# plot the elbow curve to select optimal k
plot(1:8, wss, type="b", xlab="Number of Clusters", ylab="SSE",main = "Elbow Plot")


# Applying BiSecting K-means for k= , which is the optimal number
# of clusters as derived from the elbow method
# and plotting data after Bisecting kmeans clustering
bkc <- bisectKmeans(clusterData, 4)
plot(clusterData,col = bkc$clusters,xlab = "length", ylab = "width", main = "Bisecting K-means Clustering (k = 4) Output Plot")
