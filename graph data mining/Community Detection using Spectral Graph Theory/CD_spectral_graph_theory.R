# Homework: Community Detection using Spectral Graph Theory
install.packages('KRLS')
install.packages('cccd')
install.packages('igraph')
require(KRLS) 
library('cccd')
library('igraph')
#Problem 1: Generating the data set.

#Genenrating 20 normal random (x,y) coordinate points with mean 2
x1=rnorm(n = 20,mean = 2,sd = 1)
y1=rnorm(n = 20,mean = 2,sd = 1)
#Genenrating 20 normal random (x,y) coordinate points with mean 4
x2=rnorm(n = 20,mean = 4,sd = 1)
y2=rnorm(n = 20,mean = 4,sd = 1)
#Genenrating 20 normal random (x,y) coordinate points with mean 6
x3=rnorm(n = 20,mean = 6,sd = 1)
y3=rnorm(n = 20,mean = 6,sd = 1)

# Set seed to use same random numbers again
set.seed(6)

# Creating a mixture of 3 Gaussians generated above 
X=c(x1,x2,x3) 
Y=c(y1,y2,y3)
att = c(rep(16,20),rep(17,20),rep(18,20))
attc = c(rep(1,20),rep(2,20),rep(3,20))
A = cbind(X,Y)

# (a) Plotting all the points in a single 2-dimensional space by using different shapes for each mixture.
plot(X,Y,pch=as.integer(att),col=as.integer(attc),title( "Mixture of Gaussians in 2-dimensional space"))
legend( "topleft" , legend = c("mean 2","mean 4","mean 6"), col = c(1,2,3), text.col = c(1,2,3), pch = c(16,17,18))

# (b) Plotting a histogram of all the points.
hist(X)
hist(Y)

# Problem 2: Generating the similarity graphs.

# (a) KNN: The K-nearest neighbor graph using the value of K=10.

G1 <- nng(A,k=10,mutual = TRUE)
plot(G1)

# (b) GK: The complete similarity graph using the Gaussian kernel with sigma=1 as similarity function.

GaussianMatrix <- gausskernel(X = A, sigma = 1)
G2 <- graph.adjacency(adjmatrix = GaussianMatrix, weighted = TRUE,mode = 'undirected')
G2$layout=A
plot(G2)


# Problem 3:Characterizing the graph spectra.

# Graph Laplacian Matrix for K-Nearest Neighbor graph
L1 <- graph.laplacian(G1)

#Graph Lapacian Matrix for Complete Similarity Graph with GaussianKernal as the similarity function
L2 <- graph.laplacian(G2)

# Normalized Graph Laplacian Matrix for K-Nearest Neighbor graph
L1_hat <- graph.laplacian(G1,normalized = TRUE)

# Normalized Graph Lapacian Matrix for Complete Similarity Graph with GaussianKernal as the similarity function
L2_hat <- graph.laplacian(G2,normalized = TRUE)

# Calculate and Plot graph spectra for K-Nearest Neighbor graph
Lambda1 <- eigen(L1)
plot(seq_along(Lambda1$values),Lambda1$values,xlab="Index",ylab="Eigen Values",title("KNN Graph Spectra"))

# Calculate and Plot graph spectra for Gaussian Kernel Similarity Graph
Lambda2 <- eigen(L2)
plot(seq_along(Lambda2$values),Lambda2$values,xlab="Index",ylab="Eigen Values",title("Gaussian Kernel Graph Spectra"))

# Calculate and Plot graph spectra for Normalized K-Nearest Neighbor graph
Lambda1_hat <- eigen(L1_hat)
plot(seq_along(Lambda1_hat$values),Lambda1_hat$values,xlab="Index",ylab="Eigen Values",title("Normalized KNN Graph Spectra"))

# Calculate and Plot graph spectra for Normalized Gaussian Kernel Similarity Graph
Lambda2_hat <- eigen(L2_hat)
plot(seq_along(Lambda2_hat$values),Lambda2_hat$values,xlab="Index",ylab="Eigen Values",title("Normalized Gaussian Kernel Graph Spectra"))

# (c) Plot each graph's eigenvector plot for the eigenvector u corresponding to the second smallest eigenvalue
plot(seq_along(Lambda1$vectors[,59]),Lambda1$vectors[,59],pch=as.integer(att),col=as.integer(attc),xlab="Index",ylab="Eigen Vectors",title("KNN Graph Vector with second smallest Eigen Value"))
legend( "topleft" , legend = c("mean 2","mean 4","mean 6"), col = c(1,2,3), text.col = c(1,2,3), pch = c(16,17,18))
plot(seq_along(Lambda2$vectors[,59]),Lambda2$vectors[,59],pch=as.integer(att),col=as.integer(attc),xlab="Index",ylab="Eigen Vectors",title("Gaussian Kernel Graph Vector with second smallest Eigen Value"))
legend( "topleft" , legend = c("mean 2","mean 4","mean 6"), col = c(1,2,3), text.col = c(1,2,3), pch = c(16,17,18))
plot(seq_along(Lambda1$vectors[,59]),Lambda1$vectors[,59],pch=as.integer(att),col=as.integer(attc),xlab="Index",ylab="Eigen Vectors",title("Normalized KNN Graph Vector with second smallest Eigen Value"))
legend( "topleft" , legend = c("mean 2","mean 4","mean 6"), col = c(1,2,3), text.col = c(1,2,3), pch = c(16,17,18))
plot(seq_along(Lambda2$vectors[,59]),Lambda2$vectors[,59],pch=as.integer(att),col=as.integer(attc),xlab="Index",ylab="Eigen Vectors",title("Normalized Gaussian Kernel Graph Vector with second smallest Eigen Value"))
legend( "topleft" , legend = c("mean 2","mean 4","mean 6"), col = c(1,2,3), text.col = c(1,2,3), pch = c(16,17,18))


# (D) Plot for 2-way graph partitioning into S and V-S, the points from different mixtures will end up in different partitions
attc_new = rep(0,60)

partition <- function(dataframe,threshold,titl){

  for (i in 1:59){
    if(dataframe$vectors[i,59] < threshold)
    {
      attc_new[i] = 1
    }
    else{
      attc_new[i] = 2
    }
  }

  plot(seq_along(dataframe$vectors[,59]),dataframe$vectors[,59],pch=as.integer(att),col=as.integer(attc_new),xlab="Index",ylab="Eigen Vectors",title(titl))
  legend( "topleft" , legend = c("Red - (S)","Black - (V - S)"), col = c(2,1), text.col = c(2,1))
}

knnT = 0.1
gkT = 0.05
knnNT = 0.09
gkNT = 0.04

# Performing 2-way partitioning on normalized and unnormalized KNN Graph 

partition(Lambda1,knnT,"KNN Graph Vector with Threshold 0.1")
partition(Lambda2,gkT,"Gaussian Kernel Graph Vector with Threshold 0.05")
partition(Lambda1_hat,knnNT,"Normalized KNN Graph Vector with Threshold 0.09")
partition(Lambda2_hat,gkNT,"Normalized Gaussian Kernel Graph Vector with Threshold 0.04")

# (E),(F) Calculate the conductance (write the script) for each of the identified partitions, S and V-S for the KNN graph using both the normalized and unnormalized Laplacian.
part = rep("VS",60)
Normal_part = rep("VS",60)

for (k in 1:59){
  if(Lambda1$vectors[k,59] > knnT)
  {
    part[k] = "S"
  }
}

for (k in 1:59){
  if(Lambda1_hat$vectors[k,59] > knnNT)
  {
    Normal_part[k] = "S"
  }
}

AdjKNN <- get.adjacency(G1)

calculateConductance <- function(vector,partType){

  Cs = 0.0
  Ms = 0.0
  
  for(i in 1:nrow(AdjKNN)){
    for(j in 1:ncol(AdjKNN)){
      if(AdjKNN[i,j] == 1){
        if(vector[j] == partType && vector[i] == partType)
          Ms <- Ms + 1
        else if(vector[j] != partType && vector[i] == partType)
          Cs <- Cs + 1
      }
    }
  }
  print(Cs)
  print(Ms)
  return(Cs/(Ms + Cs))
}

Lower = (Lambda1$values[59]/2)
Lower

ConductanceS = calculateConductance(part,"S")
ConductanceS

ConductanceVS = calculateConductance(part,"VS")
ConductanceVS

Upper = sqrt(2*Lambda1$values[59])
Upper

Lower_Normal = (Lambda1_hat$values[59]/2)
Lower_Normal

ConductanceNormalS = calculateConductance(Normal_part,"S")
ConductanceNormalS

ConductanceNormalVS = calculateConductance(Normal_part,"VS")
ConductanceNormalVS

Upper_Normal = sqrt(2*Lambda1_hat$values[59])
Upper_Normal

# Problem 4 : Spectral graph clustering

# (A) k-means clustering algorithm provided by R on the data set in Exercise 1, using the Euclidean distance as the dissimilarity metric, and the value of k=3.
library(cluster)
library(fpc)
kc <- kmeans(A, 3)
table(kc$cluster)
plot(A, pch=kc$cluster,col=kc$cluster)
title( "Mixture of Gaussians K-Means Clustering")
legend( "topleft" , legend = c("Cluster 1","Cluster 2","Cluster 3"), col = c(1,2,3), text.col = c(1,2,3), pch = c(1,2,3))

# (B) Spectral graph clustering and plot the corresponding points in Ex.1 with the shapes based on the identified cluster 
UKnn <- cbind(Lambda1_hat$vectors[,60],Lambda1_hat$vectors[,59],Lambda1_hat$vectors[,58])
UGk <- cbind(Lambda2_hat$vectors[,60],Lambda2_hat$vectors[,59],Lambda2_hat$vectors[,58])

# Spectral graph clustering for Normalized KNN Graph Laplacian
kc1 <- kmeans(UKnn,3)
table(kc1$cluster)
plot(A,pch=kc1$cluster,col=kc1$cluster)
title( "Normalized KNN Graph Spectral Clustering")
legend( "topleft" , legend = c("Cluster 1","Cluster 2","Cluster 3"), col = c(1,2,3), text.col = c(1,2,3), pch = c(1,2,3))

# Spectral graph clustering for Normalized Gaussian Kernel Graph Laplacian
kc2 <- kmeans(UGk,3)
table(kc2$cluster)
plot(A,pch=kc2$cluster,col=kc2$cluster)
title( "Normalized Gaussian Kernel Graph Spectral Clustering")
legend( "topleft" , legend = c("Cluster 1","Cluster 2","Cluster 3"), col = c(1,2,3), text.col = c(1,2,3), pch = c(1,2,3))
