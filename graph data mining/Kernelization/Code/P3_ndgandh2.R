#############################################
# Name: Nisarg Gandhi                       #
# Unity Id: ndgandh2                        #
# Project 3: Kernelization, Kernel Tricks   #
#############################################

library('kernlab')
library('e1071')

source('bad_kmeans.R')
source('bad_pca.R')
source('bad_svm.R')
source('pipelining.R')
source('contingencyTableMetrics.R')
source('twoCrossConfusionMatrixMetrics.R')

set.seed(64)
cluster <- as.matrix(c(rep(1,500),rep(2,700)))

# This function generates a circle with given inner and outer radius, 
# center and no of points in a 2-D space
generateCircle <- function(n,innerRadius,outerRadius,centre=c(0,0)){ 
     x0 <- centre[1] ; y0 <- centre[2] 
     theta <- 2*pi*runif(n) 
     r <- sqrt(runif(n,min =innerRadius,max = outerRadius )) 
     cbind(x=(outerRadius - innerRadius)*r*cos(theta)+x0, y=(outerRadius - innerRadius)*r*sin(theta)+y0) 
}

# This function generates a circle with given inner and outer radius, 
# center and no of points in a 7-D space
generate7DCircle <- function(n,innerRadius,outerRadius,centre=c(0,0,0,0,0,0,0)){ 
  x0 <- centre[1] ; y0 <- centre[2] ; z0 <- centre[3] ; a0 <- centre[4] ; b0 <- centre[5] ; c0 <- centre[6] ; d0 <- centre[7]
  theta <- 2*pi*runif(n) 
  r <- sqrt(runif(n,min =innerRadius,max = outerRadius )) 
  cbind(x=(outerRadius - innerRadius)*r*cos(theta)+x0, y=(outerRadius - innerRadius)*r*sin(theta)+y0,z=(outerRadius - innerRadius)*r*sin(theta)+z0,a=(outerRadius - innerRadius)*r*sin(theta)+a0,b=(outerRadius - innerRadius)*r*sin(theta)+b0,c=(outerRadius - innerRadius)*r*sin(theta)+c0,d=(outerRadius - innerRadius)*r*sin(theta)+d0) 
}

# This function generates a Sphere with given radius, center and 
# no of points in a 2-D space
generate2DSphere <- function(n,radius,centre=c(0,0)){ 
  x0 <- centre[1] ; y0 <- centre[2] 
  theta <- 2*pi*runif(n) 
  r <- sqrt(runif(n)) 
  cbind(x=radius*r*cos(theta)+x0, y=radius*r*sin(theta)+y0) 
}

# This function generates a Sphere with given radius, center and 
# no of points in a 7-D space
generate7DSphere <- function(n,radius,centre=c(0,0,0,0,0,0,0)){ 
  x0 <- centre[1] ; y0 <- centre[2] ; z0 <- centre[3] ; a0 <- centre[4] ; b0 <- centre[5] ; c0 <- centre[6] ; d0 <- centre[7]
  theta <- 2*pi*runif(n) 
  r <- sqrt(runif(n)) 
  cbind(x=radius*r*cos(theta)+x0, y=radius*r*sin(theta)+y0,z=radius*r*sin(theta)+z0,a=radius*r*sin(theta)+a0,b=radius*r*sin(theta)+b0,c=radius*r*sin(theta)+c0,d=radius*r*sin(theta)+d0) 
}

# This function generates a confusion matrix for given two clusters
getConfusionMatrix <- function(cluster1,cluster2,labels){
  
  if(cluster2[1] == 2)
    x = 1
  else
    x = 2
  
  if(cluster1[1] != cluster2[1]){
    cluster1 <- as.matrix(c(rep(cluster2[1],500),rep(x,700)))
  }
  
  
  tp <- 0
  fp <- 0
  tn <- 0
  fn <- 0
  
 
  for(i in 1:length(cluster1)){
    if(cluster1[i] == labels[1] && cluster2[i] == labels[1])
      tp <- tp + 1
    else if(cluster1[i] == labels[2] && cluster2[i] == labels[1])
      fp <- fp + 1
    else if(cluster1[i] == labels[2] && cluster2[i] == labels[2])
      fn <- fn + 1
    else if(cluster1[i] == labels[1] && cluster2[i] == labels[2])
      tn <- tn + 1
  }
  
  matrix(c(tp,tn,fp,fn),nrow = 2,ncol = 2)
  
}

kernelize <- function(){

# Generates bad_kmeans dataset, performs un-kernelized and kernelized kmeans 
# on the bad_kmeans data and outputs the corresponding performance metrics
bad_kmeans()

# Generates bad_pca dataset, performs un-kernelized and kernelized pca 
# on the bad_pca data and outputs the corresponding performance metrics
bad_pca()

# Generates bad_svm dataset, performs un-kernelized and kernelized svm 
# on the bad_svm data and outputs the corresponding performance metrics
bad_svm()

# Generates bad_kmeans dataset for 7-D space, performs un-kernelized kmeans, 
# applies kernel pca on the bad_kmeans data and performs kmenas and kernelized kmeans 
# on the bad_kmeans data and outputs the corresponding performance metrics
pipelining()


}

