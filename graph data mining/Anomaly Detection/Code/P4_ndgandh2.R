#########################################################################
# Name : Nisarg Gandhi                                                  #
# Unity Id: ndgandh2                                                    #
# Project 4: Anomaly Detection                                          #
# Paper : Event Detection in Time Series of Mobile Communication Graphs #
#########################################################################

rm(list=ls(all = T))

# Checks if the required packages are installed. I f not, install the packages listed
packages <- c("igraph","geometry","gtools","egonet")

for(pkg in packages){
  if(!is.element(pkg, installed.packages()[,1]))
  {install.packages(pkg, repos="http://cran.fhcrc.org")
  } else {print(paste(pkg, " library already installed"))}
}

library('egonet')
library('igraph')
library('geometry')
library('gtools')

findAnomaly <- function(graphDir, graphName){
  
  # This function generates a correlation matrix between every pair of nodes for a given feature 
  # for each time window, calculates largest eigenvectors of each correlation matrix and finally calculates
  # a similarity score Z-Score, to detect if a particular node in the timeseries is an anomaly or not.
  anomalyDetection <- function(feature,graph,feature_name){
    # T * N matrix for a given feature
    feature <- feature[,2:ncol(feature)]
    featureTN <- t(as.matrix(feature))
    
    # Generates correlation matrix between a pair of nodes for the first time window of 7 days
    C1 <- cor(x = featureTN[1:7,])
    C1[which(is.na(C1))] <- 0
    
    # Calculating eigenvectors for the above correlation matrix
    eigC1 <- eigen(C1)
    R <- eigC1$vectors[,1]
    
    # Runs the loop for the remaining time windows of 7 days
    i <- 2
    z_score <- 1
    while((i+6) <= nrow(featureTN)){
      # Correlation matrix
      C2 <- cor(x = as.matrix(featureTN[i:(i+6),]))
      C2[which(is.na(C2))] <- 0
      
      # Largest eigenvector for the above Correlation  matrix
      eigC2 <- eigen(C2)
      U_t <- eigC2$vectors[,1]
      
      # Calculating typical eigen behavior of the past time windows
      if(ncol(as.matrix(R)) > 1){
        if(i < 7)
          R_t <- rowMeans(R[,1:(i-1)])
        else
          R_t <- rowMeans(R[,(i-7):(i-1)])
      }
      else
        R_t <- R
      
      # Calculates a Z-Score(similarity score) for the current time window based on the 
      # past eigen behavior
      Z <- 1 - (t(as.matrix(U_t)) %*% as.matrix(R_t))
      z_score <- c(z_score,Z)
      R <- cbind(R,U_t)
      i <- i + 1
    }
    
    # Calculation moving range average of the timeseries
    MR <- abs(diff(z_score))
    SigmaMR <- mean(MR)
    
    # Calculating median value of timeseries
    Median <- median(z_score)
    
    # Calculating upper threshold which is given by "median + (3 * MR)"
    upper_threshold <- Median + (3 * SigmaMR)
    lower_threshold <- Median - (3 * SigmaMR)
    
    # Plotting the timeseries, with anomalies plotted in red and other normal nodes in black
    col <- rep("black",nrow(featureTN))
    col[which(z_score > upper_threshold)] <- "red"
    plot(z_score, type = "h",col = col,ylim = c(0,2),xlab = "Timeseries[day]", ylab = "z-score similarity",main = paste(graph, feature_name,sep = " : "))
    
    # Horizontal red line denotes the upper threshold
    abline(h = upper_threshold,col=2,lty = 3)
    text(x = 500,y=1.8,labels = paste("th = ",round(upper_threshold,digits = 3),sep = ""))
    # Returning the list of anomaly detailed information
    list(z_score = z_score,col = col, threshold = upper_threshold, moving_range = MR, moving_range_avg = SigmaMR, median = Median)
  }
  
  # This function writes the anomalies in a files in the required format
  # First line in the file gives the total anomalies detected, followed by the anomalies detected.
  # If anomalies are less than 10, all anomalies are written, else if more than 10 but less than 100 than top 10 anomalies are written
  # else if more than 100, top 10% anomalies are written.
  writeAnomaly <- function(filename, anomalyData){
    noOfAnomalies <- length(which(anomalyData$z_score > anomalyData$threshold))
    seq <- seq(from = 1,to = length(anomalyData$z_score))
    z_frame <- cbind(seq,anomalyData$z_score)
    o <- order(z_frame[,2],z_frame[,1],decreasing = T)
    z_frame <- z_frame[o,]
    anomalies <- z_frame[which(z_frame[,2] > anomalyData$threshold),1]
    sink(filename)
    cat(noOfAnomalies)
    cat("\n")
    if(noOfAnomalies > 10 && noOfAnomalies <= 100)
      printAnomaly <- anomalies[1:10]
    else if(noOfAnomalies > 100)
      printAnomaly <- anomalies[1:(0.1 * noOfAnomalies)]
    else
      printAnomaly <- anomalies
    for(i in 1:length(printAnomaly)){
      cat(printAnomaly[i])
      cat("\n")
    }
    sink()
  }
  
  # Read graphs in the grapg directory
  TSGraphDir <- graphDir
  graphsByDay <- mixedsort(list.files(path = TSGraphDir,pattern = "*",full.names = T))
  
  n <- length(graphsByDay)
  
  # Generate graphs and extract features for each nodes in the graph
  graphData <- read.table(graphsByDay[1])
  for(i in 2:n){
    graphData <- rbind(graphData, read.table(graphsByDay[i]))
  }
  
  vertices <- unique(c(unique(graphData[,1]),unique(graphData[,2])))
  
  nodeDegree <- c(NA,length(vertices))
  nodeEgonetEdges <- c(NA,length(vertices))
  nodeClusterCoeff <- c(NA,length(vertices))
  for(i in 1:n){
    graphData <- read.table(graphsByDay[i]) + 1
    vertices <- unique(c(unique(graphData[,1]),unique(graphData[,2])))
    max_vertices <- c(1:max(vertices))
    G <- graph.empty() + vertices(max_vertices)
    G <- add.edges(G,t(graphData))
    G <- as.directed(G) # Directed graph
    
    #G_dash <- G - max_vertices[-vertices]
    #plot(G_dash,main = i)
    
    # Adjacency matrix of the generated graph
    AdjacencyMatrix <- as.matrix(get.adjacency(G))
    
    # Extract features for the given nodes
    
    # Degree of the given N nodes
    nodeDegree <- cbind(nodeDegree,degree(G))
    
    # Clustering Coefficient of the given nodes
    nodeClusterCoeff <- cbind(nodeClusterCoeff,transitivity(graph = G,type = "local", vids = V(G),isolates = "zero"))
    
    # Number of edges in the egonet of the node
    nodeEgonetEdges <- cbind(nodeEgonetEdges,sapply(ego.extract(AdjacencyMatrix,ego = V(G)),sum)/2)
    
  }
  anomaly1 <- anomalyDetection(feature = nodeDegree,graphName, "Degree of the node")
  writeAnomaly("degreeofnode_anomaly.txt",anomaly1)
  anomaly2 <- anomalyDetection(feature = nodeClusterCoeff,graphName, "Clustering Coefficient of the node")
  writeAnomaly("clusterCoeff_anomaly.txt",anomaly2)
  anomaly3 <- anomalyDetection(feature = nodeEgonetEdges,graphName, "Number of edges in the egonet of the node")
  writeAnomaly("egonet_anomaly.txt",anomaly3)
}