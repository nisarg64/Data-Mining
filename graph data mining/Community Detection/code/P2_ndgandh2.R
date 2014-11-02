# Author : Nisarg Gandhi
# Unity Id: ndgandh2
# Paper 3: Large Scale Spectral Clustering on Graphs

rm(list = ls())

# This function implements the Efficient Spectral Clustering on Graphs with Regeneration (ESCG-R)
# arguments: 
#           graphFile: Name of file whose communities are to be detected
#           supernodes: numerical value which defines the number of supernodes to be sampled for ESCG-R algorithm
#           iterations: numerical value which defines the number of times the supernodes are to be regenerated for efficiency
#           communities: numerical value which defines the number of communities to be detected
ESCGR <- function(graphFile, supernodes, iterations, communities){
  
  library('igraph')
  library('cluster')
  library('fpc')
  library('MASS')
  
  # This function generates a graph from the vertex id pair vectors along with assigning weights to the edges
  getGraph <- function(args){   
    graph_sub <- graphData[2:nrow(graphData),] + 1	# Adding 1 to move vertex id, because of 0 vertex-id
    G <- graph.empty() + vertices(1:NV)
    G <- add.edges(G, t(graph_sub))
    G <- as.undirected(G)      # making the graph undirected
    
    #Assign weight = 1 to the edges
    E(G)$weight <- runif(length(E(G)),max = 1,min = 1)
    G
  }
  
  # This function samples d supernodes and assigns regular nodes to the super nodes based on the shortest path 
  # computed by the Dijkstra's algorithm. Then d disjoints subsets are connected tothe corresponding supernodes
  # and R binary matrix is generated which defines their connectivity
  algo1GenerateRMatrix <- function(graph){
    
    # Sampling Random d supernodes from the vertices
    getSuperNodes <- function(){  
      SuperNodes <- as.integer(sample(V(graph),d))
      SuperNodes
    }
    
    # Computing the shortest path of all regular nodes from super nodes
    getShortestPathMatrix <- function(){
      SPMat <- shortest.paths(graph,v = SuperNodes, to = as.integer(V(graph)), weights = E(graph)$weight,algorithm = "dijkstra")
      SPMat
    }
    
    SuperNodes <- getSuperNodes()
    SPMat <- getShortestPathMatrix()
    
    # Partition nodes into d disjoint subsets based on shortest paths from supernodes toregular nodes
    NodeClass <- rep(0,NV)
    for (i in 1:NV){
      min = Inf
      min_d <- 0
      for(j in 1:d){
        if(SPMat[j,i] < min){
          min = SPMat[j,i]
          min_d = j
        }
      }
      NodeClass[i] = min_d
    }
    
    # Computes R binary matrix d*n using d disjoint sets of d supernodes
    R <- matrix(nrow = as.numeric(d),ncol = as.numeric(NV))
    for(i in 1:d){
      for(j in 1:NV){
        if(NodeClass[j] == i && SuperNodes[i] != j)
          R[i,j] = 1
        else
          R[i,j] = 0
      }
    }
    R
  }
  
  # This function regenerates (2k - 2) supernodes for optimal supernode selection efficient clustering of large graphs
  algo2GenerateRMatrix <- function(U){
    
    # Compute R matrix d*n using d disjoint sets of d supernodes
    R <- matrix(0,nrow = as.numeric(as.numeric(2 * as.numeric(k)) - 2),ncol = as.numeric(NV))
    for(i in 2:k){
      U_bar <- mean(U[,i])
      R[ 2 * i - 3, which(U[,i] <= U_bar)] = 1
      R[ 2 * i - 2, which(U[,i] > U_bar)] = 1
    }
    R
  }
  
  # This function generates the Right Singualar Vector from the Singualar Value Decomposition of R binary matrix 
  # for further clustering the large graphs and reduce the graph size
  generateUMatrix <- function(R){
    
    # Computer matrix W_hat indicating the adjacency for transformed bipartite graph
    W_hat <- R %*% W
    
    # Computing Eigen Value Decomposition of Bipartite Graph without computing L' and D' matrices thereby reducing the computation
    
    # D1 is 1/sqrt of diagonal matrix which contains the column sums of the transformed bipartite graph adjacency matrix
    D1 <- 1/sqrt(colSums(as.matrix(W_hat[,1:ncol(W_hat)])) * diag(ncol(W_hat)))
    D1[D1 == Inf] = 0    # Handling the infinite values generated due to disconnected components
    
    # D2 is 1/sqrt of diagonal matrix which contains the row sums of the transformed bipartite graph adjacency matrix
    D2 <- 1/sqrt(rowSums(as.matrix(W_hat[1:nrow(W_hat),])) * diag(nrow(W_hat)))
    D2[D2 == Inf] = 0    # Handling the infinite values generated due to disconnected components
    
    Z <- D2 %*% W_hat %*% D1
    Z_t <- t(as.matrix(Z))  # Transpose of Z

    
    # Singular Value Decomposition
    # Computing the Left Singular Vectors of Z
    ZZT <- Z %*% Z_t
    ZZT_eig <- eigen(ZZT)
    Sigma <- ZZT_eig$values[1:k]
    Y <- ZZT_eig$vectors[,1:k]
    
    # Computing the Right Singular Vectors of Z
    X <- t(as.matrix((as.vector(as.matrix(ginv(Sigma))) * diag(nrow(as.matrix(Sigma)))) %*% t(Y) %*% Z))
    
    # U contains top k columns of X which are used to apply k-means algorithm for clustering
    U <- D1 %*% X
    U
  }
  
  # This function writes the communities detected into a file in a describes format
  writeCommunities <- function(kc){
    
    # Detecting communities in the clusters and writing to a file
    sink("communities.txt")
    for(j in 1:k){
      for(i in 1:NV){
        if(kc$cluster[i] == j ){
          cat(i-1)
          cat(" ")
        }
      }
      cat("\n")
    }
    sink()
  }
  
  #Implementation of EFCG-R algorithm
  # Read the graph containing vertex id pair
  graphData <- read.table(graphFile)
  NV <- graphData$V1[1] # Store the total number of vertices
  d <- supernodes
  iter <- iterations
  k <- communities
  
  # Get graph from the graphdata file
  graph <- getGraph(graphData)
  
  # Generate adjacency matrix for the given input graph
  W <- get.adjacency(graph)
  
  # Generates supernodes for iternation number of times and partition nodes into d disjoint sets according to the 
  # shortest path from supernodes for efficient community detection
  for (i in 1:iter){
    if(i == 1)
      R <- algo1GenerateRMatrix(graph)
    else
      R <- algo2GenerateRMatrix(U)
    
    # Generate right singular vector from the binary matrix and the original graph for further clustering
    U <- generateUMatrix(R)
  }
  
  # Applying K-means clustering to cluster top k column vectors of Right Singular Vector
  kc <- kmeans(U,centers = k)
  
  # Write communities to a file
  writeCommunities(kc)
}