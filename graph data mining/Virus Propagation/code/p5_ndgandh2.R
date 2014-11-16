#########################################################################
# Name : Nisarg Gandhi                                                  #
# Unity Id: ndgandh2                                                    #
# Project 4: Virus Propagation                                          # 
#########################################################################

rm(list = ls())
library('igraph')

virusPropagation <- function(graphFile){
  
  # This function generates a graph from the vertex id pair vectors along with assigning weights to the edges
  getGraph <- function(graphData,NV){   
    graph_sub <- graphData[2:nrow(graphData),] + 1  # Adding 1 to move vertex id, because of 0 vertex-id
    G <- graph.empty() + vertices(1:NV)
    G <- add.edges(G, t(graph_sub))
    G <- as.undirected(G)      # making the graph undirected
    
    #Assign weight = 1 to the edges
    E(G)$weight <- runif(length(E(G)),max = 1,min = 1)
    G
  }
  
  # This function calculates the effective strength of the virus propagation model for given beta and delta values
  VPMEffectiveStrength <- function(beta, theta, graph1){
    
    # Calculate the largest eigen value of the graph
    lambda1 <- graph.eigen(graph = graph1,which = list(pos = "LA", howmany = 1))
    # Effective strength
    s <- lambda1$values * (beta / theta)
    list(s = s, lambda1 = lambda1$values)
  }
  
  # This function runs a grid search on a range of beta values to find the minimum beta for which epidemic can be prevented
  betaGridSearch <- function(from, to, by,theta,Epidemic_details){
    s1 <- c()
    Sequence <- seq(from = from,to = to,by = by)
    for(i in 1:length(Sequence)){
      s1 <- c(s1,VPMEffectiveStrength(Sequence[i], theta, graph = graph)$s)
    }
    # Plot the beta variation graph
    plot(Sequence,s1,type = "p", col = "red", main = "Transmission Probability Variation", xlab = "Transmission Probability(Beta)", ylab = "Effective Strength(S)")
    lines(Sequence,rep(1,length(Sequence)),col = "black")
    
    # Calculating the mininmum beta for which epidemic can be prevented
    minBeta <- (min(s1[which(s1 > 1)]) * theta)/Epidemic_details$lambda1
    lines(rep(minBeta,length(s1)),s1)
    text(to - 0.01, 1.2,labels = paste("Min. Beta",minBeta,sep = ":"))
    cat("Minimum transmission probability (β) that results in a network-wide epidemic:",minBeta,"\n")
  }
  
  # This function runs a grid search on a range of delta values to find the maximum delta for which epidemic can be prevented
  thetaGridSearch <- function(from,to,by,beta,Epidemic_details){
    s1 <- c()
    Sequence <- seq(from = from,to = to,by = by)
    for(i in 1:length(Sequence)){
      s1 <- c(s1,VPMEffectiveStrength(beta, Sequence[i], graph = graph)$s)
    }
    # Plot the beta variation graph
    plot(Sequence,s1,type = "p", col = "red", main = "Healing Probability Variation", xlab = "Healing Probability(Delta)", ylab = "Effective Strength(S)")
    lines(Sequence,rep(1,length(Sequence)),col = "black")
    # Calculating the maximum delta for which epidemic can be prevented
    maxDelta <- (Epidemic_details$lambda1 * beta)/min(s1[which(s1 > 1)])
    lines(rep(maxDelta,length(s1)),s1)
    cat("Maximum healing probability (δ) that results in a network-wide epidemic:",maxDelta,"\n")
    maxDelta
  }
  
  # This function tell whether the epidemic will grow or die based on the effective strength of the network
  epidemicStatus <- function(s){
    if(s >= 1)
      cat("Infection will grow\n")
    else
      cat("Infection will die\n")
  }
  
  # Reading and intrepreting the graph
  graphData <- read.table(graphFile)
  NV <- graphData$V1[1]
  # Get graph from the graphdata file
  graph <- getGraph(graphData,NV)
  
  
  # Problem 1 : Effective Strength
  
  beta1 <- 0.20
  theta1 <- 0.70
  
  # Calculating effective strength of the contact network
  Epidemic_details1 <- VPMEffectiveStrength(beta1,theta1,graph = graph)
  cat("Largest Eigen Value:",Epidemic_details1$lambda1,"\n")
  cat("Effective Strength:",Epidemic_details1$s,"\n")
  epidemicStatus(Epidemic_details1$s)
  
  # Analysing the effect of beta and delta on the epidemic growth over the contact network
  betaGridSearch(from = 0,to = 0.05,by = 0.002,theta = theta1,Epidemic_details = Epidemic_details1)
  maxDelta1 <- thetaGridSearch(from = 8,to = 9,by = 0.03,beta = beta1,Epidemic_details = Epidemic_details1)
  text(8.2, 1.01,labels = paste("Max. Delta",round(maxDelta1,digits = 3),sep = ":"))
  
  beta2 <- 0.01
  theta2 <- 0.60
  
  # Calculating effective strength of the contact network
  Epidemic_details2 <- VPMEffectiveStrength(beta2,theta2,graph = graph)
  cat("Largest Eigen Value:",Epidemic_details2$lambda1,"\n")
  cat("Effective Strength:",Epidemic_details2$s,"\n")
  epidemicStatus(Epidemic_details2$s)
  
  # Analysing the effect of beta and delta on the epidemic growth over the contact network
  betaGridSearch(from = beta2,to = 0.04,by = 0.001,theta = theta2,Epidemic_details = Epidemic_details2)
  maxDelta2 <- thetaGridSearch(from = 0.25,to = 0.5,by = 0.01,beta = beta2,Epidemic_details = Epidemic_details2)
  text(0.3, 1.1,labels = paste("Max. Delta",round(maxDelta2,digits = 3),sep = ":"))
  
  # Problem 2 : Simulation
  
  # This function simulates the virus propagation model for 100 time steps based on beta and delta values and given initially infected nodes
  simulation <- function(infected, beta, theta, simulations,graph1){
    j <- 1
    infectedMatrix <- c()
    while(j <= simulations){
      avgInfection <- c()
      for(k in 1:100){
        # Susceptible nodes
        susceptible <- which(infected == "green")
        # Infected nodes
        unsusceptible <- which(infected == "red")
        if(length(unsusceptible) > 0){
          # Finding the neighbors of the infected nodes which are susceptible
          susceptNeighbors <- c()
          for(i in 1:length(unsusceptible)){
            susceptNeighbors <- c(susceptNeighbors,neighbors(graph = graph1,v = unsusceptible[i]))
          }
          susceptNeighbors <- unique(susceptNeighbors[which(infected == "green")])
          # Infecting suscept neighbors with a beta probability
          infected[sample(susceptNeighbors,size = round(beta * length(susceptNeighbors)))] <- "red"
          # Healing infected nodes with a delta probability
          infected[sample(unsusceptible,size = round(theta * length(unsusceptible)))] <- "green"
          # Calculating average infected nodes over a time tick
          avgInfection <- append(avgInfection,length(infected[infected == "red"])/length(V(graph1)))
        }else{
          avgInfection <- append(avgInfection,0)
          
        }
        #cat("Initial Infected:",length(unsusceptible)," | Healed:",round(theta * length(unsusceptible))," | After 1 timetick:",length(infected[infected == "red"]),"\n")
      }
      infectedMatrix <- rbind(infectedMatrix, avgInfection)
      j <- j + 1
    }
    infectedFraction <- colMeans(infectedMatrix)
    # Plotting the fraction of infected nodes over 100 time ticks
    plot(infectedFraction,col = "red",log = "xy",type = "b", ylab = "Fraction of infected nodes", xlab = "Time Ticks",main = "Simulation ")
    cat("Virus Infection is spread across the network ","\n")
  }
  
  # Infecting n/10 nodes initially
  infected <- rep("green",length(V(graph)))
  infected[sample(V(graph),size = length(V(graph))/10)] <- "red"
  cat("Number of initially infected nodes:",floor(length(V(graph))/10),"\n")
  # Simulation for beta1 and delta1 values
  simulation(infected = infected,beta = beta1,theta = theta1,simulations = 10,graph)
  # Simulation for beta2 and delta2 values
  simulation(infected = infected,beta = beta2,theta = theta2,simulations = 10,graph)
  
  
  # Problem 3 Immunization Policy
  
  set.seed(64)
  k1 <- 200
  
  # This function calculates the effective strength of network after immunization
  ImmunizationEffect <- function(G, beta, theta){
    Epidemic_Effect <- VPMEffectiveStrength(beta,theta,graph = G)
    Epidemic_Effect
  }
  
  # This function carries out the grid search to find the minimum number of vaccines required to prevent the epidemic from spreading in the conatct network
  kGridSearch <- function(from,to,by,G,beta,theta, pType){
    Sequence <- seq(from = from,to = to,by = by)
    minK <- Inf
    EpidemicStat <- c()
    # Calculating the effective strength of the contact network after immunizing the network for different values of k
    for(i in 1:length(Sequence)){
      G_dash <- getImmunizedGraph(policyType = pType,k = Sequence[i],G = G)
      s <- ImmunizationEffect(G_dash, beta1, theta1)
      EpidemicStat <- append(EpidemicStat, s$s)
    }
    # Calculating the minimum number of vaccines required to immunize the network from epidemic
    minK <- min(Sequence[EpidemicStat < 1])
    # Plotting the k variation
    plot(Sequence,EpidemicStat,type = "b", col = "red",xlab = "No of Vaccines",ylab = "Effective Strength(s)",main = "Variation of No. of Vaccines")
    lines(Sequence,rep(1,length(Sequence)))
    lines(rep(minK,length(EpidemicStat)),EpidemicStat)
    minK
  }
  
  # Policy 1 Implementation
  p1 <- function(k , G){
    ImmunizeNodes <- sample(V(G),size = k)
    G <- delete.vertices(graph = G,v = ImmunizeNodes)
    G
  }
  
  # Policy 2 Implementation
  p2 <- function(k, G){
    nodeDegree <- degree(graph = G,v = V(G))
    nodeDetails <- cbind(nodeDegree,V(G))
    nodeDetails <- nodeDetails[order(-nodeDetails[,1]),]
    ImmunizeNodes <- nodeDetails[1:k,2]
    G <- delete.vertices(graph = G,v = ImmunizeNodes)
    G
  }
  
  # Policy 3 Implementation
  p3 <- function(k, G){
    for(i in 1:k){
      nodeDeg <- which.max(degree(graph = G,v = V(G)))
      G <- delete.vertices(graph = G,v = nodeDeg)
    }
    G
  }
  
  # Policy 4 Implementation
  p4 <- function(k, G){
    eigVector <- abs(graph.eigen(graph = G,which = list(pos = "LA", howmany = 1))$vectors)
    eigDetails <- cbind(eigVector,V(G))
    eigDetails <- eigDetails[order(-eigDetails[,1]),]
    ImmunizeNodes <- eigDetails[1:k,2]
    G <- delete.vertices(graph = G,v = ImmunizeNodes)
    G
  }
  
  # Generates Immunized graph based on the input policy type
  getImmunizedGraph <- function(policyType, k, G){
    switch(policyType,
           p1 = p1(k ,G),
           p2 = p2(k, G),
           p3 = p3(k, G),
           p4 = p4(k, G))
  }
  
  
  # Policy 1
  G1 <- graph
  # Efective strength of the graph after immunizing using policy 1
  G11 <- getImmunizedGraph(policyType = "p1",k = k1,G = G1)
  s1 <- ImmunizationEffect(G11, beta1, theta1)
  cat("Largest Eigen Value:",s1$lambda1,"\n")
  cat("Effective Strength:",s1$s,"\n")
  epidemicStatus(s1$s)
  
  # Carrying out the grid search to get the minimum number of vaccines required to prevent the epidemic 
  kGridSearch(5000,5500,50,G1,beta1,theta1,pType = "p1")
  
  # Carrying out the simulation of virus propagation model after immunizing k nodes
  infected1 <- rep("green",length(V(G11)))
  infected1[sample(V(G11),size = length(V(G11))/10)] <- "red"
  simulation(infected = infected1,beta = beta1,theta = theta1,simulations = 10,graph1 = G11)
  
  
  
  # Policy 2
  G2 <- graph
  # Efective strength of the graph after immunizing using policy 2
  G22 <- getImmunizedGraph(policyType = "p2", k = k1,G = G2)
  s2 <- ImmunizationEffect(G22, beta1, theta1)
  cat("Largest Eigen Value:",s2$lambda1,"\n")
  cat("Effective Strength:",s2$s,"\n")
  epidemicStatus(s2$s)
  
  # Carrying out the grid search to get the minimum number of vaccines required to prevent the epidemic
  kGridSearch(150,400,50,G2,beta1,theta1,pType = "p2")
  
  # Carrying out the simulation of virus propagation model after immunizing k nodes
  infected2 <- rep("green",length(V(G22)))
  infected2[sample(V(G22),size = length(V(G22))/10)] <- "red"
  simulation(infected = infected2,beta = beta1,theta = theta1,simulations = 10,graph1 = G22)
  
  
  # Policy 3
  G3 <- graph
  # Efective strength of the graph after immunizing using policy 3
  G33 <- getImmunizedGraph(policyType = "p3", k = k1,G = G3)
  s3 <- ImmunizationEffect(G33, beta1, theta1)
  cat("Largest Eigen Value:",s3$lambda1,"\n")
  cat("Effective Strength:",s3$s,"\n")
  epidemicStatus(s3$s)
  
  # Carrying out the grid search to get the minimum number of vaccines required to prevent the epidemic
  kGridSearch(160,240,10,G3,beta1,theta1,pType = "p3")
  
  # Carrying out the simulation of virus propagation model after immunizing k nodes
  infected3 <- rep("green",length(V(G33)))
  infected3[sample(V(G33),size = length(V(G33))/10)] <- "red"
  simulation(infected = infected3,beta = beta1,theta = theta1,simulations = 10,graph1 = G33)
  
  
  
  # Policy 4
  G4 <- graph
  # Efective strength of the graph after immunizing using policy 4
  G44 <- getImmunizedGraph(policyType = "p4", k = k1,G = G4)
  s4 <- ImmunizationEffect(G44, beta1, theta1)
  cat("Largest Eigen Value:",s4$lambda1,"\n")
  cat("Effective Strength:",s4$s,"\n")
  epidemicStatus(s4$s)
  
  # Carrying out the grid search to get the minimum number of vaccines required to prevent the epidemic
  kGridSearch(3500,5500,50,G4,beta1,theta1,pType = "p4")
  
  # Carrying out the simulation of virus propagation model after immunizing k nodes
  infected4 <- rep("green",length(V(G44)))
  infected4[sample(V(G44),size = length(V(G44))/10)] <- "red"
  simulation(infected = infected4,beta = beta1,theta = theta1,simulations = 10,graph1 = G44)
  
}