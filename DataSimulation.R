library("IsingSampler")
library("Formula")
library("Matrix")
library("party")
library("GeneNet")
library("mvtnorm")
library("MASS")

#### Simulate Data ####
####Ising ####

IsingSimulationSplit <- function(n, p, prob, 
                                 delta_interaction, # c(0.15, 0.3, 0.6) # values taken from Haslbeck 2020
                                 thresh = 0){
  
  nedges <- p*(p-1)/2 # number of edges
  
  selection <- sample(0:1, size = nedges, prob = c(1-prob, prob), replace  = TRUE)
  sel_par <- runif(sum(selection), min = -.6, max = 1.2) # min and max values taken from Haslbeck 2020
  selection[selection==1] <- sel_par
  
  threshold <- rep(thresh, p) # Threshold values
  #Interaction Matrix
  graph <- matrix(0, p, p)
  graph[lower.tri(graph)] <- selection
  graph[upper.tri(graph)] <- t(graph)[upper.tri(graph)] 
  
  
  
  if(delta_interaction == 0) {
    
    # Sample Data
    data <- IsingSampler(n, graph, thresholds = threshold)
    
  } else {
    
    # Sample Data
    data1 <- IsingSampler(n/2, graph, thresholds = threshold)
    
    # Add the change in interaction
    change <- sample(0:1, size = nedges, prob = c(1-prob, prob), replace  = TRUE)
    change[change == 1] <- sample(c(+ delta_interaction, - delta_interaction), size = sum(change), replace = TRUE)
    
    #Interaction Matrix
    graph2 <- matrix(0, p, p)
    graph2[lower.tri(graph2)] <- change
    graph2[upper.tri(graph2)] <- t(graph2)[upper.tri(graph2)] 
    graph2 <- graph + graph2
    
    # Sample Data
    data <- rbind(data1, IsingSampler(n/2, graph2, thresholds = threshold))
  }
  
  
  # Generate a splitting variable/moderator variable
  splitVar <- c(rep(0, times = n/2), rep(1, times = n/2))
  
  data <- cbind(splitVar, data)
  
  return(data)
}

# data <- IsingSimulationSplit(1000, 10, 0.2, 0)


#### GGM ####

GGMSimulationSplit <- function(p, n, prob, delta_interaction){
  
  nedges <- p*(p-1)/2 # number of edges
  
  posdef <- FALSE
  
  while(posdef == FALSE) {
    # Step A:
    # generate random network with x nodes, simulate pcor
    graph <- ggm.simulate.pcor(p, etaA=prob)
    diag(graph) <- 1
    
    if (delta_interaction == 0) {
      
      # Simulate only one big dataset in case of no correlation difference
      pcor.inv <- pseudoinverse(graph)
      mean <- numeric(p)
      
      d <- mvtnorm::rmvnorm(n, mean = mean, sigma = pcor.inv)
      posdef <- TRUE
    } else {
      
      # If in variable correlation setting, modify the partial correlation matrix
      change <- sample(0:1, size = nedges, prob = c(1-prob, prob), replace  = TRUE)
      change[change == 1] <- sample(c(+ delta_interaction, - delta_interaction), size = sum(change), replace = TRUE)
      
      #Interaction Matrix
      graph2 <- matrix(0, p, p)
      graph2[lower.tri(graph2)] <- change
      graph2[upper.tri(graph2)] <- t(graph2)[upper.tri(graph2)] 
      graph2 <- graph + graph2
      
      
      if (!is.positive.definite(graph2, tol=1e-8)){
        posdef <- FALSE
      } else {
        posdef <- TRUE
      }
      
      pcor.inv <- pseudoinverse(graph)
      pcor.inv2 <- pseudoinverse(graph2)
      mean <- numeric(p)
      
      # Sample Data 
      d <- rbind(
        mvtnorm::rmvnorm(n/2, mean = mean,
                         sigma = pcor.inv),
        mvtnorm::rmvnorm(n/2, mean = mean,
                         sigma = pcor.inv2)
      )
    }
    
  }
  
  # Generate a splitting variable/moderator variable
  z1 <- c(rep(0, times = n/2), rep(1, times = n/2))
  
  
  data <- cbind(as.data.frame(z1) , d)
  
  colnames(data)[2:(1+p)] <- paste0("y", 1:p)
  
  return(data)
}

# data <- GGMSimulationSplit(15, 1000, delta_interaction = 0, prob = .5)

