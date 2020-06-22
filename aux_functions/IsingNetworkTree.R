library("party")
library("GeneNet")
library("mvtnorm")
library("MASS")
library("corpcor")
library("networktree")
library("strucchange")
library("Formula")
library("tidyverse")
library("IsingSampler")
library("profvis")

# -----------------------------------------------------------------------------------------

source("~/Documents/PhD/Programming/IsingNetworkTree/R:/B Ising - EMVS - Functions.R")


# -----------------------------------------------------------------------------------------

#### STEP 1: Set up Ising fit function 
####### - needs to have objects 
####### - (1) coef (coefficients, estimated parameters), 
####### - (2) logLik (maximized log-likelihood),
####### - (3) emprical estimating function (score function)


IsingFit <- function(y, x = NULL, start = NULL, weights = NULL,
                     offset = NULL, model = c("correlation", "main"), ...,
                     estfun = TRUE, object = FALSE) {
  
  nodevars <- y 
  
  p <- ncol(nodevars) # number of nodes
  n <- nrow(nodevars) # number of observations
  varnames <- colnames(nodevars)
  
  nodevars <- as.matrix(nodevars)
  
  ### put dots in a list
  dotlist <- list(...)
  
  ### check if correlation matrix is identified
  if(n <= p*(p-1)/2) {
    stop("Isingfit: n < k*(k-1)/2, correlation matrix is not identified.")
  }
  
  
  # --- STEP 1: 
  
  fit <- tryCatch(maximum.pseudoposterior(x = as.matrix(nodevars), prior.variance = Inf),
                  error = function(e) e, warning = function(war) war)
  val <- colSums(nodevars) # Checking for low variablity leading to errors
  if( any(val < 2) ){
    if(inherits(fit, c("error", "warning"))) {
      
      loglik <- -Inf # Ising model can't be estimated, thus Log-Likelihood set to infinity;
      
      # Just open some empty vectors, otherwise R produces errors
      scores <- c()
      coef <- c()
      vc <- NULL
      
      warning('Ising model cannot be estimated in (sub-)model, presumably because of low variability')
      # Sacha's Function
      # EIfit <- EstimateIsing(nodevars, method = "pl")
      # 
      # main <- EIfit$thresholds
      # 
      # R <- EIfit$graph
      # inter <- as.vector(EIfit$graph)[as.vector(lower.tri(EIfit$graph))]
    }
    if(!inherits(fit, c("error", "warning"))) {
      
      warning('Structural change test cannot be computed for (sub-)model, because of low variability')
      
      # Setup objects to avoid errors 
      loglik <- -Inf
      scores <- c()
      vc <- NULL
      
      # Get parameter estimates for submodel
      coef <- c(diag(fit$sigma), as.vector(fit$sigma)[as.vector(lower.tri(fit$sigma))])
      ynam <- if(is.null(varnames)) 1L:p else varnames
      objnames <- c(paste0("main_", ynam), combn(p, 2, function(x) paste0("cor_", x[1], "_", x[2])))
      
      id <- NULL
      if(any("main"        == model)) id <- c(id, 1:p)
      if(any("correlation" == model)) id <- c(id, 1:(p*(p-1)/2) + p)
      
      coef <- coef[id]
      objnames <- objnames[id]
      names(coef) <- objnames
    }
    
  } else {
    # Maartens Function 
    # main effect of variables on diagonal
    main <- diag(fit$sigma)
    
    # interactions on off-diagonal
    R <- fit$sigma
    inter <- as.vector(fit$sigma)[as.vector(lower.tri(fit$sigma))]
    diag(R) <- 0
    
    # --- STEP 2: Compute matrixes to make actual computation easier
    nodevars <- t(nodevars)
    
    M <- R%*%nodevars # compute the sums 
    A <- exp(main + M)
    D <- 1 + A
    E <- A/ D
    
    # --- STEP 3: Compute Log Likelihood of Ising function
    S <- main + M
    loglik <- sum(nodevars*S) - sum(log(D))
    
    # --- STEP 4: Score Functions
    
    # Score Function for main effect
    score_main <- nodevars - E # 0 up to the 3rd decimal
    
    
    # Score Function for interaction 
    score_rho <- combn(p, 2,
                       function(x) (2*(nodevars[x[1] ,] * nodevars[x[2], ]) - (E[x[1], ] * nodevars[x[2], ]) - (E[x[2], ] * nodevars[x[1], ])) )
    
    # --- STEP 5: Store objects depending which measures to consider
    
    # Set-Up
    coef <- c(main, inter)
    scores <- cbind(t(score_main), score_rho)
    ynam <- if(is.null(varnames)) 1L:p else varnames
    objnames <- c(paste0("main_", ynam), combn(p, 2, function(x) paste0("cor_", x[1], "_", x[2])))
    
    id <- NULL
    if(any("main"        == model)) id <- c(id, 1:p)
    if(any("correlation" == model)) id <- c(id, 1:(p*(p-1)/2) + p)
    
    coef <- coef[id]
    scores   <- scores[, id]
    objnames <- objnames[id]
    
    # Naming the coefficients and scores
    names(coef) <- objnames
    colnames(scores) <- objnames
    
    # --- STEP 6: Compute inverse of Hessian
    index <- rbind(1: (p * (p - 1)/2), combn(p, 2))
    inverse.hessian <- invert.hessian(sigma.map = fit$sigma, index = t(index), 
                                      x = t(nodevars), prior.variance = Inf)
    vc <- - inverse.hessian
  }
  
  if(estfun == FALSE) {
    scores <- NULL
  }
  
  
  res <- list(coefficients = coef,
              objfun = -loglik,
              estfun = scores,
              object = vc)
  
  return(res)
}

#fit <- IsingFit(y = nodevars, estfun = TRUE, model = c("main", "correlation"))

#-----------------------------------------------------

#### STEP 2: Feed function to mob 
####### - fit is the self-defined fit function 
####### - control is pre-defined control function 
####### - 

IsingNetworkTree <- function(nodevars, splitvars, #Data 
                             type = c("cor", "pcor", "glasso"),
                             model = c("main", "correlation"), # Define which parameters to split by
                             alpha = .05, bonferroni = TRUE, minsplit = "simple" ,  # Control for MOB  
                             verbose = TRUE) {
  
  # Rename column names of the node and split variables
  nodevars <- nodevars
  splitvars <- splitvars
  
  if(is.null(colnames(nodevars))){
    if(ncol(nodevars) == 1){
      colnames(nodevars) <- "nodevars1"
    } else {
      colnames(nodevars) <- paste('nodevars',1:ncol(nodevars),sep="")
    }
  }
  
  if(is.null(colnames(splitvars))){
    if(ncol(splitvars) == 1){
      colnames(splitvars) <- "splitvars1"
    } else {
      colnames(splitvars) <- paste('splitvars',1:ncol(splitvars),sep="")
    }
  }
  
  # Create dataframe used by the MOB function 
  d <- data.frame(nodevars,splitvars)
  
  # Write the formula; everything behind the ~ is used as a partitioning variable
  form <- paste(paste(colnames(nodevars), collapse=" + "), "~",paste(colnames(splitvars), 
                                                                     collapse=" + "))
  form <- as.formula(form)
  
  # Define control variables for the MOB
  if(minsplit == "simple") {
    p <- ncol(nodevars)
    minsplit <- (p*(p-1)) # !!!! Further define what the minimum sample size needed is for a reliable estimate
  }
  control <- partykit::mob_control(alpha = alpha, bonferroni = bonferroni, 
                                   minsplit = minsplit, verbose = verbose)
  control$ytype <- "matrix"
  
  # Run MOB function
  res <- partykit::mob(formula = form, data = d, fit = IsingFit,  
                       model = model, control = control)
  
  # Change class of 
  class(res) <- c("networktree", "mob_networktree", type[1], class(res))
  
  # Print results
  return(res)
}

# data <- IsingDataSimulation(n.nds = 10, n.obs = 1000, n.zvar = 1, delta.cor = 0.3,
#                             delta.all = FALSE, mean.dif = FALSE, mean.change = .6,
#                             part.type = "binary")
# # 
# # 
# nodevars <- as.data.frame(data[, 2:11])
# splitvars <- as.data.frame(data[, 1])
# 
# test <- IsingNetworkTree(nodevars = nodevars, splitvars = splitvars, model = c("correlation"))
# test


