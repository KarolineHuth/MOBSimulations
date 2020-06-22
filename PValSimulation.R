source("~/Documents/PhD/Programming/IsingNetworkTree/R:/IsingNetworkTree.R")
source("~/Documents/PhD/Programming/MOBSimulation/R/DataSimulation.R")
library("readr")
library("networktree")
library("MASS")

##### CHECK P-VAL ####

# Setting values
n <- c(250, 500, 750, 1000)
p <- c(5, 10, 15)
prob <- .2
delta_interaction <- 0 
res_Ising <- data.frame(iter = numeric(1000*length(n) *length(p)))
counter <- 0

for(pi in 1:length(p)) {
  for(ni in 1:length(n)) {
    for(i in 1:1000) {
      counter <- counter + 1
      #### IsingMOB ####
      # Simulate Data
      data <- IsingSimulationSplit(n[ni], p[pi], prob, delta_interaction)
      
      nodevars <- as.data.frame(data[, 2:(p[pi]+1)])
      splitvars <- as.data.frame(data[, 1])
      
      # putting output in text-file
      capture.output(IsingNetworkTree(nodevars, splitvars, model = c("correlation"), verbose = TRUE, minsplit = (n[ni]-50)/2), file = "out.txt")
      out <- readr::read_delim("~/Documents/PhD/Programming/MOBSimulation/out.txt", delim = " ")
      index <- which(grepl("p.value", t(out)))
      
      p.val <- as.vector(t(out))[index[1]+1]
      
      res_Ising$iter[counter] <- i
      res_Ising$n[counter] <- n[ni]
      res_Ising$p[counter] <- p[pi]
      res_Ising$model[counter] <- "Ising"
      res_Ising$pval[counter] <- p.val
    }
  }
}

res_Ising$pval <- as.numeric(res_Ising$pval)

res_Ising %>%
  ggplot(aes(x = pval, color = as.factor(n))) +
  geom_density() + 
  facet_wrap(n ~ p )


par(mfrow = c(3,4))
for(pi in p){
  for(ni in n) {
    sub <- res_Ising %>%
      filter(n == ni) %>%
      filter(p == pi)
    plot(ecdf(sub$pval), main = paste0("n=", ni, "p= ", pi), xlim = c(0, 1))
    lines(c(0,1), c(0, 1), col = "red")
  }
}


write_csv(res_Ising, "res_Ising_sim1000.csv")

#### GGM-MOB ####

n <-  c(250, 500, 750, 1000, 5000)
p <-  c(5, 10, 15)
prob <- .2
delta_interaction <- 0 
res_GGM <- data.frame(iter = numeric(1000*length(n) *length(p)))
counter <- 0

for(pi in 1:length(p)) {
  for(ni in 1:length(n)) {
    for(i in 1:1000) {
      counter <- counter + 1
      
      data <- GGMSimulationSplit(p[pi], n[ni], prob, delta_interaction)
      nodevars <- as.data.frame(data[, 2:(p[pi]+1)])
      splitvars <- as.data.frame(data[, 1])
      
      # putting output in text-file
      capture.output(networktree(nodevars, splitvars, method = c("mob"), model = "correlation", verbose = TRUE, minsplit = (n[ni]-50)/2), file = "out.txt")
      out <- readr::read_delim("~/Documents/PhD/Programming/MOBSimulation/out.txt", delim = " ")
      index <- which(grepl("p.value", t(out)))
      
      p.val <- as.vector(t(out))[index+1]
      
      res_GGM$iter[counter] <- i
      res_GGM$n[counter] <- n[ni]
      res_GGM$p[counter] <- p[pi]
      res_GGM$model[counter] <- "GGM"
      res_GGM$pval[counter] <- p.val
    }
  }
}

# Quick Visualization of the density of the p-values
res_GGM$pval <- as.numeric(res_GGM$pval)

res_GGM %>%
  ggplot(aes(x = pval, color = as.factor(n))) +
  geom_density() + 
  facet_wrap(n ~ p )

# plot cumulative distribution function
par(mfrow = c(3,5))
for(pi in p){
  for(ni in n) {
    sub <- res_GGM_prelim %>%
      filter(n == ni) %>%
      filter(p == pi)
    plot(ecdf(sub$pval), main = paste0("n=", ni, "p= ", pi))
    lines(c(0,1), c(0, 1), col = "red")
  }
}

res_GGM_prelim <- rbind(res_GGM %>% filter(iter != 0), res_GGM_prelim )
res_GGM <- res_GGM_prelim
# not p 15 n 250 
# not p 

write_csv(res_GGM, "res_GGM_sim1000.csv")

#### GGM-MOB - Jonas ####

n <- 250 #c(250, 500, 750, 1000)
p <- 5#c(5, 10, 15)
prob <- .2
delta_interaction <- 0 
res_GGM_Jonas <- data.frame(iter = numeric(1000*length(n) *length(p)))
counter <- 0

for(pi in 1:length(p)) {
  for(ni in 1:length(n)) {
    for(i in 1:1000) {
      counter <- counter + 1
      
      data <- GGMSimulationSplitJonas(p[pi], n[ni], prob, delta_interaction)
      nodevars <- as.data.frame(data[, 2:(p[pi]+1)])
      splitvars <- as.data.frame(data[, 1])
      
      # putting output in text-file
      capture.output(networktree(nodevars, splitvars, method = c("mob"), model = "correlation", verbose = TRUE, minsplit = (n[ni]-50)/2), file = "out.txt")
      out <- readr::read_delim("~/Documents/PhD/Programming/MOBSimulation/out.txt", delim = " ")
      index <- which(grepl("p.value", t(out)))
      
      p.val <- as.vector(t(out))[index+1]
      
      res_GGM_Jonas$iter[counter] <- i
      res_GGM_Jonas$n[counter] <- n[ni]
      res_GGM_Jonas$p[counter] <- p[pi]
      res_GGM_Jonas$model[counter] <- "GGM"
      res_GGM_Jonas$pval[counter] <- p.val
    }
  }
}

# Quick Visualization of the density of the p-values
res_GGM_Jonas$pval <- as.numeric(res_GGM_Jonas$pval)

res_GGM_Jonas %>%
  ggplot(aes(x = pval, color = as.factor(n))) +
  geom_density() + 
  facet_wrap(n ~ p )

# plot cumulative distribution function
par(mfrow = c(3,4))
for(pi in p){
  for(ni in n) {
    sub <- res_GGM_Jonas %>%
      filter(n == ni) %>%
      filter(p == pi)
    plot(ecdf(sub$pval), main = paste0("n=", ni, "p= ", pi), xlim = c(0,1))
    lines(c(0,1), c(0, 1), col = "red")
  }
}



#### LinearModel-MOB ####

n <- c(50, 100, 250, 500, 1000)
p <- c(3, 5, 8)
prob <- .2
delta_interaction <- 0 
res_reg <- data.frame(iter = numeric(1000* length(n) *length(p)))
counter <- 0

for(pi in 1:length(p)){
  for(ni in 1:length(n)) {
    for(i in 1:1000) {
      counter <- counter + 1
      
      data <- GGMSimulationSplit(p[pi], n[ni], prob, delta_interaction)
      nodevars <- as.data.frame(data[, 2:(p[pi]+1)])
      
      # putting output in text-file
      form <- paste("y1", "~",paste(colnames(nodevars[,2:p[pi]]), collapse=" + "), "| z1")
      form <- as.formula(form)
      capture.output(mob(form,
                         control = mob_control(minsplit = 10, verbose = TRUE), data = data,
                         model = linearModel), file = "out.txt")
      out <- readr::read_delim("~/Documents/PhD/Programming/MOBSimulation/out.txt", delim = "-") # does not produce same output as previous
      out <- strsplit(as.character(out$X1), ' ') # to get into right output
      
      p.val <- out[[4]]
      p.val <- p.val %>% as.numeric() 
      index <- which(!is.na(p.val))
      
      res_reg$iter[counter] <- i
      res_reg$n[counter] <- n[ni]
      res_reg$p[counter] <- p[pi]
      res_reg$model[counter] <- "Regression"
      res_reg$pval[counter] <- p.val[index]
    }
  }
}

# Quick Visualization of the density of the p-values
res_reg$pval <- as.numeric(res_reg$pval)

# plot cumulative distribution function
par(mfrow = c(3,5))
for(pi in p){
  for(ni in n) {
    sub <- res_reg %>%
      filter(n == ni) %>%
      filter(p == pi)
    plot(ecdf(sub$pval), main = paste0("n=", ni, "pred= ", pi-1), xlim = c(0,1))
    lines(c(0,1), c(0, 1), col = "red")
  }
}

test <- res_reg %>% filter(p > 3)



# EXAMPLE
mob(y1 ~ y2 + y3 | z1,
    control = mob_control(minsplit = 10, verbose = TRUE), data = data,
    model = linearModel)
summary(fit)
