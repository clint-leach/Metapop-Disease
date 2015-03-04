
###########################################################################################
#Script to run SPOM over different habitat quality variances and environmental longevities
#and write relevant data to csv.
###########################################################################################

library(doSNOW)
library(foreach)
library(rlecuyer)

source("Code/highlow.runs.r")
source("Code/contSPOM(rates).r")

#Set the number of replicates
replicates <- 100

#Setting longevity
#longevity <- seq(-1, 0.5, length.out = 10)
longevity <- c(-0.5, 0, 0.5)
longevity <- 10 ^ longevity

# Setting transmission
delta <- seq(0, 0.9, by = 0.1)

# Setting disease-induced mortality
nu <- seq(0.1, 1, by = 0.1)

#Generates a new matrix with each long-var combination repeated for the set number of replicates
treatment <- expand.grid(rep(c("+high","+low"), replicates), longevity, delta, nu)

#Defines fixed inputs for diseaseSPOM
parms<-data.frame("xi_im" = 0.5, 
                  "xi_em" = 0.5,      
                  "D" = 5, 
                  "alpha" = 1,
                  "es" = 0.1,
                  "gamma0" = 0.5)

# Connectivity matrix

# Lattice
# library(igraph)
# distance<-get.adjacency(graph.lattice(c(10,10),directed=F,circular=T))
# distance<-as.matrix(distance)

#Fully connected
distance <- matrix(1, nrow=100, ncol=100)
diag(distance) <- 0


#Initial conditions
initial <- c(1:100)
initial[1:50] <- "S"
initial[51:100] <- "E"

timesteps <- 5000

# Setting up parallelization
w <- makeCluster(11, type="SOCK")

clusterSetupRNG(w, seed = c(8294, 49867, 71531,  50191, 30331, 13590))

registerDoSNOW(w)

getDoParWorkers()

#Looping through each parameter combo and calling modelrun
out <- foreach(i = 1:dim(treatment)[1], .verbose=TRUE, .combine="rbind") %dopar% 
  modelrun(treatment[i, 2], treatment[i, 3], treatment[i, 4], parms, distance, initial, timesteps, treatment[i, 1])

stopCluster(w)










