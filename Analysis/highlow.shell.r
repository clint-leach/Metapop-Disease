
# Script to explore the effect of habitat quality distribution and pathogen 
# traits on epidemic dynamics and host population size.
#
# Runs replicated SPOM simulations over all combinations of:
#    habitat quality distribution (high or low quality)
#    pathogen longevities (3 values)
#    direct transmission probabilities (10 values)
#    disease-induced mortalities (10 values)
#
# The connectivity matrix (either fully connected or lattice) must be set manually
#
# Produces results stored in "Output/pathgrid.RData" and "Output/pathgrid(lattice).RData" 
# that are used to generate Figures 1 - 3, and Supplemental Figures 1 - 2


library(doSNOW)
library(foreach)
library(rlecuyer)

source("Code/highlow.runs.r")
source("Code/contSPOM(rates).r")

#===============================================================================
# Setting up parameters for full-factorial simulations

# Set the number of replicates
replicates <- 100

# Setting longevity
longevity <- c(-0.5, 0, 0.5)
longevity <- 10 ^ longevity

# Setting transmission
delta <- seq(0, 0.9, by = 0.1)

# Setting disease-induced mortality
nu <- seq(0.1, 1, by = 0.1)

# Generates full factorial matrix with each combination repeated for the set number of replicates
treatment <- expand.grid(rep(c("+high","+low"), replicates), longevity, delta, nu)

#===============================================================================
# Defines fixed inputs for diseaseSPOM
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

# Time steps
timesteps <- 5000

#===============================================================================
# Setting up parallelization

w <- makeCluster(11, type="SOCK")

clusterSetupRNG(w, seed = c(8294, 49867, 71531,  50191, 30331, 13590))

registerDoSNOW(w)

getDoParWorkers()

#===============================================================================
#Looping through each parameter combo and calling modelrun
out <- foreach(i = 1:dim(treatment)[1], .verbose=TRUE, .combine="rbind") %dopar% 
  modelrun(treatment[i, 2], treatment[i, 3], treatment[i, 4], parms, distance, initial, timesteps, treatment[i, 1])

stopCluster(w)
