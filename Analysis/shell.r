
###########################################################################################
#Script to run SPOM over different habitat quality variances and environmental longevities
#and write relevant data to csv.
###########################################################################################

library(foreach)
library(doSNOW)
library(rlecuyer)

source("Code/modelrun3.r")
source("Code/contSPOM(rates).r")

#Set the number of replicates
replicates <- 100

#Generating vectors of longevity and variance values to use
longevity <- seq(-1, 0.5, length.out = 10)
longevity <- 10 ^ longevity

range <- c(0.2, 1.8)

xi_im = c(0, 0.5)
xi_em = c(0, 0.5)

#Generates nx2 matrix of all pairwise longevity variance combinations
parspace<-expand.grid(longevity, xi_im, xi_em)
colnames(parspace) <- c("longevity", "xi_im", "xi_em")

#Generates a new matrix with each long-var combination repeated for the set number of replicates
par.reps<-parspace[rep(1:length(parspace[,1]),replicates),]

#Assigns a unique simulation ID and random number seed for each replicate
simID<-c(1:length(par.reps[,1]))
par.reps<-cbind(par.reps,simID)


#Defines fixed inputs for diseaseSPOM
parms<-data.frame("D" = 5, 
                  "alpha" = 1,
                  "es" = 0.1,
                  "nu" = 0.2,
                  "delta" = 0.3,
                  "gamma0" = 0.5)
  
#Fully connected matrix
distance<-matrix(1,nrow=100,ncol=100)
diag(distance)<-0

#Lattice matrix (on a torus)
# library(igraph)
# distance<-get.adjacency(graph.lattice(c(10,10),directed=F,circular=T))
# distance<-as.matrix(distance)


initial<-c(1:100)
initial[1:50]<-"S"
initial[51:100]<-"E"

timesteps<-5000

#Setting up "parallel backend"
w<-makeCluster(2,type="SOCK")

clusterSetupRNG(w,seed = c(8293, 49866, 71530,  50190, 30330, 13589))

registerDoSNOW(w)

#Checks that the number of workers is set up correctly.
getDoParWorkers()

#Looping through each parameter combo and calling modelrun
system.time(out<-foreach(i = 1:length(par.reps[,1]),.verbose=TRUE,.combine="rbind") %dopar% {
  parms["xi_im"] <- par.reps[i, 2]
  parms["xi_em"] <- par.reps[i, 3]
  modelrun(par.reps[i,1],range,par.reps[i,4],parms,distance,initial,timesteps)
  })

stopCluster(w)






		
		
	
	
