
###########################################################################################
#Script to run SPOM over different habitat quality variances and environmental longevities
#and write relevant data to csv.
###########################################################################################

#Set the number of replicates
replicates<-100

#Generating vectors of longevity and variance values to use
longevity<-seq(20,200,by=20)
variance<-seq(0.02,0.2,by=0.02)

#Generates nx2 matrix of all pairwise longevity variance combinations
parspace<-expand.grid(longevity,variance)
colnames(parspace)[1]<-"longevity"
colnames(parspace)[2]<-"variance"


#Generates a new matrix with each long-var combination repeated for the set number of replicates
par.reps<-parspace[rep(1:length(parspace[,1]),replicates),]

#Assigns a unique simulation ID and random number seed for each replicate
simID<-c(1:length(par.reps[,1]))
par.reps<-cbind(par.reps,simID)



#Defines fixed inputs for diseaseSPOM
  parms<-data.frame("im"=0.5, #Controls how a patch imports colonists
                  "b"=0.5,  #Controls how a patch exports colonists      
                  "D"=5, 
                  "y"=10,
                  "x"=1,
				          "es"=0.1,
                  "ei"=0.5,
                  "delta"=0.5,
                  "gamma0"=0.5,
                  "alpha"=0)
  
  #Fully connected matrix
	distance<-matrix(1,nrow=100,ncol=100)
	diag(distance)<-0

  #Lattice matrix (on a torus)
#   distance<-get.adjacency(graph.lattice(c(10,10),directed=F,circular=T))
#   distance<-as.matrix(distance)


	initial<-c(1:100)
	initial[1:50]<-"S"
	initial[51:100]<-"E"

	timesteps<-5000

#Calling requisite libraries for parallel computing
library(foreach)
library(doSNOW)

#Setting up "parallel backend"
w<-makeCluster(2,type="SOCK")

clusterSetupRNG(w,seed = c(8293, 49866, 71530,  50190, 30330, 13589))

registerDoSNOW(w)

#Checks that the number of workers is set up correctly.
getDoParWorkers()

#Looping through each parameter combo and calling modelrun
system.time(out<-foreach(i = 1:length(par.reps[,1]),.verbose=TRUE,.combine="rbind") %dopar% 
  modelrun(par.reps[i,1],par.reps[i,2],par.reps[i,3],parms,distance,initial,timesteps))

stopCluster(w)






		
		
	
	
