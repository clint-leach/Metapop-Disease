
###########################################################################################
#Script to run SPOM over different habitat quality variances and environmental longevities
#and write relevant data to csv.
###########################################################################################

#Set the number of replicates
replicates<-100

#Setting longevity
longevity<-100


#Generates a new matrix with each long-var combination repeated for the set number of replicates
treatment<-rep(c("+high","+low","full","low.var"), replicates)


#Defines fixed inputs for diseaseSPOM
parms<-data.frame("im"=0.5, 
                  "b"=0.5, 
                  "D"=5, 
                  "x"=1,
                  "es"=0.1,
                  "ei"=0.5,
                  "delta"=0.5,
                  "gamma0"=0.5,
                  "xi"=1,
                  "alpha"=0)

#Connectivity matrix

#Lattice
# library(igraph)
# distance<-get.adjacency(graph.lattice(c(10,10),directed=F,circular=T))
# distance<-as.matrix(distance)

#Fully connected
distance<-matrix(1,nrow=100,ncol=100)
diag(distance)<-0


#Initial conditions
initial<-c(1:100)
initial[1:50]<-"S"
initial[51:100]<-"E"

timesteps<-5000

#Calling requisite libraries for parallel computing
library(doSNOW)

#Setting up "parallel backend"
w<-makeCluster(2,type="SOCK")

clusterSetupRNG(w,seed = c(8294, 49867, 71531,  50191, 30331, 13590))

registerDoSNOW(w)


#Checks that the number of workers is set up correctly.
getDoParWorkers()

#Looping through each parameter combo and calling modelrun
out<-foreach(i = 1:length(treatment),.verbose=TRUE,.combine="rbind") %dopar% modelrun(longevity,parms,distance,initial,timesteps,treatment[i])
stopCluster(w)

names(out)<-c("treatment","longevity","S","I","R","avgS","avgI","avgR","maxI","tmaxI","quality0")








