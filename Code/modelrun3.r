#############################################################################################
#modelrun function that takes in longevity and variance values, runs the SPOM without infection
#until steady state, then introduces infection on a randomly selected occupied patch and runs SPOM for
#an additional 2000 time steps.  Collects the relevant occupancy data from diseaseSPOM output
#and writes to csv file. Same as modelrun, but with number of replicates as an input and steady
#state threshold.

#Changes as of 01/2012
  #Proportion of time in each state calculations corrected 12/20/11
  #Number of state transitions for each patch added
  #Unique simulation ID added to both metapop and patch-level data
  #Patch quality drawn from uniform instead of gamma distribution
  #Initial infected patch has to be occupied
#############################################################################################

modelrun<-function(longevity, variance, k, parms, distance, initial, timesteps){
  
	#Generates quality vector with desired variance
	#quality<-rgamma(100,shape=1/variance,scale=variance
  max<-0.5*sqrt(12*variance)+1
  min<-2-max
  quality<-runif(100,min=min,max=max)
	
	#Converts longevity (half-life) to decay rate
	r<-log(2)/longevity
	
  #Runs simulation without infection until it reaches steady state
	Sonly<-diseaseSPOM(distance, quality, initial, parms, r, timesteps,0.15)
	Sonly<-na.omit(Sonly)
	initial.inf<-Sonly[length(Sonly[,1]),2:101]
  
  #Randomly infects an occupied patch
	initial.inf[sample(which(initial.inf=="S"),1)]<-"I"
	
	#Stores quality of initial infected patch
	initial.quality<-quality[which(initial.inf=="I")]
	
	#Calls SPOM to run with infection
	output<-diseaseSPOM(distance, quality, initial.inf, parms, r, timesteps,0.05)
	output<-na.omit(output)
	
	#Generates occupancy vectors for each state
	S<-vector(mode="numeric",length=length(output[,1]))
	I<-vector(mode="numeric",length=length(output[,1]))
	R<-vector(mode="numeric",length=length(output[,1]))
	
	for(i in 1:length(S)){
		S[i]<-length(which(output[i,]=="S"))/100}

	for(i in 1:length(I)){
		I[i]<-length(which(output[i,]=="I"))/100}

	for(i in 1:length(R)){
		R[i]<-length(which(output[i,]=="R"))/100}

	
	#Metaopulation-level output
	#Occupancy at the final time step
	final.occupancy<-c(S[length(S)],I[length(I)],R[length(R)])
	#Average occupancy over last 100 time steps
	avg.occupancy<-c(mean(S[(length(S)-100):length(S)]),mean(I[(length(I)-100):length(I)]),mean(R[(length(I)-100):length(I)]))
	#Maximum infection occupanc and time at which reached (as a measure of R0)
	max.infected<-c(max(I), output[which.max(I),1])
	
	#Writes all metapop data to a single vector
	metapop.data<-c(variance, longevity, k, final.occupancy, avg.occupancy, max.infected,initial.quality)
	
	#Patch-level output
	patch.data<-matrix(nrow=length(output[1,])-1,ncol=55)
	
	timesteps<-length(output[,1])
	
	for(i in 1:length(output[1,])-1){
		#Patch quality
		patch.data[i,1]<-quality[i]
		#Time to infection (can be inf)
		patch.data[i,2]<-output[min(which(output[,i+1]=="I"),na.rm=TRUE),1]
		#Proportion of time steps in each state
		patch.data[i,3]<-length(which(output[,i+1]=="S"))/(timesteps)
		patch.data[i,4]<-length(which(output[,i+1]=="I"))/(timesteps)
		patch.data[i,5]<-length(which(output[,i+1]=="R"))/(timesteps)
		patch.data[i,6]<-length(which(output[,i+1]=="E"))/(timesteps)
    
		#Proportion of time steps i each state, broken into 10 intervals
		#Susceptible
		patch.data[i,7]<-length(which(output[1:(floor(timesteps/10)),i+1]=="S"))/floor(timesteps/10)
		patch.data[i,8]<-length(which(output[(floor(timesteps/10)+1):(2*floor(timesteps/10)),i+1]=="S"))/floor(timesteps/10)
		patch.data[i,9]<-length(which(output[((2*floor(timesteps/10))+1):(3*floor(timesteps/10)),i+1]=="S"))/floor(timesteps/10)
		patch.data[i,10]<-length(which(output[((3*floor(floor(timesteps/10)))+1):(4*(floor(timesteps/10))),i+1]=="S"))/(floor(timesteps/10))
		patch.data[i,11]<-length(which(output[((4*floor(timesteps/10))+1):(5*(floor(timesteps/10))),i+1]=="S"))/(floor(timesteps/10))
		patch.data[i,12]<-length(which(output[((5*floor(timesteps/10))+1):(6*(floor(timesteps/10))),i+1]=="S"))/(floor(timesteps/10))
		patch.data[i,13]<-length(which(output[((6*floor(timesteps/10))+1):(7*(floor(timesteps/10))),i+1]=="S"))/(floor(timesteps/10))
		patch.data[i,14]<-length(which(output[((7*floor(timesteps/10))+1):(8*(floor(timesteps/10))),i+1]=="S"))/(floor(timesteps/10))
		patch.data[i,15]<-length(which(output[((8*floor(timesteps/10))+1):(9*(floor(timesteps/10))),i+1]=="S"))/(floor(timesteps/10))
		patch.data[i,16]<-length(which(output[((9*floor(timesteps/10))+1):timesteps,i+1]=="S"))/(timesteps-((9*floor(timesteps/10))))
		
		#Infected
		patch.data[i,17]<-length(which(output[1:(floor(timesteps/10)),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,18]<-length(which(output[(floor(timesteps/10)+1):(2*(floor(timesteps/10))),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,19]<-length(which(output[((2*floor(timesteps/10))+1):(3*(floor(timesteps/10))),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,20]<-length(which(output[((3*floor(timesteps/10))+1):(4*(floor(timesteps/10))),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,21]<-length(which(output[((4*floor(timesteps/10))+1):(5*(floor(timesteps/10))),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,22]<-length(which(output[((5*floor(timesteps/10))+1):(6*(floor(timesteps/10))),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,23]<-length(which(output[((6*floor(timesteps/10))+1):(7*(floor(timesteps/10))),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,24]<-length(which(output[((7*floor(timesteps/10))+1):(8*(floor(timesteps/10))),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,25]<-length(which(output[((8*floor(timesteps/10))+1):(9*(floor(timesteps/10))),i+1]=="I"))/(floor(timesteps/10))
		patch.data[i,26]<-length(which(output[((9*floor(timesteps/10))+1):timesteps,i+1]=="I"))/(timesteps-((9*floor(timesteps/10))))

		#Reservoir
		patch.data[i,27]<-length(which(output[1:(floor(timesteps/10)),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,28]<-length(which(output[(floor(timesteps/10)+1):(2*(floor(timesteps/10))),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,29]<-length(which(output[((2*floor(timesteps/10))+1):(3*(floor(timesteps/10))),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,30]<-length(which(output[((3*floor(timesteps/10))+1):(4*(floor(timesteps/10))),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,31]<-length(which(output[((4*floor(timesteps/10))+1):(5*(floor(timesteps/10))),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,32]<-length(which(output[((5*floor(timesteps/10))+1):(6*(floor(timesteps/10))),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,33]<-length(which(output[((6*floor(timesteps/10))+1):(7*(floor(timesteps/10))),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,34]<-length(which(output[((7*floor(timesteps/10))+1):(8*(floor(timesteps/10))),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,35]<-length(which(output[((8*floor(timesteps/10))+1):(9*(floor(timesteps/10))),i+1]=="R"))/(floor(timesteps/10))
		patch.data[i,36]<-length(which(output[((9*floor(timesteps/10))+1):timesteps,i+1]=="R"))/(timesteps-((9*floor(timesteps/10))))
		
		#Empty
		patch.data[i,37]<-length(which(output[1:(floor(timesteps/10)),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,38]<-length(which(output[(floor(timesteps/10)+1):(2*(floor(timesteps/10))),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,39]<-length(which(output[((2*floor(timesteps/10))+1):(3*(floor(timesteps/10))),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,40]<-length(which(output[((3*floor(timesteps/10))+1):(4*(floor(timesteps/10))),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,41]<-length(which(output[((4*floor(timesteps/10))+1):(5*(floor(timesteps/10))),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,42]<-length(which(output[((5*floor(timesteps/10))+1):(6*(floor(timesteps/10))),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,43]<-length(which(output[((6*floor(timesteps/10))+1):(7*(floor(timesteps/10))),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,44]<-length(which(output[((7*floor(timesteps/10))+1):(8*(floor(timesteps/10))),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,45]<-length(which(output[((8*floor(timesteps/10))+1):(9*(floor(timesteps/10))),i+1]=="E"))/(floor(timesteps/10))
		patch.data[i,46]<-length(which(output[((9*floor(timesteps/10))+1):timesteps,i+1]=="E"))/(timesteps-((9*floor(timesteps/10))))
		
    #Number of transitions for each patch
    transitions<-matrix(0,nrow=1,ncol=3)

      for(j in 1:(length(output[,1])-1)){
        if(output[j,i+1]!=output[j+1,i+1]){
          if(output[j+1,i+1]=="E"|output[j+1,i+1]=="R"){
            transitions[1,1]<-transitions[1,1]+1}
          if((output[j,i+1]=="E"|output[j,i+1]=="R")&(output[j+1,i+1]=="S"|output[j+1,i+1]=="I")){
            transitions[1,2]<-transitions[1,2]+1}
          if(output[j,i+1]=="S"& output[j+1,i+1]=="I"){
            transitions[1,3]<-transitions[1,3]+1}
        }
      }
    
    #Stores number of extinctions, colonizations, and infections for each patch
    patch.data[i,47:49]<-transitions
    
    #Replicate ID and steady state
    patch.data[i,50]<-k
    patch.data[i,51:53]<-final.occupancy
    
    #Longevity, variance
		patch.data[i,54]<-longevity
		patch.data[i,55]<-variance
    
    
		}
	
  #Convert patch data from a character to a numeric matrix
  patch.data<-apply(patch.data,1,as.character)
  patch.data<-apply(patch.data,1,as.numeric)
  
	#Write metapop.data and patch.data to different csv files and figure out how to "append"
 	write.table(t(metapop.data),file="metapop(final).csv",append=TRUE,col.names=FALSE,sep=",")
	
  return(patch.data)
	}