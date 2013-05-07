###############################################################################
#Stochastic patch occupancy model (continuous time using Gillespie algorithm)
###############################################################################

#Same model as contSPOM.r with the following changes:
  #Accepts a steady state threshold
  #Able to store and output connectivity and colonization probabilities 
  #Altered to replace patch level loop with vector operations


diseaseSPOM<-function(distance, quality, initial, parms, r, timesteps,ss.threshold){

#Inputs:
#distance = nxn matrix of distances between patches
#quality = nx1 vector of patch qualities
#initial = nx1 vector of initial patch states
#parms = dataframe of parameter values
          #im = scaling parameter for effect of target patch quality in immigration
          #b = scaling parameter for how patch quality affects the number of colonists produces
          #D = inverse of mean dispersal distance
          #y = minimum connectivity required for colonization prob to equal 0.5 = strength of allee effect
          #e = extinction probability of unit patch
          #x = how strongly extinction risk depends on patch quality
          #theta = factor scaling extinction rate of infected patches
          #delta = probability of direct infection
          #gamma0 = probability of reservoir infection immediately following extinction of an infected patch
          #r = decay rate of reservoirs
          #alpha = probability of transition from reservoir to empty patch
          
im<-parms$im
b<-parms$b
D<-parms$D
x<-parms$x
es<-parms$es
ei<-parms$ei
delta<-parms$delta
gamma0<-parms$gamma0
alpha<-parms$alpha

  
n<-length(initial)

#State matrix
state<-matrix(nrow=timesteps+1, ncol=n)
state[1,]<-initial

#Event times vector
event.times<-vector(mode="numeric",length=timesteps+1)

#Reservoir creation times
ti<-vector(mode="numeric", length=n)

#Infectivity vector
gammas<-vector(mode="numeric",length=n)

#Pairwise connectivity matrix
connectivity<-(quality^im%*%t(quality^b))*exp(-D*distance)
connectivity<-connectivity*distance

#Extinction vectors
extinctionS<-vector(mode="numeric",length=n)
extinctionI<-vector(mode="numeric",length=n)

#Calculates extinction rates
for(i in 1:n){
  extinctionS[i]<-es/quality[i]^x}
for(i in 1:n){
  extinctionI[i]<-ei/quality[i]^x}


#Generates occupancy vector for S and I
occ.S<-vector(mode="numeric",length=timesteps)
occ.S[1]<-length(which(initial=="S"))/n

occ.I<-vector(mode="numeric",length=timesteps)
occ.I[1]<-length(which(initial=="I"))/n

#Sets the steady state index to a large value
ssS.index<-100
ssI.index<-0

#Time stepping process
t<-1

#Runs model until steady state or timesteps is reached
while((ssS.index>ss.threshold|ssI.index>ss.threshold) & t<=timesteps){

  #Calculates connectivities to susceptible and infected patches
  connectivityS<-connectivity%*%(state[t,]=="S")
  connectivityI<-connectivity%*%(state[t,]=="I")
    
  
  #Updates infectivity of reservoirs 
  gammas[which(ti>0)]<-gamma0*exp(-r*(t-ti[which(ti>0)]))
  
  #Creates empty rates data frame
  rates<-matrix(nrow=8*n,ncol=6)
  rates<-as.data.frame(rates)
  
  #Fills first column with patch IDs
  rates[,1]<-c(1:n)
  
  #Rate of S to I from direct contact
  rates[1:n,2]<-"S"
  rates[1:n,3]<-"I"
  rates[1:n,4]<-(delta*connectivityI)*(state[t,]=="S")
  
  #S to E/R
  EorR<-c(1:n)
  EorR[which(gammas==0)]<-"E"
  EorR[which(gammas>0)]<-"R"
  rates[(n+1):(2*n),2]<-"S"
  rates[(n+1):(2*n),3]<-EorR
  rates[(n+1):(2*n),4]<-extinctionS*(state[t,]=="S")

  #I to R
  rates[(2*n+1):(3*n),2]<-"I"
  rates[(2*n+1):(3*n),3]<-"R"
  rates[(2*n+1):(3*n),4]<-extinctionI*(state[t,]=="I")

  #R to S
  rates[(3*n+1):(4*n),2]<-"R"
  rates[(3*n+1):(4*n),3]<-"S"
  rates[(3*n+1):(4*n),4]<-connectivityS*(state[t,]=="R")
  
  #R to I
  rates[(4*n+1):(5*n),2]<-"R"
  rates[(4*n+1):(5*n),3]<-"I"
  rates[(4*n+1):(5*n),4]<-connectivityI*(state[t,]=="R")
  
  #E to S
  rates[(5*n+1):(6*n),2]<-"E"
  rates[(5*n+1):(6*n),3]<-"S"
  rates[(5*n+1):(6*n),4]<-connectivityS*(state[t,]=="E")
  
  #E to I
  rates[(6*n+1):(7*n),2]<-"E"
  rates[(6*n+1):(7*n),3]<-"I"
  rates[(6*n+1):(7*n),4]<-connectivityI*(state[t,]=="E")

  #S to I via reservoir
  rates[(7*n+1):(8*n),2]<-"S"
  rates[(7*n+1):(8*n),3]<-"I"
  rates[(7*n+1):(8*n),4]<-gammas*(state[t,]=="S")
  
  
  #Trims unused rows and rows with 0 probability events from rates matrix
  rates<-rates[-which(rates[,4]==0),]
    
  #############################################
  #Can do with sapply
  #Draws a random waiting time from an exp for each event
  for(i in 1:length(rates[,1])){
    rates[i,6]<-rexp(1,rate=rates[i,4])
    }
     
  #Finds the time until the first event 
  time.advance<-min(rates[,6],na.rm=TRUE)
    
  #Advances the time by that ammount
  event.times[t+1]<-event.times[t]+time.advance
     
  #Identifies which event happens first 
  event<-which.min(rates[,6])
    
  #Copies state at time t into state at time t+1
  state[t+1,]<-state[t,]
    
  #Updates state at t+1 with the event that occurs
  state[t+1,rates[event,1]]<-rates[event,3]
    
  #Ends script if extinction has occurred.
  if(length(which(state[t+1,]=="S"))==0 & length(which(state[t+1,]=="I"))==0){
    return(cbind(event.times,state))}
    
  #If the event creates a reservoir, drops the event time into ti, the reservoir
  #creation time vector, for that patch  
  if(rates[event,2]=="I" & rates[event,3]=="R") 
     ti[rates[event,1]]<-t
    
		
  #If susceptibles colonize a reservoir without infection, erases the patch's 
  #infectivity (optional)
  #if(rates[event,3]=="S" & rates[event,2]=="R") 
  #  ti[rates[event,1]]<-0
  #  gammas[rates[event,1]]<-0
    
  #If reservoir becomes an empty patch, erases the patch's infectivity (optional)
  #if(rates[event,3]=="E" & rates[event,2]=="R") 
  #  ti[rates[event,1]]<-0
  #  gammas[rates[event,1]]<-0
    
  #Updates occupancy vector
  occ.S[t+1]<-length(which(state[t+1,]=="S"))/n
	occ.I[t+1]<-length(which(state[t+1,]=="I"))/n
	
	#Updates time index
	t<-t+1
	
	#Updates steady state index = (max - min) over last 200 time steps
	if(t>500){
		ssS.index<-max(occ.S[(t-200):t])-min(occ.S[(t-200):t])
		ssI.index<-max(occ.I[(t-200):t])-min(occ.I[(t-200):t])
		}
	
	#if(t>500){
		#ssS.index<-abs(mean(occ.S[(t-200):(t-100)])-mean(occ.S[(t-100):t]))
		#ssI.index<-abs(mean(occ.I[(t-200):(t-100)])-mean(occ.I[(t-100):t]))
	#}

	#print(c(t,ssS.index,ssI.index,occ.S[t],occ.I[t]))
	
}#End of time loop

output<-cbind(event.times,state)
return(output)

}#End of function
