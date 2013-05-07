#############################################################################################
#Code to run SPOM model with just high, just low, or all quality patches to evaluate the relative
#contribution of each quality class to pandemic dynamics.
#############################################################################################

modelrun<-function(longevity, parms, distance, initial, timesteps, treatment){
  
  out.max<-0.5*sqrt(12*0.2)+1
  out.min<-2-out.max
  
  in.max<-0.5*sqrt(12*0.02)+1
  in.min<-2-in.max
  
  #Generates quality vector with desired variance
  if(treatment=="+high"){
    quality<-runif(100,min=in.min,max=out.max)}
  
  if(treatment=="+low"){
    quality<-runif(100,min=out.min,max=in.max)}
  
  if(treatment=="full"){
    quality<-runif(100,min=out.min,max=out.max)}
  
  if(treatment=="low.var"){
    quality<-runif(100,min=in.min,max=in.max)
  }

  
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
  
  for(i in (length(S)-100):length(S)){
    S[i]<-length(which(output[i,]=="S"))/100}
  
  for(i in (length(I)-100):length(I)){
    I[i]<-length(which(output[i,]=="I"))/100}
  
  for(i in (length(R)-100):length(R)){
    R[i]<-length(which(output[i,]=="R"))/100}
  
  
  #Generating and collecting output
  
  #Occupancy at the final time step
  final.occupancy<-c(S[length(S)],I[length(I)],R[length(R)])
  
  #Average occupancy over last 100 time steps
  avg.occupancy<-c(mean(S[(length(S)-100):length(S)]),mean(I[(length(I)-100):length(I)]),mean(R[(length(I)-100):length(I)]))
  
  #Maximum infection occupanc and time at which reached (as a measure of R0)
  max.infected<-c(max(I), output[which.max(I),1])
  
  #Writes all metapop data to a single vector
  metapop.data<-data.frame(treatment, longevity, S[length(S)],I[length(I)],R[length(R)],
                           mean(S[(length(S)-100):length(S)]),mean(I[(length(I)-100):length(I)]),mean(R[(length(I)-100):length(I)]),
                           max(I), output[which.max(I),1],initial.quality)
  

  return(metapop.data)
}