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
  max <- 0.5 * sqrt(12 * variance) + 1
  min <- 2 - max
  quality <- runif(100, min=min, max=max)
	
	#Converts longevity (half-life) to decay rate
	r <- log(2) / longevity
	
  #Runs simulation without infection until it reaches steady state
	Sonly <- diseaseSPOM(distance, quality, initial, parms, r, timesteps, 0.15)
	Sonly <- na.omit(Sonly)
	initial.inf <- Sonly[length(Sonly[, 1]), 2:101]
  
  #Randomly infects an occupied patch
	initial.inf[sample(which(initial.inf=="S"), 1)] <- "I"
	
	#Stores quality of initial infected patch
	initial.quality <- quality[which(initial.inf == "I")]
	
	#Calls SPOM to run with infection
	output <- diseaseSPOM(distance, quality, initial.inf, parms, r, timesteps, 0.05)
	output <- na.omit(output)
  
  # Processing output
  
  state <- output[, -1]
  susc <- state == "S"
  inf <- state == "I"
  empty <- state == "E" | state == "R"
  
  S <- rowSums(susc) / 100
  I <- rowSums(inf) / 100
  E <- rowSums(empty) / 100
  
  inf.events <- susc[-(timesteps + 1), ] * inf[-1, ]
  susc.col <- empty[-(timesteps + 1), ] * susc[-1, ]
  inf.col <- empty[-(timesteps + 1), ] * inf[-1, ]
  susc.ex <- susc[-(timesteps + 1), ] * empty[-1, ]
  inf.ex <- inf[-(timesteps + 1), ] * empty[-1, ]
  
  output <- data.frame(quality, 
                       tinf = apply(inf, 2, which.max),
                       I = colSums(inf) / (timesteps + 1),
                       S = colSums(susc) / (timesteps + 1),
                       E = colSums(empty) / (timesteps + 1),
                       inf.events.early = colSums(inf.events[c(1:(timesteps / 5)), ]),
                       inf.events.tot = colSums(inf.events),
                       susc.col = colSums(susc.col),
                       inf.col = colSums(inf.col),
                       susc.ex = colSums(susc.ex),
                       inf.ex = colSums(inf.ex),
                       Sfin = S[timesteps + 1],
                       Ifin = I[timesteps + 1],
                       Efin = E[timesteps + 1],
                       maxI = max(I),
                       quality0 = initial.quality,
                       repID = k,
                       longevity = longevity,
                       variance = variance
  )
  

  return(output)
  
	}