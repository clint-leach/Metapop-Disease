
diseaseSPOM <- function(distance, quality, initial, parms, r, timesteps, ss.threshold){
  # Runs stochastic patch occupancy model
  #
  # Args:
  #   distance = nxn matrix of distances between patches
  #   quality = nx1 numeric vector of patch qualities
  #   initial = nx1 character vector of initial patch states
  #   parms = dataframe of parameter values
  #     xi_im = scaling parameter for effect of target patch quality in immigration
  #     xi_em = scaling parameter for how patch quality affects the number of colonists produced
  #     D = inverse of mean dispersal distance
  #     es = extinction probability of unit susceptible patch
  #     ei = extinction probability of unit infectious patch
  #     alpha = how strongly extinction risk depends on patch quality
  #     delta = probability of direct infection
  #     gamma0 = probability of reservoir infection immediately following extinction of an infected patch
  #   r = decay rate of reservoirs
  #   timesteps = number of time steps to run
  #   ss.threshold = range within which occupancy can vary over 200 time steps to determine "steady state"
  #                  has been reached
  #
  # Returns:
  #   matrix of dimension timesteps X (n+1) giving the state of each patch at each timestep
  
  xi_im <- parms$xi_im
  xi_em <- parms$xi_em
  D <- parms$D
  alpha <- parms$alpha
  es <- parms$es
  ei <- parms$ei
  delta <- parms$delta
  gamma0 <- parms$gamma0
    
  n <- length(initial)
  
  # State matrix
  state <- matrix(nrow=timesteps+1, ncol=n)
  state[1, ] <- initial
  
  # Event times vector
  event.times <- vector(mode="numeric", length=timesteps+1)
  
  # Reservoir creation times
  ti <- vector(mode="numeric", length=n)
  
  # Infectivity vector
  gammas <- vector(mode="numeric", length=n)
  
  # Pairwise connectivity matrix
  connectivity <- (quality ^ xi_im %*% t(quality ^ xi_em)) * exp(-D * distance)
  connectivity <- connectivity * distance
  
  #Extinction vectors
  extinctionS <- es / quality ^ alpha
  extinctionI <- ei / quality ^ alpha
  
  # Generates occupancy vector for S and I
  occ.S <- vector(mode="numeric", length=timesteps)
  occ.S[1] <- length(which(initial=="S")) / n
  
  occ.I <- vector(mode="numeric", length=timesteps)
  occ.I[1] <- length(which(initial=="I")) / n
  
  # Sets the steady state index to a large value
  ssS.index <- 100
  ssI.index <- 0
  
  # Time stepping process
  t <- 1
  
  # Runs model until steady state or timesteps is reached
  while((ssS.index > ss.threshold | ssI.index > ss.threshold) & t <= timesteps){
    
    # Calculates connectivities to susceptible and infected patches
    connectivityS <- connectivity %*% (state[t, ] == "S")
    connectivityI <- connectivity %*% (state[t ,] == "I")
    
    # Updates infectivity of reservoirs 
    gammas[which(ti > 0)] <- gamma0 * exp(-r * (event.times[t] - ti[which(ti > 0)]))
    
    # Creates empty rates data frame
    rates <- matrix(nrow=8*n, ncol=5)
    rates <- as.data.frame(rates)
    
    # Fills first column with patch IDs
    rates[, 1] <- c(1:n)
    
    # Rate of S to I from direct contact
    rates[1:n, 2] <- "S"
    rates[1:n, 3] <- "I"
    rates[1:n, 4] <- (delta * connectivityI) * (state[t, ] == "S")
    
    # S to E/R
    EorR <- c(1:n)
    EorR[which(gammas == 0)] <- "E"
    EorR[which(gammas > 0)] <- "R"
    rates[(n+1):(2*n), 2] <- "S"
    rates[(n+1):(2*n), 3] <- EorR
    rates[(n+1):(2*n), 4] <- extinctionS * (state[t, ] == "S")
    
    # I to R
    rates[(2*n+1):(3*n), 2] <- "I"
    rates[(2*n+1):(3*n), 3] <- "R"
    rates[(2*n+1):(3*n), 4] <- extinctionI*(state[t, ] == "I")
    
    # R to S
    rates[(3*n+1):(4*n), 2] <- "R"
    rates[(3*n+1):(4*n), 3] <- "S"
    rates[(3*n+1):(4*n), 4] <- connectivityS * (state[t, ] == "R")
    
    # R to I
    rates[(4*n+1):(5*n), 2] <- "R"
    rates[(4*n+1):(5*n), 3] <- "I"
    rates[(4*n+1):(5*n), 4] <- connectivityI * (state[t, ] == "R")
    
    # E to S
    rates[(5*n+1):(6*n), 2] <- "E"
    rates[(5*n+1):(6*n), 3] <- "S"
    rates[(5*n+1):(6*n), 4] <- connectivityS * (state[t, ] == "E")
    
    # E to I
    rates[(6*n+1):(7*n), 2] <- "E"
    rates[(6*n+1):(7*n), 3] <- "I"
    rates[(6*n+1):(7*n), 4] <- connectivityI * (state[t, ] == "E")
    
    # S to I via reservoir
    rates[(7*n+1):(8*n), 2] <- "S"
    rates[(7*n+1):(8*n), 3] <- "I"
    rates[(7*n+1):(8*n), 4] <- gammas * (state[t, ] == "S")
    
    rates[, 5] <- rexp(dim(rates)[1], rate = rates[, 4])
    
    # Finds the time until the first event 
    time.advance <- min(rates[, 5], na.rm=TRUE)
    
    # Advances the time by that ammount
    event.times[t+1] <- event.times[t] + time.advance
    
    # Identifies which event happens first 
    event <- which.min(rates[, 5])
    
    # Copies state at time t into state at time t+1
    state[t+1, ] <- state[t, ]
    
    # Updates state at t+1 with the event that occurs
    state[t+1, rates[event, 1]] <- rates[event, 3]
    
    # Ends script if extinction has occurred.
    if(length(which(state[t+1, ] == "S")) == 0 & length(which(state[t+1, ] == "I")) == 0){
      return(cbind(event.times, state))}
    
    # If the event creates a reservoir, drops the event time into ti, the reservoir
    # creation time vector, for that patch  
    if(rates[event, 2] == "I" & rates[event, 3] == "R") 
      ti[rates[event, 1]] <- event.times[t+1]
        
    # Updates occupancy vector
    occ.S[t+1] <- length(which(state[t+1, ] == "S")) / n
    occ.I[t+1] <- length(which(state[t+1, ] == "I")) / n
    
    # Updates time index
    t <- t + 1
    
    # Updates steady state index = (max - min) over last 200 time steps
    if(t > 500){
      ssS.index <- max(occ.S[(t-200):t]) - min(occ.S[(t-200):t])
      ssI.index <- max(occ.I[(t-200):t]) - min(occ.I[(t-200):t])
    }
       
  }#End of time loop
  
  output<-cbind(event.times,state)
  return(output)

}