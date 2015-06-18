
modelrun <- function(longevity, delta, nu, parms, distance, initial, timesteps, treatment){
  # Function that runs SPOM model for given pathogen parameters and habitat quality
  # distribution and collects metapopulation-level simulation results.
  #
  # Takes in pathogen parameters and type of patch quality distribution, 
  # runs the SPOM without infection until steady state, then introduces infection
  # on a randomly selected occupied patch, calls SPOM model and generates output.
  #
  # Args:
  #   longevity: value for pathogen longevity in environment, half-life of infectivity
  #   delta: value for transmission probability
  #   nu: value for disease-induced mortality (reduction in population size from infection)
  #   parms: named dataframe of parameter values
  #   distance: nxn between patch distance matrix
  #   initial: character vector giving initial state of each patch ("S", "I", or "E")
  #   timesteps: number of timesteps to run SPOM for
  #   treatment: character value giving type of quality distribution (one of "low.var", "+high", "+low", "full")
  #
  # Returns:
  #   1 x 11 dataframe giving processed results from simulation
  #     treatment: patch quality distribution used
  #     longevity: pathogen environmental longevity
  #     S: final susceptible occupancy
  #     I: final infectious occupancy
  #     R: final reservoir occupancy
  #     S.pop: mean effective susceptible population size (sum(A_i * p_i(S)))
  #     I.pop: mean effective infectious population size (sum(nu * A_i * p_i(I)))
  #     maxI: maximum infectious occupancy
  #     tmaxI: time to maximum infectious occupancy
  #     quality0: quality of initially infected patch
  
  # Setting boundaries for patch quality distributions
  out.max <- 0.5 * sqrt(12 * 0.2) + 1
  out.min <- 2 - out.max
  
  in.max <- 0.5 * sqrt(12 * 0.02) + 1
  in.min <- 2 - in.max
  
  # Generates quality vector with desired range
  
  # High quality (0.75 to 1.77)
  if(treatment == "+high"){
    quality <- runif(100, min=in.min, max=out.max)}
  
  # Low quality (0.23 to 1.25)
  if(treatment == "+low"){
    quality <- runif(100, min=out.min, max=in.max)}
  
  # Outer range (0.23 to 1.77)
  if(treatment == "full"){
    quality <- runif(100, min=out.min, max=out.max)}
  
  # Inner range (0.75 to 1.25)
  if(treatment == "low.var"){
    quality <- runif(100, min=in.min, max=in.max)
  }

  # Adds nu and delta to parameter vector
  parms["nu"] <- nu
  parms["delta"] <- delta
  
  # Converts longevity (half-life in units of occupancy time) to decay rate in natural time
  long.nat <- longevity / parms$es
  r <- log(2) / long.nat
  
  # Runs simulation without infection until it reaches steady state
  Sonly <- diseaseSPOM(distance, quality, initial, parms, r, timesteps, 0.15)
  Sonly <- na.omit(Sonly)
  initial.inf <- Sonly[length(Sonly[, 1]), 2:101]
   
  # Randomly infects an occupied patch
  initial.inf[sample(which(initial.inf == "S"), 1)] <- "I"
  
  # Stores quality of initial infected patch
  initial.quality <- quality[which(initial.inf == "I")]
  
  # Calls SPOM to run with infection
  output <- diseaseSPOM(distance, quality, initial.inf, parms, r, timesteps, 0.05)
  output <- na.omit(output)
  
  runtime <- dim(output)[1]
  
  # Processing output
  state <- output[(runtime - 499):runtime, -1]
  susc <- state == "S"
  inf <- state == "I"
  empty <- state == "E" | state == "R"
  
  S <- rowSums(susc) / 100
  I <- rowSums(inf) / 100
  E <- rowSums(empty) / 100
  
  # Effective population size (assuming pop proportional to quality)
  S.pop <- (colSums(susc) / (dim(state)[1])) %*% quality
  I.pop <- (colSums(inf) / (dim(state)[1])) %*% (parms$nu * quality)
  
  # Generating and collecting output
  output <- data.frame(treatment, longevity, delta, nu, S[length(S)], I[length(I)], E[length(E)],
                           S.pop, I.pop, max(I), output[which.max(I),1], initial.quality)
  
  names(output) <- c("treatment", "longevity", "delta", "nu", "S", "I", "R", "S.pop", "I.pop", "maxI", "tmaxI", "quality0")
  
  return(output)
  
}