
modelrun <- function(longevity, parms, distance, initial, timesteps, treatment){
  # Takes in longevity and  type of patch quality distribution, runs the SPOM without infection
  # until steady state, then introduces infection on a randomly selected occupied patch,
  # calls SPOM model and generates output.
  #
  # Args:
  #   longevity: value for pathogen longevity in environment, half-life of infectivity
  #   parms: named numeric vector of parameter values
  #   distance: nxn between patch distance matrix
  #   initial: character vector giving initial state of each patch
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
  #     avgS: susceptible occupancy averaged over last 100 time steps
  #     avgI: infectious occupancy averaged over last 100 time steps
  #     avgR: reservoir occupancy averaged over last 100 time steps
  #     maxI: maximum infectious occupancy
  #     tmaxI: time to maximum infectious occupancy
  #     quality0: quality of initially infected patch
  
  
  out.max <- 0.5 * sqrt(12 * 0.2) + 1
  out.min <- 2 - out.max
  
  in.max <- 0.5 * sqrt(12 * 0.02) + 1
  in.min <- 2 - in.max
  
  #Generates quality vector with desired variance
  if(treatment == "+high"){
    quality <- runif(100, min=in.min, max=out.max)}
  
  if(treatment == "+low"){
    quality <- runif(100, min=out.min, max=in.max)}
  
  if(treatment == "full"){
    quality <- runif(100, min=out.min, max=out.max)}
  
  if(treatment == "low.var"){
    quality <- runif(100, min=in.min, max=in.max)
  }

  
  #Converts longevity (half-life in units of occupancy time) to decay rate in natural time
  long.nat <- longevity / parms$es
  r <- log(2) / long.nat
  
  # Runs simulation without infection until it reaches steady state
  Sonly <- diseaseSPOM(distance, quality, initial, parms, r, timesteps, 0.15)
  Sonly <- na.omit(Sonly)
  initial.inf <- Sonly[length(Sonly[, 1]), 2:101]
   
  #Randomly infects an occupied patch
  initial.inf[sample(which(initial.inf == "S"), 1)] <- "I"
  
  # Stores quality of initial infected patch
  initial.quality <- quality[which(initial.inf == "I")]
  
  # Calls SPOM to run with infection
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
  
  # Generating and collecting output
  output <- data.frame(treatment, longevity, S[length(S)], I[length(I)], E[length(E)],
                           mean(S[(length(S)-100):length(S)]), mean(I[(length(I)-100):length(I)]), mean(E[(length(I)-100):length(I)]),
                           max(I), output[which.max(I),1], initial.quality)
  
  names(output) <- c("treatment", "longevity", "S", "I", "R", "avgS", "avgI", "avgR", "maxI", "tmaxI", "quality0")
  
  return(output)
  
}