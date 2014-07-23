
modelrun <- function(longevity, variance, k, parms, distance, initial, timesteps){
  # Takes in longevity and variance values, runs the SPOM without infection
  # until steady state, then introduces infection on a randomly selected occupied patch,
  # calls SPOM model and generates output.
  #
  # Args:
  #   longevity: value for pathogen longevity in environment, half-life of infectivity
  #   variance: variance of patch quality distribution (uniform with mean 1)
  #   k: replicate ID
  #   parms: named numeric vector of parameter values
  #   distance: nxn between patch distance matrix
  #   initial: character vector giving initial state of each patch
  #   timesteps: number of timesteps to run SPOM for
  #
  # Returns:
  #   n x 19 dataframe giving patch-level results from simulation
  #     quality: quality of patch
  #     tinf: time to first infection
  #     I: proportion of time patch infected
  #     S: proportion of time patch susceptible
  #     E: proportion of time patch empty
  #     inf.events.early: number of infection events in first fifth of sim
  #     inf.events.tot: total number of infection events
  #     susc.col: number of susceptible colonization events
  #     inf.col: number of infected colonization events
  #     susc.ex: number of susceptible extinction events
  #     inf.ex: number of infected extinction events
  #     Sfin: final susceptible occupancy
  #     Ifin: final infectious occupancy
  #     Efin: final empty occupancy
  #     maxI: maximum infectious occupancy
  #     quality0: quality of initially infected patch
  #     repID: replicate ID
  #     longevity: pathogen environemental longevity
  #     variance: variance of patch quality distribution
  
  
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
  
  inf.events <- susc[-(dim(output)[1]), ] * inf[-1, ]
  susc.col <- empty[-(dim(output)[1]), ] * susc[-1, ]
  inf.col <- empty[-(dim(output)[1]), ] * inf[-1, ]
  susc.ex <- susc[-(dim(output)[1]), ] * empty[-1, ]
  inf.ex <- inf[-(dim(output)[1]), ] * empty[-1, ]
  
  data <- data.frame(quality, 
                       tinf = apply(inf, 2, which.max),
                       I = colSums(inf) / (dim(output)[1]),
                       S = colSums(susc) / (dim(output)[1]),
                       E = colSums(empty) / (dim(output)[1]),
                       inf.events.early = colSums(inf.events[c(1:(timesteps / 5)), ]),
                       inf.events.tot = colSums(inf.events),
                       susc.col = colSums(susc.col),
                       inf.col = colSums(inf.col),
                       susc.ex = colSums(susc.ex),
                       inf.ex = colSums(inf.ex),
                       Sfin = S[dim(output)[1]],
                       Ifin = I[dim(output)[1]],
                       Efin = E[dim(output)[1]],
                       maxI = max(I),
                       quality0 = initial.quality,
                       repID = k,
                       longevity = longevity,
                       variance = variance
  )
  

  return(data)
  
	}