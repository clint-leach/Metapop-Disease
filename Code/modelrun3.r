
patch.collect <- function(longevity, range, k, parms, distance, initial, timesteps){
  # Function that runs the SPOM model for fixed pathogen parameters and 
  # habitat quality distribution and collects patch-level simulation results.
  #
  # Takes in longevity and habitat quality range, runs the SPOM without infection
  # until steady state, then introduces infection on a randomly selected occupied patch,
  # calls SPOM model and generates output.
  #
  # Args:
  #   longevity: value for pathogen longevity in environment, half-life of infectivity
  #   range: vector containing the min and max of the quality distribution
  #   k: replicate ID
  #   parms: named dataframe of parameter values
  #   distance: nxn between patch distance matrix
  #   initial: character vector giving initial state of each patch ("S", "I", or "E")
  #   timesteps: number of timesteps to run SPOM for
  #
  # Returns:
  #   (number of patches) x 19 dataframe giving patch-level results from simulation
  #     quality: quality of patch
  #     tinf: time to first infection
  #     I: proportion of time patch infected over last 500 time steps
  #     S: proportion of time patch susceptible over last 500 time steps
  #     E: proportion of time patch empty over last 500 time steps
  #     inf.events: number of infection events
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
  #     xi_im: the value of xi_im used
  #     xi_em: the value of xi_em used
  
	# Generates quality vector with desired range
  max <- range[2]
  min <- range[1]
  quality <- seq(min, max, length.out = 100)
	
	# Converts longevity (half-life in units of occupancy time) to decay rate in natural time
  long.nat <- longevity / parms$es
	r <- log(2) / long.nat
	
  # Runs simulation without infection until it reaches steady state
	Sonly <- diseaseSPOM(distance, quality, initial, parms, r, timesteps, 0.15)
	Sonly <- na.omit(Sonly)
	initial.inf <- Sonly[length(Sonly[, 1]), 2:101]
  
  # Randomly infects an occupied patch
	initial.inf[sample(which(initial.inf=="S"), 1)] <- "I"
	
	# Stores quality of initial infected patch
	initial.quality <- quality[which(initial.inf == "I")]
	
	# Calls SPOM to run with infection
	output <- diseaseSPOM(distance, quality, initial.inf, parms, r, timesteps, 0.05)
	output <- na.omit(output)
  
  # Processing output
  runtime <- dim(output)[1]
  
  state <- output[(runtime - 499):runtime, -1]
  susc <- state == "S"
  inf <- state == "I"
  empty <- state == "E" | state == "R"
  
  S <- rowSums(susc) / 100
  I <- rowSums(inf) / 100
  E <- rowSums(empty) / 100
  
  inf.events <- susc[-(dim(state)[1]), ] * inf[-1, ]
  susc.col <- empty[-(dim(state)[1]), ] * susc[-1, ]
  inf.col <- empty[-(dim(state)[1]), ] * inf[-1, ]
  susc.ex <- susc[-(dim(state)[1]), ] * empty[-1, ]
  inf.ex <- inf[-(dim(state)[1]), ] * empty[-1, ]
  
  data <- data.frame(quality, 
                       tinf = apply(inf, 2, which.max),
                       I = colSums(inf) / (dim(state)[1]),
                       S = colSums(susc) / (dim(state)[1]),
                       E = colSums(empty) / (dim(state)[1]),
                       inf.events = colSums(inf.events),
                       susc.col = colSums(susc.col),
                       inf.col = colSums(inf.col),
                       susc.ex = colSums(susc.ex),
                       inf.ex = colSums(inf.ex),
                       Sfin = S[dim(state)[1]],
                       Ifin = I[dim(state)[1]],
                       Efin = E[dim(state)[1]],
                       maxI = max(I),
                       quality0 = initial.quality,
                       repID = k,
                       longevity = longevity,
                       xi_im = parms["xi_im"],
                       xi_em = parms["xi_em"]
  )
  

  return(data)
  
	}