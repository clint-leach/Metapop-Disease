# Patch level results -- scripts to analyze patch-level data from full simulation
# model.

# Select model parameterization to analyze
conn <- "(full)"
conn <- "(lattice)"
conn <- "(alpha0)"
conn <- "(delta)"

load(paste(getwd(), "/Output/out", conn, ".Rdata", sep = ""))

#Filter patch data to include only endemic cases
patch<-out[which(out$Ifin>0 & out$Sfin > 0), ]

#Filter patch data to include only high variance runs
patch<-patch[which(abs(patch$variance - 0.2) < 0.01), ]

#===============================================================================
#Plotting scatter cloud ofproportion of time occupied by susceptibles and infecteds as a function of quality

par(mfrow=c(1,2),cex=1)

smoothScatter(patch$quality,patch$S,ylab="Proportion of time occupied by susceptibles",xlab="Patch Quality")
lines(loess.smooth(patch$quality,patch$S),lwd=2,col="red")
smoothScatter(patch$quality,patch$I,ylab="Proportion of time occupied by infecteds",xlab="Patch Quality")
lines(loess.smooth(patch$quality,patch$I),lwd=2,col="red")


# Plotting scatter cloud of number of infection events
par(mfrow=c(1,1),cex=1)
smoothScatter(patch$quality,patch$inf.events,ylab="Number of infection events",xlab="Patch Quality")
lines(loess.smooth(patch$quality,patch$inf.events),lwd=3,col="red")

#===============================================================================
# Figures to show how longevity affects role of quality

# Plotting infection events against quality for range of longevities
plot(seq(0.25, 1.75, by = 0.1), c(1:16), ylim = c(0, 10), type = "n", xlab = "Quality", ylab = "Infections")
for(i in c(40, 60, 80, 100)){
  lines(loess.smooth(patch$quality[patch$longevity == i], patch$inf.events[patch$longevity == i]), type = "b", pch = as.character(i))  
}


#===============================================================================
# Showing simulation-level variability

# Number of infection events

library(scales)

reps <- unique(patch$repID)

par(mfrow = c(1, 1))
plot(patch$quality, patch$inf.events.early, type = "n", ylab = "Number of infection events", xlab = "Patch quality", bty = "l")
for(i in 1:length(reps)){
  sim <- patch[patch$repID == reps[i], ]
  try(lines(loess.smooth(sim$quality, sim$inf.events.tot), type = "l", col = alpha("gray", 0.3)))
}
lines(loess.smooth(patch$quality, patch$inf.events), lwd = 2, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 40], patch$inf.events[patch$longevity == 40]), lwd = 1, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 60], patch$inf.events[patch$longevity == 60]), lwd = 1, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 80], patch$inf.events[patch$longevity == 80]), lwd = 1, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 100], patch$inf.events[patch$longevity == 100]), lwd = 1, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 120], patch$inf.events[patch$longevity == 120]), lwd = 1, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 140], patch$inf.events[patch$longevity == 140]), lwd = 1, col = "black")


par(mfrow = c(2, 2))

# Number of infected extinction events
plot(patch$quality, patch$inf.ex, type = "n", ylab = "Number of infected extinctions", xlab = "Patch quality")
for(i in 1:length(reps)){
  sim <- patch[patch$repID == reps[i], ]
  lines(loess.smooth(sim$quality, sim$inf.ex), type = "l", col = alpha("gray", 0.3))
}
lines(loess.smooth(patch$quality, patch$inf.ex), lwd = 2, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 40], patch$inf.ex[patch$longevity == 40]), lwd = 1.5, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 100], patch$inf.ex[patch$longevity == 100]), lwd = 1.5, col = "black")


# Number of susceptible extinction events
plot(patch$quality, patch$susc.ex, type = "n", ylab = "Number of susceptible extinctions", xlab = "Patch quality")
for(i in 1:length(reps)){
  sim <- patch[patch$repID == reps[i], ]
  lines(loess.smooth(sim$quality, sim$susc.ex), type = "l", col = alpha("gray", 0.3))
}
lines(loess.smooth(patch$quality, patch$susc.ex), lwd = 2, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 40], patch$susc.ex[patch$longevity == 40]), lwd = 1.5, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 100], patch$susc.ex[patch$longevity == 100]), lwd = 1.5, col = "black")


# Infected colonizations
plot(patch$quality, patch$inf.col, type = "n", ylab = "Number of infected colonizations", xlab = "Patch quality")
for(i in 1:length(reps)){
  sim <- patch[patch$repID == reps[i], ]
  lines(loess.smooth(sim$quality, sim$inf.col), type = "l", col = alpha("gray", 0.3))
}
lines(loess.smooth(patch$quality, patch$inf.col), lwd = 2, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 40], patch$inf.col[patch$longevity == 40]), lwd = 1.5, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 100], patch$inf.col[patch$longevity == 100]), lwd = 1.5, col = "black")


# Susceptible colonizations
plot(patch$quality, patch$inf.col, type = "n", ylab = "Number of susceptible colonizations", xlab = "Patch quality")
for(i in 1:length(reps)){
  sim <- patch[patch$repID == reps[i], ]
  lines(loess.smooth(sim$quality, sim$susc.col), type = "l", col = alpha("gray", 0.3))
}
lines(loess.smooth(patch$quality, patch$susc.col), lwd = 2, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 40], patch$susc.col[patch$longevity == 40]), lwd = 1.5, col = "black")
lines(loess.smooth(patch$quality[patch$longevity == 100], patch$susc.col[patch$longevity == 100]), lwd = 1.5, col = "black")
