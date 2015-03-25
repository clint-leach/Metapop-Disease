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
for(i in unique(patch$longevity)){
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

#===============================================================================
# Analyzing trap effect when both xi > 0

quality <- unique(out$quality)
longevity <- unique(out$longevity)

pinf <- tapply(out$I, list(out$xi_im, out$xi_em, out$longevity, out$quality), median)

par(mfrow = c(1, 2))
plot(out$quality, out$I, type = "n", ylim = c(0, 0.4), main = "xi_em = 0")
lines(supsmu(quality, pinf[1, 1, 10, ]), col = "black", pch = 20)
lines(supsmu(quality, pinf[2, 1, 10, ]), col = "red", pch = 20)

plot(out$quality, out$I, type = "n", ylim = c(0, 0.4), main = "xi_em = 0.5")
lines(supsmu(quality, pinf[1, 2, 10, ]), col = "black", pch = 20)
lines(supsmu(quality, pinf[2, 2, 10, ]), col = "red", pch = 20)

#===============================================================================
# Consequences of trap effect (with xi_em = 0.5)

trap <- out[out$xi_em == 0.5, ]

# Population size
Spop <- tapply(trap$S * trap$quality, list(trap$longevity, trap$xi_im), sum) / 100
Ipop <- tapply(trap$I * parms$nu *trap$quality, list(trap$longevity, trap$xi_im), sum) / 100
pop <- Spop + Ipop

par(mfrow = c(1, 1))
plot(log10(longevity), pop[, 1], type = "b", col = "black", pch = 20, ylim = c(0, 85))
lines(log10(longevity), pop[, 2], type = "b", col = "red", pch = 20)

# Pathogen persistence
ppersist <- tapply(trap$Ifin > 0, list(trap$longevity, trap$xi_im), sum) / 100 / 100

# Epidemiological outcome 
nd <- tapply(trap$Sfin > 0 & trap$Ifin == 0, list(trap$longevity, trap$xi_im), sum) / 100 / 100
end <- tapply(trap$Sfin > 0 & trap$Ifin > 0, list(trap$longevity, trap$xi_im), sum) / 100 / 100
pan <- tapply(trap$Sfin == 0 & trap$Ifin > 0, list(trap$longevity, trap$xi_im), sum) / 100 / 100
ext <- tapply(trap$Sfin == 0 & trap$Ifin == 0, list(trap$longevity, trap$xi_im), sum) / 100 / 100

par(mfrow = c(1, 4))
plot(log10(longevity), nd[, 1], type = "b", pch = 20, ylim = c(0, 1), main = "no disease")
lines(log10(longevity), nd[, 2], type = "b", col = "red", pch = 20)

plot(log10(longevity), end[,1], type = "b", pch = 20, ylim = c(0, 1), main = "endemic")
lines(log10(longevity), end[, 2], type = "b", col = "red", pch = 20)

plot(log10(longevity), pan[,1], type = "b", pch = 20, ylim = c(0, 1), main = "pandemic")
lines(log10(longevity), pan[, 2], type = "b", col = "red", pch = 20)

plot(log10(longevity), ext[,1], type = "b", pch = 20, ylim = c(0, 1), main = "extinction")
lines(log10(longevity), ext[, 2], type = "b", col = "red", pch = 20)


