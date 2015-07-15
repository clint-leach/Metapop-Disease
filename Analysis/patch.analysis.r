# Patch level results -- scripts to analyze patch-level data 

load(paste(getwd(), "/Output/trap.RData", sep = ""))

quality <- unique(out$quality)
longevity <- unique(out$longevity)

# Plotting probability of infection (i.e. time spent infected) as function of
# the xi's, longevity, and quality
pinf <- tapply(out$I, list(out$xi_im, out$xi_em, out$longevity, out$quality), median)

# Plots for highest longevity
par(mfrow = c(1, 2))
plot(out$quality, out$I, type = "n", ylim = c(0, 0.4), main = "xi_em = 0")
lines(supsmu(quality, pinf[1, 1, 10, ]), col = "black", pch = 20)
lines(supsmu(quality, pinf[2, 1, 10, ]), col = "red", pch = 20)

plot(out$quality, out$I, type = "n", ylim = c(0, 0.4), main = "xi_em = 0.5")
lines(supsmu(quality, pinf[1, 2, 10, ]), col = "black", pch = 20)
lines(supsmu(quality, pinf[2, 2, 10, ]), col = "red", pch = 20)

# Scatter plots for all longevities
par(mfrow = c(2, 5))
for(i in 1:length(longevity)){
  plot(out$quality, out$I, type = "n", ylim = c(0, 0.4), main = longevity[i])
  points(quality, pinf[1, 2, i, ], col = "black", pch = 20)
  points(quality, pinf[2, 2, i, ], col = "red", pch = 20)
}


###

# Boxplots for all longevities and xi_im = 0.5
par(mfrow = c(2, 5))
for(i in longevity){
  sub <- subset(out, xi_em == 0.5 & xi_im == 0.5 & longevity == i)
  boxplot(I ~ quality, data = sub)
}

# Boxplots for highest longevity and  both values of xi_im 
par(mfrow = c(1, 2))
boxplot(I ~ quality, data = out, subset = xi_em == 0.5 & xi_im == 0 & longevity > 3)
boxplot(I ~ quality, data = out, subset = xi_em == 0.5 & xi_im == 0.5 & longevity > 3)

###

# Tufte-style boxplot for highest longevity and both xi_im side-by-side
sub <- subset(out, xi_em == 0.5 & longevity > 3)

p <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p <- p + xlab("Quality") + ylab("Probability of infection") 
p <- p + geom_tufteboxplot(aes(colour = factor(xi_im)), outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_colour_manual(values = c("black", "red"))
p

###

# Multipanel figure comparing xi_im for three longevities
longs <- longevity[c(4, 7, 10)]

labels <- vector(length = 100, mode = "character")
labels[c(1, 50, 100)] <- c("0.2", "1", "1.8")

sub <- subset(out, xi_em == 0.5 & xi_im == 0 & longevity == longs[1])
p1 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p1 <- p1 + xlab("Quality") + ylab("Probability of infection") 
p1 <- p1 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels)

sub <- subset(out, xi_em == 0.5 & xi_im == 0.5 & longevity == longs[1])
p2 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p2 <- p2 + xlab("Quality") + ylab("Probability of infection") 
p2 <- p2 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels)


sub <- subset(out, xi_em == 0.5 & xi_im == 0 & longevity == longs[2])
p3 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p3 <- p3 + xlab("Quality") + ylab("Probability of infection") 
p3 <- p3 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels)

sub <- subset(out, xi_em == 0.5 & xi_im == 0.5 & longevity == longs[2])
p4 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p4 <- p4 + xlab("Quality") + ylab("Probability of infection") 
p4 <- p4 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels)


sub <- subset(out, xi_em == 0.5 & xi_im == 0 & longevity == longs[3])
p5 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p5 <- p5 + xlab("Quality") + ylab("Probability of infection") 
p5 <- p5 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels)

sub <- subset(out, xi_em == 0.5 & xi_im == 0.5 & longevity == longs[3])
p6 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p6 <- p6 + xlab("Quality") + ylab("Probability of infection") 
p6 <- p6 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels)

grid.arrange(p1, p3, p5, p2, p4, p6, nrow = 2)

#===============================================================================
# Multipanel figure comparing extinctions across longevity and xi_im

ex <- tapply(out$inf.ex + out$susc.ex, list(out$xi_em, out$xi_im, out$longevity, out$quality), mean)

par(mfrow = c(2, 3))

plot(quality, ex[2, 1, 4, ], ylim = c(0, 7), ylab = "Extinctions", xlab = "Quality", pch = 20)
plot(quality, ex[2, 1, 7, ], ylim = c(0, 7), ylab = "Extinctions", xlab = "Quality", pch = 20)
plot(quality, ex[2, 1, 10, ], ylim = c(0, 7), ylab = "Extinctions", xlab = "Quality", pch = 20)
lines(supsmu(quality, ex[2, 1, 4, ]))

plot(quality, ex[2, 2, 4, ], ylim = c(0, 7), ylab = "Extinctions", xlab = "Quality", pch = 20)
plot(quality, ex[2, 2, 7, ], ylim = c(0, 7), ylab = "Extinctions", xlab = "Quality", pch = 20)
plot(quality, ex[2, 2, 10, ], ylim = c(0, 7), ylab = "Extinctions", xlab = "Quality", pch = 20)
lines(supsmu(quality, ex[2, 2, 4, ]))

#===============================================================================
# Multipanel figure comparing S and I across longevity and xi_im

S <- tapply(out$S, list(out$xi_em, out$xi_im, out$longevity, out$quality), median)

par(mfrow = c(2, 3))

plot(quality, S[2, 1, 4, ], ylim = c(0, 1), ylab = "Probability susceptible", xlab = "Quality", pch = 4)
points(quality, S[2, 2, 4, ], col = "brown", pch = 4)
plot(quality, S[2, 1, 7, ], ylim = c(0, 1), ylab = "Probability susceptible", xlab = "Quality", pch = 4)
points(quality, S[2, 2, 7, ], col = "brown", pch = 4)
plot(quality, S[2, 1, 10, ], ylim = c(0, 1), ylab = "Probability susceptible", xlab = "Quality", pch = 4)
points(quality, S[2, 2, 10, ], col = "brown", pch = 4)

I <- tapply(out$I, list(out$xi_em, out$xi_im, out$longevity, out$quality), median)

plot(quality, I[2, 1, 4, ], ylim = c(0, 1), ylab = "Probability infectious", xlab = "Quality", pch = 4)
points(quality, I[2, 2, 4, ], col = "brown", pch = 4)
plot(quality, I[2, 1, 7, ], ylim = c(0, 1), ylab = "Probability infectious", xlab = "Quality", pch = 4)
points(quality, I[2, 2, 7, ], col = "brown", pch = 4)
plot(quality, I[2, 1, 10, ], ylim = c(0, 1), ylab = "Probability infectious", xlab = "Quality", pch = 4)
points(quality, I[2, 2, 10, ], col = "brown", pch = 4)

par(mfrow = c(1, 3))

occ <- tapply(out$I + out$S, list(out$xi_em, out$xi_im, out$longevity, out$quality), median)

plot(quality, occ[2, 1, 4, ], ylim = c(0, 1), ylab = "Probability occupied", xlab = "Quality", pch = 4)
points(quality, occ[2, 2, 4, ], col = "brown", pch = 4)
plot(quality, occ[2, 1, 7, ], ylim = c(0, 1), ylab = "Probability occupied", xlab = "Quality", pch = 4)
points(quality, occ[2, 2, 7, ], col = "brown", pch = 4)
plot(quality, occ[2, 1, 10, ], ylim = c(0, 1), ylab = "Probability occupied", xlab = "Quality", pch = 4)
points(quality, occ[2, 2, 10, ], col = "brown", pch = 4)
#===============================================================================
# Consequences of preference (with xi_em = 0.5)

trap <- out[out$xi_em == 0.5, ]

# Population size
sims <- tapply((trap$S + trap$I * parms$nu) *trap$quality, list(trap$longevity, trap$xi_im), matrix, nrow = 100, byrow = T)
pop <- lapply(sims, rowSums)
pop <- lapply(pop, median)
pop <- matrix(unlist(pop), nrow = 10, byrow = F)

par(mfrow = c(1, 1))
plot(log10(longevity), pop[, 1], type = "b", col = "black", pch = 20, ylim = c(0, 85))
lines(log10(longevity), pop[, 2], type = "b", col = "red", pch = 20)

# Pathogen persistence
ppersist <- tapply(trap$maxI > 0, list(trap$longevity, trap$xi_im), sum) / 100 / 100

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


###

# Stand-alone plot of extinction probability

par(mfrow = c(1, 1))
plot(log10(longevity[c(4, 7, 10)]), ext[c(4, 7, 10), 1], type = "b", pch = 20, ylim = c(0, 0.4),
     xlab = "log(longevity)", ylab = "Probability of extinction", bty = "l")
lines(log10(longevity[c(4, 7, 10)]), ext[c(4, 7, 10), 2], type = "b", pch = 20, col = "red")

