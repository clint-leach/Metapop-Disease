
library(gridExtra)
library(lattice)
library(RColorBrewer)
library(ggplot)

#===============================================================================
# Plots of difference in total population size between high and low quality habitat

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))

dat <- tapply(out$S.pop + out$I.pop, list(out$longevity, out$treatment, out$delta, out$nu), median)

cols <- colorRampPalette(brewer.pal(11, "RdBu"))(140)

low <- levelplot(dat[1, 2, , ] - dat[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
               at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med <- levelplot(dat[2, 2, , ] - dat[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high <- levelplot(dat[3, 2, , ] - dat[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))


load(paste(getwd(), "/Output/pathgrid(lattice).RData", sep = ""))

dat.lat <- tapply(out$S.pop + out$I.pop, list(out$longevity, out$treatment, out$delta, out$nu), median)

low.lattice <- levelplot(dat.lat[1, 2, , ] - dat.lat[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med.lattice <- levelplot(dat.lat[2, 2, , ] - dat.lat[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high.lattice <- levelplot(dat.lat[3, 2, , ] - dat.lat[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                  at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))


key <- draw.colorkey(list(col = cols, 
                          at = seq(-70, 70, by = 1),
                          labels = list(at = seq(-50, 50, by = 25)),
                          height = 0.7), 
                     draw = F)

grid.arrange(arrangeGrob(low, med, high, low.lattice, med.lattice, high.lattice, ncol = 3), 
             key, ncol = 2, widths = c(0.9, 0.1))

#===============================================================================
# Plots of susceptible and infected population sizes

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))

dat.S <- tapply(out$S.pop, list(out$longevity, out$treatment, out$delta, out$nu), median)
dat.I <- tapply(out$I.pop, list(out$longevity, out$treatment, out$delta, out$nu), median)

cols <- colorRampPalette(brewer.pal(9, "Greys"))(120)

low.h <- levelplot(dat.S[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med.h <- levelplot(dat.S[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high.h <- levelplot(dat.S[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                  at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))

low.l <- levelplot(dat.S[1, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med.l <- levelplot(dat.S[2, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high.l <- levelplot(dat.S[3, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                    at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))


key <- draw.colorkey(list(col = cols, 
                          at = seq(0, 120, by = 1),
                          labels = list(at = seq(0, 120, by = 20)),
                          height = 0.7), 
                     draw = F)

grid.arrange(arrangeGrob(low.h, med.h, high.h, low.l, med.l, high.l, ncol = 3), 
             key, ncol = 2, widths = c(0.9, 0.1))

###

cols <- colorRampPalette(brewer.pal(9, "Greys"))(120)

low.h <- levelplot(dat.I[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med.h <- levelplot(dat.I[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high.h <- levelplot(dat.I[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                    at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))

low.l <- levelplot(dat.I[1, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med.l <- levelplot(dat.I[2, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high.l <- levelplot(dat.I[3, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                    at = seq(0, 120, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))


key <- draw.colorkey(list(col = cols, 
                          at = seq(0, 120, by = 1),
                          labels = list(at = seq(0, 120, by = 20)),
                          height = 0.7), 
                     draw = F)

grid.arrange(arrangeGrob(low.h, med.h, high.h, low.l, med.l, high.l, ncol = 3), 
             key, ncol = 2, widths = c(0.9, 0.1))
#===============================================================================
# Plots of probability of host extinction

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))

dat <- tapply(out$S == 0 & out$I == 0, list(out$longevity, out$treatment, out$delta, out$nu), sum)

cols <- colorRampPalette(brewer.pal(11, "RdBu"))(200)

low <- levelplot(dat[1, 2, , ] - dat[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med <- levelplot(dat[2, 2, , ] - dat[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high <- levelplot(dat[3, 2, , ] - dat[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                  at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))


load(paste(getwd(), "/Output/pathgrid(lattice).RData", sep = ""))

dat.lat <- tapply(out$S == 0 & out$I == 0, list(out$longevity, out$treatment, out$delta, out$nu), sum)

low.lattice <- levelplot(dat.lat[1, 2, , ] - dat.lat[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                         at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med.lattice <- levelplot(dat.lat[2, 2, , ] - dat.lat[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                         at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high.lattice <- levelplot(dat.lat[3, 2, , ] - dat.lat[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                          at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))


key <- draw.colorkey(list(col = cols, 
                          at = seq(-100, 100, by = 1),
                          labels = list(at = seq(-50, 50, by = 25)),
                          height = 0.7), 
                     draw = F)

grid.arrange(arrangeGrob(low, med, high, low.lattice, med.lattice, high.lattice, ncol = 3), 
             key, ncol = 2, widths = c(0.9, 0.1))


#===============================================================================
# Plots of probability of pathogen persistence

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))

dat <- tapply(out$maxI > 0, list(out$longevity, out$treatment, out$delta, out$nu), sum)

cols <- colorRampPalette(brewer.pal(11, "RdBu"))(200)

low <- levelplot(dat[1, 2, , ] - dat[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med <- levelplot(dat[2, 2, , ] - dat[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high <- levelplot(dat[3, 2, , ] - dat[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                  at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))


load(paste(getwd(), "/Output/pathgrid(lattice).RData", sep = ""))

dat.lat <- tapply(out$maxI > 0, list(out$longevity, out$treatment, out$delta, out$nu), sum)

low.lattice <- levelplot(dat.lat[1, 2, , ] - dat.lat[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                         at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
med.lattice <- levelplot(dat.lat[2, 2, , ] - dat.lat[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                         at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
high.lattice <- levelplot(dat.lat[3, 2, , ] - dat.lat[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                          at = seq(-100, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))


key <- draw.colorkey(list(col = cols, 
                          at = seq(-100, 100, by = 1),
                          labels = list(at = seq(-50, 50, by = 25)),
                          height = 0.7), 
                     draw = F)

grid.arrange(arrangeGrob(low, med, high, low.lattice, med.lattice, high.lattice, ncol = 3), 
             key, ncol = 2, widths = c(0.9, 0.1))

#===============================================================================
# Plots of effect of longevity for a fixed delta and nu

# delta = 0.3
# nu = 0.2

sub <- out[abs(out$delta - 0.3) < 0.001 & out$nu == 0.2, ]

p.S <- ggplot(sub, aes(x = factor(log10(longevity)), y = S.pop)) + theme_classic()
p.S <- p.S + xlab("log(longevity)") + ylab("S population") + theme(legend.position = "none")
p.S <- p.S + geom_tufteboxplot(aes(colour = factor(treatment)), outlier.colour = "grey80", position = position_dodge(width = 0.4)) +
  scale_y_continuous(expand = c(0, 0.1))


p.I <- ggplot(sub, aes(x = factor(log10(longevity)), y = I.pop)) + theme_classic()
p.I <- p.I + xlab("log(longevity)") + ylab("I population") + theme(legend.position = "none") 
p.I <- p.I + geom_tufteboxplot(aes(colour = factor(treatment)), outlier.colour = "grey80", position = position_dodge(width = 0.4)) +
  scale_y_continuous(expand = c(0, 0.1))

p.tot <- ggplot(sub, aes(x = factor(log10(longevity)), y = I.pop + S.pop)) + theme_classic()
p.tot <- p.tot + xlab("log(longevity)") + ylab("Total population") + theme(legend.position = "none")
p.tot <- p.tot + geom_tufteboxplot(aes(colour = factor(treatment)), outlier.colour = "grey80", position = position_dodge(width = 0.4)) + 
  scale_y_continuous(expand = c(0, 0.1))

grid.arrange(p.S, p.I, p.tot, ncol = 3)
