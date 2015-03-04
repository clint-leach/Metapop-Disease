


library(gridExtra)
library(lattice)
library(RColorBrewer)

#===============================================================================
# Plots of totol population size

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

