conn <- "(full)"
conn <- "(lattice)"
conn <- "(x0)"
conn <- "(delta)"


# Metapopulation-level results
#===============================================================

load(paste(getwd(), "/Data/out", conn, ".Rdata", sep = ""))

library(gridExtra)
library(lattice)
library(RColorBrewer)

metapop.S <- tapply(out$Sfin, list(out$longevity, out$variance), mean)
metapop.I <- tapply(out$Ifin, list(out$longevity, out$variance), mean)
metapop.occ <- tapply(out$Sfin + out$Ifin, list(out$longevity, out$variance), mean)
maxI <- tapply(out$maxI, list(out$longevity, out$variance), mean)

sd.S <- tapply(out$Sfin, list(out$longevity, out$variance), sd)
sd.I <- tapply(out$Ifin, list(out$longevity, out$variance), sd)
sd.maxI <- tapply(out$maxI, list(out$longevity, out$variance), sd)

metapop.pan <- tapply(out$Sfin == 0 & out$Ifin > 0, list(out$longevity, out$variance), sum) / 100
metapop.end <- tapply(out$Sfin > 0 & out$Ifin > 0, list(out$longevity, out$variance), sum) / 100
metapop.ext <- tapply(out$Sfin == 0 & out$Ifin == 0, list(out$longevity, out$variance), sum) / 100
metapop.nd <- tapply(out$Sfin > 0 & out$Ifin == 0, list(out$longevity, out$variance), sum) / 100

cols <- colorRampPalette(brewer.pal(9, "Reds"))(100)
p.pan <- levelplot(metapop.pan, col.regions = cols, xlab = "Longevity", ylab = "Variance", 
                   at = seq(0, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
p.end <- levelplot(metapop.end, col.regions = cols, xlab = "Longevity", ylab = "Variance", 
                   at = seq(0, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
p.ext <- levelplot(metapop.ext, col.regions = cols, xlab = "Longevity", ylab = "Variance",
                   at = seq(0, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
p.nd <- levelplot(metapop.nd, col.regions = cols, xlab = "Longevity", ylab = "Variance",
                  at = seq(0, 100, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))

key <- draw.colorkey(list(col = cols, 
                          at = seq(0, 100, by = 1),
                          labels = list(at = seq(0, 100, by = 25)),
                          height = 0.7), 
                     draw = F)

grid.arrange(p.end, p.pan, key, ncol = 3, widths = c(0.45, 0.45, 0.1))

cols <- colorRampPalette(brewer.pal(9, "Reds"))(100)

S <- levelplot(metapop.S, col.regions = cols, xlab = "Longevity", ylab = "Variance",
               at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
I <- levelplot(metapop.I, col.regions = cols, xlab = "Longevity", ylab = "Variance",
               at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
occ <- levelplot(metapop.occ, col.regions = cols, xlab = "Longevity", ylab = "Variance",
                 at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))

key <- draw.colorkey(list(col = cols, 
                          at = seq(0, 1, by = 0.01),
                          labels = list(at = seq(0, 1, by = 0.25)),
                          height = 0.7), 
                     draw = F)


grid.arrange(S, I, occ, key, ncol = 4, widths = c(0.3, 0.3, 0.3, 0.1))


#===============================================================================

# Plotting logistic regression surface for prob. of pathogen persistence
# (no longer in use -- needless level of abstraction from actual results)

# Generating pathogen persistence and pandemic indicator vectors
persists <- (metapop$I>0)
pandemic <- (metapop$I>0 & metapop$S==0)

logistic.per <- glm(formula = persists ~ metapop$var*metapop$longevity, family = binomial(link = "logit"))
logistic.pan <- glm(formula = pandemic[persists] ~ metapop$var[persists]*metapop$longevity[persists], 
                    family = binomial(link = "logit"))

p.per <- predict(logistic.per, type="response")
p.pan <- predict(logistic.pan, type="response")

plot.data <- data.frame(c(p.per, p.pan), rbind(metapop[c("var", "longevity")], metapop[persists, c("var", "longevity")]),
                        rep(c(1,2), times = c(length(persists), sum(persists))))
names(plot.data) <- c("p", "var", "longevity", "perpan")

trellis.par.set("axis.line", list(col="transparent"), las=0)
wireframe(p ~ var + longevity | perpan, data=plot.data,
          scales = list(arrows=FALSE, cex=1.2, col="black", tck=1.4, z = list(distance=1.2)),
          xlab = list(label="variance", cex=1.2, font=2, rot=30),
          ylab = list(label="longevity", cex=1.2, font=2, rot=320),
          zlab = list(label="Probability of outcome", cex=1.2, font=2, rot=93),
          layout = c(2, 1), strip=F,
          shade=F, col="gray50")

