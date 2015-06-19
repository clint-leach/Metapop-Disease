# Script to generate all plots for manuscript

library(gridExtra)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(plyr)

#===============================================================================
### Figure 1

pdf("Manuscript/figure/figure_1.pdf", width = 10, height = 6)

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))

out$endemic <- out$S > 0 & out$I > 0
dat <- ddply(out, c("longevity", "treatment", "delta", "nu"), summarise, sum = sum(endemic) / 100)

cols <- colorRampPalette(brewer.pal(9, "Greys"))(100)
labels <- c("d", "e", "f", "a", "b", "c")

plot <- levelplot(sum ~ delta * nu | longevity * treatment, data = dat,
          col.regions = cols,
          at = seq(0, 1, by = 0.01),
          xlab = list(label = expression("Probability of direct transmission"~(delta)), cex = 1),
          ylab = list(label = expression("Infectious survival"~(nu)), cex = 1),
          strip = F,
          scales = list(alternating = F, cex = 1),
          ylab.right = list(label = "Probability of endemic", cex = 1),
          par.settings = list(layout.widths = list(axis.key.padding = 0, ylab.right = 3)),
          panel=function(...){
            panel.levelplot(...)
            panel.text(0, 1, labels[panel.number()], cex = 1)
          })

col.lab <- textGrob(expression(Longevity ~ symbol("\256")))
row.low <- textGrob("Low quality", rot = 90)
row.high <- textGrob("High quality", rot = 90)

grid.arrange(arrangeGrob(row.low, row.high, ncol = 1), plot, ncol = 2, 
             widths = c(0.05, 0.95), sub = col.lab)

dev.off()

#===============================================================================
### Figure 2

rm(list = ls())

pdf("Manuscript/figure/figure_2.pdf", width = 10, height = 4)

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))

out$pop<- out$S.pop + out$I.pop

high <- ddply(subset(out, treatment == "+high"), c("treatment", "longevity", "delta", "nu"), 
              summarise, high = median(pop))
low <- ddply(subset(out, treatment == "+low"), c("treatment", "longevity", "delta", "nu"), 
              summarise, low = median(pop))
joint <- join(high, low, by = c("longevity", "delta", "nu"))
joint$diff <- joint$low - joint$high

cols <- colorRampPalette(brewer.pal(11, "RdBu"))(140)

labels = c("a", "b", "c")

plot <- levelplot(diff ~ delta * nu | longevity, data = joint,
                  col.regions = cols,
                  at = seq(-70, 70, by = 1),
                  xlab = list(label = expression("Probability of direct transmission"~(delta)), cex = 1),
                  ylab = list(label = expression("Infectious survival"~(nu)), cex = 1),
                  strip = F,
                  scales = list(alternating = F, cex = 1),
                  ylab.right = list(label = "Pop on low quality - pop on high quality", cex = 1),
                  par.settings = list(layout.widths = list(axis.key.padding = 0, ylab.right = 3)),
                  panel=function(...){
                    panel.levelplot(main = labels[panel.number()], ...)
                    panel.text(0, 1, labels[panel.number()], cex = 1.5, col = "white")
                  })

col.lab <- textGrob(expression(Longevity ~ symbol("\256")))
grid.arrange(plot, sub = col.lab)

dev.off()

#===============================================================================
### Figure 3

rm(list = ls())

pdf("Manuscript/figure/figure_3.pdf", width = 10, height = 5)

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))
sub <- out[abs(out$delta - 0.3) < 0.001 & abs(out$nu - 0.2) < 0.001, ]

S <- ddply(sub, c("longevity", "treatment"), summarise, 
           med = median(S.pop), 
           lowq = quantile(S.pop, 0.25),
           highq = quantile(S.pop, 0.75),
           min = min(S.pop),
           max = max(S.pop))

I <- ddply(sub, c("longevity", "treatment"), summarise, 
           med = median(I.pop), 
           lowq = quantile(I.pop, 0.25),
           highq = quantile(I.pop, 0.75),
           min = min(I.pop),
           max = max(I.pop))

tot <- ddply(sub, c("longevity", "treatment"), summarise, 
           med = median(S.pop + I.pop), 
           lowq = quantile(S.pop + I.pop, 0.25),
           highq = quantile(S.pop + I.pop, 0.75),
           min = min(S.pop + I.pop),
           max = max(S.pop + I.pop))

pS <- ggplot(S, aes(x = factor(log10(longevity)), y = med))
pS <- pS + geom_point(aes(colour = factor(treatment)), position = position_dodge(width = 0.4)) +
             geom_point(aes(colour = factor(treatment), y = min), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
             geom_point(aes(colour = factor(treatment), y = max), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
             geom_errorbar(aes(colour = factor(treatment), ymin = lowq, ymax = highq), 
                           position = position_dodge(width = 0.4), width = 0)
pS <- pS + scale_y_continuous(limits = c(0, 125), expand = c(0, 0.1)) +
             ylab("S population") + xlab("log(longevity)") + ggtitle("a") +
             theme_classic() + theme(legend.position = "none")


pI <- ggplot(I, aes(x = factor(log10(longevity)), y = med))
pI <- pI + geom_point(aes(colour = factor(treatment)), position = position_dodge(width = 0.4)) +
           geom_point(aes(colour = factor(treatment), y = min), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
           geom_point(aes(colour = factor(treatment), y = max), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
           geom_errorbar(aes(colour = factor(treatment), ymin = lowq, ymax = highq), 
                         position = position_dodge(width = 0.4), width = 0)
pI <- pI + scale_y_continuous(limits = c(0, 12), expand = c(0, 0.1)) +
           ylab("I population") + xlab("log(longevity)") + ggtitle("b") +
           theme_classic() + theme(legend.position = "none")


ptot <- ggplot(tot, aes(x = factor(log10(longevity)), y = med))
ptot <- ptot + geom_point(aes(colour = factor(treatment)), position = position_dodge(width = 0.4)) +
               geom_point(aes(colour = factor(treatment), y = min), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
               geom_point(aes(colour = factor(treatment), y = max), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
               geom_errorbar(aes(colour = factor(treatment), ymin = lowq, ymax = highq), 
                             position = position_dodge(width = 0.4), width = 0)
ptot <- ptot + scale_y_continuous(limits = c(0, 125), expand = c(0, 0.1)) +
               ylab("Total population") + xlab("log(longevity)") + ggtitle("c") +
               theme_classic() + theme(legend.position = "none")

grid.arrange(pS, pI, ptot, ncol = 3)

dev.off()

#===============================================================================
### Figure 4

# See code in simvis.r

#===============================================================================
### Figure 5

rm(list = ls())

pdf("Manuscript/figure/figure_5.pdf", width = 12, height = 6)

load(paste(getwd(), "/Output/trap.RData", sep = ""))

quality <- unique(out$quality)
longevity <- unique(out$longevity)

longs <- longevity[c(4, 7, 10)]

labels <- vector(length = 100, mode = "character")
labels[c(1, 50, 100)] <- c("0.2", "1", "1.8")

sub <- subset(out, xi_em == 0.5 & xi_im == 0 & longevity == longs[1])
p1 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p1 <- p1 + xlab("Quality") + ylab("Probability of infection") + ggtitle("a.")
p1 <- p1 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels) + theme(axis.ticks.x = element_blank())

sub <- subset(out, xi_em == 0.5 & xi_im == 0.5 & longevity == longs[1])
p2 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p2 <- p2 + xlab("Quality") + ylab("Probability of infection") + ggtitle("d.")
p2 <- p2 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels) + theme(axis.ticks.x = element_blank())


sub <- subset(out, xi_em == 0.5 & xi_im == 0 & longevity == longs[2])
p3 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p3 <- p3 + xlab("Quality") + ylab("Probability of infection") + ggtitle("b.") 
p3 <- p3 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels) + theme(axis.ticks.x = element_blank())

sub <- subset(out, xi_em == 0.5 & xi_im == 0.5 & longevity == longs[2])
p4 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p4 <- p4 + xlab("Quality") + ylab("Probability of infection") + ggtitle("e.")
p4 <- p4 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels) + theme(axis.ticks.x = element_blank())


sub <- subset(out, xi_em == 0.5 & xi_im == 0 & longevity == longs[3])
p5 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p5 <- p5 + xlab("Quality") + ylab("Probability of infection") + ggtitle("c.")
p5 <- p5 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels) + theme(axis.ticks.x = element_blank())

sub <- subset(out, xi_em == 0.5 & xi_im == 0.5 & longevity == longs[3])
p6 <- ggplot(sub, aes(x = as.factor(quality), y = I)) + theme_classic()
p6 <- p6 + xlab("Quality") + ylab("Probability of infection") + ggtitle("f.")
p6 <- p6 + geom_tufteboxplot(outlier.colour = "grey80", position = position_dodge(width = 1)) +
  scale_y_continuous(expand = c(0, 0.01)) + scale_x_discrete(labels = labels) + theme(axis.ticks.x = element_blank())

grid.arrange(p1, p3, p5, p2, p4, p6, nrow = 2)

dev.off()
  
#===============================================================================
### Figure 6

rm(list = ls())

pdf("Manuscript/figure/figure_6.pdf", width = 5, height = 5)

load(paste(getwd(), "/Output/trap.RData", sep = ""))

quality <- unique(out$quality)
longevity <- unique(out$longevity)

trap <- out[out$xi_em == 0.5, ]
ext <- tapply(trap$Sfin == 0 & trap$Ifin == 0, list(trap$longevity, trap$xi_im), sum) / 100 / 100

par(mfrow = c(1, 1))
plot(log10(longevity[c(4, 7, 10)]), ext[c(4, 7, 10), 1], type = "b", pch = 20, ylim = c(0, 0.4),
     xlab = "log(longevity)", ylab = "Probability of extinction", bty = "l")
lines(log10(longevity[c(4, 7, 10)]), ext[c(4, 7, 10), 2], type = "b", pch = 20, col = "red")

dev.off()

#===============================================================================
### Supplemental figure 1

rm(list = ls())

pdf("Manuscript/figure/supplement_1.pdf", width = 10, height = 6)

load(paste(getwd(), "/Output/pathgrid(lattice).RData", sep = ""))

dat <- tapply(out$I > 0 & out$S > 0, list(out$longevity, out$treatment, out$delta, out$nu), sum) / 100

cols <- colorRampPalette(brewer.pal(9, "Greys"))(100)

low.l <- levelplot(dat[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "a.")
med.l <- levelplot(dat[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "b.")
high.l <- levelplot(dat[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                    at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "c.")

low.h <- levelplot(dat[1, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "d.")
med.h <- levelplot(dat[2, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                   at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "e.")
high.h <- levelplot(dat[3, 2, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                    at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "f.")

key <- draw.colorkey(list(col = cols, 
                          at = seq(0, 1, by = 0.01),
                          labels = list(at = seq(0, 1, by = 0.25)),
                          height = 0.7), 
                     draw = F)

grid.arrange(arrangeGrob(low.l, med.l, high.l, low.h, med.h, high.h, ncol = 3), 
             key, ncol = 2, widths = c(0.9, 0.1))

dev.off()

#===============================================================================
### Supplemental figure 2

rm(list = ls())

pdf("Manuscript/figure/supplement_2.pdf", width = 10, height = 5)

load(paste(getwd(), "/Output/pathgrid(lattice).RData", sep = ""))

dat <- tapply(out$S.pop + out$I.pop, list(out$longevity, out$treatment, out$delta, out$nu), median)

cols <- colorRampPalette(brewer.pal(11, "RdBu"))(140)

low <- levelplot(dat[1, 2, , ] - dat[1, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "a.")
med <- levelplot(dat[2, 2, , ] - dat[2, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                 at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "b.")
high <- levelplot(dat[3, 2, , ] - dat[3, 1, , ], col.regions = cols, xlab = expression(delta), ylab = expression(nu),
                  at = seq(-70, 70, by = 1), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)), main = "c.")

key <- draw.colorkey(list(col = cols, 
                          at = seq(-70, 70, by = 1),
                          labels = list(at = seq(-50, 50, by = 25)),
                          height = 0.7), 
                     draw = F)

grid.arrange(low, med, high, key, ncol = 4, widths = c(0.3, 0.3, 0.3, 0.1))

dev.off()

#===============================================================================
### Supplemental figure 3

rm(list = ls())

pdf("Manuscript/figure/supplement_3.pdf", width = 10, height = 10)

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))

dat <- tapply(out$S.pop + out$I.pop, list(out$longevity, out$treatment, out$delta, out$nu), median)

cols <- colorRampPalette(brewer.pal(11, "RdBu"))(140)

diff <- dat[, 2, , ] - dat[, 1, , ]
signdiff <- sign(diff)
absdiff <- abs(diff)

par(mfcol = c(10, 10), mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
for(delta in 1:10){
  for(nu in 10:1){
    plot(log10(longevity), diff[, delta, nu], type = "n", axes = F, xlab = "", ylab = "")
    gradient.rect(-10, -70, 10, 70, col = cols, gradient = "y")
    lines(log10(longevity), diff[, delta, nu], lwd = 2)
  }
}
mtext(expression(nu), side = 2, outer = T, line = 2.5, cex = 1.2)
mtext(expression(delta), side = 1, outer = T, line = 2.5, cex = 1.2)

xlabs = as.character(seq(0, 0.9, by = 0.1))
mtext(xlabs, side = 1, outer = T, line = 0.5, at = seq(0.05, 0.95, by = 0.1))

ylabs = as.character(seq(0.1, 1, by = 0.1))
mtext(ylabs, side = 2, outer = T, line = 0.5, at = seq(0.05, 0.95, by = 0.1), las = 1)

dev.off()
