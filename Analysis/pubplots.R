# Script to generate all plots for manuscript

library(gridExtra)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(plyr)

#===============================================================================
### Figure 1

postscript("Manuscript/figure/figure_1.eps", width = 10, height = 6)

load(paste(getwd(), "/Output/pathgrid(lattice).RData", sep = ""))

dat <- ddply(out, c("longevity", "treatment", "delta", "nu"), summarise, susc = median(S.pop))

cols <- colorRampPalette(brewer.pal(9, "Greys"))(120)
labels <- c("d", "e", "f", "a", "b", "c")

plot <- levelplot(susc ~ delta * nu | longevity * treatment, data = dat,
          col.regions = cols,
          at = seq(0, 120, by = 1),
          xlab = list(label = expression("Probability of direct transmission"~(delta)), cex = 1),
          ylab = list(label = expression("Infectious survival"~(nu)), cex = 1),
          strip = F,
          scales = list(alternating = F, cex = 1),
          ylab.right = list(label = "Susceptible population size", cex = 1),
          par.settings = list(layout.widths = list(axis.key.padding = 0, ylab.right = 3)),
          panel=function(...){
            panel.levelplot(...)
            panel.text(0, 1, labels[panel.number()], cex = 1)
          })

row.low <- textGrob("Low quality", rot = 90)
row.high <- textGrob("High quality", rot = 90)

grid.arrange(arrangeGrob(row.low, row.high, ncol = 1), plot, ncol = 2, 
             widths = c(0.05, 0.95))

grid.text("Longevity = 0.32", x = unit(0.26, units = "npc"), y = unit(0.96, units = "npc"), draw = T)
grid.text("Longevity = 1.0", x = unit(0.50, units = "npc"), y = unit(0.96, units = "npc"), draw = T)
grid.text("Longevity = 3.2", x = unit(0.75, units = "npc"), y = unit(0.96, units = "npc"), draw = T)

dev.off()

#===============================================================================
### Figure 2

rm(list = ls())

postscript("Manuscript/figure/figure_2.eps", width = 10, height = 4)

load(paste(getwd(), "/Output/pathgrid(lattice).RData", sep = ""))

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
                    panel.levelplot(...)
                    panel.text(0, 1, labels[panel.number()], cex = 1.5, col = "white")
                  })

plot

grid.text("Longevity = 0.32", x = unit(0.22, units = "npc"), y = unit(0.96, units = "npc"), draw = T)
grid.text("Longevity = 1.0", x = unit(0.48, units = "npc"), y = unit(0.96, units = "npc"), draw = T)
grid.text("Longevity = 3.2", x = unit(0.74, units = "npc"), y = unit(0.96, units = "npc"), draw = T)

dev.off()

#===============================================================================
### Figure 3

rm(list = ls())

postscript("Manuscript/figure/figure_3.eps", width = 10, height = 5)

load(paste(getwd(), "/Output/pathgrid(lattice).RData", sep = ""))
sub <- out[abs(out$delta - 0.5) < 0.001 & abs(out$nu - 0.2) < 0.001, ]
sub$longevity <- round(sub$longevity, 1)

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

pS <- ggplot(S, aes(x = factor(longevity), y = med))
pS <- pS + geom_point(aes(colour = factor(treatment)), position = position_dodge(width = 0.4)) +
             geom_point(aes(colour = factor(treatment), y = min), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
             geom_point(aes(colour = factor(treatment), y = max), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
             geom_errorbar(aes(colour = factor(treatment), ymin = lowq, ymax = highq), 
                           position = position_dodge(width = 0.4), width = 0)
pS <- pS + scale_y_continuous(limits = c(0, 125), expand = c(0, 0.1)) +
             ylab("S population") + xlab("longevity") + ggtitle("a") +
             theme_classic() + theme(legend.position = "none")


pI <- ggplot(I, aes(x = factor(longevity), y = med))
pI <- pI + geom_point(aes(colour = factor(treatment)), position = position_dodge(width = 0.4)) +
           geom_point(aes(colour = factor(treatment), y = min), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
           geom_point(aes(colour = factor(treatment), y = max), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
           geom_errorbar(aes(colour = factor(treatment), ymin = lowq, ymax = highq), 
                         position = position_dodge(width = 0.4), width = 0)
pI <- pI + scale_y_continuous(limits = c(0, 12), expand = c(0, 0.1)) +
           ylab("I population") + xlab("longevity") + ggtitle("b") +
           theme_classic() + theme(legend.position = "none")


ptot <- ggplot(tot, aes(x = factor(longevity), y = med))
ptot <- ptot + geom_point(aes(colour = factor(treatment)), position = position_dodge(width = 0.4)) +
               geom_point(aes(colour = factor(treatment), y = min), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
               geom_point(aes(colour = factor(treatment), y = max), shape = 45, size = 5, position = position_dodge(width = 0.4)) +
               geom_errorbar(aes(colour = factor(treatment), ymin = lowq, ymax = highq), 
                             position = position_dodge(width = 0.4), width = 0)
ptot <- ptot + scale_y_continuous(limits = c(0, 125), expand = c(0, 0.1)) +
               ylab("Total population") + xlab("longevity") + ggtitle("c") +
               theme_classic() + theme(legend.position = "none")

grid.arrange(pS, pI, ptot, ncol = 3)

dev.off()

#===============================================================================
### Figure 4

# See code in simvis.r

#===============================================================================
### Supplemental figure 1

rm(list = ls())

postscript("Manuscript/figure/supplement_1.eps", width = 10, height = 6)

load(paste(getwd(), "/Output/pathgrid.RData", sep = ""))

dat <- ddply(out, c("longevity", "treatment", "delta", "nu"), summarise, susc = median(S.pop))

cols <- colorRampPalette(brewer.pal(9, "Greys"))(120)
labels <- c("d", "e", "f", "a", "b", "c")

plot <- levelplot(susc ~ delta * nu | longevity * treatment, data = dat,
                  col.regions = cols,
                  at = seq(0, 120, by = 1),
                  xlab = list(label = expression("Probability of direct transmission"~(delta)), cex = 1),
                  ylab = list(label = expression("Infectious survival"~(nu)), cex = 1),
                  strip = F,
                  scales = list(alternating = F, cex = 1),
                  ylab.right = list(label = "Susceptible population size", cex = 1),
                  par.settings = list(layout.widths = list(axis.key.padding = 0, ylab.right = 3)),
                  panel=function(...){
                    panel.levelplot(...)
                    panel.text(0, 1, labels[panel.number()], cex = 1)
                  })

row.low <- textGrob("Low quality", rot = 90)
row.high <- textGrob("High quality", rot = 90)

grid.arrange(arrangeGrob(row.low, row.high, ncol = 1), plot, ncol = 2, 
             widths = c(0.05, 0.95))

grid.text("Longevity = 0.32", x = unit(0.26, units = "npc"), y = unit(0.96, units = "npc"), draw = T)
grid.text("Longevity = 1.0", x = unit(0.50, units = "npc"), y = unit(0.96, units = "npc"), draw = T)
grid.text("Longevity = 3.2", x = unit(0.75, units = "npc"), y = unit(0.96, units = "npc"), draw = T)


dev.off()

#===============================================================================
### Supplemental figure 2

rm(list = ls())

postscript("Manuscript/figure/supplement_2.eps", width = 10, height = 4)

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


plot

grid.text("Longevity = 0.32", x = unit(0.22, units = "npc"), y = unit(0.96, units = "npc"), draw = T)
grid.text("Longevity = 1.0", x = unit(0.48, units = "npc"), y = unit(0.96, units = "npc"), draw = T)
grid.text("Longevity = 3.2", x = unit(0.74, units = "npc"), y = unit(0.96, units = "npc"), draw = T)

dev.off()

#===============================================================================
### Supplemental figure 3

rm(list = ls())

postscript("Manuscript/figure/supplement_3.eps", width = 8, height = 4)

load(paste(getwd(), "/Output/trap(lattice).RData", sep = ""))

longs <- unique(out$longevity)[c(4, 7, 10)]
labels <- c("a", "b", "c")

dat <- subset(out, xi_em == 0.5 & longevity %in% longs)

S <- ddply(dat, .(xi_im, longevity, quality), summarise, occ = median(S))

xyplot(occ ~ quality | longevity, groups = xi_im, data = S,
        xlab = list(label = "Quality", cex = 1.2),
        ylab = list(label = "Probability susceptible", cex = 1.2),
        scales = list(alternating = F, cex = 1.2),
        strip = F,
        pch = 4,
        col = c("black", "brown"),
        panel=function(...){
          panel.xyplot(...)
          panel.text(0.15, 0.97, labels[panel.number()], cex = 1.2, col = "black")
        }
      )

dev.off()

#===============================================================================
### Supplemental figure 4

rm(list = ls())

postscript("Manuscript/figure/supplement_4.eps", width = 8, height = 4)

load(paste(getwd(), "/Output/trap(lattice).RData", sep = ""))

longs <- unique(out$longevity)[c(4, 7, 10)]
labels <- c("a", "b", "c")

dat <- subset(out, xi_em == 0.5 & longevity %in% longs)

I <- ddply(dat, .(xi_im, longevity, quality), summarise, occ = median(I))

xyplot(occ ~ quality | longevity, groups = xi_im, data = I,
       xlab = list(label = "Quality", cex = 1.2),
       ylab = list(label = "Probability infectious", cex = 1.2),
       scales = list(alternating = F, cex = 1.2),
       strip = F,
       pch = 4,
       col = c("black", "brown"),
       panel=function(...){
         panel.xyplot(...)
         panel.text(0.15, 0.33, labels[panel.number()], cex = 1.2, col = "black")
       }
)

dev.off()


