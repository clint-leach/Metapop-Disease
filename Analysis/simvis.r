
# Code to generate visualization of a single simulation run

source("Code/contSPOM(rates).r")

set.seed(8689)

# Defines fixed inputs for diseaseSPOM
parms<-data.frame("xi_im" = 0.5, 
                  "xi_em" = 0.5,      
                  "D" = 2, 
                  "alpha" = 1,
                  "es" = 0.1,
                  "nu" = 0.2,
                  "delta" = 0.5,
                  "gamma0" = 0.5)

#Fully connected matrix
# distance<-matrix(1,nrow=100,ncol=100)
# diag(distance)<-0

# Lattice matrix (on a torus)
library(igraph)
distance<-get.adjacency(graph.lattice(c(10,10),directed=F,circular=T))
distance<-as.matrix(distance)

# Initial conditions
initial<-c(1:100)
initial[1:50]<-"S"
initial[51:100]<-"E"

timesteps<-5000

longevity <- 10 ^ 0.5
variance <- 0.2

#===============================================================================

#Generates quality vector with desired variance
max <- 0.5 * sqrt(12 * variance) + 1
min <- 2 - max
quality <- runif(100, min=min, max=max)

#Converts longevity (half-life) to decay rate
long.nat <- longevity / parms$es
r <- log(2) / long.nat

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

#===============================================================================

state <- output[, -1]

susc <- state == "S"
inf <- state == "I"
empty <- state == "E" | state == "R"

S <- rowSums(susc) / 100
I <- rowSums(inf) / 100
E <- rowSums(empty) / 100

#===============================================================================

library(ggplot2)
library(reshape2)
library(gridExtra)

state <- state[, order(quality)]
rownames(state) <- c(1:(timesteps + 1))
colnames(state) <- c(1:100)

state <- melt(state)
names(state) <- c("Time", "Patch", "color")
levels(state$color) <- c("white", "#ee2c2c", "white", "#1874cd")

p1 <- ggplot(state, aes(x = Time, y = Patch)) + geom_tile(fill = state$color) + 
              scale_fill_identity() + 
              scale_x_continuous(expand = c(0, 0)) +
              scale_y_continuous(expression(Patch ~ Quality ~ symbol('\256')), expand = c(0, 0)) +
              theme_classic()


tseries <- data.frame(Time = rep(c(1:(timesteps + 1)), 2), Occupancy = c(S, I), class = rep(c("S", "I"), each = (timesteps + 1)))

p2 <- ggplot(tseries, aes(x = Time, y = Occupancy, group = class, color = class)) + geom_line() + 
              scale_color_manual(values = c("#ee2c2c", "#1874cd")) +
              guides(color = FALSE, alpha = FALSE) +
              theme_classic() + labs(x = NULL) + 
              scale_x_continuous(labels = NULL, expand = c(0, 0)) +
              scale_y_continuous(expand = c(0, 0))

postscript("Manuscript/figure/figure_4.eps", width = 12, height = 8)

grid.arrange(p2, p1, ncol = 1, heights = c(0.25, 0.75))

dev.off()

# MoMA version
p1 + theme(line = element_blank(), text = element_blank())
