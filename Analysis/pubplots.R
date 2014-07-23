# Script to generate all plots for manuscript

models <- c("(full)", "(lattice)", "(alpha0)", "(delta)")

source("Code/metapop.plot.r")
source("Code/highlow.plot.r")
source("Code/patch.plot.r")


for(i in 1:length(models)){
  metapop.plot(models[i])
}

for(i in 1:length(models)){
  highlow.plot(models[i])
}

for(i in 1:length(models)){
  patch.plot(models[i])
}




