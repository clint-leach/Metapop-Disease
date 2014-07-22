

models <- c("(full)", "(lattice)", "(x0)", "(delta)")

source(paste(getwd(), "/Code/metapop.plot.r", sep = ""))
source(paste(getwd(), "/Code/highlow.plot.r", sep = ""))
source(paste(getwd(), "/Code/patch.plot.r", sep = ""))


for(i in 1:length(models)){
  metapop.plot(models[i])
}

for(i in 1:length(models)){
  highlow.plot(models[i])
}

for(i in 1:length(models)){
  patch.plot(models[i])
}




