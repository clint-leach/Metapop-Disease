# Code to analyze results from highlow scripts, in which the patch quality distribution
# varies asymmetrically around 1.  The four distributions explored are:
#   low var: patch qualities between 0.75 and 1.25
#   intermediate var, high quality: 0.75 - 1.75
#   intermediate var, low quality: 0.25 - 1.25
#   high var: 0.25 - 1.75

# Choose which parameterization to explore
conn <- "(full)"      
conn <- "(lattice)"
conn <- "(alpha0)"
conn <- "(delta)"

# Read in model output
highlow<-read.csv(paste(getwd(), "/Output/highlow", conn, ".csv", sep = ""), header=TRUE)

highlow<-highlow[,-1]

longevity <- sort(unique(highlow$longevity))

#===============================================================================
# Plotting mean occupancy as a function of longevity and treatment

par(mfrow = c(1, 3), bty = "l")

S.occ <- tapply(highlow$S, list(highlow$longevity, highlow$treatment), mean)
S.occ <- S.occ - S.occ[, 4]

plot(log10(longevity), S.occ[, 1], type = "b", lty = 1, ylim = c(-0.5, 0.5), pch = 19, col = "red", ylab = "S occupancy", xlab = "Longevity")
lines(log10(longevity), S.occ[, 2], type = "b", lty = 2, pch = 19, col = "blue")
lines(log10(longevity), S.occ[, 3], type = "b", lty = 3, pch = 19, col = "black")
lines(log10(longevity), S.occ[, 4], type = "l", lty = 4, pch = 19, col = "grey")

I.occ <- tapply(highlow$I, list(highlow$longevity, highlow$treatment), mean)
I.occ <- I.occ - I.occ[, 4]

plot(log10(longevity), I.occ[, 1], type = "b", lty = 1, ylim = c(-0.5, 0.5), pch = 19, col = "red", ylab = "I occupancy", xlab = "Longevity")
lines(log10(longevity), I.occ[, 2], type = "b", lty = 2, pch = 19, col = "blue")
lines(log10(longevity), I.occ[, 3], type = "b", lty = 3, pch = 19, col = "black")
lines(log10(longevity), I.occ[, 4], type = "l", lty = 4, pch = 19, col = "grey")

occ <- tapply(highlow$I + highlow$S, list(highlow$longevity, highlow$treatment), mean)
occ <- occ - occ[, 4]

plot(log10(longevity), occ[, 1], type = "b", lty = 1, ylim = c(-0.5, 0.5), pch = 19, col = "red", ylab = "Total occupancy", xlab = "Longevity")
lines(log10(longevity), occ[, 2], type = "b", lty = 2, pch = 19, col = "blue")
lines(log10(longevity), occ[, 3], type = "b", lty = 3, pch = 19, col = "black")
lines(log10(longevity), occ[, 4], type = "l", lty = 4, pch = 19, col = "grey")

#===============================================================================
# Plotting proportion of different epidemiological outcomes

par(mfrow = c(1, 4), bty = "l")

pandemic <- tapply(highlow$S == 0 & highlow$I > 0, list(highlow$longevity, highlow$treatment), sum)

plot(log10(longevity), pandemic[, 1], type = "b", lty = 1, ylim = c(0, 100), pch = 20, col = "red", xlab = "Longevity", ylab = "Pandemics")
lines(log10(longevity), pandemic[, 2], type = "b", lty = 2, pch = 20, col = "blue")
lines(log10(longevity), pandemic[, 3], type = "b", lty = 3, pch = 20, col = "black")
lines(log10(longevity), pandemic[, 4], type = "b", lty = 4, pch = 20, col = "grey")

endemic <- tapply(highlow$S > 0 & highlow$I > 0, list(highlow$longevity, highlow$treatment), sum)

plot(log10(longevity), endemic[, 1], type = "b", lty = 1, ylim = c(0, 100), pch = 20, col = "red", xlab = "Longevity", ylab = "Endemics")
lines(log10(longevity), endemic[, 2], type = "b", lty = 2, pch = 20, col = "blue")
lines(log10(longevity), endemic[, 3], type = "b", lty = 3, pch = 20, col = "black")
lines(log10(longevity), endemic[, 4], type = "b", lty = 4, pch = 20, col = "grey")

nodisease <- tapply(highlow$S > 0 & highlow$I == 0, list(highlow$longevity, highlow$treatment), sum)

plot(log10(longevity), nodisease[, 1], type = "b", lty = 1, ylim = c(0, 100), pch = 20, col = "red", xlab = "Longevity", ylab = "No disease")
lines(log10(longevity), nodisease[, 2], type = "b", lty = 2, pch = 20, col = "blue")
lines(log10(longevity), nodisease[, 3], type = "b", lty = 3, pch = 20, col = "black")
lines(log10(longevity), nodisease[, 4], type = "b", lty = 4, pch = 20, col = "grey")

extinct <- tapply(highlow$S == 0 & highlow$I == 0, list(highlow$longevity, highlow$treatment), sum)

plot(log10(longevity), extinct[, 1], type = "b", lty = 1, ylim = c(0, 100), pch = 20, col = "red", xlab = "Longevity", ylab = "Extinctions")
lines(log10(longevity), extinct[, 2], type = "b", lty = 2, pch = 20, col = "blue")
lines(log10(longevity), extinct[, 3], type = "b", lty = 3, pch = 20, col = "black")
lines(log10(longevity), extinct[, 4], type = "b", lty = 4, pch = 20, col = "grey")



#===============================================================================
# Analysis and plots for a single longevity

highlow <- highlow[highlow$longevity == 100, ]

low<-highlow[which(highlow$treatment=="+low"),]
high<-highlow[which(highlow$treatment=="+high"),]
full<-highlow[which(highlow$treatment=="full"),]
low.var<-highlow[which(highlow$treatment=="low.var"),]


outcome<-function(data){
  extinction<-length(which(data$S==0 & data$I==0))
  nodisease<-length(which(data$S>0 & data$I==0))
  endemic<-length(which(data$S>0 & data$I>0))
  pandemic<-length(which(data$S==0 & data$I>0))
  return(c(extinction,nodisease,endemic,pandemic))
}

outcomes<-data.frame(cbind(outcome(low.var), outcome(low), outcome(high), outcome(full)))

colnames(outcomes)<-c("low.var","low","high","high.var")

outcomes<-cbind(factor(c("Extinction","Disease-Free","Endemic","Pandemic"), levels=c("Extinction","Disease-Free", "Endemic", "Pandemic")),stack(outcomes))
names(outcomes)<-c("outcome", "number", "trt")

outcomes$trt<-factor(outcomes$trt, levels = c("low.var", "low", "high", "high.var"))

ggplot(outcomes, aes(x=outcome, y=number)) + geom_bar(stat="identity") + facet_grid(.~trt)


par(mfrow=c(2,4))
hist(mid$I,xlab="Infected Occupancy",main=("Low Variance"))
hist(high$I,xlab="Infected Occupancy",main=("+High quality"))
hist(low$I,xlab="Infected Occupancy",main=("+Low quality"))
hist(both$I,xlab="Infected Occupancy",main=("High variance"))

hist(mid$S,xlab="Infected Occupancy",main=("Low Variance"))
hist(high$S,xlab="Infected Occupancy",main=("+High quality"))
hist(low$S,xlab="Infected Occupancy",main=("+Low quality"))
hist(both$S,xlab="Infected Occupancy",main=("High variance"))

#===============================================================================
#Quantitative analysis

mean.S<-tapply(highlow$S,highlow$treatment,mean)
sd.S<-tapply(highlow$S,highlow$treatment,sd)

mean.I<-tapply(highlow$I,highlow$treatment,mean)
sd.I<-tapply(highlow$I,highlow$treatment,sd)

mean.R<-tapply(highlow$R,highlow$treatment,mean)
sd.R<-tapply(highlow$R,highlow$treatment,sd)

mean.PO<-tapply(highlow$S+highlow$I,highlow$treatment,mean)
sd.PO<-tapply(highlow$S+highlow$I,highlow$treatment,sd)

mean.maxI<-tapply(highlow$maxI,highlow$treatment,mean)
sd.maxI<-tapply(highlow$maxI,highlow$treatment,sd)

mean.tmaxI<-tapply(highlow$tmaxI,highlow$treatment,mean)
sd.tmaxI<-tapply(highlow$tmaxI,highlow$treatment,sd)


#Calculating Scheffe simultaneous CI

lm.S<-lm(S~as.factor(treatment),data=highlow)
MSE.S<-sum(lm.S$residuals^2)/396
spread.S<-sqrt(MSE.S*0.02*3*qf(0.95,3,396))

lm.I<-lm(I~as.factor(treatment),data=highlow)
MSE.I<-sum(lm.I$residuals^2)/396
spread.I<-sqrt(MSE.I*0.02*3*qf(0.95,3,396))


#================================================================================
#ggplot2 figures

library(ggplot2)

#Plotting susceptible occupancy

df.S<-data.frame(
  trt<-factor(c("low var","+ high/low","+ high/low","+ both")),
  S<-c(mean.S["low.var"],mean.S["+high"],mean.S["+low"],mean.S["full"]),
  group<-factor(c("control","+high","+low","control")),
  se<-rep(spread.S,4)
)

names(df.S)<-c("treatment","S","group","se")

limits<-aes(ymax=S+se,ymin=S-se)

p<-ggplot(df.S,aes(colour=group,y=S,x=trt))
p<-p+geom_point(size=5)
p<-p+geom_errorbar(limits,width=0.2,size=1)
p<-p+opts(axis.title.x=theme_text(colour="black",size=20),
       axis.text.x=theme_text(colour="black",size=20),
       axis.text.y=theme_text(colour="black",size=20),
       axis.title.y=theme_text(size=20,angle=90,vjust=0.25),
       panel.grid.minor=theme_blank(),
       panel.grid.major=theme_blank(),
       panel.background=theme_rect(fill="white"),
       legend.text=theme_text(size=20))
p<-p+xlim(c("low var","+ high/low","+ both"))+xlab("")+ylab("Susceptible Occupancy") 
p+scale_colour_manual(values=c("control"="black","+high"="blue","+low"="red"),name="",breaks=c("+high","+low"))


#Plotting infected occupancy

df.I<-data.frame(
  trt<-factor(c("low var","+ high/low","+ high/low","+ both")),
  S<-c(mean.I["low.var"],mean.I["+high"],mean.I["+low"],mean.I["full"]),
  group<-factor(c("control","+high","+low","control")),
  se<-rep(spread.I,4)
  )

names(df.I)<-c("treatment","I","group","se")

limits<-aes(ymax=I+se,ymin=I-se)

p<-ggplot(df.I,aes(colour=group,y=I,x=trt))
p<-p+geom_point(size=5)
p<-p+geom_errorbar(limits,width=0.2,size=1)
p<-p+opts(axis.title.x=theme_text(colour="black",size=20),
          axis.text.x=theme_text(colour="black",size=20),
          axis.text.y=theme_text(colour="black",size=20),
          axis.title.y=theme_text(size=20,angle=90,vjust=0.25),
          panel.grid.minor=theme_blank(),
          panel.grid.major=theme_blank(),
          panel.background=theme_rect(fill="white"),
          legend.text=theme_text(size=20))
p<-p+xlim(c("low var","+ high/low","+ both"))+xlab("") +ylab("Infected Occupancy")
p+scale_colour_manual(values=c("control"="black","+high"="blue","+low"="red"),name="",breaks=c("+high","+low"))


#Plotting total occupancy

df.PO<-data.frame(
  trt<-factor(c("low var","+ high/low","+ high/low","+ both")),
  PO<-c(mean.PO["mid"],mean.PO["high"],mean.PO["low"],mean.PO["both"]),
  group<-factor(c("control","+high","+low","control")),
  se<-c(sd.PO["mid"]/10,sd.PO["high"]/10,sd.PO["low"]/10,sd.PO["both"]/10)
)

names(df.PO)<-c("treatment","PO","group","se")

limits<-aes(ymax=PO+se,ymin=PO-se)

p<-ggplot(df.PO,aes(colour=group,y=PO,x=trt))
p<-p+geom_point(size=5)
p<-p+geom_errorbar(limits,width=0.2,size=1)
p<-p+opts(axis.title.x=theme_text(colour="black",size=20),
          axis.text.x=theme_text(colour="black",size=20),
          axis.text.y=theme_text(colour="black",size=20),
          axis.title.y=theme_text(size=20,angle=90,vjust=0.25),
          panel.grid.minor=theme_blank(),
          panel.grid.major=theme_blank(),
          panel.background=theme_rect(fill="white"),
          legend.text=theme_text(size=20))
p<-p+xlim(c("low var","+ high/low","+ both"))+xlab("") +ylab("Total Occupancy")
p+scale_colour_manual(values=c("control"="black","+high"="blue","+low"="red"),name="",breaks=c("+high","+low"))



#Plotting susceptible and infected occupancy side-by-side

df<-data.frame(
  trt<-factor(c("low var","+ high/low","+ high/low","+ both")),
  y<-c(mean.S["low.var"],mean.S["+high"],mean.S["+low"],mean.S["full"],
       mean.I["low.var"],mean.I["+high"],mean.I["+low"],mean.I["full"]),
  group<-factor(rep(c("control","+high","+low","control"),2)),
  se<-c(rep(spread.S,4),rep(spread.I,4)),
  case<-rep(c("Susceptible","Infected"),each=4)
)

names(df)<-c("treatment","S","group","se","case")

limits<-aes(ymax=S+se,ymin=S-se)

p<-ggplot(df,aes(colour=group,y=S,x=treatment))
p<-p+geom_point(size=5)
p<-p+geom_errorbar(limits,width=0.2,size=1)
p<-p+opts(axis.title.x=theme_text(colour="black",size=20),
          axis.text.x=theme_text(colour="black",size=20),
          axis.text.y=theme_text(colour="black",size=20),
          axis.title.y=theme_text(size=20,angle=90,vjust=0.25),
          panel.grid.minor=theme_blank(),
          panel.grid.major=theme_blank(),
          panel.background=theme_rect(fill="white"),
          legend.text=theme_text(size=20),
          strip.text.x=theme_text(colour="black",size=20),
          strip.background=theme_rect(fill="white",colour="white"))
p<-p+xlim(c("low var","+ high/low","+ both"))+xlab("")+ylab("Occupancy") 
p<-p+facet_grid(.~case)
p+scale_colour_manual(values=c("control"="black","+high"="blue","+low"="red"),name="",breaks=c("+high","+low"))


#===============================================================================

#Digging into the role of initial quality

plot(both$quality0,both$S)
lines(loess.smooth(both$quality0,both$S))

plot(both$quality0,both$I)
lines(loess.smooth(both$quality0,both$I))

#===============================================================================

par(mfrow = c(1, 5))

for(i in c("(full)", "(lattice)", "(alpha0)", "(delta)", "(ei)")){
  
  highlow<-read.csv(paste(getwd(), "/Output/highlow", i, ".csv", sep = ""), header=TRUE)
  
  S.occ <- tapply(highlow$S, list(highlow$longevity, highlow$treatment), mean)
  I.occ <- tapply(highlow$I, list(highlow$longevity, highlow$treatment), mean)
  occ <- tapply(highlow$I + highlow$S, list(highlow$longevity, highlow$treatment), mean)
  
  plot(log10(longevity), occ[, 1] - occ[, 2], type = "b", pch = 19, main = i, bty = "l",
       ylim = c(-0.3, 0.5), xlab = "longevity", ylab = "High quality occupancy - low quality occupancy")
  abline(h = 0, col = "red")

}


