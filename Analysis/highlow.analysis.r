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
conn <- "(ei)"
conn <- "(midei)"
conn <- "(lowei)"

# Read in model output
highlow<-read.csv(paste(getwd(), "/Output/highlow", conn, ".csv", sep = ""), header=TRUE)

highlow<-highlow[,-1]

longevity <- sort(unique(highlow$longevity))

#===============================================================================
# Plotting mean occupancy as a function of longevity and treatment

par(mfrow = c(1, 3), bty = "l")

S.occ <- tapply(highlow$S.pop, list(highlow$longevity, highlow$treatment), median)

plot(log10(longevity), S.occ[, "+high"], type = "b", lty = 1, ylim = c(0, 120), pch = 19, col = "red", ylab = "S occupancy", xlab = "Longevity")
lines(log10(longevity), S.occ[, "+low"], type = "b", lty = 2, pch = 19, col = "blue")
lines(log10(longevity), S.occ[, "full"], type = "b", lty = 3, pch = 19, col = "black")
lines(log10(longevity), S.occ[, "low.var"], type = "l", lty = 4, pch = 19, col = "grey")

I.occ <- tapply(highlow$I.pop, list(highlow$longevity, highlow$treatment), median)

plot(log10(longevity), I.occ[, "+high"], type = "b", lty = 1, ylim = c(0, 120), pch = 19, col = "red", ylab = "I occupancy", xlab = "Longevity")
lines(log10(longevity), I.occ[, "+low"], type = "b", lty = 2, pch = 19, col = "blue")
lines(log10(longevity), I.occ[, "full"], type = "b", lty = 3, pch = 19, col = "black")
lines(log10(longevity), I.occ[, "low.var"], type = "l", lty = 4, pch = 19, col = "grey")

occ <- tapply(highlow$I.pop + highlow$S.pop, list(highlow$longevity, highlow$treatment), median)

plot(log10(longevity), occ[, "+high"], type = "b", lty = 1, ylim = c(0, 120), pch = 19, col = "red", ylab = "Total occupancy", xlab = "Longevity")
lines(log10(longevity), occ[, "+low"], type = "b", lty = 2, pch = 19, col = "blue")
lines(log10(longevity), occ[, "full"], type = "b", lty = 3, pch = 19, col = "black")
lines(log10(longevity), occ[, "low.var"], type = "l", lty = 4, pch = 19, col = "grey")

#===============================================================================
# Plotting proportion of different epidemiological outcomes

par(mfrow = c(1, 4), bty = "l")

pandemic <- tapply(highlow$S == 0 & highlow$I > 0, list(highlow$longevity, highlow$treatment), sum)

plot(log10(longevity), pandemic[, "+high"], type = "b", lty = 1, ylim = c(0, 100), pch = 20, col = "red", xlab = "Longevity", ylab = "Pandemics")
lines(log10(longevity), pandemic[, "+low"], type = "b", lty = 2, pch = 20, col = "blue")
lines(log10(longevity), pandemic[, "full"], type = "b", lty = 3, pch = 20, col = "black")
lines(log10(longevity), pandemic[, "low.var"], type = "b", lty = 4, pch = 20, col = "grey")

endemic <- tapply(highlow$S > 0 & highlow$I > 0, list(highlow$longevity, highlow$treatment), sum)

plot(log10(longevity), endemic[, "+high"], type = "b", lty = 1, ylim = c(0, 100), pch = 20, col = "red", xlab = "Longevity", ylab = "Endemics")
lines(log10(longevity), endemic[, "+low"], type = "b", lty = 2, pch = 20, col = "blue")
lines(log10(longevity), endemic[, "full"], type = "b", lty = 3, pch = 20, col = "black")
lines(log10(longevity), endemic[, "low.var"], type = "b", lty = 4, pch = 20, col = "grey")

nodisease <- tapply(highlow$S > 0 & highlow$I == 0, list(highlow$longevity, highlow$treatment), sum)

plot(log10(longevity), nodisease[, "+high"], type = "b", lty = 1, ylim = c(0, 100), pch = 20, col = "red", xlab = "Longevity", ylab = "No disease")
lines(log10(longevity), nodisease[, "+low"], type = "b", lty = 2, pch = 20, col = "blue")
lines(log10(longevity), nodisease[, "full"], type = "b", lty = 3, pch = 20, col = "black")
lines(log10(longevity), nodisease[, "low.var"], type = "b", lty = 4, pch = 20, col = "grey")

extinct <- tapply(highlow$S == 0 & highlow$I == 0, list(highlow$longevity, highlow$treatment), sum)

plot(log10(longevity), extinct[, "+high"], type = "b", lty = 1, ylim = c(0, 100), pch = 20, col = "red", xlab = "Longevity", ylab = "Extinctions")
lines(log10(longevity), extinct[, "+low"], type = "b", lty = 2, pch = 20, col = "blue")
lines(log10(longevity), extinct[, "full"], type = "b", lty = 3, pch = 20, col = "black")
lines(log10(longevity), extinct[, "low.var"], type = "b", lty = 4, pch = 20, col = "grey")



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
sim <- c("(full)", "(lattice)", "(alpha0)", "(delta)", "(ei)")
label <- c("a", "b", "c", "d","e")

highlow<-read.csv(paste(getwd(), "/Output/highlow",sim[1], ".csv", sep = ""), header=TRUE)

longevity <- sort(unique(highlow$longevity))

S.occ <- tapply(highlow$S.pop, list(highlow$longevity, highlow$treatment), median)
I.occ <- tapply(highlow$I.pop, list(highlow$longevity, highlow$treatment), median)
occ <- tapply(highlow$I.pop + highlow$S.pop, list(highlow$longevity, highlow$treatment), median)

plot(log10(longevity), occ[, "+high"] - occ[, "+low"], type = "b", pch = 19, bty = "l",
    ylim = c(-60, 60), xlab = "log(longevity)", ylab = "High quality occupancy - low quality occupancy",
    main = title(label[1], adj = 0),
    cex.lab = 1.2)

abline(h = 0, col = "red")

for(i in 2:5){
  
  highlow<-read.csv(paste(getwd(), "/Output/highlow",sim[i], ".csv", sep = ""), header=TRUE)
  
  longevity <- sort(unique(highlow$longevity))
  
  S.occ <- tapply(highlow$S.pop, list(highlow$longevity, highlow$treatment), median)
  I.occ <- tapply(highlow$I.pop, list(highlow$longevity, highlow$treatment), median)
  occ <- tapply(highlow$I.pop + highlow$S.pop, list(highlow$longevity, highlow$treatment), median)
  
  plot(log10(longevity), occ[, "+high"] - occ[, "+low"], type = "b", pch = 19, bty = "l",
       ylim = c(-60, 60), xlab = "log(longevity)", ylab = "", yaxt = "n",
       main = title(label[i], adj = 0),
       cex.lab = 1.2)
  Axis(side= 2, labels = F)
  
  abline(h = 0, col = "red")

}

#===============================================================================
# Occupancy boxplots

library(ggplot2)
library(gridExtra)

sub <- highlow[highlow$treatment %in% c("+high", "+low"), ]

p.S <- ggplot(sub, aes(x = factor(log10(longevity)), y = S.pop))
p.S <- p.S + xlab("log(longevity)") + ylab("S occupancy") + theme(legend.position = "none") + theme_bw()
p.S <- p.S + geom_boxplot(aes(fill = factor(treatment)))

p.I <- ggplot(sub, aes(x = factor(log10(longevity)), y = I.pop))
p.I <- p.I + xlab("log(longevity)") + ylab("I occupancy") + theme(legend.position = "none")
p.I <- p.I + geom_boxplot(aes(fill = factor(treatment)))

p.occ <- ggplot(sub, aes(x = factor(log10(longevity)), y = S.pop + I.pop))
p.occ <- p.occ + xlab("log(longevity)") + ylab("Total occupancy") + theme(legend.position = "none")
p.occ <- p.occ + geom_boxplot(aes(fill = factor(treatment)))

grid.arrange(p.S, p.I, p.occ, ncol = 3)

