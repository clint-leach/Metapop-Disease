#=================================================================================
#
#Code to analyze results from highlow scripts -- sensitivity analyses the remove
#high or low quality patches.  The treatment names are as follows:
#  low.var = missing highest and lowest 10 patches (replaced with 1)
#   low = low quality patches added to mid 
#   high = high quality patches added to mid 
#   both = high quality variance (low.var+low+high)
#=================================================================================

highlow<-read.csv("highlow(lattice).csv",header=TRUE)

highlow<-highlow[,-1]

#Outcome analysis

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

ggplot(outcomes, aes(x=trt, y=number)) + geom_bar(stat="identity") + facet_grid(.~outcome)


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


#==================================================================================================

#Digging into the role of initial quality

plot(both$quality0,both$S)
lines(loess.smooth(both$quality0,both$S))

plot(both$quality0,both$I)
lines(loess.smooth(both$quality0,both$I))
