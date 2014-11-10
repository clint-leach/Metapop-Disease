patch.plot <- function(conn){
  
  library(scales)

  load(paste(getwd(), "/Output/out", conn, ".RData", sep = ""))
  
  #Filter patch data to include only endemic cases
  patch<-out[which(out$Ifin>0 & out$Sfin > 0), ]
  
  #Filter patch data to include only high variance runs
  patch<-patch[which(abs(patch$variance - 0.2) < 0.01),]
  
  reps <- unique(patch$repID)
  
  pdf(paste(getwd(), "/Manuscript", "/infevents", conn, ".pdf", sep=""), height=6, width=8)
  
  par(mfrow = c(1, 1))
  
  plot(patch$quality, patch$inf.events.early, type = "n", ylab = "Number of infection events", xlab = "Patch quality", bty = "l")
  for(i in 1:length(reps)){
    sim <- patch[patch$repID == reps[i], ]
    try(lines(loess.smooth(sim$quality, sim$inf.events.tot), type = "l", col = alpha("gray", 0.3)))
  }
  
  lines(loess.smooth(patch$quality, patch$inf.events.tot), lwd = 2, col = "black")
 
  dev.off()
}