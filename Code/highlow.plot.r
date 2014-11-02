highlow.plot <- function(conn){
  
  highlow<-read.csv(paste(getwd(), "/Output/highlow", conn, ".csv", sep = ""), header=TRUE)
  
  longevity <- sort(unique(highlow$longevity))
    
  S.occ <- tapply(highlow$S, list(highlow$longevity, highlow$treatment), median)
  I.occ <- tapply(highlow$I, list(highlow$longevity, highlow$treatment), median)
  occ <- tapply(highlow$I + highlow$S, list(highlow$longevity, highlow$treatment), median)
  
  pdf(paste(getwd(), "/Manuscript", "/highlow", conn, ".pdf", sep=""), height=5, width=10)
  
  par(mfrow = c(1, 3), bty = "l")
  plot(log10(longevity), S.occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red", 
       ylab = "S occupancy", xlab = "log(longevity)", main = title("a", adj = 0))
  lines(log10(longevity), S.occ[, 2], type = "b", lty = 1, pch = 19, col = "blue")
  
  plot(log10(longevity), I.occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red", 
       ylab = "I occupancy", xlab = "log(longevity)", main = title("b", adj = 0))
  lines(log10(longevity), I.occ[, 2], type = "b", lty = 1, pch = 19, col = "blue")
  
  occ <- tapply(highlow$I + highlow$S, list(highlow$longevity, highlow$treatment), mean)
  
  plot(log10(longevity), occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red",
       ylab = "Total occupancy", xlab = "log(longevity)", main = title("c", adj = 0))
  lines(log10(longevity), occ[, 2], type = "b", lty = 1, pch = 19, col = "blue")
  
  dev.off()
}