highlow.plot <- function(conn){
  
  highlow<-read.csv(paste(getwd(), "/Data/highlow", conn, ".csv", sep = ""), header=TRUE)
    
  S.occ <- tapply(highlow$S, list(highlow$longevity, highlow$treatment), mean)
  I.occ <- tapply(highlow$I, list(highlow$longevity, highlow$treatment), mean)
  occ <- tapply(highlow$I + highlow$S, list(highlow$longevity, highlow$treatment), mean)
  
  pdf(paste(getwd(), "/Manuscript", "/highlow", conn, ".pdf", sep=""), height=5, width=10)
  
  par(mfrow = c(1, 3), bty = "l")
  plot(seq(20, 200, by = 20), S.occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red", ylab = "S occupancy", xlab = "Longevity")
  lines(seq(20, 200, by = 20), S.occ[, 2], type = "b", lty = 2, pch = 19, col = "blue")
  lines(seq(20, 200, by = 20), S.occ[, 3], type = "b", lty = 3, pch = 19, col = "black")
  lines(seq(20, 200, by = 20), S.occ[, 4], type = "l", lty = 4, pch = 19, col = "grey")
  
  plot(seq(20, 200, by = 20), I.occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red", ylab = "I occupancy", xlab = "Longevity")
  lines(seq(20, 200, by = 20), I.occ[, 2], type = "b", lty = 2, pch = 19, col = "blue")
  lines(seq(20, 200, by = 20), I.occ[, 3], type = "b", lty = 3, pch = 19, col = "black")
  lines(seq(20, 200, by = 20), I.occ[, 4], type = "l", lty = 4, pch = 19, col = "grey")
  
  occ <- tapply(highlow$I + highlow$S, list(highlow$longevity, highlow$treatment), mean)
  
  plot(seq(20, 200, by = 20), occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red", ylab = "Total occupancy", xlab = "Longevity")
  lines(seq(20, 200, by = 20), occ[, 2], type = "b", lty = 2, pch = 19, col = "blue")
  lines(seq(20, 200, by = 20), occ[, 3], type = "b", lty = 3, pch = 19, col = "black")
  lines(seq(20, 200, by = 20), occ[, 4], type = "l", lty = 4, pch = 19, col = "grey")
  
  dev.off()
}