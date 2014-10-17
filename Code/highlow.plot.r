highlow.plot <- function(conn){
  
  highlow<-read.csv(paste(getwd(), "/Output/highlow", conn, ".csv", sep = ""), header=TRUE)
    
  S.occ <- tapply(highlow$S, list(highlow$longevity, highlow$treatment), mean)
  I.occ <- tapply(highlow$I, list(highlow$longevity, highlow$treatment), mean)
  occ <- tapply(highlow$I + highlow$S, list(highlow$longevity, highlow$treatment), mean)
  
  pdf(paste(getwd(), "/Manuscript", "/highlow", conn, ".pdf", sep=""), height=5, width=10)
  
  par(mfrow = c(1, 3), bty = "l")
  plot(rownames(S.occ), S.occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red", 
       ylab = "S occupancy", xlab = "Longevity", main = title("a", adj = 0))
  lines(rownames(S.occ), S.occ[, 2], type = "b", lty = 1, pch = 19, col = "blue")
  lines(rownames(S.occ), S.occ[, 3], type = "b", lty = 1, pch = 19, col = "black")
  lines(rownames(S.occ), S.occ[, 4], type = "b", lty = 1, pch = 19, col = "grey")
  
  plot(rownames(I.occ), I.occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red", 
       ylab = "I occupancy", xlab = "Longevity", main = title("b", adj = 0))
  lines(rownames(I.occ), I.occ[, 2], type = "b", lty = 1, pch = 19, col = "blue")
  lines(rownames(I.occ), I.occ[, 3], type = "b", lty = 1, pch = 19, col = "black")
  lines(rownames(I.occ), I.occ[, 4], type = "b", lty = 1, pch = 19, col = "grey")
  
  occ <- tapply(highlow$I + highlow$S, list(highlow$longevity, highlow$treatment), mean)
  
  plot(rownames(occ), occ[, 1], type = "b", lty = 1, ylim = c(0, 1), pch = 19, col = "red",
       ylab = "Total occupancy", xlab = "Longevity", main = title("c", adj = 0))
  lines(rownames(occ), occ[, 2], type = "b", lty = 1, pch = 19, col = "blue")
  lines(rownames(occ), occ[, 3], type = "b", lty = 1, pch = 19, col = "black")
  lines(rownames(occ), occ[, 4], type = "b", lty = 1, pch = 19, col = "grey")
  
  dev.off()
}