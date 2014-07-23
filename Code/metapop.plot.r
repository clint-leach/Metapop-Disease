metapop.plot <- function(conn){
  
  library(gridExtra)
  library(lattice)
  library(RColorBrewer)
  
  load(paste(getwd(), "/Output/out", conn, ".Rdata", sep = ""))
  
  metapop.S <- tapply(out$Sfin, list(out$longevity, out$variance), mean)
  metapop.I <- tapply(out$Ifin, list(out$longevity, out$variance), mean)
  metapop.occ <- tapply(out$Sfin + out$Ifin, list(out$longevity, out$variance), mean)
  
  cols <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  
  S <- levelplot(metapop.S, col.regions = cols, xlab = "Longevity", ylab = "Variance",
                 at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
  I <- levelplot(metapop.I, col.regions = cols, xlab = "Longevity", ylab = "Variance",
                 at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
  occ <- levelplot(metapop.occ, col.regions = cols, xlab = "Longevity", ylab = "Variance",
                   at = seq(0, 1, by = 0.01), colorkey = F, scales = list(at = c(1, 3, 5, 7, 9)))
  
  key <- draw.colorkey(list(col = cols, 
                            at = seq(0, 1, by = 0.01),
                            labels = list(at = seq(0, 1, by = 0.25)),
                            height = 0.7), 
                       draw = F)
  
  pdf(paste(getwd(), "/Manuscript", "/metapop", conn, ".pdf", sep=""), height=5, width=14)
  
  grid.arrange(S, I, occ, key, ncol = 4, widths = c(0.3, 0.3, 0.3, 0.1))
  
  dev.off()
}