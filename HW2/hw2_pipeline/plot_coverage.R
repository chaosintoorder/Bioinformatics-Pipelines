cat("Loading coverage.txt...\n")
data <- read.table("coverage.txt", header=FALSE, col.names=c("contig", "pos", "cov"))

data$cumulative_pos <- 0
current_offset <- 0
last_contig <- data$contig[1]

for (i in 1:nrow(data)) {
  if (data$contig[i] != last_contig) {
    current_offset <- current_offset + max(data$pos[data$contig == last_contig]) + 100
    last_contig <- data$contig[i]
  }
  data$cumulative_pos[i] <- current_offset + data$pos[i]
}

cat("Loaded", nrow(data), "positions\n")
cat("Number of contigs:", length(unique(data$contig)), "\n")

png("coverage_r.png", width=14, height=6, units="in", res=150)

plot(data$cumulative_pos, data$cov, 
     type="l", 
     col="blue",
     lwd=0.5,
     xlab="Genome Position (bp)",
     ylab="Coverage Depth",
     main="Coverage Plot",
     cex.lab=1.2)

avg_cov <- mean(data$cov)
abline(h=avg_cov, col="red", lty=2, lwd=1.5)

legend("topright", 
       legend=c(paste("Mean coverage:", round(avg_cov, 1), "x")),
       col="red", lty=2, lwd=1.5)

grid(col="gray", lty=3)

dev.off()

cat("Saved: coverage_r.png\n")
cat(paste("\nStatistics:\n"))
cat(paste("  Mean coverage:", round(avg_cov, 1), "x\n"))
cat(paste("  Min coverage:", min(data$cov), "x\n"))
cat(paste("  Max coverage:", max(data$cov), "x\n"))
cat(paste("  Total positions:", format(nrow(data), big.mark=","), "\n"))
cat(paste("  Number of contigs:", length(unique(data$contig)), "\n"))
