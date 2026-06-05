#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
depth_file <- paste0(sample_id, ".coverage.tsv")
output_png <- paste0(sample_id, ".coverage.png")
summary_tsv <- paste0(sample_id, ".coverage_summary.tsv")

data <- read.table(depth_file, header=FALSE, 
                  col.names=c("contig", "pos", "cov"))

png(output_png, width=14, height=6, units="in", res=150)

plot(data$pos, data$cov, 
     type="l", 
     col="blue",
     lwd=0.5,
     xlab="Genome Position (bp)", 
     ylab="Coverage Depth", 
     main=paste("Coverage Plot for", sample_id),
     cex.lab=1.2)

avg_cov <- mean(data$cov)
abline(h=avg_cov, col="red", lty=2, lwd=1.5)

legend("topright", 
       legend=c(paste("Mean coverage:", round(avg_cov, 1), "x")),
       col="red", lty=2, lwd=1.5)

grid(col="gray", lty=3)

dev.off()

# Write summary
write.table(data.frame(
    sample = sample_id,
    contig = data$contig[1],
    length = nrow(data),
    mean_depth = round(avg_cov, 3),
    covered_bases = sum(data$cov > 0),
    breadth = round(sum(data$cov > 0) / nrow(data), 5)
), summary_tsv, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

cat(paste("\nStatistics for", sample_id, ":\n"))
cat(paste("  Mean coverage:", round(avg_cov, 1), "x\n"))
cat(paste("  Min coverage:", min(data$cov), "x\n"))
cat(paste("  Max coverage:", max(data$cov), "x\n"))
