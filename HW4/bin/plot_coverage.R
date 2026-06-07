#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
depth_file <- paste0(sample_id, "_depth.txt")
output_png <- paste0(sample_id, "_coverage.png")

# Загрузка данных
data <- read.table(depth_file, header=FALSE, 
                  col.names=c("contig", "pos", "cov"))

# Создание графика
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

# Вывод статистики
cat(paste("\nStatistics for", sample_id, ":\n"))
cat(paste("  Mean coverage:", round(avg_cov, 1), "x\n"))
cat(paste("  Min coverage:", min(data$cov), "x\n"))
cat(paste("  Max coverage:", max(data$cov), "x\n"))