require(MaizeGeneFamilies)

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/dosageEffectPlots.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly=TRUE)


#'##############
#' ALL TISSUES #
#' Paralags    #
#'##############

# Scatter plot with regression line:
pdf(file.path(input.args[[1]], "inst", "copyNoRpkmScatter.pdf"))
p.df <- gene.copy.no.rpkm.corr.df
plot(x = p.df$Copies, y = p.df$RPKM, type = "p", pch = 20, xlab = "Copies", ylab = "mean RPKM of all tissues")
abline(gene.copy.no.rpkm.corr$lm, col = "red", lwd = 3)
text(x = max(p.df$Copies) * 0.9, y = max(p.df$RPKM) * 0.9, labels = paste("R^2", 
    "=", sprintf("%.2f", summary(gene.copy.no.rpkm.corr$lm)$r.squared)), col = "blue")
dev.off()
# Log2 scaled:
p.df$Copies <- log2(p.df$Copies)
p.df$RPKM <- log2(p.df$RPKM)
log.lm <- lm(formula = RPKM ~ Copies, data = p.df)
pdf(file.path(input.args[[1]], "inst", "copyNoRpkmLog2Scatter.pdf"))
plot(x = p.df$Copies, y = p.df$RPKM, type = "p", pch = 20, xlab = "log2( Copies )", 
    ylab = "log2( mean RPKM of all tissues )")
abline(log.lm, col = "red", lwd = 3)
text(x = max(p.df$Copies) * 0.9, y = max(p.df$RPKM) * 0.9, labels = paste("R^2", 
    "=", sprintf("%.2f", summary(log.lm)$r.squared)), col = "blue")
dev.off()

#'##############
#' ALL TISSUES #
#' TANDEMS     #
#'##############

# Scatter plot with regression line:
pdf(file.path(input.args[[1]], "inst", "tandNoRpkmScatter.pdf"))
p.df <- gene.tand.no.rpkm.corr.df
plot(x = p.df$Copies, y = p.df$RPKM, type = "p", pch = 20, xlab = "n Tandems", ylab = "mean RPKM of all tissues")
abline(gene.copy.no.rpkm.corr$lm, col = "red", lwd = 3)
text(x = max(p.df$Copies) * 0.9, y = max(p.df$RPKM) * 0.9, labels = paste("R^2", 
    "=", sprintf("%.2f", summary(gene.copy.no.rpkm.corr$lm)$r.squared)), col = "blue")
dev.off()
# Log2 scaled:
p.df$Copies <- log2(p.df$Copies)
p.df$RPKM <- log2(p.df$RPKM)
log.lm <- lm(formula = RPKM ~ Copies, data = p.df)
pdf(file.path(input.args[[1]], "inst", "tandNoRpkmLog2Scatter.pdf"))
plot(x = p.df$Copies, y = p.df$RPKM, type = "p", pch = 20, xlab = "log2( n Tandems )", 
    ylab = "log2( mean RPKM of all tissues )")
abline(log.lm, col = "red", lwd = 3)
text(x = max(p.df$Copies) * 0.9, y = max(p.df$RPKM) * 0.9, labels = paste("R^2", 
    "=", sprintf("%.2f", summary(log.lm)$r.squared)), col = "blue")
dev.off()


message("DONE")
