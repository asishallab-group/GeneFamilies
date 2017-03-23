#' Scatter-Plots the pecies specific annotation (function) diversity for each
#' single gene group found in 'entropy.df'.
#'
#' @param plot.pdf The PDF File into which to save the scatter-plot
#' @param entropy.df A data.frame with species specific annotation entropies
#' per gene group. Gene groups are rows and columns need to include the species
#' in argument 'species'.
#' @param species A character of length 2 holding the names of the two species
#' to use for the plot. The species names need to occur as column names in
#' 'entropy.df'.
#'
#' @export
#' @return TRUE if and only if the plot can be generated successfully 
plotSpeciesSpecificAnnotationDiversity <- function(plot.pdf, entropy.df, species) {
    axis.rng <- c(0, max(unlist(entropy.df[, species]), na.rm = TRUE))
    pdf(plot.pdf)
    plot(x = entropy.df[, species[[1]]], entropy.df[, species[[2]]], pch = 20, xlab = species[[1]], 
        ylab = species[[2]], xlim = axis.rng, ylim = axis.rng)
    abline(0, 1, col = "brown")
    dev.off()
    TRUE
}

#' Boxplots the distributions of species specific annotation (function)
#' diversity, measured for each single gene group found in 'entropy.df'.
#'
#' @param plot.pdf The PDF File into which to save the distributions
#' @param entropy.df A data.frame with species specific annotation entropies
#' per gene group. Gene groups are rows and columns need to include the species
#' in argument 'species'.
#' @param species A character of length 2 or greater holding the names of the species
#' to use for the plot. The species names need to occur as column names in
#' 'entropy.df'.
#' @param gene.groups.type The type of gene groups: e.g. 'Tandem Clusters'.
#' Will be used in the main title of the boxplot.
#' @param star.species The names of those species above whose boxplots to print
#' a '*'. Default is none (\code{c()}).
#' @param colors A character vector holding the colors for each species'
#' density. Default is the brewer palette 'Dark2' for length(species).
#'
#' @export
#' @return TRUE if and only if the plot can be generated successfully 
plotSpeciesSpecificAnnotationDiversityDists <- function(plot.pdf, entropy.df, species, 
    gene.groups.type, star.species = c(), colors = brewer.pal(length(species), "Dark2")) {
    plot.df <- do.call("rbind", lapply(species, function(x) data.frame(entropy = unlist(entropy.df[, 
        x]), species = x, stringsAsFactors = FALSE)))
    y.max <- max(plot.df$entropy, na.rm = TRUE)
    pdf(plot.pdf)
    boxplot(entropy ~ species, data = plot.df, border = colors, main = paste("Species specific diversity of compound annotations per", 
        gene.groups.type, sep = "\n"), ylab = "Shannon-Entropy", xlab = "Species", 
        ylim = c(0, (y.max + 0.25)), pch = "-", cex = 0.25)
    txt.specs <- rep( "", length(species) )
    txt.specs.inds <- which( sort( species ) %in% star.species )
    txt.specs[txt.specs.inds] <- "*"
    text(x = 1:length(species), y = y.max + 0.25, txt.specs, cex=3, col='brown')
    dev.off()
}
