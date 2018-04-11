require(MaizeGeneFamilies)
options(mc.cores = getMcCores())

message("USAGE: Rscript path/2/MaizeGeneFamilies/exec/domainBasedFunctionEvolutionAfterDuplication.R path/2/MaizeGeneFamilies")

input.args <- commandArgs(trailingOnly = TRUE)


#' Gene Sets have been defined in ./exec/defineGeneSets.R and stored in
#' ./data/GeneGroups.RData


#' Generate protein domain vector spaces:
dupl.w.orths.dom.vec.spaces <- lapply(dupl.w.orths, constructAnnotationVectorSpace)
tands.w.orths.dom.vec.spaces <- lapply(tands.w.orths, constructAnnotationVectorSpace)


#' Generate gene class specific vector clouds:
tands.w.orths.vec.clouds <- lapply(tands.w.orths.dom.vec.spaces, geneClassesVectorClouds, 
    gene.sets = tand.classifier)
dupl.w.orths.vec.clouds <- lapply(dupl.w.orths.dom.vec.spaces, geneClassesVectorClouds, 
    gene.sets = dupl.classifier)


#' Investigate number of annotation (domain) gains, losses, conserved (shared),
#' after duplication:
tands.w.orths.dom.evol <- Reduce(rbind, mclapply(tands.w.orths, annotationEvolution, 
    gene.sets = tand.classifier))
dupl.w.orths.dom.evol <- Reduce(rbind, mclapply(dupl.w.orths, annotationEvolution, 
    gene.sets = dupl.classifier))


#' Compare gene class specific vector clouds -
#' in order to investigate domain architecture evolution after duplication:
#' - For Tandems
t.c.i <- which(!as.logical(unlist(lapply(tands.w.orths.vec.clouds, function(x) any(is.na(x))))))
tands.w.orths.vec.clds.comp <- Reduce(rbind, lapply(tands.w.orths.vec.clouds[t.c.i], 
    function(t.c) {
        # Discard cloud comparisons that are meaningless, i.e. cause errors:
        suppressWarnings(tryCatch({
            cl.dist <- distVectorClouds(t.c[[1]], t.c[[2]])
            mean.orth.dom.vers <- (1 - cosDiag(t.c$ortholog$stat.vec)/sqrt(2))
            mean.tand.dom.vers <- (1 - cosDiag(t.c$tandem$stat.vec)/sqrt(2))
            diff.dom.vers <- mean.tand.dom.vers - mean.orth.dom.vers
            diff.dom.specif <- rad2deg(acos(cosAngleVec(t.c[[1]]$orth.on.diag.2.stat.vec, 
                t.c[[2]]$orth.on.diag.2.stat.vec)))
            data.frame(cloud.dist = cl.dist, mean.orth.dom.vers = mean.orth.dom.vers, 
                mean.tand.dom.vers = mean.tand.dom.vers, diff.dom.vers = diff.dom.vers, 
                diff.dom.specif = diff.dom.specif, stringsAsFactors = FALSE)
        }, error = function(e) {
        }))
    }))
#' - For Duplicated
d.c.i <- which(!as.logical(unlist(lapply(dupl.w.orths.vec.clouds, function(x) any(is.na(x))))))
dupl.w.orths.vec.clds.comp <- Reduce(rbind, lapply(dupl.w.orths.vec.clouds[d.c.i], 
    function(t.c) {
        # Discard cloud comparisons that are meaningless, i.e. cause errors:
        suppressWarnings(tryCatch({
            cl.dist <- distVectorClouds(t.c[[1]], t.c[[2]])
            mean.orth.dom.vers <- (1 - cosDiag(t.c$ortholog$stat.vec)/sqrt(2))
            mean.tand.dom.vers <- (1 - cosDiag(t.c$duplicated$stat.vec)/sqrt(2))
            diff.dom.vers <- mean.tand.dom.vers - mean.orth.dom.vers
            diff.dom.specif <- rad2deg(acos(cosAngleVec(t.c[[1]]$orth.on.diag.2.stat.vec, 
                t.c[[2]]$orth.on.diag.2.stat.vec)))
            data.frame(cloud.dist = cl.dist, mean.orth.dom.vers = mean.orth.dom.vers, 
                mean.tand.dom.vers = mean.tand.dom.vers, diff.dom.vers = diff.dom.vers, 
                diff.dom.specif = diff.dom.specif, stringsAsFactors = FALSE)
        }, error = function(e) {
        }))
    }))

message("DONE")
