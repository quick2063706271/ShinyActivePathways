library(igraph)
library(ActivePathways)
library(data.table)

# getDataFromGMT
#' Create a list containing id, name and length of the genset from gmt file
#'
#' A function that creates a list containing id, name and length of the genset 
#' from gmt file
#'
#' @param gmt A list object reading from gmt file
#' tree
#'
#' @return Returns a list object 
#'
#' @examples
#' # Example 1
#' library(ActivePathways)
#' gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
#' getDataFromGMT(gmt.file)
#'
#' @export

getDataFromGMT <- function(gmt) {
  lens <- numeric(length = length(names(gmt)))
  ids <- character(length = length(names(gmt)))
  pathwayNames <- character(length = length(names(gmt)))
  # add data to different vectors
  for (i in seq_along(names(gmt))) {
    name <- names(gmt)[i]
    lens[i] <- length(gmt[[name]]$genes)
    ids[i] <- gmt[[name]]$id
    pathwayNames[i] <- gmt[[name]]$name
  }
  return(list(ids = ids, 
              pathwayNames = pathwayNames, 
              lengths = lens))
}


# computeSimilarityCoefficient
#' Compute the similarity coefficient from two sets
#'
#' A function that computes the similarity coefficient from two sets providing 
#' algorithm
#'
#' @param metric A string representing metric for computing similarity coefficient
#' @param intersectionSet The intersection of genes1 and genes2
#' @param unionSet The union of genes1 and genes2
#' @param genes1 One of the sets
#' @param genes2 The other one of the sets
#' @param k When metric is neither jaccard nor overlap, using combined method with
#' coefficient k
#'
#' @return Returns a similarity score
#'
#' @examples
#' # Example 1
#' set1 <- c("a", "b", "c")
#' set2 <- c("a", "b", "d")
#' computeSimilarityCoefficient("jaccard", intersect(set1, set2), 
#'  union(set1, set2), 
#'  set1, 
#'  set2)
#'
#' @export

computeSimilarityCoefficient <- function(metric, intersectionSet, unionSet, genes1, genes2, k) {
  if (metric == "jaccard") {
    return(length(intersectionSet) / length(unionSet))
  } else if (metric == "overlap") {
    return(length(intersectionSet) / min(length(unionSet), length(genes2)))
  } else {
    intersectionLength <- length(intersectionSet)
    jaccard <- intersectionLength / length(unionSet)
    overlap <- intersectionLength / min(length(unionSet), length(genes2))
    return(k * overlap) + ((1 - k) * jaccard)
  }
}

# computeSCMatrix
#' Compute the similarity coefficient between ids
#'
#' A function that computes the similarity coefficient between each geneset 
#' selected from ids
#'
#' @param algorithm A string representing metric for computing similarity coefficient
#' @param gmt The list represents gmt 
#' @param ids The ids of genesets for which you want to compute similarity coefficient
#' @param k When metric is neither jaccard nor overlap, using combined method with
#' coefficient k
#'
#' @return Returns a matrix of similarity score
#'
#' @examples
#' # Example 1
#' library(ActivePathways)
#' gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
#' gmt <- read.GMT(gmt.file)
#' enrich <- ActivePathways(scores, gmt.file,cutoff = 0.1, significant = 0.05, merge.method = "Brown")
#' pathways <- enrich[, c("term.id", "term.name", "adjusted.p.val")]
#' matchesColumns <- match(names(gmt), pathways$term.id) > 0
#' matchesColumns[is.na(matchesColumns)] = FALSE 
#' gmt <- gmt[matchesColumns]
#' gmtData <- getDataFromGMT(gmt)
#' scMatrix <- computeSCMatrix(algorithm = "jaccard", 
#'                             gmt = gmt, 
#'                             ids = gmtData$ids, 
#'                             k = 0.5)
#'
#' @export

computeSCMatrix <- function(algorithm, gmt, ids, k) {
  simiScoreMatrix <- matrix(ncol = length(ids), 
                            nrow = length(ids), 
                            dimnames = list(ids, 
                                            ids)
  )
  # initiate Matrix storing similarity score between gene sets
  simiScoreMatrix <- matrix(ncol = length(ids), 
                            nrow = length(ids), 
                            dimnames = list(ids, 
                                            ids)
  )
  
  # calculate similarity score and store them in the gene sets
  for (i in seq_along(ids)) {
    for (j in seq_along(ids)) {
      genes1 <- gmt[[ids[i]]]$genes
      genes2 <- gmt[[ids[j]]]$genes
      simiScoreMatrix[ids[i], ids[j]] <- computeSimilarityCoefficient(algorithm, 
                                                                      intersect(genes1, genes2), 
                                                                      union(genes1, genes2), 
                                                                      genes1, 
                                                                      genes2, 
                                                                      k = k)
    }
  }
  print("simiScoreMatrix1")
  # print(simiScoreMatrix)
  return(simiScoreMatrix)
}

# plotSimiScoreMatrix
#' Plot the SimiScoreMatrix as an enrichmentmap igraph object
#'
#' A function that plots the SimiScoreMatrix as an enrichmentmap igraph object
#'
#' @param simiScoreMatrix A matrix containing similarity coefficient
#' @param similarityCutoff A numeric value. If a similarity coefficient between 2 genesets is
#'  less than this value, this similarity is ignored (set to 0)
#' @param pvalueCutoff A numeric value. If a pvalue of a geneset is
#'  less than this value, this geneset is removed
#' @param lens An integer vector representing the sizes of genesets
#' @param pathwayNames A character vector representing the names of genesets
#' @param pathways A data.frame that holds information about genesets from
#' ActivePathways result
#'
#' @return Plot a matrix of similarity score
#'
#' @examples
#' # Example 1
#' library(ActivePathways)
#' gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
#' gmt <- read.GMT(gmt.file)
#' enrich <- ActivePathways(scores, gmt.file,cutoff = 0.1, significant = 0.05, merge.method = "Brown")
#' pathways <- enrich[, c("term.id", "term.name", "adjusted.p.val")]
#' matchesColumns <- match(names(gmt), pathways$term.id) > 0
#' matchesColumns[is.na(matchesColumns)] = FALSE 
#' gmt <- gmt[matchesColumns]
#' gmtData <- getDataFromGMT(gmt)
#' scMatrix <- computeSCMatrix(algorithm = "jaccard", 
#'                             gmt = gmt, 
#'                             ids = gmtData$ids, 
#'                             k = 0.5)
#' plotSimiScoreMatrix(scMatrix, 1, 0.05, gmtData$lengths, gmtData$pathwayNames, pathways)
#'
#' @export
plotSimiScoreMatrix <- function(simiScoreMatrix, similarityCutoff, pvalueCutoff, lens, pathwayNames, pathways) {
  # remove similaity Score less than the cutoff
  simiScoreMatrix[simiScoreMatrix < similarityCutoff] = 0
  # remove self-loop
  for (i in seq_along(rownames(simiScoreMatrix))) {
    print(simiScoreMatrix[i, i])
    print(rownames(simiScoreMatrix)[i])
    simiScoreMatrix[i, i] = 0
  }
  # simiScoreMatrix[simiScoreMatrix == 1] = 0
  print("simiScoreMatrix2")
  # print(simiScoreMatrix)
  
  # make colnames, rownames to pathway name instead of id
  colnames(simiScoreMatrix) <- pathwayNames
  rownames(simiScoreMatrix) <- pathwayNames
  # Create igraph object
  mode(simiScoreMatrix) <- "numeric"
  g2 <- igraph::graph.adjacency(simiScoreMatrix, mode = "undirected", weighted = TRUE)
  
  # set vertex size as the size of gene sets
  V(g2)$size <- lens * 0.04 + 4
  # set edge thickness (width) as the weight of edge (similarity score)
  E(g2)$width <- E(g2)$weight * 2
  
  # color
  # get a vector of pathway names ordered by p.val
  orderedPathways <- pathways[order(adjusted.p.val)]$term.name
  # get position of vertex in ordered Pathways
  position <- match(V(g2)$name, orderedPathways)
  # get color based on position
  heatColors <- heat.colors(length(V(g2)))[position]
  V(g2)$color <- heatColors
  
  # apply force-directed layout algorithm by Fruchterman and Reingold
  layout_with_fr(g2) 
  # plot
  return(g2)
} 


# plotEnrichmentMap
#' Plot enrichment map from result of ActivePathways
#'
#' A function that plots enrichment map from result of ActivePathways
#'
#' @param gmt The list represents gmt 
#' @param enrichment The output of ActivePathways
#' @param algorithm A string representing metric for computing similarity coefficient
#' @param similarityCutoff A numeric value. If a similarity coefficient between 2 genesets is
#'  less than this value, this similarity is ignored (set to 0)
#' @param pvalueCutoff A numeric value. If a pvalue of a geneset is
#'  less than this value, this geneset is removed
#' When metric is neither jaccard nor overlap, using combined method with
#' coefficient k
#'
#' @return An igraph object that plots enrichment map from result of ActivePathways
#'
#' @examples
#' scores <- read.table(
#'   system.file('extdata', 'Adenocarcinoma_scores_subset.tsv', package = 'ActivePathways'),
#'  header = TRUE, sep = '\t', row.names = 'Gene')
#' scores <- as.matrix(scores)
#' scores[is.na(scores)] <- 1
#' gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
#' enrich <- ActivePathways(scores, gmt.file)
#' gmt <- read.GMT(gmt.file)
#' g <- plotEnrichmentMap(gmt = gmt, 
#'                   enrichment = enrich, 
#'                   algorithm = "jaccard", 
#'                   similarityCutoff = 0.25, 
#'                   pvalueCutoff = NULL, 
#'                   k = 0.5)
#' plot(g, vertex.label.cex = 0.6)
#' @export
#' 
plotEnrichmentMap <- function(gmt, enrichment, algorithm, similarityCutoff, pvalueCutoff, k) {
  pathways <- enrichment[, c("term.id", "term.name", "adjusted.p.val")]
  matchesColumns <- match(names(gmt), pathways$term.id) > 0
  matchesColumns[is.na(matchesColumns)] = FALSE 
  gmt <- gmt[matchesColumns]
  gmtData <- getDataFromGMT(gmt)
  scMatrix <- computeSCMatrix(algorithm = algorithm, 
                              gmt = gmt, 
                              ids = gmtData$ids, 
                              k = k)
  enrichPlot <- plotSimiScoreMatrix(simiScoreMatrix = scMatrix, 
                                    similarityCutoff = similarityCutoff, 
                                    pvalueCutoff = pvalueCutoff, 
                                    lens = gmtData$lengths, 
                                    pathwayNames = gmtData$pathwayNames,
                                    pathways = pathways)
  return(enrichPlot)
}

