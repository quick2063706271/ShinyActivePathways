library(grDevices)

# Set color schemes
MAX_PHENOTYPE_1 = rgb(178, 24, 43, max = 255)
LIGHTER_PHENOTYPE_1 = rgb(214, 96, 77, max = 255)
LIGHTEST_PHENOTYPE_1 = rgb(244, 165, 130, max = 255)
OVER_COLOR = rgb(247, 247, 247, max = 255)
MAX_PHENOTYPE_2 = rgb(33, 102, 172, max = 255)
LIGHTER_PHENOTYPE_2 = rgb(67, 147, 195, max = 255)
LIGHTEST_PHENOTYPE_2 = rgb(146, 197, 222, max = 255)

colors = c(OVER_COLOR, OVER_COLOR, LIGHTEST_PHENOTYPE_1, LIGHTER_PHENOTYPE_1, MAX_PHENOTYPE_1)
colorRanges <- c(0.0, 0.9, 0.95, 0.995, 1.0)

# findColorRange
#' Find the color range of color point value
#'
#' A function that returns a vector containing the color range of color point values
#'
#' @param colorPoint A numeric value ranging from 0 to 1
#' @return A vector containing the color range of color point values
#'
#' @examples
#' findColorRange(0.2)
#'
#' @export
findColorRange <- function(colorPoint) {
  colorRanges <- c(0.0, 0.9, 0.95, 0.995, 1.0)
  i <- 1
  while (i <= length(colorRanges)) {
    if (colorPoint <= colorRanges[i]) {
      return(i - 1)
    }
    i = i + 1
  }
}
# getRangeValue
#' Find the color interpolation of two colors
#'
#' A function that returns a new color interpolated from lowerColor and upperColor
#'
#' @param frac A numeric value ranging from 0 to 1
#' @param lowerColor A hex color value
#' @param upperColor A hex color value
#' @return A hex color value
#'
#' @examples
#' getRangeValue(0.5, rgb(178, 24, 43, max = 255), rgb(214, 96, 77, max = 255))
#'
#' @export
getRangeValue <- function(frac, lowerColor, upperColor) {
  lowerRange <- col2rgb(lowerColor)
  upperRange <- col2rgb(upperColor)
  red <- lowerRange[1] + frac * (upperRange[1] - lowerRange[1])
  green <- lowerRange[2] + frac * (upperRange[2] - lowerRange[2])
  blue <- lowerRange[3] + frac * (upperRange[3] - lowerRange[3])
  return(rgb(red, green, blue, max = 255))
}

# getColors
#' Compute colors for an enrichment result
#'
#' A function that returns a vector of colors for an enrichment result
#'
#' @param enrich The output of enrichment map
#' @return A vector of colors for an enrichment result
#'
#' @examples
#' library(ActivePathways)
#' scores <- read.table(
#' system.file('extdata', 'Adenocarcinoma_scores_subset.tsv', package = 'ActivePathways'),
#' header = TRUE, sep = '\t', row.names = 1)
#' scores <- as.matrix(scores)
#' scores[is.na(scores)] <- 1
#' gmt.file <- system.file('extdata', 'hsapiens_REAC_subset.gmt', package = 'ActivePathways')
#' enrich <- ActivePathways(scores, gmt.file,cutoff = 0.1, significant = 0.05, merge.method = "Brown")
#' getColors(enrich)
#'
#' @export
getColors <- function(enrich) {
  colorsPoints <- 1 - enrich$adjusted.p.val
  enrich$colorsPoints <- 1 - enrich$adjusted.p.val
  enrichColors <- c()
  for (i in seq_along(enrich$colorsPoints)) {
    colorRange <- findColorRange(enrich$colorsPoints[i])
    print(colorRange)
    frac <- (enrich$colorsPoints[i] - colorRanges[colorRange])/(colorRanges[colorRange + 1] - colorRanges[colorRange])
    print(enrich$term.name[i])
    print(colors[colorRange])
    cc <- getRangeValue(frac, colors[colorRange-1], colors[colorRange])
    enrichColors <- c(enrichColors, cc)
  }
  return(enrichColors)
}