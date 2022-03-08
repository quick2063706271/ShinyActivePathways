MAX_PHENOTYPE_1 = rgb(178, 24, 43, max = 255)
LIGHTER_PHENOTYPE_1 = rgb(214, 96, 77, max = 255)
LIGHTEST_PHENOTYPE_1 = rgb(244, 165, 130, max = 255)
OVER_COLOR = rgb(247, 247, 247, max = 255)
MAX_PHENOTYPE_2 = rgb(33, 102, 172, max = 255)
LIGHTER_PHENOTYPE_2 = rgb(67, 147, 195, max = 255)
LIGHTEST_PHENOTYPE_2 = rgb(146, 197, 222, max = 255)



colors = c(OVER_COLOR, OVER_COLOR, LIGHTEST_PHENOTYPE_1, LIGHTER_PHENOTYPE_1, MAX_PHENOTYPE_1)
colorRanges <- c(0.0, 0.9, 0.95, 0.995, 1.0)
colorsPoints <- 1 - enrich$adjusted.p.val
enrich$colorsPoints <- 1 - enrich$adjusted.p.val

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

getRangeValue <- function(frac, lowerColor, upperColor) {
  lowerRange <- col2rgb(lowerColor)
  upperRange <- col2rgb(upperColor)
  red <- lowerRange[1] + frac * (upperRange[1] - lowerRange[1])
  green <- lowerRange[2] + frac * (upperRange[2] - lowerRange[2])
  blue <- lowerRange[3] + frac * (upperRange[3] - lowerRange[3])
  return(rgb(red, green, blue, max = 255))
}

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
# enrichColors <- c()
# for (i in seq_along(enrich$colorsPoints)) {
#   colorRange <- findColorRange(enrich$colorsPoints[i])
#   print(colorRange)
#   frac <- (enrich$colorsPoints[i] - colorRanges[colorRange])/(colorRanges[colorRange + 1] - colorRanges[colorRange])
#   print(enrich$term.name[i])
#   print(colors[colorRange])
#   cc <- getRangeValue(frac, colors[colorRange-1], colors[colorRange])
#   enrichColors <- c(enrichColors, cc)
# }

visg$x$nodes$color <- enrichColors
visg

# calculateColor <- function(colorRange, colorPoint) {
#   ranges <- colorRanges[colorRange: colorRange + 1]
#   
# # }
# axonColorPoints <- enrich[enrich$term.name=="Axon guidance", ]$colorsPoints
# # 0.6602298
# DescTools::MixColor(OVER_COLOR, LIGHTEST_PHENOTYPE_1, 0.6602298)




getRangeValue(0.6602298, OVER_COLOR, LIGHTEST_PHENOTYPE_1)

