#Packages
library(tidyverse)
library(reshape2)
library(fGarch)
library(rugarch)
library(depmixS4)
library(HiddenMarkov)

#color-blinded friendly
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
# scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
# scale_colour_manual(values=cbPalette)

# plotl <- function(x, ...){
#   #naive lines plot
#   #...: pass arguments to the function
#   plot(x, type = "l", ...)
# }

plotl <- function(x, main = NULL){
  if (is.null(main)) {
    main1 <- deparse(substitute(x))
    main <- paste("Line plot of ", main1)
  }
  #quickly draw a line plot 
  plot(x, type = "l")
}


qqnorm2 <- function(x){
  #quickly draw a qqnorm plot 
  #or use qqnormPlot
  qqnorm(x)
  qqline(x, col = 2)
}
