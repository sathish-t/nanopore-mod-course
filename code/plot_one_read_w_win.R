#!/usr/bin/env Rscript

# Usage: Rscript ./plot_one_read_w_win.R plot_data plot.png

# Script: plot_one_read_w_win.R

# Input file: plot_data

# Input file description:
#  must have five columns: id, start, end, val, label.
#  columns must be separated by tabs.
#  must not have a header.
#  comment lines begin with '#'.
#  id is the readID. it must contain only one value since we are plotting only one read.  
#  start and end are positions on the reference genome.
#  val is the probability of the modification for raw or windowed data.
#  label is the the type of data. allowed labels are 'rawVal' for raw data and 'winVal' for windowed data.

# Ouput file: plot.png

# load R packages
require(ggplot2)
require(ggthemes)
library("dplyr")
options(bitmapType = "cairo")

# set dpi of output plot
dpi <- 50
width <- 18
height <- 12

# set default colour of the plot of windowed detect data
colour_win_detect <- "#C61F16"

# load command line arguments if available
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args) < 2) {
  stop("Usage: Rscript ./<script_name>.R plot_data plot.png", call.= FALSE)
}

# load data
read_data <- read.table(args[1], header = TRUE, comment.char = "#", sep = ",")
colnames(read_data) <- c("id", "start", "end", "val", "label")

# check that the id column has only one value
# if it has multiple values, then stop the script
if(length(unique(read_data$id)) > 1){
  stop("id column has more than one value! Script behaviour is undefined in this case", call.= FALSE)
}

# convert b to kb.
read_data$start <- read_data$start/1000
read_data$end <- read_data$end/1000

# set up plot choices
break_vector <- c('rawVal', 'winVal')
color_order_vector <- c('rawVal', 'winVal')
label_vector <- c("rawVal" = "Raw data", "winVal" = "Windowed data")
linetype_vector <- c("blank", "solid")
shape_vector <- c("circle", ".")

# mark which labels are present in input
indices <- vector()

if('rawVal' %in% read_data$label){
  indices <- append(indices,1)
}

if('winVal' %in% read_data$label){
  indices <- append(indices,2)
}

# plot the curve
plot1 <-  ggplot() +
      geom_point(data = subset(read_data, label == "rawVal"),
        aes(x = start, y = val, colour = label), shape='circle', alpha = 0.2, show.legend = TRUE) +
      geom_segment(data = subset(read_data, label == "winVal"),
        aes(x = start, y = val, xend = end, yend = val, colour = label), size = 2, show.legend = FALSE) +
      geom_step(data = subset(read_data, label == "winVal"),
        aes(x = start, y = val, colour = label), size = 2, show.legend = TRUE) +
      xlab("Coordinate (kb)") +
      ylab("Probability of base modification") +
      ylim(c(0, 1)) +
      scale_colour_manual(name = NULL,
                              values = c("rawVal" = "#888888", "winVal" = colour_win_detect),
			                  limits = color_order_vector,
                              breaks = break_vector[indices],
                              labels = label_vector[indices],
                              guide = guide_legend(override.aes = list(
                                                 linetype = linetype_vector[indices],
                                                 shape = shape_vector[indices]))
                              ) +
      theme_bw(base_size = 60) +
      theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

# save plot
ggsave(args[2], plot = plot1, dpi = dpi, width = width, height = height)