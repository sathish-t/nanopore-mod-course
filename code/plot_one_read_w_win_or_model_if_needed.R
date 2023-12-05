#!/usr/bin/env Rscript
# usage Rscript ./program_name.R plot_data plot.png plot_elements_this_is_optional

# Input file description
#
# * input file plot_data must have five columns: id start end val label
# * No header and columns separated by spaces.
# * comment lines begin with #
# * allowed labels are 'rawDetect', 'winDetect', 'model'. any of them can be missing.
# * id column must contain only one value since we are plotting only one read.
# * for rawDetect and model, only the start column is used in plotting, whereas
#   for winDetect both are used.
# * input file plot_elements_this_is_optional is optional. it has three necessary columns and
#   one optional column, separated by spaces and no header.
#   the necessary columns are start, end, label; the optional column is line style.
#   start and end are positions on the reference genome.
#   start must always be less than or equal to end i.e. don't switch them for left/right fork.
#   label can be 'origin', 'leftFork', 'rightFork', 'termination', or 'pause'.
#   line style is a (binary) column for line style in plot for a fork feature.
#   set it to zero or one if you want to distinguish between fork features.

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

# flag to check if the colour of the windowed detect data has been set
is_colour_set <- FALSE

# flag to check if annotations have been set in the plot_data file
is_annotation_set_in_plot_data <- FALSE

# load command line arguments if available
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args) < 2) {
  stop("Usage: Rscript ./<script_name>.R plot_data plot.png plot_elements_this_is_optional", call.= FALSE)
}

# process the annotation file if it exists
if (length(args) == 3 && file.size(args[3]) > 0) {

  # see if the annotation file has some data in it
  n_comment_lines <- 0
  n_total_lines <- 0

  for(line in readLines(args[3])){
    n_comment_lines <- n_comment_lines + (substr(line, 0, 1) == "#")
    n_total_lines <- n_total_lines + 1
  }

  if(n_comment_lines < n_total_lines){
    annotations <- read.table(args[3], header = FALSE, comment.char = "#")

    if(ncol(annotations) == 3){
      colnames(annotations) <- c("start", "end", "label")
      annotations$size <- 1
    } else if(ncol(annotations) == 4){
      colnames(annotations) <- c("start", "end", "label", "size")
    }  else {
      stop("Usage: plot_elements_this_is_optional must have three or four columns", call.= FALSE)
    }

    annotations$start <- annotations$start/1000
    annotations$end <- annotations$end/1000

    annotations$y_start <- 1.00
    annotations$y_end <- 1.00

    annotations[annotations$label == "origin", "y_start"] <- 1.00
    annotations[annotations$label == "origin", "y_end"] <- 1.00

    annotations[annotations$label == "rightFork", "y_start"] <- 0.99
    annotations[annotations$label == "rightFork", "y_end"] <- 0.98

    annotations[annotations$label == "leftFork", "y_start"] <- 0.98
    annotations[annotations$label == "leftFork", "y_end"] <- 0.99

    annotations[annotations$label == "termination", "y_start"] <- 0.97
    annotations[annotations$label == "termination", "y_end"] <- 0.97

    annotations[annotations$label == "pause", "y_start"] <- 1.00
    annotations[annotations$label == "pause", "y_end"] <- 0.97
  }

}

# load data
read_data <- read.table(args[1], header = FALSE, comment.char = "#")
colnames(read_data) <- c("id", "start", "end", "val", "label")

# check that the id column has only one value
# if it has multiple values, then stop the script
if(length(unique(read_data$id)) > 1){
  stop("id column has more than one value! Script behaviour is undefined in this case", call.= FALSE)
}

# read comments from input file
comments <- readLines(args[1], n = 1000)
comments <- comments[substr(comments, 0, 1) == "#"]

# remove leading and trailing spaces from comments
comments <- gsub("^\\s+|\\s+$", "", comments)

# if any entry in comments is "orientation: +", then set the colour of the windowed detect data to blue
if(any(grepl("orientation: \\+", comments))){
  colour_win_detect <- "#1F68C4"
  is_colour_set <- TRUE
}

# if any entry in comments is "orientation: -", then set the colour of the windowed detect data to red
# if the colour has already been set, then stop the script
if(any(grepl("orientation: -", comments))){
  colour_win_detect <- "#C61F16"
  if(is_colour_set){
    stop("orientation: + and orientation: - both present in comments. Script behaviour is undefined in this case",
         call.= FALSE)
  }
  is_colour_set <- TRUE
}


# read all entries in comments that start with "annotate BED6+1: " or "annotate BED3+1: "
annotate_bed6plus1_comments <- comments[grepl("^# annotate BED6\\+1: ", comments)]
annotate_bed3plus1_comments <- comments[grepl("^# annotate BED3\\+1: ", comments)]

# remove the "annotate BED6+1: " or "annotate BED3+1: " part from the annotation comments
annotate_bed6plus1_comments <- gsub("^# annotate BED6\\+1: ", "", annotate_bed6plus1_comments)
annotate_bed3plus1_comments <- gsub("^# annotate BED3\\+1: ", "", annotate_bed3plus1_comments)

# split the entries in annotate_bed6plus1_comments on tab and store in a data frame with appropriate column names
annotate_bed6plus1_comments_data <- data.frame(do.call("rbind", strsplit(annotate_bed6plus1_comments, "\\t")),
                                               stringsAsFactors = FALSE, fix.empty.names = TRUE)
if(nrow(annotate_bed6plus1_comments_data) > 0){
  colnames(annotate_bed6plus1_comments_data) <- c("contig", "start", "end", "name", "score", "strand", "label")
  # if any label starts with "plot_", then set the is_annotation_set_in_plot_data flag to TRUE
    if(any(grepl("^plot_", annotate_bed6plus1_comments_data$label))){
        is_annotation_set_in_plot_data <- TRUE
    }
}

# split the entries in annotate_bed3plus1_comments on tab and store in a data frame with appropriate column names
annotate_bed3plus1_comments_data <- data.frame(do.call("rbind", strsplit(annotate_bed3plus1_comments, "\\t")),
                                               stringsAsFactors = FALSE, fix.empty.names = TRUE)
if(nrow(annotate_bed3plus1_comments_data) > 0){
  colnames(annotate_bed3plus1_comments_data) <- c("contig", "start", "end", "label")
  # if any label starts with "plot_", then set the is_annotation_set_in_plot_data flag to TRUE
  if(any(grepl("^plot_", annotate_bed3plus1_comments_data$label))){
      is_annotation_set_in_plot_data <- TRUE
  }
}

# if any annotation data was read, then remove any rows that don't have the string "plot_" in the label column
# and then bind the data frames together, convert positions to kb, set line type, and assign y coordinates
if(is_annotation_set_in_plot_data){

  annotate_comments_data <- bind_rows(annotate_bed6plus1_comments_data, annotate_bed3plus1_comments_data)
  annotate_comments_data <- annotate_comments_data[grepl("^plot_", annotate_comments_data$label),]

  if(nrow(annotate_comments_data) > 0){

    # check that the 'strand' column exists, otherwise create it and set it to "."
    if(!"strand" %in% colnames(annotate_comments_data)){
      annotate_comments_data$strand <- "."
    }

    # if strand is NA, then set it to "."
    annotate_comments_data$strand[is.na(annotate_comments_data$strand)] <- "."

    # convert positions to kb
    annotate_comments_data$start <- as.integer(annotate_comments_data$start)/1000
    annotate_comments_data$end <- as.integer(annotate_comments_data$end)/1000

    annotate_comments_data$y_start <- 1.02
    annotate_comments_data$y_end <- 1.02

    # set both y_start and y_end to 1.03 if the strand is "+"
    annotate_comments_data$y_start[annotate_comments_data$strand == "+"] <- 1.03
    annotate_comments_data$y_end[annotate_comments_data$strand == "+"] <- 1.03

    # create a lineend column that is set to 'butt' by default
    annotate_comments_data$lineend <- "butt"

    # create a colour column that is set to 'black' by default,
    # "#C61F16" if the strand is "-", and "#1F68C4" if the strand is "+"
    annotate_comments_data$colour <- "black"
    annotate_comments_data$colour[annotate_comments_data$strand == "-"] <- "#C61F16"
    annotate_comments_data$colour[annotate_comments_data$strand == "+"] <- "#1F68C4"

  } else {
    is_annotation_set_in_plot_data <- FALSE
  }
}

# convert b to kb.
read_data$start <- read_data$start/1000
read_data$end <- read_data$end/1000

# set up plot choices
break_vector <- c('rawDetect', 'winDetect', 'model')
color_order_vector <- c('model', 'rawDetect', 'winDetect')
label_vector <- c("rawDetect" = "Raw data", "winDetect" = "Windowed data", "model" = "Best-fit model")
linetype_vector <- c("blank", "solid", "solid")
shape_vector <- c("circle", ".", ".")

# mark which labels are present in input
indices <- vector()

if('rawDetect' %in% read_data$label){
  indices <- append(indices,1)
}

if('winDetect' %in% read_data$label){
  indices <- append(indices,2)
}

if('model' %in% read_data$label){
  indices <- append(indices,3)
}

# plot the curve
plot1 <-  ggplot() +
      geom_point(data = subset(read_data, label == "rawDetect"),
        aes(x = start, y = val, colour = label), shape='circle', alpha = 0.2, show.legend = TRUE) +
      geom_segment(data = subset(read_data, label == "winDetect"),
        aes(x = start, y = val, xend = end, yend = val, colour = label), size = 2, show.legend = FALSE) +
      geom_step(data = subset(read_data, label == "winDetect"),
        aes(x = start, y = val, colour = label), size = 2, show.legend = TRUE) +
      geom_path(data = subset(read_data, label == "model"),
        aes(x = start, y = val, colour = label), size = 2, show.legend = TRUE) +
      xlab("Reference coordinate (kb)") +
      ylab("Probability of BrdU") +
      ylim(c(0, 1)) +
      scale_colour_manual(name = NULL,
                              values = c("model" = "#000000", "rawDetect" = "#888888", "winDetect" = colour_win_detect),
			                  limits = color_order_vector,
                              breaks = break_vector[indices],
                              labels = label_vector[indices],
                              guide = guide_legend(override.aes = list(
                                                 linetype = linetype_vector[indices],
                                                 shape = shape_vector[indices]))
                              ) +
      theme_bw(base_size = 60) +
      theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank())

if (length(args) == 3 && exists('annotations')) {
  plot1 <- plot1 +
    geom_segment(
      aes(x = as.numeric(start), y = as.numeric(y_start),
          xend = as.numeric(end), yend = as.numeric(y_end), linetype = as.factor(size)),
      data = annotations,
      colour = "olivedrab4", alpha = 1, size = 1
    )
}

if(is_annotation_set_in_plot_data){

  # get min value of read_data$start and max value of read_data$end, and go 1kb either side
  min_start <- min(read_data$start) - 1
  max_end <- max(read_data$end) + 1

  # if min_start is less than 0, then set it to 0
  if(min_start < 0){
    min_start <- 0
  }

  # plot annotations
  plot1 <- plot1 +
    geom_segment(
      aes(x = as.numeric(start), y = as.numeric(y_start),
          xend = as.numeric(end), yend = as.numeric(y_end)),
      data = annotate_comments_data,
      lineend = annotate_comments_data$lineend,
      colour = annotate_comments_data$colour, alpha = 0.4, size = 4
    ) + ylim(c(0, 1.04)) + coord_cartesian(xlim=c(min_start, max_end))
}

ggsave(args[2], plot = plot1, dpi = dpi, width = width, height = height)
