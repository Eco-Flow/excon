# Load necessary libraries
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
if (!requireNamespace("ggtree", quietly = TRUE)) {
  BiocManager::install("ggtree")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

library(ggtree)
library(ggplot2)
library(patchwork)
library(argparse)
library(dplyr)

# Parse command-line arguments
parser <- ArgumentParser(description = 'Plot phylogenetic tree with statistics and true/false data')
parser$add_argument('tree_file', type = 'character', help = 'Path to the Newick formatted tree file')
parser$add_argument('data_file', type = 'character', help = 'Path to the CSV file with statistic and true/false data')
parser$add_argument('phylogeny_width', type = 'integer', help = 'Width percentage of the phylogeny plot (0-100)')
args <- parser$parse_args()

# Read the Newick tree from the file
tree <- read.tree(args$tree_file)

# Read the data table from the file, ensuring species column is read as character
data <- read.csv(args$data_file, sep = "\t", colClasses = c("species" = "character"))

# Extract the column headers and the second row to determine plot types
column_headers <- colnames(data)
plot_types <- as.character(data[1, ])
data <- data[-1, ]  # Remove the second row used for plot types

# Debugging: Print species names from the tree and the data
cat("Species names in the tree:\n")
print(tree$tip.label)
cat("\nSpecies names in the data:\n")
print(data$species)

# Debugging: Print plot types
cat("\nPlot types:\n")
print(plot_types)

# Ensure species names match the tree tips
if (!all(data$species %in% tree$tip.label)) {
  stop("Species names in the data do not match the tree tips")
}

# Ensure the 'statistic' column is numeric
data$statistic <- as.numeric(data$statistic)

# Plot the phylogenetic tree
tree_plot <- ggtree(tree) + 
  geom_tiplab(size=3)  # Add labels to the tips of the tree

# Create plots based on the plot types
plots <- list()
for (i in 2:length(column_headers)) {
  column_name <- column_headers[i]
  plot_type <- plot_types[i]
  
  if (plot_type == "bar") {
    bar_plot <- ggplot(data, aes(y = reorder(species, !!sym(column_name)), x = !!sym(column_name))) + 
      geom_bar(stat = "identity", fill = "darkblue") + 
      theme_minimal() + 
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),  # Remove horizontal grid lines
            panel.grid.minor.y = element_blank()) + 
      labs(x = column_name)
    plots[[length(plots) + 1]] <- bar_plot
  } else if (plot_type == "text") {
    text_plot <- ggplot(data, aes(y = reorder(species, statistic), x = 1, label = !!sym(column_name))) + 
      geom_text(size = 3, hjust = 0) + 
      theme_void() + 
      labs(x = NULL, y = NULL)
    plots[[length(plots) + 1]] <- text_plot
  } else {
    cat("Unknown plot type:", plot_type, "for column:", column_name, "\n")
  }
}

# Debugging: Print the list of plots
cat("List of plots:\n")
print(plots)

# Calculate the widths for the plot layout
phylogeny_width <- args$phylogeny_width / 100
other_plots_width <- (1 - phylogeny_width) / length(plots)

# Combine the plots if there are any
if (length(plots) > 0) {
  combined_plot <- tree_plot + 
    wrap_plots(plots, ncol = length(plots)) + 
    plot_layout(ncol = length(plots) + 1, widths = c(phylogeny_width, rep(other_plots_width, length(plots))))  # Adjust widths based on input

  # Save the plot to a PDF file
  ggsave("Phyloplot.pdf", plot = combined_plot, width = 10, height = 8)
} else {
  cat("No plots to combine.\n")
}

# Print warnings
warnings()
