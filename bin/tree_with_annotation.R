# Load necessary libraries
library(ggtree)
library(ggplot2)
library(patchwork)  # For combining plots

# Sample tree data (replace this with your actual tree)
tree <- rtree(10)  # Generates a random tree with 10 species

# Sample statistics data (replace this with your actual data)
data <- data.frame(
  species = tree$tip.label,  # Ensure species names match the tree tips
  statistic = runif(10, 0, 1)  # Generate random statistics
)

# Define a threshold for True/False
threshold <- 0.5

# Add a True/False column based on the threshold
data$above_threshold <- data$statistic > threshold

# Plot the phylogenetic tree
tree_plot <- ggtree(tree) + 
  geom_tiplab(size=3)  # Add labels to the tips of the tree

# Create a horizontal bar plot
bar_plot <- ggplot(data, aes(y = reorder(species, statistic), x = statistic)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  labs(x = "Statistic")

# Create a plot for the True/False values
true_false_plot <- ggplot(data, aes(y = reorder(species, statistic), x = 1, label = above_threshold)) + 
  geom_text(size = 3, hjust = 0) + 
  theme_void() + 
  labs(x = NULL, y = NULL)

# Combine the plots
combined_plot <- tree_plot + 
  (bar_plot | true_false_plot) + 
  plot_layout(ncol = 2, widths = c(3, 1))  # Adjust widths to your preference

# Display the plot
print(combined_plot)
