library(tidyverse)
# Assuming 'df' is your data frame with columns 'Time', 'Value', 'Group', and 'LineID'
# Example data frame
df <- data.frame(
  Time = rep(c(0,0.5,1,1.5,2,4,6), times = 8),
  Value = c(0,1.8,2.8,3.2,3.0,1.7,1.6,0,0.3,1.5,2.3,2.8,2.1,1.5,
            0,1.0,1.4,2.1,2.5,2.9,2.85,0,0.2,0.4,0.65,1.0,2.3,2.85,
            0,-1.3,-1.8,-1.5,-1.4,-0.3,-0.1,0,0,-0.3,-0.8,-1.3,-0.6,-0.2,
            0,-0.8,-1.0,-1.15,-1.4,-1.55,-1.5,0,0,-0.1,-0.25,-0.4,-1.25,-1.6),
  Group = rep(c("Pattern1", "Pattern2","Pattern3","Pattern4"), each = 14),
  LineID = rep(c(rep(1,7),rep(2,7)), times = 4)
)


# Plotting the data with points and a legend, all plots in one row with their own y-axis
p1 <- ggplot(df %>% filter(Group == "Pattern1"), aes(x = Time, y = Value, color = as.factor(LineID))) +
  geom_line() +                                # Add lines
  geom_point(size = 3) +                       # Add points with a size of 3
  # facet_grid(. ~ Group) +   # Arrange plots in one row with free y-axis
  theme_classic() +                            # Use a classic theme
  labs(x = "Time (hours)", y = "Expression") +  # Add labels
  scale_color_manual(values = c("Blue", "Gray"),
                     name = "Data",            # Define the legend title
                     labels = c("Enhancer", "Gene")) +  # Define the legend labels
  theme(legend.position = "none")                            # Set the same y-axis limits for all facets

p2 <- ggplot(df %>% filter(Group == "Pattern2"), aes(x = Time, y = Value, color = as.factor(LineID))) +
  geom_line() +                                # Add lines
  geom_point(size = 3) +                       # Add points with a size of 3
  # facet_grid(. ~ Group) +   # Arrange plots in one row with free y-axis
  theme_classic() +                            # Use a classic theme
  labs(x = "Time (hours)", y = "Expression") +  # Add labels
  scale_color_manual(values = c("Blue", "Gray"),
                     name = "Data",            # Define the legend title
                     labels = c("Enhancer", "Gene")) +  # Define the legend labels
  theme(legend.position = "none")                            # Set the same y-axis limits for all facets

p3 <- ggplot(df %>% filter(Group == "Pattern3"), aes(x = Time, y = Value, color = as.factor(LineID))) +
  geom_line() +                                # Add lines
  geom_point(size = 3) +                       # Add points with a size of 3
  # facet_grid(. ~ Group) +   # Arrange plots in one row with free y-axis
  theme_classic() +                            # Use a classic theme
  labs(x = "Time (hours)", y = "Expression") +  # Add labels
  scale_color_manual(values = c("Blue", "Gray"),
                     name = "Data",            # Define the legend title
                     labels = c("Enhancer", "Gene")) +  # Define the legend labels
  theme(legend.position = "none")                            # Set the same y-axis limits for all facets

p4 <- ggplot(df %>% filter(Group == "Pattern4"), aes(x = Time, y = Value, color = as.factor(LineID))) +
  geom_line() +                                # Add lines
  geom_point(size = 3) +                       # Add points with a size of 3
  # facet_grid(. ~ Group) +   # Arrange plots in one row with free y-axis
  theme_classic() +                            # Use a classic theme
  labs(x = "Time (hours)", y = "Expression") +  # Add labels
  scale_color_manual(values = c("Blue", "Gray"),
                     name = "Data",            # Define the legend title
                     labels = c("Enhancer", "Gene")) +  # Define the legend labels
  theme(legend.position = "none")                            # Set the same y-axis limits for all facets
combined_plot <- p1 + p2 + p3+p4+ plot_layout(ncol = 4)
combined_plot
jpeg(file="/Users/wancen/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/GRA/lima/plots/introduction.jpg",width = 15, height = 3.5,units = "in",res=450)
p1 + p2 + p3+p4+ plot_layout(ncol = 4)&
  theme(text = element_text(size = 22),
        plot.tag = element_text(size = 16)
  )
dev.off()
