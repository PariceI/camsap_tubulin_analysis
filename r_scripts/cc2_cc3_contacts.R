library (ggplot2) # plotting
library (dplyr) # data manipulation

# This script plots efficiency of CC2-CC3 from different constructs across
# all three CAMSAPs

# DATA

data <- data.frame(
  Paralog = c("CAMSAP1", "CAMSAP1", "CAMSAP1", "CAMSAP1",
              "CAMSAP2", "CAMSAP2", "CAMSAP2", "CAMSAP2",
              "CAMSAP3", "CAMSAP3", "CAMSAP3", "CAMSAP3"),
  Construct = c("Full", "delCKK", "delCH", "Isolated",
                "Full", "delCKK", "delCH", "Isolated",
                "Full", "delCKK", "delCH", "Isolated"),
  Efficiency = c(0.88, 1.08, 1.37, 3.1,
                 0.64, 0.51, 0.5, 3.1,
                 1.3, 1.58, 2, 4.51))

# Factor levels so constructs appear in the correct order

data$Construct <- factor(data$Construct, levels = c("Full", "delCKK", "delCH", "Isolated"))

# Position dodge groups bars

plot <- ggplot(data, aes(x = Paralog, y = Efficiency, fill = Construct)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.8)) +
  
  # Value labels above bars and postion_dodge2 to ensure bars don't overlap
  
  geom_text(aes(label = format(Efficiency, big.mark = ",")),
            position = position_dodge2(width = 0.8, preserve = "single"),
            vjust = -0.3,
            size = 3.5,
            colour = "black") +
  
  # y axis scale
  
  scale_y_continuous(
    limits = c(0, 5),
    breaks = seq(0, 5, 0.5)
  ) +
  
  # Axis labels and title
  
  labs(
    title = "CC2-CC3 Linker Efficiency Across CAMSAP Paralogs",
    x = "",
    y = "Contacts Per Residue",
    fill = "Construct"
  ) +
  
  # COLOURING
  
  scale_fill_manual(values = c("Full" = 'lightblue',
                               "delCKK" = 'lightgreen',
                               "delCH" = 'pink',
                               "Isolated" = 'yellow')) +
  
  # THEME
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, colour = "black", face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 11, colour = "black"),
    legend.title = element_text(size = 11, colour = "black"),
    legend.text = element_text(size = 10, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# SAVE AND PRINT

ggsave("mbd_efficiency_barplot.png", plot, width = 8, height = 6, dpi = 300)
print(plot)

# MBD interface occupancy

library (ggplot2) # plotting
library (dplyr) # data manipulation

# This script plots occupancy of cc2-cc3 from different constructs across
# all three CAMSAPs

# DATA

data <- data.frame(
  Paralog = c("CAMSAP1", "CAMSAP1", "CAMSAP1", "CAMSAP1",
              "CAMSAP2", "CAMSAP2", "CAMSAP2", "CAMSAP2",
              "CAMSAP3", "CAMSAP3", "CAMSAP3", "CAMSAP3"),
  Construct = c("Full", "delCKK", "delCH", "Isolated",
                "Full", "delCKK", "delCH", "Isolated",
                "Full", "delCKK", "delCH", "Isolated"),
  Occupancy = c(8.26, 5.79, 7.02, 24.79,
                7.95, 6.28, 5.44, 24.27,
                18.1, 22.29, 22.89, 40.36)
) 

# Factor levels so constructs appear in the correct order

data$Construct <- factor(data$Construct, levels = c("Full", "delCKK", "delCH", "Isolated"))

# Position dodge groups bars

plot <- ggplot(data, aes(x = Paralog, y = Occupancy, fill = Construct)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.8)) +
  
  # Value labels above bars and postion_dodge2 to ensure bars don't overlap
  
  geom_text(aes(label = format(Occupancy, big.mark = ",")),
            position = position_dodge2(width = 0.8, preserve = "single"),
            vjust = -0.3,
            size = 3.5,
            colour = "black") +
  
  # y axis scale
  
  scale_y_continuous(
    limits = c(0, 42),
    breaks = seq(0, 42, 2)
  ) +
  
  # Axis labels and title
  
  labs(
    title = "CC2-CC3 Linker Occupancy Across CAMSAP Paralogs",
    x = "",
    y = "Occupancy (%)",
    fill = "Construct"
  ) +
  
  # COLOURING
  
  scale_fill_manual(values = c("Full" = 'lightblue',
                               "delCKK" = 'lightgreen',
                               "delCH" = 'pink',
                               "Isolated" = 'yellow')) +
  
  # THEME
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, colour = "black", face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 11, colour = "black"),
    legend.title = element_text(size = 11, colour = "black"),
    legend.text = element_text(size = 10, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# SAVE AND PRINT

ggsave("mbd_occupancy_barplot.png", plot, width = 8, height = 6, dpi = 300)
print(plot)