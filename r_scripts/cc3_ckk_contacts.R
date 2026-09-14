library (ggplot2) # plotting
library (dplyr) # data manipulation

# This script plots efficiency and occupancy of the CC3-CKK region 
# from different constructs across all three CAMSAPs


# DATA

data <- data.frame(
  Paralog = c("CAMSAP1", "CAMSAP1", "CAMSAP1", "CAMSAP1", "CAMSAP1",
              "CAMSAP2", "CAMSAP2", "CAMSAP2", "CAMSAP2", "CAMSAP2",
              "CAMSAP3", "CAMSAP3", "CAMSAP3", "CAMSAP3", "CAMSAP3"),
  Construct = c("Full", "delCKK", "delCH", "delCC", "Isolated",
                "Full", "delCKK", "delCH", "delCC", "Isolated",
                "Full", "delCKK", "delCH", "delCC", "Isolated"),
  Efficiency = c(0.68, 0.84, 0.55, 0, 2.9,
                 0.14, 5.19, 0.64, 1, 2.65,
                 0.44, 5.72, 0.29, 0.35, 4.02))

# Factor levels so constructs appear in the correct order

data$Construct <- factor(data$Construct, levels = c("Full", "delCKK", "delCH", "delCC", "Isolated"))

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
    limits = c(0, 6),
    breaks = seq(0, 6, 0.5)
  ) +
  
  # Axis labels and title
  
  labs(
    title = "CC3-CKK Linker Efficiency Across CAMSAP Paralogs",
    x = "",
    y = "Contacts Per Residue",
    fill = "Construct"
  ) +
  
  # COLOURING
  
  scale_fill_manual(values = c("Full" = 'lightblue',
                               "delCKK" = 'lightgreen',
                               "delCH" = 'pink',
                               "delCC" = 'orange',
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

ggsave("d2_efficiency_barplot.png", plot, width = 8, height = 6, dpi = 300)
print(plot)

# CC3-CKK occupancy

library (ggplot2) # plotting
library (dplyr) # data manipulation

# This script plots occupancy of CC3-CKK from different constructs across
# all three CAMSAPs

# DATA

data <- data.frame(
  Paralog = c("CAMSAP1", "CAMSAP1", "CAMSAP1", "CAMSAP1", "CAMSAP1",
              "CAMSAP2", "CAMSAP2", "CAMSAP2", "CAMSAP2", "CAMSAP2",
              "CAMSAP3", "CAMSAP3", "CAMSAP3", "CAMSAP3", "CAMSAP3"),
  Construct = c("Full", "delCKK", "delCH", "delCC", "Isolated",
                "Full", "delCKK", "delCH", "delCC", "Isolated",
                "Full", "delCKK", "delCH", "delCC", "Isolated"),
  Occupancy = c(12.61, 11.76, 8.4, 0, 34.45,
                2.72, 40.91, 11.82, 7.27, 36.36,
                6.98, 20.93, 2.91, 8.14, 40.7)
) 

# Factor levels so constructs appear in the correct order

data$Construct <- factor(data$Construct, levels = c("Full", "delCKK", "delCH", "delCC", "Isolated"))

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
    title = "CC3-CKK Linker Occupancy Across CAMSAP Paralogs",
    x = "",
    y = "Occupancy (%)",
    fill = "Construct"
  ) +
  
  # COLOURING
  
  scale_fill_manual(values = c("Full" = 'lightblue',
                               "delCKK" = 'lightgreen',
                               "delCH" = 'pink',
                               "delCC" = 'orange',
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

ggsave("D2_occupancy_barplot.png", plot, width = 8, height = 6, dpi = 300)
print(plot)