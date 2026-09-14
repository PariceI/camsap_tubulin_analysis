library (ggplot2) # plotting
library (dplyr) # data manipulation

# This script plots BSA of CC2-CC3 from different constructs across
# all three CAMSAPs

# DATA

data <- data.frame(
  Paralog = c("CAMSAP1", "CAMSAP1", "CAMSAP1", "CAMSAP1",
              "CAMSAP2", "CAMSAP2", "CAMSAP2", "CAMSAP2",
              "CAMSAP3", "CAMSAP3", "CAMSAP3", "CAMSAP3"),
  Construct = c("Full", "delCKK", "delCH", "Isolated",
                "Full", "delCKK", "delCH", "Isolated",
                "Full", "delCKK", "delCH", "Isolated"),
  BSA = c(1326, 1126, 1435, 4747,
          1498, 979, 917, 4597,
          2126, 2457, 2584, 5849)
)

# Factor levels so constructs appear in the correct order

data$Construct <- factor(data$Construct, levels = c("Full", "delCKK", "delCH", "Isolated"))

# Position dodge groups bars

plot <- ggplot(data, aes(x = Paralog, y = BSA, fill = Construct)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.8)) +
  
  # Value labels above bars and postion_dodge2 to ensure bars don't overlap
  
  geom_text(aes(label = format(BSA, big.mark = ",")),
            position = position_dodge2(width = 0.8, preserve = "single"),
            vjust = -0.3,
            size = 3.5,
            colour = "black") +
  
  # Axis labels and title
  
  labs(
    title = "CC2_CC3_Linker Region: Buried Surface Area by Construct",
    x = "",
    y = "Buried Surface Area (Å)",
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

ggsave("mbd_bsa_barplot.png", plot, width = 8, height = 6, dpi = 300)
print(plot)