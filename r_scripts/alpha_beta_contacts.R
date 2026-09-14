library (ggplot2)
library (dplyr)
library (tidyr)

# MBD DATA

mbd_data <- data.frame(
  Paralog = c(rep("C1", 4), rep("C2", 4), rep("C3", 4)),
  Construct = c("Full", "delCKK", "delCH", "Isolated",
                "Full", "delCKK", "delCH", "Isolated",
                "Full", "delCKK", "delCH", "Isolated"),
  Length = c(242, 242, 242, 242,
             239, 239, 239, 239,
             166, 166, 166, 166),
  Contacts = c(151, 132, 169, 551,
               153, 122, 120, 742,
               215, 262, 332, 749),
  Alpha = c(151, 132, 169, 229,
            108, 111, 113, 178,
            210, 231, 283, 418),
  Beta = c(0, 0, 0, 332,
           45, 11, 7, 564,
           5, 31, 49, 331),
  Interface_Residues = c(20, 14, 17, 60,
                         19, 15, 13, 58,
                         30, 37, 38, 67)
)

# D2 DATA

d2_data <- data.frame(
  Paralog = c(rep("C1", 5), rep("C2", 5), rep("C3", 5)),
  Construct = c("Full", "delCKK", "delCH", "delCC", "Isolated",
                "Full", "delCKK", "delCH", "delCC", "Isolated",
                "Full", "delCKK", "delCH", "delCC", "Isolated"),
  Length = c(119, 119, 119, 119, 119,
             110, 110, 110, 110, 110,
             172, 172, 172, 172, 172),
  Contacts = c(81, 100, 65, 0, 350,
               16, 571, 70, 110, 292,
               76, 984, 50, 61, 691),
  Alpha = c(50, 68, 25, 0, 206,
            9, 224, 53, 19, 112,
            63, 775, 23, 0, 317),
  Beta = c(31, 32, 40, 0, 114,
           7, 347, 17, 91, 180,
           13, 209, 27, 61, 374),
  Interface_Residues = c(15, 14, 10, 0, 41,
                         3, 45, 13, 8, 40,
                         12, 36, 5, 14, 70)
)

# ALPHA / BETA CONTACTS IN MBD REGION

# Prepare data
mbd_data$Total <- mbd_data$Alpha + mbd_data$Beta

alpha_beta_data <- mbd_data %>%
  select(Paralog, Construct, Alpha, Beta) %>%
  pivot_longer(cols = c(Alpha, Beta), names_to = "Tubulin", values_to = "Contacts")

# Set stacking order (Alpha on bottom)
alpha_beta_data$Tubulin <- factor(alpha_beta_data$Tubulin,
                                  levels = c("Alpha", "Beta"))

alpha_beta_plot <- ggplot(alpha_beta_data, aes(x = Construct, y = Contacts, fill = Tubulin)) +
  
  # Bars
  geom_bar(stat = "identity", position = "stack") +
  
  # Total contacts label above each bar
  geom_text(data = mbd_data,
            aes(x = Construct, y = Total + 15, label = Total),
            inherit.aes = FALSE,
            size = 3,
            colour = "black",
            show.legend = FALSE) +
  
  facet_wrap(
    ~ Paralog,
    labeller = labeller(
      Paralog = c(
        "C1" = "CAMSAP1",
        "C2" = "CAMSAP2",
        "C3" = "CAMSAP3"
      )
    )
  ) +
  
  # y axis scale
  scale_y_continuous(
    limits = c(0, 800),
    breaks = seq(0, 800, 100)
  ) +
  
  labs(title = "Alpha vs Beta Tubulin Contacts (CC2_CC3_Linker Region)",
       x = "Construct", y = "Contacts", fill = "Tubulin") +
  
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  
  scale_fill_manual(values = c("Alpha" = "lightblue", "Beta" = "pink"))

print(alpha_beta_plot)
ggsave("Alpha_Beta_contacts.png", alpha_beta_plot, dpi = 300, width = 10, height = 6)


# ALPHA BETA CONSTRUCTS IN D2 REGION

# Alpha vs Beta for D2

# Prepare data
d2_data$Total <- d2_data$Alpha + d2_data$Beta
d2_data$Construct <- factor(d2_data$Construct,
                            levels = c("Full", "delCKK", "delCH", "delCC", "Isolated"))

alpha_beta_data <- d2_data %>%
  select(Paralog, Construct, Alpha, Beta) %>%
  pivot_longer(cols = c(Alpha, Beta), names_to = "Tubulin", values_to = "Contacts")

# Set stacking order (Alpha on bottom)
alpha_beta_data$Tubulin <- factor(alpha_beta_data$Tubulin,
                                  levels = c("Alpha", "Beta"))

alpha_beta_plot <- ggplot(alpha_beta_data, aes(x = Construct, y = Contacts, fill = Tubulin)) +
  
  # Bars
  geom_bar(stat = "identity", position = "stack") +
  
  # Total contacts label above each bar
  geom_text(data = d2_data,
            aes(x = Construct, y = Total + 15, label = Total),
            inherit.aes = FALSE,
            size = 3,
            colour = "black",
            show.legend = FALSE) +
  
  facet_wrap(
    ~ Paralog,
    labeller = labeller(
      Paralog = c(
        "C1" = "CAMSAP1",
        "C2" = "CAMSAP2",
        "C3" = "CAMSAP3"
      )
    )
  ) +
  
  # y axis scale
  scale_y_continuous(
    limits = c(0, 1000),
    breaks = seq(0, 1000, 100)
  ) +
  
  labs(title = "Alpha vs Beta Tubulin Contacts (CC3_CKK_Linker Region)",
       x = "Construct", y = "Contacts", fill = "Tubulin") +
  
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top",
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  
  scale_fill_manual(values = c("Alpha" = "lightblue", "Beta" = "pink"))

print(alpha_beta_plot)
ggsave("Alpha_Beta_D2_contacts.png", alpha_beta_plot, dpi = 300, width = 10, height = 6)