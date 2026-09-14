# Load libraries
library(Biostrings)
library(ggplot2)

# Set working directory
setwd("C:/Users/paric/Documents/msc_project_new/blast")

# Read alignment
alignment <- readAAStringSet("maftt_alignment.fasta", format = "fasta")
seq_matrix <- as.matrix(alignment)
alignment_length <- ncol(seq_matrix)

# Calculate percent identity per position
window <- 30
half <- floor(window / 2)
identity_scores <- numeric(alignment_length)

for (i in 1:alignment_length) {
  residues <- seq_matrix[, i]
  residues <- residues[residues != "-"]
  if (length(residues) < 2) {
    identity_scores[i] <- 0
    next
  }
  total_pairs <- 0
  identical_pairs <- 0
  for (j in 1:(length(residues) - 1)) {
    for (k in (j + 1):length(residues)) {
      if (residues[j] == residues[k]) {
        identical_pairs <- identical_pairs + 1
      }
      total_pairs <- total_pairs + 1
    }
  }
  identity_scores[i] <- if (total_pairs > 0) (identical_pairs / total_pairs) * 100 else 0
}

# Apply sliding window smoothing
smoothed <- numeric(alignment_length)
for (i in 1:alignment_length) {
  start <- max(1, i - half)
  end <- min(alignment_length, i + half)
  smoothed[i] <- mean(identity_scores[start:end])
}

data <- data.frame(Position = 0:(alignment_length - 1), Identity = smoothed)

# Map domains
idx <- grep("Q5T5Y3", names(alignment))[1]
camsap1_seq <- as.character(alignment[[idx]])
aln_positions <- which(strsplit(camsap1_seq, "")[[1]] != "-")

domains <- data.frame(
  Domain = c("CH", "CC1", "CC2", "CC3", "CKK"),
  Uniprot_Start = c(216, 873, 1016, 1291, 1463),
  Uniprot_End = c(331, 909, 1048, 1343, 1597)
)

for (i in 1:nrow(domains)) {
  domains$Aln_Start[i] <- aln_positions[domains$Uniprot_Start[i]]
  domains$Aln_End[i] <- aln_positions[domains$Uniprot_End[i]]
}
domains$Aln_Start <- domains$Aln_Start - 1
domains$Aln_End <- domains$Aln_End - 1

# Plot
ggplot(data, aes(x = Position, y = Identity)) +
  geom_rect(data = domains,
            aes(xmin = Aln_Start, xmax = Aln_End, ymin = -Inf, ymax = Inf, fill = Domain),
            alpha = 0.15, inherit.aes = FALSE) +
  geom_text(data = domains,
            aes(x = (Aln_Start + Aln_End) / 2, y = 95, label = Domain),
            size = 3.5, inherit.aes = FALSE) +
  geom_line(colour = "black", linewidth = 1) +
  labs(title = "Percentage Identity of CAMSAP Homologs",
       x = "Alignment Position",
       y = "Identity (%)") +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, max(data$Position)),
    breaks = seq(0, max(data$Position), by = 200)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Save

ggsave("identity_plot.png", dpi = 300, width = 10, height = 5, units = "in")