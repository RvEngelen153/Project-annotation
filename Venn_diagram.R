# === Load required packages ===
library(readxl)   # To read Excel files
library(dplyr)    # For data manipulation
library(grid)     # For low-level drawing of graphics (used here for manual Venn diagrams)

# === Load the data from Excel ===
file_path <- "plot table.xlsx"  # Path to the Excel file
data <- read_excel(file_path, sheet = "Blad1", skip = 1)
# Reads sheet "Blad1" from the Excel file, skipping the first row (e.g., header notes)

# === Rename columns to standardized names ===
colnames(data) <- c("VirusGroup", "Sample", "DIAMOND", "VIRify", "VIGA")

# === Replace NA values (missing data) with zero ===
data$DIAMOND[is.na(data$DIAMOND)] <- 0
data$VIRify[is.na(data$VIRify)] <- 0
data$VIGA[is.na(data$VIGA)] <- 0

# === Summarize hits per tool per sample and virus group ===
grouped_hits <- data %>%
  group_by(Sample, VirusGroup) %>%
  summarise(
    DIAMOND_hits = sum(DIAMOND),  # Sum of DIAMOND hits in each group
    VIRify_hits = sum(VIRify),    # Sum of VIRify hits
    VIGA_hits = sum(VIGA),        # Sum of VIGA hits
    .groups = "drop"              # Do not keep grouping structure
  ) %>%
  mutate(
    DIAMOND = DIAMOND_hits > 0,   # Logical: was there a DIAMOND hit?
    VIRify = VIRify_hits > 0,     # Logical: was there a VIRify hit?
    VIGA = VIGA_hits > 0          # Logical: was there a VIGA hit?
  )

# === Create output folder for plots ===
dir.create("venn_plots", showWarnings = FALSE)
# Creates the folder only if it doesn't exist; suppress warnings if it already exists

# === Get the list of unique samples ===
samples <- unique(grouped_hits$Sample)

# === Loop over each sample and create a custom Venn diagram ===
for (s in samples) {
  df <- subset(grouped_hits, Sample == s)  # Filter data for this sample
  
  d <- df$DIAMOND  # Logical vector: TRUE if DIAMOND hit exists
  v <- df$VIRify   # Logical vector: TRUE if VIRify hit exists
  g <- df$VIGA     # Logical vector: TRUE if VIGA hit exists

  # === Count hits for each region of the Venn diagram ===
  only_d_hits <- sum(df$DIAMOND_hits[d & !v & !g])
  only_v_hits <- sum(df$VIRify_hits[!d & v & !g])
  only_g_hits <- sum(df$VIGA_hits[!d & !v & g])
  
  d_v_hits     <- sum(df$DIAMOND_hits[d & v & !g] + df$VIRify_hits[d & v & !g])
  d_g_hits     <- sum(df$DIAMOND_hits[d & !v & g] + df$VIGA_hits[d & !v & g])
  v_g_hits     <- sum(df$VIRify_hits[!d & v & g] + df$VIGA_hits[!d & v & g])
  d_v_g_hits   <- sum(df$DIAMOND_hits[d & v & g] + df$VIRify_hits[d & v & g] + df$VIGA_hits[d & v & g])

  # === Start a new PNG file for plotting ===
  png(paste0("venn_plots/venn_", make.names(s), ".png"), width = 2000, height = 2000, res = 300)
  grid.newpage()  # Start fresh page for the plot

  # === Draw the three overlapping circles ===
  grid.circle(x = 0.4, y = 0.6, r = 0.25, gp = gpar(col = NA, fill = "skyblue", alpha = 0.5))    # DIAMOND
  grid.circle(x = 0.6, y = 0.6, r = 0.25, gp = gpar(col = NA, fill = "orange", alpha = 0.5))     # VIRify
  grid.circle(x = 0.5, y = 0.4, r = 0.25, gp = gpar(col = NA, fill = "lightgreen", alpha = 0.5)) # VIGA

  # === Add titles for each tool ===
  grid.text("DIAMOND", x = 0.25, y = 0.85, gp = gpar(fontsize = 12, fontface = "bold"))
  grid.text("VIRify", x = 0.75, y = 0.85, gp = gpar(fontsize = 12, fontface = "bold"))
  grid.text("VIGA", x = 0.5, y = 0.13, gp = gpar(fontsize = 12, fontface = "bold"))

  # === Add main title with sample name ===
  grid.text(paste("Venn Diagram -", s), y = 0.95, gp = gpar(fontsize = 16, fontface = "bold"))

  # === Add counts in the respective circle areas ===
  grid.text(only_d_hits, x = 0.25, y = 0.6, gp = gpar(fontsize = 20, fontface = "bold"))
  grid.text(only_v_hits, x = 0.75, y = 0.6, gp = gpar(fontsize = 20, fontface = "bold"))
  grid.text(only_g_hits, x = 0.5, y = 0.22, gp = gpar(fontsize = 20, fontface = "bold"))

  grid.text(d_v_hits, x = 0.5, y = 0.7, gp = gpar(fontsize = 20, fontface = "bold"))
  grid.text(d_g_hits, x = 0.37, y = 0.42, gp = gpar(fontsize = 20, fontface = "bold"))
  grid.text(v_g_hits, x = 0.63, y = 0.42, gp = gpar(fontsize = 20, fontface = "bold"))

  grid.text(d_v_g_hits, x = 0.5, y = 0.5, gp = gpar(fontsize = 20, fontface = "bold"))

  dev.off()  # Close and save the PNG file
}
