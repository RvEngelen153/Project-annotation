# Vereiste packages
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# 1. Data inlezen
df <- read_excel("plot_table.xlsx", skip = 1)
colnames(df) <- c("Virus", "Sample", "DIAMOND", "VIRify", "VIGA")

# 2. Schoonmaken en voorbereiden
df <- df %>%
  filter(!is.na(Virus)) %>%
  mutate(across(c(DIAMOND, VIRify, VIGA), as.numeric)) %>%
  mutate(
    total_hits = DIAMOND + VIRify + VIGA
  ) %>%
  arrange(desc(total_hits)) %>%
  mutate(virus_id = factor(Virus, levels = unique(Virus)))  # x-as volgorde

# 3. Matrixdata: lang formaat, en bepaal of er minstens 1 hit was per tool
matrix_data <- df %>%
  select(virus_id, DIAMOND, VIRify, VIGA) %>%
  pivot_longer(cols = c(DIAMOND, VIRify, VIGA), names_to = "tool", values_to = "hits") %>%
  mutate(hit = hits > 0)

# 4. Barplot: total hits per virus
bar <- ggplot(df, aes(x = virus_id, y = total_hits)) +
  geom_col(fill = "steelblue") +
  labs(y = "Totale hits", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

# 5. Matrixplot eronder
matrix <- ggplot(matrix_data, aes(x = virus_id, y = tool)) +
  geom_point(aes(color = hit), size = 4) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# 6. Combineer beide plots
(bar / matrix) +
  plot_annotation(title = "Totaal aantal hits per virus met tooldetectie (bolletjes)")
