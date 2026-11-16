# Матриця пріоритетів: дані з «Результати тестування_Дідух.docx».

packages <- c("ggplot2", "ggrepel", "viridisLite")
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(ggplot2)
library(ggrepel)

df <- data.frame(
  Biotope = c("R15.1", "R15.2", "R15.3", "R15.4", "R15.5", "R16b.2", "R13.2"),
  Z  = c(0.27, 0.31, 0.40, 0.39, 0.39, 0.33, 0.57),
  T  = c(0.93, 0.73, 1.00, 1.00, 0.93, 0.68, 0.47),
  V  = c(0.60, 0.40, 0.60, 0.50, 0.50, 0.40, 0.23),
  ST = c(0.45, 0.53, 0.53, 0.53, 0.53, 0.53, 0.48)
)

p <- ggplot(df, aes(x = Z, y = T)) +
  geom_vline(xintercept = c(0.3, 0.6),
             linetype = "dashed",
             colour = "grey70") +
  geom_hline(yintercept = c(0.3, 0.6),
             linetype = "dashed",
             colour = "grey70") +
  geom_point(aes(colour = V, size = ST),
             alpha = 0.9) +
  geom_text_repel(aes(label = Biotope),
                  size = 5,
                  fontface = "bold",
                  max.overlaps = Inf) +
  scale_colour_viridis_c(
    option = "plasma",
    direction = 1,
    name = "V"
  ) +
  scale_size(
    range = c(5, 12),
    name = "ST"
  ) +
  scale_x_reverse(
    limits = c(1, 0),
    breaks = c(1, 0.6, 0.3, 0.1),
    name = "Z"
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = c(0.1, 0.3, 0.6, 1),
    name = "T"
  ) +
  coord_fixed() +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.title = element_text(face = "bold")
  )

print(p)

# TIFF буде збережено в поточній робочій папці R (getwd()).
ggsave(
  "Priority_matrix.tiff",
  p,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
