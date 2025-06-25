library(ggplot2)
library(reshape2)

corr <- read.table("corr.txt", header = TRUE, row.names = 1, check.names = FALSE)
shrinkage <- read.table("shrinkage_weight.txt", header = TRUE, row.names = 1, check.names = FALSE)

corr_long <- melt(as.matrix(corr))
shrinkage_long <- melt(as.matrix(shrinkage))
colnames(corr_long) <- c("row", "col", "value")
colnames(shrinkage_long) <- c("row", "col", "shrinkage")

df <- merge(corr_long, shrinkage_long, by = c("row", "col"))

df$base_color <- NA
df$base_color[grepl("_sp$", df$row)] <- "blue"
df$base_color[grepl("_ex$", df$row)] <- "red"

df$fill_color <- mapply(function(color, alpha) {
  scales::alpha(color, ifelse(alpha >= 0.5, 1, 0.1))
}, df$base_color, df$shrinkage)

df$text_color <- ifelse(df$shrinkage >= 0.5, "white", "black")

df$row <- factor(df$row, levels = rev(rownames(corr)))
df$col <- factor(df$col, levels = colnames(corr))

p <- ggplot(df, aes(x = col, y = row)) +
  geom_tile(aes(fill = fill_color), color = "grey90", linewidth = 0.5, width = 0.95, height = 0.95) +
  geom_text(aes(label = round(value, 2), color = text_color), size = 2, fontface = "bold") +
  scale_fill_identity() +
  scale_color_identity() +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

plot(p)

ggsave("MBD.tiff", plot = p, width = 2000, height = 2000 * 12 / 22, dpi = 300, unit = "px")
