load("../simus/resData/PBD/rejection_sets_length_SC.RData")

library(ggplot2)
library(tidyr)
library(dplyr)

df_long <- final_results_df %>%
  pivot_longer(
    cols = starts_with("length_set_"),
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(
    Method = factor(
      gsub("length_set_", "", Method),
      levels = c("SC_cHBY", "SC_cDBY", "DBY", "cBY")
    )
  )

fig1 <- ggplot(df_long, aes(x = Method, y = Rejections, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "darkred", outlier.size = 1.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Comparison of rejection counts across procedures",
    subtitle = paste("Based on", nrow(final_results_df), "simulation replicates"),
    x = "Procedure",
    y = "Number of rejected hypotheses"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )


ggsave("plots/figure1.pdf", plot = fig1, units = "mm",width = 297,height = 210)



