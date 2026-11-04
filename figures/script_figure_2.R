load("../simus/resData/PBD/rejection_sets_length_SC2.RData")

library(ggplot2)
library(tidyr)
library(dplyr)

final_results_df[10] <- final_results_df[2] -final_results_df[4]
final_results_df[11] <- final_results_df[3] -final_results_df[5]
final_results_df[12] <- final_results_df[2] -final_results_df[3]
names(final_results_df)[10] <-"diff_(SC_cHBY-cBY)"
names(final_results_df)[11] <-"diff_(SC_cDBY-DBY)"
names(final_results_df)[12] <-"diff_(SC_cHBY-SC_cDBY)"

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
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )
fig1

ggsave("plots/figure3_eclosure.pdf", plot = fig1, units = "mm",width = 297,height = 210)

df_long <- final_results_df %>%
  pivot_longer(
    cols = starts_with("diff_"),
    names_to = "Method",
    values_to = "Rejections"
  ) %>%
  mutate(
    Method = factor(
      gsub("diff_", "", Method),
      levels = c("(SC_cHBY-cBY)", "(SC_cDBY-DBY)", "(SC_cHBY-SC_cDBY)")
    )
  )

fig2 <- ggplot(df_long, aes(x = Method, y = Rejections, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "darkred", outlier.size = 1.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Difference of rejection counts across procedures",
    subtitle = paste("Based on", nrow(final_results_df), "simulation replicates"),
    x = "Difference",
    y = "Number of rejected hypotheses"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave("plots/figure4_eclosure.pdf", plot = fig2, units = "mm",width = 297,height = 210)
fig2

