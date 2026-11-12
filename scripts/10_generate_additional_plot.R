# 10_generate_additional_plot.R
library(ggplot2)
library(tibble)
library(viridis)
library(forcats)
library(dplyr)

results <- read.csv("./output/02_generate_mps/performance/aggregated_test_results.csv")
main_out_dir <- file.path("./output/10_generate_additional_plot")

ranking_method_labels <- c(
  `mutual_information` = "Mutual Information",
  `lasso_stability_selection` = "Lasso Stability Selection"
)

dir.create(main_out_dir, showWarnings = FALSE)
results <- dplyr::distinct(results, outer_fold, ranking_method, training_method, n_probes, Rsquared, .keep_all = TRUE)

results <- results[results$Rsquared_ss >= 0,]
top_results <- results %>%
  arrange(desc(Rsquared_ss)) %>%
  slice_head(n = 20) %>%
  mutate(
    model_label = paste0(
      "F", outer_fold,
      " | ", training_method,
      " | ", n_probes, " probes"
    ),
    model_label = fct_reorder(model_label, Rsquared_ss)
  )

p_top <- ggplot(top_results, aes(x = Rsquared_ss, y = model_label, fill = ranking_method)) +

  geom_col(alpha = 1) +
  geom_text(aes(label = sprintf("%.2f", Rsquared_ss)), hjust = -0.2, size = 3.5, color = "black") +
  labs(
    x = expression("Performance" ~ (italic(R)^2)),
    y = "Model Config",
    fill = "Ranking Method"
  ) +
  scale_x_continuous(limits = c(0, max(top_results$Rsquared_ss) * 1.1)) +
  scale_fill_viridis_d(option = "B", begin = 0.15, end = 0.85, labels = ranking_method_labels) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(margin = margin(b = 5)),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size = 10)
  )

filename_top <- file.path(main_out_dir, "top_performance.pdf")
ggsave(filename_top, plot = p_top, width = 9, height = 6)