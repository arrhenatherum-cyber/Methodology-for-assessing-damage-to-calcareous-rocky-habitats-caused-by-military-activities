# Аналіз чутливості методики: оновлені дані станом на 30.10.2026.
# Джерело: «Для валідації методики.xlsx», аркуш «Аркуш1».
# Первинні 161 оцінка перенесені нижче; вибирати Excel-файл не потрібно.
# Відкрийте цей файл у RStudio та натисніть Source.
# Результати: папка Habitat_sensitivity_20261030 у поточній папці getwd().
# Це аналіз внутрішньої чутливості, а не зовнішня польова валідація.

packages <- c("ggplot2", "patchwork")
missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
library(ggplot2)
library(patchwork)

seed <- 20261030L
n_sim <- 10000L
p_base_values <- c(0.60, 0.80, 0.90)
primary_p <- 0.60
out_dir <- file.path(getwd(), "Habitat_sensitivity_20261030")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. ПЕРВИННІ ОЦІНКИ ----------------------------------------------------
# Z: вища оцінка = краща збереженість.
# T: вища оцінка = вища созологічна цінність.
# V: вища оцінка = складніше / триваліше відновлення.
# ST: вища оцінка = сприятливіша просторово-часова характеристика.
# Чотири індекси аналізуються окремо; сукупний індекс не обчислюється.

df <- read.table(header = TRUE, stringsAsFactors = FALSE, text = '
Biotope z1 z2 z3 z4 z5 z6 z7 z8 z9 z10 t1 t2 t3 t4 t5 t6 v1 v2 v3 st1 st2 st3 st4
R15.1 0.1 0.1 0.3 0.1 0.3 0.6 0.3 0.3 0.3 0.3 1.0 1.0 0.6 1.0 1.0 1.0 0.6 0.6 0.6 0.3 0.6 0.6 0.3
R15.2 0.1 0.1 0.1 0.1 0.3 0.6 0.3 0.6 0.3 0.6 1.0 0.6 0.6 0.6 1.0 0.6 0.3 0.6 0.3 0.3 0.6 0.6 0.6
R15.3 0.1 0.1 0.1 0.3 0.3 1.0 0.6 0.3 0.6 0.6 1.0 1.0 1.0 1.0 1.0 1.0 0.6 0.6 0.6 0.3 0.6 0.6 0.6
R15.4 0.1 0.1 0.1 0.6 0.6 0.6 0.3 0.3 0.6 0.6 1.0 1.0 1.0 1.0 1.0 1.0 0.6 0.6 0.3 0.3 0.6 0.6 0.6
R15.5 0.1 0.1 0.1 0.6 0.3 0.6 0.6 0.3 0.6 0.6 1.0 1.0 0.6 1.0 1.0 1.0 0.6 0.6 0.3 0.3 0.6 0.6 0.6
R16b.2 0.1 0.1 0.1 0.3 0.6 0.6 0.3 0.3 0.3 0.6 1.0 0.6 0.3 0.6 1.0 0.6 0.3 0.6 0.3 0.3 0.6 0.6 0.6
R13.2 0.6 0.6 0.3 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.3 0.6 0.6 0.1 0.6 0.3 0.1 0.3 0.1 0.6 0.6 0.6
')

groups <- list(Z = paste0("z", 1:10), T = paste0("t", 1:6),
               V = paste0("v", 1:3), ST = paste0("st", 1:4))
index_names <- names(groups)
criterion_names <- unlist(groups, use.names = FALSE)
habitats <- df$Biotope
scale_integer <- c(1L, 3L, 6L, 10L)

stopifnot(nrow(df) == 7L, !anyDuplicated(habitats),
          identical(names(df)[-1], criterion_names))
raw <- as.matrix(df[, criterion_names])
stopifnot(!anyNA(raw), all(is.finite(raw)))
# Суми цілих балів запобігають хибному розділенню однакових рангів.
scores <- round(raw * 10)
stopifnot(all(abs(raw * 10 - scores) < 1e-10),
          all(scores %in% scale_integer))
storage.mode(scores) <- "integer"
colnames(scores) <- criterion_names
rownames(scores) <- habitats

baseline <- data.frame(Biotope = habitats)
for (ind in index_names) {
  baseline[[ind]] <- rowSums(scores[, groups[[ind]], drop = FALSE]) /
    (10 * length(groups[[ind]]))
}

# Контроль проти актуальної таблиці, з точними середніми до округлення.
expected <- cbind(
  Z = c(.27, .31, .40, .39, .39, .33, .57),
  T = c(5.6, 4.4, 6, 6, 5.6, 4.1, 2.8) / 6,
  V = c(1.8, 1.2, 1.8, 1.5, 1.5, 1.2, .7) / 3,
  ST = c(.45, .525, .525, .525, .525, .525, .475)
)
stopifnot(max(abs(as.matrix(baseline[, index_names]) - expected)) < 1e-12)
# Якщо первинні оцінки згодом зміняться, оновіть також expected.

write_csv <- function(x, name) {
  write.csv(x, file.path(out_dir, name), row.names = FALSE,
            fileEncoding = "UTF-8", na = "NA")
}
safe_cor <- function(x, y, method) {
  if (length(unique(x)) < 2L || length(unique(y)) < 2L) return(NA_real_)
  unname(cor(x, y, method = method))
}
q95 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  unname(quantile(x, c(.025, .975), type = 7))
}

# 2. LOCO: ВИЛУЧЕННЯ ОДНОГО КРИТЕРІЮ -----------------------------------
# Решта критеріїв усереднюються з однаковими вагами.
# Кореляції обчислюються між сімома біотопами всередині кожного індексу.
# Однаковим значенням надаються середні ранги; Kendall враховує зв'язки (tau-b).

loco_rows <- list()
loco_summary_rows <- list()
rank_direction <- c(Z = 1, T = -1, V = -1, ST = 1)
for (ind in index_names) {
  b <- baseline[[ind]]
  base_rank <- rank(rank_direction[[ind]] * b, ties.method = "average")
  for (removed in groups[[ind]]) {
    kept <- setdiff(groups[[ind]], removed)
    changed <- rowSums(scores[, kept, drop = FALSE]) / (10 * length(kept))
    abs_change <- abs(changed - b)
    rel_change <- 100 * abs_change / b
    new_rank <- rank(rank_direction[[ind]] * changed, ties.method = "average")
    loco_rows[[length(loco_rows) + 1L]] <- data.frame(
      Biotope = habitats, Index = ind, Removed_criterion = removed,
      Baseline = b, Modified = changed, Signed_change = changed - b,
      Absolute_change = abs_change, Relative_change_percent = rel_change,
      Baseline_rank = base_rank, Modified_rank = new_rank,
      Absolute_rank_change = abs(new_rank - base_rank)
    )
    loco_summary_rows[[length(loco_summary_rows) + 1L]] <- data.frame(
      Index = ind, Removed_criterion = removed,
      Mean_absolute_change = mean(abs_change),
      Maximum_absolute_change = max(abs_change),
      Mean_relative_change_percent = mean(rel_change),
      Maximum_relative_change_percent = max(rel_change),
      Spearman_rho = safe_cor(b, changed, "spearman"),
      Kendall_tau_b = safe_cor(b, changed, "kendall"),
      Number_changed_ranks = sum(new_rank != base_rank),
      Maximum_rank_change = max(abs(new_rank - base_rank))
    )
  }
}
loco <- do.call(rbind, loco_rows)
loco_summary <- do.call(rbind, loco_summary_rows)

# 3. МОНТЕ-КАРЛО ---------------------------------------------------------
# Кожна оцінка незалежно залишається незмінною з імовірністю p_base.
# Інакше вона змінюється на ОДИН сусідній рівень шкали {0.1, 0.3, 0.6, 1.0}.
# Для 0.3 і 0.6 решта ймовірності порівну розподілена між двома сусідами.
# На краях 0.1 і 1.0 вся ймовірність зміни переходить до єдиного сусіда.
# p_base = 0.60 — сценарій з попереднього коду, а не виміряна точність експертів.
# 0.80 і 0.90 — додаткові, м'якші сценарії. Не обирайте сценарій за результатом.
# Через нерівні кроки та межі шкали ця модель НЕ є незміщеним шумом.
# Інтервали — умовні квантілі симуляцій, не емпіричні довірчі інтервали.

rank_correlations <- function(b, simulations) {
  rb <- rank(b, ties.method = "average")
  rs <- t(apply(simulations, 1, rank, ties.method = "average"))
  cb <- rb - mean(rb)
  cs <- sweep(rs, 1, rowMeans(rs), "-")
  den <- sqrt(sum(cb^2) * rowSums(cs^2))
  rho <- as.vector(cs %*% cb) / den
  pairs <- combn(length(b), 2)
  sb <- sign(b[pairs[1, ]] - b[pairs[2, ]])
  ss <- vapply(seq_len(ncol(pairs)), function(k) {
    sign(simulations[, pairs[1, k]] - simulations[, pairs[2, k]])
  }, numeric(nrow(simulations)))
  den_tau <- sqrt(sum(sb != 0) * rowSums(ss != 0))
  tau <- as.vector(ss %*% sb) / den_tau
  rho[den == 0] <- NA_real_
  tau[den_tau == 0] <- NA_real_
  # Перевірка векторизованих обчислень стандартною функцією R.
  for (k in c(1L, 18L, 100L)) {
    stopifnot(isTRUE(all.equal(rho[k], safe_cor(b, simulations[k, ], "spearman"),
                               tolerance = 1e-12)),
              isTRUE(all.equal(tau[k], safe_cor(b, simulations[k, ], "kendall"),
                               tolerance = 1e-12)))
  }
  data.frame(Spearman = rho, Kendall = tau)
}

set.seed(seed)
# Спільні випадкові числа для порівняння сценаріїв; клітинки незалежні.
uniforms <- array(runif(n_sim * nrow(scores) * ncol(scores)),
                  dim = c(n_sim, nrow(scores), ncol(scores)))
mc_rows <- list()
mc_rank_rows <- list()
primary_samples <- list()

for (p_base in p_base_values) {
  for (ind in index_names) {
    columns <- match(groups[[ind]], criterion_names)
    simulated_sums <- matrix(0L, nrow = n_sim, ncol = nrow(scores))
    expected_sums <- numeric(nrow(scores))
    for (i in seq_len(nrow(scores))) {
      for (j in columns) {
        pos <- match(scores[i, j], scale_integer)
        u <- uniforms[, i, j]
        if (pos == 1L) {
          values <- scale_integer[1:2]
          probs <- c(p_base, 1 - p_base)
          draw <- ifelse(u < p_base, values[1], values[2])
        } else if (pos == 4L) {
          values <- scale_integer[3:4]
          probs <- c(1 - p_base, p_base)
          draw <- ifelse(u < 1 - p_base, values[1], values[2])
        } else {
          values <- scale_integer[(pos - 1L):(pos + 1L)]
          probs <- c((1 - p_base) / 2, p_base, (1 - p_base) / 2)
          draw <- ifelse(u < probs[1], values[1],
                         ifelse(u < sum(probs[1:2]), values[2], values[3]))
        }
        simulated_sums[, i] <- simulated_sums[, i] + draw
        expected_sums[i] <- expected_sums[i] + sum(values * probs)
      }
    }
    samples <- simulated_sums / (10 * length(columns))
    colnames(samples) <- habitats
    theoretical_mean <- expected_sums / (10 * length(columns))
    if (p_base == primary_p) primary_samples[[ind]] <- samples
    for (i in seq_along(habitats)) {
      limits <- q95(samples[, i])
      mc_rows[[length(mc_rows) + 1L]] <- data.frame(
        P_base = p_base, Biotope = habitats[i], Index = ind,
        Baseline = baseline[[ind]][i], Mean = mean(samples[, i]),
        Median = median(samples[, i]), SD = sd(samples[, i]),
        Lower_95 = limits[1], Upper_95 = limits[2],
        Theoretical_mean = theoretical_mean[i],
        Expected_shift = theoretical_mean[i] - baseline[[ind]][i],
        Observed_shift = mean(samples[, i]) - baseline[[ind]][i]
      )
    }
    correlations <- rank_correlations(baseline[[ind]], samples)
    for (metric in names(correlations)) {
      vals <- correlations[[metric]]
      limits <- q95(vals)
      mc_rank_rows[[length(mc_rank_rows) + 1L]] <- data.frame(
        P_base = p_base, Index = ind, Metric = metric, N_sim = n_sim,
        N_undefined = sum(!is.finite(vals)),
        Median = median(vals[is.finite(vals)]),
        Lower_95 = limits[1], Upper_95 = limits[2]
      )
    }
  }
}
mc <- do.call(rbind, mc_rows)
mc_rank <- do.call(rbind, mc_rank_rows)
mc_primary <- mc[mc$P_base == primary_p, ]

# 4. ЗБЕРЕЖЕННЯ ТАБЛИЦЬ ------------------------------------------------
write_csv(df, "00_primary_criteria.csv")
write_csv(baseline, "01_baseline_indices_exact.csv")
write_csv(loco, "02_LOCO_full.csv")
write_csv(loco_summary, "03_LOCO_summary.csv")
write_csv(mc, "04_MC_intervals_all_scenarios.csv")
write_csv(mc_primary, "05_MC_intervals_primary.csv")
write_csv(mc_rank, "06_MC_rank_stability.csv")
saveRDS(primary_samples, file.path(out_dir, "MC_primary_samples.rds"))
writeLines(c(
  paste("seed:", seed), paste("n_sim:", n_sim),
  paste("p_base:", paste(p_base_values, collapse = ", ")),
  "Source: Для валідації методики.xlsx; Аркуш1; 30.10.2026",
  "Independent adjacent-level perturbations; conditional sensitivity scenarios.",
  "Boundary rule: all change probability goes to the sole neighbour.",
  "Full precision; equal scores have average ranks; Kendall tau-b.",
  "95% simulation ranges are not empirical confidence intervals.",
  "Identical input profiles can have slightly different MC summaries due to sampling.",
  "No combined priority index or tests of statistical significance."
), file.path(out_dir, "analysis_settings.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))

# 5. РИСУНКИ ------------------------------------------------------------
index_colours <- c(Z = "#0072B2", T = "#009E73", V = "#D55E00", ST = "#7B3294")
plot_theme <- theme_bw(base_size = 11) + theme(
  panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
  strip.text = element_text(face = "bold"),
  plot.title = element_text(face = "bold", size = 12),
  legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 9)
)
loco_summary$Index <- factor(loco_summary$Index, levels = index_names)
loco_summary$Removed_criterion <- factor(loco_summary$Removed_criterion,
                                         levels = criterion_names)

p_a <- ggplot(loco_summary, aes(Removed_criterion,
                                Mean_relative_change_percent, fill = Index)) +
  geom_col(width = .72, show.legend = FALSE) +
  geom_point(aes(y = Maximum_relative_change_percent), shape = 21,
             fill = "white", colour = "black", size = 2, show.legend = FALSE) +
  facet_wrap(~Index, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = index_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "A. Leave-one-criterion-out sensitivity",
       x = "Removed criterion", y = "Absolute relative change (%)",
       caption = "Bars: mean across seven habitats. Circles: maximum across habitats.") +
  plot_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

rank_plot <- rbind(
  data.frame(loco_summary[, c("Index", "Removed_criterion")],
             Metric = "Spearman rho", Value = loco_summary$Spearman_rho),
  data.frame(loco_summary[, c("Index", "Removed_criterion")],
             Metric = "Kendall tau-b", Value = loco_summary$Kendall_tau_b)
)
p_b <- ggplot(rank_plot, aes(Removed_criterion, Value,
                             colour = Metric, shape = Metric, group = Metric)) +
  geom_hline(yintercept = 0, colour = "grey75", linetype = "dashed") +
  geom_line(linewidth = .5) + geom_point(size = 2) +
  facet_wrap(~Index, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("Spearman rho" = "#0072B2", "Kendall tau-b" = "#D55E00")) +
  scale_shape_manual(values = c("Spearman rho" = 16, "Kendall tau-b" = 17)) +
  scale_y_continuous(limits = c(-1, 1.02), breaks = seq(-1, 1, .5)) +
  labs(title = "B. Rank stability after removing one criterion",
       x = "Removed criterion", y = "Rank correlation", colour = NULL, shape = NULL) +
  plot_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

mc_primary$Index <- factor(mc_primary$Index, levels = index_names)
mc_primary$Biotope <- factor(mc_primary$Biotope, levels = rev(habitats))
p_c <- ggplot(mc_primary, aes(y = Biotope, colour = Index)) +
  geom_segment(aes(x = Lower_95, xend = Upper_95, yend = Biotope),
               linewidth = .8, show.legend = FALSE) +
  geom_point(aes(x = Mean, shape = "Simulation mean"), size = 2.4) +
  geom_point(aes(x = Baseline, shape = "Baseline"), size = 2.6, colour = "black") +
  facet_wrap(~Index, nrow = 1) +
  scale_colour_manual(values = index_colours, guide = "none") +
  scale_shape_manual(values = c("Simulation mean" = 16, "Baseline" = 4)) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, .5, 1)) +
  labs(title = paste0("C. Monte Carlo sensitivity (N = ", n_sim,
                      "; p_base = ", primary_p, ")"),
       x = "Index value", y = NULL, shape = NULL,
       caption = "Lines: conditional 2.5th-97.5th percentiles; they are not empirical confidence intervals.") +
  plot_theme

combined <- (p_a / p_b / p_c) + plot_layout(heights = c(1, 1.1, 1.1))
print(combined)
ggsave(file.path(out_dir, "Sensitivity_combined.tiff"), combined,
       width = 10, height = 11, units = "in", dpi = 600, compression = "lzw",
       bg = "white")
ggsave(file.path(out_dir, "Sensitivity_combined.png"), combined,
       width = 10, height = 11, units = "in", dpi = 200, bg = "white")
for (name in c("A_LOCO", "B_Rank_stability", "C_Monte_Carlo")) {
  plot <- switch(name, A_LOCO = p_a, B_Rank_stability = p_b, C_Monte_Carlo = p_c)
  ggsave(file.path(out_dir, paste0(name, ".tiff")), plot,
         width = 10, height = 3.8, units = "in", dpi = 600,
         compression = "lzw", bg = "white")
}

cat("\nБазові індекси (розрахунки без округлення):\n")
print(baseline, digits = 6, row.names = FALSE)
cat("\nLOCO: мінімальна рангова кореляція для кожного індексу:\n")
print(aggregate(cbind(Spearman_rho, Kendall_tau_b) ~ Index,
                data = loco_summary, FUN = min), row.names = FALSE)
cat("\nМонте-Карло: стабільність рангів у трьох сценаріях:\n")
print(mc_rank, digits = 4, row.names = FALSE)
cat("\nЗбережено у:", normalizePath(out_dir), "\n")
