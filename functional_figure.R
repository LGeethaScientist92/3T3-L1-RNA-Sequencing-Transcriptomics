library(tidyverse)
library(patchwork)

# -----------------------------
# 1. Load data
# -----------------------------
df <- read_csv("functional_assays.csv")

# Order factors
df <- df %>%
  mutate(
    condition = factor(condition, levels = c("Control", "R130", "Q130")),
    dose = factor(dose, levels = c(0, 25, 50, 100)),
    stimulus = factor(stimulus, levels = c("None", "Basal", "Isoproterenol", "Insulin"))
  )

# Publication colors
group_cols <- c(
  "Control" = "#9E9E9E",
  "R130"    = "#2C7FB8",
  "Q130"    = "#D95F02"
)

# Theme
theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Summary function
sum_df <- function(data) {
  data %>%
    group_by(condition, dose, stimulus) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}

# -----------------------------
# 2. Figure panel B: Adipogenesis
# -----------------------------
adipo <- df %>% filter(assay == "Adipogenesis")
adipo_sum <- sum_df(adipo)

p_adipo <- ggplot(adipo_sum, aes(x = dose, y = mean, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_point(
    data = adipo,
    aes(x = dose, y = value, color = condition),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
    size = 2,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "B  Lipid accumulation",
    x = "IL-13 concentration (ng/mL)",
    y = "Lipid accumulation (% control)"
  ) +
  theme_pub

# -----------------------------
# 3. Figure panels C and D: Lipolysis
# -----------------------------
lipo <- df %>% filter(assay == "Lipolysis")
lipo_sum <- sum_df(lipo)

p_lipo_basal <- lipo_sum %>%
  filter(stimulus == "Basal") %>%
  ggplot(aes(x = dose, y = mean, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_point(
    data = lipo %>% filter(stimulus == "Basal"),
    aes(x = dose, y = value, color = condition),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
    size = 2,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "C  Basal lipolysis",
    x = "IL-13 concentration (ng/mL)",
    y = "Glycerol release"
  ) +
  theme_pub

p_lipo_iso <- lipo_sum %>%
  filter(stimulus == "Isoproterenol") %>%
  ggplot(aes(x = dose, y = mean, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_point(
    data = lipo %>% filter(stimulus == "Isoproterenol"),
    aes(x = dose, y = value, color = condition),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
    size = 2,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "D  Isoproterenol-stimulated lipolysis",
    x = "IL-13 concentration (ng/mL)",
    y = "Glycerol release"
  ) +
  theme_pub

# Net response = stimulated - basal
lipo_wide <- lipo %>%
  select(condition, dose, replicate, stimulus, value) %>%
  pivot_wider(names_from = stimulus, values_from = value) %>%
  mutate(delta = Isoproterenol - Basal)

lipo_delta_sum <- lipo_wide %>%
  group_by(condition, dose) %>%
  summarise(
    mean = mean(delta, na.rm = TRUE),
    sd = sd(delta, na.rm = TRUE),
    .groups = "drop"
  )

p_lipo_delta <- ggplot(lipo_delta_sum, aes(x = dose, y = mean, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_point(
    data = lipo_wide,
    aes(x = dose, y = delta, color = condition),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
    size = 2,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "E  Net lipolytic response",
    x = "IL-13 concentration (ng/mL)",
    y = expression(Delta*" glycerol release")
  ) +
  theme_pub

# -----------------------------
# 4. Figure panels F, G, H: Glucose uptake
# -----------------------------
gluc <- df %>% filter(assay == "Glucose")
gluc_sum <- sum_df(gluc)

p_gluc_basal <- gluc_sum %>%
  filter(stimulus == "Basal") %>%
  ggplot(aes(x = dose, y = mean, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_point(
    data = gluc %>% filter(stimulus == "Basal"),
    aes(x = dose, y = value, color = condition),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
    size = 2,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "F  Basal glucose uptake",
    x = "IL-13 concentration (ng/mL)",
    y = "2-DG uptake"
  ) +
  theme_pub

p_gluc_ins <- gluc_sum %>%
  filter(stimulus == "Insulin") %>%
  ggplot(aes(x = dose, y = mean, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_point(
    data = gluc %>% filter(stimulus == "Insulin"),
    aes(x = dose, y = value, color = condition),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
    size = 2,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "G  Insulin-stimulated glucose uptake",
    x = "IL-13 concentration (ng/mL)",
    y = "2-DG uptake"
  ) +
  theme_pub

# Fold induction = insulin / basal
gluc_wide <- gluc %>%
  select(condition, dose, replicate, stimulus, value) %>%
  pivot_wider(names_from = stimulus, values_from = value) %>%
  mutate(fold_change = Insulin / Basal)

gluc_fc_sum <- gluc_wide %>%
  group_by(condition, dose) %>%
  summarise(
    mean = mean(fold_change, na.rm = TRUE),
    sd = sd(fold_change, na.rm = TRUE),
    .groups = "drop"
  )

p_gluc_fc <- ggplot(gluc_fc_sum, aes(x = dose, y = mean, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  geom_point(
    data = gluc_wide,
    aes(x = dose, y = fold_change, color = condition),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
    size = 2,
    alpha = 0.9,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(
    title = "H  Fold induction by insulin",
    x = "IL-13 concentration (ng/mL)",
    y = "Fold change"
  ) +
  theme_pub

# -----------------------------
# 5. Optional placeholder for microscopy panel A
# -----------------------------
p_placeholder <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "A  Microscopy images\n(insert separately in Illustrator/PowerPoint)", size = 5) +
  theme_void()

# -----------------------------
# 6. Assemble multi-panel figure
# -----------------------------
final_plot <- (
  p_placeholder + p_adipo
) /
(
  p_lipo_basal + p_lipo_iso + p_lipo_delta
) /
(
  p_gluc_basal + p_gluc_ins + p_gluc_fc
)

# Show
print(final_plot)

# Save
ggsave(
  filename = "Figure1_functional_characterisation.png",
  plot = final_plot,
  width = 14,
  height = 16,
  dpi = 300
)