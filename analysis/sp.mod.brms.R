############################################################
## FULL WORKFLOW: from `net` to Fig S5-style panels (A–C)
## Replaces your OLS decomposition + partial residual plots
## with a BRMS measurement-error model:
##   degC | se(degC.se, sigma=TRUE) ~ 0 + me(estimate, SE) + ngen_centered
##
## What you get (mirrors your original figure structure):
##   A) Decomposition stacked bars: modeled contributions of
##      - Population size (density change)
##      - Symbiont genera
##      - Unexplained (residual-like)
##      + observed point with ±SE
##   B) “Partial residual” style plot for density (adjusted for ngen=1)
##   C) “Partial residual” style plot for ngen (adjusted for density change=0)
##
## Notes / Interpretation:
## - The PaBa line you used earlier is still fine for the bivariate plot.
## - This BRMS model replaces ONLY the subsequent multivariable OLS step.
## - Decomposition bars use posterior medians for each component (keeps your style).
##   Component medians won’t sum EXACTLY to degC, but should be very close.
## - Partial residual analog uses posterior-median slopes + a residual-like term
##   defined from the posterior-median conditional mean; this preserves the
##   species-level scatter behavior you liked in your original plots.
############################################################

## ----------------------------
## 0) Libraries
## ----------------------------
library(tidyverse)
library(forcats)
library(ggrepel)
library(cowplot)
library(brms)
library(tidybayes)

## ----------------------------
## 1) Starting point: build `net` exactly like you did
##    (Assumes you already have sp.degC, densdiffs, symsumm.ag in memory,
##     and abbrev_names, myticks2 defined elsewhere as in your script.)
## ----------------------------

# densdiffs <- readRDS("data/processed/dens20142023.rds")  # uncomment if needed

net <- sp.degC %>%
  left_join(select(as_tibble(densdiffs), Species, estimate, SE), by = "Species") %>%
  left_join(symsumm.ag, by = "Species") %>%
  drop_na(estimate)

## ----------------------------
## 2) Clean data for the model (same columns as your OLS workflow)
## ----------------------------
net_clean <- na.omit(net[, c("Species", "degC", "estimate", "ngen", "degC.se", "SE", "D")])

# Center ngen so that ngen_centered = 0 corresponds to ngen = 1
net_clean$ngen_centered <- net_clean$ngen - 1

# Strip any names/attributes and coerce to numeric (prevents annoying failures)
net_clean <- net_clean %>%
  mutate(
    degC.se = unname(as.numeric(degC.se)),
    SE      = unname(as.numeric(SE)),
    estimate = as.numeric(estimate),
    ngen_centered = as.numeric(ngen_centered)
  )

# Basic checks (brms will fail if SEs are missing or <= 0)
stopifnot(all(net_clean$degC.se > 0), all(net_clean$SE > 0))

## ----------------------------
## 3) Fit the BRMS measurement-error model
##    - se(degC.se, sigma=TRUE): your degC estimates have known SE, plus extra residual scatter
##    - me(estimate, SE): accounts for measurement error in your density-change predictor
##    - 0 + ... : keeps intercept fixed at 0 when estimate = 0 and ngen_centered = 0 (ngen = 1)
## ----------------------------

# Weakly-informative priors (adjust if you have stronger expectations)
pri <- c(
  prior(normal(0, 1), class = "b"),
  prior(exponential(2), class = "sigma")
)

fit_me_add <- brm(
  bf(degC | se(degC.se, sigma = TRUE) ~ 0 + me(estimate, SE) + ngen_centered),
  data = net_clean,
  family = stats::gaussian(),
  prior = pri,
  chains = 4, cores = 4, iter = 4000, warmup = 1500,
  seed = 8,
  control = list(adapt_delta = 0.99)
)

# Optional quick checks
print(fit_me_add)
# pp_check(fit_me_add)   # uncomment if you want a quick PPC
# plot(fit_me_add)       # trace plots

## ----------------------------
## 4) Extract posterior draws for slopes (robust to brms version naming)
## ----------------------------
draws_beta <- as_draws_df(fit_me_add)
b_est  <- draws_beta[["bsp_meestimateSE"]]
b_ngen <- draws_beta[["b_ngen_centered"]]
n_draw <- nrow(draws_beta)

## ----------------------------
## 5) BRMS analog of your OLS decomposition (posterior contributions)
##    Compute contributions for each posterior draw × species, then summarize
## ----------------------------

# vectors by species (row order = net_clean rows)
est_i <- net_clean$estimate
ng_i  <- net_clean$ngen_centered
y_i   <- net_clean$degC

# Contribution matrices: [draw x species]
Est_contrib_mat  <- outer(b_est,  est_i)  # population size component
Ngen_contrib_mat <- outer(b_ngen, ng_i)   # symbiont genera component

# Conditional mean per draw/species
mu_mat <- Est_contrib_mat + Ngen_contrib_mat

# Residual-like term per draw/species: observed - predicted mean
Resid_mat <- matrix(rep(y_i, each = n_draw), nrow = n_draw) - mu_mat

# summarize across draws for each species (column)
summ_one <- function(mat) {
  tibble(
    med = apply(mat, 2, median),
    lwr = apply(mat, 2, quantile, probs = 0.05),
    upr = apply(mat, 2, quantile, probs = 0.95)
  )
}

Est_sum   <- summ_one(Est_contrib_mat)
Ngen_sum  <- summ_one(Ngen_contrib_mat)
Resid_sum <- summ_one(Resid_mat)

net_brms_comp <- net_clean %>%
  transmute(Species, degC, degC.se, estimate, SE, ngen, ngen_centered) %>%
  bind_cols(
    Est_sum   %>% rename(Estimate_Contribution = med,
                         Estimate_LWR = lwr,
                         Estimate_UPR = upr),
    Ngen_sum  %>% rename(Ngen_Contribution = med,
                         Ngen_LWR = lwr,
                         Ngen_UPR = upr),
    Resid_sum %>% rename(Unexplained_Variation = med,
                         Unexplained_LWR = lwr,
                         Unexplained_UPR = upr)
  ) %>%
  mutate(
    # Componentwise medians won't sum exactly to degC; that's expected.
    Check_Sum = Estimate_Contribution + Ngen_Contribution + Unexplained_Variation
  )

message("Summary of (degC - Check_Sum) using component medians:")
print(summary(net_brms_comp$degC - net_brms_comp$Check_Sum))

## ----------------------------
## 6) Figure A: Decomposition stacked bars (same style as yours)
## ----------------------------
out_brms <- net_brms_comp %>%
  select(Species, degC, degC.se,
         Estimate_Contribution, Ngen_Contribution, Unexplained_Variation) %>%
  pivot_longer(
    cols = c(Estimate_Contribution, Ngen_Contribution, Unexplained_Variation),
    names_to = "name", values_to = "value"
  ) %>%
  arrange(degC) %>%
  mutate(Species = fct_reorder(Species, degC))

decomp_brms <- ggplot(out_brms, aes(x = Species, y = value, fill = name)) +
  geom_col(position = "stack", alpha = 1) +
  geom_point(aes(y = degC, shape = "Net change (±SE)"), size = 2) +
  geom_errorbar(aes(ymin = degC - degC.se, ymax = degC + degC.se), width = 0) +
  scale_x_discrete(labels = abbrev_names[as.character(out_brms$Species)]) +
  scale_fill_manual(
    name = "Modeled contributions",
    breaks = c("Estimate_Contribution", "Ngen_Contribution", "Unexplained_Variation"),
    labels = c("Population size", "Symbiont genera", "Unexplained"),
    values = c("Estimate_Contribution" = "#67a9cf",
               "Ngen_Contribution" = "#ef8a62",
               "Unexplained_Variation" = "#C0C0C0")
  ) +
  scale_shape_manual(name = "", values = 16) +
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, size = 10),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(-0.2, "cm"),
    legend.position = c(0.12, 0.86),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  labs(x = "", y = "Change in bleaching threshold (°C)") +
  guides(
    fill = guide_legend(ncol = 1, override.aes = list(shape = NA)),
    shape = guide_legend(order = 1, label.theme = element_text(size = 10))
  )

## ----------------------------
## 7) Figures B & C: BRMS analog partial-residual plots (mirrors your approach)
##
## We use posterior median slopes to define:
##   mu_med = b_est_med * estimate + b_ngen_med * ngen_centered
##   resid  = degC - mu_med
## Then:
##   Adj_degC_Estimate = b_est_med * estimate + resid   (adjusted for ngen=1)
##   Adj_degC_Ngen     = b_ngen_med * ngen_centered + resid (adjusted for estimate=0)
##
## This preserves the same “total residual added back” visual logic you used.
## ----------------------------

b_est_med  <- median(b_est)
b_ngen_med <- median(b_ngen)

net_pr <- net_clean %>%
  mutate(
    mu_med = b_est_med * estimate + b_ngen_med * ngen_centered,
    resid_med = degC - mu_med,
    Adj_degC_Estimate = (b_est_med * estimate) + resid_med,
    Adj_degC_Ngen     = (b_ngen_med * ngen_centered) + resid_med
  )

# jitter ngen to prevent overlap
set.seed(42)
net_pr$Jittered_Ngen <- jitter(net_pr$ngen, amount = 0.15)

# regression line data
estimate_seq <- seq(min(net_pr$estimate), max(net_pr$estimate), length.out = 100)
ngen_seq     <- seq(min(net_pr$ngen),     max(net_pr$ngen),     length.out = 100)

regression_line_estimate_brms <- data.frame(
  estimate = estimate_seq,
  predicted_degC = b_est_med * estimate_seq
)

regression_line_ngen_brms <- data.frame(
  ngen = ngen_seq,
  predicted_degC = b_ngen_med * (ngen_seq - 1)
)

# Plot B: effect of Estimate on adjusted degC (holding ngen=1 baseline)
p1_brms <- ggplot(net_pr, aes(x = estimate, y = Adj_degC_Estimate)) +
  geom_line(
    data = regression_line_estimate_brms,
    aes(x = estimate, y = predicted_degC),
    color = "#67a9cf", linewidth = 1.5, alpha = 1
  ) +
  geom_point(size = 3, alpha = 0.7, color = "#7F7F7F", stroke = 0) +
  geom_errorbar(
    aes(ymin = Adj_degC_Estimate - degC.se, ymax = Adj_degC_Estimate + degC.se),
    width = 0, alpha = 0.7, lwd = 0.3
  ) +
  geom_errorbar(
    aes(xmin = estimate - SE, xmax = estimate + SE),
    width = 0, alpha = 0.7, lwd = 0.3
  ) +
  ggrepel::geom_text_repel(
    aes(label = abbrev_names[as.character(Species)]),
    fontface = "italic", size = 3, seed = 8
  ) +
  scale_x_continuous(
    labels = ~ paste0(round(100 * (exp(.x) - 1), 0), "%"),
    breaks = myticks2
  ) +
  labs(
    x = "Change in population size (2014-2023)",
    y = "Change in bleaching threshold (°C)\n(adjusted for symbiont genera = 1)"
  ) +
  theme_classic()

# Plot C: effect of Ngen on adjusted degC (holding estimate=0 baseline)
p2_brms <- ggplot(net_pr, aes(x = Jittered_Ngen, y = Adj_degC_Ngen)) +
  geom_line(
    data = regression_line_ngen_brms,
    aes(x = ngen, y = predicted_degC),
    color = "#ef8a62", linewidth = 1.5, alpha = 1
  ) +
  geom_point(size = 3, alpha = 0.7, color = "#7F7F7F", stroke = 0) +
  geom_errorbar(
    aes(ymin = Adj_degC_Ngen - degC.se, ymax = Adj_degC_Ngen + degC.se),
    width = 0, alpha = 0.7, lwd = 0.3
  ) +
  ggrepel::geom_text_repel(
    aes(label = abbrev_names[as.character(Species)]),
    fontface = "italic", size = 3, seed = 8
  ) +
  labs(
    x = "Number of symbiont genera",
    y = "Change in bleaching threshold (°C)\n(adjusted for density change = 0)"
  ) +
  theme_classic()

ps_brms <- cowplot::plot_grid(p1_brms, p2_brms, nrow = 1, labels = c("B", "C"))

## ----------------------------
## 8) Combine into the full panel and save (mirrors your s5 object)
## ----------------------------
s5_brms <- cowplot::plot_grid(
  decomp_brms, ps_brms,
  nrow = 2, rel_heights = c(0.55, 0.45),
  labels = "A"
)

s5_brms

# Save
ggsave(filename = "figures/FigS5with24_brms.png",
       plot = s5_brms, width = 190, height = 200, units = "mm")

bayes_R2(fit_me_add)




############################################################
## VARIANCE EXPLAINED BY EACH PREDICTOR (BRMS)
############################################################

draws_beta <- as_draws_df(fit_me_add)

# your slope names (from earlier)
b_est  <- draws_beta[["bsp_meestimateSE"]]
b_ngen <- draws_beta[["b_ngen_centered"]]

# species vectors
x <- net_clean$estimate
z <- net_clean$ngen_centered

n_draw <- length(b_est)

var_density <- numeric(n_draw)
var_ngen    <- numeric(n_draw)
var_total   <- numeric(n_draw)

for(i in seq_len(n_draw)){
  
  d_comp <- b_est[i]  * x
  n_comp <- b_ngen[i] * z
  mu     <- d_comp + n_comp
  
  var_density[i] <- var(d_comp)
  var_ngen[i]    <- var(n_comp)
  var_total[i]   <- var(mu)
}

# proportion of explained variance attributable to each predictor
prop_density <- var_density / var_total
prop_ngen    <- var_ngen / var_total

summary_df <- tibble(
  density = prop_density,
  ngen = prop_ngen
)

summary_df %>%
  summarise(
    density_med = median(density),
    density_lwr = quantile(density,.05),
    density_upr = quantile(density,.95),
    ngen_med    = median(ngen),
    ngen_lwr    = quantile(ngen,.05),
    ngen_upr    = quantile(ngen,.95)
  )



############################################################
## OLS COMPARISON
############################################################

ols_model_adjusted <- lm(
  degC ~ 0 + estimate + ngen_centered,
  data = net_clean,
  weights = 1/degC.se^2
)

summary(ols_model_adjusted)$r.squared

# OLS variance partition
b_est_ols  <- coef(ols_model_adjusted)["estimate"]
b_ngen_ols <- coef(ols_model_adjusted)["ngen_centered"]

d_comp_ols <- b_est_ols  * x
n_comp_ols <- b_ngen_ols * z
mu_ols     <- d_comp_ols + n_comp_ols

var_density_ols <- var(d_comp_ols)
var_ngen_ols    <- var(n_comp_ols)
var_total_ols   <- var(mu_ols)

c(
  density = var_density_ols / var_total_ols,
  ngen    = var_ngen_ols    / var_total_ols
)






############################################################
## POSTERIOR UNIQUE vs SHARED VARIANCE (BRMS)
############################################################

draws_beta <- as_draws_df(fit_me_add)

b_est  <- draws_beta[["bsp_meestimateSE"]]
b_ngen <- draws_beta[["b_ngen_centered"]]

x <- net_clean$estimate
z <- net_clean$ngen_centered
y <- net_clean$degC

n_draw <- length(b_est)

# storage
r2_full   <- numeric(n_draw)
r2_donly  <- numeric(n_draw)
r2_nonly  <- numeric(n_draw)

for(i in seq_len(n_draw)){
  
  mu_full  <- b_est[i]*x + b_ngen[i]*z
  mu_donly <- b_est[i]*x
  mu_nonly <- b_ngen[i]*z
  
  # R2 relative to observed y variance
  r2_full[i]  <- var(mu_full)  / var(y)
  r2_donly[i] <- var(mu_donly) / var(y)
  r2_nonly[i] <- var(mu_nonly) / var(y)
}

# unique contributions (semi-partial R2)
unique_density <- r2_full - r2_nonly
unique_ngen    <- r2_full - r2_donly

# shared variance
shared <- r2_full - unique_density - unique_ngen

variance_partition <- tibble(
  total   = r2_full,
  density_unique = unique_density,
  ngen_unique    = unique_ngen,
  shared  = shared
)

variance_partition %>%
  summarise(
    total_med   = median(total),
    total_lwr   = quantile(total,.05),
    total_upr   = quantile(total,.95),
    
    density_unique_med = median(density_unique),
    density_unique_lwr = quantile(density_unique,.05),
    density_unique_upr = quantile(density_unique,.95),
    
    ngen_unique_med = median(ngen_unique),
    ngen_unique_lwr = quantile(ngen_unique,.05),
    ngen_unique_upr = quantile(ngen_unique,.95),
    
    shared_med = median(shared),
    shared_lwr = quantile(shared,.05),
    shared_upr = quantile(shared,.95)
  ) %>% t()
