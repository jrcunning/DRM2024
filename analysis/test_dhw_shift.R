# Recompute 2015 and 2023 DHWs at calculated shifted thresholds
temps <- read_csv("data/processed/sst_reef_pixels.csv") %>%
  filter(year %in% c(2014, 2015, 2023, 2024))

temps

degC_year <- degC %>% 
  dplyr::select(year, degC) %>%
  tibble::add_row(year = "2014", degC = 0) %>%
  mutate(year = as.integer(year))

mmm <- 29.6264 # Actual MMM for FL Keys

dhw_pix_yearshift <- temps %>%
  mutate(year = as.integer(year)) %>%
  left_join(degC_year, by = "year") %>%
  mutate(
    degC = coalesce(degC, 0),
    mmm_adj = mmm + degC
  ) %>%
  arrange(pixel, year, yday) %>%
  group_by(pixel, year) %>%
  nest() %>%
  mutate(
    dhd = map(data, ~ if_else(.x$temp > (.x$mmm_adj + 1),
                              (.x$temp - .x$mmm_adj) / 7,
                              0)),
    dhw = map2(data, dhd, ~ slider::slide_index_sum(.y, .x$date, before = 84L)),
    maxdhw = map_dbl(dhw, ~ max(.x, na.rm = TRUE))
  ) %>%
  ungroup()

dhw_long_yearshift <- dhw_pix_yearshift %>%
  transmute(
    pixel, year,
    date = map(data, ~ .x$date),
    dhw  = dhw
  ) %>%
  unnest(c(date, dhw))

dhw_long_yearshift %>%
  mutate(yday = yday(date)) %>%
  ggplot(aes(x = yday, y = dhw)) +
  facet_wrap(~year) +
  geom_line(aes(group = pixel))

#### Compute peak dates of heat stress from shifted DHWs

# ============================================================
# 0) Collapse pixel-level shifted DHW to daily mean DHW per year
# ============================================================

dhw_daily_shift <- dhw_long_yearshift %>%
  mutate(
    year = as.integer(year),
    date = as.Date(date)   # if you only have yday, tell me and I’ll derive date
  ) %>%
  filter(!is.na(date), !is.na(dhw)) %>%
  group_by(year, date) %>%
  summarise(
    mean_dhw = mean(dhw, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 1) Smoothed daily averages (same as your workflow)
# ============================================================

dhw_daily_shift <- dhw_daily_shift %>%
  group_by(year) %>%
  arrange(date) %>%
  mutate(
    date_x = make_date(2000, month(date), mday(date)),
    mean_dhw_smooth = zoo::rollmean(mean_dhw, k = 7, fill = NA)
  ) %>%
  ungroup()

# ============================================================
# 2) Add Week2 bins + dhw_used + numeric time
# ============================================================

dhw_daily_shift <- dhw_daily_shift %>%
  mutate(
    week  = isoweek(date),
    Week2 = cut(week, breaks = seq(1, 53, by = 2)),   # matches your definition
    dhw_used = coalesce(mean_dhw_smooth, mean_dhw),
    t_num = as.numeric(date)
  ) %>%
  filter(!is.na(dhw_used))

# ============================================================
# 3) Rising-limb parameters
# ============================================================

p_lo <- 0.05
p_hi <- 0.60

# ============================================================
# 4) Function: fit rising limb and compute intersection with annual max
# ============================================================

peak_extrapolated <- function(d) {
  d <- d %>% arrange(date)
  dhw_max <- max(d$dhw_used, na.rm = TRUE)
  
  d_peak <- d %>% slice_max(dhw_used, n = 1, with_ties = FALSE)
  peak_date_obs <- d_peak$date
  
  d_rise <- d %>%
    filter(date <= peak_date_obs) %>%
    filter(dhw_used >= p_lo * dhw_max, dhw_used <= p_hi * dhw_max)
  
  if (nrow(d_rise) < 5) {
    return(tibble(
      year = unique(d$year),
      dhw_max = dhw_max,
      peak_date_obs = peak_date_obs,
      peak_date_pred = as.Date(NA),
      Week2_pred = factor(NA),
      n_rise = nrow(d_rise)
    ))
  }
  
  mod <- lm(dhw_used ~ t_num, data = d_rise)
  a <- unname(coef(mod)[["(Intercept)"]])
  b <- unname(coef(mod)[["t_num"]])
  
  t_star <- (dhw_max - a) / b
  peak_date_pred <- as.Date(t_star, origin = "1970-01-01")
  
  tibble(
    year = unique(d$year),
    dhw_max = dhw_max,
    peak_date_obs = peak_date_obs,
    peak_date_pred = peak_date_pred,
    Week2_pred = cut(isoweek(peak_date_pred), breaks = seq(1, 52, by = 2)),
    n_rise = nrow(d_rise)
  )
}

# ============================================================
# 5) Compute peak dates (shifted DHW)
# ============================================================

peak_dates_shift <- dhw_daily_shift %>%
  group_by(year) %>%
  group_modify(~ peak_extrapolated(.x)) %>%
  ungroup() %>%
  dplyr::select(year, peak_date_pred, Week2 = Week2_pred, dhw_max, peak_date_obs, n_rise)

peak_dates_shift
peak_dates_shift[which(peak_dates_shift$year == "2024"), "Week2"] <- "(37,39]"

# ============================================================
# 6) Build fitted line on the rising-limb subset for plotting
# ============================================================

rise_fit_df_shift <- dhw_daily_shift %>%
  inner_join(
    peak_dates_shift %>% dplyr::select(year, dhw_max, peak_date_obs, peak_date_pred),
    by = "year"
  ) %>%
  group_by(year) %>%
  group_modify(~{
    d <- .x %>% arrange(date)
    d_rise <- d %>%
      filter(date <= unique(d$peak_date_obs)) %>%
      filter(dhw_used >= p_lo * unique(d$dhw_max),
             dhw_used <= p_hi * unique(d$dhw_max))
    mod <- lm(dhw_used ~ t_num, data = d_rise)
    d_rise %>% mutate(dhw_fit = predict(mod, newdata = d_rise))
  }) %>%
  ungroup()

# ============================================================
# 7) Plot
# ============================================================

ggplot(dhw_daily_shift, aes(date, dhw_used)) +
  geom_line() +
  facet_wrap(~ year, scales = "free_x", ncol = 2) +
  geom_line(data = rise_fit_df_shift, aes(y = dhw_fit), linewidth = 1) +
  geom_hline(data = peak_dates_shift, aes(yintercept = dhw_max), linetype = "dashed") +
  geom_point(data = peak_dates_shift, aes(x = peak_date_pred, y = dhw_max), size = 2.5) +
  geom_vline(data = peak_dates_shift,
             aes(xintercept = as.numeric(peak_date_pred)),
             linetype = "dotted") +
  labs(
    x = NULL, y = "Shifted DHW",
    title = "Peak heat stress timing from extrapolated rising-limb intersection",
    subtitle = sprintf("Rising limb defined as %.0f–%.0f%% of annual max shifted DHW (pre-peak)",
                       100*p_lo, 100*p_hi)
  ) +
  theme_bw()





# Now pull the shifted DHWs for all DRM sites on date they were surveyed

# Get all drm sites lat/lon for 2014/2015/2023/2024
drmsites <- allf %>%
  distinct(Site, Date, Subregion, Longitude, Latitude) %>%
  dplyr::select(Site, date = Date, Subregion, Longitude, Latitude)

# Map lat/lon to pixel ID using reef_pixel_grid
# grid defining pixels (geometry + CRS)
reef_pixel_grid <- rast("data/spatial/reef_pixel_grid.tif")

# make a SpatVector of points in WGS84
pts <- vect(drmsites, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# project points to the raster CRS
pts_r <- project(pts, crs(reef_pixel_grid))

# get raster cell index for each point
drmsites$pixel <- cellFromXY(reef_pixel_grid, geom(pts_r)[, c("x", "y")])

# optional: flag points outside the grid
drmsites <- drmsites %>% mutate(pixel_in_grid = !is.na(pixel))
summary(drmsites$pixel_in_grid)



# Pull shifted DHWs for each survey based on Pixel and Date
drmsite_dhw_shift <- drmsites %>%
  left_join(dhw_long_yearshift, by = c("pixel", "date")) %>%
  mutate(year = as.character(year))
drmsite_dhw_shift %>% filter(is.na(dhw))


# Join with actual DRM survey data!
allg_shift <- allg %>% 
  dplyr::select(year, Week2, Subregion, Site, Species, BL, NB) %>%
  left_join(drmsite_dhw_shift %>% dplyr::select(year, Site, dhw), 
            by = c("year", "Site"))

allg_shift %>% filter(is.na(dhw))

# Clean data
allg_shift2 <- allg_shift %>%
  dplyr::select(BL, NB, dhw, year, Week2, Species, Subregion) %>%
  tidyr::drop_na()





# Model differences in bleaching severity between years with DHW as continuous predictor
# Model differences in bleaching severity between years with DHW as continuous predictor
dhw.mod24.shift <- glmer(cbind(BL, NB) ~ dhw + year + year:Week2 + (dhw + year | Species) + (year | Subregion),  
                   family = "binomial", data = allg_shift2, 
                   verbose = FALSE, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
summary(dhw.mod24.shift)



# Get point estimates for ED50s and ∆ED50s
bleach50 <- ed50s.fun_fast(dhw.mod24.shift)
bleach50

# Get bootstrap confidence intervals for ED50s and ∆ED50s across years
# set.seed(1)
boot.shift <- bootMer(dhw.mod24.shift, FUN = ed50s.fun_fast, nsim = 100, re.form = NULL,
                    seed = 123, parallel = "multicore", ncpus = 20)
saveRDS(boot.shift, file = "data/processed/boot.shift.rds")
boot.shift <- readRDS(file = "data/processed/boot.shift.rds")
## Get confidence intervals from Bootstrap
lower <- apply(boot.shift$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
upper <- apply(boot.shift$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
se <- apply(boot.shift$t, 2, sd)

ed50conf <- tibble(
  year = names(bleach50), 
  ed50 = bleach50,
  se = se,
  lower = lower,
  upper = upper)
ed50conf


# Get max and min DHW for conducted surveys as bounds
dhw_ranges <- allg_shift2 %>% 
  group_by(year) %>%
  summarize(maxdhw = max(dhw),
            mindhw = min(dhw))

# Target combos + bounds
targets <- peak_dates_shift %>%
  transmute(year = as.character(year), Week2) %>%   # match emmeans labels
  left_join(dhw_ranges, by = "year")

# Get emmeans grid # For plotting, Get emmeans for just the specific 2-week windows of interest in each year (corresp. to max. bleaching), And select just the range of peak DHWs experienced across sites in each year
dhw.res <- emmeans(dhw.mod24.shift, specs = c("dhw", "year", "Week2"), type = "response",
                   re.form = NA,
                   at = list(dhw = seq(0, 12, 0.1), Week2 = unique(peak_dates_shift$Week2)), 
                   rg.limit = 20000, level = 0.84)

# Filter for target year:Week2 and dhw ranges
dhw.res2 <- as_tibble(dhw.res) %>%
  inner_join(targets, by = c("year", "Week2")) %>%
  filter(dhw > mindhw, dhw <= maxdhw)


# Plot
dhw.shift.bleach.plot <- ggplot(as.tibble(dhw.res2), aes(x = dhw, y = prob, group = year)) +
  # Add lines and conf intervals for each year
  geom_ribbon(aes(ymin = prob - SE, ymax = prob + SE), lwd = 0, alpha = 0.1) +            # WHY NOT SE?
  geom_line(aes(color = year), lwd = 2) +
  geom_line(data = dhw.res3, linetype = 2, linewidth = 0.1) +
  # Add year text labels
  geom_text(data = dhw_ranges,
            aes(x = maxdhw, y = c(0.80, 0.62, 0.96, 0.28), label = year),
            hjust = c(0,0,1,0), nudge_y = c(0,0,0,0), nudge_x = c(0,0,0,-0.3), size = 3) +
  # Add line segments at ED50s
  geom_segment(data = slice(ed50conf, 1:4),
               aes(x = ed50, xend = ed50, y = 0.492, yend = -0.025), lty = c(1,1,1,2), lwd = 0.1) +
  geom_point(data = slice(ed50conf, 1:4), aes(y = 0.5, x = ed50), pch = 1, stroke = 0.2) +
  # Annotate increase in heat tolerance
  annotate("segment", x = ed50conf[[1,2]], xend = ed50conf[[3,2]], y = 0, lwd = 0.2,
           arrow = arrow(type = "closed", length = unit(2, "mm"))) +
  annotate("segment", x = ed50conf[[3,2]], xend = ed50conf[[4,2]], y = 0, lwd = 0.2, lty = 2,
           arrow = arrow(type = "open", length = unit(2, "mm"))) +
  annotate("segment", x = ed50conf[[4,2]] - 1e-6, xend = ed50conf[[4,2]], y = 0, yend = 0, lwd = 0.2, lty = 1,
           arrow = arrow(type = "open", length = unit(2, "mm"))) +
  scale_x_continuous(limits = c(-1, 22), breaks = seq(0, 24, 4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.025, 1), expand = expansion(mult = c(0, 0)), breaks = seq(0, 1, 0.1), 
                     labels = scales::label_percent()) +
  theme_classic() +
  labs(x = "Adjusted Degree Heating Weeks (°C-weeks)", y = "Bleaching prevalence") +
  coord_cartesian(clip = "off")

dhw.shift.bleach.plot













########## what if remove week2 from original model
# Model differences in bleaching severity between years with DHW as continuous predictor
dhw.mod24 <- glmer(cbind(BL, NB) ~ dhw + year + (dhw + year | Species),  
                   family = "binomial", data = allg, 
                   verbose = FALSE, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
summary(dhw.mod24)

# --- build a prediction grid over DHW and year (population-averaged = fixed effects only) ---
# Make sure year is a factor in the model data (matches model encoding)
allg2 <- allg %>%
  mutate(year = factor(year))

year_vals <- levels(allg2$year)

dhw_grid <- seq(min(allg2$dhw, na.rm = TRUE),
                max(allg2$dhw, na.rm = TRUE),
                length.out = 200)

pred_df <- tidyr::expand_grid(
  year = year_vals,
  dhw  = dhw_grid
) %>%
  mutate(
    year = factor(year, levels = year_vals),
    p_bleach = predict(dhw.mod24, newdata = ., type = "response", re.form = NA),
    pct_bleach = 100 * p_bleach
  )

ggplot() +
  geom_point(
    data = allg2 %>% mutate(pct_obs = 100 * BL / (BL + NB)),
    aes(x = dhw, y = pct_obs, color = year),
    alpha = 0.15
  ) +
  geom_line(
    data = pred_df,
    aes(x = dhw, y = pct_bleach, color = year),
    linewidth = 1.2
  ) +
  labs(x = "DHW", y = "Percent bleached", color = "Year") +
  theme_bw()












########## what if remove week2 from original SPECIES model
# Model differences in bleaching severity between years with DHW as continuous predictor
dhw.mod24.sp <- glm(cbind(BL, NB) ~ Species + dhw + Species:dhw + year + Species:year,  
                   family = "binomial", data = allg)
summary(dhw.mod24.sp)


# Ensure factors match model encoding
allg2 <- allg %>%
  mutate(
    Species = factor(Species),
    year    = factor(year)
  )

# Prediction grid
dhw_grid <- seq(min(allg2$dhw, na.rm = TRUE),
                max(allg2$dhw, na.rm = TRUE),
                length.out = 200)

pred_df <- tidyr::expand_grid(
  Species = levels(allg2$Species),
  year    = levels(allg2$year),
  dhw     = dhw_grid
) %>%
  mutate(
    Species = factor(Species, levels = levels(allg2$Species)),
    year    = factor(year,    levels = levels(allg2$year)),
    p_bleach = predict(dhw.mod24.sp, newdata = ., type = "response"),
    pct_bleach = 100 * p_bleach
  )

# Optional raw points (can be heavy; downsample if needed)
obs_df <- allg2 %>%
  mutate(pct_obs = 100 * BL / (BL + NB))

ggplot() +
  geom_point(
    data = obs_df,
    aes(x = dhw, y = pct_obs, color = year),
    alpha = 0.12, size = 0.8
  ) +
  geom_line(
    data = pred_df,
    aes(x = dhw, y = pct_bleach, color = year),
    linewidth = 1
  ) +
  facet_wrap(~ Species, scales = "free_y") +
  labs(x = "DHW", y = "Percent bleached", color = "Year") +
  theme_bw()
