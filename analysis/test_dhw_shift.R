##### 1. ORIGINAL MODEL
# Model differences in bleaching severity between years with DHW as continuous predictor
dhw.mod24 <- glmer(cbind(BL, NB) ~ dhw + year + year:Week2 + (dhw + year | Species) + (year | Subregion),  
                   family = "binomial", data = allg, 
                   verbose = FALSE, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
summary(dhw.mod24)

# dhw.mod24.depth <- glmer(cbind(BL, NB) ~ siteDepth + dhw + year + year:Week2 + 
#                            (dhw + year | Species) + (year | Subregion),  
#                  family = "binomial", data = allg, 
#                  verbose = FALSE, nAGQ = 0, control = glmerControl(optimizer = "nloptwrap"))
# 
# anova(dhw.mod24, dhw.mod24.depth)
# summary(dhw.mod24.depth)
# emmeans(dhw.mod24, specs = c("year", "dhw"), type = "response")
# emmeans(dhw.mod24.depth, specs = c("year", "dhw"), type = "response")

# Function to get ED50s and ∆ED50s
ed50s.fun_fast <- function(mod, peak_dates = NULL) {
  if (is.null(peak_dates)) {
    if (!exists("peak_dates", envir = .GlobalEnv, inherits = FALSE)) {
      stop("No peak_dates provided and no object named 'peak_dates' found in the global environment.")
    }
    peak_dates <- get("peak_dates", envir = .GlobalEnv)
  }
  
  combos <- peak_dates
  b <- fixef(mod)
  
  mf <- model.frame(mod)
  
  # ensure year/Week2 are factors with the model's levels (this will also catch fakes)
  combos$year  <- factor(as.character(combos$year),  levels = levels(mf$year))
  combos$Week2 <- factor(as.character(combos$Week2), levels = levels(mf$Week2))
  
  # explicit checks (friendlier error messages than "new levels" downstream)
  bad_year  <- setdiff(as.character(peak_dates$year),  levels(mf$year))
  bad_week2 <- setdiff(as.character(peak_dates$Week2), levels(mf$Week2))
  if (length(bad_year))  stop("peak_dates has year(s) not in model: ",  paste(bad_year, collapse = ", "))
  if (length(bad_week2)) stop("peak_dates has Week2 level(s) not in model: ", paste(bad_week2, collapse = ", "))
  
  # build a fixed-effect design matrix for dhw = 0 rows
  combos$dhw <- 0
  tt <- delete.response(terms(mod))
  X <- model.matrix(tt, combos)
  
  # handle rank-deficiency: columns in X that were dropped from fixef(mod)
  miss <- setdiff(colnames(X), names(b))
  if (length(miss)) {
    warning("Some fixed-effect columns are not in fixef(mod) (rank deficiency); treating as 0: ",
            paste(miss, collapse = ", "))
    b <- c(b, setNames(rep(0, length(miss)), miss))
  }
  
  eta0 <- drop(X %*% b[colnames(X)])
  
  # ED50 where eta = 0 on logit scale: dhw = -eta0 / beta_dhw
  if (!("dhw" %in% names(b)) || abs(unname(b["dhw"])) < 1e-12) {
    ed50 <- rep(NA_real_, length(eta0))
  } else {
    ed50 <- -eta0 / unname(b["dhw"])
  }
  
  ed50_names <- as.character(combos$year)
  delta <- ed50[-1] - ed50[1]
  delta_names <- paste0(as.character(combos$year[-1]), "-", as.character(combos$year[1]))
  
  res <- c(ed50, delta)
  setNames(res, c(ed50_names, delta_names))
}

# Get point estimates for ED50s and ∆ED50s
bleach50 <- ed50s.fun_fast(dhw.mod24)
bleach50

# Get bootstrap confidence intervals for ED50s and ∆ED50s across years
# set.seed(1)
# boot.outnull <- bootMer(dhw.mod24, FUN = ed50s.fun_fast, nsim = 100, re.form = NULL,
#                     seed = 123, parallel = "multicore", ncpus = 10) 
# saveRDS(boot.outnull, file = "data/processed/boot.outnull.rds")
boot.outnull <- readRDS(file = "data/processed/boot.outnull.rds")
## Get confidence intervals from Bootstrap
lower <- apply(boot.outnull$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
upper <- apply(boot.outnull$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
se <- apply(boot.outnull$t, 2, sd)

ed50conf <- tibble(
  year = names(bleach50), 
  ed50 = bleach50,
  se = se,
  lower = lower,
  upper = upper)


# Get max and min DHW for conducted surveys as bounds
dhw_ranges <- all0 %>% 
  group_by(year) %>%
  summarize(maxdhw = max(dhw),
            mindhw = min(dhw))

# Target combos + bounds
targets <- peak_dates %>%
  transmute(year = as.character(year), Week2) %>%   # match emmeans labels
  left_join(dhw_ranges, by = "year")

# Get emmeans grid # For plotting, Get emmeans for just the specific 2-week windows of interest in each year (corresp. to max. bleaching), And select just the range of peak DHWs experienced across sites in each year
dhw.res <- emmeans(dhw.mod24, specs = c("dhw", "year", "Week2"), type = "response",
                   re.form = NA,
                   at = list(dhw = seq(0, 23, 0.1), Week2 = unique(peak_dates$Week2)), 
                   rg.limit = 20000, level = 0.84)

# Filter for target year:Week2 and dhw ranges
dhw.res2 <- as_tibble(dhw.res) %>%
  inner_join(targets, by = c("year", "Week2")) %>%
  filter(dhw > mindhw, dhw <= maxdhw)

dhw.res3 <- as_tibble(dhw.res) %>%
  inner_join(targets, by = c("year", "Week2")) %>%
  filter(year == 2024, dhw > maxdhw, prob <= 0.5)


# Plot
dhw.bleach.plot <- ggplot(as.tibble(dhw.res2), aes(x = dhw, y = prob, group = year)) +
  # Add lines and conf intervals for each year
  geom_ribbon(aes(ymin = prob - SE, ymax = prob + SE), lwd = 0, alpha = 0.1) +            # WHY NOT SE?
  geom_line(aes(color = prob), lwd = 2) +
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
  scale_color_gradient2(low = "forestgreen", mid = "yellow", high = "firebrick1", 
                        limits = c(0, 1), midpoint = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Degree Heating Weeks (°C-weeks)", y = "Bleaching prevalence") +
  coord_cartesian(clip = "off")

dhw.bleach.plot



######## delta ED50s from ORIGINAL MODEL
ed50conf




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
  left_join(drmsite_dhw_shift %>% dplyr::select(year, Site, dhw, pixel, date), 
            by = c("year", "Site"))

allg_shift %>% filter(is.na(dhw))

# Clean data
allg_shift2 <- allg_shift %>%
  dplyr::select(BL, NB, dhw, year, Week2, Species, Subregion, pixel, date) %>%
  tidyr::drop_na()




# Validation model
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
# boot.shift <- bootMer(dhw.mod24.shift, FUN = ed50s.fun_fast, nsim = 100, re.form = NULL,
#                     seed = 123, parallel = "multicore", ncpus = 20)
# saveRDS(boot.shift, file = "data/processed/boot.shift.rds")
boot.shift <- readRDS(file = "data/processed/boot.shift.rds")
## Get confidence intervals from Bootstrap
lower <- apply(boot.shift$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
upper <- apply(boot.shift$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
se <- apply(boot.shift$t, 2, sd)

ed50conf_shift <- tibble(
  year = names(bleach50), 
  ed50 = bleach50,
  se = se,
  lower = lower,
  upper = upper)
ed50conf_shift


# Get max and min DHW for conducted surveys as bounds
dhw_ranges_shift <- allg_shift2 %>% 
  group_by(year) %>%
  summarize(maxdhw = max(dhw),
            mindhw = min(dhw))

# Target combos + bounds
targets <- peak_dates_shift %>%
  transmute(year = as.character(year), Week2) %>%   # match emmeans labels
  left_join(dhw_ranges_shift, by = "year")

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
  geom_line(data = dhw.res2, linetype = 2, linewidth = 0.1) +
  # Add year text labels
  geom_text(data = dhw_ranges_shift,
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



# 2024 has no DHWS! remove?
allg_shift_no2024 <- allg_shift2 %>% filter(year != "2024")

dhw.mod24.shift_no2024 <- glmer(
  cbind(BL, NB) ~ dhw + year + year:Week2 +
    (dhw + year | Species) + (year | Subregion),
  family = "binomial",
  data = allg_shift_no2024,
  control = glmerControl(optimizer = "nloptwrap"),
  nAGQ = 0
)
summary(dhw.mod24.shift_no2024)

peak_dates_no2024 <- peak_dates_shift %>% filter(year != 2024)

bleach50_no2024 <- ed50s.fun_fast(dhw.mod24.shift_no2024, peak_dates = peak_dates_no2024)
bleach50_no2024


set.seed(123)
boot_no2024 <- bootMer(
  dhw.mod24.shift_no2024,
  FUN = function(m) ed50s.fun_fast(m, peak_dates = peak_dates_no2024),
  nsim = 200,
  re.form = NULL,
  parallel = "multicore",
  ncpus = 20
)

ci_no2024 <- tibble(
  term  = names(bleach50_no2024),
  est   = as.numeric(bleach50_no2024),
  lower = apply(boot_no2024$t, 2, quantile, probs = 0.025, na.rm = TRUE),
  upper = apply(boot_no2024$t, 2, quantile, probs = 0.975, na.rm = TRUE),
  se    = apply(boot_no2024$t, 2, sd, na.rm = TRUE)
)
ci_no2024



### ANOTHER TEST 

# survey observations used in the bleaching model
# needs: pixel, date, year, and baseline dhw used in model (call it dhw0)
obs <- allg_shift2 %>%
  transmute(
    year  = as.integer(as.character(year)),
    pixel,
    date,
    dhw0 = dhw
  )

adds <- seq(0, 2.5, 0.1)

obs_expanded <- obs %>%
  tidyr::crossing(add = adds)

obs_diffs <- obs_expanded %>%
  left_join(dhw_by_add, by = c("pixel","year","date","add")) %>%
  filter(!is.na(dhw)) %>%
  mutate(diff_obs = dhw0 - dhw)

summary(obs_diffs$diff_obs)
sum(obs_diffs$diff_obs < -1e-8, na.rm = TRUE)

# temps must have: pixel, year, date (or yday), temp
# mmm is scalar
recompute_dhw_for_add <- function(add) {
  temps %>%
    mutate(year = as.integer(year),
           mmm_adj = mmm + add) %>%
    arrange(pixel, year, date) %>%
    group_by(pixel, year) %>%
    mutate(
      dhd = if_else(temp > (mmm_adj + 1), (temp - mmm_adj) / 7, 0),
      dhw = slide_index_sum(dhd, date, before = 84L, complete = FALSE)
    ) %>%
    ungroup() %>%
    dplyr::select(pixel, year, date, dhw) %>%
    mutate(add = add)
}

dhw_by_add <- map_dfr(adds, recompute_dhw_for_add)

# Join recomputed DHW to the survey observations at the survey date/pixel
obs_diffs <- obs %>%
  left_join(dhw_by_add, by = c("pixel", "year", "date")) %>%
  filter(!is.na(dhw)) %>%
  mutate(diff_obs = dhw0 - dhw)  # decrease in DHW at survey obs

# Build year-specific curve of mean decrease in "model-used DHW"
curve_year <- obs_diffs %>%
  group_by(year, add) %>%
  summarise(diff_obs = mean(diff_obs, na.rm = TRUE), .groups = "drop")


invert_add <- function(df_year, target) {
  df_year <- df_year %>% arrange(add)
  # solve for add such that diff_obs(add) == target
  approx(x = df_year$diff_obs, y = df_year$add, xout = target, rule = 1)$y
}

target_2015 <- 5.0
target_2023 <- 8.5

degC_2015 <- invert_add(filter(curve_year, year == 2015), target_2015)
degC_2023 <- invert_add(filter(curve_year, year == 2023), target_2023)

degC_2015; degC_2023
