ed50s.fun_fast_sp <- function(mod) {
  
  combos <- peak_dates  # same object you used before: year, Week2 rows (e.g., 4 rows)
  b <- fixef(mod)
  
  # helper to safely grab coefficients that may not exist (rank deficiency)
  getb <- function(nm) if (nm %in% names(b)) unname(b[nm]) else 0
  
  # levels in the model (to identify reference levels robustly)
  mf <- model.frame(mod)
  ref_year <- levels(mf$year)[1]
  ref_sp   <- levels(mf$Species)[1]
  
  # all species we will compute
  sp_levels <- levels(mf$Species)
  
  # storage
  out_list <- vector("list", length(sp_levels))
  
  for (s in seq_along(sp_levels)) {
    sp <- sp_levels[s]
    
    ## slope for this species (dhw + Species:dhw if not reference)
    slope <- getb("dhw")
    if (sp != ref_sp) {
      # interaction name is usually "Species<sp>:dhw" (not "dhw:Species<sp>")
      slope <- slope + getb(paste0("Species", sp, ":dhw"))
    }
    
    ## eta0 for each year/week combo at dhw = 0
    eta0 <- numeric(nrow(combos))
    for (i in seq_len(nrow(combos))) {
      yr <- as.character(combos$year[i])
      wk <- as.character(combos$Week2[i])
      
      eta <- getb("(Intercept)")
      
      # Species main effect (intercept offset) if not reference
      if (sp != ref_sp) eta <- eta + getb(paste0("Species", sp))
      
      # year main effect if not reference year
      if (yr != ref_year) eta <- eta + getb(paste0("year", yr))
      
      # Species:year interaction (treatment-coded; present only for non-ref sp and non-ref year)
      if (sp != ref_sp && yr != ref_year) {
        eta <- eta + getb(paste0("Species", sp, ":year", yr))
      }
      
      # year:Week2 interaction (you previously noted even ref year can appear here)
      eta <- eta + getb(paste0("year", yr, ":Week2", wk))
      
      eta0[i] <- eta
    }
    
    # ED50 per combo
    ed50 <- if (abs(slope) < 1e-12) rep(NA_real_, length(eta0)) else -eta0 / slope
    
    # name ED50s by year (assumes each year appears once in combos)
    yrs <- as.character(combos$year)
    names(ed50) <- yrs
    
    # deltas vs first combo/year (exactly your prior pattern)
    ed50_ref <- ed50[1]
    d1 <- ed50[2] - ed50_ref
    d2 <- ed50[3] - ed50_ref
    d3 <- ed50[4] - ed50_ref
    
    res <- c(ed50, d1, d2, d3)
    
    nm <- c(yrs,
            paste0(yrs[2], "-", yrs[1]),
            paste0(yrs[3], "-", yrs[1]),
            paste0(yrs[4], "-", yrs[1]))
    
    out_list[[s]] <- setNames(res, paste0(sp, "|", nm))
  }
  
  unlist(out_list)
}

bleach50_sp <- ed50s.fun_fast_sp(dhw.sp.mod)
bleach50_sp



set.seed(1)
boot.out_sp <- bootMer(
  dhw.sp.mod,
  FUN = ed50s.fun_fast_sp,
  nsim = 100,
  re.form = NULL,
  seed = 123,
  parallel = "multicore",
  ncpus = 10
)

# CIs
lower <- apply(boot.out_sp$t, 2, quantile, probs = 0.025, na.rm = TRUE)
upper <- apply(boot.out_sp$t, 2, quantile, probs = 0.975, na.rm = TRUE)
se    <- apply(boot.out_sp$t, 2, sd, na.rm = TRUE)

ci_sp <- tibble::tibble(
  term  = names(boot.out_sp$t0),
  est   = unname(boot.out_sp$t0),
  se    = unname(se),
  lower = unname(lower),
  upper = unname(upper)
)

sp.ed50conf <- ci_sp %>%
  filter(grepl("2023-2014", term)) %>%
  separate(term, into = c("Species", "year"), sep = "\\|") %>%
  mutate(species = fct_reorder(Species, est))

ggplot(sp.ed50conf, aes(x = species, y = est)) +
  geom_point()
