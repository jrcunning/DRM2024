ed50s.fun_fast_sp <- function(mod) {
  peak_dates <- get("peak_dates", envir = .GlobalEnv)
  
  # drop 2024
  peak_dates <- peak_dates[peak_dates$year != 2024, ]
  
  b <- fixef(mod)
  mf <- model.frame(mod)
  depth_value <- mean(mf$siteDepth)
  
  ref_year <- levels(mf$year)[1]
  ref_sp   <- levels(mf$Species)[1]
  sp_levels <- levels(mf$Species)
  
  getb <- function(nm) unname(b[nm])
  
  get_interaction <- function(a, b2) {
    n1 <- paste0(a, ":", b2)
    n2 <- paste0(b2, ":", a)
    if (n1 %in% names(b)) unname(b[n1]) else unname(b[n2])
  }
  
  out_list <- vector("list", length(sp_levels))
  
  for (s in seq_along(sp_levels)) {
    sp <- sp_levels[s]
    
    slope <- getb("dhw")
    if (sp != ref_sp) {
      slope <- slope + get_interaction(paste0("Species", sp), "dhw")
    }
    
    eta0 <- numeric(nrow(peak_dates))
    
    for (i in seq_len(nrow(peak_dates))) {
      yr <- as.character(peak_dates$year[i])
      wk <- as.character(peak_dates$Week2[i])
      
      eta <- getb("(Intercept)") + getb("siteDepth") * depth_value
      
      if (sp != ref_sp) {
        eta <- eta + getb(paste0("Species", sp))
      }
      
      if (yr != ref_year) {
        eta <- eta + getb(paste0("year", yr))
      }
      
      if (sp != ref_sp && yr != ref_year) {
        eta <- eta + get_interaction(paste0("Species", sp), paste0("year", yr))
      }
      
      eta <- eta + get_interaction(paste0("year", yr), paste0("Week2", wk))
      
      eta0[i] <- eta
    }
    
    #ed50 <- -eta0 / slope
    # censor at zero
    ed50 <- pmax(-eta0 / slope, 0)
    
    yrs <- as.character(peak_dates$year)
    names(ed50) <- yrs
    
    deltas <- ed50[-1] - ed50[1]
    names(deltas) <- paste0(yrs[-1], "-", yrs[1])
    
    res <- c(ed50, deltas)
    names(res) <- paste0(sp, "|", c(yrs, names(deltas)))
    
    out_list[[s]] <- res
  }
  
  unlist(out_list)
}


bleach50new_sp <- ed50s.fun_fast_sp(dhw.sp.mod.24.depth)
bleach50new_sp
sp.bleach50.24


set.seed(1)
boot.out_sp <- bootMer(
  dhw.sp.mod.24.depth,
  FUN = ed50s.fun_fast_sp,
  nsim = 24,
  re.form = NULL,
  seed = 123,
  parallel = "multicore", 
  ncpus = 8
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

sp.ed50conf2 <- ci_sp %>%
  filter(grepl("2023-2014", term)) %>%
  separate(term, into = c("Species", "year"), sep = "\\|") %>%
  mutate(Species = fct_reorder(Species, est),
         ed50 = est)

ggplot(sp.ed50conf2, aes(x = Species, y = ed50)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) + ylim(0,20)

ggplot(sp.ed50conf.24 %>% filter(grepl("2023-2014", year)) %>% mutate(Species = fct_reorder(Species, ed50)), 
       aes(x = Species, y = ed50)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper))+ ylim(0,20)
