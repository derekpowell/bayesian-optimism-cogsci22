# Reliability of person-level updating parameters in the actuarial risk task,
# as a function of the number of items per participant.
#
# Data are generated from the two-part model fit to the 2021 data (M5 in irt_optimism.Rmd):
#   - update at all?  logit P = 0.797 - 0.178*g + u_i + v_j
#   - size of update  logit(post) = logit(pre) + (alpha + a_i + a_j + (beta + b_i + b_j) * g) * E
# with person and item effects drawn from M5's estimated SDs and correlations.
# Each simulated item borrows the empirical (pre, E, g) distribution of one of the 2021 items,
# so good/bad news stays item-driven unless the "balanced" design reassigns it.
#
# Reliability = squared correlation between true and estimated person parameters (classical test theory).
# Models are fit with glmmTMB (maximum likelihood) for speed; estimates are conditional modes.

suppressMessages({
  library(tidyverse)
  library(glmmTMB)
  library(MASS, include.only = "mvrnorm")
  library(parallel)
})

load("processed_data.Rdata")
trials_2021 <- df_reg %>%
  ungroup() %>%
  transmute(item_src = item_num, lp = logit_pre, E = logit_truebase - logit_base, g = good_news)

informative_items <- trials_2021 %>%
  group_by(item_src) %>% summarize(mE2 = mean(E^2)) %>%
  slice_max(mE2, n = 10) %>% pull(item_src)

# --- generating parameters (posterior means from M5) ---
pars <- list(
  upd0 = 0.797, upd_g = -0.178, alpha = 0.716, beta = -0.083, phi = 17,
  person_sd = c(1.287, 0.208, 0.214),
  person_cor = matrix(c(1, .18, .14, .18, 1, -.69, .14, -.69, 1), 3),
  item_sd = c(1.048, 0.118, 0.103),
  item_cor = matrix(c(1, .24, .29, .24, 1, -.52, .29, -.52, 1), 3)
)
covmat <- function(sd, cor) diag(sd) %*% cor %*% diag(sd)
clamp <- function(p, lo = .005) pmin(pmax(p, lo), 1 - lo)

simulate_data <- function(n_subj, n_items, design) {
  people <- mvrnorm(n_subj, rep(0, 3), covmat(pars$person_sd, pars$person_cor))
  items <- mvrnorm(n_items, rep(0, 3), covmat(pars$item_sd, pars$item_cor))
  # "informative" draws only from the 10 items with the largest mean squared evidence term in 2021
  src_pool <- if (design == "informative") informative_items else unique(trials_2021$item_src)
  src <- sample(src_pool, n_items, replace = n_items > length(src_pool))

  d <- expand_grid(subj = seq_len(n_subj), item = seq_len(n_items)) %>%
    mutate(item_src = src[item])
  # draw each trial's (pre, E, g) from the 2021 trials for that source item
  pool <- split(seq_len(nrow(trials_2021)), trials_2021$item_src)
  rows <- vapply(as.character(d$item_src), function(s) pool[[s]][sample.int(length(pool[[s]]), 1)], integer(1))
  d <- bind_cols(d, trials_2021[rows, c("lp", "E", "g")])
  if (design == "balanced") d$g <- rbinom(nrow(d), 1, .5)

  u <- people[d$subj, 1]; a_i <- people[d$subj, 2]; b_i <- people[d$subj, 3]
  v <- items[d$item, 1]; a_j <- items[d$item, 2]; b_j <- items[d$item, 3]
  d %>% mutate(
    upd = rbinom(n(), 1, plogis(pars$upd0 + pars$upd_g * g + u + v)),
    mu = plogis(lp + (pars$alpha + a_i + a_j + (pars$beta + b_i + b_j) * g) * E),
    y = if_else(upd == 1, clamp(rbeta(n(), mu * pars$phi, (1 - mu) * pars$phi)), clamp(plogis(lp))),
    subj = factor(subj), item = factor(item)
  ) %>%
    structure(truth = tibble(subj = factor(seq_len(n_subj)), u = people[, 1],
                             alpha = pars$alpha + people[, 2], beta = pars$beta + people[, 3],
                             # overall evidence weight including no-update trials, which one-part models estimate
                             alpha_overall = (pars$alpha + people[, 2]) * plogis(pars$upd0 + people[, 1])))
}

person_estimates <- function(fit, terms) {
  fe <- fixef(fit)$cond
  re <- ranef(fit)$cond$subj
  as_tibble(sapply(terms, function(t) unname(fe[t] + re[[t]])), .name_repair = "minimal") %>%
    setNames(terms) %>%
    mutate(subj = factor(rownames(re)))
}

safe_fit <- function(...) tryCatch(suppressWarnings(glmmTMB(...)), error = function(e) NULL)

rel <- function(est, truth, est_col, true_col) {
  j <- inner_join(est, truth, by = "subj")
  cor(j[[est_col]], j[[true_col]])^2
}

run_one <- function(n_items, design, rep, n_subj = 300) {
  set.seed(1e5 * n_items + 1e3 * (design == "balanced") + rep)
  d <- simulate_data(n_subj, n_items, design)
  truth <- attr(d, "truth")
  out <- list()

  # A: person effects only, all trials (the 2022 paper's approach)
  f <- safe_fit(y ~ 0 + E + E:g + offset(lp) + (0 + E + E:g | subj), data = d, family = beta_family())
  if (!is.null(f)) { e <- person_estimates(f, c("E", "E:g"))
    out$A <- c(alpha = rel(e, truth, "E", "alpha"), alpha_overall = rel(e, truth, "E", "alpha_overall"),
                 beta = rel(e, truth, "E:g", "beta")) }

  # B: person and item effects, all trials (M3)
  f <- safe_fit(y ~ 0 + E + E:g + offset(lp) + (0 + E + E:g | subj) + (0 + E + E:g | item), data = d, family = beta_family())
  if (!is.null(f)) { e <- person_estimates(f, c("E", "E:g"))
    out$B <- c(alpha = rel(e, truth, "E", "alpha"), alpha_overall = rel(e, truth, "E", "alpha_overall"),
                 beta = rel(e, truth, "E:g", "beta")) }

  # C: two-part model: size of update on changed trials, plus whether people update (M5)
  f <- safe_fit(y ~ 0 + E + E:g + offset(lp) + (0 + E + E:g | subj) + (0 + E + E:g | item),
                data = filter(d, upd == 1), family = beta_family())
  fu <- safe_fit(upd ~ 1 + g + (1 | subj) + (1 | item), data = d, family = binomial())
  if (!is.null(f) && !is.null(fu)) {
    e <- person_estimates(f, c("E", "E:g"))
    eu <- tibble(subj = factor(rownames(ranef(fu)$cond$subj)), u_hat = ranef(fu)$cond$subj[[1]])
    out$C <- c(alpha = rel(e, truth, "E", "alpha"), beta = rel(e, truth, "E:g", "beta"),
               update = rel(eu, truth, "u_hat", "u"))
  }

  imap_dfr(out, ~ tibble(model = .y, parameter = names(.x), reliability = unname(.x))) %>%
    mutate(n_items = n_items, design = design, rep = rep)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  # parallel workers; each large run peaks near 1.7 GB, so keep this modest
  n_cores <- as.integer(Sys.getenv("SIM_CORES", "4"))
  reps <- if (length(args)) as.integer(args[1]) else 20
  designs <- if (length(args) > 1) strsplit(args[2], ",")[[1]] else c("2021", "balanced", "informative")
  out_file <- if (length(args) > 2) args[3] else "local/sims/reliability_results.csv"
  grid <- expand_grid(n_items = c(6, 8, 10, 12, 15, 20, 30, 40), design = designs, rep = seq_len(reps))
  res <- mclapply(seq_len(nrow(grid)), function(k) {
    tryCatch(run_one(grid$n_items[k], grid$design[k], grid$rep[k]), error = function(e) NULL)
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  results <- bind_rows(res)
  dir.create("local/sims", recursive = TRUE, showWarnings = FALSE)
  write_csv(results, out_file)
  results %>%
    group_by(model, parameter, design, n_items) %>%
    summarize(rel = mean(reliability, na.rm = TRUE), lo = quantile(reliability, .1, na.rm = TRUE),
              hi = quantile(reliability, .9, na.rm = TRUE), n_ok = sum(!is.na(reliability)), .groups = "drop") %>%
    print(n = Inf)
}
