# How precisely can a cross-task correlation of updating parameters be estimated?
#
# Each simulated participant completes two tasks generated from the 2021 two-part model
# (see reliability_sim.R), each with its own items. Their evidence weight (alpha) and
# good/bad-news bias (beta) correlate rho across the two tasks.
#
# Two estimators of rho, both using the size-of-update part (trials where the estimate changed):
#   two-stage:    fit each task separately, correlate the person-level estimates (attenuated)
#   single-stage: fit both tasks jointly with correlated person effects, read rho off the
#                 estimated covariance matrix (the latent correlation)

source("sims/reliability_sim.R")

simulate_two_tasks <- function(n_subj, n_items, rho) {
  # person effects: alpha1, beta1, alpha2, beta2 (within-task correlation from M5)
  sd4 <- rep(pars$person_sd[2:3], 2)
  r_ab <- pars$person_cor[2, 3]
  cor4 <- matrix(c(
    1,          r_ab,       rho,        r_ab * rho,
    r_ab,       1,          r_ab * rho, rho,
    rho,        r_ab * rho, 1,          r_ab,
    r_ab * rho, rho,        r_ab,       1), 4)
  ab <- mvrnorm(n_subj, rep(0, 4), covmat(sd4, cor4))
  u <- mvrnorm(n_subj, c(0, 0), covmat(rep(pars$person_sd[1], 2), matrix(c(1, .6, .6, 1), 2)))

  pool <- split(seq_len(nrow(trials_2021)), trials_2021$item_src)
  map_dfr(1:2, function(tk) {
    items <- mvrnorm(n_items, rep(0, 3), covmat(pars$item_sd, pars$item_cor))
    src <- sample(unique(trials_2021$item_src), n_items, replace = n_items > 20)
    d <- expand_grid(subj = seq_len(n_subj), item = seq_len(n_items)) %>% mutate(item_src = src[item])
    rows <- vapply(as.character(d$item_src), function(s) pool[[s]][sample.int(length(pool[[s]]), 1)], integer(1))
    d <- bind_cols(d, trials_2021[rows, c("lp", "E", "g")])
    # index person effects by row here, outside mutate(): inside it, `task` would be a whole column
    a_i <- ab[d$subj, 2 * tk - 1]; b_i <- ab[d$subj, 2 * tk]; u_i <- u[d$subj, tk]
    d %>% mutate(
      task = tk,
      upd = rbinom(n(), 1, plogis(pars$upd0 + pars$upd_g * g + u_i + items[item, 1])),
      mu = plogis(lp + (pars$alpha + a_i + items[item, 2] + (pars$beta + b_i + items[item, 3]) * g) * E),
      y = if_else(upd == 1, clamp(rbeta(n(), mu * pars$phi, (1 - mu) * pars$phi)), clamp(plogis(lp))),
      item = paste(tk, item, sep = "_")
    )
  }) %>%
    mutate(subj = factor(subj), item = factor(item), E1 = E * (task == 1), E2 = E * (task == 2),
         E1g = E1 * g, E2g = E2 * g)
}

# the joint fit only needs its covariance estimates, so skip glmmTMB's standard-error step
# (the per-task fits keep it: ranef() needs it)
fit_light <- function(...) tryCatch(suppressWarnings(glmmTMB(..., se = FALSE)), error = function(e) NULL)

run_latent <- function(n_subj, n_items, rep, rho = .4) {
  set.seed(7e5 + 1e4 * n_items + 10 * n_subj + rep)
  d <- simulate_two_tasks(n_subj, n_items, rho) %>% filter(upd == 1)
  out <- tibble(n_subj = n_subj, n_items = n_items, rep = rep, rho = rho)

  # two-stage
  est <- map(1:2, function(t) {
    f <- safe_fit(y ~ 0 + E + E:g + offset(lp) + (0 + E + E:g | subj) + (0 + E + E:g | item),
                  data = filter(d, task == t), family = beta_family())
    if (is.null(f)) NULL else person_estimates(f, c("E", "E:g"))
  })
  if (!any(map_lgl(est, is.null))) {
    j <- inner_join(est[[1]], est[[2]], by = "subj")
    out$alpha_2stage <- cor(j$E.x, j$E.y)
    out$beta_2stage <- cor(j$`E:g.x`, j$`E:g.y`)
  }

  # single-stage: joint model with a 4 x 4 person covariance matrix
  f <- fit_light(y ~ 0 + E1 + E1g + E2 + E2g + offset(lp) + (0 + E1 + E1g + E2 + E2g | subj) + (0 + E + E:g | item),
                data = d, family = beta_family())
  if (!is.null(f)) {
    cm <- attr(VarCorr(f)$cond$subj, "correlation")
    out$alpha_1stage <- cm["E1", "E2"]
    out$beta_1stage <- cm["E1g", "E2g"]
    out$converged <- f$fit$convergence == 0
  }
  out
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  # parallel workers; each large run peaks near 1.7 GB, so keep this modest
  n_cores <- as.integer(Sys.getenv("SIM_CORES", "4"))
  reps <- if (length(args)) as.integer(args[1]) else 30
  grid <- expand_grid(n_subj = c(300, 600), n_items = c(10, 20, 30), rep = seq_len(reps))
  # each run is saved as it finishes, so an interrupted grid resumes where it stopped
  run_dir <- "local/sims/latent_runs"
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_file <- function(k) file.path(run_dir, sprintf("n%d_k%d_r%02d.rds", grid$n_subj[k], grid$n_items[k], grid$rep[k]))
  todo <- which(!file.exists(map_chr(seq_len(nrow(grid)), run_file)))
  message(length(todo), " of ", nrow(grid), " runs to do")
  invisible(mclapply(todo, function(k) {
    r <- tryCatch(run_latent(grid$n_subj[k], grid$n_items[k], grid$rep[k]), error = function(e) NULL)
    if (!is.null(r)) saveRDS(r, run_file(k))
  }, mc.cores = n_cores, mc.preschedule = FALSE))
  results <- map_dfr(list.files(run_dir, full.names = TRUE), readRDS)
  write_csv(results, "local/sims/latent_cor_results.csv")
}
