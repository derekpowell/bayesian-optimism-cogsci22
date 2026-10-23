# The proposed Stage 1 design with an urn block in place of the 4 neutral actuarial items:
# 6 anchors + 10 new valenced actuarial items (balanced news) + 12 or 24 urn trials.
#
# Urn trials: a stated prior (.3, .5 or .7), an urn of known composition (60/40 or 70/30) and
# 3 or 5 draws, so the evidence (log likelihood ratio) is known exactly. Responses follow
#   logit(post) = logit(prior) + alpha_urn * logLR
# with alpha_urn = .55 + person effect (SD .25), correlated .5 with actuarial (neutral) updating
# and unrelated to the good/bad-news shifts. 90% of urn trials are updated.
# The urn values are assumptions; there are no urn data from the 2021 study.
#
# Without neutral actuarial items, actuarial updating is identified for bad news (alpha_bad) and
# the good - bad difference (beta). The urn block adds alpha_urn and lets a joint model estimate
# the cross-paradigm correlation of conservatism (H1b), here between alpha_bad and alpha_urn.

source("sims/design_mix_sim.R")

urn_pars <- list(alpha = .55, sd = .25, cor_with_alpha = .5, p_update = .9)

simulate_urn <- function(n_subj, n_urn, alpha_act_dev) {
  # alpha_urn person effects correlated with the actuarial alpha deviations
  z <- rnorm(n_subj)
  r <- urn_pars$cor_with_alpha
  dev <- urn_pars$sd * (r * alpha_act_dev / sd(alpha_act_dev) + sqrt(1 - r^2) * z)
  d <- expand_grid(subj = seq_len(n_subj), trial = seq_len(n_urn)) %>%
    mutate(
      prior = sample(c(.3, .5, .7), n(), replace = TRUE),
      p = sample(c(.6, .7), n(), replace = TRUE),
      n_draws = sample(c(3, 5), n(), replace = TRUE),
      red_urn = rbinom(n(), 1, prior),
      reds = rbinom(n(), n_draws, if_else(red_urn == 1, p, 1 - p)),
      Lu = (2 * reds - n_draws) * qlogis(p),
      lp = qlogis(prior),
      w = urn_pars$alpha + dev[subj],
      upd = rbinom(n(), 1, urn_pars$p_update),
      mu = plogis(lp + w * Lu),
      y = if_else(upd == 1, clamp(rbeta(n(), mu * pars$phi, (1 - mu) * pars$phi)), clamp(plogis(lp)))
    )
  structure(d, alpha_urn = urn_pars$alpha + dev)
}

run_urn <- function(n_urn, rep, n_subj = 300) {
  set.seed(1.1e6 + 1e3 * n_urn + rep)
  # actuarial part: D2 without its neutral items
  act <- simulate_mix(n_subj, "D2") %>% filter(type == "valenced")
  truth <- attr(act, "truth")
  urn <- simulate_urn(n_subj, n_urn, truth$alpha - pars$alpha)
  truth$alpha_urn <- attr(urn, "alpha_urn")

  stacked <- bind_rows(
    act %>% filter(upd == 1) %>% transmute(subj, item = as.character(item), y, lp, Ea = E, Eag = E * good, Lu = 0),
    urn %>% filter(upd == 1) %>% transmute(subj = factor(subj), item = "urn", y, lp, Ea = 0, Eag = 0, Lu)
  ) %>% mutate(item = factor(item), act = as.numeric(item != "urn"))

  # item effects apply to actuarial items only (urn trials have exactly known evidence)
  f <- safe_fit(y ~ 0 + Ea + Eag + Lu + offset(lp) + (0 + Ea + Eag + Lu | subj) + (0 + Ea + Eag | item),
                data = stacked, family = beta_family())
  if (is.null(f)) return(NULL)
  e <- person_estimates(f, c("Ea", "Eag", "Lu"))
  cm <- attr(VarCorr(f)$cond$subj, "correlation")
  tibble(
    n_urn = n_urn, rep = rep,
    parameter = c("alpha_bad", "beta", "alpha_urn", "latent_cor_alpha", "true_cor_alpha"),
    value = c(rel_mix(e, truth, "Ea", "alpha_bad"), rel_mix(e, truth, "Eag", "beta"),
              rel_mix(e, truth, "Lu", "alpha_urn"), cm["Ea", "Lu"],
              cor(truth$alpha_bad, truth$alpha_urn))
  )
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_cores <- as.integer(Sys.getenv("SIM_CORES", "4"))
  reps <- if (length(args)) as.integer(args[1]) else 20
  grid <- expand_grid(n_urn = c(12, 24), rep = seq_len(reps))
  run_dir <- "local/sims/urn_runs"
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_file <- function(k) file.path(run_dir, sprintf("u%d_r%02d.rds", grid$n_urn[k], grid$rep[k]))
  todo <- which(!file.exists(map_chr(seq_len(nrow(grid)), run_file)))
  message(length(todo), " of ", nrow(grid), " runs to do")
  invisible(mclapply(todo, function(k) {
    r <- tryCatch(run_urn(grid$n_urn[k], grid$rep[k]), error = function(e) NULL)
    if (!is.null(r)) saveRDS(r, run_file(k))
  }, mc.cores = n_cores, mc.preschedule = FALSE))
  results <- map_dfr(list.files(run_dir, full.names = TRUE), readRDS)
  write_csv(results, "local/sims/urn_block_results.csv")
}
