# Reliability for a proposed Stage 1 actuarial design: 6 anchor items from the 2021 study,
# 10 new valenced items and 4 neutral items (20 per participant), against alternatives.
#
# Valence has three levels here, so updating on good and bad news can each be compared with
# neutral updating:
#   logit(post) = logit(pre) + (alpha + dg * good + db * bad) * E
# 2021 identifies only beta = dg - db (about -.08), so the split is assumed symmetric
# (dg = -.04, db = +.04) and person-level dg and db are independent with SD .15, which
# reproduces the 2021 person SD of beta (about .21).
#
# Designs (all 20 items):
#   D0  20 valenced items, news as in 2021 (baseline)
#   D1  6 anchors + 10 new valenced (news as in 2021) + 4 neutral
#   D2  6 anchors + 10 new valenced (balanced news) + 4 neutral
#   D3  6 anchors + 14 new valenced (balanced news), no neutral items

source("sims/reliability_sim.R")

anchors <- c(2, 4, 5, 9, 12, 20)  # die before 80, die before 70, stroke, cataracts, live past 90, no cold
others <- setdiff(unique(trials_2021$item_src), anchors)

mix_pars <- list(
  dg = -.04, db = .04,
  # person effects: update propensity, alpha, dg, db
  person_sd = c(pars$person_sd[1], pars$person_sd[2], .15, .15),
  person_cor = matrix(c(
    1,   .18, .1,  -.1,
    .18, 1,   -.5, .5,
    .1,  -.5, 1,   0,
    -.1, .5,  0,   1), 4),
  # item effects: update propensity, alpha, dg, db
  item_sd = c(pars$item_sd[1], pars$item_sd[2], .073, .073)
)

design_spec <- function(design) {
  switch(design,
    D0 = tibble(src = sample(unique(trials_2021$item_src), 20), type = "valenced", news = "2021"),
    D1 = bind_rows(tibble(src = anchors, type = "valenced", news = "2021"),
                   tibble(src = sample(others, 10), type = "valenced", news = "2021"),
                   tibble(src = sample(others, 4), type = "neutral", news = "none")),
    D2 = bind_rows(tibble(src = anchors, type = "valenced", news = "2021"),
                   tibble(src = sample(others, 10), type = "valenced", news = "balanced"),
                   tibble(src = sample(others, 4), type = "neutral", news = "none")),
    D3 = bind_rows(tibble(src = anchors, type = "valenced", news = "2021"),
                   tibble(src = sample(others, 14), type = "valenced", news = "balanced"))
  ) %>% mutate(item = row_number())
}

simulate_mix <- function(n_subj, design) {
  spec <- design_spec(design)
  people <- mvrnorm(n_subj, rep(0, 4), covmat(mix_pars$person_sd, mix_pars$person_cor))
  items <- sapply(mix_pars$item_sd, function(s) rnorm(nrow(spec), 0, s))

  pool <- split(seq_len(nrow(trials_2021)), trials_2021$item_src)
  d <- expand_grid(subj = seq_len(n_subj), item = spec$item) %>% left_join(spec, by = "item")
  rows <- vapply(as.character(d$src), function(s) pool[[s]][sample.int(length(pool[[s]]), 1)], integer(1))
  d <- bind_cols(d, trials_2021[rows, c("lp", "E", "g")])
  d$g[d$news == "balanced"] <- rbinom(sum(d$news == "balanced"), 1, .5)

  # person and item effects indexed outside mutate()
  pe <- people[d$subj, ]; ie <- items[d$item, ]
  d %>% mutate(
    good = if_else(type == "neutral", 0, g),
    bad = if_else(type == "neutral", 0, 1 - g),
    upd = rbinom(n(), 1, plogis(pars$upd0 + u_shift(good) + pe[, 1] + ie[, 1])),
    w = pars$alpha + pe[, 2] + ie[, 2] + (mix_pars$dg + pe[, 3] + ie[, 3]) * good + (mix_pars$db + pe[, 4] + ie[, 4]) * bad,
    mu = plogis(lp + w * E),
    y = if_else(upd == 1, clamp(rbeta(n(), mu * pars$phi, (1 - mu) * pars$phi)), clamp(plogis(lp))),
    subj = factor(subj), item = factor(item),
    Eg = E * good, Eb = E * bad
  ) %>%
    structure(truth = tibble(
      subj = factor(seq_len(n_subj)),
      alpha = pars$alpha + people[, 2],                       # neutral updating
      dg = mix_pars$dg + people[, 3], db = mix_pars$db + people[, 4],
      beta = (mix_pars$dg + people[, 3]) - (mix_pars$db + people[, 4]),
      alpha_bad = pars$alpha + people[, 2] + mix_pars$db + people[, 4]))
}
u_shift <- function(good) pars$upd_g * good

# like rel() in reliability_sim.R, but safe when estimate and truth share a column name, and
# scoring a collapsed (constant) estimate as 0: the model then carries no person-level information
rel_mix <- function(est, truth, est_col, true_col) {
  j <- inner_join(est %>% select(subj, x = all_of(est_col)), truth %>% select(subj, y = all_of(true_col)), by = "subj")
  if (sd(j$x) < 1e-6) return(0)
  cor(j$x, j$y)^2
}

run_mix <- function(design, rep, n_subj = 300) {
  set.seed(9e5 + 1e3 * match(design, c("D0", "D1", "D2", "D3")) + rep)
  d <- simulate_mix(n_subj, design)
  truth <- attr(d, "truth")
  dc <- filter(d, upd == 1)
  has_neutral <- any(d$good == 0 & d$bad == 0)

  if (has_neutral) {
    f <- safe_fit(y ~ 0 + E + Eg + Eb + offset(lp) + (0 + E + Eg + Eb | subj) + (0 + E + Eg + Eb | item),
                  data = dc, family = beta_family())
    if (is.null(f)) return(NULL)
    e <- person_estimates(f, c("E", "Eg", "Eb")) %>% mutate(beta = Eg - Eb)
    out <- c(alpha_neutral = rel_mix(e, truth, "E", "alpha"), good_shift = rel_mix(e, truth, "Eg", "dg"),
             bad_shift = rel_mix(e, truth, "Eb", "db"), beta = rel_mix(e, truth, "beta", "beta"))
  } else {
    # without neutral items only bad-news updating and the good-bad difference are identified
    f <- safe_fit(y ~ 0 + E + Eg + offset(lp) + (0 + E + Eg | subj) + (0 + E + Eg | item),
                  data = mutate(dc, Eg = E * good), family = beta_family())
    if (is.null(f)) return(NULL)
    e <- person_estimates(f, c("E", "Eg"))
    out <- c(alpha_bad = rel_mix(e, truth, "E", "alpha_bad"), beta = rel_mix(e, truth, "Eg", "beta"))
  }
  tibble(design = design, rep = rep, parameter = names(out), reliability = unname(out))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_cores <- as.integer(Sys.getenv("SIM_CORES", "4"))
  reps <- if (length(args)) as.integer(args[1]) else 20
  grid <- expand_grid(design = c("D0", "D1", "D2", "D3"), rep = seq_len(reps))
  run_dir <- "local/sims/mix_runs"
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_file <- function(k) file.path(run_dir, sprintf("%s_r%02d.rds", grid$design[k], grid$rep[k]))
  todo <- which(!file.exists(map_chr(seq_len(nrow(grid)), run_file)))
  message(length(todo), " of ", nrow(grid), " runs to do")
  invisible(mclapply(todo, function(k) {
    r <- tryCatch(run_mix(grid$design[k], grid$rep[k]), error = function(e) NULL)
    if (!is.null(r)) saveRDS(r, run_file(k))
  }, mc.cores = n_cores, mc.preschedule = FALSE))
  results <- map_dfr(list.files(run_dir, full.names = TRUE), readRDS)
  write_csv(results, "local/sims/design_mix_results.csv")
}
