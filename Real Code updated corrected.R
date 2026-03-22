# ============================================================
# BE-COD REAL DATA ANALYSIS: CORRECTED AND STABLE VERSION
# ToothGrowth illustration for the paper
# ============================================================

rm(list = ls())
gc()

# ----------------------------
# 0. Packages
# ----------------------------
needed_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "forcats", "tibble", "scales", "broom", "patchwork"
)

to_install <- needed_pkgs[!vapply(needed_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

invisible(lapply(needed_pkgs, library, character.only = TRUE))

# ----------------------------
# 1. Output folders
# ----------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
base_dir  <- file.path(getwd(), paste0("BE_COD_ToothGrowth_", timestamp))
fig_dir   <- file.path(base_dir, "figures")
tab_dir   <- file.path(base_dir, "tables")
obj_dir   <- file.path(base_dir, "objects")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir,  recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("Output directory:\n", base_dir, "\n")
cat("============================================================\n\n")

# ----------------------------
# 2. Helper functions
# ----------------------------
save_plot_both <- function(plot_obj, filename, width = 8, height = 5, dpi = 320) {
  ggplot2::ggsave(
    filename = file.path(fig_dir, paste0(filename, ".png")),
    plot = plot_obj, width = width, height = height, dpi = dpi
  )
  ggplot2::ggsave(
    filename = file.path(fig_dir, paste0(filename, ".pdf")),
    plot = plot_obj, width = width, height = height
  )
}

save_table_csv <- function(df, filename) {
  readr::write_csv(df, file.path(tab_dir, paste0(filename, ".csv")))
}

safe_solve <- function(A, ridge = 1e-8) {
  A <- as.matrix(A)
  p <- nrow(A)
  solve(A + diag(ridge, p))
}

stable_logdet <- function(A, eps = 1e-8) {
  A <- as.matrix(A)
  A <- 0.5 * (A + t(A))
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  ev[ev < eps] <- eps
  sum(log(ev))
}

row_quad_form <- function(X, V) {
  rowSums((X %*% V) * X)
}

mvnorm_logdens <- function(y, mean, Sigma, eps = 1e-8) {
  Sigma <- as.matrix(Sigma)
  Sigma <- 0.5 * (Sigma + t(Sigma)) + diag(eps, nrow(Sigma))
  k <- length(y)
  cholS <- chol(Sigma)
  z <- backsolve(cholS, y - mean, transpose = TRUE)
  quad <- sum(z^2)
  logdetS <- 2 * sum(log(diag(cholS)))
  -0.5 * (k * log(2 * pi) + logdetS + quad)
}

# ----------------------------
# 3. Load and preprocess ToothGrowth
# ----------------------------
data("ToothGrowth", package = "datasets")

tg <- ToothGrowth %>%
  dplyr::mutate(
    supp   = factor(supp, levels = c("VC", "OJ")),
    dose   = as.numeric(dose),
    dose_f = factor(dose, levels = sort(unique(dose))),
    s      = ifelse(supp == "OJ", 1, 0),
    z      = log(dose),
    z2     = z^2,
    cell   = paste0(supp, "_", dose)
  ) %>%
  dplyr::arrange(supp, dose)

make_X <- function(df) {
  cbind(
    `(Intercept)` = 1,
    s   = df$s,
    z   = df$z,
    z2  = df$z2,
    sz  = df$s * df$z,
    sz2 = df$s * df$z2
  )
}

X_all <- make_X(tg)
y_all <- tg$len

# ----------------------------
# 4. Descriptive summaries
# ----------------------------
cell_summary_all <- tg %>%
  dplyr::group_by(supp, dose) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    mean   = mean(len),
    sd     = sd(len),
    se     = sd(len) / sqrt(dplyr::n()),
    median = median(len),
    min    = min(len),
    max    = max(len),
    .groups = "drop"
  )

mean_profile <- tg %>%
  dplyr::group_by(supp, dose) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(len),
    sd = sd(len),
    se = sd / sqrt(n),
    lower = mean - qt(0.975, n - 1) * se,
    upper = mean + qt(0.975, n - 1) * se,
    .groups = "drop"
  )

heat_df <- tg %>%
  dplyr::group_by(supp, dose_f) %>%
  dplyr::summarise(mean_len = mean(len), .groups = "drop")

cat("\n==================== FULL DATA SUMMARY ====================\n")
print(cell_summary_all)

save_table_csv(cell_summary_all, "table_01_cell_summary_full")
save_table_csv(mean_profile,      "table_02_mean_profile_full")
save_table_csv(heat_df,           "table_03_heatmap_means_full")

# ----------------------------
# 5. Historical / pseudo-current split
#    5 per cell in D0, 5 per cell in D1
# ----------------------------
set.seed(12345)

tg_split <- tg %>%
  dplyr::group_by(supp, dose) %>%
  dplyr::mutate(rand_order = sample(seq_len(dplyr::n()))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sample_role = ifelse(rand_order <= 5, "D0_historical", "D1_current"))

D0 <- tg_split %>%
  dplyr::filter(sample_role == "D0_historical") %>%
  dplyr::arrange(supp, dose)

D1 <- tg_split %>%
  dplyr::filter(sample_role == "D1_current") %>%
  dplyr::arrange(supp, dose)

X0 <- make_X(D0)
y0 <- D0$len
X1 <- make_X(D1)
y1 <- D1$len

cell_summary_D0 <- D0 %>%
  dplyr::group_by(supp, dose) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(len),
    sd = sd(len),
    .groups = "drop"
  )

cell_summary_D1 <- D1 %>%
  dplyr::group_by(supp, dose) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(len),
    sd = sd(len),
    .groups = "drop"
  )

cat("\n==================== D0 SUMMARY ====================\n")
print(cell_summary_D0)

cat("\n==================== D1 SUMMARY ====================\n")
print(cell_summary_D1)

save_table_csv(cell_summary_D0, "table_04_cell_summary_D0")
save_table_csv(cell_summary_D1, "table_05_cell_summary_D1")

# ----------------------------
# 6. Historical model fit
# ----------------------------
fit_D0 <- lm(len ~ s + z + z2 + s:z + s:z2, data = D0)

coef_D0 <- broom::tidy(fit_D0, conf.int = TRUE)
glance_D0 <- broom::glance(fit_D0)
anova_D0 <- broom::tidy(anova(fit_D0))

cat("\n==================== HISTORICAL MODEL COEFFICIENTS ====================\n")
print(coef_D0)

cat("\n==================== HISTORICAL MODEL GLANCE ====================\n")
print(glance_D0)

cat("\n==================== HISTORICAL MODEL ANOVA ====================\n")
print(anova_D0)

save_table_csv(coef_D0,  "table_06_historical_coefficients")
save_table_csv(glance_D0, "table_07_historical_glance")
save_table_csv(anova_D0, "table_08_historical_anova")

beta0_hat <- coef(fit_D0)
V0_hat <- vcov(fit_D0)
sigma2_hat <- summary(fit_D0)$sigma^2
sigma_hat <- sqrt(sigma2_hat)

cat("\nEstimated sigma^2 from D0 =", round(sigma2_hat, 6), "\n")

# ----------------------------
# 7. EB commensurate prior and robust mixture prior
# ----------------------------
p <- ncol(X0)
param_names <- colnames(X0)

m_comm <- as.numeric(beta0_hat)
V_comm <- as.matrix(V0_hat) + diag(1e-8, p)

kappa2 <- 1e4
m_vague <- rep(0, p)
V_vague <- diag(kappa2, p)

omega0 <- 0.70

eb_prior_table <- tibble::tibble(
  parameter = param_names,
  prior_mean_comm = m_comm,
  prior_sd_comm   = sqrt(diag(V_comm)),
  tau_proxy       = 1 / diag(V_comm)
)

cat("\n==================== EB PRIOR SUMMARY ====================\n")
print(eb_prior_table)

save_table_csv(eb_prior_table, "table_09_eb_prior_summary")

# ----------------------------
# 8. Target contrasts
# ----------------------------
make_delta_row <- function(d) {
  z <- log(d)
  c(0, 1, 0, 0, z, z^2)
}

L_targets <- rbind(
  "Delta_0.5" = make_delta_row(0.5),
  "Delta_1.0" = make_delta_row(1.0),
  "Delta_2.0" = make_delta_row(2.0),
  "beta_S1"   = c(0, 0, 0, 0, 1, 0),
  "beta_S2"   = c(0, 0, 0, 0, 0, 1)
)

# ----------------------------
# 9. Posterior functions
# ----------------------------
posterior_normal_linear <- function(X, y, m_prior, V_prior, sigma2) {
  Prec_prior <- safe_solve(V_prior)
  Prec_post  <- Prec_prior + crossprod(X) / sigma2
  V_post     <- safe_solve(Prec_post)
  m_post     <- V_post %*% (Prec_prior %*% m_prior + crossprod(X, y) / sigma2)
  
  list(
    m_post = as.numeric(m_post),
    V_post = as.matrix(V_post),
    Prec_post = as.matrix(Prec_post)
  )
}

marginal_loglik_linear <- function(X, y, m_prior, V_prior, sigma2) {
  mu_y <- as.numeric(X %*% m_prior)
  Sigma_y <- sigma2 * diag(nrow(X)) + X %*% V_prior %*% t(X)
  mvnorm_logdens(y = y, mean = mu_y, Sigma = Sigma_y)
}

posterior_mixture_linear <- function(X, y, sigma2, m1, V1, m2, V2, omega = 0.7) {
  post1 <- posterior_normal_linear(X, y, m1, V1, sigma2)
  post2 <- posterior_normal_linear(X, y, m2, V2, sigma2)
  
  logml1 <- marginal_loglik_linear(X, y, m1, V1, sigma2)
  logml2 <- marginal_loglik_linear(X, y, m2, V2, sigma2)
  
  a1 <- log(omega) + logml1
  a2 <- log(1 - omega) + logml2
  mmax <- max(a1, a2)
  
  w1 <- exp(a1 - mmax) / (exp(a1 - mmax) + exp(a2 - mmax))
  w2 <- 1 - w1
  
  list(
    post1 = post1,
    post2 = post2,
    w1 = w1,
    w2 = w2,
    logml1 = logml1,
    logml2 = logml2
  )
}

contrast_summary_single <- function(m_post, V_post, L, method_name, level = 0.95) {
  alpha <- 1 - level
  zcrit <- qnorm(1 - alpha / 2)
  
  est <- as.numeric(L %*% m_post)
  se  <- sqrt(diag(L %*% V_post %*% t(L)))
  
  tibble::tibble(
    method = method_name,
    estimand = rownames(L),
    mean = est,
    sd = se,
    lower = est - zcrit * se,
    upper = est + zcrit * se
  )
}

contrast_summary_mixture <- function(mix_post, L, method_name = "Robust Mixture", level = 0.95) {
  alpha <- 1 - level
  zcrit <- qnorm(1 - alpha / 2)
  
  mu1 <- as.numeric(L %*% mix_post$post1$m_post)
  mu2 <- as.numeric(L %*% mix_post$post2$m_post)
  
  S1 <- L %*% mix_post$post1$V_post %*% t(L)
  S2 <- L %*% mix_post$post2$V_post %*% t(L)
  
  w1 <- mix_post$w1
  w2 <- mix_post$w2
  
  mu_mix <- w1 * mu1 + w2 * mu2
  
  S_mix <- w1 * (S1 + tcrossprod(mu1 - mu_mix)) +
    w2 * (S2 + tcrossprod(mu2 - mu_mix))
  
  se_mix <- sqrt(diag(S_mix))
  
  tibble::tibble(
    method = method_name,
    estimand = rownames(L),
    mean = mu_mix,
    sd = se_mix,
    lower = mu_mix - zcrit * se_mix,
    upper = mu_mix + zcrit * se_mix
  )
}

# ----------------------------
# 10. Posterior inference on D1
# ----------------------------
post_comm <- posterior_normal_linear(X1, y1, m_comm, V_comm, sigma2_hat)
post_vague <- posterior_normal_linear(X1, y1, m_vague, V_vague, sigma2_hat)
post_mix <- posterior_mixture_linear(
  X = X1, y = y1, sigma2 = sigma2_hat,
  m1 = m_comm, V1 = V_comm,
  m2 = m_vague, V2 = V_vague,
  omega = omega0
)

contr_comm <- contrast_summary_single(
  m_post = post_comm$m_post,
  V_post = post_comm$V_post,
  L = L_targets,
  method_name = "Commensurate EB"
)

contr_vague <- contrast_summary_single(
  m_post = post_vague$m_post,
  V_post = post_vague$V_post,
  L = L_targets,
  method_name = "Vague Bayes"
)

contr_mix <- contrast_summary_mixture(
  mix_post = post_mix,
  L = L_targets,
  method_name = "Robust Mixture"
)

contrast_table <- dplyr::bind_rows(contr_comm, contr_vague, contr_mix) %>%
  dplyr::mutate(
    estimand = factor(estimand, levels = rownames(L_targets)),
    method = factor(method, levels = c("Vague Bayes", "Commensurate EB", "Robust Mixture"))
  ) %>%
  dplyr::arrange(estimand, method)

posterior_weight_table <- tibble::tibble(
  prior_weight_commensurate = omega0,
  posterior_weight_commensurate = post_mix$w1,
  posterior_weight_vague = post_mix$w2,
  log_marglik_commensurate = post_mix$logml1,
  log_marglik_vague = post_mix$logml2
)

cat("\n==================== POSTERIOR MIXTURE WEIGHTS ====================\n")
print(posterior_weight_table)

cat("\n==================== TARGET CONTRAST SUMMARY ====================\n")
print(contrast_table)

save_table_csv(contrast_table,       "table_10_target_contrast_summary")
save_table_csv(posterior_weight_table, "table_11_mixture_weight_summary")

# ----------------------------
# 11. Candidate design space for future exact design
#     For real data, use the actual 2 x 3 factorial grid
# ----------------------------
cand <- expand.grid(
  supp = factor(c("VC", "OJ"), levels = c("VC", "OJ")),
  dose = sort(unique(tg$dose))
) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    s = ifelse(supp == "OJ", 1, 0),
    z = log(dose),
    z2 = z^2,
    dose_f = factor(dose, levels = sort(unique(dose))),
    cell = paste0(supp, "_", dose)
  )

# ----------------------------
# 12. Design utility
#     IMPORTANT:
#     The 5 target estimands are linearly redundant.
#     So we use a regularized logdet:
#       U = -log det( Cov(L beta | xi) + eps I )
#     This keeps the criterion finite and stable.
# ----------------------------
design_from_indices <- function(idx, cand_df) {
  cand_df[idx, , drop = FALSE]
}

counts_from_indices <- function(idx, cand_df) {
  design_from_indices(idx, cand_df) %>%
    dplyr::count(supp, dose, name = "count") %>%
    tidyr::complete(supp, dose, fill = list(count = 0)) %>%
    dplyr::arrange(supp, dose)
}

posterior_cov_for_design <- function(idx, cand_df, V_prior, sigma2) {
  Xd <- make_X(design_from_indices(idx, cand_df))
  Prec_prior <- safe_solve(V_prior)
  Prec_post <- Prec_prior + crossprod(Xd) / sigma2
  safe_solve(Prec_post)
}

targeted_utility_regularized <- function(V_post, L, eps = 1e-8) {
  S <- L %*% V_post %*% t(L)
  S_reg <- 0.5 * (S + t(S)) + diag(eps, nrow(S))
  -stable_logdet(S_reg, eps = eps)
}

utility_comm <- function(idx, cand_df, L, sigma2, V_prior, eps = 1e-8) {
  V_post <- posterior_cov_for_design(idx, cand_df, V_prior, sigma2)
  targeted_utility_regularized(V_post, L, eps = eps)
}

utility_robust <- function(idx, cand_df, L, sigma2, V1, V2, omega = 0.7, eps = 1e-8) {
  V_post1 <- posterior_cov_for_design(idx, cand_df, V1, sigma2)
  V_post2 <- posterior_cov_for_design(idx, cand_df, V2, sigma2)
  
  omega * targeted_utility_regularized(V_post1, L, eps = eps) +
    (1 - omega) * targeted_utility_regularized(V_post2, L, eps = eps)
}

# ----------------------------
# 13. Coordinate exchange
# ----------------------------
coordinate_exchange <- function(
    n_runs,
    cand_df,
    utility_fun,
    init_idx = NULL,
    max_iter = 100,
    tol = 1e-10,
    verbose = TRUE
) {
  n_cand <- nrow(cand_df)
  
  if (is.null(init_idx)) {
    init_idx <- sample(seq_len(n_cand), size = n_runs, replace = TRUE)
  }
  
  idx <- init_idx
  current_utility <- utility_fun(idx)
  trace <- current_utility
  
  if (verbose) {
    cat("\nInitial utility =", round(current_utility, 6), "\n")
  }
  
  for (iter in seq_len(max_iter)) {
    improved <- FALSE
    
    for (i in seq_len(n_runs)) {
      util_vec <- vapply(seq_len(n_cand), function(j) {
        idx_try <- idx
        idx_try[i] <- j
        val <- utility_fun(idx_try)
        if (!is.finite(val)) val <- -1e12
        val
      }, numeric(1))
      
      j_best <- which.max(util_vec)
      best_val <- util_vec[j_best]
      
      if (best_val > current_utility + tol) {
        idx[i] <- j_best
        current_utility <- best_val
        improved <- TRUE
      }
    }
    
    trace <- c(trace, current_utility)
    
    if (verbose) {
      cat("Iteration", iter, "utility =", round(current_utility, 6), "\n")
    }
    
    if (!improved) break
  }
  
  list(
    idx = idx,
    utility = current_utility,
    trace = trace,
    counts = counts_from_indices(idx, cand_df)
  )
}

multi_start_coordinate_exchange <- function(
    n_starts,
    n_runs,
    cand_df,
    utility_fun,
    balanced_idx = NULL,
    max_iter = 60,
    verbose = FALSE
) {
  res_list <- vector("list", n_starts)
  
  for (s in seq_len(n_starts)) {
    if (!is.null(balanced_idx) && s == 1) {
      init_idx <- sample(balanced_idx)
    } else {
      init_idx <- sample(seq_len(nrow(cand_df)), size = n_runs, replace = TRUE)
    }
    
    res_list[[s]] <- coordinate_exchange(
      n_runs = n_runs,
      cand_df = cand_df,
      utility_fun = utility_fun,
      init_idx = init_idx,
      max_iter = max_iter,
      verbose = verbose
    )
  }
  
  best_id <- which.max(vapply(res_list, function(x) x$utility, numeric(1)))
  res_list[[best_id]]
}

# ----------------------------
# 14. Future design comparison
# ----------------------------
n_future <- 30

if (n_future %% nrow(cand) != 0) {
  stop("n_future must be divisible by the number of candidate cells (= 6 here).")
}

balanced_idx <- rep(seq_len(nrow(cand)), each = n_future / nrow(cand))

u_becod <- function(idx) {
  utility_robust(
    idx = idx, cand_df = cand, L = L_targets, sigma2 = sigma2_hat,
    V1 = V_comm, V2 = V_vague, omega = omega0, eps = 1e-8
  )
}

u_nibod <- function(idx) {
  utility_comm(
    idx = idx, cand_df = cand, L = L_targets, sigma2 = sigma2_hat,
    V_prior = V_vague, eps = 1e-8
  )
}

set.seed(2026)
becod_design <- multi_start_coordinate_exchange(
  n_starts = 10,
  n_runs = n_future,
  cand_df = cand,
  utility_fun = u_becod,
  balanced_idx = balanced_idx,
  max_iter = 60,
  verbose = FALSE
)

set.seed(2027)
nibod_design <- multi_start_coordinate_exchange(
  n_starts = 10,
  n_runs = n_future,
  cand_df = cand,
  utility_fun = u_nibod,
  balanced_idx = balanced_idx,
  max_iter = 60,
  verbose = FALSE
)

bal_counts <- counts_from_indices(balanced_idx, cand)
bal_utility_becod <- u_becod(balanced_idx)
bal_utility_nibod <- u_nibod(balanced_idx)

design_comparison <- dplyr::bind_rows(
  becod_design$counts %>% dplyr::mutate(method = "BE-COD"),
  nibod_design$counts %>% dplyr::mutate(method = "NI-BOD"),
  bal_counts %>% dplyr::mutate(method = "Balanced")
) %>%
  dplyr::select(method, supp, dose, count) %>%
  dplyr::mutate(
    method = factor(method, levels = c("Balanced", "NI-BOD", "BE-COD"))
  ) %>%
  dplyr::arrange(method, supp, dose)

utility_comparison <- tibble::tibble(
  method = c("BE-COD", "NI-BOD", "Balanced"),
  utility_under_BECOD = c(
    becod_design$utility,
    u_becod(nibod_design$idx),
    bal_utility_becod
  ),
  utility_under_NIBOD = c(
    u_nibod(becod_design$idx),
    nibod_design$utility,
    bal_utility_nibod
  )
)

design_counts_wide <- design_comparison %>%
  tidyr::pivot_wider(names_from = method, values_from = count)

cat("\n==================== FUTURE EXACT DESIGN ALLOCATIONS ====================\n")
print(design_comparison)

cat("\n==================== DESIGN UTILITY COMPARISON ====================\n")
print(utility_comparison)

save_table_csv(design_comparison,  "table_12_future_design_allocations")
save_table_csv(design_counts_wide, "table_13_future_design_allocations_wide")
save_table_csv(utility_comparison, "table_14_design_utility_comparison")

# ----------------------------
# 15. Posterior target SDs under future designs
# ----------------------------
future_target_sd <- function(idx, cand_df, V_prior, sigma2, L, method_name) {
  V_post <- posterior_cov_for_design(idx, cand_df, V_prior, sigma2)
  S <- L %*% V_post %*% t(L)
  tibble::tibble(
    method = method_name,
    estimand = rownames(L),
    posterior_sd = sqrt(diag(S))
  )
}

future_sd_table <- dplyr::bind_rows(
  future_target_sd(balanced_idx, cand, V_vague, sigma2_hat, L_targets, "Balanced"),
  future_target_sd(nibod_design$idx, cand, V_vague, sigma2_hat, L_targets, "NI-BOD"),
  future_target_sd(becod_design$idx, cand, V_comm,  sigma2_hat, L_targets, "BE-COD (commensurate core)")
)

cat("\n==================== FUTURE TARGET POSTERIOR SD TABLE ====================\n")
print(future_sd_table)

save_table_csv(future_sd_table, "table_15_future_target_posterior_sd")

# ----------------------------
# 16. Fitted curves on D1
# ----------------------------
dose_grid <- seq(min(tg$dose), max(tg$dose), length.out = 200)

pred_grid <- expand.grid(
  supp = factor(c("VC", "OJ"), levels = c("VC", "OJ")),
  dose = dose_grid
) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    s = ifelse(supp == "OJ", 1, 0),
    z = log(dose),
    z2 = z^2
  )

Xg <- make_X(pred_grid)

curve_summary_single <- function(Xnew, m_post, V_post, method_name) {
  mu <- as.numeric(Xnew %*% m_post)
  v  <- row_quad_form(Xnew, V_post)
  tibble::tibble(
    method = method_name,
    mean = mu,
    sd = sqrt(v),
    lower = mu - 1.96 * sqrt(v),
    upper = mu + 1.96 * sqrt(v)
  )
}

curve_summary_mixture <- function(Xnew, mix_post, method_name = "Robust Mixture") {
  mu1 <- as.numeric(Xnew %*% mix_post$post1$m_post)
  mu2 <- as.numeric(Xnew %*% mix_post$post2$m_post)
  
  v1 <- row_quad_form(Xnew, mix_post$post1$V_post)
  v2 <- row_quad_form(Xnew, mix_post$post2$V_post)
  
  w1 <- mix_post$w1
  w2 <- mix_post$w2
  
  mu_mix <- w1 * mu1 + w2 * mu2
  var_mix <- w1 * (v1 + (mu1 - mu_mix)^2) +
    w2 * (v2 + (mu2 - mu_mix)^2)
  
  tibble::tibble(
    method = method_name,
    mean = mu_mix,
    sd = sqrt(var_mix),
    lower = mu_mix - 1.96 * sqrt(var_mix),
    upper = mu_mix + 1.96 * sqrt(var_mix)
  )
}

curve_comm <- curve_summary_single(Xg, post_comm$m_post, post_comm$V_post, "Commensurate EB")
curve_vague <- curve_summary_single(Xg, post_vague$m_post, post_vague$V_post, "Vague Bayes")
curve_mix <- curve_summary_mixture(Xg, post_mix, "Robust Mixture")

curve_df <- dplyr::bind_cols(pred_grid, curve_comm %>% dplyr::select(-method)) %>%
  dplyr::mutate(method = "Commensurate EB") %>%
  dplyr::bind_rows(
    dplyr::bind_cols(pred_grid, curve_vague %>% dplyr::select(-method)) %>%
      dplyr::mutate(method = "Vague Bayes"),
    dplyr::bind_cols(pred_grid, curve_mix %>% dplyr::select(-method)) %>%
      dplyr::mutate(method = "Robust Mixture")
  ) %>%
  dplyr::mutate(
    method = factor(method, levels = c("Vague Bayes", "Commensurate EB", "Robust Mixture"))
  )

save_table_csv(curve_df, "table_16_fitted_curve_grid")

# ----------------------------
# 17. Plot theme
# ----------------------------
theme_be <- ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  )

# ----------------------------
# 18. Plots
# ----------------------------

# Plot 1: Raw data boxplot + jitter
p1 <- ggplot2::ggplot(tg, ggplot2::aes(x = dose_f, y = len, fill = supp)) +
  ggplot2::geom_boxplot(alpha = 0.45, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.8)) +
  ggplot2::geom_point(
    ggplot2::aes(color = supp),
    position = ggplot2::position_jitterdodge(jitter.width = 0.10, dodge.width = 0.8),
    alpha = 0.78, size = 2
  ) +
  ggplot2::labs(
    title = "ToothGrowth: Raw odontoblast length by dose and supplement",
    x = "Dose (mg/day)",
    y = "Tooth length"
  ) +
  theme_be

print(p1)
save_plot_both(p1, "fig_01_raw_box_jitter", width = 9, height = 5.5)

# Plot 2: Observed mean profiles
p2 <- ggplot2::ggplot(mean_profile, ggplot2::aes(x = dose, y = mean, color = supp, group = supp)) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.08, linewidth = 0.8) +
  ggplot2::scale_x_continuous(breaks = sort(unique(tg$dose))) +
  ggplot2::labs(
    title = "Observed mean dose-response profiles with 95% confidence intervals",
    x = "Dose (mg/day)",
    y = "Mean tooth length"
  ) +
  theme_be

print(p2)
save_plot_both(p2, "fig_02_mean_profile_ci", width = 8.5, height = 5.2)

# Plot 3: D0 / D1 split
p3 <- ggplot2::ggplot(tg_split, ggplot2::aes(x = dose_f, y = len, color = sample_role, shape = sample_role)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0), size = 2.4, alpha = 0.85) +
  ggplot2::facet_wrap(~ supp, nrow = 1) +
  ggplot2::labs(
    title = "Historical (D0) and pseudo-current (D1) split",
    x = "Dose (mg/day)",
    y = "Tooth length",
    color = "Sample role",
    shape = "Sample role"
  ) +
  theme_be

print(p3)
save_plot_both(p3, "fig_03_split_D0_D1", width = 9.5, height = 4.8)

# Plot 4: Historical coefficient forest
coef_plot_df <- coef_D0 %>%
  dplyr::mutate(term = factor(term, levels = rev(term)))

p4 <- ggplot2::ggplot(coef_plot_df, ggplot2::aes(x = estimate, y = term)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = conf.low, xmax = conf.high),
    orientation = "y",
    width = 0.18,
    linewidth = 0.9
  ) +
  ggplot2::labs(
    title = "Historical model coefficients (D0) with 95% confidence intervals",
    x = "Estimate",
    y = NULL
  ) +
  theme_be

print(p4)
save_plot_both(p4, "fig_04_historical_coef_forest", width = 8, height = 4.8)

# Plot 5: Posterior target intervals
p5 <- ggplot2::ggplot(contrast_table, ggplot2::aes(x = mean, y = estimand, color = method)) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
  ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.55), size = 2.6) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower, xmax = upper),
    orientation = "y",
    position = ggplot2::position_dodge(width = 0.55),
    width = 0.18,
    linewidth = 0.9
  ) +
  ggplot2::labs(
    title = "Posterior summaries for target estimands on D1",
    x = "Posterior mean with 95% interval",
    y = NULL,
    color = "Method"
  ) +
  theme_be

print(p5)
save_plot_both(p5, "fig_05_target_contrast_posterior", width = 9.5, height = 5.5)

# Plot 6: Fitted curves by method
p6 <- ggplot2::ggplot() +
  ggplot2::geom_point(data = D1, ggplot2::aes(x = dose, y = len, color = supp), alpha = 0.70, size = 2) +
  ggplot2::geom_ribbon(
    data = curve_df,
    ggplot2::aes(x = dose, ymin = lower, ymax = upper, fill = supp),
    alpha = 0.18
  ) +
  ggplot2::geom_line(
    data = curve_df,
    ggplot2::aes(x = dose, y = mean, color = supp),
    linewidth = 1.15
  ) +
  ggplot2::facet_wrap(~ method, nrow = 1) +
  ggplot2::labs(
    title = "Posterior fitted mean curves on D1",
    subtitle = "Points = D1 observations; ribbons = 95% intervals for mean curves",
    x = "Dose (mg/day)",
    y = "Tooth length"
  ) +
  theme_be

print(p6)
save_plot_both(p6, "fig_06_posterior_fitted_curves", width = 13, height = 4.8)

# Plot 7: Heatmap of observed cell means
p7 <- ggplot2::ggplot(heat_df, ggplot2::aes(x = dose_f, y = supp, fill = mean_len)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean_len)), size = 4.2, fontface = "bold") +
  ggplot2::labs(
    title = "Observed mean tooth length by factorial cell",
    x = "Dose (mg/day)",
    y = "Supplement",
    fill = "Mean length"
  ) +
  theme_be

print(p7)
save_plot_both(p7, "fig_07_cell_mean_heatmap", width = 7.5, height = 4.5)

# Plot 8: Future design allocation heatmap
alloc_plot_df <- design_comparison %>%
  dplyr::mutate(
    dose_f = factor(dose, levels = sort(unique(dose))),
    method = factor(method, levels = c("Balanced", "NI-BOD", "BE-COD"))
  )

p8 <- ggplot2::ggplot(alloc_plot_df, ggplot2::aes(x = dose_f, y = supp, fill = count)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = count), size = 4.8, fontface = "bold") +
  ggplot2::facet_wrap(~ method, nrow = 1) +
  ggplot2::labs(
    title = paste0("Recommended future exact design allocations (n = ", n_future, ")"),
    x = "Dose (mg/day)",
    y = "Supplement",
    fill = "Count"
  ) +
  theme_be

print(p8)
save_plot_both(p8, "fig_08_future_design_heatmap", width = 12, height = 4.6)

# Plot 9: Utility comparison
utility_plot_df <- utility_comparison %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("utility_"),
    names_to = "criterion",
    values_to = "utility"
  ) %>%
  dplyr::mutate(
    criterion = dplyr::recode(
      criterion,
      utility_under_BECOD = "BE-COD criterion",
      utility_under_NIBOD = "NI-BOD criterion"
    )
  )

p9 <- ggplot2::ggplot(utility_plot_df, ggplot2::aes(x = method, y = utility, fill = criterion)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65, alpha = 0.9) +
  ggplot2::labs(
    title = "Utility comparison across future exact designs",
    x = NULL,
    y = "Utility value",
    fill = "Criterion"
  ) +
  theme_be

print(p9)
save_plot_both(p9, "fig_09_design_utility_barplot", width = 8.5, height = 5)

# Plot 10: Coordinate exchange traces
trace_df <- dplyr::bind_rows(
  tibble::tibble(iter = seq_along(becod_design$trace), utility = becod_design$trace, method = "BE-COD"),
  tibble::tibble(iter = seq_along(nibod_design$trace), utility = nibod_design$trace, method = "NI-BOD")
)

p10 <- ggplot2::ggplot(trace_df, ggplot2::aes(x = iter, y = utility, color = method)) +
  ggplot2::geom_line(linewidth = 1.15) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::labs(
    title = "Coordinate exchange utility trace",
    x = "Iteration",
    y = "Utility"
  ) +
  theme_be

print(p10)
save_plot_both(p10, "fig_10_coordinate_exchange_trace", width = 8.3, height = 4.8)

# Plot 11: Residual diagnostics
D0_diag <- broom::augment(fit_D0)

p11a <- ggplot2::ggplot(D0_diag, ggplot2::aes(x = .fitted, y = .resid)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_point(size = 2.3, alpha = 0.8) +
  ggplot2::labs(
    title = "Residuals vs fitted (historical model)",
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_be

p11b <- ggplot2::ggplot(D0_diag, ggplot2::aes(sample = .std.resid)) +
  ggplot2::stat_qq(size = 1.8, alpha = 0.8) +
  ggplot2::stat_qq_line() +
  ggplot2::labs(
    title = "Normal Q-Q plot (historical model)",
    x = "Theoretical quantiles",
    y = "Standardized residuals"
  ) +
  theme_be

print(p11a)
print(p11b)
save_plot_both(p11a, "fig_11_residuals_vs_fitted_D0", width = 7.2, height = 4.8)
save_plot_both(p11b, "fig_12_qqplot_D0", width = 7.2, height = 4.8)

# Plot 12: Future posterior SD comparison
p12 <- ggplot2::ggplot(future_sd_table, ggplot2::aes(x = estimand, y = posterior_sd, fill = method)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.65) +
  ggplot2::labs(
    title = "Posterior SD of target estimands under future design options",
    x = NULL,
    y = "Posterior SD",
    fill = "Design"
  ) +
  theme_be +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

print(p12)
save_plot_both(p12, "fig_13_future_target_sd_barplot", width = 10, height = 5.2)

# ----------------------------
# 19. Save all main objects
# ----------------------------
saveRDS(
  list(
    tg = tg,
    tg_split = tg_split,
    D0 = D0,
    D1 = D1,
    fit_D0 = fit_D0,
    sigma2_hat = sigma2_hat,
    sigma_hat = sigma_hat,
    m_comm = m_comm,
    V_comm = V_comm,
    m_vague = m_vague,
    V_vague = V_vague,
    omega0 = omega0,
    L_targets = L_targets,
    post_comm = post_comm,
    post_vague = post_vague,
    post_mix = post_mix,
    contrast_table = contrast_table,
    posterior_weight_table = posterior_weight_table,
    becod_design = becod_design,
    nibod_design = nibod_design,
    design_comparison = design_comparison,
    utility_comparison = utility_comparison,
    future_sd_table = future_sd_table,
    curve_df = curve_df
  ),
  file = file.path(obj_dir, "becod_toothgrowth_analysis_objects.rds")
)

writeLines(capture.output(sessionInfo()), con = file.path(obj_dir, "sessionInfo.txt"))

# ----------------------------
# 20. Final console summary
# ----------------------------
cat("\n============================================================\n")
cat("FINAL SUMMARY\n")
cat("============================================================\n")
cat("Historical sigma^2 estimate:", round(sigma2_hat, 6), "\n")
cat("Posterior robust weight on commensurate component:", round(post_mix$w1, 6), "\n")
cat("Posterior robust weight on vague component:", round(post_mix$w2, 6), "\n\n")

cat("Target contrast summary:\n")
print(contrast_table)

cat("\nFuture design utility comparison:\n")
print(utility_comparison)

cat("\nFuture design allocations:\n")
print(design_comparison)

cat("\nFuture target posterior SDs:\n")
print(future_sd_table)

cat("\nAll tables saved to:\n", tab_dir, "\n")
cat("All figures saved to:\n", fig_dir, "\n")
cat("All objects saved to:\n", obj_dir, "\n")
cat("============================================================\n")