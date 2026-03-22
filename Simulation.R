############################################################
# BE-COD Simulation Verification (R)
# - Compares BE-COD (Empirical-Bayes prior + targeted Bayes design)
#   vs NI-BOD (same targeted Bayes design but vague prior)
#   vs BAL (balanced allocation across candidate points)
# - Saves ALL outputs (plots + tables) to an auto-created folder
# - Prints key tables to console for copy/paste feedback
############################################################

# ---------------------------
# 0) Packages (auto-install)
# ---------------------------
pkgs <- c(
  "dplyr","tidyr","purrr","tibble","ggplot2","Matrix","MASS","knitr"
)
to_install <- setdiff(pkgs, rownames(installed.packages()))
if(length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(Matrix); library(MASS); library(knitr)
})

# ---------------------------
# 1) Output folder setup
# ---------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTDIR <- file.path("BE_COD_sim_outputs", paste0("run_", timestamp))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "figs"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "tables"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "rds"), showWarnings = FALSE)

cat("\n[Output folder]\n", OUTDIR, "\n\n")

# ---------------------------
# 2) Core model specification
# ---------------------------
# Candidate design points: supplement S in {VC, OJ} and dose d in {0.25,0.5,1,2,4}
# Model: mu(s,z)=β0 + βS*s + β1*z + β2*z^2 + βS1*(s*z) + βS2*(s*z^2)
# where s=1 for OJ and 0 for VC, z=log(dose)

make_candidates <- function(doses = c(0.25, 0.5, 1, 2, 4)) {
  expand.grid(
    supp = c("VC","OJ"),
    dose = doses,
    stringsAsFactors = FALSE
  ) %>% as_tibble() %>%
    mutate(s = ifelse(supp == "OJ", 1, 0),
           z = log(dose))
}

f_row <- function(s, z) {
  # β = (β0, βS, β1, β2, βS1, βS2)
  c(1, s, z, z^2, s*z, s*z^2)
}

design_matrix_from_alloc <- function(cand, alloc) {
  # alloc is integer vector length n with indices into cand (1..N)
  X <- do.call(rbind, lapply(alloc, function(j) {
    with(cand[j,], f_row(s, z))
  }))
  colnames(X) <- c("b0","bS","b1","b2","bS1","bS2")
  X
}

# Target estimands (contrasts) as L * beta:
# Interaction at d=0.5,1,2: Δ(d)=mu(OJ,d)-mu(VC,d)=βS + βS1*z + βS2*z^2
# Plus interaction-shape parameters: βS1 and βS2 directly (optional but informative)
make_L <- function(d_targets = c(0.5, 1, 2)) {
  L_list <- lapply(d_targets, function(d) {
    z <- log(d)
    c(0, 1, 0, 0, z, z^2)  # picks βS + βS1 z + βS2 z^2
  })
  L <- do.call(rbind, L_list)
  rownames(L) <- paste0("Delta_d", d_targets)
  # add rows for shape (optional)
  L2 <- rbind(L,
              "bS1" = c(0,0,0,0,1,0),
              "bS2" = c(0,0,0,0,0,1))
  L2
}

# ---------------------------
# 3) Bayesian updating (Gaussian, sigma known)
# ---------------------------
posterior_cov <- function(V0, X, sigma2) {
  # Vpost = (V0^{-1} + X'X/sigma2)^{-1}
  V0i <- solve(V0)
  XtX <- crossprod(X)
  Vpost <- solve(V0i + XtX / sigma2)
  Vpost
}

posterior_mean <- function(m0, V0, X, y, sigma2) {
  # mpost = Vpost * (V0^{-1} m0 + X'y/sigma2)
  V0i <- solve(V0)
  Vpost <- solve(V0i + crossprod(X)/sigma2)
  rhs <- V0i %*% m0 + crossprod(X, y) / sigma2
  as.vector(Vpost %*% rhs)
}

utility_targeted_logdet <- function(m0, V0, sigma2, cand, alloc, L) {
  # design utility uses posterior covariance only (mean irrelevant for covariance)
  X <- design_matrix_from_alloc(cand, alloc)
  Vpost <- posterior_cov(V0, X, sigma2)
  S <- L %*% Vpost %*% t(L)
  # numerical stabilization
  S <- as.matrix(nearPD(S, corr = FALSE)$mat)
  val <- -as.numeric(determinant(S, logarithm = TRUE)$modulus)
  val
}

# ---------------------------
# 4) Exact design optimization (coordinate exchange)
# ---------------------------
coord_exchange <- function(cand, n, m0, V0, sigma2, L,
                           n_iter = 15, n_restarts = 10, seed = 1, verbose = FALSE) {
  set.seed(seed)
  N <- nrow(cand)
  
  best_alloc <- sample.int(N, size = n, replace = TRUE)
  best_u <- utility_targeted_logdet(m0, V0, sigma2, cand, best_alloc, L)
  
  for(r in seq_len(n_restarts)) {
    alloc <- sample.int(N, size = n, replace = TRUE)
    u_curr <- utility_targeted_logdet(m0, V0, sigma2, cand, alloc, L)
    
    for(it in seq_len(n_iter)) {
      improved <- FALSE
      for(i in seq_len(n)) {
        u_best_i <- u_curr
        best_j <- alloc[i]
        # try all candidate points for run i
        for(j in seq_len(N)) {
          if(j == alloc[i]) next
          alloc_try <- alloc
          alloc_try[i] <- j
          u_try <- utility_targeted_logdet(m0, V0, sigma2, cand, alloc_try, L)
          if(u_try > u_best_i + 1e-10) {
            u_best_i <- u_try
            best_j <- j
          }
        }
        if(best_j != alloc[i]) {
          alloc[i] <- best_j
          u_curr <- u_best_i
          improved <- TRUE
        }
      }
      if(verbose) cat(sprintf(" restart=%d iter=%d utility=%.4f improved=%s\n", r, it, u_curr, improved))
      if(!improved) break
    }
    
    if(u_curr > best_u) {
      best_u <- u_curr
      best_alloc <- alloc
    }
  }
  
  list(alloc = best_alloc, utility = best_u)
}

# Balanced design (as close as possible)
balanced_alloc <- function(cand, n, seed = 1) {
  set.seed(seed)
  N <- nrow(cand)
  base <- rep(1:N, length.out = n)
  sample(base, size = n, replace = FALSE)
}

# ---------------------------
# 5) Empirical Bayes prior from K historical studies
# ---------------------------
simulate_study <- function(cand, alloc, beta_true, sigma2) {
  X <- design_matrix_from_alloc(cand, alloc)
  mu <- as.vector(X %*% beta_true)
  y <- rnorm(nrow(X), mean = mu, sd = sqrt(sigma2))
  list(X = X, y = y)
}

ols_fit <- function(X, y) {
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  beta_hat <- solve(XtX, Xty)
  # residual variance estimate
  n <- nrow(X); p <- ncol(X)
  res <- y - as.vector(X %*% beta_hat)
  s2 <- sum(res^2) / max(1, (n - p))
  V_hat <- s2 * solve(XtX)
  list(beta_hat = as.vector(beta_hat), V_hat = V_hat, s2 = s2)
}

estimate_EB_prior <- function(hist_fits, omega = 0.7, kappa2 = 1e4) {
  # hist_fits: list of OLS fits from K historical studies
  B <- do.call(rbind, lapply(hist_fits, `[[`, "beta_hat"))
  beta_bar <- colMeans(B)
  # between-study covariance (regularize)
  S_between <- cov(B)
  if(any(!is.finite(S_between))) S_between <- diag(rep(0, ncol(B)))
  # average within-study covariance
  V_within <- Reduce(`+`, lapply(hist_fits, `[[`, "V_hat")) / length(hist_fits)
  
  V_comm <- S_between + V_within + 1e-8 * diag(ncol(B))
  
  # robust moment-matching (mixture with vague component)
  m0 <- omega * beta_bar + (1 - omega) * rep(0, length(beta_bar))
  V0 <- omega * V_comm + (1 - omega) * (kappa2 * diag(length(beta_bar)))
  
  list(m0 = m0, V0 = V0, beta_bar = beta_bar, V_comm = V_comm)
}

# ---------------------------
# 6) Simulation engine
# ---------------------------
set.seed(123)

cand <- make_candidates()
L <- make_L()

# True parameters (historical)
beta_hist <- c(
  b0  = 18,   # baseline (VC at dose=1 since z=0)
  bS  = 3.0,  # OJ vs VC baseline diff at z=0
  b1  = 6.0,  # dose effect (log-dose)
  b2  = -1.0, # curvature
  bS1 = 1.0,  # interaction slope (OJ vs VC changes with dose)
  bS2 = -0.2  # interaction curvature
)

# Current truth: allow drift (prior-data conflict scenario)
make_beta_current <- function(beta_hist, drift = c(0,0,0,0,0,0)) {
  beta_hist + drift
}

# Two scenarios: "commensurate" and "conflict"
scenarios <- tibble(
  scenario = c("commensurate", "conflict"),
  drift_bS  = c(0.0, -2.0),   # flip/attenuate supplement effect
  drift_bS1 = c(0.0, -1.2),   # alter interaction slope
  drift_bS2 = c(0.0,  0.4)    # alter curvature
) %>%
  mutate(drift = pmap(list(drift_bS, drift_bS1, drift_bS2), function(dbS, dbS1, dbS2) {
    c(0, dbS, 0, 0, dbS1, dbS2)
  }))

# Global simulation controls
M <- 250          # Monte Carlo replicates per scenario (increase later if you want)
K_hist <- 8       # number of historical/pilot studies to learn EB prior
n0 <- 36          # per historical study sample size
n1 <- 40          # current sample size (to design)
sigma2_true <- 9  # observation variance
omega_rob <- 0.7  # robust weight for commensurate component
kappa2_vague <- 1e4

# Methods:
# - BE_COD: design with EB prior (m0,V0) learned from historical fits
# - NI_BOD: design with vague prior (0, kappa2 I) (non-informative Bayes design)
# - BAL: balanced design (no optimization; allocations equal-ish)
methods <- c("BE_COD","NI_BOD","BAL")

# Helper: summarize design allocations
alloc_table <- function(cand, alloc) {
  cand[alloc, ] %>%
    count(supp, dose, name = "n") %>%
    complete(supp, dose, fill = list(n = 0)) %>%
    arrange(supp, dose)
}

# Main run for one replicate
run_one <- function(beta_current, seed_base) {
  
  # 1) simulate K historical studies (balanced allocations)
  hist_fits <- vector("list", K_hist)
  for(k in 1:K_hist) {
    alloc0 <- balanced_alloc(cand, n0, seed = seed_base + 1000*k)
    dat0 <- simulate_study(cand, alloc0, beta_true = beta_hist, sigma2 = sigma2_true)
    hist_fits[[k]] <- ols_fit(dat0$X, dat0$y)
  }
  EB <- estimate_EB_prior(hist_fits, omega = omega_rob, kappa2 = kappa2_vague)
  
  # For design, assume sigma2 known (use sigma2_true; you can plug-in mean(hist s2))
  sigma2_design <- sigma2_true
  
  # 2) build designs
  # 2a) BE-COD design (uses EB prior)
  des_BE <- coord_exchange(
    cand = cand, n = n1, m0 = EB$m0, V0 = EB$V0, sigma2 = sigma2_design, L = L,
    n_iter = 15, n_restarts = 10, seed = seed_base + 1, verbose = FALSE
  )
  
  # 2b) NI-BOD design (same utility but vague prior)
  m0_v <- rep(0, length(EB$m0))
  V0_v <- kappa2_vague * diag(length(EB$m0))
  des_NI <- coord_exchange(
    cand = cand, n = n1, m0 = m0_v, V0 = V0_v, sigma2 = sigma2_design, L = L,
    n_iter = 15, n_restarts = 10, seed = seed_base + 2, verbose = FALSE
  )
  
  # 2c) Balanced design
  alloc_BAL <- balanced_alloc(cand, n1, seed = seed_base + 3)
  
  designs <- list(
    BE_COD = des_BE$alloc,
    NI_BOD = des_NI$alloc,
    BAL    = alloc_BAL
  )
  
  # 3) simulate current data under each design and compute posterior for estimands
  truth_contrasts <- function(beta) as.vector(L %*% beta)
  
  truth <- truth_contrasts(beta_current)
  
  eval_one_method <- function(method_name) {
    alloc <- designs[[method_name]]
    dat1 <- simulate_study(cand, alloc, beta_true = beta_current, sigma2 = sigma2_true)
    
    if(method_name == "BE_COD") {
      m0 <- EB$m0; V0 <- EB$V0
    } else {
      m0 <- m0_v;  V0 <- V0_v
    }
    
    mpost <- posterior_mean(m0, V0, dat1$X, dat1$y, sigma2_true)
    Vpost <- posterior_cov(V0, dat1$X, sigma2_true)
    
    est  <- as.vector(L %*% mpost)
    Vc   <- L %*% Vpost %*% t(L)
    se   <- sqrt(pmax(0, diag(Vc)))
    
    # 95% credible interval (normal approx, since conjugate normal w/ known sigma^2)
    lo <- est - 1.96 * se
    hi <- est + 1.96 * se
    
    tibble(
      method = method_name,
      estimand = rownames(L),
      truth = truth,
      est = est,
      se = se,
      cover = (lo <= truth) & (truth <= hi),
      sqerr = (est - truth)^2
    )
  }
  
  res <- bind_rows(lapply(names(designs), eval_one_method))
  
  # attach design summaries too
  des_summ <- bind_rows(lapply(names(designs), function(mn) {
    tab <- alloc_table(cand, designs[[mn]])
    tab$method <- mn
    tab
  }))
  
  list(res = res, des_summ = des_summ, EB = EB)
}

# ---------------------------
# 7) Run simulation over scenarios
# ---------------------------
all_res <- list()
all_des <- list()

pb <- txtProgressBar(min = 0, max = nrow(scenarios)*M, style = 3)
ctr <- 0

for(ss in seq_len(nrow(scenarios))) {
  sc <- scenarios[ss, ]
  beta_cur <- make_beta_current(beta_hist, drift = sc$drift[[1]])
  
  for(m in 1:M) {
    ctr <- ctr + 1
    setTxtProgressBar(pb, ctr)
    
    out <- run_one(beta_current = beta_cur, seed_base = 10000*ss + m)
    
    all_res[[length(all_res)+1]] <- out$res %>% mutate(scenario = sc$scenario)
    all_des[[length(all_des)+1]] <- out$des_summ %>% mutate(scenario = sc$scenario)
  }
}
close(pb)

df <- bind_rows(all_res)
df_des <- bind_rows(all_des)

saveRDS(df, file.path(OUTDIR, "rds", "simulation_results.rds"))
saveRDS(df_des, file.path(OUTDIR, "rds", "design_allocations.rds"))

# ---------------------------
# 8) Summaries (tables printed + saved)
# ---------------------------
# Overall metrics by scenario/method/estimand
summ_estimand <- df %>%
  group_by(scenario, method, estimand) %>%
  summarise(
    RMSE = sqrt(mean(sqerr)),
    mean_SE = mean(se),
    coverage = mean(cover),
    .groups = "drop"
  )

summ_method <- df %>%
  group_by(scenario, method) %>%
  summarise(
    RMSE_avg = sqrt(mean(sqerr)),
    mean_SE_avg = mean(se),
    coverage_avg = mean(cover),
    .groups = "drop"
  )

# Print to console (copy/paste friendly)
cat("\n==============================\n")
cat("SUMMARY BY METHOD (avg over estimands)\n")
cat("==============================\n")
print(summ_method %>% arrange(scenario, method))
cat("\n\n")

cat("\n==============================\n")
cat("SUMMARY BY METHOD x ESTIMAND\n")
cat("==============================\n")
print(summ_estimand %>% arrange(scenario, estimand, method))
cat("\n\n")

# Also pretty kable (still console-friendly)
cat("\n--- kable: Summary by method (avg) ---\n")
print(knitr::kable(summ_method %>% arrange(scenario, method), digits = 4))
cat("\n\n")

cat("\n--- kable: Summary by method x estimand ---\n")
print(knitr::kable(summ_estimand %>% arrange(scenario, estimand, method), digits = 4))
cat("\n\n")

# Save tables
write.csv(summ_method, file.path(OUTDIR, "tables", "summary_by_method.csv"), row.names = FALSE)
write.csv(summ_estimand, file.path(OUTDIR, "tables", "summary_by_method_estimand.csv"), row.names = FALSE)

# Design allocation summaries (average counts)
des_avg <- df_des %>%
  group_by(scenario, method, supp, dose) %>%
  summarise(n_avg = mean(n), .groups = "drop")

cat("\n==============================\n")
cat("AVERAGE DESIGN ALLOCATIONS (counts per cell)\n")
cat("==============================\n")
print(des_avg %>% arrange(scenario, method, supp, dose))
cat("\n\n")

write.csv(des_avg, file.path(OUTDIR, "tables", "avg_design_allocations.csv"), row.names = FALSE)

# ---------------------------
# 9) Plots (saved + also printed)
# ---------------------------
# 9a) RMSE by estimand (boxplots over replicates)
df_rmse_rep <- df %>%
  group_by(scenario, method, estimand, .drop = FALSE) %>%
  mutate(err = sqrt(sqerr)) %>% ungroup()

p_rmse <- ggplot(df_rmse_rep, aes(x = method, y = err)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_grid(scenario ~ estimand, scales = "free_y") +
  labs(title = "Replicate-wise absolute error (sqrt squared error) by method",
       x = "Method", y = "Abs. error (|estimate - truth|)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(p_rmse)
ggsave(file.path(OUTDIR, "figs", "rmse_boxplots.png"), p_rmse, width = 13, height = 6, dpi = 180)
ggsave(file.path(OUTDIR, "figs", "rmse_boxplots.pdf"), p_rmse, width = 13, height = 6)

# 9b) Mean posterior SE (bar/point)
p_se <- summ_estimand %>%
  ggplot(aes(x = method, y = mean_SE, group = method)) +
  geom_point(size = 2) +
  facet_grid(scenario ~ estimand, scales = "free_y") +
  labs(title = "Mean posterior SD for each estimand (lower is better)",
       x = "Method", y = "Mean posterior SD") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(p_se)
ggsave(file.path(OUTDIR, "figs", "mean_se_by_estimand.png"), p_se, width = 13, height = 6, dpi = 180)
ggsave(file.path(OUTDIR, "figs", "mean_se_by_estimand.pdf"), p_se, width = 13, height = 6)

# 9c) Coverage (should be near 0.95 ideally)
p_cov <- summ_estimand %>%
  ggplot(aes(x = method, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  geom_col() +
  facet_grid(scenario ~ estimand) +
  ylim(0, 1) +
  labs(title = "Empirical 95% credible interval coverage (target ~0.95)",
       x = "Method", y = "Coverage") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(p_cov)
ggsave(file.path(OUTDIR, "figs", "coverage_by_estimand.png"), p_cov, width = 13, height = 6, dpi = 180)
ggsave(file.path(OUTDIR, "figs", "coverage_by_estimand.pdf"), p_cov, width = 13, height = 6)

# 9d) Average design allocation heatmaps
p_alloc <- des_avg %>%
  mutate(dose = factor(dose, levels = sort(unique(dose)))) %>%
  ggplot(aes(x = dose, y = supp, fill = n_avg)) +
  geom_tile() +
  facet_grid(scenario ~ method) +
  labs(title = "Average allocation counts per design cell (heatmap)",
       x = "Dose", y = "Supplement", fill = "Avg n") +
  theme_bw()

print(p_alloc)
ggsave(file.path(OUTDIR, "figs", "avg_alloc_heatmap.png"), p_alloc, width = 10, height = 5, dpi = 180)
ggsave(file.path(OUTDIR, "figs", "avg_alloc_heatmap.pdf"), p_alloc, width = 10, height = 5)

# ---------------------------
# 10) A quick “better or not” check printed clearly
# ---------------------------
# Compare BE_COD vs NI_BOD in each scenario using RMSE_avg
comp <- summ_method %>%
  select(scenario, method, RMSE_avg, mean_SE_avg, coverage_avg) %>%
  pivot_wider(names_from = method, values_from = c(RMSE_avg, mean_SE_avg, coverage_avg))

cat("\n==============================\n")
cat("DIRECT COMPARISON (wide table)\n")
cat("==============================\n")
print(comp)

write.csv(comp, file.path(OUTDIR, "tables", "direct_comparison_wide.csv"), row.names = FALSE)

cat("\n\n[Done] All plots + tables saved under:\n", OUTDIR, "\n")

############################################################
# Notes (edit if desired):
# - Increase M to 500 or 1000 for a tighter Monte Carlo estimate.
# - Increase n_restarts / n_iter if you want more aggressive optimization.
# - If you want the *analysis* prior to be identical across methods
#   (pure design comparison), set the posterior prior choice equal for all.
############################################################
