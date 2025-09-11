psi_alpha <- function(uis, mu_alpha, Sigma_alpha){
  lin <- -uis %*% mu_alpha
  quad <- 0.5 * rowSums((uis %*% Sigma_alpha) * uis)
  res <- exp(lin + quad)
  
  as.vector(res)
}

psi_tau <- function(wijs, mu_tau, Sigma_tau){
  lin <- -wijs %*% mu_tau
  quad <- 0.5 * rowSums((wijs %*% Sigma_tau) * wijs)
  res <- exp(lin + quad)
  
  as.vector(res)
}

psi_omega <- function(mu_omegas, Sigma_omegas){
  lin <- -mu_omegas
  quad <- 0.5 * Sigma_omegas
  res <- exp(lin + quad)
  
  as.vector(res)
}

beta_update_ncr <- function(X_mat, y_vec, id, uq_id, mu_nu, psi_taus, psi_omegas, beta_prior){
  psi_omega_ij <- psi_omegas[match(id, uq_id)]
  X_weights <- psi_taus * psi_omega_ij
  nuij <- mu_nu[match(id, uq_id)]
  sum_xx <- t(X_mat) %*% (X_weights * X_mat)
  
  weight_resid <- (y_vec - nuij) * X_weights
  sum_yx <- t(X_mat) %*% weight_resid
  
  prior_mat <- diag(1 / beta_prior, ncol(X_mat))
  
  Sigma_beta_new <- solve(sum_xx + prior_mat)
  mu_beta_new <- Sigma_beta_new %*% sum_yx
  
  list(mu_beta_new = mu_beta_new, Sigma_beta_new = Sigma_beta_new)
}

nui_grad_hes_ncr <- function(idx, X_mat, y_vec, uq_id, id_list, mu_beta, psi_taus, psi_omegas, psi_alphas, mu_nu){
  idx_char <- as.character(uq_id[idx]) 
  cur_rows <- id_list[[idx_char]]
  
  yi <- y_vec[cur_rows]
  xi <- X_mat[cur_rows, ]
  psi_taui <- psi_taus[cur_rows]
  psi_omegai <- psi_omegas[idx]
  psi_alphai <- psi_alphas[idx]
  mu_nui <- mu_nu[idx]
  
  psi_i <- psi_taui * psi_omegai
  resid_i <- yi - xi %*% mu_beta - mu_nui
  
  lik_grad <- sum(psi_i * resid_i)
  prior_grad <- -psi_alphai * mu_nui
  tol_grad <- lik_grad + prior_grad
  
  lik_hess <- -sum(psi_i)
  prior_hess <- -psi_alphai
  tol_hess <- lik_hess + prior_hess
  
  c(grad = tol_grad, hessian = tol_hess)
}

nu_update_ncr <- function(X_mat, y_vec, uq_id, id_list, mu_beta, psi_taus, psi_omegas, psi_alphas, mu_nu){
  num_ids <- length(uq_id)
  nu_grad_hes <- sapply(1:num_ids, nui_grad_hes_ncr, X_mat = X_mat, y_vec = y_vec, uq_id = uq_id, id_list = id_list, 
                        mu_beta = mu_beta, psi_taus = psi_taus, psi_omegas = psi_omegas, 
                        psi_alphas = psi_alphas, mu_nu = mu_nu)
  nu_grad <- nu_grad_hes[1,]
  nu_hes <- nu_grad_hes[2,]
  
  Sigma_nu_new <- -1 / nu_hes
  mu_nu_new <- mu_nu + Sigma_nu_new * nu_grad
  
  list(mu_nu_new = mu_nu_new, Sigma_nu_new = Sigma_nu_new)
}

alpha_update_ncr <- function(uis, psi_alphas, mu_nu, Sigma_nu, mu_alpha, alpha_prior){
  nu_sq <- mu_nu^2 + Sigma_nu
  
  Sigma_alpha_w <- 0.5 * psi_alphas * nu_sq
  Sigma_alpha_sum <- t(uis) %*% (Sigma_alpha_w * uis)
  prior_mat <- diag(1 / alpha_prior, ncol(uis))
  Sigma_alpha_new <- solve(Sigma_alpha_sum + prior_mat)
  
  mu_alpha_w <- 0.5 * (psi_alphas * nu_sq - 1)
  mu_alpha_sum <- t(uis) %*% mu_alpha_w
  prior_grad <- -(1 / alpha_prior) * mu_alpha
  tol_grad <- mu_alpha_sum + prior_grad
  mu_alpha_new <- mu_alpha + Sigma_alpha_new %*% tol_grad
  
  list(mu_alpha_new = mu_alpha_new, Sigma_alpha_new = Sigma_alpha_new)
}

tau_update_ncr <- function(X_mat, y_vec, w_mat, id, uq_id, mu_nu, Sigma_nu, mu_beta, Sigma_beta, psi_taus, 
                           psi_omegas, mu_tau, tau_prior){
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  psi_omega_ij <- psi_omegas[match(id, uq_id)]
  
  fixed <- y_vec - (X_mat %*% mu_beta) - mu_nuis
  var_beta <- rowSums((X_mat %*% Sigma_beta) * X_mat)
  err_sq <- fixed^2 + var_beta + Sigma_nuis
  hij <- as.vector(err_sq * psi_omega_ij)
  
  hess_w <- 0.5 * psi_taus * hij
  sum_ww <- t(w_mat) %*% (hess_w * w_mat)
  prior_mat <- diag(1/tau_prior, ncol(w_mat))
  Sigma_tau_new <- solve(sum_ww + prior_mat)
  
  grad_w <- 0.5 * (psi_taus * hij - 1)
  lik_grad <- t(w_mat) %*% grad_w
  prior_grad <- -(1/tau_prior) * mu_tau
  tol_grad <- lik_grad + prior_grad
  mu_tau_new <- mu_tau + Sigma_tau_new %*% tol_grad
  
  list(mu_tau_new = mu_tau_new, Sigma_tau_new = Sigma_tau_new)
}

omega_grad_hess_ncr <- function(idx, X_mat, y_vec, w_mat, uq_id, id_list, mu_beta, Sigma_beta, mu_nu, Sigma_nu, 
                                psi_taus, mu_inv_sigma, psi_omegas, mu_omegas){
  idx_char <- as.character(uq_id[idx])
  cur_rows <- id_list[[idx_char]]
  
  yi <- y_vec[cur_rows]
  xi <- X_mat[cur_rows, ]
  mu_nui <- mu_nu[idx]
  Sigma_nui <- Sigma_nu[idx]
  mu_omegai <- mu_omegas[idx]
  psi_omegai <- psi_omegas[idx]
  psi_taui <- psi_taus[cur_rows]
  
  fixed_i <- yi - xi %*% mu_beta - mu_nui
  var_beta_i <- rowSums((xi %*% Sigma_beta) * xi)
  err_sq_i <- fixed_i^2 + var_beta_i + Sigma_nui
  hij <- err_sq_i * psi_taui
  
  lik_grad <- 0.5 * psi_omegai * sum(hij) - 0.5 * length(yi)
  prior_grad <- -mu_inv_sigma * mu_omegai
  tol_grad <- lik_grad + prior_grad
  
  lik_hess <- -0.5 * psi_omegai * sum(hij)
  prior_hess <- -mu_inv_sigma
  tol_hess <- lik_hess + prior_hess
  
  c(grad = tol_grad, hessian = tol_hess)
}

omega_update_ncr <- function(X_mat, y_vec, w_mat, uq_id, id_list, mu_beta, Sigma_beta, mu_nu, Sigma_nu, mu_inv_sigma, 
                             psi_taus, psi_omegas, mu_omegas){
  num_ids <- length(uq_id)
  omega_grad_hess <- sapply(1:num_ids, omega_grad_hess_ncr, X_mat = X_mat, y_vec = y_vec, w_mat = w_mat, 
                            uq_id = uq_id, id_list = id_list, mu_beta = mu_beta, Sigma_beta = Sigma_beta, 
                            mu_nu = mu_nu, Sigma_nu = Sigma_nu, psi_taus = psi_taus, mu_inv_sigma = mu_inv_sigma, 
                            psi_omegas = psi_omegas, mu_omegas = mu_omegas)
  
  omega_grad <- omega_grad_hess[1, ]
  omega_hess <- omega_grad_hess[2, ]
  
  Sigma_omega_new <- -1 / omega_hess
  mu_omega_new <- mu_omegas + Sigma_omega_new * omega_grad
  
  list(mu_omega_new = as.vector(mu_omega_new), Sigma_omega_new = as.vector(Sigma_omega_new))
}

sigma_inv_update_ncr <- function(nids, mu_omega, Sigma_omega, mu_inv_a){
  sigma_shape <- (nids+1)/2
  sum_omega_sq <- sum(Sigma_omega + mu_omega^2)
  rate_sigma <- mu_inv_a + 0.5 * sum_omega_sq
  
  mu_inv_sigma_new <- sigma_shape / rate_sigma
  mu_log_sigma <- log(rate_sigma) - digamma(sigma_shape)
  
  list(mu_inv_sigma_new = mu_inv_sigma_new, mu_log_sigma = mu_log_sigma)
}

a_inv_update_ncr <- function(mu_inv_sigma, A_omega){
  rate_a <- (1 / A_omega^2) + mu_inv_sigma
  mu_inv_a_new <- 1 / rate_a
  mu_log_a <- log(rate_a) - digamma(1)
  
  list(mu_inv_a_new = mu_inv_a_new, mu_log_a = mu_log_a)
}

elbo_cal_ncr <- function(X_mat, y_vec, u_mat, w_mat, mu_beta, Sigma_beta, mu_nu, Sigma_nu, 
                         mu_alpha, Sigma_alpha, mu_tau, Sigma_tau, mu_omega, Sigma_omega, 
                         mu_inv_sigma, mu_log_sigma, mu_inv_a, mu_log_a, 
                         psi_ij_tau, psi_i_alpha, psi_i_omega,
                         beta_prior, alpha_prior, tau_prior, A_prior,
                         id, uq_id){
  num_ids <- length(uq_id)
  
  mu_omega_expanded <- mu_omega[match(id, uq_id)]
  psi_omega_expanded <- psi_i_omega[match(id, uq_id)]
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  
  fixed <- y_vec - (X_mat %*% mu_beta) - mu_nuis
  var_beta <- rowSums((X_mat %*% Sigma_beta) * X_mat)
  err_sq <- fixed^2 + var_beta + Sigma_nuis
  hij <- as.vector(err_sq * psi_ij_tau)
  
  log_var_mean_expanded <- as.vector(w_mat %*% mu_tau) + mu_omega_expanded
  log_p_y <- -0.5 * sum(log_var_mean_expanded) - 0.5 * sum(psi_omega_expanded * hij)
  
  log_p_beta <- -0.5 * (1 / beta_prior) * (sum(mu_beta^2) + sum(diag(Sigma_beta)))
  log_p_tau <- -0.5 * (1 / tau_prior) * (sum(mu_tau^2) + sum(diag(Sigma_tau)))
  log_p_alpha <- -0.5 * (1 / alpha_prior) * (sum(mu_alpha^2) + sum(diag(Sigma_alpha)))
  
  log_p_nu <- -0.5 * sum(u_mat %*% mu_alpha + psi_i_alpha * (mu_nu^2 + Sigma_nu))
  log_p_omega <- -0.5 * mu_inv_sigma * sum(mu_omega^2 + Sigma_omega) 
  
  log_p_sigma <- -0.5 * mu_log_a - 1.5 * mu_log_sigma - mu_inv_a * mu_inv_sigma
  log_p_a <- -1.5 * mu_log_a - mu_inv_a * (1 / A_prior^2)
  
  log_q_beta <- 0.5 * determinant(Sigma_beta, logarithm = TRUE)$modulus
  log_q_tau <- 0.5 * determinant(Sigma_tau, logarithm = TRUE)$modulus
  log_q_alpha <- 0.5 * determinant(Sigma_alpha, logarithm = TRUE)$modulus
  log_q_nu <- 0.5 * sum(log(Sigma_nu))
  log_q_omega <- 0.5 * sum(log(Sigma_omega)) 
  
  rate_sigma <- mu_inv_a + 0.5 * sum(Sigma_omega + mu_omega^2) 
  log_q_sigma <- -lgamma((num_ids+1)/2) + ((num_ids+1)/2)*log(rate_sigma)
  
  rate_a <- (1 / A_prior^2) + mu_inv_sigma
  log_q_a <- -lgamma(1) + log(rate_a)
  
  p_terms <- log_p_y + log_p_beta + log_p_tau + log_p_alpha + log_p_nu + log_p_omega + log_p_sigma + log_p_a
  q_terms <- log_q_beta + log_q_alpha + log_q_tau + log_q_nu + log_q_omega + log_q_sigma + log_q_a
  
  elbo <- p_terms - q_terms
  elbo
}

# Riesby data set
st_time <- Sys.time()
riesby <- read.table("RIESBY.DAT.txt", na.strings = ".")
colnames(riesby) <- c("id", "hamd", "intcpt", "week", "endog", "endweek")
riesby <- na.omit(riesby)

X_matrix <- model.matrix(~ week + endog + endweek, data = riesby)
y_vector <- as.matrix(riesby$hamd)
w_matrix <- model.matrix(~ week + endog, data = riesby)

riesby_id <- riesby[!duplicated(riesby$id), ]
ids <- riesby$id
u_matrix <- model.matrix(~ endog, data = riesby_id)
unique_ids <- riesby_id$id
num_ids <- length(unique_ids)

id_indices <- split(1:nrow(riesby), riesby$id)

# Priors
beta_prior <- 100000; alpha_prior <- 100000; tau_prior <- 100000; A_prior <- 100000

library(lme4)
lme_fit <- lmer(hamd ~ week + endog + endweek + (1 | id), data = riesby)
beta_init <- fixef(lme_fit)
nu_init <- as.vector(ranef(lme_fit)$id[, "(Intercept)"])
resid_init <- residuals(lme_fit) |> as.vector()

# Initialize variational parameters
p_beta <- ncol(X_matrix); mu_beta_q <- beta_init; Sigma_beta_q <- diag(0.1, p_beta)
mu_nu_i_q <- nu_init; Sigma_nu_i_q <- rep(0.1, num_ids)

p_alpha <- ncol(u_matrix); mu_alpha_q <- matrix(0, nrow = p_alpha, ncol = 1); Sigma_alpha_q <- diag(0.1, p_alpha)

p_tau <- ncol(w_matrix); lm_tau <- lm(log(resid_init^2 + 1e-4) ~ w_matrix - 1); mu_tau_q <- coef(lm_tau); Sigma_tau_q <- diag(0.1, p_tau)

# Correctly initialize omega_i at the group level
mu_omega_q <- tapply(residuals(lm_tau), ids, mean); Sigma_omega_q <- rep(0.1, num_ids)
mu_inv_sigma_omega_q <- 1; mu_log_sigma_omega_q <- 0
mu_inv_a_omega_q <- 1; mu_log_a_omega_q <- 0

# Initialize psi's
psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
psi_i_omega <- psi_omega(mu_omega_q, Sigma_omega_q)

max_iter <- 1000
tolerance <- 1e-7
elbo_history <- c()
for (i in 1:max_iter){
  beta_new <- beta_update_ncr(X_matrix, y_vector, ids, unique_ids, mu_nu_i_q, psi_ij_tau, psi_i_omega, beta_prior)
  mu_beta_q <- beta_new$mu_beta_new; Sigma_beta_q <- beta_new$Sigma_beta_new
  
  nu_new <- nu_update_ncr(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, psi_ij_tau, psi_i_omega, psi_i_alpha, mu_nu_i_q)
  mu_nu_i_q <- nu_new$mu_nu_new; Sigma_nu_i_q <- nu_new$Sigma_nu_new
  
  alpha_new <- alpha_update_ncr(u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)
  mu_alpha_q <- alpha_new$mu_alpha_new; Sigma_alpha_q <- alpha_new$Sigma_alpha_new
  psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  
  tau_new <- tau_update_ncr(X_matrix, y_vector, w_matrix, ids, unique_ids, mu_nu_i_q, Sigma_nu_i_q, 
                            mu_beta_q, Sigma_beta_q, psi_ij_tau, psi_i_omega, mu_tau_q, tau_prior)
  mu_tau_q <- tau_new$mu_tau_new; Sigma_tau_q <- tau_new$Sigma_tau_new
  psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  
  omega_new <- omega_update_ncr(X_matrix, y_vector, w_matrix, unique_ids, id_indices, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                                mu_inv_sigma_omega_q, psi_ij_tau, psi_i_omega, mu_omega_q)
  mu_omega_q <- omega_new$mu_omega_new; Sigma_omega_q <- omega_new$Sigma_omega_new
  psi_i_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  sigma_omega_new <- sigma_inv_update_ncr(num_ids, mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q)
  mu_inv_sigma_omega_q <- sigma_omega_new$mu_inv_sigma_new; mu_log_sigma_omega_q <- sigma_omega_new$mu_log_sigma
  
  a_omega_new <- a_inv_update_ncr(mu_inv_sigma_omega_q, A_prior)
  mu_inv_a_omega_q <- a_omega_new$mu_inv_a_new; mu_log_a_omega_q <- a_omega_new$mu_log_a
  
  elbo <- elbo_cal_ncr(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                       mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q, 
                       mu_inv_sigma_omega_q, mu_log_sigma_omega_q, mu_inv_a_omega_q, mu_log_a_omega_q, 
                       psi_ij_tau, psi_i_alpha, psi_i_omega,
                       beta_prior, alpha_prior, tau_prior, A_prior,
                       ids, unique_ids)
  elbo_history <- c(elbo_history, elbo)
  if ((i > 1) && (abs(elbo_history[i] - elbo_history[i-1]) < tolerance * abs(elbo_history[i]))){
    cat("CAVI converges at iteration ", i, "\n")
    break
  }
}
ed_time <- Sys.time()

print(beta_new)
print(alpha_new)
print(tau_new)
sigma_omega_rate <- mu_inv_a_omega_q + 0.5 * sum(Sigma_omega_q + mu_omega_q^2)
# FIX #4: Use num_ids, not nobs, for the denominator
sigma_omega_mean <- sigma_omega_rate / ((num_ids + 1)/2 - 1)
print(paste("Std Dev of omega_i:", sqrt(sigma_omega_mean)))
print(ed_time - st_time)

# Harvard data set
st_time <- Sys.time()
indat=read.table("Dataset_HealthBehavAcadPerfAffect.dat", header = FALSE, 
                 col.names=c("id", "day", "sex", "age", "sem", "sq", "physact", "pa", "na", "lga", 
                             "exam", "hsg", "bdi", "day_c"), na.strings="-99")
indat <- subset(indat, !is.na(indat$pa))
indat <- subset(indat, !is.na(indat$sex))

X_matrix <- model.matrix(~ day_c, data = indat)
y_vector <- as.matrix(indat$pa)
w_matrix <- model.matrix(~ day_c, data = indat)

indat_id <- indat[!duplicated(indat$id), ]
ids <- indat$id
u_matrix <- model.matrix(~ sex, data = indat_id)
unique_ids <- indat_id$id
num_ids <- length(unique_ids)

id_indices <- split(1:nrow(indat), indat$id)

# Priors
beta_prior <- 100000; alpha_prior <- 100000; tau_prior <- 100000; A_prior <- 100000

lme_fit <- lmer(pa ~ day_c + (1 | id), data = indat)
beta_init <- fixef(lme_fit)
nu_init <- as.vector(ranef(lme_fit)$id[, "(Intercept)"])
resid_init <- residuals(lme_fit) |> as.vector()

# Initialize variational parameters
p_beta <- ncol(X_matrix); mu_beta_q <- beta_init; Sigma_beta_q <- diag(0.1, p_beta)
mu_nu_i_q <- nu_init; Sigma_nu_i_q <- rep(0.1, num_ids)

p_alpha <- ncol(u_matrix); mu_alpha_q <- matrix(0, nrow=p_alpha, ncol=1); Sigma_alpha_q <- diag(0.1, p_alpha)
p_tau <- ncol(w_matrix); lm_tau <- lm(log(resid_init^2 + 1e-4) ~ w_matrix - 1); mu_tau_q <- coef(lm_tau); Sigma_tau_q <- diag(0.1, p_tau)

mu_omega_q <- tapply(residuals(lm_tau), ids, mean); Sigma_omega_q <- rep(0.1, num_ids)

mu_inv_sigma_omega_q <- 1; mu_log_sigma_omega_q <- 0
mu_inv_a_omega_q <- 1; mu_log_a_omega_q <- 0

# Initialize psi's
psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
psi_i_omega <- psi_omega(mu_omega_q, Sigma_omega_q)

max_iter <- 1000
tolerance <- 1e-7
elbo_history <- c()
for (i in 1:max_iter){
  beta_new <- beta_update_ncr(X_matrix, y_vector, ids, unique_ids, mu_nu_i_q, psi_ij_tau, psi_i_omega, beta_prior)
  mu_beta_q <- beta_new$mu_beta_new; Sigma_beta_q <- beta_new$Sigma_beta_new
  
  nu_new <- nu_update_ncr(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, psi_ij_tau, psi_i_omega, psi_i_alpha, mu_nu_i_q)
  mu_nu_i_q <- nu_new$mu_nu_new; Sigma_nu_i_q <- nu_new$Sigma_nu_new
  
  alpha_new <- alpha_update_ncr(u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)
  mu_alpha_q <- alpha_new$mu_alpha_new; Sigma_alpha_q <- alpha_new$Sigma_alpha_new
  psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  
  tau_new <- tau_update_ncr(X_matrix, y_vector, w_matrix, ids, unique_ids, mu_nu_i_q, Sigma_nu_i_q, 
                            mu_beta_q, Sigma_beta_q, psi_ij_tau, psi_i_omega, mu_tau_q, tau_prior)
  mu_tau_q <- tau_new$mu_tau_new; Sigma_tau_q <- tau_new$Sigma_tau_new
  psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  
  omega_new <- omega_update_ncr(X_matrix, y_vector, w_matrix, unique_ids, id_indices, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                                mu_inv_sigma_omega_q, psi_ij_tau, psi_i_omega, mu_omega_q)
  mu_omega_q <- omega_new$mu_omega_new; Sigma_omega_q <- omega_new$Sigma_omega_new
  psi_i_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  sigma_omega_new <- sigma_inv_update_ncr(num_ids, mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q)
  mu_inv_sigma_omega_q <- sigma_omega_new$mu_inv_sigma_new; mu_log_sigma_omega_q <- sigma_omega_new$mu_log_sigma
  
  a_omega_new <- a_inv_update_ncr(mu_inv_sigma_omega_q, A_prior)
  mu_inv_a_omega_q <- a_omega_new$mu_inv_a_new; mu_log_a_omega_q <- a_omega_new$mu_log_a
  
  elbo <- elbo_cal_ncr(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                       mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q, 
                       mu_inv_sigma_omega_q, mu_log_sigma_omega_q, 
                       mu_inv_a_omega_q, mu_log_a_omega_q, 
                       psi_ij_tau, psi_i_alpha, psi_i_omega, beta_prior, alpha_prior, tau_prior, A_prior,
                       ids, unique_ids)
  elbo_history <- c(elbo_history, elbo)
  if ((i > 1) && (abs(elbo_history[i] - elbo_history[i-1]) < tolerance * abs(elbo_history[i]))){
    cat("CAVI converges at iteration ", i, "\n")
    break
  }
}
ed_time <- Sys.time()

print(beta_new)
print(alpha_new)
print(tau_new)
sigma_omega_rate <- mu_inv_a_omega_q + 0.5 * sum(Sigma_omega_q + mu_omega_q^2)
sigma_omega_mean <- sigma_omega_rate / ((num_ids + 1)/2 - 1)
print(paste("Std Dev of omega_i:", sqrt(sigma_omega_mean)))
print(ed_time - st_time)

# Simulated data set
st_time <- Sys.time()
sim_dat = read.table("simulated_data.dat", header = TRUE)

X_matrix <- model.matrix(~ time + group, data = sim_dat)
y_vector <- as.matrix(sim_dat$y)
w_matrix <- model.matrix(~ time, data = sim_dat)

sim_dat_id <- sim_dat[!duplicated(sim_dat$id), ]
ids <- sim_dat$id
u_matrix <- model.matrix(~ group, data = sim_dat_id)
unique_ids <- sim_dat_id$id
num_ids <- length(unique_ids)

id_indices <- split(1:nrow(sim_dat), sim_dat$id)

# Priors
beta_prior <- 100000; alpha_prior <- 100000; tau_prior <- 100000; A_prior <- 100000

lme_fit <- lmer(y ~ time + group + (1 | id), data = sim_dat)
beta_init <- fixef(lme_fit)
nu_init <- as.vector(ranef(lme_fit)$id[, "(Intercept)"])
resid_init <- residuals(lme_fit) |> as.vector()

# Initialize variational parameters
p_beta <- ncol(X_matrix); mu_beta_q <- beta_init; Sigma_beta_q <- diag(0.1, p_beta)
mu_nu_i_q <- nu_init; Sigma_nu_i_q <- rep(0.1, num_ids)

p_alpha <- ncol(u_matrix); mu_alpha_q <- matrix(0, nrow=p_alpha, ncol=1); Sigma_alpha_q <- diag(0.1, p_alpha)
p_tau <- ncol(w_matrix); lm_tau <- lm(log(resid_init^2 + 1e-4) ~ w_matrix - 1); mu_tau_q <- coef(lm_tau); Sigma_tau_q <- diag(0.1, p_tau)

mu_omega_q <- tapply(residuals(lm_tau), ids, mean); Sigma_omega_q <- rep(0.1, num_ids)

mu_inv_sigma_omega_q <- 1; mu_log_sigma_omega_q <- 0
mu_inv_a_omega_q <- 1; mu_log_a_omega_q <- 0

# Initialize psi's
psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
psi_i_omega <- psi_omega(mu_omega_q, Sigma_omega_q)

max_iter <- 10000
tolerance <- 1e-7
elbo_history <- c()
for (i in 1:max_iter){
  beta_new <- beta_update_ncr(X_matrix, y_vector, ids, unique_ids, mu_nu_i_q, psi_ij_tau, psi_i_omega, beta_prior)
  mu_beta_q <- beta_new$mu_beta_new; Sigma_beta_q <- beta_new$Sigma_beta_new
  
  nu_new <- nu_update_ncr(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, psi_ij_tau, psi_i_omega, psi_i_alpha, mu_nu_i_q)
  mu_nu_i_q <- nu_new$mu_nu_new; Sigma_nu_i_q <- nu_new$Sigma_nu_new
  
  alpha_new <- alpha_update_ncr(u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)
  mu_alpha_q <- alpha_new$mu_alpha_new; Sigma_alpha_q <- alpha_new$Sigma_alpha_new
  psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  
  tau_new <- tau_update_ncr(X_matrix, y_vector, w_matrix, ids, unique_ids, mu_nu_i_q, Sigma_nu_i_q, 
                            mu_beta_q, Sigma_beta_q, psi_ij_tau, psi_i_omega, mu_tau_q, tau_prior)
  mu_tau_q <- tau_new$mu_tau_new; Sigma_tau_q <- tau_new$Sigma_tau_new
  psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  
  omega_new <- omega_update_ncr(X_matrix, y_vector, w_matrix, unique_ids, id_indices, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                                mu_inv_sigma_omega_q, psi_ij_tau, psi_i_omega, mu_omega_q)
  mu_omega_q <- omega_new$mu_omega_new; Sigma_omega_q <- omega_new$Sigma_omega_new
  psi_i_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  sigma_omega_new <- sigma_inv_update_ncr(num_ids, mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q)
  mu_inv_sigma_omega_q <- sigma_omega_new$mu_inv_sigma_new; mu_log_sigma_omega_q <- sigma_omega_new$mu_log_sigma
  
  a_omega_new <- a_inv_update_ncr(mu_inv_sigma_omega_q, A_prior)
  mu_inv_a_omega_q <- a_omega_new$mu_inv_a_new; mu_log_a_omega_q <- a_omega_new$mu_log_a
  
  elbo <- elbo_cal_ncr(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                       mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q, 
                       mu_inv_sigma_omega_q, mu_log_sigma_omega_q, 
                       mu_inv_a_omega_q, mu_log_a_omega_q, 
                       psi_ij_tau, psi_i_alpha, psi_i_omega, beta_prior, alpha_prior, tau_prior, A_prior,
                       ids, unique_ids)
  elbo_history <- c(elbo_history, elbo)
  if ((i > 1) && (abs(elbo_history[i] - elbo_history[i-1]) < tolerance * abs(elbo_history[i]))){
    cat("CAVI converges at iteration ", i, "\n")
    break
  }
}
ed_time <- Sys.time()

print(beta_new)
print(alpha_new)
print(tau_new)
sigma_omega_rate <- mu_inv_a_omega_q + 0.5 * sum(Sigma_omega_q + mu_omega_q^2)
sigma_omega_mean <- sigma_omega_rate / ((num_ids + 1)/2 - 1)
print(paste("Std Dev of omega_i:", sqrt(sigma_omega_mean)))
print(ed_time - st_time)