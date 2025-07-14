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

psi_r_nu <- function(mu_r, Sigma_r, mu_nus, Sigma_nus){
  lin <- -mu_r * mu_nus
  quad <- mu_nus^2 * Sigma_r + mu_r^2 * Sigma_nus
  res <- exp(lin + 0.5 * quad)
  as.vector(res)
}

beta_update <- function(X_mat, y_vec, id, uq_id, mu_nu, psi_taus, psi_omegas, psi_r_nus, beta_prior){
  nuij <- mu_nu[match(id, uq_id)]
  psi_r_nu_ij <- psi_r_nus[match(id, uq_id)]
  X_weights <- psi_taus * psi_omegas * psi_r_nu_ij
  sum_xx <- t(X_mat) %*% (X_weights * X_mat)
  
  weight_resid <- (y_vec - nuij) * X_weights
  sum_yx <- t(X_mat) %*% weight_resid
  
  prior_mat <- diag(1 / beta_prior, ncol(X_mat))
  Sigma_beta_new <- solve(sum_xx + prior_mat)
  mu_beta_new <- Sigma_beta_new %*% sum_yx
  
  list(mu_beta_new = mu_beta_new, Sigma_beta_new = Sigma_beta_new)
}

cal_Bij <- function(y_vec, X_mat, mu_beta, Sigma_beta, mu_nu, Sigma_nu, mu_r, Sigma_r){
  rss <- y_vec - (X_mat %*% mu_beta)
  t1 <- rss^2 + rowSums((X_mat %*% Sigma_beta) * X_mat)
  
  t2 <- -2 * rss * (mu_nu - mu_r * Sigma_nu)
  
  t3 <- Sigma_nu + (mu_nu - mu_r * Sigma_nu)^2

  Bij <- t1 + t2 + t3
  as.vector(Bij)
}

cal_Cij <- function(y_vec, X_mat, Bij, mu_beta, mu_nu, Sigma_nu, mu_r, Sigma_r){
  t1 <- (-mu_r + mu_nu * Sigma_r) * Bij
  t2 <- -2 * (y_vec - (X_mat %*% mu_beta))
  t3 <- 2 * (mu_nu - mu_r * Sigma_nu)
  Cij <- t1 + t2 + t3
  as.vector(Cij)
}

cal_Cij_deriv <- function(y_vec, X_mat, Bij, mu_beta, mu_nu, Sigma_nu, mu_r, Sigma_r) {
  
  deriv_multiplier <- Sigma_r
  deriv_B <- -2 * (y_vec - (X_mat %*% mu_beta)) + 2 * (mu_nu - mu_r * Sigma_nu)
  
  product_rule_part <- (deriv_multiplier * Bij) + ((-mu_r + mu_nu * Sigma_r) * deriv_B)
  
  rest_deriv <- 2
  
  return(as.vector(product_rule_part + rest_deriv))
}

nui_grad_hes <- function(idx, X_mat, y_vec, uq_id, id_list, mu_beta, Sigma_beta, mu_nu, Sigma_nu, mu_r, Sigma_r,
                         psi_taus, psi_omegas, psi_alphas, psi_nu_rs){
  idx_char <- as.character(uq_id[idx])
  cur_rows <- id_list[[idx_char]]
  
  yi <- y_vec[cur_rows]
  xi <- X_mat[cur_rows, ]
  
  psi_taui <- psi_taus[cur_rows]
  psi_omegai <- psi_omegas[cur_rows]
  psi_alphai <- psi_alphas[idx]
  psi_nu_ri <- psi_nu_rs[idx]
  
  mu_nui <- mu_nu[idx]
  Sigma_nui <- Sigma_nu[idx]
  Bi <- cal_Bij(yi, xi, mu_beta, Sigma_beta, mu_nui, Sigma_nui, mu_r, Sigma_r)
  Ci <- cal_Cij(yi, xi, Bi, mu_beta, mu_nui, Sigma_nui, mu_r, Sigma_r)
  ni <- length(yi)
  
  psi_i <- psi_taui * psi_omegai * psi_nu_ri
  
  logvar_grad <- -0.5 * ni * mu_r
  prior_grad <- -psi_alphai * mu_nui
  lik_grad <- -0.5 * sum(psi_i * Ci)
  tol_grad <- logvar_grad + prior_grad + lik_grad
  
  DCi <- cal_Cij_deriv(yi, xi, Bi, mu_beta, mu_nui, Sigma_nui, mu_r, Sigma_r)
  hess_in <- ((-mu_r + mu_nui * Sigma_r) * Ci) + DCi
  lik_hess <- -0.5 * sum(psi_i * hess_in)
  prior_hess <- -psi_alphai
  tol_hess <- lik_hess + prior_hess
  
  c(grad = tol_grad, hessian = tol_hess)
}

nu_update <- function(X_mat, y_vec, uq_id, id_list, mu_beta, Sigma_beta, mu_nu, Sigma_nu, mu_r, Sigma_r,
                      psi_taus, psi_omegas, psi_alphas, psi_nu_rs){
  num_ids <- length(uq_id)
  nu_grad_hes <- sapply(1:num_ids, nui_grad_hes, X_mat = X_mat, y_vec = y_vec, uq_id = uq_id, id_list = id_list, 
                        mu_beta = mu_beta, Sigma_beta = Sigma_beta, mu_nu = mu_nu, Sigma_nu = Sigma_nu, 
                        mu_r = mu_r, Sigma_r = Sigma_r, psi_taus = psi_taus, psi_omegas = psi_omegas, 
                        psi_alphas = psi_alphas, psi_nu_rs = psi_nu_rs)
  nu_grad <- nu_grad_hes[1,]
  nu_hes <- nu_grad_hes[2,]
  
  Sigma_nu_new <- -1 / nu_hes
  mu_nu_new <- mu_nu + Sigma_nu_new * nu_grad
  
  list(mu_nu_new = mu_nu_new, Sigma_nu_new = Sigma_nu_new)
}

alpha_update <- function(uis, psi_alphas, mu_nu, Sigma_nu, mu_alpha, alpha_prior){
  p <- ncol(uis)
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

cal_hij <- function(y_vec, X_mat, id, uq_id, mu_beta, Sigma_beta, mu_nu, Sigma_nu, mu_r, Sigma_r, psi_r_nus){
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  psi_r_nuis <- psi_r_nus[match(id, uq_id)]
  
  t1 <- (y_vec - X_mat %*% mu_beta)^2
  t2 <- rowSums((X_mat %*% Sigma_beta) * X_mat)
  t3 <- -2 * (mu_nuis - mu_r * Sigma_nuis) * (y_vec - X_mat %*% mu_beta)
  t4 <- Sigma_nuis + (mu_nuis - mu_r * Sigma_nuis)^2
  
  hij <- psi_r_nuis * (t1 + t2 + t3 + t4)
  as.vector(hij)
}

tau_update <- function(X_mat, y_vec, w_mat, id, uq_id, mu_nu, Sigma_nu, mu_beta, Sigma_beta, mu_r, Sigma_r, 
                       psi_taus, psi_omegas, psi_r_nus, mu_tau, tau_prior){
  p <- ncol(w_mat)
  
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  psi_r_nuis <- psi_r_nus[match(id, uq_id)]
  
  hij <- cal_hij(y_vec, X_mat, id, uq_id, mu_beta, Sigma_beta, mu_nuis, Sigma_nuis, mu_r, Sigma_r, psi_r_nuis)
  
  hess_w <- 0.5 * psi_taus * psi_omegas * hij
  sum_ww <- t(w_mat) %*% (hess_w * w_mat)
  prior_mat <- diag(1 / tau_prior, ncol(w_mat))
  Sigma_tau_new <- solve(sum_ww + prior_mat)
  
  grad_w <- 0.5 * (psi_taus * psi_omegas * hij - 1)
  lik_grad <- t(w_mat) %*% grad_w
  prior_grad <- -(1 / tau_prior) * mu_tau
  tol_grad <- lik_grad + prior_grad
  mu_tau_new <- mu_tau + Sigma_tau_new %*% tol_grad
  
  list(mu_tau_new = mu_tau_new, Sigma_tau_new = Sigma_tau_new)
}

omega_update <- function(X_mat, y_vec, id, uq_id, mu_beta, Sigma_beta, mu_nu, Sigma_nu, mu_r, Sigma_r, mu_inv_sigma, 
                         psi_taus, psi_omegas, psi_r_nus, mu_omegas){
  
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  psi_r_nuis <- psi_r_nus[match(id, uq_id)]
  
  hij <- cal_hij(y_vec, X_mat, id, uq_id, mu_beta, Sigma_beta, mu_nuis, Sigma_nuis, mu_r, Sigma_r, psi_r_nuis)
  
  tol_grad <- -0.5 + 0.5 * psi_omegas * psi_taus * hij - mu_inv_sigma * mu_omegas
  tol_hess <- -0.5 * psi_omegas * psi_taus * hij - mu_inv_sigma
  
  Sigma_omega_new <- -1 / tol_hess
  mu_omega_new <- mu_omegas + Sigma_omega_new * tol_grad
  
  list(mu_omega_new = as.vector(mu_omega_new), Sigma_omega_new = as.vector(Sigma_omega_new))
}

sigma_inv_update <- function(nobs, mu_omega, Sigma_omega, mu_inv_a, A_omega){
  sigma_shape <- (nobs+1)/2
  sum_omega_sq <- sum(Sigma_omega + mu_omega^2)
  rate_sigma <- mu_inv_a + 0.5 * sum_omega_sq
  
  mu_inv_sigma_new <- sigma_shape / rate_sigma
  mu_log_sigma <- log(rate_sigma) - digamma(sigma_shape)
  
  list(mu_inv_sigma_new = mu_inv_sigma_new, mu_log_sigma = mu_log_sigma)
}

a_inv_update <- function(mu_inv_sigma, A_omega){
  rate_a <- (1 / A_omega^2) + mu_inv_sigma
  mu_inv_a_new <- 1 / rate_a
  mu_log_a <- log(rate_a) - digamma(1)
  
  list(mu_inv_a_new = mu_inv_a_new, mu_log_a = mu_log_a)
}

r_update <- function(X_mat, y_vec, id, uq_id, mu_beta, Sigma_beta, mu_nu, Sigma_nu, mu_r, Sigma_r, 
                     psi_taus, psi_omegas, psi_r_nus, r_prior){
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  psi_r_nuis <- psi_r_nus[match(id, uq_id)]

  Bij <- cal_Bij(y_vec, X_mat, mu_beta, Sigma_beta, mu_nuis, Sigma_nuis, mu_r, Sigma_r)
  mul_B <- (-mu_nuis + mu_r * Sigma_nuis)
  DBij <- 2 * Sigma_nuis * (y_vec - X_mat %*% mu_beta) - 2 * Sigma_nuis * (mu_nuis - mu_r * Sigma_nuis)
  Eij <- mul_B * Bij + DBij
  DEij <- Sigma_nuis * Bij + mul_B * DBij + 2 * Sigma_nuis^2
  
  grad_prior <- -(1 / r_prior) * mu_r
  grad_var <- -0.5 * sum(mu_nuis)
  grad_lik <- -0.5 * sum(psi_taus * psi_omegas * psi_r_nuis * Eij)
  tol_grad <- grad_prior + grad_var + grad_lik
  
  hess_prior <- -1 / r_prior
  hess_in <- mul_B * Eij + DEij
  hess_lik <- -0.5 * sum(psi_taus * psi_omegas * psi_r_nuis * hess_in)
  tol_hess <- hess_prior + hess_lik
  
  Sigma_r_new <- -1 / tol_hess
  mu_r_new <- mu_r + Sigma_r_new * tol_grad
  
  list(mu_r_new = mu_r_new, Sigma_r_new = Sigma_r_new)
}

elbo_cal <- function(X_mat, y_vec, u_mat, w_mat, mu_beta, Sigma_beta, mu_nu, Sigma_nu, 
                     mu_alpha, Sigma_alpha, mu_tau, Sigma_tau, mu_omega, Sigma_omega, 
                     mu_inv_sigma, mu_log_sigma, mu_inv_a, mu_log_a, mu_r, Sigma_r, 
                     psi_taus, psi_alphas, psi_omegas, psi_r_nus, 
                     beta_prior, alpha_prior, tau_prior, A_prior, r_prior, 
                     id, uq_id){
  nobs <- length(y_vec)
  
  mu_nuis <- mu_nu[match(id, uq_id)]
  Sigma_nuis <- Sigma_nu[match(id, uq_id)]
  psi_r_nuis <- psi_r_nus[match(id, uq_id)]
  
  hij <- cal_hij(y_vec, X_mat, id, uq_id, mu_beta, Sigma_beta, mu_nuis, Sigma_nuis, mu_r, Sigma_r, psi_r_nuis)
  log_p_y <- -0.5 * sum(w_mat %*% mu_tau + mu_r * mu_nuis + mu_omega) - 0.5 * sum(psi_taus * psi_omegas * hij)
  
  log_p_beta <- -0.5 * (1 / beta_prior) * (sum(mu_beta^2) + sum(diag(Sigma_beta)))
  
  log_p_tau <- -0.5 * (1 / tau_prior) * (sum(mu_tau^2) + sum(diag(Sigma_tau)))
  
  log_p_alpha <- -0.5 * (1 / alpha_prior) * (sum(mu_alpha^2) + sum(diag(Sigma_alpha)))
  
  log_p_r <- -0.5 * (1 / r_prior) * (mu_r^2 + Sigma_r)
  
  log_p_nu <- -0.5 * sum(u_mat %*% mu_alpha + psi_alphas * (mu_nu^2 + Sigma_nu))
  
  log_p_omega <- -0.5 * mu_inv_sigma * sum(mu_omega^2 + Sigma_omega)
  
  log_p_sigma <- -0.5 * mu_log_a - 1.5 * mu_log_sigma - mu_inv_a * mu_inv_sigma
  
  log_p_a <- -1.5 * mu_log_a - (mu_inv_a / A_prior^2)
  
  log_q_beta <- 0.5 * log(det(Sigma_beta))
  
  log_q_tau <- 0.5 * log(det(Sigma_tau))
  
  log_q_alpha <- 0.5 * log(det(Sigma_alpha))
  
  log_q_r <- 0.5 * log(Sigma_r)
  
  log_q_nu <- 0.5 * sum(log(Sigma_nu))
  
  log_q_omega <- 0.5 * sum(log(Sigma_omega))
  
  rate_sigma <- mu_inv_a + 0.5 * sum(Sigma_omega + mu_omega^2)
  log_q_sigma <- log(rate_sigma)
  
  rate_a <- (1 / A_prior^2) + mu_inv_sigma
  log_q_a <- log(rate_a)
  
  elbo <- log_p_y + log_p_beta + log_p_tau + log_p_alpha + log_p_r + log_p_nu + log_p_omega + log_p_sigma + log_p_a + 
    log_q_beta + log_q_alpha + log_q_tau + log_q_r + log_q_nu + log_q_omega + log_q_sigma + log_q_a
  
  elbo
}

beta_update(X_matrix, y_vector, riesby$id, unique_ids, mu_nu_i_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu, beta_prior)

nu_update(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, Sigma_r_q,
          psi_ij_tau, psi_ij_omega, psi_i_alpha, psi_i_r_nu)

alpha_update(uis = u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)

tau_update(X_matrix, y_vector, w_matrix, riesby$id, unique_ids, mu_nu_i_q, Sigma_nu_i_q, 
           mu_beta_q, Sigma_beta_q, mu_r_q, Sigma_r_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu,
           mu_tau_q, tau_prior)

omega_update(X_matrix, y_vector, riesby$id, unique_ids, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, Sigma_r_q,
             mu_inv_sigma_omega_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu, mu_omega_q)

sigma_inv_update(length(y_vector), mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q, A_prior)

a_inv_update(mu_inv_sigma_omega_q, A_prior)

r_update(X_matrix, y_vector, riesby$id, unique_ids, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, Sigma_r_q, 
         psi_ij_tau, psi_ij_omega, psi_i_r_nu, r_prior)

elbo_cal(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
         mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q, 
         mu_inv_sigma_omega_q, mu_log_sigma_omega_q, 
         mu_inv_a_omega_q, mu_log_a_omega_q, mu_r_q, Sigma_r_q, 
         psi_ij_tau, psi_i_alpha, psi_ij_omega, psi_i_r_nu, 
         beta_prior, alpha_prior, tau_prior, A_prior, r_prior, 
         riesby$id, unique_ids)

###########################

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

# Initialize the priors
beta_prior <- 100000
alpha_prior <- 100000
tau_prior <- 100000
A_prior <- 100000
r_prior <- 100000

library(lme4)
lme_fit <- lmer(hamd ~ week + endog + endweek + (1 | id), data = riesby)
beta_init <- fixef(lme_fit)
nu_init <- as.vector(ranef(lme_fit)$id[, "(Intercept)"])
resid_init <- residuals(lme_fit) |> as.vector()

# Initialize beta variational parameters
p_beta <- ncol(X_matrix)
mu_beta_q <- beta_init
Sigma_beta_q <- diag(0.1, p_beta)

# Initialize nu_i's variational parameters
mu_nu_i_q <- nu_init
Sigma_nu_i_q <- rep(0.1, num_ids)

# Initialize alpha variational parameters
p_alpha <- ncol(u_matrix)
log_var_nu <- log(nu_init^2 + 1e-6)
lm_alpha <- lm(log_var_nu ~ u_matrix - 1)
mu_alpha_q <- matrix(coef(lm_alpha), nrow = p_alpha, ncol = 1)
Sigma_alpha_q <- diag(0.1, p_alpha)

# Initialize tau variational parameters
p_tau <- ncol(w_matrix)
log_var_resid <- log(resid_init^2 + 1e-6)
lm_tau <- lm(log_var_resid ~ w_matrix - 1)
mu_tau_q <- matrix(coef(lm_tau), nrow = p_tau, ncol = 1)
Sigma_tau_q <- diag(0.1, p_tau)

# Initialize omega_ij's variational parameters
nobs <- nrow(X_matrix)
mu_omega_q <- residuals(lm_tau)
Sigma_omega_q <- rep(0.1, nobs)

# Initialize sigma_omega^2's variational parameters
mu_log_sigma_omega_q <- 0
mu_inv_sigma_omega_q <- 1

# Initialize a_omega's variational parameters
mu_log_a_omega_q <- 0
mu_inv_a_omega_q <- 1

# Initialize r's variational parameters
mu_r_q <- 0
Sigma_r_q <- 0.01

# Initialize psi's
psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
psi_ij_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
psi_i_r_nu <- psi_r_nu(mu_r_q, Sigma_r_q, mu_nu_i_q, Sigma_nu_i_q)

max_iter <- 1000
tolerance <- 1e-6
elbo_history <- c()
#st_time <- Sys.time()

for (i in 1:max_iter){
  beta_new <- beta_update_ncr(X_matrix, y_vector, ids, unique_ids, mu_nu_i_q, psi_ij_tau, psi_ij_omega, beta_prior)
  mu_beta_q <- beta_new$mu_beta_new
  Sigma_beta_q <- beta_new$Sigma_beta_new
  
  nu_new <- nu_update_ncr(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, 
                          psi_ij_tau, psi_ij_omega, psi_i_alpha, mu_nu_i_q)
  mu_nu_i_q <- nu_new$mu_nu_new
  Sigma_nu_i_q <- nu_new$Sigma_nu_new
  
  alpha_new <- alpha_update(uis = u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)
  mu_alpha_q <- alpha_new$mu_alpha_new
  Sigma_alpha_q <- alpha_new$Sigma_alpha_new
  psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  
  tau_new <- tau_update_ncr(X_matrix, y_vector, w_matrix, ids, unique_ids, mu_nu_i_q, Sigma_nu_i_q, 
                            mu_beta_q, Sigma_beta_q, psi_ij_tau, psi_ij_omega, mu_tau_q, tau_prior)
  mu_tau_q <- tau_new$mu_tau_new
  Sigma_tau_q <- tau_new$Sigma_tau_new
  psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  
  omega_new <- omega_update_ncr(X_matrix, y_vector, ids, unique_ids, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                                mu_inv_sigma_omega_q, psi_ij_tau, psi_ij_omega, mu_omega_q)
  mu_omega_q <- omega_new$mu_omega_new
  Sigma_omega_q <- omega_new$Sigma_omega_new
  psi_ij_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  sigma_omega_new <- sigma_inv_update(length(y_vector), mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q)
  mu_inv_sigma_omega_q <- sigma_omega_new$mu_inv_sigma_new
  mu_log_sigma_omega_q <- sigma_omega_new$mu_log_sigma
  
  a_omega_new <- a_inv_update(mu_inv_sigma_omega_q, A_prior)
  mu_inv_a_omega_q <- a_omega_new$mu_inv_a_new
  mu_log_a_omega_q <- a_omega_new$mu_log_a
  
  elbo <- elbo_cal_ncr(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                       mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q, 
                       mu_inv_sigma_omega_q, mu_log_sigma_omega_q, 
                       mu_inv_a_omega_q, mu_log_a_omega_q, 
                       psi_ij_tau, psi_i_alpha, psi_ij_omega, beta_prior, alpha_prior, tau_prior, A_prior,
                       ids, unique_ids)
  elbo_history <- c(elbo_history, elbo)
  if ((i > 1) && (abs(elbo_history[i] - elbo_history[i-1]) < tolerance)){
    cat("CAVI converges at iteration ", i)
    break
  }
}

for (i in 1:max_iter){
  beta_new <- beta_update(X_matrix, y_vector, ids, unique_ids, mu_nu_i_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu, beta_prior)
  mu_beta_q <- beta_new$mu_beta_new
  Sigma_beta_q <- beta_new$Sigma_beta_new
  #print(beta_new)
  
  nu_new <- nu_update(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, 
                      Sigma_r_q, psi_ij_tau, psi_ij_omega, psi_i_alpha, psi_i_r_nu)
  mu_nu_i_q <- nu_new$mu_nu_new
  Sigma_nu_i_q <- nu_new$Sigma_nu_new
  psi_i_r_nu <- psi_r_nu(mu_r_q, Sigma_r_q, mu_nu_i_q, Sigma_nu_i_q)
  #print(nu_new)
  
  alpha_new <- alpha_update(uis = u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)
  mu_alpha_q <- alpha_new$mu_alpha_new
  Sigma_alpha_q <- alpha_new$Sigma_alpha_new
  psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  #print(alpha_new)
  
  tau_new <- tau_update(X_matrix, y_vector, w_matrix, ids, unique_ids, mu_nu_i_q, Sigma_nu_i_q, 
                        mu_beta_q, Sigma_beta_q, mu_r_q, Sigma_r_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu,
                        mu_tau_q, tau_prior)
  mu_tau_q <- tau_new$mu_tau_new
  Sigma_tau_q <- tau_new$Sigma_tau_new
  psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  #print(tau_new)
  
  omega_new <- omega_update(X_matrix, y_vector, ids, unique_ids, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, Sigma_r_q,
                            mu_inv_sigma_omega_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu, mu_omega_q)
  mu_omega_q <- omega_new$mu_omega_new
  Sigma_omega_q <- omega_new$Sigma_omega_new
  psi_ij_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  sigma_omega_new <- sigma_inv_update(length(y_vector), mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q)
  mu_inv_sigma_omega_q <- sigma_omega_new$mu_inv_sigma_new
  mu_log_sigma_omega_q <- sigma_omega_new$mu_log_sigma
  
  a_omega_new <- a_inv_update(mu_inv_sigma_omega_q, A_prior)
  mu_inv_a_omega_q <- a_omega_new$mu_inv_a_new
  mu_log_a_omega_q <- a_omega_new$mu_log_a
  
  r_new <- r_update(X_matrix, y_vector, ids, unique_ids, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, Sigma_r_q, 
                    psi_ij_tau, psi_ij_omega, psi_i_r_nu, r_prior)
  mu_r_q <- r_new$mu_r_new
  Sigma_r_q <- r_new$Sigma_r_new
  psi_i_r_nu <- psi_r_nu(mu_r_q, Sigma_r_q, mu_nu_i_q, Sigma_nu_i_q)
  
  elbo <- elbo_cal(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                   mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q, 
                   mu_inv_sigma_omega_q, mu_log_sigma_omega_q, 
                   mu_inv_a_omega_q, mu_log_a_omega_q, mu_r_q, Sigma_r_q, 
                   psi_ij_tau, psi_i_alpha, psi_ij_omega, psi_i_r_nu, 
                   beta_prior, alpha_prior, tau_prior, A_prior, r_prior, 
                   ids, unique_ids)
  elbo_history <- c(elbo_history, elbo)
  if ((i > 1) && (abs(elbo_history[i] - elbo_history[i-1]) < tolerance)){
    cat("CAVI converges at iteration ", i)
    break
  }
}
ed_time <- Sys.time()

print(beta_new)
print(alpha_new)
print(tau_new)
print(r_new)
sigma_omega_rate <- mu_inv_a_omega_q + 0.5 * sum(Sigma_omega_q + mu_omega_q^2)
sigma_omega_mean <- sigma_omega_rate / ((length(y_vector) + 1)/2 - 1)
print(sigma_omega_mean |> sqrt())
print((ed_time - st_time) |> round(4))

############################################

# Harvard data set
st_time <- Sys.time()
indat=read.table("Dataset_HealthBehavAcadPerfAffect.dat", header = FALSE, 
                 col.names=c("id", "day", "sex", "age", "sem", "sq", "physact", "pa", "na", "lga", 
                             "exam", "hsg", "bdi", "day_c"), na.strings="-99")
summary(indat)
indat <- subset(indat, !is.na(indat$pa))

X_matrix <- model.matrix(~ day_c, data = indat)
y_vector <- as.matrix(indat$pa)
w_matrix <- model.matrix(~ day_c, data = indat)

indat_id <- indat[!duplicated(indat$id), ]
ids <- indat$id
u_matrix <- model.matrix(~ day_c, data = indat_id)
unique_ids <- indat_id$id
num_ids <- length(unique_ids)

id_indices <- split(1:nrow(indat), indat$id)

# Initialize the priors
beta_prior <- 100000
alpha_prior <- 100000
tau_prior <- 100000
A_prior <- 100000
r_prior <- 100000

library(lme4)
lme_fit <- lmer(pa ~ day_c + (1 | id), data = indat)
beta_init <- fixef(lme_fit)
nu_init <- as.vector(ranef(lme_fit)$id[, "(Intercept)"])
resid_init <- residuals(lme_fit) |> as.vector()

# Initialize beta variational parameters
p_beta <- ncol(X_matrix)
mu_beta_q <- beta_init
Sigma_beta_q <- diag(0.1, p_beta)

# Initialize nu_i's variational parameters
mu_nu_i_q <- nu_init
Sigma_nu_i_q <- rep(0.1, num_ids)

# Initialize alpha variational parameters
p_alpha <- ncol(u_matrix)
log_var_nu <- log(nu_init^2 + 1e-6)
lm_alpha <- lm(log_var_nu ~ u_matrix - 1)
mu_alpha_q <- matrix(coef(lm_alpha), nrow = p_alpha, ncol = 1)
Sigma_alpha_q <- diag(0.1, p_alpha)

# Initialize tau variational parameters
p_tau <- ncol(w_matrix)
log_var_resid <- log(resid_init^2 + 1e-6)
lm_tau <- lm(log_var_resid ~ w_matrix - 1)
mu_tau_q <- matrix(coef(lm_tau), nrow = p_tau, ncol = 1)
Sigma_tau_q <- diag(0.1, p_tau)

# Initialize omega_ij's variational parameters
nobs <- nrow(X_matrix)
mu_omega_q <- residuals(lm_tau)
Sigma_omega_q <- rep(0.1, nobs)

# Initialize sigma_omega^2's variational parameters
mu_log_sigma_omega_q <- 0
mu_inv_sigma_omega_q <- 1

# Initialize a_omega's variational parameters
mu_log_a_omega_q <- 0
mu_inv_a_omega_q <- 1

# Initialize r's variational parameters
mu_r_q <- 0
Sigma_r_q <- 0.1

# Initialize psi's
psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
psi_ij_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
psi_i_r_nu <- psi_r_nu(mu_r_q, Sigma_r_q, mu_nu_i_q, Sigma_nu_i_q)

max_iter <- 100000
tolerance <- 1e-10
elbo_history <- c()

for (i in 1:max_iter){
  beta_new <- beta_update_ncr(X_matrix, y_vector, ids, unique_ids, mu_nu_i_q, psi_ij_tau, psi_ij_omega, beta_prior)
  mu_beta_q <- beta_new$mu_beta_new
  Sigma_beta_q <- beta_new$Sigma_beta_new
  
  nu_new <- nu_update_ncr(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, 
                          psi_ij_tau, psi_ij_omega, psi_i_alpha, mu_nu_i_q)
  mu_nu_i_q <- nu_new$mu_nu_new
  Sigma_nu_i_q <- nu_new$Sigma_nu_new
  
  alpha_new <- alpha_update(uis = u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)
  mu_alpha_q <- alpha_new$mu_alpha_new
  Sigma_alpha_q <- alpha_new$Sigma_alpha_new
  psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  
  tau_new <- tau_update_ncr(X_matrix, y_vector, w_matrix, ids, unique_ids, mu_nu_i_q, Sigma_nu_i_q, 
                            mu_beta_q, Sigma_beta_q, psi_ij_tau, psi_ij_omega, mu_tau_q, tau_prior)
  mu_tau_q <- tau_new$mu_tau_new
  Sigma_tau_q <- tau_new$Sigma_tau_new
  psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  
  omega_new <- omega_update_ncr(X_matrix, y_vector, ids, unique_ids, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                                mu_inv_sigma_omega_q, psi_ij_tau, psi_ij_omega, mu_omega_q)
  mu_omega_q <- omega_new$mu_omega_new
  Sigma_omega_q <- omega_new$Sigma_omega_new
  psi_ij_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  sigma_omega_new <- sigma_inv_update(length(y_vector), mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q)
  mu_inv_sigma_omega_q <- sigma_omega_new$mu_inv_sigma_new
  mu_log_sigma_omega_q <- sigma_omega_new$mu_log_sigma
  
  a_omega_new <- a_inv_update(mu_inv_sigma_omega_q, A_prior)
  mu_inv_a_omega_q <- a_omega_new$mu_inv_a_new
  mu_log_a_omega_q <- a_omega_new$mu_log_a
  
  elbo <- elbo_cal_ncr(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                       mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q, 
                       mu_inv_sigma_omega_q, mu_log_sigma_omega_q, 
                       mu_inv_a_omega_q, mu_log_a_omega_q, 
                       psi_ij_tau, psi_i_alpha, psi_ij_omega, beta_prior, alpha_prior, tau_prior, A_prior,
                       ids, unique_ids)
  elbo_history <- c(elbo_history, elbo)
  if ((i > 1) && (abs(elbo_history[i] - elbo_history[i-1]) < tolerance)){
    cat("CAVI converges at iteration ", i)
    break
  }
}

for (i in 1:max_iter){
  beta_new <- beta_update(X_matrix, y_vector, ids, unique_ids, mu_nu_i_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu, beta_prior)
  mu_beta_q <- beta_new$mu_beta_new
  Sigma_beta_q <- beta_new$Sigma_beta_new
  #print(beta_new)
  
  nu_new <- nu_update(X_matrix, y_vector, unique_ids, id_indices, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, 
                      Sigma_r_q, psi_ij_tau, psi_ij_omega, psi_i_alpha, psi_i_r_nu)
  mu_nu_i_q <- nu_new$mu_nu_new
  Sigma_nu_i_q <- nu_new$Sigma_nu_new
  psi_i_r_nu <- psi_r_nu(mu_r_q, Sigma_r_q, mu_nu_i_q, Sigma_nu_i_q)
  #print(nu_new)
  
  alpha_new <- alpha_update(uis = u_matrix, psi_i_alpha, mu_nu_i_q, Sigma_nu_i_q, mu_alpha_q, alpha_prior)
  mu_alpha_q <- alpha_new$mu_alpha_new
  Sigma_alpha_q <- alpha_new$Sigma_alpha_new
  psi_i_alpha <- psi_alpha(u_matrix, mu_alpha_q, Sigma_alpha_q)
  #print(alpha_new)
  
  tau_new <- tau_update(X_matrix, y_vector, w_matrix, ids, unique_ids, mu_nu_i_q, Sigma_nu_i_q, 
                        mu_beta_q, Sigma_beta_q, mu_r_q, Sigma_r_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu,
                        mu_tau_q, tau_prior)
  mu_tau_q <- tau_new$mu_tau_new
  Sigma_tau_q <- tau_new$Sigma_tau_new
  psi_ij_tau <- psi_tau(w_matrix, mu_tau_q, Sigma_tau_q)
  #print(tau_new)
  
  omega_new <- omega_update(X_matrix, y_vector, ids, unique_ids, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, Sigma_r_q,
                            mu_inv_sigma_omega_q, psi_ij_tau, psi_ij_omega, psi_i_r_nu, mu_omega_q)
  mu_omega_q <- omega_new$mu_omega_new
  Sigma_omega_q <- omega_new$Sigma_omega_new
  psi_ij_omega <- psi_omega(mu_omega_q, Sigma_omega_q)
  
  sigma_omega_new <- sigma_inv_update(length(y_vector), mu_omega_q, Sigma_omega_q, mu_inv_a_omega_q)
  mu_inv_sigma_omega_q <- sigma_omega_new$mu_inv_sigma_new
  mu_log_sigma_omega_q <- sigma_omega_new$mu_log_sigma
  
  a_omega_new <- a_inv_update(mu_inv_sigma_omega_q, A_prior)
  mu_inv_a_omega_q <- a_omega_new$mu_inv_a_new
  mu_log_a_omega_q <- a_omega_new$mu_log_a
  
  r_new <- r_update(X_matrix, y_vector, ids, unique_ids, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, mu_r_q, Sigma_r_q, 
                    psi_ij_tau, psi_ij_omega, psi_i_r_nu, r_prior)
  mu_r_q <- r_new$mu_r_new
  Sigma_r_q <- r_new$Sigma_r_new
  # print(r_new)
  psi_i_r_nu <- psi_r_nu(mu_r_q, Sigma_r_q, mu_nu_i_q, Sigma_nu_i_q)
  
  elbo <- elbo_cal(X_matrix, y_vector, u_matrix, w_matrix, mu_beta_q, Sigma_beta_q, mu_nu_i_q, Sigma_nu_i_q, 
                   mu_alpha_q, Sigma_alpha_q, mu_tau_q, Sigma_tau_q, mu_omega_q, Sigma_omega_q, 
                   mu_inv_sigma_omega_q, mu_log_sigma_omega_q, 
                   mu_inv_a_omega_q, mu_log_a_omega_q, mu_r_q, Sigma_r_q, 
                   psi_ij_tau, psi_i_alpha, psi_ij_omega, psi_i_r_nu, 
                   beta_prior, alpha_prior, tau_prior, A_prior, r_prior, 
                   ids, unique_ids)
  elbo_history <- c(elbo_history, elbo)
  if ((i > 1) && (abs(elbo_history[i] - elbo_history[i-1]) < tolerance)){
    cat("CAVI converges at iteration ", i)
    break
  }
}
ed_time <- Sys.time()

print(beta_new)
print(alpha_new)
print(tau_new)
sigma_omega_rate <- mu_inv_a_omega_q + 0.5 * sum(Sigma_omega_q + mu_omega_q^2)
sigma_omega_mean <- sigma_omega_rate / ((length(y_vector) + 1)/2 - 1)
print(sigma_omega_mean)
print(r_new)
print((ed_time - st_time) |> round(4))
