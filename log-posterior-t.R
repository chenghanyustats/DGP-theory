################################################################################
# Log posterior of t in the theory paper (given hyperparameters)
# with Beta prior
################################################################################
log_post_t <- function(t, y, x, Kff, A, sig2, lambda, h,
                            shape1 = 1, shape2 = 1, a = 0, b = 2) {

    Kdf <- computeCovDer1(t, x, ker_der1_fcn = se_ker_der1_no_tau, h = h)
    # B <- get_B(Kdf = Kdf, h = h)    
    a_mat <- get_a_mat(Kdf = Kdf, h = h)
    at_Ainv <- t(a_mat) %*% chol2inv(A)
    mu_t <- get_mu_t(y, A, a_mat, h, at_Ainv = at_Ainv)
    sig2_t <-  get_sig2_t(sig2, lambda, A, a_mat, h, at_Ainv = at_Ainv)
    
    (-1 / 2) * (log(sig2_t) + mu_t ^ 2 / sig2_t) + 
        log_dgbeta(t, shape1 = shape1, shape2 = shape2, a = a, b = b)
}



log_post_t_new <- function(t, y, x, Kff, A, Ainv, sig2, lambda, h, tau2, 
                       shape1 = 1, shape2 = 1, a = 0, b = 2) {
    n <- length(y)
    Kdf <- computeCovDer1(t, x, ker_der1_fcn = se_ker_der1_no_tau, h = h)
    # B <- get_B(Kdf = Kdf, h = h)    
    # a_mat <- get_a_mat(Kdf = Kdf, h = h)
    # at_Ainv <- t(a_mat) %*% chol2inv(A)
    mu_t <- get_mu_t_new(y, Kdf = Kdf, Ainv = Ainv)
    sig2_t <-  get_sig2_t_new(sig2, lambda, tau2, h, Kdf = Kdf, Ainv = Ainv)
    
    log_t <- (-1 / 2) * (log(sig2_t) + mu_t ^ 2 / sig2_t) + 
        log_dgbeta(t, shape1 = shape1, shape2 = shape2, a = a, b = b)
    
    log_C <- (-n / 2) * log(2 * pi * tau2) - 
        (1 / 2) * determinant(A, logarithm = TRUE)$modulus - 
        quad.form(Ainv, y) - (2 * tau2) + (1 / 2) * log(tau2) - log(h)
    
    return(log_t + log_C)
}

# k <- 55
# t = grid_t[i]
# y = YY[[k]]$y
# x = YY[[k]]$x
# Kff = Kff_gp
# A = A_gp
# lambda = lambda_gp
# h = h_gp
# sig2 = sig_gp ^ 2
# shape1 = 1
# shape2 = 1
# a = x_a 
# b = x_b


log_post_t_theory <- function(t, y, x, Kff, A, lambda, h, sig2,
                              shape1 = 1, shape2 = 1, a = 0, b = 2) {
    n <- length(x)
    tau2 <- sig2 / (n * lambda)
    Kdf <- computeCovDer1(idx1 = t, idx2 = x, 
                          tau = 1, h = h) 
    
    a_mat <- t(Kdf) * h
    mu_t <- t(a_mat) %*% solve(A, y) / h
    dd <- 1 - quad.form.inv(A, a_mat)
    sig2_t <- tau2 * dd / (h ^ 2)
    cholA <- chol(A)
    log_det_A <- 2 * sum(log(diag(cholA)))
    log_C <- (-1 / 2) * (n * log(2 * pi * tau2) + log_det_A +
                             quad.form.inv(A, y) / tau2 - log(tau2) + 
                             2 * log(h))
    # log_C <- (-n / 2) * log(2 * pi * tau2) -
    #     (1 / 2) * log(det(A)) - quad.form.inv(A, y) / (2 * tau2) + 
    #     (1 / 2) * log(tau2) - log(h)
    log_den_t <- (-1 / 2) * (log(sig2_t) + crossprod(mu_t) / sig2_t) 
    log_prior <- log_dgbeta(t, shape1, shape2, a, b)
    
    return(as.vector(log_C + log_den_t + log_prior))
}


log_post_t_theory_matern <- function(t, y, x, Kff, A, lambda, l, sig2, nu = 1.5,
                                     shape1 = 1, shape2 = 1, a = 0, b = 2) {
    n <- length(x)
    tau2 <- sig2 / (n * lambda)
    Kdf <- computeCovDer1(idx1 = t, idx2 = x, ker_der1_fcn = matern_ker_der1,
                          tau = 1, l = l, nu = nu) 
    if (nu == 1.5) {
        Kdd <- 3 / (l ^ 2)
    } else {
        Kdd <- 5 / (3 * l ^ 2)
    }
    
    a_mat <- t(Kdf) * (1 / sqrt(Kdd))
    mu_t <- t(a_mat) %*% solve(A, y) * sqrt(Kdd)
    dd <- 1 - quad.form.inv(A, a_mat)
    sig2_t <- tau2 * dd * Kdd
    cholA <- chol(A)
    log_det_A <- 2 * sum(log(diag(cholA)))

    # log_C <- (-n / 2) * log(2 * pi * tau2) -
    #     (1 / 2) * log(det(A)) - quad.form.inv(A, y) / (2 * tau2) +
    #     (1 / 2) * (log(tau2) + log(Kdd))
    
    log_C <- (-1 / 2) * (n * log(2 * pi * tau2) + log_det_A + 
                             quad.form.inv(A, y) / tau2 - 
                             (log(tau2) + log(Kdd)))
    
    log_den_t <- (-1 / 2) * (log(sig2_t) + crossprod(mu_t) / sig2_t) 
    log_prior <- log_dgbeta(t, shape1, shape2, a, b)
    
    return(log_C + log_den_t + log_prior)
}


get_mu_t_new <- function(y, Kdf, Ainv) {
    Kdf %*% Ainv %*% y    
}

get_sig2_t_new <- function(sig2, lambda, tau2 = NULL, h, Kdf, Ainv) {
    if(!is.null(tau2)) {
        return(tau2 * (1 / h ^ 2 - quad.tform(Ainv, Kdf)))
    } else {
        sig2 / (nrow(A) * lambda) * (1 / h ^ 2 - quad.tform(Ainv, Kdf))
    }
}


log_post_diff_t <- function(t_star, t, y, x, Kff, A, sig2, lambda, h,
                            shape1 = 1, shape2 = 1, a, b) {
    
    Kdf_star <- computeCovDer1(t_star, x, ker_der1_fcn = se_ker_der1_no_tau, h = h)
    Kdf <- computeCovDer1(t, x, ker_der1_fcn = se_ker_der1_no_tau, h = h)
    # B <- get_B(Kdf = Kdf, h = h)    
    a_mat_star <- get_a_mat(Kdf = Kdf_star, h = h)
    a_mat <- get_a_mat(Kdf = Kdf, h = h)
    
    A_inv <- chol2inv(A)
    at_Ainv_star <- t(a_mat_star) %*% A_inv
    at_Ainv <- t(a_mat) %*% A_inv
    
    mu_t_star <- get_mu_t(y, A, a_mat_star, h, at_Ainv = at_Ainv_star)
    mu_t <- get_mu_t(y, A, a_mat, h, at_Ainv = at_Ainv)
    
    sig2_t_star <-  get_sig2_t(sig2, lambda, A, a_mat_star, h, 
                               at_Ainv = at_Ainv_star)
    sig2_t <-  get_sig2_t(sig2, lambda, A, a_mat, h, at_Ainv = at_Ainv)
    
    (-1 / 2) * ((sig2_t_star / sig2_t) - 
        (mu_t_star ^ 2 / sig2_t_star -  mu_t ^ 2 / sig2_t)) + 
        log_dgbeta(t_star, shape1 = shape1, shape2 = shape2, a = a, b = b) -
        log_dgbeta(t, shape1 = shape1, shape2 = shape2, a = a, b = b)
}




log_mar_lik_gp <- function(theta, y, H0) {
    
    # K <- compute_cov_1d(idx1 = x_vec, tau = theta[2], h = theta[3])
    K <- se_ker(H0 = H0, tau = theta[2], h = theta[3])
    n <- length(y)
    return(-mvnfast::dmvn(y, rep(0L, n), K + diag(theta[1] ^ 2, n), log = TRUE))
}





