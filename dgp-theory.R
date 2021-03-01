# ##############################################################################
# Semiparametric Bayesian inference for local extrema of functions in the presence of noise
# Simulation Study
# Cheng-Han Yu
# ##############################################################################


## ============================================================================
## load data and packages
## ============================================================================
load("./data/sim_data_n100.RData", verbose = TRUE)

library(ftnonpar)
library(matrixcalc)
library(emulator)
library(R.utils)
library(parallel)
library(doParallel)

## ============================================================================
## Source functions
## ============================================================================
source("./song_fcns.R")
source("./log-posterior-t.R")
sourceDirectory("./Rcode/R/", modifiedOnly = FALSE)

## get hpd interval from density 
get_hpd_interval_from_den <- function(post_prob, grid_t,
                                      max_iter = 1000, tol = 0.01,
                                      step_size = 0.001,
                                      target_prob = 0.95) {
    sum_density <- sum(post_prob)
    max_density <- max(post_prob)
    current_density <- max_density
    prob_value <- 0
    count <- 1
    # sample_size <- length(samples)
    while(abs(prob_value - target_prob) > tol) {
        if (prob_value > target_prob) {
            current_density <- current_density + step_size
        } else {
            current_density <- current_density / sqrt(1 + count)
        }
        
        high_density_idx <- which(post_prob > current_density)
        clu_idx <- which(diff(high_density_idx) > 2)
        (prob_value <- sum(post_prob[high_density_idx]) / sum_density)
        
        ci_lower_b <- c(grid_t[high_density_idx][1], 
                        grid_t[high_density_idx][clu_idx + 1])
        ci_upper_b <- c(grid_t[high_density_idx][clu_idx],
                        grid_t[high_density_idx][length(high_density_idx)])
        no_clu <- length(ci_lower_b)
        # samples[samples > ci_lower & samples < ci_upper]
        count <- count + 1
        # print(prob_value)
        if (count == max_iter) break
    }
    return(list(no_cluster = no_clu, 
                ci_lower = ci_lower_b, ci_upper = ci_upper_b, 
                prob_value = prob_value,
                den_value = current_density))
}


## smooth taut string method
smqreg <- function (y, thr.const = 2.5, verbose = FALSE, bandwidth = -1, 
                    sigma = -1, localsqueezing = TRUE, squeezing.factor = 0.5, 
                    DYADIC = TRUE, firstlambda = 100, smqeps = 1/length(y), 
                    fsign = double(0), 
                    gensign = TRUE, tolerance = 1e-12, ...) {
    n <- length(y)
    lambda <- rep(firstlambda, n - 1)
    sigma <- mad((y[-1] - y[-n])/sqrt(2))
    if (verbose) 
        print(c("sigma is ", sigma))
    if (gensign) 
        fsign <- gensign(y, thr.const = thr.const, extrema.mean = TRUE, 
                         sigma = sigma, localsqueezing = localsqueezing, squeezing.factor = squeezing.factor)
    repeat {
        f <- smqnew(y = y, lambda = lambda, eps = smqeps, fsign = fsign, 
                    tolerance = tolerance)
        if (bandwidth < 0) {
            residuals <- y - f
            residuals <- residuals - mean(residuals)
            if (DYADIC) 
                residuals.wr <- multiwdwr(residuals, sqrt(thr.const * 
                                                              log(n)) * sigma)
            else residuals.wr <- nondymwdr(residuals, sqrt(thr.const * 
                                                               log(n)) * sigma)
        }
        if (verbose) {
            par(mfrow = c(2, 1))
            plot(y, col = "grey", ...)
            lines(f, col = "red")
            lines(residuals.wr, type = "l", col = "green")
            lines(fsign - 2, col = "blue")
            plot(lambda, ty = "b", ...)
            print("Press Enter")
            dum <- readline()
        }
        if (bandwidth > 0) 
            break
        ind <- (abs(residuals.wr) > 1e-10)
        ind2 <- ind[-1] | ind[-n]
        if (sum(ind) == 0) 
            break
        if (localsqueezing) 
            lambda[ind2] <- lambda[ind2] * squeezing.factor
        else lambda <- lambda * squeezing.factor
        if (min(lambda) < 1e-08) {
            par(mfrow = c(2, 1))
            plot(y, col = "grey", ...)
            lines(f, col = "red")
            lines(residuals.wr, type = "l", col = "green")
            lines(fsign - 2, col = "blue")
            plot(lambda, ty = "b", ...)
            print("ERROR, This should not happen")
            break
        }
    }
    list(y = f, nmax = findmod(f)$mod, loc = findmod(f)$x, sigma = sigma)
}

## get maximum a posteriori
get_map_new <- function(post_den, grid_t, hpdi_1_2_lst) {
    no_data <- nrow(post_den)
    grid_len <- length(grid_t)
    
    map_est <- vector("list", length = no_data)
    for (i in 1:no_data) {
        num_map <- nrow(hpdi_1_2_lst[[i]])
        for (k in 1:num_map) {
            lower <- which.min(abs(grid_t - hpdi_1_2_lst[[i]][k, 1]))
            upper <- which.min(abs(grid_t - hpdi_1_2_lst[[i]][k, 2]))
            
            map_est[[i]][[k]] <- grid_t[lower + which.max(post_den[i, lower:upper]) - 1]
            # print(map_est[[i]][[k]])
        }
        
        # map_est[i, 1] <- grid_t[which.max(post_den[i, 1:313])]
        # map_est[i, 2] <- grid_t[which.max(post_den[i, 314:1000]) + 313]
    }
    map_est
}



## ============================================================================
## Objects/variables to be used
## ============================================================================
H0_diff <- lapply(YY, function(d) {
    outer(as.vector(d$x), as.vector(d$x),
          FUN = function(x1, x2) (x1 - x2))
})

idx <- seq(x_a, x_b, length.out = 500)
x_test <- sort(seq(x_a, x_b, length.out = 100))
n_test <- length(x_test)
n <- length(YY[[1]]$x)
no_data <- 100
grid_t <- seq(x_a, x_b, length.out = 400)
type1err <- 0.05


cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)
EB_gp <- foreach(k = 1:no_data) %dopar% {
    res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = YY[[k]]$y, H0 = H0_diff[[k]])
    res$par
}
stopCluster(cl)

## ============================================================================
## DGP beta(1, 1) by Li et al. (2021)
## ============================================================================
log_post_gp_theory <- matrix(0, no_data, length(grid_t))
post_prob_t_gp_theory <- matrix(0, no_data, length(grid_t))
post_prob_gp_hist <- matrix(0, no_data, length(grid_t))

for (k in 1:no.data) { 
    cat("k =", k, "\r")
    sig_gp <- EB_gp[[k]][1]
    tau_gp <- EB_gp[[k]][2]
    h_gp <- EB_gp[[k]][3]
    lambda_gp <- sig_gp ^ 2 / (n * tau_gp ^ 2)
    Kff_gp <- se_ker(H0 = H0_diff[[k]], tau = 1, h = h_gp)
    A_gp <- Kff_gp + diag((n * lambda_gp), n)
    for (i in 1:length(grid_t)) {
        log_post_gp_theory[k, i] <- log_post_t_theory(t = grid_t[i],
                                                      y = YY[[k]]$y,
                                                      x = YY[[k]]$x,
                                                      Kff = Kff_gp,
                                                      A = A_gp,
                                                      lambda = lambda_gp,
                                                      h = h_gp,
                                                      sig2 = sig_gp ^ 2,
                                                      shape1 = 1, shape2 = 1,
                                                      a = x_a, b = x_b)
    }
}


for (k in 1:no_data) {
    post_prob_t_gp_theory[k, ] <- exp(log_post_gp_theory[k, ] -
                                              max(log_post_gp_theory[k, ]))
    post_prob_gp_hist[k, ] <- post_prob_t_gp_theory[k, ] * 2 * diff(grid_t)[1] /
        sum(post_prob_t_gp_theory[k, ] * 2 * diff(grid_t)[1])
}



hpdi_lst_n100 <- apply(post_prob_t_gp_theory, 1, get_hpd_interval_from_den, 
                       grid_t = grid_t, target_prob = 1 - type1err)
hpdi_1_lst <- lapply(hpdi_lst_n100, function(x) x$ci_lower)
hpdi_2_lst <- lapply(hpdi_lst_n100, function(x) x$ci_upper)

hpdi_1_2_lst <- list()
for (i in 1:no_data) {
    hpdi_1_2_lst[[i]] <- cbind(hpdi_1_lst[[i]], hpdi_2_lst[[i]])
}


## ============================================================================
## DGP beta(2, 3) by Li et al. (2021)
## ============================================================================
log_post_gp_theory_2_3 <- matrix(0, no_data, length(grid_t))
post_prob_t_gp_theory_2_3 <- matrix(0, no_data, length(grid_t))
post_prob_gp_hist_2_3 <- matrix(0, no_data, length(grid_t))

for (k in 1:no_data) { 
    cat("k =", k, "\r")
    sig_gp <- EB_gp[[k]][1]
    tau_gp <- EB_gp[[k]][2]
    h_gp <- EB_gp[[k]][3]
    lambda_gp <- sig_gp ^ 2 / (n * tau_gp ^ 2)
    Kff_gp <- se_ker(H0 = H0_diff[[k]], tau = 1, h = h_gp)
    A_gp <- Kff_gp + diag((n * lambda_gp), n)
    for (i in 1:length(grid_t)) {
        log_post_gp_theory_2_3[k, i] <- log_post_t_theory(t = grid_t[i],
                                                          y = YY[[k]]$y,
                                                          x = YY[[k]]$x,
                                                          Kff = Kff_gp,
                                                          A = A_gp,
                                                          lambda = lambda_gp,
                                                          h = h_gp,
                                                          sig2 = sig_gp ^ 2,
                                                          shape1 = 2, 
                                                          shape2 = 3,
                                                          a = x_a, b = x_b)
    }
}


for (k in 1:no_data) {
    post_prob_t_gp_theory_2_3[k, ] <- exp(log_post_gp_theory_2_3[k, ] -
                                              max(log_post_gp_theory_2_3[k, ]))
    post_prob_gp_hist_2_3[k, ] <- post_prob_t_gp_theory_2_3[k, ] * 2 * 
        diff(grid_t)[1] / sum(post_prob_t_gp_theory_2_3[k, ] * 
                                  2 * diff(grid_t)[1])
}

hpdi_lst_n100_beta_2_3 <- apply(post_prob_t_gp_theory_2_3, 1, 
                                get_hpd_interval_from_den, grid_t = grid_t)

hpdi_1_lst_beta_2_3_n100 <- lapply(hpdi_lst_n100_beta_2_3, 
                                   function(x) x$ci_lower)
hpdi_2_lst_beta_2_3_n100 <- lapply(hpdi_lst_n100_beta_2_3, 
                                   function(x) x$ci_upper)
hpdi_1_2_lst_beta_2_3_n100 <- list()
for (i in 1:no_data) {
    hpdi_1_2_lst_beta_2_3_n100[[i]] <- cbind(hpdi_1_lst_beta_2_3_n100[[i]], 
                                             hpdi_2_lst_beta_2_3_n100[[i]])
}


## ============================================================================
## STS by Davies, P. L. and Kovac, A. (2001) and Kovac (2006)
## ============================================================================
smqreg_fit_lst <- list()
kovac_est <- list()

for (i in 1:no_data) {
    smqreg_fit_lst[[i]] <- smqreg(YY[[i]]$y)
    kovac_est[[i]] <- YY[[i]]$x[smqreg_fit_lst[[i]]$loc+1]
}

## ============================================================================
## NKS by Song (2006)
## ============================================================================
bw_vec <- rep(0, no_data)
beta_lst <- list()
ci_song_lst_bon <- list()


for (i in 1:no_data) {
    bw_vec[i] <- dpill(x = YY[[i]]$x, y = YY[[i]]$y)
    beta_lst[[i]] <- LQfit(x = YY[[i]]$x, y = YY[[i]]$y, h = bw_vec[i])
    ci_song_lst_bon[[i]] <- AddCI(x = YY[[i]]$x, y = YY[[i]]$y, h = bw_vec[i], 
                                  beta = beta_lst[[i]], alpha = type1err / 3)
    all_info_lst[[i]] <- cbind(pos = YY[[i]]$x, beta_lst[[i]], 
                               ci_song_lst_bon[[i]])
}

song_der_zero_idx_lst <- find_song_der_zero_idx(all_info_lst, no_data, 
                                                is.print = TRUE)
song_der_zero_pt_lst <- song_der_zero_idx_lst

## Obtain CI (most annoying part! Should wrap up as functions)
ci_pts_song_lst <- list()
for (i in 1:no_data) {
    # print(paste("i=", i))
    stationary_pt <- all_info_lst[[i]][song_der_zero_pt_lst[[i]], ]
    if(length(song_der_zero_pt_lst[[i]]) == 1) {
        stationary_pt_upper <- stationary_pt[4]
        stationary_pt_lower <- stationary_pt[5]
    } else {
        stationary_pt_upper <- stationary_pt[, 4]
        stationary_pt_lower <- stationary_pt[, 5]
    }
    
    upper_vec <- c()
    lower_vec <- c()
    for(j in 1:length(song_der_zero_pt_lst[[i]])) {
        idx <- song_der_zero_pt_lst[[i]][j]
        if (all_info_lst[[i]][idx, 3] - 
            all_info_lst[[i]][idx-1, 3] < 0) {
            
            
            upper_idx_search <- which(diff(rev(all_info_lst[[i]][1:(idx-1), 3])) < 0)
            
            if (length(upper_idx_search) != 0) {
                upper_idx <- all_info_lst[[i]][(idx-1):(idx-1-upper_idx_search[1]), 3][all_info_lst[[i]][(idx-1):(idx-1-upper_idx_search[1]), 3] > 
                                                                                           stationary_pt_upper[j]][1]
            } else {
                upper_idx <- rev(all_info_lst[[i]][1:(idx-1), 3][all_info_lst[[i]][1:(idx-1), 3] > 
                                                                     stationary_pt_upper[j]])[1]
            }
            
            if (!is.na(upper_idx)) {
                lower_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == upper_idx), 1]
            }
            
            
            
            lower_idx_search <-  which(diff(all_info_lst[[i]][(idx+1):n, 3]) > 0)
            
            
            if(length(lower_idx_search) != 0) {
                lower_idx <- all_info_lst[[i]][(idx+1):(idx+lower_idx_search[1]), 3][all_info_lst[[i]][(idx+1):(idx+lower_idx_search[1]), 3] < 
                                                                                         stationary_pt_lower[j]][1]
            } else{
                lower_idx <- all_info_lst[[i]][(idx+1):n, 3][all_info_lst[[i]][(idx+1):n, 3] < 
                                                                 stationary_pt_lower[j]][1]
            }
            
            
            if (!is.na(lower_idx)) {
                upper_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == lower_idx), 1]
            }
            
            
            if(is.na(upper_idx) & !is.na(lower_idx)) {
                
                lower_pt <- 2 * stationary_pt[j, 1] - upper_pt
            }
            
            if(!is.na(upper_idx) & is.na(lower_idx)) {
                upper_pt <- 2 * stationary_pt[j, 1] - lower_pt
            }
            
            if(is.na(upper_idx) & is.na(lower_idx)) {
                lower_pt <- NA
                upper_pt <- NA
            }
            
            
        } else {
            upper_idx_search <- which(diff(all_info_lst[[i]][(idx+1):n, 3]) < 0)
            if (length(upper_idx_search) != 0) {
                upper_idx <- all_info_lst[[i]][(idx+1):(idx+upper_idx_search[1]), 3][all_info_lst[[i]][(idx+1):(idx+upper_idx_search[1]), 3] > 
                                                                                         stationary_pt_upper[j]][1]
            } else {
                upper_idx <- all_info_lst[[i]][(idx+1):n, 3][all_info_lst[[i]][(idx+1):n, 3] >
                                                                 stationary_pt_upper[j]][1]
            }
            
            if (!is.na(upper_idx)) {
                upper_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == upper_idx), 1]
            }
            
            lower_idx_search <- which(diff(rev(all_info_lst[[i]][1:(idx-1), 3])) > 0)
            if(length(lower_idx_search) != 0) {
                lower_idx <- all_info_lst[[i]][(idx-1):(idx-lower_idx_search[1]), 3][all_info_lst[[i]][(idx-1):(idx-lower_idx_search[1]), 3] < 
                                                                                         stationary_pt_lower[j]][1]
            } else{
                lower_idx <- rev(all_info_lst[[i]][1:(idx-1), 3][all_info_lst[[i]][1:(idx-1), 3] <
                                                                     stationary_pt_lower[j]])[1]
            }
            
            if (!is.na(lower_idx)) {
                lower_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == lower_idx), 1]
            }
            
            
            if(is.na(upper_idx) & !is.na(lower_idx)) {
                upper_pt <- 2 * stationary_pt[j, 1] - lower_pt
            }
            
            if(!is.na(upper_idx) & is.na(lower_idx)) {
                lower_pt <- 2 * stationary_pt[j, 1] - upper_pt
            }
            
            if(is.na(upper_idx) & is.na(lower_idx)) {
                lower_pt <- NA
                upper_pt <- NA
            }
        }
        upper_vec <- c(upper_vec, upper_pt)
        lower_vec <- c(lower_vec, lower_pt)
    }
    
    
    ci_pts_song_lst[[i]] <- cbind(lower_vec, upper_vec)
}

imp_lst <- list()
for(i in 1:no_data) {
    position <- song_der_zero_idx_lst[[i]]
    x_val <- YY[[i]]$x[position]
    lCI <- ci_pts_song_lst[[i]][, 1]
    uCI <- ci_pts_song_lst[[i]][, 2]
    imp_lst[[i]] <- list(position = position,
                         x = x_val,
                         lCI = lCI,
                         uCI = uCI)
}


## ============================================================================
## Number of local extrema
## ============================================================================
# -----------------
## DGP-beta(1, 1)
# -----------------
map_est_n100 <- get_map_new(post_prob_t_gp_theory, 
                            grid_t, hpdi_1_2_lst_n100)
table(unlist(lapply(map_est_n100, length)))

# -----------------
## DGP-beta(2, 3)
# -----------------
map_est_beta_2_3_n100 <- get_map_new(post_prob_t_gp_theory_2_3, grid_t, 
                                     hpdi_1_2_lst_beta_2_3_n100)

table(unlist(lapply(map_est_beta_2_3_n100, length)))

# -----------------
## STS
# -----------------
n_sp_sts <- sapply(smqreg_fit_lst, function(x) x$nmax)
table(sapply(smqreg_fit_lst, function(x) x$nmax))


# -----------------
## NKS
# -----------------
table(unlist(lapply(song_der_zero_idx_lst, function(x) length(x))))


## ============================================================================
## RMSE 
## ============================================================================

b0 <- 0
b3 <- 1
b1 <- (cri_pts[1] + cri_pts[2]) / 2
b2 <- (cri_pts[2] + cri_pts[3]) / 2

# -----------------
## DGP-beta(1, 1)
# -----------------

map_mat_m2_n100 <- matrix(0, nrow = no_data, 3)
for (i in 1:no_data) {
    x_1 <- map_est_n100[[i]][map_est_n100[[i]] > b0 & map_est_n100[[i]] < b1]
    map_mat_m2_n100[i, 1] <- mean(unlist(x_1))
    
    x_2 <- map_est_n100[[i]][map_est_n100[[i]] > b1 & map_est_n100[[i]] < b2]
    map_mat_m2_n100[i, 2] <- mean(unlist(x_2))
    
    x_3 <- map_est_n100[[i]][map_est_n100[[i]] > b2 & map_est_n100[[i]] < b3]
    map_mat_m2_n100[i, 3] <- mean(unlist(x_3))
}

sqrt(mean((map_mat_m2_n100[, 1] - cri_pts[1]) ^ 2)) * 100
sqrt(mean((map_mat_m2_n100[, 2] - cri_pts[2]) ^ 2)) * 100
sqrt(mean((map_mat_m2_n100[, 3] - cri_pts[3]) ^ 2)) * 100

# -----------------
## DGP-beta(2, 3)
# -----------------
map_mat_m2_beta_2_3_n100 <- matrix(0, nrow = no_data, 3)
for (i in 1:no_data) {
    x_1 <- map_est_beta_2_3_n100[[i]][map_est_beta_2_3_n100[[i]] > b0 & 
                                          map_est_beta_2_3_n100[[i]] < b1]
    map_mat_m2_beta_2_3_n100[i, 1] <- mean(unlist(x_1))
    
    x_2 <- map_est_beta_2_3_n100[[i]][map_est_beta_2_3_n100[[i]] > b1 & 
                                          map_est_beta_2_3_n100[[i]] < b2]
    map_mat_m2_beta_2_3_n100[i, 2] <- mean(unlist(x_2))
    
    x_3 <- map_est_beta_2_3_n100[[i]][map_est_beta_2_3_n100[[i]] > b2 & 
                                          map_est_beta_2_3_n100[[i]] < b3]
    map_mat_m2_beta_2_3_n100[i, 3] <- mean(unlist(x_3))
}

sqrt(mean((map_mat_m2_beta_2_3_n100[, 1] - cri_pts[1]) ^ 2)) * 100
sqrt(mean((map_mat_m2_beta_2_3_n100[, 2] - cri_pts[2]) ^ 2)) * 100
sqrt(mean((map_mat_m2_beta_2_3_n100[, 3] - cri_pts[3]) ^ 2)) * 100



# -----------------
## STS
# -----------------
kovac_est_mat_m2 <- matrix(0, nrow = no_data, 3)
na_count_1_kovac <- 0
na_count_2_kovac <- 0
na_count_3_kovac <- 0
multi_count_1_kovac <- 0
multi_count_2_kovac <- 0
multi_count_3_kovac <- 0

for (i in 1:100) {
    x_1 <- kovac_est[[i]][kovac_est[[i]] > b0 & kovac_est[[i]] < b1]
    if (length(x_1) == 0) na_count_1_kovac <- na_count_1_kovac + 1
    if (length(x_1) > 1) multi_count_1_kovac <- multi_count_1_kovac + 1
    kovac_est_mat_m2[i, 1] <- mean(x_1, na.rm = TRUE)
    
    x_2 <- kovac_est[[i]][kovac_est[[i]] > b1 & kovac_est[[i]] < b2]
    if (length(x_2) == 0) na_count_2_kovac <- na_count_2_kovac + 1
    if (length(x_2) > 1) multi_count_2_kovac <- multi_count_2_kovac + 1
    kovac_est_mat_m2[i, 2] <- mean(x_2, na.rm = TRUE)
    
    x_3 <- kovac_est[[i]][kovac_est[[i]] > b2 & kovac_est[[i]] < b3]
    if (length(x_3) == 0) na_count_3_kovac <- na_count_3_kovac + 1
    if (length(x_3) > 1) multi_count_3_kovac <- multi_count_3_kovac + 1
    kovac_est_mat_m2[i, 3] <- mean(x_3, na.rm = TRUE)
}

sqrt(mean((kovac_est_mat_m2[, 1][!is.nan(kovac_est_mat_m2[, 1])] - cri_pts[1]) ^ 2))
sqrt(mean((kovac_est_mat_m2[, 2] - cri_pts[2]) ^ 2))
sqrt(mean((kovac_est_mat_m2[, 3] - cri_pts[3]) ^ 2))


# -----------------
## NKS
# -----------------
t_mat_m2 <- matrix(0, nrow = no_data, 3)
multi_count_1_nks <- 0
multi_count_2_nks <- 0
multi_count_3_nks <- 0

for (i in 1:no_data) {
    x_1 <- imp_lst[[i]]$x[imp_lst[[i]]$x > b0 & imp_lst[[i]]$x < b1]
    if (length(x_1) > 1) multi_count_1_nks <- multi_count_1_nks + 1
    t_mat_m2[i, 1] <- mean(x_1)
    
    x_2 <- imp_lst[[i]]$x[imp_lst[[i]]$x > b1 & imp_lst[[i]]$x < b2]
    if (length(x_2) > 1) multi_count_2_nks <- multi_count_2_nks + 1
    t_mat_m2[i, 2] <- mean(x_2)
    
    x_3 <- imp_lst[[i]]$x[imp_lst[[i]]$x > b2 & imp_lst[[i]]$x < b3]
    if (length(x_3) > 1) multi_count_3_nks <- multi_count_3_nks + 1
    t_mat_m2[i, 3] <- mean(x_3)
}

c(sqrt(mean((t_mat_m2[, 1] - cri_pts[1]) ^ 2)),
  sqrt(mean((t_mat_m2[, 2] - cri_pts[2]) ^ 2)),
  sqrt(mean((t_mat_m2[, 3] - cri_pts[3]) ^ 2))) * 100




## ============================================================================
## Number of simulations with missing or multiple estimated local extrema in each interval
## ============================================================================

# -----------------
## DGP-beta(1, 1)
# -----------------
sum(unlist(lapply(map_est_n100, function(x) sum(b0 < x & x < b1))) > 1)
sum(unlist(lapply(map_est_n100, function(x) sum(b1 < x & x < b2))) > 1)
sum(unlist(lapply(map_est_n100, function(x) sum(b2 < x & x < b3))) > 1)


# -----------------
## DGP-beta(2, 3)
# -----------------
sum(unlist(lapply(map_est_beta_2_3_n100, function(x) sum(b0 < x & x < b1))) > 1)
sum(unlist(lapply(map_est_beta_2_3_n100, function(x) sum(b1 < x & x < b2))) > 1)
sum(unlist(lapply(map_est_beta_2_3_n100, function(x) sum(b2 < x & x < b3))) > 1)


# -----------------
## STS
# -----------------
multi_count_1_kovac
multi_count_2_kovac
multi_count_3_kovac


# -----------------
## NKS
# -----------------
multi_count_1_nks
multi_count_2_nks
multi_count_3_nks
