################################################################################
## Simulation of confidence interval for t in Theorem 3 (ii)
################################################################################
## load packages and functions
library(emulator)
library(ftnonpar)
library(matrixcalc)
library(R.utils)
library(parallel)
library(doParallel)
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


## get muf
get_mufprime <- function(y, x, t, H0, tau, h, lambda) {
    n <- length(y)
    Kff_gp <- se_ker(H0 = H0, tau = tau, h = h)
    A_gp <- Kff_gp + diag((n * lambda), n)
    A_inv <- solve(A_gp)
    Kdf <- computeCovDer1(idx1 = t, idx2 = x, 
                          tau = tau, h = h) 
    A_inv_y <- A_inv %*% y
    mufprime <- Kdf %*% A_inv_y
}


get_t <- function(b, a) {
    for (j in 1:length(b)) {
        if (j == 1) {
            idx <- which.min(abs(b[j] - a))
            b[idx] <- b[j] 
            idx_old <- idx
        } else {
            idx <- which.min(abs(b[j] - a))
            if (idx_old == idx) {
                if (abs(b[j] - a[idx_old]) < abs(b[idx_old] - a[idx_old])) {
                    b[idx] <- b[j] 
                    idx_old <- idx
                } else {
                    next
                }
            } else {
                b[idx] <- b[j] 
                idx_old <- idx
            }
        }
    }
    return(b[1:length(a)])
}


## load replicates of data
################################################################################
load("./data/sim_data_n1000.RData", verbose = TRUE)

# -----------------
## create needed objects
# -----------------
H0_diff <- lapply(YY, function(d) {
    outer(as.vector(d$x), as.vector(d$x),
          FUN = function(x1, x2) (x1 - x2))
})

idx <- seq(x_a, x_b, length.out = 500)
x_test <- sort(seq(x_a, x_b, length.out = 100))
n_test <- length(x_test)
n <- length(YY[[1]]$x)
grid_t <- seq(x_a, x_b, length.out = 100000)


## estimate hyperparameters
cl <- makeCluster(6, type = "FORK")
registerDoParallel(cl)
system.time(
    EB_gp <- foreach(k = 1:no_data) %dopar% {
        res <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
                             LB = c(0.0001, 0.0001, 0.0001),
                             UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                             control = list(TOL = 1e-5, trace = 0),
                             y = YY[[k]]$y, H0 = H0_diff[[k]])
        res$par
    }
)
stopCluster(cl)



## Confidence interval construction
################################################################################

ci_mat <- matrix(0, nrow = no_data, ncol = 2)
lst <- list("ci90" = ci_mat, "ci95" = ci_mat, "ci99" = ci_mat)
ci_lst <- list("t1" = lst, "t2" = lst, "t3" = lst, "t4" = lst, "t5" = lst)

ci_lst_bon <- list(lst, lst, lst, lst, lst)

names(ci_lst_bon) <- c("t1", "t2", "t3", "t4", "t5")
names(ci_lst_bon$t1) <- names(ci_lst_bon$t2) <- names(ci_lst_bon$t3) <- 
    names(ci_lst_bon$t4) <- names(ci_lst_bon$t5) <- c("alpha_01", "alpha_005", "alpha_001")

t_est_mat <- matrix(0, no_data, 3)

## margin of error
for (k in 1:no_data) { 
    cat("k =", k, "\r")
    x <- YY[[k]]$x
    y <- YY[[k]]$y
    sig_gp <- EB_gp[[k]][1]
    tau_gp <- EB_gp[[k]][2]
    h_gp <- EB_gp[[k]][3]
    lambda_gp <- sig_gp ^ 2 / (n * tau_gp ^ 2)
    Kff_gp <- se_ker(H0 = H0_diff[[k]], tau = 1, h = h_gp)
    A_gp <- Kff_gp + diag((n * lambda_gp), n)
    A_inv <- solve(A_gp)
    A_inv_y <- A_inv %*% y
    muf <- Kff_gp %*% A_inv_y
    
    ## --------
    ## find t
    ## --------
    ## mu_fprime
    Kdf_grid <- computeCovDer1(idx1 = grid_t, idx2 = x, 
                               tau = 1, h = h_gp)
    mufprime <- Kdf_grid %*% A_inv_y

    idx0 <- which(diff(sign(mufprime)) != 0)
    idx1 <- which(diff(sign(mufprime)) != 0) + 1
    closeone <- apply(abs(cbind(mufprime[idx0], mufprime[idx1])), 1, which.min)
    idx <- idx0 + closeone - 1
    t_all <- grid_t[idx]
    t_est <- get_t(t_all, cri_pts)
    t_est_mat[k, ] <- t_est
    
    ### prepare terms for CI
    Kdf <- computeCovDer1(idx1 = t_est, idx2 = x, 
                          tau = 1, h = h_gp) 
    del <- Kdf %*% A_inv %*% muf  ## 3 by 1 
    Kddf <- computeCovDer2(idx1 = t_est, idx2 = x, tau = 1, h = h_gp) 
    mu_fprime_1 <- Kddf %*% A_inv_y ## 3 by 1 

    ## sqrt term
    A_inv_d <- A_inv %*% t(Kdf)
    B <- sqrt(diag(crossprod(A_inv_d)))
    
    ## confidence interval
    for (j in 1:length(t_est)) {
        t_adj <- t_est[j] + del[j] / abs(mu_fprime_1[j])
        E <- sig * B[j] / abs(mu_fprime_1[j])
        mar_err <- qnorm(c(0.95, 0.975, 0.995)) * E
        mar_err_bon <- qnorm(1 - (1 - c(0.95, 0.975, 0.995)) / 3) * E
        
        outer(mar_err, c(-1, 1))
        ci <- t_adj + outer(mar_err, c(-1, 1))
        ci_bon <-  t_adj + outer(mar_err_bon, c(-1, 1))
        
        ci_lst[[j]][[1]][k, ] <- ci[1, ]
        ci_lst[[j]][[2]][k, ] <- ci[2, ]
        ci_lst[[j]][[3]][k, ] <- ci[3, ]
        
        ci_lst_bon[[j]][[1]][k, ] <- ci_bon[1, ]
        ci_lst_bon[[j]][[2]][k, ] <- ci_bon[2, ]
        ci_lst_bon[[j]][[3]][k, ] <- ci_bon[3, ]
        
    }
}


## CI summary and plotting

out_ci <- matrix(0, 5, 3, 
                 dimnames = list(c("t1", "t2", "t3", "t4", "t5"),
                                 c("alpha_01", "alpha_005", "alpha_001")))

out_ci_bon <- lapply(vector(length = 5, mode = "list"), 
                     function(x) {setNames(as.list(rep(0, 3)),
                                           c("alpha_01", "alpha_005", "alpha_001"))})
names(out_ci_bon) <- paste0("t", 1:5)
    
    
par(mar = c(4, 4, 2, 1), mgp = c(2.5, 1, 0))
par(mfrow = c(1, 3))
alpha <- c(0.1, 0.05, 0.01)
for (al in 1:3) {
    for (j in 1:length(t_est)) {
        idx <- which.min(abs(t_est[j] - cri_pts))
        true_t <- cri_pts[idx]
        ci_min <- min(ci_lst[[j]][[al]])
        ci_max <- max(ci_lst[[j]][[al]])
        ci_lower <- ci_lst[[j]][[al]][, 1]
        ci_upper <- ci_lst[[j]][[al]][, 2]
        plot(rep(true_t, no_data), seq(no_data), type = 'l', 
             xlim = c(ci_min, ci_max), 
             col = "yellow", las = 1,
             xlab = "t", ylab = "No. Data", 
             main = paste0((1-alpha[al])*100, "% CI for t", j, " (abs)"))
        out <- (true_t < ci_lst[[j]][[al]][, 1] | true_t > ci_lst[[j]][[al]][, 2])
        out_bon <- (true_t < ci_lst_bon[[j]][[al]][, 1] | true_t > ci_lst_bon[[j]][[al]][, 2])
        out_ci_bon[[j]][[al]] <- out_bon
        out_ci[j, al] <- sum(out)
        segments(x0 = ci_lower, y0 = 1:no_data, x1 = ci_upper, y1 = 1:no_data, 
                 col = '#003366', lwd = 1)
        segments(x0 = ci_lower[out], y0 = (1:no_data)[out], x1 = ci_upper[out], 
                 y1 = (1:no_data)[out], col = 'red', lwd = 1)
        abline(v = true_t, col = "#FFCC00", lwd = 2)
    }
    
}


# out_ci_bon_which_1 <- lapply(out_ci_bon$t1, function(x) which(x == TRUE))
# 
# out_ci_bon_which_2 <- lapply(out_ci_bon$t2, function(x) which(x == TRUE))
# 
# out_ci_bon_which_3 <- lapply(out_ci_bon$t3, function(x) which(x == TRUE))
# 
# length(union(union(out_ci_bon_which_1[[1]], out_ci_bon_which_2[[1]]), 
#              out_ci_bon_which_3[[1]]))
# 
# length(union(union(out_ci_bon_which_1[[2]], out_ci_bon_which_2[[2]]), 
#              out_ci_bon_which_3[[2]]))
# 
# length(union(union(out_ci_bon_which_1[[3]], out_ci_bon_which_2[[3]]), 
#              out_ci_bon_which_3[[3]]))
# 
# 
# intersect(intersect(out_ci_bon_which_1[[1]], out_ci_bon_which_2[[1]]), 
#           out_ci_bon_which_3[[1]])
# 
# intersect(intersect(out_ci_bon_which_1[[2]], out_ci_bon_which_2[[2]]), 
#           out_ci_bon_which_3[[2]])
# 
# intersect(intersect(out_ci_bon_which_1[[3]], out_ci_bon_which_2[[3]]), 
#           out_ci_bon_which_3[[3]])







