
# dataset <- read.csv("./data/Raw_ERP.csv", header = TRUE)
# names(dataset) <- c("id", "group", "region", "electrode", "condition",
#                     "trial", seq(-100, 896, 4))
# TT <- 250
# erp_data <- dataset[, 7:(TT+6)]
# lwr_idx <- which(names(erp_data) == 0)
# upr_idx <- which(names(erp_data) == 360)
# 
# lwr_idx <- which(names(erp_data) == -100)
# upr_idx <- which(names(erp_data) == 896)

lwr_idx <- which(names(erp_data) == 100)
upr_idx <- which(names(erp_data) == 348)
msec_idx <- seq(-100, 896, 4)
## 21 and 43/44

# number_trials <- seq(2, nrow(erp_data), by = 2)

number_trials <- 1:72


n_erp <- length(seq(-100, 896, 4)[lwr_idx:upr_idx])
avg_erp <- apply(erp_data, 2, mean)

avg_erp[which.max(avg_erp)] - mean(avg_erp)
avg_erp[which.min(avg_erp)] - mean(avg_erp)

erp_x <- seq(0.004, 1, length = 250)


x_a_erp <- min(erp_x[lwr_idx:upr_idx])
x_b_erp <- max(erp_x[lwr_idx:upr_idx])






load("./erp_avg_trial_hpd_result.RData", verbose = TRUE)

# ==========

plot(seq(length(number_trials)), rep(10, length(number_trials)), type = "l", 
     xlab = "Number of trials averaged", yaxt='n', xaxt='n',
     ylab = "msec", cex.axis = 1, tcl = -0.4, ylim = c(x_a_erp, x_b_erp),
     main = "95% HPD interval beta(1, 1) # of trials averaged")
axis(2, at = erp_x[lwr_idx:upr_idx], labels = seq(-100, 896, 4)[lwr_idx:upr_idx], 
     las = 1, tcl = -0.3, cex.axis = 0.9)
axis(1, at = 1:length(number_trials), labels = number_trials, las = 1, tcl = -0.3, cex.axis = 1)
for (k in 1: length(number_trials)) {
    n_clu <- length(erp_hpdi_1_2_lst_1_1[[k]][, 1])
    for (j in 1:n_clu) {
        segments(x0 = k, 
                 y0 = erp_hpdi_1_2_lst_1_1[[k]][j, 1], 
                 x1 = k,
                 y1 = erp_hpdi_1_2_lst_1_1[[k]][j, 2], col = "blue", lwd = 3)
    }
    points(x = rep(k, length(erp_map_est_1_1[[k]])), y = erp_map_est_1_1[[k]],
           col = "red", cex = 0.5, pch = 19)
}
abline(h = erp_map_est_1_1[[72]][2:3], lwd = 2)

# ==========

plot(seq(length(number_trials)), rep(10, length(number_trials)), type = "l", ylab = "msec", yaxt='n', xaxt='n',
     xlab = "Number of trials averaged", cex.axis = 1, tcl = -0.4, ylim = c(x_a_erp, x_b_erp),
     main = "95% HPD interval beta(3, 3) # of trials averaged")
axis(2, at = erp_x[lwr_idx:upr_idx], labels = seq(-100, 896, 4)[lwr_idx:upr_idx], 
     las = 1, tcl = -0.3, cex.axis = 0.9)
axis(1, at = 1:length(number_trials), labels = number_trials, las = 1, tcl = -0.3, cex.axis = 0.8)
for (k in 1: length(number_trials)) {
    n_clu <- length(erp_hpdi_1_2_lst_3_3[[k]][, 1])
    for (j in 1:n_clu) {
        segments(y0 = erp_hpdi_1_2_lst_3_3[[k]][j, 1], 
                 x0 = k, 
                 y1 = erp_hpdi_1_2_lst_3_3[[k]][j, 2], 
                 x1 = k, col = "blue", lwd = 3)
    }
    points(y = erp_map_est_3_3[[k]], x = rep(k, length(erp_map_est_3_3[[k]])), 
           col = "red", cex = 0.5, pch = 19)
}
abline(h = erp_map_est_3_3[[72]], lwd = 2)

# ==========

plot(seq(length(number_trials)), rep(10, length(number_trials)), type = "l", ylab = "msec", yaxt='n', xaxt='n',
     xlab = "Number of trials averaged", cex.axis = 1, tcl = -0.4, ylim = c(x_a_erp, x_b_erp),
     main = "95% HPD interval beta(5, 5) # of trials averaged")
axis(2, at = erp_x[lwr_idx:upr_idx], labels = seq(-100, 896, 4)[lwr_idx:upr_idx], 
     las = 1, tcl = -0.3, cex.axis = 0.9)
axis(1, at = 1:length(number_trials), labels = number_trials, las = 1, tcl = -0.3, cex.axis = 1)
for (k in 1: length(number_trials)) {
    n_clu <- length(erp_hpdi_1_2_lst_5_5[[k]][, 1])
    for (j in 1:n_clu) {
        segments(y0 = erp_hpdi_1_2_lst_5_5[[k]][j, 1], 
                 x0 = k, 
                 y1 = erp_hpdi_1_2_lst_5_5[[k]][j, 2], 
                 x1 = k, col = "blue", lwd = 3)
    }
    points(y = erp_map_est_5_5[[k]], x = rep(k, length(erp_map_est_5_5[[k]])), 
           col = "red", cex = 0.5, pch = 19)
}
abline(h = erp_map_est_5_5[[72]], lwd = 2)



