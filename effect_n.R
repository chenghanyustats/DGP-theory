## change path
load(file = "./data/sim_data_n1000.RData", verbose = TRUE)

# Loading objects:
#     YY
# sig
# x_a
# x_b
# cri_pts
# no_data

# first replicate
x_1000 <- YY[[1]]$x
y_1000 <- YY[[1]]$y
# x_1000 <- lapply(YY, function(i) i$x)
# y_1000 <- lapply(YY, function(i) i$y)

## downsampling
n_vec <- c(20, 50, 100, 200, 500, 1000)
down_sample_lst <- vector("list", length = 6)
for (i in 1:5) {
    # sample_idx <- sample(1:1000, n_vec[i], replace = FALSE)
    sample_idx <- round(seq(1, 1000, length.out = n_vec[i]), 0)
    down_sample_lst[[i]]$x <- x_1000[sample_idx]
    down_sample_lst[[i]]$y <- y_1000[sample_idx]
}
down_sample_lst[[6]]$x <- x_1000
down_sample_lst[[6]]$y <- y_1000


## load results
load("./post_t_effect_n.RData", verbose = TRUE)


grid_t <- seq(x_a, x_b, length.out = 400)


par(mfrow = c(1, 1))
par(mar = c(4, 4, 2, 1))
plot(grid_t, post_prob_gp_hist_lst[[6]], ylab = "p(t | y)", xlab = "t", type = 'l',
     main = "Effect of n on posteiror distribution of t beta(1, 1)",
     ylim = c(0, 0.06), lwd = 2)
for (i in 1:5) {
    lines(grid_t, post_prob_gp_hist_lst[[i]], col = i+1, lwd = 2)
}
abline(v = cri_pts, lty = 2)
legend("topright", paste("n =", n_vec), col = c(2:6, 1), lwd = rep(2, 6),
       bty = "n")
# ------------------------
plot(grid_t, post_prob_gp_hist_lst_beta_2_3[[6]], ylab = "p(t | y)", xlab = "t", type = 'l',
     main = "Effect of n on posteiror distribution of t beta(2, 3)",
     ylim = c(0, 0.125), lwd = 2)
for (i in 1:5) {
    lines(grid_t, post_prob_gp_hist_lst_beta_2_3[[i]], col = i+1, lwd = 2)
}
abline(v = cri_pts, lty = 2)
legend("topright", paste("n =", n_vec), col = c(2:6, 1), lwd = rep(2, 6),
       bty = "n")

