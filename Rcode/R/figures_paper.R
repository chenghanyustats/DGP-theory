# ################################################################################
# # Reproduce Figures in the paper                                               #
# # Cheng-Han Yu                                                                 #
# ################################################################################
# # ==============================================================================
# # Figure 1: Sample paths from a Gaussian process with (right) and without (left) 
# # first order derivative information. A squared exponential kernel with 
# # parameters $\tau=1$ and $h=1$ was used. Input values are indicated with 
# # plus signs in the plots, and stationary points $t_1 = -4, t_2 = 0,$ and 
# # $t_3 = 4$ are indicated with purple vertical lines in the right-hand plot.
# # ==============================================================================
# # noise_free_path.png
# # ===================
# n_path <- 5
# idx_obs <- sort(runif(5, -5, 5))
# idx_der <- seq(-4, 4, length = 3)
# idx_test <- seq(-5, 5, length = 101)
# y_der <- rep(0, length(idx_der))
# tau <- 1
# phi <- 1
# h <- sqrt(1 / phi)
# sig <- 0
# set.seed(10000)
# # --------------
# K_X <- compute_cov_1d(idx1 = idx_obs, tau = tau, h = h) + 
#     diag(sig ^ 2, length(idx_obs))
# K_XtestX <- compute_cov_1d(idx1 = idx_test, idx2 = idx_obs, tau = tau, h = h)
# K_XtestXtest <- compute_cov_1d(idx1 = idx_test, idx2 = idx_test, tau = tau, 
#                                h = h)
# y_obs <- mvnfast::rmvn(1, rep(0, length(idx_obs)), K_X)
# 
# predict_mean_var <- cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0, 
#                                        obs_2 = y_obs[1, ],
#                                        cov_mat_1 = K_XtestXtest, 
#                                        cov_mat_2 = K_X, cov_mat_12 = K_XtestX)
# 
# predict_val <- mvnfast::rmvn(n_path, predict_mean_var$condMean,
#                              predict_mean_var$condVar)
# # --------------
# y_joint <- c(y_obs[1, ], y_der)
# K_joint <- compute_joint_cov(idx_obs = idx_obs, idx_der = idx_der, sig = sig, 
#                              tau = tau, h = h, w = 0)
# Kffnew <- compute_cov_1d(idx1 = idx_obs, idx2 = idx_test, tau = tau, h = h)
# Kdfnew <- computeCovDer1(idx1 = idx_der, idx2 = idx_test, tau = tau, h = h) 
# Kbar <- rbind(Kffnew, Kdfnew)
# Kfnewfnew <- compute_cov_1d(idx1 = idx_test, idx2 = idx_test, tau = tau, h = h)
# 
# predict_mean_var_der <- cond_norm_mean_var(mean_vec_1 = 0, mean_vec_2 = 0,
#                                            obs_2 = y_joint, 
#                                            cov_mat_1 = Kfnewfnew,
#                                            cov_mat_2 = K_joint,
#                                            cov_mat_12 = t(Kbar))
# 
# predict_val_der <- mvnfast::rmvn(n_path, predict_mean_var_der$condMean, 
#                                  predict_mean_var_der$condVar)
# # --------------
# par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
# gp_plot(idx = idx_obs, idx_test = idx_test, obs = y_obs[1, ], 
#         predict_val = predict_val, ylim = c(-3, 3), is.multiplepath = TRUE,
#         n_path = 5, is.title = FALSE, is.ci = FALSE)
# title(substitute(paste("Paths without derivative: (", tau, ", ", h, ") = (", 
#                        v2, ", ", v3, ")"), list(v2 = tau, v3 = phi)),
#       cex.main = 1)
# 
# gp_plot(idx = idx_obs, idx_test = idx_test, obs = y_obs[1, ], 
#         predict_val = predict_val_der, idx_der = idx_der, 
#         ylim = c(-3, 3), is.multiplepath = TRUE, is.ci = FALSE,
#         n_path = 5, is.title = FALSE, is.derline = TRUE)
# title(substitute(
#     paste("Paths with 1st derivative: (", tau, ", ", h, ") = (", v2, ", ", 
#           v3, ")"), list(v2 = tau, v3 = phi)), cex.main = 1)
# 
# # ==============================================================================
# # Figure 2: Simulated data. {\it Top two rows:} Predictive curves for one 
# # simulated dataset under different models: A standard GPR model and three 
# # different DGP models (oracle, multiple and single), as described in the text. 
# # The true regression function is shown as a red-solid line and the estimated 
# # predictive curves as blue-dashed lines. Dots indicate input values. 
# # {\it Bottom row:} Posterior distributions of $t$, for single and multiple 
# # DGPs, with vertical dark lines indicating the locations of the true 
# # stationary points.
# # ==============================================================================
# 
# ## load simulated datasets for comparing methods
# ## change paths!
# load("./Data/sim_compare_data.rda", verbose = TRUE)
# # load("./results/sim_new.RData", verbose = TRUE)
# 
# ## difference matrix
# H0_diff <- lapply(YY, function(d) {
#     outer(as.vector(d$x), as.vector(d$x),
#           FUN = function(x1, x2) (x1 - x2))
# })
# 
# ## select the second simulated dataset
# k <- 2
# 
# ## empirical Bayes estimmates of kernel parameters for GPR
# eb_gp <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
#                      LB = c(0.0001, 0.0001, 0.0001), 
#                      UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001), 
#                      control = list(TOL = 1e-5, trace = 0),
#                      y = YY[[k]]$y, H0 = H0_diff[[k]])
# res$par
# 
# 
# 
# 
# 
# # pred_hist.png
# # =============
# par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
# plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
#                  mu_test = pred_no_der_lst_new[[k]]$mu_test,
#                  CI_Low_f = pred_no_der_lst_new[[k]]$ci_low,
#                  CI_High_f = pred_no_der_lst_new[[k]]$ci_high,
#                  xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                  ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
#                  plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
#                  pred_lwd = 2, title = paste("GPR"), is.legend = FALSE)
# 
# plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
#                  mu_test = pred_known_a_lst_new[[k]]$mu_test,
#                  CI_Low_f = pred_known_a_lst_new[[k]]$ci_low,
#                  CI_High_f = pred_known_a_lst_new[[k]]$ci_high,
#                  xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                  ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
#                  plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
#                  pred_lwd = 2, title = paste("oracle-DGP"), 
#                  is.legend = FALSE, legend.loc = "bottomleft")
# 
# plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
#                  mu_test = pred_EM_multi_lst_new[[k]]$mu_test,
#                  CI_Low_f = pred_EM_multi_lst_new[[k]]$ci_low,
#                  CI_High_f = pred_EM_multi_lst_new[[k]]$ci_high,
#                  xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                  ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
#                  plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
#                  pred_lwd = 2, title = paste("multiple-DGP"), is.legend = FALSE)
# 
# plot_pred_gp_f_y(x = YY[[k]]$x, y = YY[[k]]$y, idx = idx, x_test = x_test,
#                  mu_test = pred_EM_one_lst_new[[k]]$mu_test,
#                  CI_Low_f = pred_EM_one_lst_new[[k]]$ci_low,
#                  CI_High_f = pred_EM_one_lst_new[[k]]$ci_high,
#                  xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                  ylim = c(0.3, 2.4), is.true.fcn = TRUE, cex = 1, 
#                  plot.type = "p", col.poly = rgb(0, 0, 1, 0.1), 
#                  pred_lwd = 2, title = paste("single-DGP"), is.legend = FALSE)
# 
# hist_t(sample_t_EM_multi_lst_new[[k]], 0, 2, ylim = c(0, 5), col = "blue", 
#        den.line = FALSE)
# title(list(paste("Distribution of t: multiple-DGP"), cex = 1.5))
# abline(v = cri_pts, col = "darkgreen", lty = 1, lwd = 1)
# 
# hist_t(sample_t_EM_lst_new[[k]], 0, 2, ylim = c(0, 5), col = "blue", 
#        den.line = FALSE)
# title(list(paste("Distribution of t: single-DGP"), cex = 1.5))
# abline(v = cri_pts, col = "darkgreen", lty = 1, lwd = 1)
# 
# 
# 
# 
# # ==============================================================================
# # Figure 3: Simulated data: Top row: Predictive curves for n = 50 and sigma = .7. 
# # Bottom row: Predictive curves for n = 10 and sigma = 0.25. Standard GPR is 
# # shown in the plots on the left, and single-DGP in those on the right. 
# # The true regression function is shown as a red-solid line and the estimated 
# # predictive curves as blue-dashed lines. Dots indicate input values. 
# # Vertical dark lines indicate the locations of the true stationary points.
# # ==============================================================================
# load("./Analysis/Data/SimDataSizeNew.Rdata", verbose = TRUE)
# load("./results/StoEM_sim_size_new.RData", verbose = TRUE)
# 
# k <- 5
# par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
# plot_pred_gp_f_y(x = YYn_lst$n10[[k]]$x, y = YYn_lst$n10[[k]]$y, idx = idx, 
#                  x_test = x_test,
#                  mu_test = pred_no_der_lst_10_new[[k]]$mu_test,
#                  CI_Low_f = pred_no_der_lst_10_new[[k]]$ci_low,
#                  CI_High_f = pred_no_der_lst_10_new[[k]]$ci_high,
#                  is.pred.f = TRUE, ylab = "f(x)",
#                  xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                  ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
#                  plot.type = "p", col.poly = rgb(0, 0, 1, 0.2), 
#                  pred_lwd = 2, title = "GPR", 
#                  is.legend = FALSE)
# 
# plot_pred_gp_f_y(x = YYn_lst$n10[[k]]$x, y = YYn_lst$n10[[k]]$y, idx = idx, 
#                  x_test = x_test,
#                  mu_test = pred_EM_one_lst_10_new[[k]]$mu_test,
#                  CI_Low_f = pred_EM_one_lst_10_new[[k]]$ci_low,
#                  CI_High_f = pred_EM_one_lst_10_new[[k]]$ci_high,
#                  is.pred.f = TRUE, ylab = "f(x)",
#                  xlim = c(0, 2), is.der.line = TRUE, cri_pts = cri_pts,
#                  ylim = c(0.3, 2.5), is.true.fcn = TRUE, cex = 1, 
#                  plot.type = "p", col.poly = rgb(0, 0, 1, 0.2), 
#                  pred_lwd = 2, title = "single-DGP", 
#                  is.legend = FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
