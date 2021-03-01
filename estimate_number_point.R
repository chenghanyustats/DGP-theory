### Simulations
################################################################################
# Plot for the number of estimated local extrema
################################################################################
# number of points 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
library(ggplot2)
number_davis_kovac_50 <- c(22, 54, 24, 0, 0, 0, 0, 0)
number_davis_kovac_100 <- c(0, 50, 49, 0, 1, 0, 0, 0)
number_davis_kovac_500 <- c(0, 0, 99, 1, 0, 0, 0, 0)
number_davis_kovac_1000 <- c(0, 0, 100, 0, 0, 0, 0, 0)

number_kovac_50 <- c(28, 49, 23, 0, 0, 0, 0, 0)
number_kovac_100 <- c(0, 50, 49, 0, 1, 0, 0, 0)
number_kovac_500 <- c(0, 0, 98, 1, 1, 0, 0, 0)
number_kovac_1000 <- c(0, 0, 99, 0, 1, 0, 0, 0)


number_nks_50 <- c(0, 0, 51, 10, 26, 10, 2, 1, 0, 0, 1)  ## 11 elements. the 11th element is 1
number_nks_100 <- c(0, 0, 44, 10, 28, 12, 5, 1)
number_nks_500 <- c(0, 0, 45, 9, 29, 13, 3, 1)
number_nks_1000 <- c(0, 0, 47, 13, 33, 7, 0, 0)

number_dgp_50 <- c(0, 0, 40, 49, 11, 0, 0, 0)
number_dgp_100 <- c(0, 0, 47, 48, 5, 0, 0, 0)
number_dgp_500 <- c(0, 0, 95, 5, 0, 0, 0, 0)
number_dgp_1000 <- c(0, 0, 99, 1, 0, 0, 0, 0)

number_dgp_beta_2_3_50 <- c(0, 0, 91, 7, 1, 1, 0, 0)
number_dgp_beta_2_3_100 <- c(0, 0, 86, 13, 1, 0, 0, 0)
number_dgp_beta_2_3_500 <- c(0, 0, 99, 1, 0, 0, 0, 0)
number_dgp_beta_2_3_1000 <- c(0, 0, 100, 0, 0, 0, 0, 0)

number_gmm_50 <- c(0, 0, 87, 3, 10, 0, 0, 0)
number_gmm_100 <- c(0, 0, 95, 5, 0, 0, 0, 0)
number_gmm_500 <- c(0, 0, 99, 1, 0, 0, 0, 0)
number_gmm_1000 <- c(0, 0, 100, 0, 0, 0, 0, 0)

# =============================================================================

method_name <- c("TS", "STS", "NKS", 
                 "DGP beta(1, 1)", "DGP beta(2, 3)", "GMM Approximation")


# -------
number_df_50 <- data.frame("number" = rep(1:8, 6),
                            "method" = rep(method_name, each = 8),
                            "count" = c(number_davis_kovac_50,
                                        number_kovac_50,
                                        number_nks_50[1:8],
                                        number_dgp_50,
                                        number_dgp_beta_2_3_50,
                                        number_gmm_50))

ggplot(number_df_50, aes(x = number, y = count, fill = method)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "Sample size is 50",
         fill = "Method") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) +
    scale_fill_viridis_d()
# -------
number_df_100 <- data.frame("number" = rep(1:8, 6),
                            "method" = rep(method_name, each = 8),
                            "count" = c(number_davis_kovac_100,
                                        number_kovac_100,
                                        number_nks_100,
                                        number_dgp_100,
                                        number_dgp_beta_2_3_100,
                                        number_gmm_100))

ggplot(number_df_100, aes(x = number, y = count, fill = method)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "Sample size is 100",
         fill = "Method") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) +
    # scale_y_continuous(limits = c(0, 100)) + 
    scale_fill_viridis_d()
# -------
number_df_500 <- data.frame("number" = rep(1:8, 6),
                           "method" = rep(method_name, each = 8),
                           "count" = c(number_davis_kovac_500,
                                       number_kovac_500,
                                       number_nks_500[1:8],
                                       number_dgp_500,
                                       number_dgp_beta_2_3_500,
                                       number_gmm_500))

ggplot(number_df_500, aes(x = number, y = count, fill = method)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "Sample size is 500",
         fill = "Method") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) +
    scale_fill_viridis_d()
# -------

number_df_1000 <- data.frame("number" = rep(1:8, 6),
                             "method" = rep(method_name, each = 8),
                             "count" = c(number_davis_kovac_1000,
                                         number_kovac_1000,
                                         number_nks_1000,
                                         number_dgp_1000,
                                         number_dgp_beta_2_3_1000,
                                         number_gmm_1000))

ggplot(number_df_1000, aes(x = number, y = count, fill = method)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "Sample size is 1000",
         fill = "Method") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) +
    scale_fill_viridis_d()




## maybe this

## Taut String (same for other methods)
number_df_ts <- data.frame("number" = rep(1:8, 4),
                           "size" = rep(c("50", "100", "500", "1000"), each = 8),
                           "count" = c(number_davis_kovac_50,
                                       number_davis_kovac_100,
                                       number_davis_kovac_500,
                                       number_davis_kovac_1000))
ggplot(number_df_ts, aes(x = number, y = count, 
                         fill = factor(size, levels = c("50", "100", "500", "1000")))) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "Taut String",
         fill = "Size") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) + 
    scale_fill_viridis_d()

# -----

number_df_sts <- data.frame("number" = rep(1:8, 4),
                           "size" = rep(c("50", "100", "500", "1000"), each = 8),
                           "count" = c(number_kovac_50,
                                       number_kovac_100,
                                       number_kovac_500,
                                       number_kovac_1000))
ggplot(number_df_sts, aes(x = number, y = count, 
                          fill = factor(size, levels = c("50", "100", "500", "1000")))) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "Smooth Taut String",
         fill = "Size") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) + 
    scale_fill_viridis_d()


# -----


number_df_dgp <- data.frame("number" = rep(1:8, 4),
                            "size" = rep(c("50", "100", "500", "1000"), each = 8),
                            "count" = c(number_dgp_50,
                                        number_dgp_100,
                                        number_dgp_500,
                                        number_dgp_1000))
ggplot(number_df_dgp, aes(x = number, y = count, 
                          fill = factor(size, levels = c("50", "100", "500", "1000")))) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "DGP beta(1, 1)",
         fill = "Size") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) + 
    scale_fill_viridis_d()

# -----


number_df_dgp_2_3 <- data.frame("number" = rep(1:8, 4),
                            "size" = rep(c("50", "100", "500", "1000"), each = 8),
                            "count" = c(number_dgp_beta_2_3_50,
                                        number_dgp_beta_2_3_100,
                                        number_dgp_beta_2_3_500,
                                        number_dgp_beta_2_3_1000))
ggplot(number_df_dgp_2_3, aes(x = number, y = count, 
                          fill = factor(size, levels = c("50", "100", "500", "1000")))) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "DGP beta(2, 3)",
         fill = "Size") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) + 
    scale_fill_viridis_d()

# -----


number_df_nks <- data.frame("number" = rep(1:8, 4),
                                "size" = rep(c("50", "100", "500", "1000"), each = 8),
                                "count" = c(number_nks_50[1:8],
                                            number_nks_100,
                                            number_nks_500,
                                            number_nks_1000))
ggplot(number_df_nks, aes(x = number, y = count, 
                              fill = factor(size, levels = c("50", "100", "500", "1000")))) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "Nonlinear Kernel Smoothing",
         fill = "Size") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) + 
    scale_fill_viridis_d()


# -----


number_df_gmm <- data.frame("number" = rep(1:8, 4),
                                "size" = rep(c("50", "100", "500", "1000"), each = 8),
                                "count" = c(number_gmm_50,
                                            number_gmm_100,
                                            number_gmm_500,
                                            number_gmm_1000))
ggplot(number_df_gmm, aes(x = number, y = count, 
                              fill = factor(size, levels = c("50", "100", "500", "1000")))) + 
    geom_bar(stat = "identity", position = "dodge") + 
    labs(x = "Number of local extrema",
         y = "Count",
         title = "Gaussian Mixture Model Approximation",
         fill = "Size") +
    scale_x_continuous(limits = c(0.5, 8.5), breaks=1:8) + 
    scale_fill_viridis_d()








