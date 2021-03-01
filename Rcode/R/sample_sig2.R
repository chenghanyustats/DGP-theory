sample_sig2 <- function(y, x_vec, theta, Kff, Kdf, Kdd, 
                        ga_shape = 5, ga_rate = 5) {
    n <- length(y)
    
    # Kff <- se_ker(H0 = H0, tau = theta[1], h = theta[2])
    
    # Kdf <- computeCovDer1(idx1 = sample_t, idx2 = x_vec, 
    #                       tau = theta[1], h = theta[2]) 
    # Kdd <- computeCovDer2(idx1 = sample_t, tau = theta[1], h = theta[2])

    A <- Kff - emulator::quad.form.inv(Kdd, Kdf) + diag(n)
    
    # rinvgamma(n, shape, rate = 1, scale = 1/rate)

    invgamma::rinvgamma(1, shape = ga_shape + n / 2,
                        rate = ga_rate + (1 / 2) * quad.form.inv(A, y))
}
# 
# 
# sample_t <- function(y, x_vec, sample_sig2, theta, H0) {
#     
# }


log_dgbeta_prior <- function(shape1, shape2, t, 
                             a, b, mixture_prob = NULL) {
    n_mixture <- length(shape1)
    
    if (n_mixture == 1) {
        if (shape1 == 1 && shape2 == 1) {
            if (length(a) == 1) {
                return(log_dgbeta((a + b) / 2, shape1 = shape1,
                                  shape2 = shape2, 
                                  a = a, b = b))
            } else {
                val <- 0
                for (i in 1:length(a)) {
                    val <- val + log(1 / (b[i] - a[i]))
                }
                return(val)
            }

        }
    }
    
    if (n_mixture == 1) {
        if (shape1 != 1 || shape2 != 1) {
            return(log_dgbeta(t, shape1 = shape1,
                                        shape2 = shape2,
                                        a = a, b = b))
        }
    } else {
        prior_den <- sum(mixture_prob * exp(log_dgbeta(t, shape1, shape2,
                                          a = a, b = b)))
        return(log(prior_den))
    }
}



log_post_den_given_sig <- function(t, x, theta, y, Kff, sig,
                                   shape1, shape2, a, b, mixture_prob = NULL) {
    
    
    Kdf <- computeCovDer1(idx1 = t, idx2 = x, tau = theta[1], h = theta[2]) 
    
    Kdd <- computeCovDer2(idx1 = t, tau = theta[1], h = theta[2])
    
    log_den <- -log_mar_lik_gp_der_given_sig(theta = theta, y = y, 
                                             x_vec = x, 
                                             der_vec = t,
                                             Kff = Kff,
                                             Kdf = Kdf,
                                             Kdd = Kdd,
                                             sig = sig)   
    
    log_prior_den <- log_dgbeta_prior(shape1 = shape1, 
                                      shape2 = shape2, t = t, 
                                      a = a, b = b, 
                                      mixture_prob = mixture_prob)
    
    return(list(log_post_den = log_den + log_prior_den,
                Kdf = Kdf,
                Kdd = Kdd))
}


log_prior_den <- log_dgbeta_prior(shape1 = 1, 
                                  shape2 = 1, t = sample_t_mat[1, ], 
                                  a = c(0, 1), b = c(1, 2), 
                                  mixture_prob = NULL)

# log_dgbeta(c(0.5, 1.5), shape1 = 1,
#            shape2 = 1, 
#            a = c(0, 1), b = c(1, 2))
# 
# log_dgbeta

