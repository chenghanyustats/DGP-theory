## Song (2006) functions
library(KernSmooth)
LQfit <- function(x, y, h) {
    # x <- data[ , 1]
    # y <- data[ , 2]
    l <- length(x)
    beta0 <- rep(0, l)
    beta1 <- rep(0, l)
    for(i in 1:l) {
        x.reg <- x - x[i]
        w <- dnorm(x - x[i], 0, h)
        fit <- lm(y ~ x.reg + I(x.reg ^ 2), weights = w)
        beta0[i] <- fit$coe[1]
        beta1[i] <- fit$coe[2]
    }
    beta <- cbind(beta0, beta1)
    return(beta)
}

AddCI <- function(x, y, h, beta, alpha) {
    # x <- data[, 1]
    # y <- data[, 2]
    l <- length(x)
    beta0 <- beta[, 1]
    beta1 <- beta[, 2]
    ##Test for equilibruim points
    w <- rep(0, l)
    diag <- rep(0, l)
    upperCI <- rep(0, l)
    lowerCI <- rep(0, l)
    
    upperCI_f <- rep(0, l)
    lowerCI_f <- rep(0, l)
    
    se <- rep(0, l)
    se_f <- rep(0, l)
    
    z <- rep(0, l)
    p <- rep(0 ,l)
    options(object.size = 1000000000)
    
    ##Estimate sigma^2
    for(i in 1:l){
        Xi <- cbind(rep( 1, l), x - x[i], (x - x[i]) ^ 2)
        ker <- dnorm(Xi[, 2], 0, h)
        A <- matrix(0, ncol = 3, nrow = 3)
        A[1, 1] <- sum(ker)
        A[1, 2] <- ker %*% Xi[, 2]
        A[2, 1] <- A[1, 2]
        A[1, 3] <- ker %*% Xi[, 3]
        A[2, 2] <- A[1, 3]
        A[3, 1] <- A[1, 3]
        A[2, 3] <- ker %*% Xi[, 2] ^ 3
        A[3, 2] <- A[2, 3]
        A[3, 3] <- ker %*% Xi[, 3] ^ 2
        B <- solve(A)[1, ]
        C <- rbind(ker, ker * Xi[, 2], ker * Xi[,3])
        wi <- B %*% C
        diag[i] <- wi[i]
        w <- rbind(w, wi)
    }
    w <- w[2:( l + 1), ]
    second <- sum(w ^ 2)
    first <- 2 * sum(diag)
    v <- first - second
    vari <- 1 / ( l - v ) * sum((y - beta0) ^ 2)
    
    ##Calculate the 95% confidence band
    for(i in 1:l) {
        X <- cbind(rep(1, l), x - x[i], (x - x[i]) ^ 2)
        kernel <- dnorm(X[, 2], 0, h)
        An <- matrix(0, ncol = 3, nrow = 3)
        Bn <- matrix(0, ncol = 3, nrow = 3)
        An[1, 1] <- sum(kernel) / l
        An[1, 2] <- kernel %*% X[, 2] / l
        An[2, 1] <- An[1, 2]
        An[1, 3] <- kernel %*% X[, 3] / l
        An[2, 2] <- An[1, 3]
        An[3, 1] <- An[1, 3]
        An[2, 3] <- kernel %*% X[, 2] ^ 3 / l
        An[3, 2] <- An[2, 3]
        An[3, 3] <- kernel %*% X[, 3] ^ 2 / l
        kernel2 <- kernel ^ 2
        Bn[1, 1] <- sum(kernel2) / l / l
        Bn[1, 2] <- kernel2 %*% X[, 2] / l / l
        Bn[2, 1] <- Bn[1, 2]
        Bn[1, 3] <- kernel2 %*% X[, 3] / l / l
        Bn[2, 2] <- Bn[1, 3]
        Bn[3, 1] <- Bn[1, 3]
        Bn[2, 3] <- kernel2 %*% X[, 2] ^ 3 / l / l
        Bn[3, 2] <- Bn[2, 3]
        Bn[3, 3] <- kernel2 %*% X[, 3] ^ 2 / l / l
        sol <- solve(An)
        temp <- sol %*% Bn %*% sol
        temp2 <- temp[2, 2]
        temp1 <- temp[1, 1]
        se[i] <- sqrt(vari * temp2)
        se_f[i] <- sqrt(vari * temp1)
        
        z[i] <- abs(beta1[i] / se[i])
        p[i] <- (1 - pnorm(z[i])) * 2
        upperCI[i] <- beta1[i] + qnorm(1 - alpha/2) * se[i]
        lowerCI[i] <- beta1[i] - qnorm(1 - alpha/2) * se[i]
        
        upperCI_f[i] <- beta0[i] + qnorm(1 - alpha/2) * se_f[i]
        lowerCI_f[i] <- beta0[i] - qnorm(1 - alpha/2) * se_f[i]
    }
    upperCI <- round(upperCI, 5)
    lowerCI <- round(lowerCI, 5)
    
    upperCI_f <- round(upperCI_f, 5)
    lowerCI_f <- round(lowerCI_f, 5)
    
    
    p <- round(p, 5)
    CIp <- cbind(upperCI, lowerCI, p)
    return(CIp)
}

AddCI_f <- function(x, y, h, beta, alpha) {
    # x <- data[, 1]
    # y <- data[, 2]
    l <- length(x)
    beta0 <- beta[, 1]
    beta1 <- beta[, 2]
    ##Test for equilibruim points
    w <- rep(0, l)
    diag <- rep(0, l)
    upperCI <- rep(0, l)
    lowerCI <- rep(0, l)
    
    upperCI_f <- rep(0, l)
    lowerCI_f <- rep(0, l)
    
    se <- rep(0, l)
    se_f <- rep(0, l)
    
    z <- rep(0, l)
    p <- rep(0 ,l)
    options(object.size = 1000000000)
    
    ##Estimate sigma^2
    for(i in 1:l){
        Xi <- cbind(rep( 1, l), x - x[i], (x - x[i]) ^ 2)
        ker <- dnorm(Xi[, 2], 0, h)
        A <- matrix(0, ncol = 3, nrow = 3)
        A[1, 1] <- sum(ker)
        A[1, 2] <- ker %*% Xi[, 2]
        A[2, 1] <- A[1, 2]
        A[1, 3] <- ker %*% Xi[, 3]
        A[2, 2] <- A[1, 3]
        A[3, 1] <- A[1, 3]
        A[2, 3] <- ker %*% Xi[, 2] ^ 3
        A[3, 2] <- A[2, 3]
        A[3, 3] <- ker %*% Xi[, 3] ^ 2
        B <- solve(A)[1, ]
        C <- rbind(ker, ker * Xi[, 2], ker * Xi[,3])
        wi <- B %*% C
        diag[i] <- wi[i]
        w <- rbind(w, wi)
    }
    w <- w[2:( l + 1), ]
    second <- sum(w ^ 2)
    first <- 2 * sum(diag)
    v <- first - second
    vari <- 1 / ( l - v ) * sum((y - beta0) ^ 2)
    
    ##Calculate the 95% confidence band
    for(i in 1:l) {
        X <- cbind(rep(1, l), x - x[i], (x - x[i]) ^ 2)
        kernel <- dnorm(X[, 2], 0, h)
        An <- matrix(0, ncol = 3, nrow = 3)
        Bn <- matrix(0, ncol = 3, nrow = 3)
        An[1, 1] <- sum(kernel) / l
        An[1, 2] <- kernel %*% X[, 2] / l
        An[2, 1] <- An[1, 2]
        An[1, 3] <- kernel %*% X[, 3] / l
        An[2, 2] <- An[1, 3]
        An[3, 1] <- An[1, 3]
        An[2, 3] <- kernel %*% X[, 2] ^ 3 / l
        An[3, 2] <- An[2, 3]
        An[3, 3] <- kernel %*% X[, 3] ^ 2 / l
        kernel2 <- kernel ^ 2
        Bn[1, 1] <- sum(kernel2) / l / l
        Bn[1, 2] <- kernel2 %*% X[, 2] / l / l
        Bn[2, 1] <- Bn[1, 2]
        Bn[1, 3] <- kernel2 %*% X[, 3] / l / l
        Bn[2, 2] <- Bn[1, 3]
        Bn[3, 1] <- Bn[1, 3]
        Bn[2, 3] <- kernel2 %*% X[, 2] ^ 3 / l / l
        Bn[3, 2] <- Bn[2, 3]
        Bn[3, 3] <- kernel2 %*% X[, 3] ^ 2 / l / l
        sol <- solve(An)
        temp <- sol %*% Bn %*% sol
        temp2 <- temp[2, 2]
        temp1 <- temp[1, 1]
        se[i] <- sqrt(vari * temp2)
        se_f[i] <- sqrt(vari * temp1)
        
        z[i] <- abs(beta1[i] / se[i])
        p[i] <- (1 - pnorm(z[i])) * 2
        upperCI[i] <- beta1[i] + qnorm(1 - alpha/2) * se[i]
        lowerCI[i] <- beta1[i] - qnorm(1 - alpha/2) * se[i]
        
        upperCI_f[i] <- beta0[i] + qnorm(1 - alpha/2) * se_f[i]
        lowerCI_f[i] <- beta0[i] - qnorm(1 - alpha/2) * se_f[i]
    }
    upperCI <- round(upperCI, 5)
    lowerCI <- round(lowerCI, 5)
    
    upperCI_f <- round(upperCI_f, 5)
    lowerCI_f <- round(lowerCI_f, 5)
    
    
    p <- round(p, 5)
    CI_f <- cbind(upperCI_f, lowerCI_f, se_f)
    return(CI_f)
}

Curve <- function(x, y, beta, CIp, imp) {
    # x <- data[ , 1]
    # y <- data[ , 2]
    l <- length(x)
    ma <- max(x)
    beta0 <- beta[ , 1]
    beta1 <- beta[ , 2]
    upperCI <- CIp[ , 1]
    lowerCI <- CIp[ , 2]
    # par(mfrow = c( 2, 1))
    ## nonparametric fitted curve with 95% CI
    plot(x, y, xlim = c(0, ma), type = 'p', ylim = c(-.9, 0.9),
         xlab = 'x', ylab = 'y', pch = 19)
    lines(x, beta0, lty = 4, lwd = 3, col = "green4")
    # title('Nonparametric Fitted Curve')
    pos <-as.vector(imp$position)
    vx <-as.vector(imp$x)
    ll <- length(vx)
    vx1 <-as.vector(imp$lCI)
    vx2 <-as.vector(imp$uCI)
    vy <- beta0[pos]
    
    for (i in 1:ll){
        x1 <- vx1[i]
        x2 <- vx2[i]
        y1 <- vy[i]
        polygon(c(x1, x1, x2, x2), c(0, 0.02, 0.02, 0) + 0.02*i, col = i,
                border = "white")
        points(vx[i], 0.02*i, pch = 18, col = i, cex = 1.5)
    }
    
    # ##first-derivative curve with 95% CI
    # plot(x, beta1, type = 'l', xlim = c(0, ma), xlab = 'x', 
    #      ylab = 'First derivative')
    # title( 'Estimated First Derivative Curve with 95% Confidence Interval')
    # for (i in 1:( l - 1)){
    #     x1 <- x[i]
    #     x2 <- x[i + 1]
    #     y1 <- lowerCI[i]
    #     y2 <- lowerCI[i + 1]
    #     y3 <- upperCI[i]
    #     y4 <- upperCI[i + 1]
    #     polygon(c(x1, x1, x2, x2), c(y1, y3, y2, y4), col = "lightgrey")
    # }
    # cross <- rep(0, ll)
    # par(new = T)
    # plot(vx, cross, pch = "x", xlim = c( 0, ma), ylim = c( -2, 2), 
    #      xlab = '', ylab = '')
    # abline(0)
    # par(new = T)
    # plot(x, beta1, type = 'l', xlim = c( 0, ma), ylim = c( -2, 2),
    #      xlab = 'Chromosomal coordinates (in kb)', ylab = 'First derivative'
}

find_song_der_zero_idx <- function(all_info_lst, no_data, is.print = FALSE) {
    song_der_zero_idx_lst <- list()
    for(i in 1:no_data) {
        der_zero_idx <- which(diff(sign(all_info_lst[[i]][, 3]))!= 0)
        if (is.print) print(der_zero_idx)
        der_zero_idx_adj <- rep(0, length(der_zero_idx))
        
        ## choose beta1 closet to zero
        for (j in 1:length(der_zero_idx)) {
            
            if(abs(all_info_lst[[i]][der_zero_idx[j], 3]) > abs(all_info_lst[[i]][der_zero_idx[j]+1, 3])) {
                der_zero_idx_adj[j] <- der_zero_idx[j] + 1
            } else {
                der_zero_idx_adj[j] <- der_zero_idx[j]
            }
        }
        ## no end points
        end_idx <- der_zero_idx_adj %in% c(1, length(all_info_lst[[i]][, 3]))
        song_der_zero_idx_lst[[i]] <- der_zero_idx_adj[!end_idx]
    }
    song_der_zero_idx_lst
}
