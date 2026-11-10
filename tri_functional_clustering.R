# Trivariate functional clustering

library(readr)
library(idopNetwork)
library(parallel)
library(patchwork)
library(ggplot2)
library(reshape2)
library(mvtnorm)
library(dplyr)

power_equation <- function(x, power_par) {
    t(sapply(1:nrow(power_par), function(c) power_par[c, 1] * x^power_par[c, 2]))
}

power_equation_base <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    min_value = min(y[y != 0])
    lmFit <- lm(log(y + runif(1, min = 0, max = min_value)) ~ log(x))
    coefs <- coef(lmFit)
    a <- exp(coefs[1])
    b <- coefs[2]
    model <- try(nls(
        y ~ a * x^b, start = list(a = a, b = b),
        control = nls.control(maxiter = 1000, minFactor = 1e-200)
    ))
    if ("try-error" %in% class(model)) {
        result = NULL
    } else {
        result = model
    }
    return(result)
}

power_equation_all <- function(x, y) {
    maxit <- 100L
    result <- power_equation_base(x, y)
    iter <- 1
    while (is.null(result) && iter <= maxit) {
        iter <- iter + 1
        try(result <- power_equation_base(x, y))
    }
    return(result)
}

logsumexp <- function(v) {
    vm = max(v)
    log(sum(exp(v - vm))) + vm
}

get_SAD1_covmatrix <- function(par, n) {
    phi <- par[1]
    gamma <- par[2]
    sigma <- array(dim = c(n, n))

    diag(sigma) <- sapply(1:n, function(c) (1 - phi^(2 * c))/(1 - phi^2))

    sigma[lower.tri(sigma)] <- do.call(c, lapply(
        1:(n - 1), function(c) phi^seq(1:(n - c)) * diag(sigma)[c]
    ))
    sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
    return(gamma^2 * sigma)
}

biQ_function <- function(par, prob_log, omega_log, X, k, n1, n2, n3, times1, times2, times3) {
    n = dim(X)[1]
    X1 = X[, 1:n1]
    X2 = X[, (n1 + 1):(n1 + n2)]
    X3 = X[, (n1 + n2 + 1):(n1 + n2 + n3)]

    par.cov <- par[1:6]
    par.mu <- array(par[-c(1:6)], dim = c(k, 2, 3))

    cov1 = get_SAD1_covmatrix(par.cov[1:2], n1)
    cov2 = get_SAD1_covmatrix(par.cov[3:4], n2)
    cov3 = get_SAD1_covmatrix(par.cov[5:6], n3)
    mu1 <- power_equation(times1, par.mu[, , 1][1:k, ])
    mu2 <- power_equation(times2, par.mu[, , 2][1:k, ])
    mu3 <- power_equation(times3, par.mu[, , 3][1:k, ])
    mvn_log1 <- sapply(1:k, function(c) dmvnorm(X1, mu1[c, ], cov1, log = T))
    mvn_log2 <- sapply(1:k, function(c) dmvnorm(X2, mu2[c, ], cov2, log = T))
    mvn_log3 <- sapply(1:k, function(c) dmvnorm(X3, mu3[c, ], cov3, log = T))

    mvn.log = mvn_log1 + mvn_log2 + mvn_log3
    tmp = sweep(mvn.log, 2, FUN = "+", STATS = prob_log) - omega_log
    Q = -sum(tmp * exp(omega_log))
    return(Q)
}

biget_par_int <- function(X, k, times1, times2, times3, n1, n2, n3) {
    n = dim(X)[1]
    X1 = X[, 1:n1]
    X2 = X[, (n1 + 1):(n1 + n2)]
    X3 = X[, (n1 + n2 + 1):(n1 + n2 + n3)]

    cov.int = c(0.1, IQR(diag(cov(X1))), 0.1, IQR(diag(cov(X2))), 0.1, IQR(diag(cov(X3))))
    init.cluster <- kmeans(X, centers = k, iter.max = 50)
    prob <- table(init.cluster$cluster)/n

    fit1 <- lapply(1:k, function(c) power_equation_all(times1, init.cluster$centers[c, 1:n1]))
    fit2 <- lapply(1:k, function(c) power_equation_all(
        times2, init.cluster$centers[c, (n1 + 1):(n1 + n2)]
    ))
    fit3 <- lapply(1:k, function(c) power_equation_all(
        times3, init.cluster$centers[c, (n1 + n2 + 1):(n1 + n2 + n3)]
    ))

    mu.par.int1 <- t(sapply(fit1, coef))
    mu.par.int2 <- t(sapply(fit2, coef))
    mu.par.int3 <- t(sapply(fit3, coef))
    return_obj <- list(initial_cov_params = cov.int, initial_mu_params = c(mu.par.int1, mu.par.int2,
        mu.par.int3), initial_probibality = prob)
    return(return_obj)
}

bifun_clu <- function(data1, data2, data3, k) {
    iter.max <- 1000L
    parscale <- 0.001

    data1 = data1[, order(colSums(data1))]
    data2 = data2[, order(colSums(data2))]
    data3 = data3[, order(colSums(data3))]

    data = as.matrix(cbind(data1, data2, data3))
    n = dim(data)[1]
    n1 = dim(data1)[2]
    n2 = dim(data2)[2]
    n3 = dim(data3)[2]
    eplison <- 1
    iter <- 0

    times1 <- suppressWarnings(as.numeric(colnames(data1)))
    times2 <- suppressWarnings(as.numeric(colnames(data2)))
    times3 <- suppressWarnings(as.numeric(colnames(data3)))
    if (length(times1) != n1 || length(times2) != n2 || length(times3) != n3 ||
        any(!is.finite(times1)) || any(!is.finite(times2)) || any(!is.finite(times3))) {
        times1 <- as.numeric(colSums(data1))
        times2 <- as.numeric(colSums(data2))
        times3 <- as.numeric(colSums(data3))
    }
    X <- data

    initial.pars <- biget_par_int(X, k, times1, times2, times3, n1, n2, n3)
    X1 = X[, 1:n1]
    X2 = X[, (n1 + 1):(n1 + n2)]
    X3 = X[, (n1 + n2 + 1):(n1 + n2 + n3)]
    par.int <- c(initial.pars$initial_cov_params, initial.pars$initial_mu_params)
    prob_log <- log(initial.pars$initial_probibality)

    while (abs(eplison) > 0.001 && iter <= iter.max) {
        par.mu <- array(par.int[-c(1:6)], dim = c(k, 2, 3))
        par.cov = par.int[1:6]
        cov1 = get_SAD1_covmatrix(par.cov[1:2], n1)
        cov2 = get_SAD1_covmatrix(par.cov[3:4], n2)
        cov3 = get_SAD1_covmatrix(par.cov[5:6], n3)
        mu1 <- power_equation(times1, par.mu[, , 1][1:k, ])
        mu2 <- power_equation(times2, par.mu[, , 2][1:k, ])
        mu3 <- power_equation(times3, par.mu[, , 3][1:k, ])
        mvn_log1 <- sapply(1:k, function(c) dmvnorm(X1, mu1[c, ], cov1, log = T))
        mvn_log2 <- sapply(1:k, function(c) dmvnorm(X2, mu2[c, ], cov2, log = T))
        mvn_log3 <- sapply(1:k, function(c) dmvnorm(X3, mu3[c, ], cov3, log = T))
        mvn_log = mvn_log1 + mvn_log2 + mvn_log3
        mvn = sweep(mvn_log, 2, FUN = "+", STATS = prob_log)
        omega_log = t(sapply(1:n, function(c) mvn[c, ] - logsumexp(mvn[c, ])))
        omega = exp(omega_log)
        LL.mem <- biQ_function(
            par = par.int, prob_log = prob_log, omega_log, X, k,
            n1, n2, n3, times1, times2, times3
        )

        prob_exp = apply(omega_log, 2, logsumexp)
        prob_log = prob_exp - log(n)
        Q.maximization <- try(optim(
            par = par.int, biQ_function, prob_log = prob_log, omega_log = omega_log,
            X = X, k = k, n1 = n1, n2 = n2, n3 = n3,
            times1 = times1, times2 = times2, times3 = times3,
            method = "L-BFGS-B",
            control = list(
                parscale = c(rep(parscale, 6), rep(1, 6 * k)),
                maxit = 100
            )
        ))
        if ("try-error" %in% class(Q.maximization))
            break

        par.hat <- Q.maximization$par
        par.int = par.hat
        LL.next <- biQ_function(
            par = par.int, prob_log = prob_log, omega_log, X, k,
            n1, n2, n3, times1, times2, times3
        )
        eplison <- LL.next - LL.mem
        LL.mem <- LL.next
        iter = iter + 1
    }

    AIC = 2 * (LL.next) + 2 * (length(par.hat) + k - 1)
    BIC = 2 * (LL.next) + log(n) * (length(par.hat) + k - 1)
    omega = exp(omega_log)
    X.clustered <- data.frame(X, apply(omega, 1, which.max), check.names = F)

    return_obj <- list(
        cluster_number = k, Log_likelihodd = LL.mem, AIC = AIC, BIC = BIC,
        cov_par = par.cov, mu_par = par.mu, probibality = exp(prob_log),
        omega = omega, cluster = X.clustered,
        cluster2 = data.frame(data, apply(omega, 1, which.max), check.names = F),
        Time1 = times1, Time2 = times2, Time3 = times3, original_data = data
    )
    return(return_obj)
}

bifun_clu_parallel <- function(data1, data2, data3, start, end) {
    object_count <- min(nrow(data1), nrow(data2), nrow(data3))
    if (!is.finite(object_count) || object_count < 3) {
        stop("Insufficient fitted objects for clustering.")
    }
    if (length(start) != 1L || length(end) != 1L || any(!is.finite(c(start, end))) ||
        any(c(start, end) != as.integer(c(start, end))) ||
        start < 2L || end < start || end >= object_count) {
        stop("Use integer cluster limits with 2 <= start <= end < the number of fitted objects.")
    }
    candidate_k <- seq.int(start, end)

    available_cores <- parallel::detectCores(logical = FALSE)
    thread <- if (is.na(available_cores)) 2L else max(1L, min(4L, available_cores - 1L))
    core.number <- thread
    cl <- makeCluster(getOption("cl.cores", core.number))
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, {
        library(mvtnorm)
    })
    clusterExport(
        cl, c(
            "bifun_clu", "biget_par_int", "biQ_function", "logsumexp",
            "power_equation_base", "power_equation_all",
            "power_equation", "get_SAD1_covmatrix"
        ),
        envir = environment()
    )
    result <- parLapply(cl = cl, candidate_k, function(c) bifun_clu(
        data1 = data1, data2 = data2, data3 = data3, k = c
    ))
    return(result)
}

normalization <- function(x) {
    z <- 0.2
    if (max(x) == min(x)) {
        return(rep(z, length(x)))
    }
    (x - min(x)) / (max(x) - min(x)) + z
}

bifun_clu_plot <- function(result) {
    n <- 30L
    degree <- 1 / 4
    show.legend <- FALSE
    color1 <- "#38E54D"
    color2 <- "#FF8787"
    color3 <- "#87CEEB"
    cluster.result = result
    times1 = cluster.result$Time1
    times2 = cluster.result$Time2
    times3 = cluster.result$Time3

    times1_new = seq(min(times1), max(times1), length = n)
    times2_new = seq(min(times2), max(times2), length = n)
    times3_new = seq(min(times3), max(times3), length = n)
    n1 = length(times1)
    n2 = length(times2)
    n3 = length(times3)
    timesall = c(min(c(times1, times2, times3)), max(c(times1, times2, times3)))
    par.mu = cluster.result$mu_par
    k = cluster.result$cluster_number
    cluster_id = as.integer(cluster.result$cluster[[ncol(cluster.result$cluster)]])
    alpha = tabulate(cluster_id, nbins = k)

    mu.fit1 = power_equation(times1_new, par.mu[, , 1][1:k, ])
    mu.fit2 = power_equation(times2_new, par.mu[, , 2][1:k, ])
    mu.fit3 = power_equation(times3_new, par.mu[, , 3][1:k, ])

    colnames(mu.fit1) = times1_new
    mu.fit1 = melt(as.matrix(mu.fit1))
    colnames(mu.fit1) = c("cluster", "x", "y")
    mu.fit1$x = as.numeric(as.character(mu.fit1$x))

    colnames(mu.fit2) = times2_new
    mu.fit2 = melt(as.matrix(mu.fit2))
    colnames(mu.fit2) = c("cluster", "x", "y")
    mu.fit2$x = as.numeric(as.character(mu.fit2$x))

    colnames(mu.fit3) = times3_new
    mu.fit3 = melt(as.matrix(mu.fit3))
    colnames(mu.fit3) = c("cluster", "x", "y")
    mu.fit3$x = as.numeric(as.character(mu.fit3$x))

    mu.fit1$type = "one"
    mu.fit2$type = "two"
    mu.fit3$type = "three"
    mu.fit = rbind(mu.fit1, mu.fit2, mu.fit3)

    plot.df1 = cluster.result$cluster[, c(1:n1, ncol(cluster.result$cluster))]
    X1 = plot.df1[, -ncol(plot.df1)]
    plot.df1$name = rownames(plot.df1)
    colnames(plot.df1) = c(times1, "cluster", "name")
    plot.df1 = melt(plot.df1, id.vars = c("cluster", "name"))
    colnames(plot.df1) = c("cluster", "name", "x", "y")
    plot.df1$x = as.numeric(as.character(plot.df1$x))
    plot.df1$name = paste0("a_", plot.df1$name)

    plot.df2 = cluster.result$cluster[, c((n1 + 1):(n1 + n2), ncol(cluster.result$cluster))]
    X2 = plot.df2[, -ncol(plot.df2)]
    plot.df2$name = rownames(plot.df2)
    colnames(plot.df2) = c(times2, "cluster", "name")
    plot.df2 = melt(plot.df2, id.vars = c("cluster", "name"))
    colnames(plot.df2) = c("cluster", "name", "x", "y")
    plot.df2$x = as.numeric(as.character(plot.df2$x))
    plot.df2$name = paste0("b_", plot.df2$name)

    plot.df3 = cluster.result$cluster[
        , c((n1 + n2 + 1):(n1 + n2 + n3), ncol(cluster.result$cluster))
    ]
    X3 = plot.df3[, -ncol(plot.df3)]
    plot.df3$name = rownames(plot.df3)
    colnames(plot.df3) = c(times3, "cluster", "name")
    plot.df3 = melt(plot.df3, id.vars = c("cluster", "name"))
    colnames(plot.df3) = c("cluster", "name", "x", "y")
    plot.df3$x = as.numeric(as.character(plot.df3$x))
    plot.df3$name = paste0("c_", plot.df3$name)

    plot.df1$type = "one"
    plot.df2$type = "two"
    plot.df3$type = "three"
    plot.df = rbind(plot.df1, plot.df2, plot.df3)
    df.alpha = data.frame(cluster = 1:k, alpha = normalization(alpha) * degree)
    plot.df = merge(plot.df, df.alpha, by = "cluster")

    name.df = data.frame(
        label = paste0("M", 1:k, "(", alpha, ")"),
        x = mean(timesall),
        y = max(c(
            as.numeric(as.matrix(X1)), as.numeric(as.matrix(X2)),
            as.numeric(as.matrix(X3))
        ) * 0.9),
        cluster = 1:k
    )

    p = ggplot() +
        geom_point(
            plot.df, mapping = aes(x = x, y = y, colour = type, group = name),
            show.legend = show.legend, shape = 1
        ) +
        geom_line(
            mu.fit, mapping = aes(x = x, y = y, colour = type, group = type),
            linewidth = 1.25, show.legend = show.legend
        ) +
        facet_wrap(~cluster) +
        xlab("Habitat Index") + ylab("Niche Index") +
        theme(axis.title = element_text(size = 18)) +
        geom_text(
            name.df, mapping = aes(x = x, y = y, label = label),
            check_overlap = TRUE, size = 5
        ) +
        theme_bw() +
        scale_color_manual(
            name = "Type", values = c("one" = color1, "two" = color2, "three" = color3)
        ) +
        guides(alpha = "none") +
        theme(
            axis.title = element_text(size = 15),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10, hjust = 0),
            panel.spacing = unit(0.0, "lines"),
            plot.margin = unit(c(1, 1, 1, 1), "lines"),
            strip.background = element_blank(),
            plot.background = element_blank(),
            strip.text = element_blank()
        )
    return(p)
}

load(file.choose())
df1 <- res1$dataset1$power_fit
df2 <- res1$dataset2$power_fit
df3 <- res1$dataset3$power_fit

# Define n as the maximum candidate cluster count before running this script.
if (!exists("n", inherits = FALSE)) stop("Define n as the maximum candidate cluster count.")
c1 <- bifun_clu_parallel(df1,df2,df3,start=2,end=n)

print(fun_clu_BIC(c1))
cluster_result <- c1[[which.min(vapply(c1, `[[`, numeric(1), "BIC"))]]
print(bifun_clu_plot(cluster_result))
save(c1, cluster_result, file = file.choose(new = TRUE))
