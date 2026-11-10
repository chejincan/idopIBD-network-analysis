# Module network reconstruction

library(ggplot2)
library(patchwork)
library(Cairo)
library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(reshape2)
library(readxl)
library(cowplot)
library(ggalt)

get_sum_line <- function(i, k) {
    df_data <- df[df[[ncol(df)]] == i, , drop = FALSE]
    if (k == 1) {
        df_data <- df_data[, seq_along(times1), drop = FALSE]
    } else if (k == 2) {
        df_data <- df_data[, length(times1) + seq_along(times2), drop = FALSE]
    } else {
        df_data <- df_data[, length(times1) + length(times2) + seq_along(times3), drop = FALSE]
    }
    df_line <- colSums(df_data)
    return(df_line)
}

power_equation <- function(par, x) {
    y <- par[1] * x^par[2]
    y
}

power_par <- function(y, times) {
    x <- as.numeric(times)
    y <- as.numeric(y)
    lmFit <- lm(log(y + runif(1, min = 0, max = 0.1)) ~ log(x))
    coefs <- coef(lmFit)
    a <- exp(coefs[1])
    b <- coefs[2]
    tmp <- c(a, b)
    par_est <- function(par, x) {
        sum((y - power_equation(par, x))^2)
    }
    r <- optim(tmp, par_est, x = as.numeric(times), method = "Nelder-Mead")
    return(r$par)
}

power_fit <- function(x, par) {
    t(sapply(1:nrow(par), function(c) par[c, 1] * x^par[c, 2]))
}

LOP_rk4 <- function(y0, h, d_LOP_fit, effect, times) {
    f <- function(x, y) {
        dy = do.call(sum, polynomial.values(polynomials = d_LOP_fit, x = x))
        dy
    }

    res_y <- c()
    res_y[1] <- y0

    for (i in 1:(length(times) - 1)) {
        x0 <- scaleX(times, u = -1, v = 1)[i]
        k1 <- f(x0, y0)
        k2 <- f(x0 + h/2, y0 + h/2 * k1)
        k3 <- f(x0 + h/2, y0 + h/2 * k2)
        k4 <- f(x0 + h, y0 + h * k3)
        res_y[i + 1] <- y0 + h/6 * (k1 + 2 * k2 + 2 * k3 + k4) * effect[i]
        y0 <- res_y[i + 1]
    }
    return(res_y)
}

get_effect <- function(pars, effect, times, order, y0) {
    LOP <- legendre.polynomials(order, normalized = F)

    d_LOP_fit <- sapply(1:length(pars), function(c) pars[c] * polynomial.derivatives(LOP)[[c + 1]])

    h <- scaleX(times, u = -1, v = 1)[2] - scaleX(times, u = -1, v = 1)[1]

    y_fit <- LOP_rk4(y0, h, d_LOP_fit, effect, times)
    return(y_fit)
}

get_legendre_par <- function(data, times) {
    order <- 2L

    get_interaction <- function(data, col) {
        clean_data <- data
        m <- clean_data[, col]
        M <- clean_data[, -col, drop = FALSE]
        x_matrix <- M
        x_matrix <- as.matrix(x_matrix)
        name <- colnames(clean_data)
        cv_folds <- max(5L, min(10L, nrow(data) - 1L))

        ridge1_cv <- cv.glmnet(
            x = x_matrix, y = m, type.measure = "mse", family = "gaussian",
            nfolds = cv_folds, alpha = 0
        )
        best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])

        fit_res <- cv.glmnet(
            x = x_matrix, y = m, type.measure = "mse", family = "gaussian",
            nfolds = cv_folds, alpha = 1,
            penalty.factor = 1/best_ridge_coef, keep = TRUE
        )
        best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)

        selected_coef <- as.matrix(best_alasso_coef1)[-1, 1]
        selected <- which(selected_coef != 0)
        return_obj <- list(
            ind.name = name[col],
            dep.name = rownames(best_alasso_coef1)[-1][selected],
            coefficient = unname(selected_coef[selected])
        )
        return(return_obj)
    }

    module_relationship <- lapply(1:ncol(data), function(c) get_interaction(data, c))

    ode_optimize <- function(pars, ind, dep, times, data, order, effect) {
        ind_pars <- matrix(pars, ncol = order)[1, ]
        dep_pars <- matrix(pars, ncol = order)[-1, ]
        inital_value <- data[, ind][1]
        ind_effect <- get_effect(ind_pars, data[, ind], times, order, inital_value)
        if (!length(dep)) {
            y <- ind_effect
        } else if (length(dep) == 1) {
            dep_effect <- get_effect(as.numeric(dep_pars), data[, dep], times, order, 0)
            y <- ind_effect + dep_effect
        } else {
            dep_effect <- sapply(seq_along(dep), function(c) get_effect(
                dep_pars[c, ], data[, dep[c]], times, order, 0
            ))
            y <- ind_effect + rowSums(dep_effect)
        }
        ssr <- sum((data[, ind] - y)^2)

        alpha = 0.005
        ridge <- sum((data[, ind] - y)^2) + alpha * (sum(pars^2) + 0 * sum(pars^2))

        if (min(ind_effect) > 0) {
            return(ridge)
        } else {
            return(ridge * 500)
        }
    }

    get_value <- function(effect, relationship, times, order) {
        ind <- relationship[[1]]
        dep <- relationship[[2]]
        ind_no <- which(colnames(effect) == ind)
        dep_no <- match(dep, colnames(effect))

        init_pars <- rep(1e-05, (length(ind_no) + length(dep_no)) * order)
        result <- optim(
            init_pars, ode_optimize, ind = ind_no, dep = dep_no,
            times = times, data = effect, order = order,
            method = "BFGS", control = list(maxit = 500)
        )
        par_after <- matrix(result$par, length(ind_no) + length(dep_no), order)
        return(par_after)
    }

    core.number <- detectCores()
    cl <- makeCluster(getOption("cl.cores", core.number))
    on.exit(stopCluster(cl), add = TRUE)
    clusterEvalQ(cl, {
        library(orthopolynom)
    })
    clusterExport(
        cl, c(
            "LOP_rk4", "get_value", "ode_optimize", "get_effect",
            "times", "get_interaction", "module_relationship", "order"
        ),
        envir = environment()
    )
    lop_par <- pblapply(
        1:ncol(data),
        function(c) get_value(data, module_relationship[[c]], times, order),
        cl = cl
    )
    return(list(lop_par, module_relationship))
}

get_output <- function(relationship, par, effect, times) {
    order <- 2L
    output <- list()
    output[[1]] <- relationship[[1]]
    output[[2]] <- relationship[[2]]
    output[[3]] <- par[1, ]
    output[[4]] <- if (nrow(par) > 1) par[-1, , drop = FALSE] else
        matrix(numeric(0), nrow = 0, ncol = ncol(par))

    ind_no <- which(colnames(effect) == output[[1]])
    dep_no <- match(output[[2]], colnames(effect))
    inital_value <- effect[, ind_no][1]
    ind_effect <- get_effect(as.numeric(output[[3]]), effect[, ind_no], times, order, inital_value)

    if (!length(dep_no)) {
        dep_effect <- matrix(numeric(0), nrow = length(times), ncol = 0)
    } else if (length(dep_no) == 1) {
        dep_effect <- get_effect(as.numeric(output[[4]]), effect[, dep_no], times, order, 0)
    } else {
        dep_effect <- sapply(seq_along(dep_no), function(c) get_effect(
            as.numeric(output[[4]][c, ]), effect[, dep_no[c]], times, order, 0
        ))
        colnames(dep_effect) <- dep_no
    }

    all_effect <- cbind(ind_effect, dep_effect)
    effect_mean <- apply(all_effect, 2, mean)
    output[[5]] <- effect_mean
    output[[6]] <- all_effect
    return(output)
}

get_figure_all <- function(k, all_net, times) {
    col <- c("#d34752", "#5161ae", "#8ed0aa")
    get_curve_data <- function(i) {
        d1 <- data.frame(all_net[[i]][[6]], check.names = FALSE)
        colnames(d1) <- c("ind_effect", all_net[[i]][[2]])
        d1$sum <- rowSums(d1)
        d1$x <- times
        return(d1)
    }

    result <- lapply(seq_len(k), function(c) get_curve_data(c))

    myplot <- function(i, data) {
        d <- data[[i]]
        p1 <- ggplot(d, aes(x = x)) +
            geom_line(d, mapping = aes(x = x, y = sum), color = col[2], linewidth = 1.5) +
            geom_line(d, mapping = aes(x = x, y = ind_effect), color = col[1], linewidth = 1.4)

        if (ncol(d) > 3) {
            for (j in seq_len(ncol(d) - 3)) {
                p1 <- p1 +
                    geom_line(
                        aes(y = .data[[colnames(d)[1 + j]]]),
                        color = col[3], linewidth = 1.2, alpha = 1
                    ) +
                    annotate(
                        "text", x = as.numeric(d$x[length(d$x)]),
                        y = d[nrow(d), 1 + j], label = colnames(d)[1 + j], size = 4
                    )
            }
        }

        p1 <- p1 +
            theme(
                panel.grid = element_blank(),
                panel.background = element_rect(color = "black", fill = "transparent"),
                legend.title = element_blank()
            ) +
            theme(plot.title = element_text(vjust = -6, hjust = 0.5)) +
            xlab(NULL) + ylab(NULL) +
            theme(plot.margin = unit(c(0, 0, 0, 0), "lines")) +
            geom_hline(yintercept = 0, linewidth = 0.6)
        return(p1)
    }

    p <- lapply(seq_len(k), function(i) myplot(i, result))
    pp <- plot_grid(plotlist = p, ncol = 6, nrow = ceiling(k / 6))
    ggsave(
        file.choose(new = TRUE), ggdraw(pp),
        width = 10, height = 6, units = "in", limitsize = FALSE
    )
}

get_net_output <- function(i) {
    get_after <- function(i) {
        temp <- data.frame(from = as.character(i[[2]]), to = rep(i[[1]], length(i[[2]])),
            dep_effect = as.numeric(i[[5]][1 + seq_along(i[[2]])]))
        return(temp)
    }

    links <- get_after(i)
    get_link_color <- function(i) {
        tmp <- links$dep_effect[i]
        if (tmp >= 0) {
            tmp2 <- "+"
        } else {
            tmp2 <- "-"
        }
        return(tmp2)
    }

    links$effect_type <- if (nrow(links))
        vapply(seq_len(nrow(links)), get_link_color, character(1)) else character(0)
    links$dep_effect <- abs(links$dep_effect)
    nodes <- data.frame(id = i[[1]], independent = i[[5]][1], stringsAsFactors = FALSE)
    nodes$received_effect <- if (nrow(links)) sum(links$dep_effect) else 0
    return(list(links, nodes))
}

load(file.choose())
df <- cluster_result[["cluster2"]]
times1 <- cluster_result[["Time1"]]
times2 <- cluster_result[["Time2"]]
times3 <- cluster_result[["Time3"]]

line1 <- t(sapply(sort(unique(df[[ncol(df)]])), function(c) get_sum_line(c, 1)))
line2 <- t(sapply(sort(unique(df[[ncol(df)]])), function(c) get_sum_line(c, 2)))
line3 <- t(sapply(sort(unique(df[[ncol(df)]])), function(c) get_sum_line(c, 3)))

par1 <- t(sapply(1:nrow(line1), function(c) power_par(line1[c, ], times1)))
par2 <- t(sapply(1:nrow(line2), function(c) power_par(line2[c, ], times2)))
par3 <- t(sapply(1:nrow(line3), function(c) power_par(line3[c, ], times3)))

times1 <- seq(min(times1), max(times1), length = 30)
times2 <- seq(min(times2), max(times2), length = 30)
times3 <- seq(min(times3), max(times3), length = 30)

line_data1 <- t(power_fit(times1, par1))
line_data2 <- t(power_fit(times2, par2))
line_data3 <- t(power_fit(times3, par3))

colnames(line_data1) <- paste("M", sep = "", sort(unique(df[[ncol(df)]])))
colnames(line_data2) <- paste("M", sep = "", sort(unique(df[[ncol(df)]])))
colnames(line_data3) <- paste("M", sep = "", sort(unique(df[[ncol(df)]])))

all_par1 <- get_legendre_par(data = line_data1, times = times1)
all_par2 <- get_legendre_par(data = line_data2, times = times2)
all_par3 <- get_legendre_par(data = line_data3, times = times3)

module_relationship1 <- all_par1[[2]]
module_relationship2 <- all_par2[[2]]
module_relationship3 <- all_par3[[2]]

all_net1 <- pblapply(1:ncol(line_data1), function(c) get_output(
    module_relationship1[[c]], all_par1[[1]][[c]], line_data1, times = times1
))
all_net2 <- pblapply(1:ncol(line_data2), function(c) get_output(
    module_relationship2[[c]], all_par2[[1]][[c]], line_data2, times = times2
))
all_net3 <- pblapply(1:ncol(line_data3), function(c) get_output(
    module_relationship3[[c]], all_par3[[1]][[c]], line_data3, times = times3
))

get_figure_all(k = ncol(line_data3), all_net = all_net3, times = times3)

result1 <- lapply(seq_along(all_net1), function(c) get_net_output(all_net1[[c]]))
result2 <- lapply(seq_along(all_net2), function(c) get_net_output(all_net2[[c]]))
result3 <- lapply(seq_along(all_net3), function(c) get_net_output(all_net3[[c]]))

for (group in c("HC", "UC", "CD")) {
    result <- switch(group, HC = result1, UC = result2, CD = result3)
    b <- data.frame()
    b2 <- data.frame()
    for (i in seq_along(result)) {
        a <- result[[i]][[1]]
        b <- rbind(b, a)
        a2 <- result[[i]][[2]]
        b2 <- rbind(b2, a2)
    }

    write.csv(b, file.choose(new = TRUE), row.names = FALSE)
    write.csv(b2, file.choose(new = TRUE), row.names = FALSE)
}
