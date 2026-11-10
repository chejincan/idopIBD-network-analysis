# Power equation fitting

library(readr)
library(idopNetwork)
library(parallel)
library(patchwork)
library(ggplot2)
library(reshape2)

power_equation_base <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    min_value = min(y[y != 0])
    lmFit <- lm(log(y + runif(1, min = 0, max = min_value)) ~ log(x))
    coefs <- coef(lmFit)
    a <- exp(coefs[1])
    b <- coefs[2]
    model <- try(nls(y ~ a * x^b, start = list(a = a, b = b), control = nls.control(maxiter = 10000,
        minFactor = 1e-200)))
    if ("try-error" %in% class(model)) {
        result = NULL
    } else {
        result = model
    }
    return(result)
}

power_equation_all <- function(x, y) {
    maxit <- 1000L
    result <- power_equation_base(x, y)
    iter <- 1
    while (is.null(result) && iter <= maxit) {
        iter <- iter + 1
        try(result <- power_equation_base(x, y))
    }
    return(result)
}

power_equation_fit <- function(data) {
    n <- 30L
    trans <- log10
    available_cores <- parallel::detectCores(logical = FALSE)
    thread <- if (is.na(available_cores)) 2L else max(1L, min(4L, available_cores - 1L))

    data = data[, order(colSums(data)), drop = FALSE]
    if (is.null(trans)) {
        X = colSums(data)
        trans_data = data
    } else {
        X = trans(colSums(data) + 1)
        trans_data = trans(data + 1)
    }
    colnames(trans_data) = X

    core.number <- thread
    cl <- makeCluster(getOption("cl.cores", core.number))
    on.exit(stopCluster(cl), add = TRUE)
    clusterExport(
        cl, c("power_equation_all", "power_equation_base", "trans_data", "X"),
        envir = environment()
    )
    all_model = parLapply(
        cl = cl, 1:nrow(data),
        function(c) power_equation_all(X, trans_data[c, ])
    )

    names(all_model) = rownames(data)
    no = which(sapply(all_model, length) >= 1)
    all_model2 = all_model[no]
    data2 = data[no, , drop = FALSE]
    trans_data2 = trans_data[no, , drop = FALSE]

    new_x = seq(min(X), max(X), length = n)
    power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), USE.NAMES = TRUE))
    power_fit = t(vapply(
        all_model2, predict, newdata = data.frame(x = new_x),
        FUN.VALUE = numeric(n), USE.NAMES = TRUE
    ))
    colnames(power_fit) = new_x

    result = list(
        original_data = data2, trans_data = trans_data2,
        power_par = power_par, power_fit = power_fit, Time = X
    )
    return(result)
}

data_match <- function(result1, result2, result3) {
    matchname = intersect(
        intersect(rownames(result1$original_data), rownames(result2$original_data)),
        rownames(result3$original_data)
    )

    new_result1 = list(
        original_data = result1$original_data[matchname, , drop = FALSE],
        trans_data = result1$trans_data[matchname, , drop = FALSE],
        power_par = result1$power_par[matchname, , drop = FALSE],
        power_fit = result1$power_fit[matchname, , drop = FALSE],
        Time = result1$Time
    )
    new_result2 = list(
        original_data = result2$original_data[matchname, , drop = FALSE],
        trans_data = result2$trans_data[matchname, , drop = FALSE],
        power_par = result2$power_par[matchname, , drop = FALSE],
        power_fit = result2$power_fit[matchname, , drop = FALSE],
        Time = result2$Time
    )
    new_result3 = list(
        original_data = result3$original_data[matchname, , drop = FALSE],
        trans_data = result3$trans_data[matchname, , drop = FALSE],
        power_par = result3$power_par[matchname, , drop = FALSE],
        power_fit = result3$power_fit[matchname, , drop = FALSE],
        Time = result3$Time
    )
    result = list(dataset1 = new_result1, dataset2 = new_result2, dataset3 = new_result3)
    return(result)
}

power_equation_plot <- function(result) {
    data1 = as.matrix(result$trans_data)
    data2 = as.matrix(result$power_fit)
    sample_size <- min(9L, nrow(data1))
    if (sample_size < 1) {
        stop("No fitted data available to plot.")
    }

    selected = sample(seq_len(nrow(data1)), sample_size)
    df_original = reshape2::melt(as.matrix(data1[selected, , drop = FALSE]))
    df_fit = reshape2::melt(as.matrix(data2[selected, , drop = FALSE]))
    df_original$time = as.numeric(as.character(df_original$Var2))
    df_fit$time = as.numeric(as.character(df_fit$Var2))
    label_data = df_fit[!duplicated(df_fit$Var1), , drop = FALSE]

    p <- ggplot() +
        geom_point(
            df_original, mapping = aes(x = time, y = value, colour = Var1, group = Var1),
            size = 2.5, show.legend = FALSE
        ) +
        geom_line(
            df_fit, mapping = aes(x = time, y = value, colour = Var1, group = Var1),
            linewidth = 0.8, show.legend = FALSE
        ) +
        geom_text(
            label_data, mapping = aes(x = time, y = value, label = Var1, colour = Var1),
            hjust = -0.1, check_overlap = TRUE, show.legend = FALSE, size = 3
        ) +
        xlab("Habitat Index") + ylab("Niche Index") +
        theme_bw() +
        theme(
            axis.title = element_text(size = 15),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10, hjust = 0),
            panel.spacing = unit(0.0, "lines"),
            plot.margin = unit(c(1, 1, 1, 1), "lines"),
            strip.background = element_blank(),
            plot.background = element_blank()
        )
    return(p)
}

df <- read.csv(file.choose())
rownames(df) <- df[, 1]
df <- df[, -1, drop = FALSE]

df1<-df[,1:20]
df2<-df[,21:68]
df3<-df[,69:93]

r1 = power_equation_fit(df1)
r2 = power_equation_fit(df2)
r3 = power_equation_fit(df3)
res1 <- data_match(result1 = r1, result2 = r2, result3 = r3)

print(power_equation_plot(r1))
save(res1, file = file.choose(new = TRUE))
