

'usethis::use_package("lmtest")
usethis::use_package("car")
usethis::use_pipe(export = TRUE)
usethis::use_package("stringr")'

'library(car)
library(car)
library(tidyverse)
lm1 <- lm(data=iris,Sepal.Length~Petal.Length+Petal.Width)'

#' @title  All-in-One function for regression diagnostics.
#'
#' @description Computes statistical tests and plots for linear regression diagnostics.
#' Plots are also generated for a multivariate normal comparison set.
#'
#' @export
#'
#' @param lm Please enter lm object with the following notation: lm(data = data, Y~Xi+Xn)
#' @param interactionYN TRUE or FALSE, are interactions included? Will calculate GVIF instead of VIF.
#' @param sig_level A value between 0 and 1. If a test's p-value is lager than sig_level, the assumption is met.
#' @param cutoff_vif Cutoof-VIF threshold for identifying multicoliniearity.

cus_reg_diag <-
  function(lm,
           summaryOnly = TRUE,
           interactionYN = FALSE,
           sig_level = .05,
           cutoff_vif = 5) {
    #Funktion für die Regressionsdiagnostik.

    #Vorbereitung des Data.Frames, welches der Output der Funktion ist.
    tests <-
      c(
        "Formula",
        "Shapiro-Wilk-Normality",
        "BP-Homoskedasticity",
        "RESET-Linearity",
        "Cook's d-Outliers",
        "VIF-Multicoliniearity"
      )
    conclusion <- c("NA", "NA", "NA", "NA", "NA", "NA")
    value <- c("", "NA", "NA", "NA", "NA", "NA")
    threshold <- c("", "<.05", "<.05", "<.05", "NA", ">5")


    n <- nrow(lm$model)
    preds <- ncol(lm$model)

    if (!interactionYN) {
      if (ncol(lm$model) > 2) {
        var_inf <- (car::vif(lm))
        value[6] <- max(var_inf)
      }
    }

    if (interactionYN) {
      if (ncol(lm$model) > 2) {
        var_inf <- (car::vif(lm, type = "predictor"))
        value[6] <- max(var_inf$`GVIF^(1/(2*Df))`) ^ 2
        warning("Using VIF type = predictor, output is GVIF^(1/(2*Df))^2")
      }
    }


    cdist <- lm %>% cooks.distance()


    cdist_cutoff <-
      qf(.5, length(coef(lm)) + 1, nobs(lm) - length(coef(lm)) - 1)

    threshold[5] <- paste(">", cdist_cutoff)

    value[5] <- max(cdist)



    sw_p <- lm %>% rstandard() %>% shapiro.test() %>% {
      .$p.value
    }
    bp_p <- lm %>% lmtest::bptest() %>% {
      .$p.value
    }
    rt_p <- lm %>% lmtest::resettest() %>% {
      .$p.value
    }


    value[2] <- sw_p
    value[3] <- bp_p
    value[4] <- rt_p






    conclusion[1] <- ""


    if (sw_p < sig_level) {
      conclusion[2] <- "Violated"
    } else{
      conclusion[2] <- "Okay"
    }

    if (bp_p < sig_level) {
      conclusion[3] <- "Violated"
    } else{
      conclusion[3] <- "Okay"
    }

    if (rt_p < sig_level) {
      conclusion[4] <- "Violated"
    } else{
      conclusion[4] <- "Okay"
    }

    if (value[5] %>% as.numeric() > cdist_cutoff) {
      conclusion[5] <- "Violated"
    } else{
      conclusion[5] <- "Okay"
    }

    if (value[6] %>% as.numeric() > cutoff_vif & value[6] != "NA") {
      conclusion[6] <- "Violated"
    } else{
      conclusion[6] <- "Okay"
    }
    value[1] <- paste("n= ", n, "k= ", length(coef(lm)) - 1)

    result <- data.frame(tests, conclusion, value, threshold)
    print(result)

    if (any(cdist > cdist_cutoff)) {
      model <-
        lm$model[which(cdist > cdist_cutoff), ] %>% cbind("cdist" = cdist[which(cdist >
                                                                                  cdist_cutoff)])
      model[order(model$cdist, decreasing = TRUE), ] %>% print()
    }

    #Simuliert multivariat normalverteilte Daten basierend auf Covarianzmatrix, Stichprobengröße und Mittelwerte des Modells.

    'if(any(sapply(lm$model,class)=="factor")){

    lm_model_dummycodiert <- list()
    factor_postion <- which(sapply(lm$model,class)=="factor")

    for (i in factor_postion){
      lm_model_dummycodiert[which(factor_postion==i)] <- car::dumm

    }
  }'

    lm_model_matrix_minus_intercept <-
      (lm$model %>% model.matrix(data = lm$model))[, -1]
    lm_dataframe_codiert <-
      data.frame(lm$model[1], lm_model_matrix_minus_intercept)
    lm_cov <- cov(lm_dataframe_codiert)

    model_data <- lm_dataframe_codiert %>% as.data.frame()
    model_data_sim <- model_data
    model_data_sim[, ] <-
      MASS::mvrnorm(
        n = nrow(lm_dataframe_codiert),
        Sigma = lm_cov,
        mu = colSums(lm_dataframe_codiert) / nrow(lm_dataframe_codiert),
        empirical = TRUE
      )

    formula <- lm$call %>% as.character() %>% stringr::str_split(pattern = "=")

    model_data_sim_new <- model_data_sim

    for (i in 1:ncol(lm_dataframe_codiert)) {
      if (lm_dataframe_codiert[, i] %>% table() %>% names() %in% c("0", "1") %>%
          any() & lm_dataframe_codiert[, i] %>% table() %>% length() < 3) {
        model_sim_mean_temp <- mean(model_data_sim[, i])
        model_data_sim_new[, i] <- model_data_sim[, i] %>%
          cut(breaks = c(-Inf, model_sim_mean_temp, +Inf),
              labels = c(0, 1))
      }

    }


    lm2 <- lm(formula[[2]], data = model_data_sim_new)
    if (length(formula) == 2) {
      stop(
        "Please use the following structure for lm(): data = data, y ~ x1 + xn..; Execution stopped."
      )
    }

    par(mfrow = c(1, 2))
    par(mar = c(2, 2, 2, 2))
    plot(lm, 1)
    plot(lm2, 1)
    car::residualPlot(lm)
    mtext("Residuals vs fitted,\n echtes Modell")
    car::residualPlot(lm2)
    mtext("Residuals vs fitted,\n simuliertes Modell")
    linearity_plot_result <-
      readline("Assume Linearity based on left plot? (y/n):")

    car::qqPlot(lm)
    mtext("QQ-Plot,\n echtes Modell")
    car:::qqPlot(lm2$residuals)
    mtext("QQ-Plot,\n simuliertes Modell")
    normality_plot_result <-
      readline("Assume Normality based on left plot? (y/n):")

    plot(lm, 4, caption = "")
    mtext("QQ-Plot,\n echtes Modell")
    plot(lm2, 4, caption = "")
    mtext("QQ-Plot,\n simuliertes Modell")

    normality_plot_result <-
      readline("Identify Outliers based on left plot. (0 = non):")

    mtext("QQ-Plot,\n simuliertes Modell")
    #Todo, debugging





    lm %>% rstandard() %>% shapiro.test() %>% print()
    lmtest::bptest(lm) %>% print()
    lmtest::resettest(lm) %>% print()
    if (preds > 2) {
      print(var_inf)
    }


  }
