


'usethis::use_package("lmtest")
usethis::use_package("car")
usethis::use_pipe(export = TRUE)
usethis::use_package("stringr")'

'library(car)
library(car)
library(tidyverse)
lm1 <- lm(data=iris,Sepal.Length~Petal.Length+Petal.Width+Species)'

#' @title  All-in-One function for regression diagnostics.
#'
#' @description Computes statistical tests and plots for linear regression diagnostics.
#' Plots are also generated for a multivariate normal comparison set.
#'
#' @export
#'
#' @param lm Please enter lm object with the following notation: lm(data = data, Y~Xi+Xn)
#' @param summaryOnly Only shows the statistical test results, omitting plots.
#' @param interactionYN TRUE or FALSE, are interactions included? Will calculate GVIF instead of VIF.
#' @param sig_level A value between 0 and 1. If a test's p-value is lager than sig_level, the assumption is met.
#' @param cutoff_vif Cutoof-VIF threshold for identifying multicoliniearity.
#' @param seed Set seed for reproducible number generation.



cus_reg_diag <-
  function(lm,
           summaryOnly = FALSE,
           interactionYN = FALSE,
           sig_level = .05,
           cutoff_vif = 5,
           seed = 123,
           rMarkdown = FALSE) {
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
        lm$model[which(cdist > cdist_cutoff),] %>% cbind("cdist" = cdist[which(cdist >
                                                                                 cdist_cutoff)])
      model[order(model$cdist, decreasing = TRUE),] %>% print()
    }

    #Simuliert multivariat normalverteilte Daten basierend auf Covarianzmatrix, Stichprobengröße und Mittelwerte des Modells.


    #Für die Simulation wird eine Dummycodierte Matrix der Y unx X Werte benötigt.
    #Daher wird die Modellmatrix erstellt und um den intercep [1,1,1...] gekürzt.
    lm_model_matrix_minus_intercept <-
      (lm$model %>% model.matrix(data = lm$model))[,-1]

    #Jetzt werden die vorhergesagten Werte hinzugefügt.
    lm_dataframe_codiert <-
      data.frame(lm$model[1], lm_model_matrix_minus_intercept)

    #Basierend auf dieser Matrix werden die Covarianzen bestimmt.
    lm_cov <- cov(lm_dataframe_codiert)

    #Die Formatierten Modelldaten werden als Dataframe gespeichert.
    model_data <- lm_dataframe_codiert %>% as.data.frame()

    #Das Dataframe des Modells wird als Vorlage für das simulierte Modellframe genutzt.
    model_data_sim <- model_data

    #Basierend auf n, der Kovarianzmatrix und den Mittelwertsvektor wird eine
    #multivariat normalverteilte Stichprobe generiert.

    set.seed(123)
    model_data_sim[,] <-
      MASS::mvrnorm(
        n = nrow(lm_dataframe_codiert),
        Sigma = lm_cov,
        mu = colSums(lm_dataframe_codiert) / nrow(lm_dataframe_codiert),
        empirical = TRUE
      )

    #Die Modellgleichung des lm  Input Objektes wird zwischengespeichert.
    formula <-
      lm$call %>% as.character() %>% stringr::str_split(pattern = "=")

    #Jetzt besteht noch die Aufgabe die kontinuierlich generierten Werte der
    #dummycodierten Prädiktoren korrekt zu transformieren.
    #Wichtig ist, dass sie nur die Werte 0 oder 1 annehmen und dummycodierte
    #Variablen eines Prädiktors nur maximal eine 1 je Reihe aufweisen.

    #Das simulierte Model-frame stellt Grundlage des neuen Frames dar.
    model_data_sim_new <- model_data_sim

    #Erst wird geprüft, ob ein Faktor als Input vorhanden war.
    #Also, ob die Tranformation überhaupt nötig ist.

    if ("factor" %in% (lm$model %>% sapply(class))) {

      #Falls ein Faktor vorhanden ist, werden alle Faktornamen in einem Vektor
      #gespeichert.

      factor_names <-
        (lm$model %>% sapply(class) == "factor") %>% which() %>% names()

      for (i in 1:ncol(lm_dataframe_codiert)) {

        #Jetzt findet die erste Transformation statt. Wenn eine Variable nur
        #Ausprägungen der Werte 0 und oder 1 hat, wird angekommen,
        #dass die Variable Dummycodiert ist.

        if (any(
          lapply(
            list(
          c("0", "1"),
          c("1", "0"),
          c("0"),
          c("1")
          ),identical,(lm_dataframe_codiert[, i] %>%
            table() %>%
            names()
            ))%>%
          unlist())){

          #Der Variablenmittelwert der simulierten Werte der dummycodierten Variable
          #in einer Variable gespeichert. Dieser stellt den cut-off punkt dar,
          #weil die normalverteilten generierten Daten symmetrisch sind.
          #Dann wird anhand des Mittelwertes die Variable dichtotomisiert.

          model_sim_mean_temp <- mean(model_data_sim[, i])
          model_data_sim_new[, i] <- model_data_sim[, i] %>%
            cut(
              breaks = c(-Inf, model_sim_mean_temp,+Inf),
              labels = c(0, 1)
            )%>%
            as.character()%>%
            as.numeric()
        }

      }

      #Jetzt gilt es noch die Fälle zu beheben, in denen die zufällige
      #Datengenerierung überdurschnittliche Werte für zwei Dummyvariablen
      #eines kategorialen Prädiktors ergeben hat. Leider war es nicht möglich
      #binominale Wertebereiche in der multivariaten Generierung zu beachten.
      #Bei den bisherigen Transformationen hätten die Variablen jetzt z. B.
      #die Werte 1 1.

      for (i in 1:ncol(lm_dataframe_codiert)) {
        dummy_var_connected <- c()

        for (j in 1:length(factor_names)) {

          #Um zu identifizieren, welche dummycodierten Variablen zusammen gehören,
          #Werden für jeden Faktor diejenigen dummyvariablen-Positionen in einem
          #Vektor gespeichert, welche den Faktornamen enthalten.

          if (grepl(factor_names[j], names(lm_dataframe_codiert)) %>% any()) {
            dummy_var_connected <-
              grepl(
                factor_names[j],
                names(lm_dataframe_codiert)) %>%
                  which()

            #Jetzt wird geprüft, ob im tranformierten, simulierten Datensatz für
            #eine Matrix aus n beobachtungen und i = verbundene Dummyvariablen mindestens
            #eine Reihensumme >= 2 beträgt. D. h., dass eine Fall wie z.B. 1 1, 1 0 1 oder 1 1 1
            #vorliegt.
            if((rowSums(model_data_sim_new[,dummy_var_connected])>=2)%>%any()){

              #Wenn das der Fall ist, dann werden alle Reihen durchlaufen, welche diese
              #Eigenschaft aufweisen.
              #Basierend auf den kontinuierlich generierten Daten für diese Variablen
              #wird derjenige Wert identifiziert, welcher maximal ist.
              #Diesem wird im transformierten Datensatz der Wert 1 zugeordnet.
              #Allen anderen Variablen der Gruppe wird der Wert 0 zugeordnet.

              for(k in (rowSums(model_data_sim_new[,dummy_var_connected])>=2)%>%which()){

                max_sim_value <- model_data_sim[k,dummy_var_connected]%>%max()

                max_sim_value_id <- (model_data_sim[k,dummy_var_connected] == max_sim_value) %>%
                  which()

                model_data_sim_new[k, dummy_var_connected[-max_sim_value_id[1]]] <- 0
              }

            }
          }


        }


      }
    }



    #Erstellt ein lm basierend auf den neu generierten Daten
    lm2 <- lm(paste((model_data_sim_new%>%names())[1],"~ ."), data = model_data_sim_new)

    if (length(formula) == 2) {
      stop(
        "Please use the following structure for lm(): data = data, y ~ x1 + xn..; Execution stopped."
      )
    }

    par(mfrow = c(1, 2))
    par(mar = c(2, 2, 2, 2))


    car::residualPlot(lm)
    mtext("Residuals vs fitted,\n echtes Modell")
    car::residualPlot(lm2)
    mtext("Residuals vs fitted,\n simuliertes Modell")

    if(!rMarkdown){linearity_plot_result <-
      readline("Assume Linearity based on left plot? (y/n):")}

    car::qqPlot(lm)
    mtext("QQ-Plot,\n echtes Modell")
    car:::qqPlot(lm2$residuals)
    mtext("QQ-Plot,\n simuliertes Modell")

    if(!rMarkdown){normality_plot_result <-
      readline("Assume Normality based on left plot? (y/n):")}

    plot(lm, 4, caption = "")
    mtext("Cook's d Plot,\n echtes Modell")
    plot(lm2, 4, caption = "")
    mtext("Cook's d Plot,\n simuliertes Modell")

    if(!rMarkdown){outliers_plot_result <-
      readline("Identify Outliers based on left plot. (0 = non):")}





    print(paste("Seed:",seed))
    lm %>% rstandard() %>% shapiro.test() %>% print()
    lmtest::bptest(lm) %>% print()
    lmtest::resettest(lm) %>% print()
    if (preds > 2& interactionYN) {
      print(var_inf)
    }


  }

