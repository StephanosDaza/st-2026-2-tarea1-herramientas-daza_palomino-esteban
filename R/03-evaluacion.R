library(ggplot2)
library(patchwork)

medidas <- function(y, yhat, naive_mae = NULL) {
  stopifnot(
    "y y yhat deben ser numericos" = is.numeric(y) && is.numeric(yhat),
    "y y yhat deben tener la misma longitud" = length(y) == length(yhat)
  )
  idx <- !is.na(y) & !is.na(yhat)
  stopifnot("No hay pares (y, yhat) validos" = any(idx))
  y_c <- y[idx]
  yhat_c <- yhat[idx]
  e <- y_c - yhat_c
  
  mse <- mean(e^2)
  mad <- mean(abs(e))
  mape <- mean(abs(e) / abs(y_c)) * 100
  mase <- if (!is.null(naive_mae) && naive_mae > 0) mad / naive_mae else NA
  
  return(c(MSE = mse, MAD = mad, MAPE = mape, MASE = mase))
}

ljung_box <- function(r, T, m, p = 0) {
  stopifnot(
    "r debe ser numerico" = is.numeric(r),
    "r debe tener al menos m autocorrelaciones" = length(r) >= m,
    "T debe ser mayor que m" = T > m,
    "m debe ser mayor que p" = m > p,
    "p no puede ser negativo" = p >= 0
  )
  Qm <- T * (T + 2) * sum((r[1:m]^2) / (T - (1:m)))
  df <- m - p
  cv <- qchisq(0.95, df)
  pval <- 1 - pchisq(Qm, df)
  
  return(list(estadistico = Qm, m = m, grados_libertad = df, valor_critico = cv, valor_p = pval))
}

jarque_bera <- function(e) {
  stopifnot("e debe ser numerico" = is.numeric(e))
  e <- e[!is.na(e)]
  N <- length(e)
  stopifnot("Se necesitan al menos 3 errores validos" = N >= 3)
  
  mu <- mean(e)
  s2 <- mean((e - mu)^2)
  m3 <- mean((e - mu)^3)
  m4 <- mean((e - mu)^4)
  
  A <- m3 / (s2^(3/2))
  K <- m4 / (s2^2)
  
  JB <- (N / 6) * (A^2 + ((K - 3)^2) / 4)
  df <- 2
  cv <- qchisq(0.95, df)
  pval <- 1 - pchisq(JB, df)
  
  return(list(estadistico = JB, N = N, asimetria = A, curtosis = K,
              grados_libertad = df, valor_critico = cv, valor_p = pval,
              advertencia = if (N < 20) "N < 20: prueba asintotica con poca potencia" else NA))
}

durbin_watson <- function(e, dL = NULL, dU = NULL) {
  stopifnot("e debe ser numerico" = is.numeric(e))
  e <- e[!is.na(e)]
  N <- length(e)
  stopifnot("Se necesitan al menos 2 errores validos" = N >= 2)
  
  diff_e <- e[2:N] - e[1:(N - 1)]
  d <- sum(diff_e^2) / sum(e^2)
  
  decision <- NA
  if (!is.null(dL) && !is.null(dU)) {
    stopifnot("Debe cumplirse 0 < dL < dU < 2" = dL > 0 && dL < dU && dU < 2)
    decision <- if (d < dL) {
      "Rechazar H0: autocorrelacion positiva"
    } else if (d <= dU) {
      "Zona no concluyente (lado positivo)"
    } else if (d < 4 - dU) {
      "No rechazar H0"
    } else if (d <= 4 - dL) {
      "Zona no concluyente (lado negativo)"
    } else {
      "Rechazar H0: autocorrelacion negativa"
    }
  }
  
  return(list(estadistico = d, N = N, r1_aprox = 1 - d / 2,
              grados_libertad = NA, valor_critico = NA, valor_p = NA,
              dL = if (is.null(dL)) NA else dL, dU = if (is.null(dU)) NA else dU,
              decision = decision))
}

validar_errores <- function(e, yhat = NULL, m = NULL, p = 0, dL = NULL, dU = NULL) {
  stopifnot(
    "e debe ser numerico" = is.numeric(e),
    "yhat debe tener la misma longitud que e" = is.null(yhat) || length(yhat) == length(e),
    "p no puede ser negativo" = p >= 0
  )
  e_clean <- e[!is.na(e)]
  N <- length(e_clean)
  stopifnot("Se necesitan al menos 4 errores validos" = N >= 4)
  if (is.null(m)) m <- min(floor(N / 4), 24)
  
  t_stat <- sqrt(N) * mean(e_clean) / sd(e_clean)
  t_test <- list(
    estadistico = t_stat,
    grados_libertad = N - 1,
    valor_critico = qt(0.975, df = N - 1),
    valor_p = 2 * (1 - pt(abs(t_stat), df = N - 1))
  )
  
  acf_vals <- acf_muestral(e_clean, m)
  banda <- banda_ruido_blanco(N)
  
  lb <- ljung_box(acf_vals, T = N, m = m, p = p)
  jb <- jarque_bera(e_clean)
  dw <- durbin_watson(e_clean, dL = dL, dU = dU)
  
  df_e <- data.frame(t = 1:N, e = e_clean)
  p_time <- ggplot(df_e, aes(x = t, y = e)) +
    geom_line(color = "darkred") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Errores en el tiempo", x = "Indice", y = "Error") + theme_minimal()
  
  p_acf_pacf <- correlograma(e_clean, m)
  
  if (!is.null(yhat)) {
    df_res <- data.frame(yhat = yhat[!is.na(e)], e = e_clean)
    p_res <- ggplot(df_res, aes(x = yhat, y = e)) +
      geom_point(color = "blue", alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Residuos vs Ajustados", x = "Valores Ajustados", y = "Residuos") + theme_minimal()
    graficos <- (p_time | p_res) / p_acf_pacf + plot_layout(heights = c(1, 2))
  } else {
    graficos <- p_time / p_acf_pacf + plot_layout(heights = c(1, 2))
  }
  
  return(list(n_errores = N, m = m, acf = acf_vals, banda = banda,
              rezagos_fuera_banda = which(abs(acf_vals) > banda),
              prueba_t_media = t_test, ljung_box = lb, jarque_bera = jb,
              durbin_watson = dw, graficos = graficos))
}
