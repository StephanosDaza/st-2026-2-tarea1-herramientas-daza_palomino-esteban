library(ggplot2)
library(patchwork)

medidas <- function(y, yhat, naive_mae = NULL) {
  idx <- !is.na(y) & !is.na(yhat)
  y_c <- y[idx]
  yhat_c <- yhat[idx]
  
  e <- y_c - yhat_c
  
  mse <- mean(e^2)
  mad <- mean(abs(e))
  mape <- mean(abs(e) / y_c) * 100
  mase <- if (!is.null(naive_mae) && naive_mae > 0) mad / naive_mae else NA
  
  return(c(MSE = mse, MAD = mad, MAPE = mape, MASE = mase))
}

ljung_box <- function(r, T_obs, m, p = 0) {
  if (m <= p) stop("m debe ser mayor que p.")
  
  Qm <- T_obs * (T_obs + 2) * sum((r[1:m]^2) / (T_obs - (1:m)))
  df <- m - p
  cv <- qchisq(0.95, df)
  pval <- 1 - pchisq(Qm, df)
  
  return(list(estadistico = Qm, grados_libertad = df, valor_critico = cv, valor_p = pval))
}

jarque_bera <- function(e) {
  e <- e[!is.na(e)]
  N <- length(e)
  if (N == 0) stop("No hay errores validos.")
  
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
  
  return(list(estadistico = JB, grados_libertad = df, valor_critico = cv, valor_p = pval))
}

durbin_watson <- function(e) {
  e <- e[!is.na(e)]
  N <- length(e)
  if (N < 2) stop("Datos insuficientes para DW.")
  
  diff_e <- e[2:N] - e[1:(N-1)]
  d <- sum(diff_e^2) / sum(e^2)
  
  return(list(estadistico = d))
}

validar_errores <- function(e, m = NULL, p = 0) {
  e_clean <- e[!is.na(e)]
  N <- length(e_clean)
  if (N < 2) stop("Datos insuficientes en validar_errores.")
  
  if (is.null(m)) m <- min(floor(N / 4), 24)
  
  t_stat <- sqrt(N) * mean(e_clean) / sd(e_clean)
  t_cv <- qt(0.975, df = N - 1)
  t_pval <- 2 * (1 - pt(abs(t_stat), df = N - 1))
  
  t_test <- list(estadistico = t_stat, grados_libertad = N - 1, valor_critico = t_cv, valor_p = t_pval)
  
  y_bar <- mean(e_clean)
  var_y <- sum((e_clean - y_bar)^2) / N 
  acf_vals <- numeric(m)
  for (h in 1:m) {
    acf_vals[h] <- sum((e_clean[1:(N-h)] - y_bar) * (e_clean[(1+h):N] - y_bar)) / N / var_y
  }
  
  lb <- ljung_box(acf_vals, N, m, p)
  jb <- jarque_bera(e_clean)
  dw <- durbin_watson(e_clean)
  
  df_e <- data.frame(t = 1:N, e = e_clean)
  p_time <- ggplot(df_e, aes(x = t, y = e)) +
    geom_line(color = "darkred") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Errores en el tiempo", x = "Indice", y = "Error") +
    theme_minimal()
  
  p_acf_pacf <- correlograma(e_clean, m)
  graficos <- p_time / p_acf_pacf + plot_layout(heights = c(1, 2))
  
  return(list(
    n_errores = N,
    prueba_t_media = t_test,
    ljung_box = lb,
    jarque_bera = jb,
    durbin_watson = dw,
    graficos = graficos
  ))
}