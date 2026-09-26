library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)

source("R/00-lectura.R")
source("R/01-graficos.R")
source("R/02-metodos.R")
source("R/03-evaluacion.R")

if (!dir.exists("figs")) dir.create("figs")

ejecutar_protocolo <- function(nombre_serie, serie_raw, fuente, unidad,
                               metodo_fn, nombre_metodo, p_params = 0,
                               prefijo, param_rejilla = NULL, es_tendencia = FALSE,
                               tipo_tend = "lineal", corregir_ses = FALSE,
                               estacional = NULL, dL = NULL, dU = NULL) {
  
  cat(sprintf("\n\n>>> %s: %s sobre %s <<<\n", prefijo, nombre_metodo, nombre_serie))
  
  df <- leer_serie(serie_raw, fuente = fuente, unidad = unidad)
  freq <- attr(df, "frecuencia")
  N <- nrow(df)
  if (is.null(estacional)) estacional <- freq > 1
  cat(sprintf("- Fuente: %s\n- Unidad: %s\n- Frecuencia: %d\n- Rango: %s a %s\n- Obs verificadas (N): %d\n",
              attr(df, "fuente"), attr(df, "unidad"), freq, min(df$fecha), max(df$fecha), N))
  
  p_ini <- graficar_serie(df, paste(nombre_serie, "- Serie Original"))
  p_corr <- correlograma(df)
  ggsave(paste0("figs/", prefijo, "_ini.png"), p_ini / p_corr, width = 8, height = 8, dpi = 150)
  
  m_ser <- min(floor(N / 4), 24)
  acf_ser <- acf_muestral(df$y, m_ser)
  banda_ser <- banda_ruido_blanco(N)
  fuera_ser <- which(abs(acf_ser) > banda_ser)
  cat(sprintf("- Banda de ruido blanco: +/- %.4f. Rezagos fuera de la banda: %s\n",
              banda_ser, if (length(fuera_ser) == 0) "ninguno" else paste(fuera_ser, collapse = ", ")))
  cat(sprintf("- Prueba Ljung-Box (Serie original, m = %d, p = 0):\n", m_ser))
  lb_ini <- ljung_box(acf_ser, T = N, m = m_ser, p = 0)
  print(lb_ini)
  
  h <- min(12, floor(0.2 * N))
  if (estacional && h < freq) h <- freq
  cat(sprintf("- Particion extrayendo h = %d observaciones\n", h))
  
  train_y <- df$y[1:(N - h)]
  test_y <- df$y[(N - h + 1):N]
  n_train <- length(train_y)
  
  if (estacional) {
    mad_naive <- mean(abs(diff(train_y, lag = freq)))
    pred_naive <- train_y[n_train - freq + ((seq_len(h) - 1) %% freq) + 1]
    nombre_ref <- "Ingenuo estacional"
  } else {
    mad_naive <- mean(abs(diff(train_y)))
    pred_naive <- rep(train_y[n_train], h)
    nombre_ref <- "Ingenuo"
  }
  
  opt <- NULL
  if (!is.null(param_rejilla)) {
    opt <- optimizar(train_y, metodo_fn, param_rejilla)
    cat("- Parametros optimizados en tramo de estimacion:\n")
    print(opt$optimo)
    p_opt <- graficar_optimizacion(opt, paste("MSE de un paso -", nombre_metodo, "-", nombre_serie))
    ggsave(paste0("figs/", prefijo, "_opt.png"), p_opt, width = 8, height = 5, dpi = 150)
    
    ajuste <- switch(metodo_fn,
                     ajustar_mm   = ajustar_mm(train_y, k = opt$optimo[[1]]),
                     ajustar_ses  = ajustar_ses(train_y, alpha = opt$optimo[[1]]),
                     ajustar_dmm  = ajustar_dmm(train_y, k = opt$optimo[[1]]),
                     ajustar_holt = ajustar_holt(train_y, alpha = opt$optimo[[1]], beta = opt$optimo[[2]])
    )
  } else {
    if (es_tendencia) ajuste <- ajustar_tendencia(train_y, tipo = tipo_tend, corregir_sesgo = corregir_ses)
    else ajuste <- ajustar_media(train_y)
  }
  
  preds <- ajuste$pronosticar(h)
  med_in <- medidas(train_y, ajuste$yhat, naive_mae = mad_naive)
  med_out <- medidas(test_y, preds, naive_mae = mad_naive)
  med_naive_out <- medidas(test_y, pred_naive, naive_mae = mad_naive)
  
  cat("- Medidas Ajuste (IN):\n"); print(med_in)
  cat("- Medidas Pronostico (OUT) [MAD del ", nombre_ref, " en estimacion = ", mad_naive, "]:\n", sep = "")
  print(med_out)
  cat("- Medidas del ", nombre_ref, " (OUT):\n", sep = "")
  print(med_naive_out)
  
  if (es_tendencia) {
    cat("- Tabla de Coeficientes OLS:\n")
    print(ajuste$parametros$tabla)
    cat(sprintf("- R2 = %.4f, sigma2 = %.6f, DW = %.4f, rezagos HAC = %d\n",
                ajuste$parametros$R2, ajuste$parametros$sigma2,
                ajuste$parametros$DW, ajuste$parametros$rezagos_HAC))
  }
  
  usa_cotas <- es_tendencia && tipo_tend != "exponencial"
  errores <- train_y - ajuste$yhat
  cat(sprintf("- Validacion de Errores (N_errores = %d, p = %d para Ljung-Box):\n",
              sum(!is.na(errores)), p_params))
  val_err <- validar_errores(errores, yhat = ajuste$yhat, p = p_params,
                             dL = if (usa_cotas) dL else NULL,
                             dU = if (usa_cotas) dU else NULL)
  ggsave(paste0("figs/", prefijo, "_err.png"), val_err$graficos, width = 10, height = 8, dpi = 150)
  print(val_err[c("rezagos_fuera_banda", "prueba_t_media", "ljung_box", "jarque_bera", "durbin_watson")])
  
  val_reg <- NULL
  if (es_tendencia && tipo_tend == "exponencial") {
    cat("- Supuestos de la regresion sobre los residuos de ln(Y):\n")
    val_reg <- validar_errores(ajuste$parametros$residuos, yhat = ajuste$parametros$ajustados,
                               p = p_params, dL = dL, dU = dU)
    ggsave(paste0("figs/", prefijo, "_reg.png"), val_reg$graficos, width = 10, height = 8, dpi = 150)
    print(val_reg[c("rezagos_fuera_banda", "prueba_t_media", "ljung_box", "jarque_bera", "durbin_watson")])
  }
  
  df_plot <- df
  df_plot$Ajustado <- c(ajuste$yhat, rep(NA, h))
  df_plot$Pronostico <- c(rep(NA, N - h), preds)
  df_plot$Referente <- c(rep(NA, N - h), pred_naive)
  
  p_fin <- ggplot(df_plot, aes(x = fecha)) +
    geom_line(aes(y = y, color = "1. Original"), linewidth = 0.8) +
    geom_line(aes(y = Ajustado, color = "2. Ajuste (In)"), linewidth = 0.8) +
    geom_line(aes(y = Pronostico, color = "3. Pronostico (Out)"), linewidth = 1) +
    geom_line(aes(y = Referente, color = "4. Referente ingenuo"), linewidth = 0.8, linetype = "dotted") +
    scale_color_manual(values = c("1. Original" = "black", "2. Ajuste (In)" = "steelblue",
                                  "3. Pronostico (Out)" = "darkred", "4. Referente ingenuo" = "gray40")) +
    labs(title = paste(nombre_metodo, "-", nombre_serie), x = "Fecha",
         y = attr(df, "unidad"), color = "Serie") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(paste0("figs/", prefijo, "_fin.png"), p_fin, width = 8, height = 5, dpi = 150)
  
  val_err$graficos <- NULL
  if (!is.null(val_reg)) val_reg$graficos <- NULL
  
  invisible(list(
    serie = nombre_serie, metodo = nombre_metodo, prefijo = prefijo,
    fuente = attr(df, "fuente"), unidad = attr(df, "unidad"), frecuencia = freq,
    fecha_ini = min(df$fecha), fecha_fin = max(df$fecha), N = N,
    m_serie = m_ser, acf_serie = acf_ser, banda_serie = banda_ser, rezagos_fuera_serie = fuera_ser,
    lb_serie = lb_ini, h = h, n_train = n_train, estacional = estacional, referente = nombre_ref,
    optimizacion = opt, parametros = ajuste$parametros[setdiff(names(ajuste$parametros),
                                                               c("residuos", "ajustados", "E_t", "beta1_t", "MM", "DMM", "L_t", "T_t"))],
    medidas_in = med_in, medidas_out = med_out, medidas_referente_out = med_naive_out,
    mad_referente_in = mad_naive, pronosticos = preds, pronostico_referente = pred_naive,
    validacion = test_y, p_params = p_params,
    errores = val_err, errores_regresion = val_reg
  ))
}

# CORRIDA DE LOS 8 EJEMPLOS Y EL CONTRAEJEMPLO

resultados <- list()

# EJEMPLO 1 - Media Simple (nivel constante: la serie no muestra tendencia ni estacionalidad)
resultados$ej01 <- ejecutar_protocolo("discoveries", datasets::discoveries,
                                      "The World Almanac and Book of Facts, 1975 Edition, pages 315-318",
                                      "Numero de descubrimientos", "ajustar_media", "Media Simple",
                                      p_params = 0, prefijo = "01")

# EJEMPLO 2 - Media Movil (nivel localmente constante, sin tendencia clara)
resultados$ej02 <- ejecutar_protocolo("discoveries", datasets::discoveries,
                                      "The World Almanac and Book of Facts, 1975 Edition, pages 315-318",
                                      "Numero de descubrimientos", "ajustar_mm", "Media Movil",
                                      p_params = 1, prefijo = "02", param_rejilla = 2:12)

# EJEMPLO 3 - SES (nivel variante, sin tendencia agresiva)
resultados$ej03 <- ejecutar_protocolo("Nile", datasets::Nile,
                                      "Durbin, J. and Koopman, S. J. (2001) Time Series Analysis by State Space Methods. Oxford University Press",
                                      "Flujo anual (10^8 m^3)", "ajustar_ses", "Suavizamiento Exponencial Simple",
                                      p_params = 1, prefijo = "03", param_rejilla = seq(0.02, 0.98, by = 0.02))

# EJEMPLO 4 - DMM (tendencia local variable)
resultados$ej04 <- ejecutar_protocolo("WWWusage", datasets::WWWusage,
                                      "Durbin, J. and Koopman, S. J. (2001) Time Series Analysis by State Space Methods. Oxford University Press",
                                      "Usuarios conectados por minuto", "ajustar_dmm", "Doble Media Movil",
                                      p_params = 1, prefijo = "04", param_rejilla = 2:12)

# EJEMPLO 5 - Tendencia Lineal (tendencia recta global persistente; sin estacionalidad)
# DW 5% Savin y White (1977): T = 77, k' = 1 (fila T = 75)
resultados$ej05 <- ejecutar_protocolo("austres", datasets::austres,
                                      "P. J. Brockwell and R. A. Davis (1996) Introduction to Time Series and Forecasting. Springer",
                                      "Miles de residentes", "ajustar_tendencia", "Tendencia Lineal",
                                      p_params = 2, prefijo = "05", es_tendencia = TRUE, tipo_tend = "lineal",
                                      estacional = FALSE, dL = 1.598, dU = 1.652)

# EJEMPLO 6 - Tendencia Cuadratica (aceleracion visible en la historia)
# DW 5% Savin y White (1977): T = 20, k' = 2
resultados$ej06 <- ejecutar_protocolo("airmiles", datasets::airmiles,
                                      "F.A.A. Statistical Handbook of Aviation",
                                      "Millas-pasajero (millones)", "ajustar_tendencia", "Tendencia Cuadratica",
                                      p_params = 3, prefijo = "06", es_tendencia = TRUE, tipo_tend = "cuadratica",
                                      dL = 1.100, dU = 1.537)

# EJEMPLO 7 - Tendencia Exponencial (la varianza aumenta con el nivel)
# DW 5% Savin y White (1977): T = 72, k' = 1 (fila T = 70)
resultados$ej07 <- ejecutar_protocolo("JohnsonJohnson", datasets::JohnsonJohnson,
                                      "Shumway, R. H. and Stoffer, D. S. (2000) Time Series Analysis and its Applications. Second Edition. Springer. Example 1.1",
                                      "Ganancia por accion (dolares)", "ajustar_tendencia", "Tendencia Exponencial",
                                      p_params = 2, prefijo = "07", es_tendencia = TRUE, tipo_tend = "exponencial",
                                      corregir_ses = TRUE, dL = 1.583, dU = 1.641)

# EJEMPLO 8 - Holt Lineal (tendencia estocastica adaptable)
resultados$ej08 <- ejecutar_protocolo("LakeHuron", datasets::LakeHuron,
                                      "Brockwell, P. J. and Davis, R. A. (1991). Time Series and Forecasting Methods. Second edition. Springer, New York. Series A, page 555",
                                      "Nivel del lago (pies)", "ajustar_holt", "Holt Lineal",
                                      p_params = 2, prefijo = "08",
                                      param_rejilla = expand.grid(alpha = seq(0.05, 0.95, by = 0.05),
                                                                  beta  = seq(0.05, 0.95, by = 0.05)))

# EJEMPLO 9 - CONTRAEJEMPLO - Media Simple sobre serie con tendencia y estacionalidad
resultados$ej09 <- ejecutar_protocolo("AirPassengers", datasets::AirPassengers,
                                      "Box, G. E. P., Jenkins, G. M. and Reinsel, G. C. (1976) Time Series Analysis, Forecasting and Control. Third Edition. Holden-Day. Series G",
                                      "Miles de pasajeros", "ajustar_media", "Media Simple (Contraejemplo)",
                                      p_params = 0, prefijo = "09")

# BLOQUES DE VERIFICACION

cat("\n\n====================================================================\n")
cat("BLOQUES DE VERIFICACION\n")
cat("====================================================================\n")

set.seed(2026)
y_verif <- rnorm(100)
n_v <- length(y_verif)
m_v <- min(floor(n_v / 4), 24)

# ACF propia contra acf()
acf_man <- acf_muestral(y_verif, m_v)
acf_r <- as.numeric(stats::acf(y_verif, lag.max = m_v, plot = FALSE)$acf)[2:(m_v + 1)]
dif_acf <- max(abs(acf_man - acf_r))
cat("- ACF maxima diferencia absoluta:", dif_acf, "| Menor que 1e-12:", dif_acf < 1e-12, "\n")

# Ljung-Box propio contra Box.test()
lb_man <- ljung_box(acf_man, T = n_v, m = m_v, p = 0)$estadistico
lb_r <- unname(stats::Box.test(y_verif, lag = m_v, type = "Ljung-Box")$statistic)
cat("- Ljung-Box Qm diferencia absoluta:", abs(lb_man - lb_r), "\n")

# SES - forma de correccion de error contra forma de promedio ponderado
alpha_v <- 0.3
ajuste_ses <- ajustar_ses(y_verif, alpha_v)$yhat
y_pon <- rep(NA_real_, n_v)
for (t in 1:(n_v - 1)) {
  j <- 0:(t - 1)
  pesos <- alpha_v * (1 - alpha_v)^j
  pesos[t] <- (1 - alpha_v)^(t - 1)
  y_pon[t + 1] <- sum(pesos * y_verif[t - j])
}
cat("- SES coincide (correccion de error vs promedio ponderado):",
    isTRUE(all.equal(ajuste_ses, y_pon)), "\n")

# Coeficientes MCO propios contra lm() en las tres tendencias
tt_v <- seq_along(y_verif)
y_pos <- exp(0.02 * tt_v + 0.1 * y_verif)
coef_lin <- ajustar_tendencia(y_verif, "lineal")$parametros$tabla$Estimacion
coef_cua <- ajustar_tendencia(y_verif, "cuadratica")$parametros$tabla$Estimacion
coef_exp <- ajustar_tendencia(y_pos, "exponencial")$parametros$tabla$Estimacion
cat("- Tendencia lineal vs lm():",
    isTRUE(all.equal(unname(coef(lm(y_verif ~ tt_v))), coef_lin)), "\n")
cat("- Tendencia cuadratica vs lm():",
    isTRUE(all.equal(unname(coef(lm(y_verif ~ tt_v + I(tt_v^2)))), coef_cua)), "\n")
cat("- Tendencia exponencial (sobre ln Y) vs lm():",
    isTRUE(all.equal(unname(coef(lm(log(y_pos) ~ tt_v))), coef_exp)), "\n")

# Holt - forma de correccion de error
beta_v <- 0.2
ajuste_h <- ajustar_holt(y_verif, alpha_v, beta_v)
L <- ajuste_h$parametros$L_t
T_val <- ajuste_h$parametros$T_t
L_c <- rep(NA_real_, n_v); T_c <- rep(NA_real_, n_v)
L_c[1] <- y_verif[1]; T_c[1] <- 0
e_v <- y_verif - ajuste_h$yhat
for (t in 2:n_v) {
  L_c[t] <- L[t - 1] + T_val[t - 1] + alpha_v * e_v[t]
  T_c[t] <- T_val[t - 1] + alpha_v * beta_v * e_v[t]
}
cat("- Holt coincide nivel:", isTRUE(all.equal(L, L_c)), "\n")
cat("- Holt coincide pendiente:", isTRUE(all.equal(T_val, T_c)), "\n")
cat("====================================================================\n")

writeLines(capture.output(sessionInfo()), "sesion-info.txt")
