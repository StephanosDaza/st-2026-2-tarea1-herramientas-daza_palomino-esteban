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

if(!dir.exists("figs")) dir.create("figs")

ejecutar_protocolo <- function(nombre_serie, serie_raw, fuente, unidad, 
                               metodo_fn, nombre_metodo, p_params = 0, 
                               prefijo, param_rejilla = NULL, es_tendencia = FALSE, 
                               tipo_tend = "lineal", corregir_ses = FALSE) {
  
  cat(sprintf("\n\n>>> %s: %s sobre %s <<<\n", prefijo, nombre_metodo, nombre_serie))
  

  df <- leer_serie(serie_raw, fuente = fuente, unidad = unidad)
  freq <- attr(df, "frecuencia")
  N <- nrow(df)
  cat(sprintf("- Fuente: %s\n- Unidad: %s\n- Frecuencia: %d\n- Rango: %s a %s\n- Obs verificadas (N): %d\n",
              attr(df, "fuente"), attr(df, "unidad"), freq, min(df$fecha), max(df$fecha), N))
  

  p_ini <- graficar_serie(df, paste(nombre_serie, "- Serie Original"))
  p_corr <- correlograma(df)
  ggsave(paste0("figs/", prefijo, "_ini.png"), p_ini / p_corr, width = 8, height = 8, dpi = 150)
  
  acf_ser <- as.numeric(stats::acf(df$y, lag.max = 24, plot = FALSE)$acf)[-1]
  cat("- Prueba Ljung-Box (Serie original, p=0):\n")
  lb_ini <- ljung_box(acf_ser, T = N, m = 24, p = 0)
  print(lb_ini)
  
 
  h <- min(12, floor(0.2 * N))
  if (freq > 1 && h < freq) h <- freq
  cat(sprintf("- Particion extrayendo h = %d observaciones\n", h))
  
  train_y <- df$y[1:(N-h)]
  test_y <- df$y[(N-h+1):N]
  
 
  if (freq > 1) {
    mad_naive <- mean(abs(diff(train_y, lag = freq)), na.rm = TRUE)
  } else {
    mad_naive <- mean(abs(diff(train_y)), na.rm = TRUE)
  }
  

  if (!is.null(param_rejilla)) {
    opt <- optimizar(train_y, metodo_fn, param_rejilla)
    cat("- Parametros optimizados en tramo de estimacion:\n")
    print(opt$optimo)
    
    if (metodo_fn == "ajustar_mm") ajuste <- ajustar_mm(train_y, k = opt$optimo[[1]])
    else if (metodo_fn == "ajustar_ses") ajuste <- ajustar_ses(train_y, alpha = opt$optimo[[1]])
    else if (metodo_fn == "ajustar_dmm") ajuste <- ajustar_dmm(train_y, k = opt$optimo[[1]])
    else if (metodo_fn == "ajustar_holt") ajuste <- ajustar_holt(train_y, alpha = opt$optimo[[1]], beta = opt$optimo[[2]])
  } else {
    if (es_tendencia) ajuste <- ajustar_tendencia(train_y, tipo = tipo_tend, corregir_sesgo = corregir_ses)
    else ajuste <- ajustar_media(train_y)
  }
  
  preds <- ajuste$pronosticar(h)
  med_in <- medidas(train_y, ajuste$yhat, naive_mae = mad_naive)
  med_out <- medidas(test_y, preds, naive_mae = mad_naive)
  
  cat("- Medidas Ajuste (IN):\n"); print(med_in)
  cat("- Medidas Pronostico (OUT) [Referente Ingenuo MAD = ", mad_naive, "]:\n", sep="")
  print(med_out)
  
  if (es_tendencia) {
    cat("- Tabla de Coeficientes OLS:\n")
    print(ajuste$parametros$tabla)
  }
  
  
  errores <- train_y - ajuste$yhat
  cat(sprintf("- Validacion de Errores (N_errores = %d, p = %d para Ljung-Box):\n", length(errores[!is.na(errores)]), p_params))
  val_err <- validar_errores(errores, yhat = ajuste$yhat, p = p_params)
  ggsave(paste0("figs/", prefijo, "_err.png"), val_err$graficos, width = 10, height = 8, dpi = 150)
  print(val_err[c("prueba_t_media", "ljung_box", "jarque_bera", "durbin_watson")])
  
 
  df_plot <- df
  df_plot$Ajustado <- c(ajuste$yhat, rep(NA, h))
  df_plot$Pronostico <- c(rep(NA, N-h), preds)
  
  p_fin <- ggplot(df_plot, aes(x = fecha)) +
    geom_line(aes(y = y, color = "1. Original"), linewidth = 0.8) +
    geom_line(aes(y = Ajustado, color = "2. Ajuste (In)"), linewidth = 0.8) +
    geom_line(aes(y = Pronostico, color = "3. Pronostico (Out)"), linewidth = 1) +
    scale_color_manual(values = c("1. Original" = "black", "2. Ajuste (In)" = "steelblue", "3. Pronostico (Out)" = "darkred")) +
    labs(title = paste(nombre_metodo, "-", nombre_serie), x = "Fecha", y = "Valor", color = "Serie") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(paste0("figs/", prefijo, "_fin.png"), p_fin, width = 8, height = 5, dpi = 150)
}



# CORRIDA DE LOS 8 EJEMPLOS 

# EJEMPLO 1 - CONTRAEJEMPLO - Media Simple asume nivel local constante, viola fuerte estacionalidad y tendencia
ejecutar_protocolo("AirPassengers", datasets::AirPassengers, "Box & Jenkins (1976)", "Miles de pasajeros", 
                   "ajustar_media", "Media Simple (Contraejemplo)", p_params = 0, prefijo = "01")

# EJEMPLO 2 - Media Movil (Nivel constante local, sin tendencia clara)
ejecutar_protocolo("discoveries", datasets::discoveries, "World Almanac", "Numero de descubrimientos", 
                   "ajustar_mm", "Media Movil", p_params = 0, prefijo = "02", param_rejilla = 2:12)

# EJEMPLO 3 - SES (Nivel variante, sin tendencia agresiva)
ejecutar_protocolo("Nile", datasets::Nile, "Durbin & Koopman", "Flujo (10^8 m^3)", 
                   "ajustar_ses", "Suavizamiento Exponencial Simple", p_params = 1, prefijo = "03", param_rejilla = seq(0.02, 0.98, by=0.02))

# EJEMPLO 4 - DMM (Tendencia local variable que rompe media movil simple)
ejecutar_protocolo("WWWusage", datasets::WWWusage, "Durbin & Koopman", "Usuarios por minuto", 
                   "ajustar_dmm", "Doble Media Movil", p_params = 0, prefijo = "04", param_rejilla = 2:12)

# EJEMPLO 5 - Tendencia Lineal (Tendencia recta global persistente)
ejecutar_protocolo("austres", datasets::austres, "Australian Bureau of Statistics", "Miles de residentes", 
                   "ajustar_tendencia", "Tendencia Lineal", p_params = 2, prefijo = "05", es_tendencia = TRUE, tipo_tend = "lineal")

# EJEMPLO 6 - Tendencia Cuadratica (Aceleracion visible en la historia)
ejecutar_protocolo("airmiles", datasets::airmiles, "Box & Jenkins", "Millas (millones)", 
                   "ajustar_tendencia", "Tendencia Cuadratica", p_params = 3, prefijo = "06", es_tendencia = TRUE, tipo_tend = "cuadratica")

# EJEMPLO 7 - Tendencia Exponencial (Varianza aumenta con el nivel)
ejecutar_protocolo("JohnsonJohnson", datasets::JohnsonJohnson, "Shumway & Stoffer", "Ganancia por accion ($)", 
                   "ajustar_tendencia", "Tendencia Exponencial", p_params = 2, prefijo = "07", es_tendencia = TRUE, tipo_tend = "exponencial", corregir_ses = TRUE)

# EJEMPLO 8 - Holt Lineal (Tendencia estocastica adaptable)
ejecutar_protocolo("LakeHuron", datasets::LakeHuron, "Brockwell & Davis", "Nivel en pies", 
                   "ajustar_holt", "Holt Lineal", p_params = 2, prefijo = "08", param_rejilla = expand.grid(seq(0.05, 0.95, by=0.05), seq(0.05, 0.95, by=0.05)))



# BLOQUES DE VERIFICACION

cat("\n\n====================================================================\n")
cat("BLOQUES DE VERIFICACION EXIGIDOS (LITERALES 2 y 3)\n")
cat("====================================================================\n")

set.seed(2026)
y_verif <- rnorm(100)
n_v <- length(y_verif)
m_v <- min(floor(n_v / 4), 24)


y_bar_v <- mean(y_verif)
var_y_v <- sum((y_verif - y_bar_v)^2) / n_v 
acf_man <- numeric(m_v)
for (i in 1:m_v) acf_man[i] <- sum((y_verif[1:(n_v-i)] - y_bar_v) * (y_verif[(1+i):n_v] - y_bar_v)) / n_v / var_y_v

acf_r <- as.numeric(stats::acf(y_verif, lag.max = m_v, plot = FALSE)$acf)[2:(m_v+1)]
cat("- ACF Maxima Diferencia Absoluta (< 1e-12):", max(abs(acf_man - acf_r)), "\n")

lb_man <- ljung_box(acf_man, T = n_v, m = m_v, p = 0)$estadistico
lb_r <- unname(stats::Box.test(y_verif, lag = m_v, type = "Ljung-Box")$statistic)
cat("- Ljung-Box Qm Diferencia Absoluta:", abs(lb_man - lb_r), "\n")

# SES Promedio Ponderado vs Correccion Error
alpha_v <- 0.3
ajuste_ses <- ajustar_ses(y_verif, alpha_v)$yhat
y_pon <- rep(NA, n_v)
y_pon[2] <- y_verif[1]
for(t in 2:(n_v-1)) y_pon[t+1] <- alpha_v * y_verif[t] + (1 - alpha_v) * y_pon[t]
cat("- (3a) SES Coincide (Promedio vs Error):", all.equal(ajuste_ses[!is.na(ajuste_ses)], y_pon[!is.na(y_pon)]), "\n")

# OLS vs lm()
fit_lm <- coef(lm(y_verif ~ seq_along(y_verif)))
fit_man <- ajustar_tendencia(y_verif, "lineal")$parametros$tabla$Estimacion
cat("- (3b) Tendencias Coeficientes OLS vs lm():", all.equal(as.numeric(fit_lm), as.numeric(fit_man)), "\n")

# Holt Ecuaciones
beta_v <- 0.2
ajuste_h <- ajustar_holt(y_verif, alpha_v, beta_v)
L <- ajuste_h$parametros$L_t; T_val <- ajuste_h$parametros$T_t
L_c <- rep(NA, n_v); T_c <- rep(NA, n_v)
L_c[1] <- y_verif[1]; T_c[1] <- 0
e_v <- y_verif - ajuste_h$yhat
for(t in 2:n_v) {
  L_c[t] <- L[t-1] + T_val[t-1] + alpha_v * e_v[t]
  T_c[t] <- T_val[t-1] + alpha_v * beta_v * e_v[t]
}
cat("- (3c) Holt Coincide Nivel:", all.equal(L[!is.na(L)], L_c[!is.na(L_c)]), "\n")
cat("- (3c) Holt Coincide Pendiente:", all.equal(T_val[!is.na(T_val)], T_c[!is.na(T_c)]), "\n")
cat("====================================================================\n")