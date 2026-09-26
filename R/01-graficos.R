library(ggplot2)
library(patchwork)

graficar_serie <- function(datos, titulo) {
  stopifnot(
    "datos debe ser el tibble generado por leer_serie" =
      is.data.frame(datos) && all(c("fecha", "y") %in% names(datos)),
    "titulo debe ser un texto" = is.character(titulo) && length(titulo) == 1
  )
  
  fuente <- attr(datos, "fuente")
  unidad <- attr(datos, "unidad")
  n <- nrow(datos)
  
  fuente_texto <- if (is.null(fuente)) "Desconocida" else fuente
  unidad_texto <- if (is.null(unidad)) "" else unidad
  
  p <- ggplot(datos, aes(x = fecha, y = y)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    labs(
      title = titulo,
      x = "Fecha",
      y = unidad_texto,
      caption = paste0("Fuente: ", fuente_texto, " | Numero de observaciones: ", n)
    ) +
    scale_x_date(date_labels = "%Y-%m") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, face = "italic", color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

acf_muestral <- function(y, m) {
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y),
    "y no puede tener NA" = !anyNA(y),
    "m debe ser un entero entre 1 y n - 1" =
      length(m) == 1 && m == round(m) && m >= 1 && m < length(y)
  )
  n <- length(y)
  y_bar <- mean(y)
  c0 <- sum((y - y_bar)^2) / n
  r <- numeric(m)
  for (h in 1:m) {
    r[h] <- sum((y[1:(n - h)] - y_bar) * (y[(1 + h):n] - y_bar)) / n / c0
  }
  return(r)
}

banda_ruido_blanco <- function(n, ci = 0.95) {
  stopifnot("n debe ser positivo" = n > 0, "ci debe estar en (0, 1)" = ci > 0 && ci < 1)
  qnorm((1 + ci) / 2) / sqrt(n)
}

correlograma <- function(datos, m = NULL) {
  if (is.data.frame(datos) && "y" %in% names(datos)) {
    y <- datos$y
  } else {
    stopifnot("datos debe ser un tibble con columna y o un vector numerico" = is.numeric(datos))
    y <- datos
  }
  
  y_clean <- y[!is.na(y)]
  n <- length(y_clean)
  
  if (is.null(m)) m <- min(floor(n / 4), 24)
  stopifnot("m esta fuera de rango" = m > 0 && m < n)
  
  acf_vals <- acf_muestral(y_clean, m)
  
  pacf_obj <- stats::pacf(y_clean, lag.max = m, plot = FALSE)
  pacf_vals <- as.numeric(pacf_obj$acf)[1:m]
  
  climo <- banda_ruido_blanco(n)
  
  df_acf <- data.frame(lag = 1:m, acf = acf_vals)
  df_pacf <- data.frame(lag = 1:m, pacf = pacf_vals)
  
  p_acf <- ggplot(df_acf, aes(x = lag, y = acf)) +
    geom_segment(aes(x = lag, xend = lag, y = 0, yend = acf), color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(climo, -climo), color = "blue", linetype = "dashed") +
    labs(title = "ACF", x = "", y = "ACF") +
    scale_x_continuous(breaks = 1:m) +
    theme_bw()
  
  p_pacf <- ggplot(df_pacf, aes(x = lag, y = pacf)) +
    geom_segment(aes(x = lag, xend = lag, y = 0, yend = pacf), color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(climo, -climo), color = "blue", linetype = "dashed") +
    labs(title = "PACF", x = "Rezago", y = "PACF") +
    scale_x_continuous(breaks = 1:m) +
    theme_bw()
  
  return(p_acf / p_pacf)
}

graficar_optimizacion <- function(opt, titulo) {
  stopifnot(
    "opt debe ser la salida de optimizar()" =
      is.list(opt) && all(c("rejilla", "optimo") %in% names(opt))
  )
  rej <- opt$rejilla
  op <- opt$optimo
  pars <- setdiff(names(rej), "mse")
  
  if (length(pars) == 1) {
    df <- data.frame(x = rej[[pars[1]]], mse = rej$mse)
    df_op <- data.frame(x = op[[pars[1]]], mse = op$mse)
    p <- ggplot(df, aes(x = x, y = mse)) +
      geom_line(color = "steelblue") +
      geom_point(color = "steelblue", size = 1.2) +
      geom_vline(xintercept = df_op$x, linetype = "dashed", color = "darkred") +
      geom_point(data = df_op, color = "darkred", size = 3) +
      labs(title = titulo,
           subtitle = sprintf("Optimo: %s = %s, MSE = %.4f", pars[1], format(df_op$x), df_op$mse),
           x = pars[1], y = "MSE de un paso (tramo de estimacion)") +
      theme_minimal()
  } else {
    df <- data.frame(a = rej[[pars[1]]], b = rej[[pars[2]]], mse = rej$mse)
    df_op <- data.frame(a = op[[pars[1]]], b = op[[pars[2]]])
    p <- ggplot(df, aes(x = a, y = b, fill = mse)) +
      geom_tile() +
      scale_fill_viridis_c(name = "MSE") +
      geom_point(data = df_op, aes(x = a, y = b), inherit.aes = FALSE,
                 shape = 4, size = 5, stroke = 1.5, color = "red") +
      labs(title = titulo,
           subtitle = sprintf("Optimo: %s = %.2f, %s = %.2f, MSE = %.4f",
                              pars[1], df_op$a, pars[2], df_op$b, op$mse),
           x = pars[1], y = pars[2]) +
      theme_minimal()
  }
  return(p)
}
