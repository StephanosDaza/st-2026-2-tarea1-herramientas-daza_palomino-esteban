library(ggplot2)
library(patchwork)

graficar_serie <- function(datos, titulo) {
  if (!is.data.frame(datos) || !all(c("fecha", "y") %in% names(datos))) {
    stop("datos debe ser el tibble generado por leer_serie.")
  }
  
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
    scale_x_date(date_labels = "%Y-%m", date_breaks = "auto") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, face = "italic", color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

correlograma <- function(datos, m = NULL) {
  if (is.data.frame(datos) && "y" %in% names(datos)) {
    y <- datos$y
  } else if (is.numeric(datos)) {
    y <- datos
  } else {
    stop("datos debe ser un tibble o un vector numerico.")
  }
  
  y_clean <- y[!is.na(y)]
  n <- length(y_clean)
  
  if (is.null(m)) m <- min(floor(n / 4), 24)
  if (m <= 0 || m >= n) stop("m esta fuera de rango.")
  
  y_bar <- mean(y_clean)
  var_y <- sum((y_clean - y_bar)^2) / n 
  
  acf_vals <- numeric(m)
  for (h in 1:m) {
    cov_h <- sum((y_clean[1:(n-h)] - y_bar) * (y_clean[(1+h):n] - y_bar)) / n
    acf_vals[h] <- cov_h / var_y
  }
  
  pacf_obj <- stats::pacf(y_clean, lag.max = m, plot = FALSE)
  pacf_vals <- as.numeric(pacf_obj$acf)[1:m]
  
  climo <- qnorm((1 + 0.95) / 2) / sqrt(n)
  
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