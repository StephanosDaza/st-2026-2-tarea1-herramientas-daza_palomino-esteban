library(tibble)

leer_serie <- function(x, fuente, unidad) {
  if (inherits(x, "ts")) {
    n <- length(x)
    freq <- frequency(x)
    inicio <- start(x)
    anio_ini <- inicio[1]
    ciclo_ini <- inicio[2]
    
    if (freq == 1) {
      fechas <- seq(as.Date(paste0(anio_ini, "-01-01")), by = "year", length.out = n)
    } else if (freq == 4) {
      mes_ini <- (ciclo_ini - 1) * 3 + 1
      fechas <- seq(as.Date(paste0(anio_ini, "-", mes_ini, "-01")), by = "quarter", length.out = n)
    } else if (freq == 12) {
      fechas <- seq(as.Date(paste0(anio_ini, "-", ciclo_ini, "-01")), by = "month", length.out = n)
    } else {
      stop("Frecuencia no soportada para conversión de ts a Date.")
    }
    
    df <- tibble(t = 1:n, fecha = fechas, y = as.numeric(x))
    
  } else if (is.character(x)) {
    if (!file.exists(x)) stop("El archivo CSV no existe en la ruta.")
    
    datos <- read.csv(x)
    if (!all(c("fecha", "valor") %in% colnames(datos))) {
      stop("El CSV debe tener exactamente las columnas 'fecha' y 'valor'.")
    }
    
    fechas <- as.Date(datos$fecha)
    n <- nrow(datos)
    
    df <- tibble(t = 1:n, fecha = fechas, y = as.numeric(datos$valor))
    
    dias_diff <- as.numeric(diff(fechas))
    mediana_dias <- median(dias_diff, na.rm = TRUE)
    
    if (mediana_dias >= 360) {
      freq <- 1
    } else if (mediana_dias >= 85 && mediana_dias <= 95) {
      freq <- 4
    } else if (mediana_dias >= 28 && mediana_dias <= 31) {
      freq <- 12
    } else if (mediana_dias == 1) {
      freq <- 365
    } else {
      stop("No se pudo inferir la frecuencia desde el CSV.")
    }
  } else {
    stop("El objeto x debe ser ts o una ruta a un csv.")
  }
  
  if (any(diff(df$fecha) <= 0)) {
    stop("Fechas no estrictamente crecientes.")
  }
  
  paso <- if (freq == 1) "year" else if (freq == 4) "quarter" else if (freq == 12) "month" else "day"
  fecha_esperada <- seq(df$fecha[1], by = paso, length.out = n)
  
  if (!all(df$fecha == fecha_esperada)) {
    stop("Fechas sin equiespaciamiento correcto.")
  }
  
  attr(df, "frecuencia") <- freq
  attr(df, "fuente") <- fuente
  attr(df, "unidad") <- unidad
  
  return(df)
}