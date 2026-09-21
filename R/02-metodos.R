ajustar_media <- function(y) {
  stopifnot(is.numeric(y), !any(is.na(y)))
  n <- length(y)
  yhat <- rep(NA, n)
  
  if (n > 1) {
    yhat[2] <- y[1]
    suma <- y[1]
    if (n > 2) {
      for (t in 2:(n - 1)) {
        suma <- suma + y[t]
        yhat[t + 1] <- suma / t
      }
    }
  }
  
  estado_final <- sum(y) / n
  
  pronosticar <- function(h) {
    stopifnot(h > 0, h == round(h))
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(estado_final = estado_final)))
}

ajustar_mm <- function(y, k) {
  stopifnot(is.numeric(y), !any(is.na(y)), k >= 2, k <= length(y), k == round(k))
  n <- length(y)
  yhat <- rep(NA, n)
  
  if (n >= k + 1) {
    suma_ventana <- sum(y[1:k])
    yhat[k + 1] <- suma_ventana / k
    if (n > k + 1) {
      for (t in (k + 1):(n - 1)) {
        suma_ventana <- suma_ventana - y[t - k + 1] + y[t]
        yhat[t + 1] <- suma_ventana / k
      }
    }
  }
  
  estado_final <- mean(y[(n - k + 1):n])
  
  pronosticar <- function(h) {
    stopifnot(h > 0, h == round(h))
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(k = k, estado_final = estado_final)))
}

ajustar_ses <- function(y, alpha) {
  stopifnot(is.numeric(y), !any(is.na(y)), alpha > 0, alpha < 1)
  n <- length(y)
  yhat <- rep(NA, n)
  
  if (n > 1) {
    yhat[2] <- y[1]
    if (n > 2) {
      for (t in 2:(n - 1)) {
        e_t <- y[t] - yhat[t]
        yhat[t + 1] <- yhat[t] + alpha * e_t
      }
    }
  }
  
  e_n <- y[n] - yhat[n]
  estado_final <- yhat[n] + alpha * e_n
  
  pronosticar <- function(h) {
    stopifnot(h > 0, h == round(h))
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(alpha = alpha, estado_final = estado_final)))
}

ajustar_dmm <- function(y, k) {
  stopifnot(is.numeric(y), !any(is.na(y)), k >= 2, k <= length(y) / 2, k == round(k))
  n <- length(y)
  yhat <- rep(NA, n)
  E <- rep(NA, n)
  B <- rep(NA, n)
  mm1 <- rep(NA, n)
  mm2 <- rep(NA, n)
  
  suma1 <- sum(y[1:k])
  mm1[k] <- suma1 / k
  for (t in (k + 1):n) {
    suma1 <- suma1 - y[t - k] + y[t]
    mm1[t] <- suma1 / k
  }
  
  if (2 * k - 1 <= n) {
    suma2 <- sum(mm1[k:(2 * k - 1)])
    mm2[2 * k - 1] <- suma2 / k
    if (2 * k <= n) {
      for (t in (2 * k):n) {
        suma2 <- suma2 - mm1[t - k] + mm1[t]
        mm2[t] <- suma2 / k
      }
    }
  }
  
  for (t in (2 * k - 1):n) {
    E[t] <- 2 * mm1[t] - mm2[t]
    B[t] <- (2 / (k - 1)) * (mm1[t] - mm2[t])
    if (t < n) {
      yhat[t + 1] <- E[t] + B[t]
    }
  }
  
  E_n <- E[n]
  B_n <- B[n]
  
  pronosticar <- function(h) {
    stopifnot(h > 0, h == round(h))
    return(E_n + B_n * (1:h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(k = k, E_t = E, B_t = B)))
}

optimizar <- function(y, metodo, rejilla) {
  stopifnot(is.numeric(y), !any(is.na(y)), is.character(metodo))
  
  if (is.numeric(rejilla) && is.vector(rejilla)) {
    rejilla <- data.frame(param1 = rejilla)
  }
  
  mse_vals <- numeric(nrow(rejilla))
  
  for (i in 1:nrow(rejilla)) {
    params <- as.list(rejilla[i, , drop = FALSE])
    
    if (metodo == "ajustar_mm") {
      ajuste <- ajustar_mm(y, k = params[[1]])
    } else if (metodo == "ajustar_ses") {
      ajuste <- ajustar_ses(y, alpha = params[[1]])
    } else if (metodo == "ajustar_dmm") {
      ajuste <- ajustar_dmm(y, k = params[[1]])
    } else if (metodo == "ajustar_holt") {
      ajuste <- ajustar_holt(y, alpha = params[[1]], beta = params[[2]])
    } else {
      stop("Metodo no soportado.")
    }
    
    e <- y - ajuste$yhat
    mse_vals[i] <- mean(e^2, na.rm = TRUE)
  }
  
  rejilla$mse <- mse_vals
  optimo <- rejilla[which.min(rejilla$mse), ]
  
  return(list(rejilla = rejilla, optimo = optimo))
}




