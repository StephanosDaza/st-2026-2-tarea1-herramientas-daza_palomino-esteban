ajustar_media <- function(y) {
  if (any(is.na(y))) stop("El vector no puede tener NAs.")
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
    if (h <= 0) stop("h debe ser positivo")
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(estado_final = estado_final)))
}

ajustar_mm <- function(y, k) {
  if (any(is.na(y))) stop("El vector no puede tener NAs.")
  n <- length(y)
  if (k > n) stop("La ventana k supera la longitud de la serie.")
  if (k < 2) stop("k debe ser >= 2.")
  
  yhat <- rep(NA, n)
  suma_ventana <- 0
  
  for (t in 1:(n - 1)) {
    suma_ventana <- suma_ventana + y[t]
    if (t > k) suma_ventana <- suma_ventana - y[t - k]
    if (t >= k) {
      yhat[t + 1] <- suma_ventana / k
    }
  }
  
  estado_final <- mean(y[(n - k + 1):n])
  
  pronosticar <- function(h) {
    if (h <= 0) stop("h debe ser positivo")
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(k = k, estado_final = estado_final)))
}

ajustar_ses <- function(y, alpha) {
  if (any(is.na(y))) stop("El vector no puede tener NAs.")
  if (alpha <= 0 || alpha >= 1) stop("alpha debe estar en (0, 1).")
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
    if (h <= 0) stop("h debe ser positivo")
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(alpha = alpha, estado_final = estado_final)))
}

ajustar_dmm <- function(y, k) {
  if (any(is.na(y))) stop("El vector no puede tener NAs.")
  n <- length(y)
  if (k > n / 2) stop("Ventana supera el limite de la doble media movil.")
  if (k < 2) stop("k debe ser >= 2.")
  
  yhat <- rep(NA, n)
  E_t <- rep(NA, n)
  beta1_t <- rep(NA, n)
  mm1 <- rep(NA, n)
  mm2 <- rep(NA, n)
  
  suma1 <- 0
  suma2 <- 0
  
  for (t in 1:n) {
    suma1 <- suma1 + y[t]
    if (t > k) suma1 <- suma1 - y[t - k]
    
    if (t >= k) {
      mm1[t] <- suma1 / k
      suma2 <- suma2 + mm1[t]
      if (t > 2 * k - 1) suma2 <- suma2 - mm1[t - k]
      
      if (t >= 2 * k - 1) {
        mm2[t] <- suma2 / k
        E_t[t] <- 2 * mm1[t] - mm2[t]
        beta1_t[t] <- (2 / (k - 1)) * (mm1[t] - mm2[t])
        if (t < n) yhat[t + 1] <- E_t[t] + beta1_t[t]
      }
    }
  }
  
  E_n <- E_t[n]
  B_n <- beta1_t[n]
  
  pronosticar <- function(h) {
    if (h <= 0) stop("h debe ser positivo")
    return(E_n + B_n * (1:h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(k = k, E_t = E_t, beta1_t = beta1_t)))
}

ajustar_tendencia <- function(y, tipo, corregir_sesgo = FALSE) {
  if (any(is.na(y))) stop("El vector no puede tener NAs.")
  if (!(tipo %in% c("lineal", "cuadratica", "exponencial"))) stop("Tipo de tendencia no valido.")
  n <- length(y)
  
  if (tipo == "exponencial" && any(y <= 0)) stop("Serie con valores no positivos en tendencia exponencial.")
  
  y_adj <- if (tipo == "exponencial") log(y) else y
  
  X <- if (tipo == "lineal") {
    cbind(1, 1:n)
  } else if (tipo == "cuadratica") {
    cbind(1, 1:n, (1:n)^2)
  } else {
    cbind(1, 1:n)
  }
  
  beta_hat <- solve(crossprod(X), crossprod(X, y_adj))
  
  y_fit <- as.numeric(X %*% beta_hat)
  e <- y_adj - y_fit
  
  p <- ncol(X)
  sigma2 <- sum(e^2) / (n - p)
  
  Q_inv <- solve(crossprod(X))
  var_ols <- as.numeric(diag(sigma2 * Q_inv))
  se_ols <- sqrt(var_ols)
  
  L_hac <- floor(4 * (n / 100)^(2/9))
  eX <- X * e
  S <- crossprod(eX)
  if (L_hac > 0 && n > L_hac) {
    for (h in 1:L_hac) {
      Gamma_h <- crossprod(eX[-(1:h), , drop = FALSE], eX[-((n-h+1):n), , drop = FALSE])
      w_h <- 1 - h / (L_hac + 1)
      S <- S + w_h * (Gamma_h + t(Gamma_h))
    }
  }
  
  var_rob_mat <- Q_inv %*% S %*% Q_inv
  se_rob <- sqrt(as.numeric(diag(var_rob_mat)))
  
  t_stat <- as.numeric(beta_hat) / se_rob
  p_val <- 2 * (1 - pt(abs(t_stat), df = n - p))
  
  sst <- sum((y_adj - mean(y_adj))^2)
  sse <- sum(e^2)
  r2 <- 1 - sse / sst
  
  diff_e <- e[2:n] - e[1:(n-1)]
  dw <- sum(diff_e^2) / sse
  
  tabla_coef <- data.frame(
    Estimacion = as.numeric(beta_hat),
    SE_OLS = se_ols,
    SE_Robusto = se_rob,
    t = t_stat,
    p_valor = p_val
  )
  rownames(tabla_coef) <- paste0("beta_", 0:(p-1))
  
  params <- list(tabla = tabla_coef, R2 = r2, sigma2 = sigma2, DW = dw)
  
  
  # El pronostico exp(a + bt) estima la mediana condicional y no la media.
  yhat <- if (tipo == "exponencial") {
    if (corregir_sesgo) exp(y_fit) * exp(sigma2 / 2) else exp(y_fit)
  } else {
    y_fit
  }
  
  pronosticar <- function(h) {
    if (h <= 0) stop("h debe ser positivo")
    t_futuro <- (n + 1):(n + h)
    X_fut <- if (tipo == "lineal") {
      cbind(1, t_futuro)
    } else if (tipo == "cuadratica") {
      cbind(1, t_futuro, t_futuro^2)
    } else {
      cbind(1, t_futuro)
    }
    pred_adj <- as.numeric(X_fut %*% beta_hat)
    
    if (tipo == "exponencial") {
      if (corregir_sesgo) exp(pred_adj) * exp(sigma2 / 2) else exp(pred_adj)
    } else {
      pred_adj
    }
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = params))
}

ajustar_holt <- function(y, alpha, beta) {
  if (any(is.na(y))) stop("El vector no puede tener NAs.")
  if (alpha <= 0 || alpha >= 1) stop("alpha debe estar en (0, 1).")
  if (beta <= 0 || beta >= 1) stop("beta debe estar en (0, 1).")
  n <- length(y)
  yhat <- rep(NA, n)
  L <- rep(NA, n)
  T_val <- rep(NA, n)
  
  L[1] <- y[1]
  T_val[1] <- 0
  
  if (n > 1) {
    yhat[2] <- L[1] + T_val[1]
    for (t in 2:n) {
      L[t] <- alpha * y[t] + (1 - alpha) * yhat[t]
      T_val[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * T_val[t - 1]
      if (t < n) {
        yhat[t + 1] <- L[t] + T_val[t]
      }
    }
  }
  
  L_n <- L[n]
  T_n <- T_val[n]
  
  pronosticar <- function(h) {
    if (h <= 0) stop("h debe ser positivo")
    return(L_n + T_n * (1:h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar, parametros = list(alpha = alpha, beta = beta, L_t = L, T_t = T_val)))
}

optimizar <- function(y, metodo, rejilla) {
  if (any(is.na(y))) stop("El vector no puede tener NAs.")
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