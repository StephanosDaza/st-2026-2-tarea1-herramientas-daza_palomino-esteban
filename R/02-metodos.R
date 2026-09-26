ajustar_media <- function(y) {
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y),
    "y no puede tener NA" = !anyNA(y),
    "y debe tener al menos 2 observaciones" = length(y) >= 2
  )
  n <- length(y)
  yhat <- rep(NA_real_, n)
  
  media <- y[1]
  yhat[2] <- media
  for (t in 2:n) {
    media <- media + (y[t] - media) / t
    if (t < n) yhat[t + 1] <- media
  }
  estado_final <- media
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h >= 1 && h == round(h))
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar,
              parametros = list(estado_final = estado_final)))
}

ajustar_mm <- function(y, k) {
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y),
    "y no puede tener NA" = !anyNA(y),
    "k debe ser un entero mayor o igual a 2" = length(k) == 1 && k == round(k) && k >= 2,
    "la ventana k no puede ser mayor que la longitud de la serie" = k <= length(y)
  )
  n <- length(y)
  yhat <- rep(NA_real_, n)
  suma_ventana <- 0
  
  for (t in 1:n) {
    suma_ventana <- suma_ventana + y[t]
    if (t > k) suma_ventana <- suma_ventana - y[t - k]
    if (t >= k && t < n) yhat[t + 1] <- suma_ventana / k
  }
  estado_final <- suma_ventana / k
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h >= 1 && h == round(h))
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar,
              parametros = list(k = k, estado_final = estado_final)))
}

ajustar_ses <- function(y, alpha) {
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y),
    "y no puede tener NA" = !anyNA(y),
    "y debe tener al menos 2 observaciones" = length(y) >= 2,
    "alpha debe estar en (0, 1)" = length(alpha) == 1 && alpha > 0 && alpha < 1
  )
  n <- length(y)
  yhat <- rep(NA_real_, n)
  
  yhat[2] <- y[1]
  if (n > 2) {
    for (t in 2:(n - 1)) {
      e_t <- y[t] - yhat[t]
      yhat[t + 1] <- yhat[t] + alpha * e_t
    }
  }
  
  e_n <- y[n] - yhat[n]
  estado_final <- yhat[n] + alpha * e_n
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h >= 1 && h == round(h))
    return(rep(estado_final, h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar,
              parametros = list(alpha = alpha, estado_final = estado_final)))
}

ajustar_dmm <- function(y, k) {
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y),
    "y no puede tener NA" = !anyNA(y),
    "k debe ser un entero mayor o igual a 2" = length(k) == 1 && k == round(k) && k >= 2,
    "la ventana k no puede ser mayor que la longitud de la serie" = k <= length(y),
    "la serie debe tener al menos 2k - 1 observaciones" = 2 * k - 1 <= length(y)
  )
  n <- length(y)
  yhat <- rep(NA_real_, n)
  E_t <- rep(NA_real_, n)
  beta1_t <- rep(NA_real_, n)
  mm1 <- rep(NA_real_, n)
  mm2 <- rep(NA_real_, n)
  
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
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h >= 1 && h == round(h))
    return(E_n + B_n * (1:h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar,
              parametros = list(k = k, E_t = E_t, beta1_t = beta1_t, MM = mm1, DMM = mm2)))
}

ajustar_tendencia <- function(y, tipo, corregir_sesgo = FALSE) {
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y),
    "y no puede tener NA" = !anyNA(y),
    "tipo debe ser 'lineal', 'cuadratica' o 'exponencial'" =
      length(tipo) == 1 && tipo %in% c("lineal", "cuadratica", "exponencial"),
    "corregir_sesgo debe ser TRUE o FALSE" = is.logical(corregir_sesgo) && length(corregir_sesgo) == 1
  )
  if (tipo == "exponencial") {
    stopifnot("La tendencia exponencial exige una serie con valores estrictamente positivos" = all(y > 0))
  }
  n <- length(y)
  p <- if (tipo == "cuadratica") 3 else 2
  stopifnot("La serie es demasiado corta para el numero de coeficientes" = n > p)
  
  y_adj <- if (tipo == "exponencial") log(y) else y
  
  tt <- 1:n
  X <- if (tipo == "cuadratica") cbind(1, tt, tt^2) else cbind(1, tt)
  
  beta_hat <- solve(crossprod(X), crossprod(X, y_adj))
  
  y_fit <- as.numeric(X %*% beta_hat)
  e <- y_adj - y_fit
  
  sigma2 <- sum(e^2) / (n - p)
  
  Q_inv <- solve(crossprod(X))
  se_ols <- sqrt(as.numeric(diag(sigma2 * Q_inv)))
  
  L_hac <- floor(4 * (n / 100)^(2/9))
  eX <- X * e
  S <- crossprod(eX)
  if (L_hac > 0 && n > L_hac) {
    for (h in 1:L_hac) {
      Gamma_h <- crossprod(eX[-(1:h), , drop = FALSE], eX[-((n - h + 1):n), , drop = FALSE])
      w_h <- 1 - h / (L_hac + 1)
      S <- S + w_h * (Gamma_h + t(Gamma_h))
    }
  }
  
  var_rob_mat <- Q_inv %*% S %*% Q_inv
  se_rob <- sqrt(as.numeric(diag(var_rob_mat)))
  
  t_ols <- as.numeric(beta_hat) / se_ols
  p_ols <- 2 * (1 - pt(abs(t_ols), df = n - p))
  t_rob <- as.numeric(beta_hat) / se_rob
  p_rob <- 2 * (1 - pt(abs(t_rob), df = n - p))
  
  sst <- sum((y_adj - mean(y_adj))^2)
  sse <- sum(e^2)
  r2 <- 1 - sse / sst
  
  diff_e <- e[2:n] - e[1:(n - 1)]
  dw <- sum(diff_e^2) / sse
  
  tabla_coef <- data.frame(
    Estimacion = as.numeric(beta_hat),
    SE_OLS = se_ols,
    SE_Robusto = se_rob,
    t_OLS = t_ols,
    p_OLS = p_ols,
    t_Robusto = t_rob,
    p_Robusto = p_rob
  )
  rownames(tabla_coef) <- if (tipo == "exponencial") c("a", "theta") else paste0("beta_", 0:(p - 1))
  
  params <- list(tipo = tipo, tabla = tabla_coef, R2 = r2, sigma2 = sigma2, DW = dw,
                 rezagos_HAC = L_hac, gl = n - p, residuos = e, ajustados = y_fit,
                 corregir_sesgo = corregir_sesgo)
  if (tipo == "exponencial") {
    params$beta0 <- exp(beta_hat[1])
    params$beta1 <- exp(beta_hat[2])
  }
  
  # El pronostico exp(a + bt) estima la mediana condicional y no la media.
  yhat <- if (tipo == "exponencial") {
    if (corregir_sesgo) exp(y_fit) * exp(sigma2 / 2) else exp(y_fit)
  } else {
    y_fit
  }
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h >= 1 && h == round(h))
    t_futuro <- (n + 1):(n + h)
    X_fut <- if (tipo == "cuadratica") cbind(1, t_futuro, t_futuro^2) else cbind(1, t_futuro)
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
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y),
    "y no puede tener NA" = !anyNA(y),
    "y debe tener al menos 2 observaciones" = length(y) >= 2,
    "alpha debe estar en (0, 1)" = length(alpha) == 1 && alpha > 0 && alpha < 1,
    "beta debe estar en (0, 1)" = length(beta) == 1 && beta > 0 && beta < 1
  )
  n <- length(y)
  yhat <- rep(NA_real_, n)
  L <- rep(NA_real_, n)
  T_val <- rep(NA_real_, n)
  
  L[1] <- y[1]
  T_val[1] <- 0
  
  yhat[2] <- L[1] + T_val[1]
  for (t in 2:n) {
    L[t] <- alpha * y[t] + (1 - alpha) * yhat[t]
    T_val[t] <- beta * (L[t] - L[t - 1]) + (1 - beta) * T_val[t - 1]
    if (t < n) yhat[t + 1] <- L[t] + T_val[t]
  }
  
  L_n <- L[n]
  T_n <- T_val[n]
  
  pronosticar <- function(h) {
    stopifnot("h debe ser un entero positivo" = length(h) == 1 && h >= 1 && h == round(h))
    return(L_n + T_n * (1:h))
  }
  
  return(list(yhat = yhat, pronosticar = pronosticar,
              parametros = list(alpha = alpha, beta = beta, L_t = L, T_t = T_val)))
}

optimizar <- function(y, metodo, rejilla) {
  metodos_validos <- c("ajustar_mm", "ajustar_ses", "ajustar_dmm", "ajustar_holt")
  stopifnot(
    "y debe ser un vector numerico" = is.numeric(y),
    "y no puede tener NA" = !anyNA(y),
    "metodo no soportado" = length(metodo) == 1 && metodo %in% metodos_validos,
    "rejilla debe ser un vector numerico o un data.frame" = is.numeric(rejilla) || is.data.frame(rejilla)
  )
  
  if (is.numeric(rejilla) && is.null(dim(rejilla))) {
    nombre <- if (metodo %in% c("ajustar_mm", "ajustar_dmm")) "k" else "alpha"
    rejilla <- setNames(data.frame(rejilla), nombre)
  }
  if (metodo == "ajustar_holt") {
    stopifnot("La rejilla de Holt debe tener dos columnas (alpha, beta)" = ncol(rejilla) == 2)
    names(rejilla) <- c("alpha", "beta")
  }
  
  mse_vals <- numeric(nrow(rejilla))
  
  for (i in 1:nrow(rejilla)) {
    params <- as.list(rejilla[i, , drop = FALSE])
    
    ajuste <- switch(metodo,
                     ajustar_mm   = ajustar_mm(y, k = params[[1]]),
                     ajustar_ses  = ajustar_ses(y, alpha = params[[1]]),
                     ajustar_dmm  = ajustar_dmm(y, k = params[[1]]),
                     ajustar_holt = ajustar_holt(y, alpha = params[[1]], beta = params[[2]])
    )
    
    e <- y - ajuste$yhat
    mse_vals[i] <- mean(e^2, na.rm = TRUE)
  }
  
  rejilla$mse <- mse_vals
  optimo <- rejilla[which.min(rejilla$mse), , drop = FALSE]
  
  return(list(rejilla = rejilla, optimo = optimo))
}
