https://github.com/StephanosDaza/st-2026-2-tarea1-herramientas-daza_palomino-esteban

# Tarea 1 - caja de herramientas de pronóstico

Series de Tiempo, Pregrado en Estadística, UNAL Medellín. Autor: Esteban Daza Palomino.

## Qué contiene el repositorio

| Archivo | Contenido | Dependencias |
|:---|:---|:---|
| `R/00-lectura.R` | `leer_serie()`: lee un objeto `ts` o un csv y devuelve un tibble con `t`, `fecha`, `y` y los atributos `frecuencia`, `fuente` y `unidad` | tibble |
| `R/01-graficos.R` | `graficar_serie()`, `acf_muestral()`, `banda_ruido_blanco()`, `correlograma()` y `graficar_optimizacion()` | ggplot2, patchwork |
| `R/02-metodos.R` | Los ocho métodos (`ajustar_media()`, `ajustar_mm()`, `ajustar_ses()`, `ajustar_dmm()`, `ajustar_tendencia()`, `ajustar_holt()`) y `optimizar()` | R base |
| `R/03-evaluacion.R` | `medidas()`, `ljung_box()`, `jarque_bera()`, `durbin_watson()` y `validar_errores()` | ggplot2, patchwork; usa `acf_muestral()` y `banda_ruido_blanco()` de `01-graficos.R` |
| `ejemplos/ejemplos.R` | Carga `R/`, ejecuta los ocho ejemplos y el contraejemplo, guarda las figuras, corre los bloques de verificación y escribe `sesion-info.txt` | dplyr, tidyr, purrr, tibble, ggplot2, patchwork |
| `informe/informe.qmd` | Análisis detallado de los ocho ejemplos y el contraejemplo; toma las cifras de la ejecución de `ejemplos.R` | Quarto, knitr |
| `informe/informe.html` | El informe renderizado | |
| `figs/` | Figuras generadas por `ejemplos.R`: `NN_ini`, `NN_opt`, `NN_err`, `NN_reg` y `NN_fin` para el ejemplo `NN` | |
| `sesion-info.txt` | Salida de `sessionInfo()` de la máquina en la que se cerró la entrega | |

## Cómo se corre

Desde la raíz del repositorio, en R 4.3 o posterior:

```r
install.packages(c("dplyr", "tidyr", "purrr", "tibble", "ggplot2", "patchwork"))
source("ejemplos/ejemplos.R")
```

Tiempo de ejecución observado: 90 segundos (R 4.5.1, Windows).

Para regenerar el informe, con Quarto instalado:

```
quarto render informe/informe.qmd
```

## Cómo se usan las funciones

Ejemplo mínimo con el suavizamiento exponencial simple sobre `Nile`:

```r
source("R/02-metodos.R")

y <- as.numeric(Nile)
ajuste <- ajustar_ses(y, alpha = 0.24)

head(ajuste$yhat)        # pronosticos de un paso; yhat[1] es NA (calentamiento)
ajuste$pronosticar(5)    # pronosticos extramuestrales Y_{T+1}, ..., Y_{T+5}
ajuste$parametros        # alpha usado y estado final

opt <- optimizar(y, "ajustar_ses", seq(0.02, 0.98, by = 0.02))
opt$optimo               # alpha con menor MSE de un paso
```

Todos los métodos devuelven la misma lista: `yhat`, `pronosticar` y `parametros`.

## Convenciones que fijan los números

- **Inicializaciones.** Media simple y suavizamiento exponencial simple: $\hat Y_2 = Y_1$. Media móvil: calentamiento de $k$ períodos. Doble media móvil: calentamiento de $2k - 1$ períodos. Tendencias: sin calentamiento, `yhat` son los valores ajustados. Holt: $L_1 = Y_1$ y $\hat T_1 = 0$.
- **ACF.** $r_h = c_h / c_0$ con divisor único $T$: $c_h = T^{-1}\sum_{t=1}^{T-h}(Y_t - \bar Y)(Y_{t+h} - \bar Y)$, hasta $m = \min\{\lfloor T/4 \rfloor, 24\}$. La banda es la de `plot.acf()`, $\pm z_{0.975}/\sqrt{n}$, con $n$ el número de observaciones de la sucesión graficada (de errores, cuando se grafican errores).
- **Partición.** El tramo de validación son las últimas $h = \min\{12, \lfloor 0.2T \rfloor\}$ observaciones, y al menos un ciclo en series estacionales.
- **Referente.** Ingenuo $\hat Y_{T+h} = Y_T$, o ingenuo estacional en series estacionales (AirPassengers y JohnsonJohnson). El MASE se escala con el MAD del ingenuo de un paso en el tramo de estimación.
- **Ljung-Box.** $p$ es el número de parámetros estimados: 0 en la media simple, 1 en la media móvil, la doble media móvil (la ventana $k$ se elige con `optimizar()`) y el suavizamiento exponencial simple, 2 en Holt y el número de coeficientes en las tendencias.
- **Tendencias.** Ecuaciones normales con `solve(crossprod(X), crossprod(X, y))`. Errores estándar robustos con núcleo de Bartlett y $\lfloor 4(T/100)^{2/9} \rfloor$ rezagos. La tendencia exponencial se estima sobre $\ln Y_t$; `corregir_sesgo = TRUE` multiplica por $e^{\hat\sigma^2/2}$.
- **Durbin-Watson.** En las tendencias se decide con las cotas $d_L$ y $d_U$ de Savin y White (1977), 5%, modelos con intercepto. En los demás métodos se usa la aproximación $z = (2 - d)\sqrt{N}/2 \approx N(0, 1)$.

## Resumen de resultados

MASE sobre el tramo de validación del método y de su referente. Ambos se escalan con el mismo MAD del ingenuo en el tramo de estimación, así que el método supera al referente cuando su MASE es menor.

| Ejemplo | Método | Serie | Parámetros | MASE método | MASE referente |
|---:|:---|:---|:---|---:|---:|
| 1 | Media simple | `discoveries` | ninguno | 0.919 | 1.136 |
| 2 | Media móvil | `discoveries` | $k = 9$ | 0.588 | 1.136 |
| 3 | Suavizamiento exponencial simple | `Nile` | $\alpha = 0.24$ | 0.806 | 0.835 |
| 4 | Doble media móvil | `WWWusage` | $k = 2$ | 1.368 | 7.625 |
| 5 | Tendencia lineal | `austres` | $\hat\beta_0 = 12963.29$, $\hat\beta_1 = 50.71$ | 3.974 | 6.091 |
| 6 | Tendencia cuadrática | `airmiles` | $\hat\beta_0 = 1420.12$, $\hat\beta_1 = -473.93$, $\hat\beta_2 = 74.33$ | 1.205 | 4.496 |
| 7 | Tendencia exponencial | `JohnsonJohnson` | $\hat a = -0.697$, $\hat\theta = 0.0428$, con corrección de sesgo | 3.586 | 6.529 |
| 8 | Holt lineal | `LakeHuron` | $\alpha = 0.95$, $\beta = 0.05$ | 2.090 | 2.164 |

El contraejemplo (media simple sobre `AirPassengers`) obtiene un MASE de 7.017 frente a 1.571 del ingenuo estacional.

## Declaración de uso de IA

Se usó un asistente de IA (Claude) como apoyo en tres frentes:

- **Lo que se pidió:** Revisar los scripts contra la rúbrica, proponer correcciones y
  estructurar el informe en Quarto a partir del análisis y los resultados propios.
- **Lo que se recibió:** Sugerencias de corrección para las funciones y para
  `ejemplos.R` (validaciones con `stopifnot()`, ACF propia en lugar de `acf()`,
  referente evaluado en validación, gráficas de optimización) y una estructura
  para el informe con las pruebas en seis elementos.
- **Lo que se verificó por cuenta propia:** La ejecución completa de `ejemplos.R`
  desde una sesión limpia, los bloques de verificación contra las funciones de R,
  las figuras generadas, las cotas de Durbin-Watson y las fuentes de las series
  con `help()`.
