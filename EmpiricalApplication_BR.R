#################################################################################################################
###   Forecasting Tail Risk in the Brazilian Stock Market: An Out-of-Sample Comparison of VaR and ES Methods  ###
###   Authors: Esther Cleveston and Carlos Trucíos
#################################################################################################################
library(dplyr)
library(readxl)
library(stringr)
library(MASS)
library(rugarch)
library(eva)
library(betategarch)
library(MSGARCH)
library(Rcpp)
sourceCpp("CAViaR.cpp")
source("Utils.R")

# Import Data
data <- read_excel("precos_b3.xlsx", skip = 3, col_types = c("date", rep("numeric", 475)), na = c("", "-", NA)) |> 
  filter(Data > "2010-01-01" & Data < "2026-10-01") |> 
  filter(!if_all(where(is.numeric), is.na)) |> 
  select_if(~ !any(is.na(.))) |> 
  rename_with(~ str_remove(.x, "^(?s).*orig\n"), -Data) |> 
  mutate(across(where(is.numeric), ~ (log(.x) - lag(log(.x))) * 100)) |> 
  slice(-1)
names_stocks <- colnames(data)[-1]


data <- read_excel("precos_b3.xlsx", skip = 3, col_types = c("date", rep("numeric", 475)), na = c("", "-", NA))  |> 
  rename_with(~ str_remove(.x, "PU\najust p/ prov\nEm moeda orig\n")) |> 
  mutate(Data = as.Date(Data)) |>
  filter(Data > as.Date("2009-12-29") & Data < as.Date("2026-08-01")) |> 
  filter(!if_all(where(is.numeric), is.na)) |> 
  select(where(~ !any(is.na(.x)))) |> 
  mutate(across(where(is.numeric), ~ (log(.x) - lag(log(.x))) * 100)) |> 
  slice(-1)



## Baixar dados de preço diário e nós mesmos passar pra log-retorno

# Setup
alpha1 <- 0.010
alpha2 <- 0.025
alpha3 <- 0.050
in_sample <- 2500
B <- 3 * in_sample

garch_spec_n <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "norm")
garch_spec_t <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "std")
gjr_spec_n <- ugarchspec(variance.model = list(model = 'gjrGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "norm")
gjr_spec_t <- ugarchspec(variance.model = list(model = 'gjrGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "std")
egarch_spec_n <- ugarchspec(variance.model = list(model = 'eGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "norm")
egarch_spec_t <- ugarchspec(variance.model = list(model = 'eGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "std")
igarch_spec_n <- ugarchspec(variance.model = list(model = 'iGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "norm")
igarch_spec_t <- ugarchspec(variance.model = list(model = 'iGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "std")
ngarch_spec_n <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1, 1), submodel = "NAGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "norm")
ngarch_spec_t <- ugarchspec(variance.model = list(model = "fGARCH", garchOrder = c(1, 1), submodel = "NAGARCH"), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "std")
figarch_spec_n <- ugarchspec(variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'norm')
figarch_spec_t <- ugarchspec(variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'std')
ms_spec_n <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("norm", "norm")))
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))

names <- c("Normal", "T", "HS", "EVT", 
  "GARCH-N", "GARCH-T", "GJR-N", "GJR-T", "EGARCH-N", "EGARCH-T", "IGARCH-N", "IGARCH-T", "NGARCH-N", "NGARCH-T", "FIGARCH-N", "FIGARCH-T", "betategarch", "betategarchl", "betategarch2", "MS-N", "MS-T", 
  "FHS-GARCH-N", "FHS-GARCH-T", "FHS-GJR-N", "FHS-GJR-T", "FHS-EGARCH-N", "FHS-EGARCH-T", "FHS-IGARCH-N", "FHS-IGARCH-T", "FHS-NGARCH-N", "FHS-NGARCH-T", "FHS-FIGARCH-N", "FHS-FIGARCH-T", "FHS-betategarch", "FHS-betategarchl", "FHS-betategarch2", "FHS-MS-N", "FHS-MS-T",
  "EVT-GARCH-N", "EVT-GARCH-T", "EVT-GJR-N", "EVT-GJR-T", "EVT-EGARCH-N", "EVT-EGARCH-T", "EVT-IGARCH-N", "EVT-IGARCH-T", "EVT-NGARCH-N", "EVT-NGARCH-T", "EVT-FIGARCH-N", "EVT-FIGARCH-T", "EVT-betategarch", "EVT-betategarchl", "EVT-betategarch2", "EVT-MS-N", "EVT-MS-T",
  "CAVIAR-SYM", "CAVIAR-ASYM", "CAVIAR-ALD-SYM", "CAVIAR-ALD-ASYM", "CAVIAR-EVT-SYM", "CAVIAR-EVT-ASYM")
risk_levels <- c(1, 2, 5)
names <- paste(rep(names, each = length(risk_levels)), rep(risk_levels, times = length(names)))


# Rolling Window
for (stock in names_stocks){
  ret <- na.omit(data[[stock]])
  ntot <- length(ret)
  n_oos <- ntot - in_sample
  var_matrix <- matrix(0, nrow = n_oos, ncol = length(names) + 1, dimnames = list(NULL, c("returns", names)))
  es_matrix <- matrix(0, nrow = n_oos, ncol = length(names), dimnames = list(NULL,  names))
  for (i in 1:n_oos){
    print(i)
    set.seed(123 + i)
    r <- ret[i:(i + in_sample - 1)]
    mu <- mean(r)
    r_c <- scale(r, center = TRUE, scale = FALSE)
    # Pre-filtering AR(p) if needed
    if(abs(acf(r_c, plot = FALSE)$acf[2]) > 2/sqrt(nrow(r_c))) {
      ar_fit <- ar.yw(r_c, 3, aic = TRUE, se.fit = FALSE)
      mu <- mu + as.numeric(predict(ar_fit, newdata = r_c, n.ahead = 1, se.fit = FALSE))
      r_c <- na.omit(ar_fit$resid)
    }
    
    # Baseline Methods
    oos_normal <- sapply(c(alpha1, alpha2, alpha3), function(alpha) variance_covariance(r_c, alpha, distri = "Normal"))
    oos_student <- sapply(c(alpha1, alpha2, alpha3), function(alpha) variance_covariance(r_c, alpha, distri = "Student-t"))
    oos_hs <- sapply(c(alpha1, alpha2, alpha3), function(alpha) historical_simulation(r_c, alpha))
    oos_evtpot <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(r_c, alpha))
    
    # Volatility-Based
    ## Fitting models
    fit_garch_n <- ugarchfit(garch_spec_n, r_c, solver = "hybrid")
    fit_garch_t <- ugarchfit(garch_spec_t, r_c, solver = "hybrid")
    fit_gjr_n <- ugarchfit(gjr_spec_n, r_c, solver = "hybrid")
    fit_gjr_t <- ugarchfit(gjr_spec_t, r_c, solver = "hybrid")
    fit_egarch_n <- ugarchfit(egarch_spec_n, r_c, solver = "hybrid")
    fit_egarch_t <- ugarchfit(egarch_spec_t, r_c, solver = "hybrid")
    fit_igarch_n <- ugarchfit(igarch_spec_n, r_c, solver = "hybrid")
    fit_igarch_t <- ugarchfit(igarch_spec_t, r_c, solver = "hybrid")
    fit_ngarch_n <- ugarchfit(ngarch_spec_n, r_c, solver = "hybrid")
    fit_ngarch_t <- ugarchfit(ngarch_spec_t, r_c, solver = "hybrid")
    fit_figarch_n <- figarch_fit(figarch_spec_n, r_c)
    fit_figarch_t <- figarch_fit(figarch_spec_t, r_c)
    fit_betategarch <- tegarch(r_c, asym = FALSE, skew = FALSE, components = 1)
    fit_betategarch_lev <- tegarch(r_c, asym = TRUE, skew = FALSE, components = 1)
    fit_betategarch_2comp <- tegarch(r_c, asym = TRUE, skew = FALSE, components = 2)
    fit_ms_n <- msgarch_fit(ms_spec_n, r_c)
    fit_ms_t <- msgarch_fit(ms_spec_t, r_c)
    ## One-step-ahead volatility
    sigma_fore_garch_n <- ugarchforecast(fit_garch_n, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_garch_t <- ugarchforecast(fit_garch_t, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_gjr_n <- ugarchforecast(fit_gjr_n, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_gjr_t <- ugarchforecast(fit_gjr_t, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_egarch_n <- ugarchforecast(fit_egarch_n, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_egarch_t <- ugarchforecast(fit_egarch_t, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_igarch_n <- suppressWarnings(ugarchforecast(fit_igarch_n, n.ahead = 1)@forecast$sigmaFor[1, 1])
    sigma_fore_igarch_t <- suppressWarnings(ugarchforecast(fit_igarch_t, n.ahead = 1)@forecast$sigmaFor[1, 1])
    sigma_fore_ngarch_n <- ugarchforecast(fit_ngarch_n, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_ngarch_t <- ugarchforecast(fit_ngarch_t, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_figarch_n <- ugarchforecast(fit_figarch_n, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_figarch_t <- ugarchforecast(fit_figarch_t, n.ahead = 1)@forecast$sigmaFor[1, 1]
    sigma_fore_betategarch <- as.numeric(predict(fit_betategarch, n.ahead = 1, verbose = TRUE)$stdev)
    sigma_fore_betategarch_lev <- as.numeric(predict(fit_betategarch_lev, n.ahead = 1, verbose = TRUE)$stdev)
    sigma_fore_betategarch_2comp <- as.numeric(predict(fit_betategarch_2comp, n.ahead = 1, verbose = TRUE)$stdev)
    sigma_fore_ms_n <- as.numeric(predict(fit_ms_n, nahead = 1)$vol)
    sigma_fore_ms_t <- as.numeric(predict(fit_ms_t, nahead = 1)$vol)
    ## VaR and ES Parametric
    oos_garch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_garch_n, alpha, distri = "Norm"))
    oos_garch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_garch_t, alpha, distri = "Student-t", nu = coef(fit_garch_t)["shape"]))  
    oos_gjr_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_gjr_n, alpha, distri = "Norm"))
    oos_gjr_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_gjr_t, alpha, distri = "Student-t", nu = coef(fit_gjr_t)["shape"]))  
    oos_egarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_egarch_n, alpha, distri = "Norm"))
    oos_egarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_egarch_t, alpha, distri = "Student-t", nu = coef(fit_egarch_t)["shape"]))  
    oos_igarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_igarch_n, alpha, distri = "Norm"))
    oos_igarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_igarch_t, alpha, distri = "Student-t", nu = coef(fit_igarch_t)["shape"]))  
    oos_ngarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_ngarch_n, alpha, distri = "Norm"))
    oos_ngarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_ngarch_t, alpha, distri = "Student-t", nu = coef(fit_ngarch_t)["shape"]))  
    oos_figarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_figarch_n, alpha, distri = "Norm"))
    oos_figarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_figarch_t, alpha, distri = "Student-t", nu = coef(fit_figarch_t)["shape"]))  
    oos_betategarch <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_betategarch, alpha, distri = "Student-t", nu = coef(fit_betategarch)["df"]))  
    oos_betategarch_lev <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_betategarch_lev, alpha, distri = "Student-t", nu = coef(fit_betategarch_lev)["df"]))  
    oos_betategarch_2comp <- sapply(c(alpha1, alpha2, alpha3), function(alpha) parametric_var_es(sigma_fore_betategarch_2comp, alpha, distri = "Student-t", nu = coef(fit_betategarch_2comp)["df"]))  
    oos_ms_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) { aux <- Risk(fit_ms_n, nahead = 1, alpha = alpha); var <- as.numeric(aux$VaR); es <- as.numeric(aux$ES);  out <- c(var, es) })  
    oos_ms_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) { aux <- Risk(fit_ms_t, nahead = 1, alpha = alpha); var <- as.numeric(aux$VaR); es <- as.numeric(aux$ES);  out <- c(var, es) })  
    
    # FHS
    e_garch_n <- as.numeric(r_c / sigma(fit_garch_n))
    e_garch_t <- as.numeric(r_c / sigma(fit_garch_t))
    e_gjr_n <- as.numeric(r_c / sigma(fit_gjr_n))
    e_gjr_t <- as.numeric(r_c / sigma(fit_gjr_t))
    e_egarch_n <- as.numeric(r_c / sigma(fit_egarch_n))
    e_egarch_t <- as.numeric(r_c / sigma(fit_egarch_t))
    e_igarch_n <- as.numeric(r_c / sigma(fit_igarch_n))
    e_igarch_t <- as.numeric(r_c / sigma(fit_igarch_t))
    e_ngarch_n <- as.numeric(r_c / sigma(fit_ngarch_n))
    e_ngarch_t <- as.numeric(r_c / sigma(fit_ngarch_t))
    e_figarch_n <- as.numeric(r_c / sigma(fit_figarch_n))
    e_figarch_t <- as.numeric(r_c / sigma(fit_figarch_t))
    e_betategarch <- as.numeric(r_c / fitted(fit_betategarch, verbose = TRUE)$stdev)
    e_betategarch_lev <- as.numeric(r_c / fitted(fit_betategarch_lev, verbose = TRUE)$stdev)
    e_betategarch_2comp <- as.numeric(r_c / fitted(fit_betategarch_2comp, verbose = TRUE)$stdev)
    e_ms_n <- r_c / Volatility(fit_ms_n)
    e_ms_t <- r_c / Volatility(fit_ms_t)
    
    oos_fhs_garch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_garch_n, sigma_fore_garch_n, B, alpha))
    oos_fhs_garch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_garch_t, sigma_fore_garch_t, B, alpha))
    oos_fhs_gjr_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_gjr_n, sigma_fore_gjr_n, B, alpha))
    oos_fhs_gjr_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_gjr_t, sigma_fore_gjr_t, B, alpha))
    oos_fhs_egarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_egarch_n, sigma_fore_egarch_n, B, alpha))
    oos_fhs_egarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_egarch_t, sigma_fore_egarch_t, B, alpha))
    oos_fhs_igarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_igarch_n, sigma_fore_igarch_n, B, alpha))
    oos_fhs_igarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_igarch_t, sigma_fore_igarch_t, B, alpha))
    oos_fhs_ngarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_ngarch_n, sigma_fore_ngarch_n, B, alpha))
    oos_fhs_ngarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_ngarch_t, sigma_fore_ngarch_t, B, alpha))
    oos_fhs_figarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_figarch_n, sigma_fore_figarch_n, B, alpha))
    oos_fhs_figarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_figarch_t, sigma_fore_figarch_t, B, alpha))
    oos_fhs_betategarch <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_betategarch, sigma_fore_betategarch, B, alpha))
    oos_fhs_betategarch_lev <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_betategarch_lev, sigma_fore_betategarch_lev, B, alpha)) 
    oos_fhs_betategarch_2comp <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_betategarch_2comp, sigma_fore_betategarch_2comp, B, alpha))
    oos_fhs_ms_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_ms_n, sigma_fore_ms_n, B, alpha))  
    oos_fhs_ms_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) fhs(e = e_ms_t, sigma_fore_ms_t, B, alpha)) 
    
    # EVT_Based
    oos_evt_garch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_garch_n, alpha)) * sigma_fore_garch_n
    oos_evt_garch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_garch_t, alpha)) * sigma_fore_garch_t
    oos_evt_gjr_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_gjr_n, alpha)) * sigma_fore_gjr_n
    oos_evt_gjr_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_gjr_t, alpha)) * sigma_fore_gjr_t
    oos_evt_egarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_egarch_n, alpha)) * sigma_fore_egarch_n
    oos_evt_egarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_egarch_t, alpha)) * sigma_fore_egarch_t
    oos_evt_igarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_igarch_n, alpha)) * sigma_fore_igarch_n
    oos_evt_igarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_igarch_t, alpha)) * sigma_fore_igarch_t
    oos_evt_ngarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_ngarch_n, alpha)) * sigma_fore_ngarch_n
    oos_evt_ngarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_ngarch_t, alpha)) * sigma_fore_ngarch_t
    oos_evt_figarch_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_figarch_n, alpha)) * sigma_fore_figarch_n
    oos_evt_figarch_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_figarch_t, alpha)) * sigma_fore_figarch_t
    oos_evt_betategarch <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_betategarch, alpha)) * sigma_fore_betategarch
    oos_evt_betategarch_lev <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_betategarch_lev, alpha)) * sigma_fore_betategarch_lev
    oos_evt_betategarch_2comp <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_betategarch_2comp, alpha)) * sigma_fore_betategarch_2comp
    oos_evt_ms_n <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_ms_n, alpha)) * sigma_fore_ms_n
    oos_evt_ms_t <- sapply(c(alpha1, alpha2, alpha3), function(alpha) evt_pot(e_ms_t, alpha)) * sigma_fore_ms_t
    
    # CAViaR-Based
    oos_caviar_sym <- sapply(c(alpha1, alpha2, alpha3), function(alpha) CAViaR_reg(r_c, risklevel = alpha, type = "sym_abs"))
    oos_caviar_asym <- sapply(c(alpha1, alpha2, alpha3), function(alpha)CAViaR_reg(r_c, risklevel = alpha, type = "asym_slope"))
    oos_caviar_ald_sym <- sapply(c(alpha1, alpha2, alpha3), function(alpha) ALD(r_c, risklevel = alpha, type = "sym_abs"))
    oos_caviar_ald_asym <- sapply(c(alpha1, alpha2, alpha3), function(alpha) ALD(r_c, risklevel = alpha, type = "asym_slope"))
    oos_caviar_evt_sym <- sapply(c(alpha1, alpha2, alpha3), function(alpha) CAViaR_EVT(r_c, risklevel = alpha, type = "sym_abs"))
    oos_caviar_evt_asym <- sapply(c(alpha1, alpha2, alpha3), function(alpha) CAViaR_EVT(r_c, risklevel = alpha, type = "asym_slope"))
    
    # Save results
    var_matrix[i, 1] <- ret[i + in_sample]
    var_matrix[i, -1] <- mu + c(oos_normal[1, ], oos_student[1, ], oos_hs[1, ], oos_evtpot[1, ],
      oos_garch_n[1, ], oos_garch_t[1, ], oos_gjr_n[1, ], oos_gjr_t[1, ], oos_egarch_n[1, ], oos_egarch_t[1, ], oos_igarch_n[1, ], oos_igarch_t[1, ], oos_ngarch_n[1, ], oos_ngarch_t[1, ], oos_figarch_n[1, ], oos_figarch_t[1, ], oos_betategarch[1, ], oos_betategarch_lev[1, ], oos_betategarch_2comp[1, ], oos_ms_n[1, ], oos_ms_t[1, ],
      oos_fhs_garch_n[1, ], oos_fhs_garch_t[1, ], oos_fhs_gjr_n[1, ], oos_fhs_gjr_t[1, ], oos_fhs_egarch_n[1, ], oos_fhs_egarch_t[1, ], oos_fhs_igarch_n[1, ], oos_fhs_igarch_t[1, ],oos_fhs_ngarch_n[1, ], oos_fhs_ngarch_t[1, ], oos_fhs_figarch_n[1, ], oos_fhs_figarch_t[1, ], oos_fhs_betategarch[1, ], oos_fhs_betategarch_lev[1, ], oos_fhs_betategarch_2comp[1, ], oos_fhs_ms_n[1, ], oos_fhs_ms_t[1, ],
      oos_evt_garch_n[1, ], oos_evt_garch_t[1, ], oos_evt_gjr_n[1, ], oos_evt_gjr_t[1, ], oos_evt_egarch_n[1, ], oos_evt_egarch_t[1, ], oos_evt_igarch_n[1, ], oos_evt_igarch_t[1, ], oos_evt_ngarch_n[1, ], oos_evt_ngarch_t[1, ], oos_evt_figarch_n[1, ], oos_evt_figarch_t[1, ], oos_evt_betategarch[1, ], oos_evt_betategarch_lev[1, ],oos_evt_betategarch_2comp[1, ], oos_evt_ms_n[1, ], oos_evt_ms_t[1, ],
      oos_caviar_sym [1, ], oos_caviar_asym[1, ], oos_caviar_ald_sym[1, ], oos_caviar_ald_asym[1, ], oos_caviar_evt_sym [1, ], oos_caviar_evt_asym[1, ])
    es_matrix[i, ] <- mu + c(oos_normal[2, ], oos_student[2, ], oos_hs[2, ], oos_evtpot[2, ],
      oos_garch_n[2, ], oos_garch_t[2, ], oos_gjr_n[2, ], oos_gjr_t[2, ], oos_egarch_n[2, ], oos_egarch_t[2, ], oos_igarch_n[2, ], oos_igarch_t[2, ], oos_ngarch_n[2, ], oos_ngarch_t[2, ], oos_figarch_n[2, ], oos_figarch_t[2, ], oos_betategarch[2, ], oos_betategarch_lev[2, ], oos_betategarch_2comp[2, ], oos_ms_n[2, ], oos_ms_t[2, ],
      oos_fhs_garch_n[2, ], oos_fhs_garch_t[2, ], oos_fhs_gjr_n[2, ], oos_fhs_gjr_t[2, ], oos_fhs_egarch_n[2, ], oos_fhs_egarch_t[2, ], oos_fhs_igarch_n[2, ], oos_fhs_igarch_t[2, ],oos_fhs_ngarch_n[2, ], oos_fhs_ngarch_t[2, ], oos_fhs_figarch_n[2, ], oos_fhs_figarch_t[2, ], oos_fhs_betategarch[2, ], oos_fhs_betategarch_lev[2, ], oos_fhs_betategarch_2comp[2, ], oos_fhs_ms_n[2, ], oos_fhs_ms_t[2, ],
      oos_evt_garch_n[2, ], oos_evt_garch_t[2, ], oos_evt_gjr_n[2, ], oos_evt_gjr_t[2, ], oos_evt_egarch_n[2, ], oos_evt_egarch_t[2, ], oos_evt_igarch_n[2, ], oos_evt_igarch_t[2, ], oos_evt_ngarch_n[2, ], oos_evt_ngarch_t[2, ], oos_evt_figarch_n[2, ], oos_evt_figarch_t[2, ], oos_evt_betategarch[2, ], oos_evt_betategarch_lev[2, ],oos_evt_betategarch_2comp[2, ], oos_evt_ms_n[2, ], oos_evt_ms_t[2, ],
      oos_caviar_sym [2, ], oos_caviar_asym[2, ], oos_caviar_ald_sym[2, ], oos_caviar_ald_asym[2, ], oos_caviar_evt_sym [2, ], oos_caviar_evt_asym[2, ])
  }
  #Exportar resultados
  write.csv(var_matrix, paste0("VaR_Forecasts_", stock, ".csv"), row.names = FALSE)
  write.csv(es_matrix, paste0("ES_Forecasts_", stock, ".csv"), row.names = FALSE)
}



