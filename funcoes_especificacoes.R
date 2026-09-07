## -----------------------------------------------------------------------------------------------------------------------------------
especificacao_gjr_n <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                   variance.model = list(model = 'gjrGARCH', garchOrder = c(1, 1)),
                   distribution = 'norm')


## -----------------------------------------------------------------------------------------------------------------------------------
especificacao_gjr_t <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), 
                   variance.model = list(model = 'gjrGARCH', garchOrder = c(1, 1)),
                   distribution = 'std')


## -----------------------------------------------------------------------------------------------------------------------------------
calc_VaR_hist <- function(x, alpha) {
  quantile(x, probs = alpha, na.rm = TRUE)
}


## -----------------------------------------------------------------------------------------------------------------------------------
garch_fit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  n <- length(data)
  opt_methods <- c("hybrid", "nloptr", "solnp", "gosolnp", "lbfgs", "nlminb")
  fit_model <- NULL
  
  while (is_error && k < length(opt_methods)) {
    k <- k + 1
    is_error <- FALSE
    expr <- tryCatch({
      fit_model <- ugarchfit(spec, data, solver = opt_methods[k])
      warns <- character(0)
      fc <- withCallingHandlers(
        ugarchforecast(fit_model, n.ahead = 1),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")}
      )
      if (fit_model@fit$convergence != 0) stop("no convergence")
      if (any(!is.finite(fc@forecast$sigmaFor))) stop("NaN or Inf in forecast")
      if (any(grepl("Positivity Contraints", warns))) stop("Positivity Constraints NOT")
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (!isTRUE(expr)) is_error <- TRUE
  }
  
  if (!is_error) {
    return(fit_model)
  } else {
    return(NULL) 
  }
}

## -----------------------------------------------------------------------------------------------------------------------------------
especificacao1  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), distribution = 'norm')
especificacao2  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), distribution = 'std')
especificacao3  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'eGARCH', garchOrder = c(1, 1)), distribution = 'norm')
especificacao4  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'eGARCH', garchOrder = c(1, 1)), distribution = 'std')
especificacao5  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'iGARCH', garchOrder = c(1, 1)), distribution = 'norm')
especificacao6  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'iGARCH', garchOrder = c(1, 1)), distribution = 'std')
especificacao7  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), distribution = 'norm')
especificacao8  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), distribution = 'std')
especificacao9  <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), distribution = 'norm')
especificacao10 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), distribution = 'std')
especificacao11 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'fGARCH', garchOrder = c(1, 1), submodel = "NAGARCH"), distribution = 'norm')
especificacao12 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), variance.model = list(model = 'fGARCH', garchOrder = c(1, 1), submodel = "NAGARCH"), distribution = 'std')


## -----------------------------------------------------------------------------------------------------------------------------------
CAViaR_loss = function(params, retornos, risklevel, type) {
  n = length(retornos)
  Qu = rep(0, n)
  Qu[1] = quantile(retornos, risklevel)
  if (type == 1) {
    for (i in 2:n) {
      Qu[i] = params[1] + params[2]*Qu[i - 1] + params[3]*max(retornos[i - 1], 0) + params[4]*min(retornos[i - 1], 0)
    }
  } else {
      for (i in 2:n) {
        Qu[i] = params[1] + params[2]*Qu[i - 1] + params[3]*abs(retornos[i - 1])
      }
  }
  loss_function = mean((risklevel - (retornos[2:n] < Qu[2:n])) * (retornos[2:n] - Qu[2:n]))
  return(loss_function)
}

CAViaR = function(retornos, risklevel = 0.05, type = "sym_abs", par_ini = NULL) {
  if (type %in% c("asym_slope", "sym_abs")) {
    n = length(retornos)
    Qu = rep(0, n + 1)
    type = ifelse(type == "asym_slope", 1, 0)
    if (type == 1) {
      pini = matrix(0, ncol = 4, nrow = 10000)
      pini[,1] = runif(10000, -0.4, 0)
      pini[,2] = runif(10000)
      pini[,3] = runif(10000, -1.5, -0.01)
      pini[,4] = 0.5*runif(10000)
      loss = apply(pini, 1, CAViaR_loss, retornos, risklevel, type)
    } else{
      pini = matrix(0, ncol = 3, nrow = 10000)
      pini[,1] = runif(10000, -0.4, 0)
      pini[,2] = runif(10000)
      pini[,3] = runif(10000, -1.5, -0.01)
      loss = apply(pini, 1, CAViaR_loss, retornos, risklevel, type)
    }
    smallest = order(loss ,decreasing = F)[1:3]
    params1 = suppressWarnings(optim(par = pini[smallest[1],], fn = CAViaR_loss, retornos = retornos, risklevel = risklevel, type = type))
    params2 = suppressWarnings(optim(par = pini[smallest[2],], fn = CAViaR_loss, retornos = retornos, risklevel = risklevel, type = type))
    params3 = suppressWarnings(optim(par = pini[smallest[3],], fn = CAViaR_loss, retornos = retornos, risklevel = risklevel, type = type))
    best_of_three = order(c(params1$value, params2$value, params3$value),decreasing = F)[1]
    if (best_of_three == 1) params = params1
    if (best_of_three == 2) params = params2
    if (best_of_three == 3) params = params3
    if (!is.null(par_ini)) {
      params4 = suppressWarnings(optim(par = par_ini, fn = CAViaR_loss, retornos = retornos, risklevel = risklevel, type = type))
      if (params4$value < params$value) params = params4
    }
    params = params$par
    Qu[1] = quantile(retornos, risklevel) 
    if (type == 1) {
      for (i in 2:(n + 1)) {
        Qu[i] = params[1] + params[2]*Qu[i - 1] + params[3]*max(retornos[i - 1], 0) + params[4]*min(retornos[i - 1], 0)
      }
    } else {
      for (i in 2:(n + 1)) {
        Qu[i] = params[1] + params[2]*Qu[i - 1] + params[3]*abs(retornos[i - 1])
      }
    }
    return(list(Qu, params))
  } else{
    return(print('type option not found, should be either asym_slope or sym_abs.'))
  }
}


## -----------------------------------------------------------------------------------------------------------------------------------
evt <- function(alpha, retornos){
  perdas <- -retornos[retornos < 0]
  u <- quantile(perdas, 0.85)
  fit_pareto <- gpdFit(perdas, threshold = u)
  parametros_pareto <- fit_pareto$par.ests
  xi <- parametros_pareto["Shape (Intercept)"]
  sigma <- parametros_pareto["Scale (Intercept)"]
  N <- length(perdas)
  Nu <- fit_pareto$n.exceed
  var_evt <- u + (sigma / xi) * ((((N / Nu) * (alpha)))^(-xi) - 1)
  return(-as.numeric(var_evt))
}


## -----------------------------------------------------------------------------------------------------------------------------------
msgarch_fit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  opt_methods <- c("BFGS", "Nelder-Mead", "CG", "SANN", "solnp")
  expr <- NULL
  while (is_error == TRUE && k < 6) {
    k <- k + 1
    expr <- tryCatch({
      if (k < 5) {
        fit_model <- FitML(spec, data, ctr = list(do.se = FALSE, do.plm = FALSE, OptimFUN = function(vPw, f_nll, spec, data, do.plm){
          out <- stats::optim(vPw, f_nll, spec = spec, data = data, do.plm = do.plm, method = opt_methods[k])}))
      } else {
        fit_model <- FitML(spec, data, ctr = list(do.se = FALSE, do.plm = TRUE, OptimFUN = function(vPw, f_nll, spec, data, do.plm) {
          fn_obj <- function(pars) {
            val <- f_nll(pars, spec = spec, data = data, do.plm = do.plm)
            return(sum(val))
          }
          out <- Rsolnp::solnp(pars = vPw, fun = fn_obj)
          out$value <- out$values[length(out$values)]
          out$convergence <- ifelse(out$convergence == 0, 0, 1)
          return(out)
        }))
      }
      TRUE
    }, error = function(e) {
      FALSE
    }, warning = function(cond) {
      FALSE
    })
    if (isTRUE(expr)) is_error <- FALSE
  }
  return(fit_model)
}


## -----------------------------------------------------------------------------------------------------------------------------------
gas_fit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  opt_methods <- c("BFGS", "Nelder-Mead", "CG", "SANN")
  fit_model <- NULL
  while (is_error == TRUE && k < length(opt_methods)) {
    k <- k + 1
    tryCatch(
      expr <- {
        fit_model <- UniGASFit(spec, data, Compute.SE = FALSE, fn.optimizer = function(par0, data, GASSpec, FUN) {
          optimizer = optim(par0, FUN, data = data, GASSpec = GASSpec, method = opt_methods[k],
            control = list(trace = 0), hessian = FALSE)
          out = list(pars = optimizer$par,
            value = optimizer$value,
            hessian = optimizer$hessian,
            convergence = optimizer$convergence)
          return(out)
        })
        if (fit_model@Estimates$optimiser$convergence != 0) stop("No convergence")
        TRUE
      },
      error = function(e) {
        FALSE
      })
    if (isTRUE(expr)) is_error <- FALSE
  }
  return(fit_model)
}


## -----------------------------------------------------------------------------------------------------------------------------------
gas_spec_n <- UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
gas_spec_t <- UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
ms_spec_n <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("norm", "norm")))
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))


## -----------------------------------------------------------------------------------------------------------------------------------
#Simulação Histórica
calc_ES_hist <- function(retornos, alpha){
  perdas_alem_do_var <- retornos[retornos <= calc_VaR_hist(retornos, alpha)]
  return(mean(perdas_alem_do_var))
}

es_gauss <- function(retornos, alpha) {
  media <- mean(retornos, na.rm = TRUE)
  desvio_padrao <- sd(retornos, na.rm = TRUE)  
  z_alpha <- qnorm(alpha, mean = 0, sd = 1)
  pdf_alpha <- dnorm(z_alpha, mean = 0, sd = 1)
  return(media - desvio_padrao * (pdf_alpha / alpha))
}

es_student <- function(retornos, alpha){
  quantil <- qt(alpha, graus_de_liberdade)
  densidade <- dt(quantil, graus_de_liberdade)
  termo_es <- -(densidade / alpha) * ((graus_de_liberdade + quantil^2) / (graus_de_liberdade - 1))
  return(mu + desvio_padrao * termo_es)
}
#EVT
evt_ES <- function(alpha, retornos){
  perdas <- -retornos[retornos < 0]
  u <- quantile(perdas, 0.85)
  fit_pareto <- gpdFit(perdas, threshold = u)
  parametros_pareto <- fit_pareto$par.ests
  xi <- parametros_pareto["Shape (Intercept)"]
  sigma <- parametros_pareto["Scale (Intercept)"]
  N <- length(perdas)
  Nu <- fit_pareto$n.exceed
  var_evt <- u + (sigma / xi) * ((((N / Nu) * (alpha)))^(-xi) - 1)
  if(xi < 1) {
    evt_es <- (var_evt / (1 - xi)) + ((sigma - xi * u) / (1 - xi))
  } else {
    evt_es <- NA 
  }
  return(-as.numeric(evt_es))
}

