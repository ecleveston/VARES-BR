################
##   Utils   ###
################


##############################
###   Auxiliary Functions  ###
##############################

dt_c <- function(x, df, scale) {
  dt(x / scale, df = df) / scale
}

ALD_loss = function(params, r, risklevel = 0.05, type = "sym_abs") {
  n = length(r)
  Qu = rep(0, n)
  Qu[1] = quantile(r, risklevel) 
  if (type == "asym_slope") {
    for (i in 2:n) {
      Qu[i] = params[1] + params[2]*Qu[i - 1] + params[3]*max(r[i - 1], 0) + params[4]*min(r[i - 1], 0)
    }
    ES = (1 + exp(params[5]))*Qu
  } else {
    for (i in 2:n) {
      Qu[i] = params[1] + params[2]*Qu[i - 1] + params[3]*abs(r[i - 1])
    }
    ES = (1 + exp(params[4]))*Qu
  }
  ll_ald = -mean(log((risklevel - 1)/ES[2:n]) + (r[2:n] - Qu[2:n])*(risklevel - (r[2:n] <= Qu[2:n]))/(risklevel*ES[2:n]))
  return(ll_ald)
}

figarch_fit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  n <- length(data)
  opt_methods <-  c("solnp", "nlminb",  "gosolnp", "nloptr", "lbfgs", "hybrid")
  while (is_error && k < length(opt_methods)) {
    k <- k + 1
    is_error <- FALSE
    expr <- tryCatch({
      fit_model <- ugarchfit(spec, data, solver = opt_methods[k], fit.control = list(trunclag = 1000))
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
  
  if (is_error) stop("All optimization methods failed to converge.")
  return(fit_model)
}

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

###############################
###         Methods         ###
###############################

historical_simulation <- function(r, alpha) {
  var <- quantile(r, alpha)
  es <- mean(r[r < var])
  out <- c(var, es)
}

variance_covariance <- function(r, alpha, distri = c("Normal", "Student-t")) {
  distri <- match.arg(distri)
  if (distri == "Normal") {
    n <- length(r)
    sigma <- sd(r) * sqrt((n - 1) / n)
    q <- qnorm(alpha)
    var <- sigma * q
    es <- -sigma * dnorm(q) / alpha
  } else {
    fit_student <- MASS::fitdistr(r, densfun = dt_c, start = list(df = 5, scale = sd(r)), lower = c(df = 2.01, scale = 1e-8)) 
    nu <- fit_student$estimate["df"]
    s  <- fit_student$estimate["scale"] # equiv sigma/{sqrt(nu/(nu -2))}
    q <- qt(alpha, df = nu)
    var <- s * q
    es <- - s * dt(q, df = nu) * (nu + q^2) / ((nu - 1) * alpha)
  }
  out <- c(var, es)
}

evt_pot <- function(r, alpha, theta = 0.85){
  x <- -r
  u <- as.numeric(quantile(x, theta))
  fit_pareto <- gpdFit(x, threshold = u)
  params <- fit_pareto$par.ests
  xi    <- unname(params["Shape (Intercept)"])
  beta <- unname(params["Scale (Intercept)"])
  N  <- length(r)
  Nu <- fit_pareto$n.exceed
  if (abs(xi) < 1e-8) {
    var <- u + beta * log(Nu / (N * alpha))
    es  <- var + beta
  } else {
    var <- u + (beta / xi) * ((Nu / (N * alpha))^xi - 1)
    if (xi >= 1) stop("xi >= 1")
    es <- (var + beta - xi * u) / (1 - xi)
  }
  out <- c(-var, -es)
}



fhs <- function(e, sigma_fore, B, alpha) {
  eb <- sample(e, B, replace = TRUE)
  var_e <- quantile(eb, prob = alpha)
  es_e  <- mean(eb[eb < var_e])
  var <- var_e * sigma_fore
  es <- es_e * sigma_fore
  out <- c(var, es)
}

parametric_var_es <- function(sigma, alpha, distri = c("Normal", "Student-t"), nu = NULL){
  distri <- match.arg(distri)
  if (distri == "Normal") {
    q <- qnorm(alpha)
    var <- sigma * q
    es <- -sigma * dnorm(q) / alpha
  } else {
    k <- sqrt(nu / (nu - 2))
    q <- qt(alpha, df = nu)
    var <- sigma / k * q
    es <- - sigma / k * dt(q, df = nu) * (nu + q^2) / ((nu - 1) * alpha) 
  }
  out <- c(var, es)
}

CAViaR <- function(retornos, risklevel = 0.05, type = "sym_abs", par_ini = NULL) {
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
    params1 = suppressWarnings(optim(par = pini[smallest[1],], fn = CAViaR_loss, r = retornos, risklevel = risklevel, type = type))
    params2 = suppressWarnings(optim(par = pini[smallest[2],], fn = CAViaR_loss, r = retornos, risklevel = risklevel, type = type))
    params3 = suppressWarnings(optim(par = pini[smallest[3],], fn = CAViaR_loss, r = retornos, risklevel = risklevel, type = type))
    best_of_three = order(c(params1$value, params2$value, params3$value),decreasing = F)[1]
    if (best_of_three == 1) params = params1
    if (best_of_three == 2) params = params2
    if (best_of_three == 3) params = params3
    if (!is.null(par_ini)) {
      params4 = suppressWarnings(optim(par = par_ini, fn = CAViaR_loss, r = retornos, risklevel = risklevel, type = type))
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

CAViaR_reg <- function(r, risklevel = 0.05, type = "sym_abs", par_ini = NULL) {
  n <- length(r)
  caviar_method <- CAViaR(r, risklevel, type, par_ini)
  var <- caviar_method[[1]]
  r_exceed <- r[r < var[1:n]] 
  var_exceed <- var[1:n][r < var[1:n]]
  delta <- as.numeric(coef(lm(r_exceed ~ var_exceed - 1)))
  es <- delta * var
  out <- c(tail(var, 1), tail(es, 1))
}

CAViaR_EVT <- function(r, risklevel = 0.05, type = "sym_abs", theta = 0.075, par_ini = NULL) {
  n <- length(r)
  caviar_method <- CAViaR(r, theta, type, par_ini)
  Qu <- caviar_method[[1]]
  u <- 0
  std_residual <- (r[1:n] - Qu[1:n]) / Qu[1:n]             # r[1:n] - Qu[1:n] are the residuals
  std_res_exceed <- std_residual[std_residual > u]         # Select only exceedances
  params_gpd <- gpdFit(std_res_exceed, threshold = 0, method = "mle")$par.sum$Estimate
  params_gpd_shape <- params_gpd[2]                        #varepsilon
  params_gpd_scale <- params_gpd[1]                        #beta
  std_res_Qu <- u + (params_gpd_scale/params_gpd_shape)*((theta/risklevel)^(params_gpd_shape) - 1)   # Eq (7.18) in McNeil, Frey and Embrechts (2005) or Quantile (1-p) in Eq (9.7) of Engle & Manganelli (2004)
  std_res_ES <- (std_res_Qu + params_gpd_scale - params_gpd_shape*u)/(1 - params_gpd_shape)          # Eq (9.19) in Engle & Manganelli (2004) [with sigma_t = q_t] or Eq (7.19) in McNeil, Frey and Embrechts (2005)
  var <- Qu + std_res_Qu * Qu   # As in Engle & Manganelli (2004) pag 132
  es <- u + std_res_ES * Qu    # As in Taylor (2019)
  out <- c(tail(var, 1), tail(es, 1))
}

ALD <- function(r, risklevel = 0.05, type = "sym_abs", caviar_params = NULL) {
    n = length(r)
    Qu = rep(0, n + 1)
    type_code = ifelse(type == "asym_slope", 1, 0)
    if (is.null(caviar_params)) caviar_params = CAViaR(r, risklevel, type = type, par_ini = caviar_params)[[2]]
    gamma_pini = ALD_grid(caviar_params, r, risklevel, type_code)
    params_ini = c(caviar_params, gamma_pini)
    ald_params = tryCatch({
      suppressWarnings(optim(par = params_ini, fn = ALD_loss, r = r, risklevel = risklevel, type = type))
    }, error = function(e) {
      NULL
    })
    # Grid
    if (type == "asym_slope") {
      pini = matrix(0, ncol = 5, nrow = 1000)
      pini[,1] = -0.2*runif(1000)
      pini[,2] = runif(1000)
      pini[,3] = -runif(1000)
      pini[,4] = 0.5*runif(1000)
      pini[,5] = -0.8 + 3*runif(1000)
      loss = apply(pini, 1, ALD_loss, r, risklevel, type)
    } else{
      pini = matrix(0, ncol = 4, nrow = 1000)
      pini[,1] = -0.2*runif(1000)
      pini[,2] = runif(1000)
      pini[,3] = -runif(1000)
      pini[,4] = -0.8 + 3*runif(1000)
      loss = apply(pini, 1, ALD_loss, r, risklevel, type)
    }
    # Optimisation
    smallest = order(loss ,decreasing = F)[1:3]
    params1 = suppressWarnings(optim(par = pini[smallest[1],], fn = ALD_loss, r = r, risklevel = risklevel, type = type))
    params2 = suppressWarnings(optim(par = pini[smallest[2],], fn = ALD_loss, r = r, risklevel = risklevel, type = type))
    params3 = suppressWarnings(optim(par = pini[smallest[3],], fn = ALD_loss, r = r, risklevel = risklevel, type = type))
    best_of_three = order(c(params1$value, params2$value, params3$value),decreasing = F)[1]
    if (best_of_three == 1) params = params1
    if (best_of_three == 2) params = params2
    if (best_of_three == 3) params = params3
    if (!is.null(ald_params)) {
      if (ald_params$value < params$value) params = ald_params
    }
    params = params$par
    Qu[1] = quantile(r, risklevel) 
    if (type == "asym_slope") {
      for (i in 2:(n + 1)) {
        Qu[i] = params[1] + params[2]*Qu[i - 1] + params[3]*max(r[i - 1], 0) + params[4]*min(r[i - 1], 0)
      }
      ES = (1 + exp(params[5]))*Qu
    } else {
      for (i in 2:(n + 1)) {
        Qu[i] = params[1] + params[2]*Qu[i - 1] + params[3]*abs(r[i - 1])
      }
      ES = (1 + exp(params[4]))*Qu
    }
    out <- c(tail(Qu, 1), tail(ES, 1))
}


























