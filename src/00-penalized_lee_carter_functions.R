# Penalized Lee Carter Model

# Configuration ---------------------------------------------------

PLCfitConfig <- function (
  # penalties
  lambda_ax = 0,
  lambda_bx = 0,
  lambda_kt = 0,
  lambda_ridge = 1e-6,
  lambda_ax_target = 0,
  lambda_bx_target = 0,
  lambda_kt_target = 0,
  ax_target = 0,
  bx_target = 0,
  kt_target = 0,
  # optimization
  init_pars = NULL,
  init_mode = 'basic',
  outer_maxit = 50,
  bfgs_maxit = 10,
  dev_stop_crit = 1e-3,
  # kt outlier detection
  kt_outlier_removal = FALSE,
  kt_outlier_zscore = 4,
  kt_outlier_runmed_k = 5,
  kt_outlier_absolute = FALSE,
  # kt adjustment
  kt_adjustment = 'poisson',
  # log
  plot_progress = TRUE
) {
  config <- as.list(environment())
  return(config)
}


# Initialization --------------------------------------------------

#' Fill NAs in Matrix
PLCfillMatrixNAs <- function (X) {
  # row-wise forward pass over columns with special rule for leading/trailing NA
  X_ <- zoo::na.approx(t(X), rule = 2) |> t()
  colnames(X_) <- colnames(X)
  return(X_)
}

#' Get Initial LC Estimates
PLCgetInitEstimates <- function (Mcx, mode) {

  if (identical(mode, 'basic')) {
    ax = log(rowMeans(Mcx, na.rm = TRUE))
    bx = rnorm(nrow(Mcx), 0, 0.1)
    kt = rep(0, ncol(Mcx))
    init <- list(
      ax = ax,
      bx = bx,
      kt = kt
    )
  }
  if (identical(mode, 'svd')) {
    Mcx <- PLCfillMatrixNAs(Mcx)
    log_Mcx <- log(Mcx+1e-10)
    ax <- rowMeans(log_Mcx)
    log_Mcx_normalized <- sweep(log_Mcx, 1, ax, FUN = '-')
    SVD <- svd(log_Mcx_normalized)
    init <- list(
      ax = ax,
      bx = SVD$u[,1],
      kt = SVD$v[,1]
    )
  }

  return(init)

}

# Fit -------------------------------------------------------------

#' Fit Penalized Lee Carter Model
PLCfit <- function(
    Dxt, Ext, label = 'unknown population', W = NULL,
    config = PLCfitConfig()
) {

  N = nrow(Dxt)
  m = ncol(Dxt)

  Mcx <- Dxt/Ext
  Mcx[Mcx==0] <- min(Mcx[Mcx!=0], na.rm = TRUE)
  Mcx[is.nan(Mcx)] <- NA

  if (is.null(W)) { W <- PLCcreateWeightMatrix(Dxt = Dxt, Ext = Ext) }
  # information about which kt's are not estimated from data due to
  # the weight for the corresponding year being 0
  kt_zerowght <- apply(W, 2, function (x) all(x == 0))

  # initialize parameters
  if (!is.null(config$init_pars)) {
    theta_ <- config$init_pars
    theta <- PLCconstraints(theta_$ax, theta_$bx, theta_$kt)
  } else {
    theta_ <- PLCgetInitEstimates(Mcx, mode = config$init_mode)
    theta <- PLCconstraints(theta_$ax, theta_$bx, theta_$kt)
  }

  # smoothing parameter
  lambda_ax <- config$lambda_ax
  lambda_bx <- config$lambda_bx
  lambda_kt <- config$lambda_kt
  lambda_ridge <- config$lambda_ridge
  lambda_ax_target <- config$lambda_ax_target
  lambda_bx_target <- config$lambda_bx_target
  lambda_kt_target <- config$lambda_kt_target
  # difference matrix
  # penalize second differences
  D_ax <- diff(diag(N), diff = 2)
  DD_ax <- t(D_ax)%*%D_ax
  D_bx <- diff(diag(N), diff = 2)
  DD_bx <- t(D_bx)%*%D_bx
  D_kt <- diff(diag(m), diff = 2)
  DD_kt <- t(D_kt)%*%D_kt

  eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
  epsilon <- log(Mcx)-eta

  maxit = config$outer_maxit
  bfgs_maxit = config$bfgs_maxit

  dev <- rep(NA, maxit+1)
  dev[1] <- PLCpoissonDeviance(eta, Dxt, Ext, W)
  cat(dev[1], ' (', sep = '')
  r2 <- rep(NA, maxit+1)
  r2[1] <- PLCr2(epsilon, eta, W)
  cat(r2[1], ') ', sep = '')

  if (isTRUE(config$plot_progress)) {
    PLCplotFitDiagnostics(dev, epsilon, eta, theta, maxit, outlier = kt_zerowght)
  }

  # the central fitting loop
  # - ax, bx, and kt are optimized sequentially for up to
  #   <maxit> iterations
  # - deviance reduction is criterion for convergence
  # - if kt outlier detection is on, the loop works until convergence
  #   and then, based on the converged kt, detects kt outliers, removes
  #   them (0 weight), and optimizes again until convergence or maxit
  kt_outlier <- rep(FALSE, m)
  detected_kt_outliers <- NULL
  outlier_removed <- 0
  for (i in 1:maxit) {

    theta_$ax <- optim(
      theta_$ax,
      PLCpoissonLL,
      bx = theta_$bx, kt = theta_$kt,
      Dcx = Dxt, Ecx = Ext, n = N, m = m,
      DD_ax = DD_ax, DD_bx = DD_bx, DD_kt = DD_kt,
      lambda_ax = lambda_ax, lambda_bx = lambda_bx,
      lambda_kt = lambda_kt,
      lambda_ridge = lambda_ridge,
      lambda_ax_target = lambda_ax_target,
      lambda_bx_target = lambda_bx_target,
      lambda_kt_target = lambda_kt_target,
      ax_target = config$ax_target,
      bx_target = config$bx_target,
      kt_target = config$kt_target,
      W = W,
      method = 'BFGS', hessian = FALSE,
      control = list(fnscale = -1, trace = 0, maxit = bfgs_maxit)
    )[['par']]
    theta_$bx <- optim(
      theta_$bx,
      PLCpoissonLL,
      ax = theta_$ax, kt = theta_$kt,
      Dcx = Dxt, Ecx = Ext, n = N, m = m,
      DD_ax = DD_ax, DD_bx = DD_bx, DD_kt = DD_kt,
      lambda_ax = lambda_ax, lambda_bx = lambda_bx,
      lambda_kt = lambda_kt,
      lambda_ridge = lambda_ridge,
      lambda_ax_target = lambda_ax_target,
      lambda_bx_target = lambda_bx_target,
      lambda_kt_target = lambda_kt_target,
      ax_target = config$ax_target,
      bx_target = config$bx_target,
      kt_target = config$kt_target,
      W = W,
      method = 'BFGS', hessian = FALSE,
      control = list(fnscale = -1, trace = 0, maxit = bfgs_maxit)
    )[['par']]
    theta_$kt <- optim(
      theta_$kt,
      PLCpoissonLL,
      ax = theta_$ax, bx = theta_$bx,
      Dcx = Dxt, Ecx = Ext, n = N, m = m,
      DD_ax = DD_ax, DD_bx = DD_bx, DD_kt = DD_kt,
      lambda_ax = lambda_ax, lambda_bx = lambda_bx,
      lambda_kt = lambda_kt,
      lambda_ridge = lambda_ridge,
      lambda_ax_target = lambda_ax_target,
      lambda_bx_target = lambda_bx_target,
      lambda_kt_target = lambda_kt_target,
      ax_target = config$ax_target,
      bx_target = config$bx_target,
      kt_target = config$kt_target,
      W = W,
      method = 'BFGS', hessian = FALSE,
      control = list(fnscale = -1, trace = 0, maxit = bfgs_maxit)
    )[['par']]

    theta <- PLCconstraints(theta_$ax, theta_$bx, theta_$kt)
    eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
    epsilon <- log(Mcx)-eta

    dev[i+1] <- PLCpoissonDeviance(eta, Dxt, Ext, W)
    cat(dev[i+1], ' (', sep = '')
    r2[i+1] <- PLCr2(epsilon, eta, W)
    cat(r2[i+1], ') ', sep = '')

    if (isTRUE(config$plot_progress)) {
      PLCplotFitDiagnostics(dev, epsilon, eta, theta, maxit,
                            outlier = kt_outlier|kt_zerowght)
    }

    deviance_reduction <- abs(log(dev[i+1])-log(dev[i]))
    deviance_convergence <- deviance_reduction < config$dev_stop_crit

    # stop if deviance converged and outliers should not be removed
    if (!config$kt_outlier_removal && deviance_convergence) {
      break
    }
    # stop if deviance converged and outliers have been removed
    if (config$kt_outlier_removal && outlier_removed && deviance_convergence) {
      break
    }
    # remove outliers if deviance converged and continue with optimization
    if (config$kt_outlier_removal && !outlier_removed && deviance_convergence) {
      outlier_removed <- 1
      detected_kt_outliers <- PLCdetectKtOutliers(
        theta$kt,
        k = config$kt_outlier_runmed_k,
        z = config$kt_outlier_zscore,
        absolute = config$kt_outlier_absolute
      )
      kt_outlier <- detected_kt_outliers$outlier_lgl
      if (!any(kt_outlier)) { break }
      W[,kt_outlier] <- 0
      # initialize parameters
      if (!is.null(config$init_pars)) {
        theta_ <- config$init_pars
      } else {
        theta_ <- PLCgetInitEstimates(Mcx)
      }
    }

  }

  # kt adjustment
  if (isTRUE(config$kt_adjustment == 'poisson')) {
    # for the kt adjustment all data is used unless
    # its Na or Ext is 0
    W_complete <- PLCcreateWeightMatrix(Dxt = Dxt, Ext = Ext)
    theta_$kt <- optim(
      theta_$kt,
      PLCpoissonLL,
      ax = theta_$ax, bx = theta_$bx,
      Dcx = Dxt, Ecx = Ext, n = N, m = m,
      DD_ax = DD_ax, DD_bx = DD_bx, DD_kt = DD_kt,
      lambda_ax = lambda_ax, lambda_bx = lambda_bx,
      lambda_kt = 0,
      lambda_ridge = 0,
      lambda_ax_target = lambda_ax_target,
      lambda_bx_target = lambda_bx_target,
      lambda_kt_target = 0,
      ax_target = config$ax_target,
      bx_target = config$bx_target,
      kt_target = config$kt_target,
      W = W_complete,
      method = 'BFGS', hessian = FALSE,
      control = list(fnscale = -1, trace = 0, maxit = 100)
    )[['par']]

    theta <- PLCconstraints(theta_$ax, theta_$bx, theta_$kt)
    eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
    epsilon <- log(Mcx)-eta

    if (isTRUE(config$plot_progress)) {
      PLCplotFitDiagnostics(dev, epsilon, eta, theta, maxit,
                            outlier = kt_outlier|kt_zerowght)
    }
  }

  return(
    list(
      data = list(
        exposure = Ext,
        deaths = Dxt
      ),
      model_parameters = theta,
      unconstrained_model_parameters = theta_,
      predicted_log_rates = eta,
      model_diagnostics = list(
        deviance = dev,
        r2 = r2,
        epsilon = epsilon
      ),
      meta = list(
        population = label,
        date_of_model_fit = Sys.time(),
        config = config,
        years = colnames(Dxt),
        ages = rownames(Dxt),
        W = W,
        # which kt's where not fitted to data
        kt_exclude = kt_outlier|kt_zerowght,
        kt_outlier = detected_kt_outliers
      )
    )
  )

}

PLCcreateWeightMatrix <- function (Dxt, Ext, nacols = NULL, narows = NULL) {
  n <- nrow(Dxt)
  m <- ncol(Dxt)
  W <- matrix(1, nrow = n, ncol = m)
  W[is.na(Dxt)|is.na(Ext)|(Ext == 0)] <- 0
  W[,nacols] <- 0
  W[narows,] <- 0
  return(W)
}

# Prediction ------------------------------------------------------

#' Predict Lee-Carter Surface Given Parameters
PLCpredict <- function (ax, bx, kt) {
  ax + outer(bx, kt)
}

#' Apply LC Constraints to Parameters
PLCconstraints <- function (ax, bx, kt) {
  c1 <- mean(kt, na.rm = TRUE)
  c2 <- sum(bx, na.rm = TRUE)
  theta <- list(
    ax = ax + c1 * bx,
    bx = bx / c2,
    kt = c2 * (kt - c1)
  )
  return(theta)
}

# Penalized likelihood --------------------------------------------

#' Lee-Carter Penalized Poisson Log-Likelihood
PLCpoissonLL <- function (
    ax, bx, kt, Dcx, Ecx, n, m,
    DD_ax, DD_bx, DD_kt,
    lambda_ax, lambda_bx, lambda_kt, lambda_ridge,
    lambda_ax_target, lambda_bx_target, lambda_kt_target,
    ax_target = 0, bx_target = 0, kt_target = 0,
    W
) {
  theta <- PLCconstraints(ax, bx, kt)
  eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
  mu <- exp(eta)
  ll <-
    sum(W*(Dcx*eta-Ecx*mu), na.rm = TRUE) -
    lambda_ax*t(theta$ax)%*%DD_ax%*%theta$ax -
    lambda_bx*t(theta$bx)%*%DD_bx%*%theta$bx -
    lambda_kt*t(theta$kt)%*%DD_kt%*%theta$kt -
    lambda_ridge*sum((theta$bx+theta$ax)^2) -
    lambda_ax_target*sum((ax-ax_target)^2) -
    lambda_bx_target*sum((bx-bx_target)^2) -
    lambda_kt_target*sum((kt-kt_target)^2)
  return(ll)
}

# Diagnostics -----------------------------------------------------

#' Lee-Carter Poisson Deviance
PLCpoissonDeviance <- function (eta, Dcx, Ecx, W) {
  sum(
    W*2*(Dcx*(log(Dcx)-(eta+log(Ecx)))-Dcx+(exp(eta)*Ecx)),
    na.rm = TRUE
  )
}

#' Lee-Carter R2
PLCr2 <- function(epsilon, eta, W) {
  rss <- sum(W*epsilon^2, na.rm = TRUE)
  observed <- W*(eta+epsilon)
  tss <- sum((observed-mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- round((1 - rss/tss)*100, 4)
}

#' Plot PLC Fit Diagnostics
PLCplotFitDiagnostics <- function(dev, epsilon, eta, theta, maxit, outlier,
                                  theme = PLCtosTheme) {

  plot_layout_matrix <- matrix(NA, 3, 6)
  plot_layout_matrix[1:2,1:3] <- 1
  plot_layout_matrix[1,4:6] <- 2
  plot_layout_matrix[2,4:6] <- 3
  plot_layout_matrix[3,1:2] <- 4
  plot_layout_matrix[3,3:4] <- 5
  plot_layout_matrix[3,5:6] <- 6
  layout(plot_layout_matrix)
  par(
    bg = theme$bg,
    col = theme$col,
    col.axis = theme$col.axis,
    col.lab = theme$col.lab,
    col.main = theme$col.main,
    fg = theme$fg,
    mar = theme$mar
  )

  N = length(theta$ax)
  m = length(theta$kt)

  plot.new()
  ldev <- log(dev)
  plot.window(xlim = c(0, maxit), y = c(ldev[1]-4, ldev[1]), log = 'y')
  axis(1); axis(2)
  title(xlab = 'Iteration', ylab = 'Deviance', main = 'Log-deviance reduction profile')
  polygon(
    x = c(0, maxit, maxit, 0),
    y = c(ldev[1]-4, ldev[1]-4, dev[1], dev[1]),
    col = theme$deviance_bg_col, border = FALSE
  )
  grid(col = 'black', lty = 1, nx = NULL, ny = NULL, lwd = 0.2, equilogs = FALSE)
  points(x = 0:maxit, y = ldev, col = theme$deviance_line_col,
         pch = theme$deviance_pch, cex = theme$deviance_cex)
  lines(x = 0:maxit, y = ldev, col = theme$deviance_line_col)
  PlotMatrix(epsilon, type = 'd',
             clip = c(-0.4, 0.4), # clip at +- 50%
             main = 'Residuals', xlab = 't', ylab = 'x',
             divcol = theme$color_scale_cols_residuals,
             N = theme$color_scale_n_residuals)
  PlotMatrix(eta, type = 'c', clip = c(-12, 0),
             main = 'Estimated log m(x,t)', xlab = 't', ylab = 'x',
             concol = theme$color_scale_cols_log_mx,
             N = theme$color_scale_n_log_mx)
  plot(x = 1:N, y = theta$ax, xlab = 'x', ylab = 'a(x)',
       main = 'a(x) estimates')
  plot(x = 1:N, y = theta$bx, xlab = 'x', ylab = 'b(x)',
       main = 'b(x) estimates')
  plot(x = 1:m, y = theta$kt, xlab = 't', ylab = 'k(t)',
       main = 'k(t) estimates', pch = ifelse(outlier, 4, 1))
}

# Outlier detection -----------------------------------------------

#' Detect Outliers in Lee-Carter kt Vector
PLCdetectKtOutliers <- function (kt, k = 5, z = 3, absolute = TRUE) {
  # dkt <- diff(kt)
  # dkt_avg <- median(dkt)
  # dkt_sdv <- 1.48*median(abs(dkt-dkt_avg))
  # dkt_z <- abs(dkt-dkt_avg)/dkt_sdv
  # outlier <- ifelse(dkt_z > z, TRUE, FALSE)
  kt_smooth <- runmed(kt, k = k, endrule = 'constant')
  res <- kt-kt_smooth
  #res_avg <- mean(res)
  res_sdv <- 1.48*median(abs(res))
  res_z <- res/res_sdv
  if (isTRUE(absolute)) { res_z <- abs(res_z) }
  outlier <- ifelse(res_z > z, TRUE, FALSE)

  out <- list(outlier_lgl = outlier,
              outlier_zsc = res_z[outlier],
              residual_sd = res_sdv)

  return(out)
}

# Forecasting -----------------------------------------------------

#' Forecast
#'
#' @param theta List of unconstrained model parameters
PLCforecast <- function (theta, h, nsim,
                         sd_estimation = 'classic',
                         drift_estimation = 'classic',
                         jumpoff_estimation = 'classic',
                         kt_exclude = NULL,
                         p_crisis = NULL,
                         m_crisis = NULL) {

  kt <- theta$kt
  # exclude kt elements from variance and drift estimation
  kt[kt_exclude] <- NA

  # estimate linear drift of kt term
  dkt <- diff(kt)
  dkt_nona <- na.exclude(dkt)
  kt_drift <- switch(
    drift_estimation,
    classic = mean(dkt_nona),
    robust = median(dkt_nona)
  )
  # estimate standard deviation of kt innovations
  dkt_sd <- switch(
    sd_estimation,
    classic = sd(dkt_nona),
    weighted = sqrt(tail(EMVar(dkt_nona, length(dkt_nona)-1), 1)),
    robust = 1.48*median(abs(dkt_nona-median(dkt_nona)))
  )

  # jump-off kt
  kt_init <- switch(
    jumpoff_estimation,
    classic = tail(na.exclude(kt), 1),
    robust = median(tail(na.exclude(kt), 5))
  )

  # simulate kc trajectories as random walk with drift
  kt_sim <- matrix(NA, nrow = nsim, ncol = h)
  for (k in 1:nsim) {
    kt_sim[k,] <-
      kt_init + cumsum(kt_drift + rnorm(h, sd = dkt_sd))
  }

  # simulate mortality crises
  # sample from a 2 state Markov-Chain
  # giving the probability of kt switching from and to the
  # excluded state. if the excluded states mark mortality crises,
  # the the MC represents the annual probability to enter such a crisis,
  # and – if the crisis has been entered – the annual probability of
  # leaving the crisis. once a crisis has been entered, sample the
  # magnitude of the crisis
  if (!is.null(p_crisis) & !is.null(m_crisis)) {
    vt_sim <- matrix(NA, nrow = nsim, ncol = h)
    for (k in 1:nsim) {
      in_crisis <- 0
      for (i in 1:h) {
        vt_sim[k,i] <- rbinom(1, 1, p_crisis[in_crisis+1])
        in_crisis <- vt_sim[k,i]
      }
      vt_sim[k,] <-
        vt_sim[k,]*sample(m_crisis, size = h, replace = TRUE)
    }
    kt_sim <- kt_sim + vt_sim
  }

  Eta_sim <- array(NA, dim = c(length(theta$ax), h, nsim))
  for (k in 1:nsim) {
    Eta_sim[,,k] <- PLCpredict(theta$ax, theta$bx, kt_sim[k,])
  }

  list(
    Eta_forecast_sim = Eta_sim,
    kt_drift = kt_drift,
    dkt_sd = dkt_sd,
    kt_sim = kt_sim
  )

}

#' Window Exponential Smoothing Variance Estimate
EMVar <- function(x, n){
  alpha <- 2/(n+1)
  # exponential moving average
  ema <- rep(NA, n-1)
  ema[n]<- mean(x[1:n])

  for (i in (n+1):length(x)){
    ema[i]<-alpha*x[i] + (1-alpha)*ema[i-1]
  }
  # exponential moving variance
  delta <- x - lag(ema)
  emvar <- rep(NA, n-1)
  emvar[n] <- ifelse(n==1,0,var(x[1:n]))
  for(i in (n+1):length(x)){
    emvar[i] <-  (1-alpha)*(emvar[i-1] + alpha*delta[i]^2)
  }
  return(emvar)
}

# Lifetables ------------------------------------------------------

LifeExpectancyFromMortality <- function (mx) {
  n <- length(mx)
  I <- diag(n)
  U <- head(rbind(0, diag(exp(-mx))),-1)
  N <- solve(I-U)
  ex <- colSums(N)
  return(ex)
}

LifespanVarianceFromMortality <- function (mx) {
  n <- length(mx)
  I <- diag(n)
  U <- head(rbind(0, diag(exp(-mx))),-1)
  N <- solve(I-U)
  ex <- colSums(N)
  n2 <- t(ex)%*%(2*N-I)
  vx <- n2 - (ex^2)
  return(vx)
}

# Misc ------------------------------------------------------------

#' Create Cohort Lee Carter Binary Weight Matrix
#'
#' Given a matrix of cohort x age training data with
#' dimensions m x n, return a binary matrix returning 0
#' for data unknown at the year of the jump-off cohort
#'
#' @param m number of cohorts in training
#' @param n number of age groups in training
#'
#' @return Binary weight matrix
#'
#' @examples
#' CreateWeightMatrix(10, 10)
# PLCcreateWeightMatrix <- function (m, n) {
#   W <- matrix(1, nrow = n, ncol = m)
#   # for (s in 0:(m-1)) {
#   #   cohort <- m-s
#   #   admissible_ages <- 1:min(s+1, n)
#   #   W[admissible_ages,cohort] <- 1
#   # }
#   return(W)
# }

PlotMatrix <- function (
    X, N = 100, type = 'c', clip = c(NA,NA),
    main = '', xlab = '', ylab = '',
    concol = 'cubehelix', divcol = 'RdBu'
) {
  require(rje)
  XX <- t(X[nrow(X):1,])
  if (!is.na(clip[1])) { XX[XX<=clip[1]] <- clip[1] }
  if (!is.na(clip[2])) { XX[XX>=clip[2]] <- clip[2] }
  XX <- RescaleToUnit(XX, clip[1], clip[2])
  col <- switch (concol[1],
                 'cubehelix' = cubeHelix(N),
                 colorRampPalette(concol)(N))
  if (type == 'd') {
    col <- switch (divcol[1],
                   'RdBu' = hcl.colors(N, 'RdBu'),
                   colorRampPalette(divcol)(N))
  }
  if (type == 'q') {
    N <- quantile(XX, probs = seq(0, 1, length.out = N),
                  na.rm = TRUE)
  }
  XX_col <- matrix(
    as.integer(cut(XX, breaks = N, labels = FALSE)),
    nrow = nrow(XX), ncol = ncol(XX)
  )
  image(XX_col, xaxt = 'n', yaxt = 'n', col = col, main = main,
        xlab = xlab, ylab = ylab, useRaster = TRUE)
}

RescaleToUnit <- function (X, a = NULL, b = NULL) {
  if (is.na(a)) {
    a <- min(X, na.rm = TRUE)
  }
  if (is.na(b)) {
    b <- max(X, na.rm = TRUE)
  }
  X_ <- (X-a)/(b-a)
  return(X_)
}

#' PLC Plot Fit TOS Theme
PLCtosTheme <- list(
  bg = 'black',
  col = '#EDAC31',
  col.axis = '#EDAC31',
  col.lab = '#E1511F',
  col.main = '#E1511F',
  fg = '#EDAC31',
  mar = c(2,2,2,2),
  deviance_bg_col = '#EDAC31',
  deviance_grid_col = 'black',
  deviance_grid_lwd = 0.2,
  deviance_line_col = '#F33826',
  deviance_pch = 16,
  deviance_cex = 1.5,
  color_scale_cols_residuals = c('#175223', 'black', '#E1511F'),
  color_scale_n_residuals = 5,
  color_scale_cols_log_mx = c('#0285D0', '#28A578', '#EDAC31', '#E1511F', '#F33826'),
  color_scale_n_log_mx = 10
)

#' PLC Plot Fit Default Theme
PLCdefaultTheme <- list(
  bg = 'black',
  col = 'white',
  col.axis = 'white',
  col.lab = 'white',
  col.main = 'white',
  fg = 'white',
  mar = c(2,2,2,2),
  deviance_bg_col = 'grey30',
  deviance_grid_col = 'black',
  deviance_grid_lwd = 0.2,
  deviance_line_col = 'white',
  deviance_pch = 16,
  deviance_cex = 1.5,
  color_scale_cols_residuals = c('#5EDF82', '#497252', 'black', '#6F6388', '#c682fa'),
  color_scale_n_residuals = 50,
  color_scale_cols_log_mx = 'cubehelix',
  color_scale_n_log_mx = 10
)
