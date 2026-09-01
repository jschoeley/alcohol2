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
  init_mode = 'svd',
  outer_maxit = 500,
  bfgs_maxit = 20,
  dev_stop_crit = 1e-5,
  # kt outlier detection
  kt_outlier_removal = FALSE,
  kt_outlier_zscore = 4,
  kt_outlier_runmed_k = 5,
  kt_outlier_absolute = FALSE,
  # kt adjustment
  kt_adjustment = 'poisson',
  # bx sensitivity repair
  bx_sensitivity_threshold = Inf,
  # plotting
  plot_progress = FALSE,
  plot_final = TRUE
) {
  config <- as.list(environment())
  return(config)
}

# Initialization --------------------------------------------------

#' Fill NAs in Matrix
PLCfillMatrixNAs <- function (X) {
  # row-wise forward pass over columns with special rule for
  # leading/trailing NA
  X_ <- zoo::na.approx(t(X), rule = 2) |> t()
  colnames(X_) <- colnames(X)
  return(X_)
}

#' Get Initial LC Estimates
PLCgetInitEstimates <- function (Mcx, mode = c('randombx', 'svd')) {

  mode <- match.arg(mode)

  if (identical(mode, 'randombx')) {
    mean_rates <- rowMeans(Mcx, na.rm = TRUE)
    fallback_rate <- mean(Mcx[is.finite(Mcx) & Mcx > 0], na.rm = TRUE)
    if (!is.finite(fallback_rate)) { fallback_rate <- 1e-10 }
    mean_rates[!is.finite(mean_rates) | mean_rates <= 0] <- fallback_rate
    ax = log(mean_rates)
    bx = rnorm(nrow(Mcx), 0, 0.1)
    kt = rep(0, ncol(Mcx))
    init <- list(
      ax = ax,
      bx = bx,
      kt = kt
    )
  }
  if (identical(mode, 'svd')) {
    if (anyNA(Mcx)) {
      warning(
        paste(
          'SVD initialization received Mcx with missing values;',
          'missing values will be interpolated.'
        ),
        call. = FALSE
      )
    }
    positive_rates_per_row <- rowSums(is.finite(Mcx) & Mcx > 0)
    insufficient_rows <- which(positive_rates_per_row < 2)
    if (length(insufficient_rows) > 0) {
      warning(
        sprintf(
          paste(
            'SVD initialization requires at least two positive rates per',
            "age row; falling back to 'randombx' initialization. Insufficient rows: %s."
          ),
          paste(insufficient_rows, collapse = ', ')
        ),
        call. = FALSE
      )
      return(PLCgetInitEstimates(Mcx, mode = 'randombx'))
    }
    Mcx <- PLCfillMatrixNAs(Mcx)
    log_Mcx <- log(Mcx+1e-10)
    ax <- rowMeans(log_Mcx)
    log_Mcx_normalized <- sweep(log_Mcx, 1, ax, FUN = '-')
    SVD <- svd(log_Mcx_normalized)
    init <- list(
      ax = ax,
      bx = SVD$u[,1],
      kt = SVD$d[1]*SVD$v[,1]
    )
  }

  return(init)

}

# Optimization ----------------------------------------------------

#' Optimize One Lee-Carter Parameter Vector
PLCoptim <- function(
    theta, parameter = c('ax', 'bx', 'kt'), objective_args, W, maxit
) {
  parameter <- match.arg(parameter)

  reserved_args <- c('ax', 'bx', 'kt', 'W')
  if (any(reserved_args %in% names(objective_args))) {
    stop(
      'objective_args must not contain ax, bx, kt, or W.',
      call. = FALSE
    )
  }

  objective <- function(par) {
    candidate <- theta
    candidate[[parameter]] <- par
    do.call(
      PLCpoissonLL,
      c(
        candidate[c('ax', 'bx', 'kt')],
        objective_args,
        list(W = W)
      )
    )
  }

  result <- optim(
    par = theta[[parameter]],
    fn = objective,
    method = 'BFGS', hessian = FALSE,
    control = list(fnscale = -1, trace = 0, maxit = maxit)
  )

  theta[[parameter]] <- result$par
  theta <- do.call(
    PLCconstraints,
    theta[c('ax', 'bx', 'kt')]
  )

  list(theta = theta, optim = result)
}

#' Repair Overly Sensitive bx Estimates
PLCrepairBx <- function(theta, threshold, kt_exclude, age_labels) {
  if (length(threshold) != 1 || is.na(threshold) || threshold < 0) {
    stop('bx_sensitivity_threshold must be a single non-negative number.',
         call. = FALSE)
  }

  calculate_sensitivities <- function(theta) {
    kt <- theta$kt
    kt[kt_exclude] <- NA
    kt_sd <- sd(na.exclude(diff(kt)))
    setNames(abs(theta$bx) * kt_sd, age_labels)
  }

  original_sensitivities <- calculate_sensitivities(theta)
  rejected_bx <-
    is.finite(original_sensitivities) & original_sensitivities > threshold
  if (any(rejected_bx)) {
    repaired_bx <- theta$bx
    repaired_bx[rejected_bx] <- 0
    theta <- PLCconstraints(theta$ax, repaired_bx, theta$kt)
  }

  list(
    theta = theta,
    bx_scaled_sensitivities = calculate_sensitivities(theta),
    rejected_bx = rejected_bx,
    diagnostics = list(
      applied = any(rejected_bx),
      rejected_ages = age_labels[rejected_bx],
      original_bx_scaled_sensitivities = original_sensitivities
    )
  )
}

# Fit -------------------------------------------------------------

#' Fit Penalized Lee Carter Model
PLCfit <- function(
    Dxt, Ext, label = 'unknown population', W = NULL,
    config = PLCfitConfig()
) {
  
  N = nrow(Dxt)
  m = ncol(Dxt)
  
  age_labels <- rownames(Dxt)
  if (is.null(age_labels)) { age_labels <- as.character(seq_len(N)) }
  
  Mcx <- Dxt/Ext
  Mcx[Mcx==0] <- min(Mcx[Mcx!=0], na.rm = TRUE)
  Mcx[is.nan(Mcx)] <- NA
  
  if (is.null(W)) { W <- PLCcreateWeightMatrix(Dxt = Dxt, Ext = Ext) }
  # information about which kt's are not estimated from data due to
  # the weight for the corresponding year being 0
  kt_zerowght <- apply(W, 2, function (x) all(x == 0))
  
  # initialize parameters
  if (!is.null(config$init_pars)) {
    theta <- config$init_pars
  } else {
    theta <- PLCgetInitEstimates(Mcx, mode = config$init_mode)
  }
  theta <- PLCconstraints(theta$ax, theta$bx, theta$kt)
  
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

  objective_args <- list(
    Dcx = Dxt,
    Ecx = Ext,
    DD_ax = DD_ax,
    DD_bx = DD_bx,
    DD_kt = DD_kt,
    lambda_ax = lambda_ax,
    lambda_bx = lambda_bx,
    lambda_kt = lambda_kt,
    lambda_ridge = lambda_ridge,
    lambda_ax_target = lambda_ax_target,
    lambda_bx_target = lambda_bx_target,
    lambda_kt_target = lambda_kt_target,
    ax_target = config$ax_target,
    bx_target = config$bx_target,
    kt_target = config$kt_target
  )
  
  eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
  epsilon <- log(Mcx)-eta
  
  maxit = config$outer_maxit
  bfgs_maxit = config$bfgs_maxit
  
  dev <- rep(NA, maxit+1)
  dev[1] <- PLCpoissonDeviance(eta, Dxt, Ext, W)
  cat(dev[1], '\n', sep = '')
  r2 <- rep(NA, maxit+1)
  r2[1] <- PLCr2(epsilon, eta, W)
  
  if (isTRUE(config$plot_progress)) {
    PLCplotFitDiagnostics(dev, epsilon, eta, theta, maxit,
                          rejected_bx = rep(FALSE, N),
                          outlier = kt_zerowght,
                          population = label)
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
  final_iter <- 0
  for (i in 1:maxit) {
    final_iter <- i
    
    # optimize ax
    ax_step <- PLCoptim(theta, 'ax', objective_args, W, bfgs_maxit)
    theta <- ax_step$theta
    ax_optim <- ax_step$optim
    # optimize bx
    bx_step <- PLCoptim(theta, 'bx', objective_args, W, bfgs_maxit)
    theta <- bx_step$theta
    bx_optim <- bx_step$optim
    # optimize kt
    kt_step <- PLCoptim(theta, 'kt', objective_args, W, bfgs_maxit)
    theta <- kt_step$theta
    kt_optim <- kt_step$optim
    
    # predict logrates
    eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
   
    # report diagnostics
    epsilon <- log(Mcx)-eta
    dev[i+1] <- PLCpoissonDeviance(eta, Dxt, Ext, W)
    cat(dev[i+1], '\n', sep = '')
    r2[i+1] <- PLCr2(epsilon, eta, W)
    if (isTRUE(config$plot_progress)) {
      PLCplotFitDiagnostics(dev, epsilon, eta, theta, maxit,
                            rejected_bx = rep(FALSE, N),
                            outlier = kt_outlier|kt_zerowght,
                            population = label)
    }
    
    # determine deviance convergence
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
        theta <- config$init_pars
      } else {
        theta <- PLCgetInitEstimates(Mcx, mode = config$init_mode)
      }
      theta <- PLCconstraints(theta$ax, theta$bx, theta$kt)
    }
    
  }

  # convergence determination
  optim_convergence <- all(c(
    ax_optim$convergence,
    bx_optim$convergence,
    kt_optim$convergence
  ) == 0)
  # total convergence is achieved when ax,bx,kt optim calls report
  # convergence and when the fina
  total_convergence <-
    isTRUE(deviance_convergence) && isTRUE(optim_convergence)
  
  # kt adjustment
  if (isTRUE(config$kt_adjustment == 'poisson')) {
    # for the kt adjustment all data is used unless
    # its Na or Ext is 0
    W_complete <- PLCcreateWeightMatrix(Dxt = Dxt, Ext = Ext)
    adjustment_args <- objective_args
    adjustment_args$lambda_kt <- 0
    adjustment_args$lambda_ridge <- 0
    adjustment_args$lambda_kt_target <- 0
    kt_step <- PLCoptim(
      theta, 'kt', adjustment_args, W_complete, maxit = 100
    )
    theta <- kt_step$theta
    eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
    epsilon <- log(Mcx)-eta
    
    if (isTRUE(config$plot_progress)) {
      PLCplotFitDiagnostics(dev, epsilon, eta, theta, maxit,
                            rejected_bx = rep(FALSE, N),
                            outlier = kt_outlier|kt_zerowght,
                            population = label)
    }
  }

  # set bx to 0 if their sensitivity to kt changes is above threshold
  bx_repair_step <- PLCrepairBx(
    theta,
    threshold = config$bx_sensitivity_threshold,
    kt_exclude = kt_outlier | kt_zerowght,
    age_labels = age_labels
  )
  theta <- bx_repair_step$theta
  bx_scaled_sensitivities <- bx_repair_step$bx_scaled_sensitivities
  rejected_bx <- bx_repair_step$rejected_bx
  bx_repair <- bx_repair_step$diagnostics
  if (bx_repair$applied) {
    eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
    epsilon <- log(Mcx)-eta
  }
  if (isTRUE(config$plot_progress) || isTRUE(config$plot_final)) {
    PLCplotFitDiagnostics(
      dev, epsilon, eta, theta, maxit,
      outlier = kt_outlier|kt_zerowght,
      rejected_bx = rejected_bx,
      population = label
    )
  }
  
  # calculate some diagnostics on 0 counts
  included_cells <- !is.na(W) & W > 0
  included_cells_per_age <- rowSums(included_cells)
  positive_count_shares <-
    rowSums(included_cells & Dxt > 0, na.rm = TRUE) /
    included_cells_per_age
  positive_count_shares[included_cells_per_age == 0] <- NA_real_
  positive_count_shares <- setNames(positive_count_shares, age_labels)
  
  return(
    list(
      data = list(
        exposure = Ext,
        deaths = Dxt
      ),
      model_parameters = theta,
      predicted_log_rates = eta,
      jumpoff_calibration_vector = setNames(epsilon[,m], age_labels),
      model_diagnostics = list(
        deviance = dev,
        r2 = r2,
        epsilon = epsilon,
        final_iter = final_iter,
        optim_convergence = optim_convergence,
        total_convergence = total_convergence,
        bx_scaled_sensitivities = bx_scaled_sensitivities,
        positive_count_shares = positive_count_shares,
        bx_repair = bx_repair
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
  if (abs(c2) <= sqrt(.Machine$double.eps)) {
    warning(paste(
      'bx sensitivity repair leaves parameters with a near-zero sum.',
      'LC constraints not applied.'),
      call. = FALSE)
    theta <- list(
      ax = ax,
      bx = bx,
      kt = kt
    )
  } else {
    theta <- list(
      ax = ax + c1 * bx,
      bx = bx / c2,
      kt = c2 * (kt - c1)
    )
  }
  return(theta)
}

# Penalized likelihood --------------------------------------------

#' Lee-Carter Penalized Poisson Log-Likelihood
PLCpoissonLL <- function (
    ax, bx, kt, Dcx, Ecx,
    DD_ax, DD_bx, DD_kt,
    lambda_ax, lambda_bx, lambda_kt, lambda_ridge,
    lambda_ax_target, lambda_bx_target, lambda_kt_target,
    ax_target = 0, bx_target = 0, kt_target = 0,
    W,
    penalties = TRUE
) {
  theta <- PLCconstraints(ax, bx, kt)
  eta <- PLCpredict(theta$ax, theta$bx, theta$kt)
  mu <- exp(eta)
  
  unpenalized_loglike <- sum(W*(Dcx*eta-Ecx*mu), na.rm = TRUE)
  if (isTRUE(penalties)) {
    smootheness_penalties <-
      lambda_ax*t(theta$ax)%*%DD_ax%*%theta$ax +
      lambda_bx*t(theta$bx)%*%DD_bx%*%theta$bx +
      lambda_kt*t(theta$kt)%*%DD_kt%*%theta$kt
    ridge_penalty <-
      lambda_ridge*sum(c(theta$bx, theta$ax)^2)
    target_deviation_penalties <-
      lambda_ax_target*sum((ax-ax_target)^2, na.rm = TRUE) +
      lambda_bx_target*sum((bx-bx_target)^2, na.rm = TRUE) +
      lambda_kt_target*sum((kt-kt_target)^2, na.rm = TRUE)
    unpenalized_loglike -
      smootheness_penalties -
      ridge_penalty -
      target_deviation_penalties
  } else {
    unpenalized_loglike
  }
}

# Diagnostics -----------------------------------------------------

#' Lee-Carter Poisson Deviance
PLCpoissonDeviance <- function (eta, Dcx, Ecx, W) {
  included <-
    !is.na(W) & W != 0 &
    !is.na(Dcx) & !is.na(Ecx) & !is.na(eta)

  d <- Dcx[included]
  mu <- Ecx[included] * exp(eta[included])
  weight <- W[included]
  contribution <- numeric(length(d))

  zero_deaths <- d == 0
  contribution[zero_deaths] <- mu[zero_deaths]

  positive_deaths <- d > 0
  valid_mu <- positive_deaths & is.finite(mu) & mu > 0
  contribution[valid_mu] <-
    d[valid_mu] * log(d[valid_mu] / mu[valid_mu]) -
    d[valid_mu] +
    mu[valid_mu]

  invalid_mu <- positive_deaths & !valid_mu
  contribution[invalid_mu] <- Inf

  2 * sum(weight * contribution)
}

#' Lee-Carter R2
PLCr2 <- function(epsilon, eta, W) {
  rss <- sum(W*epsilon^2, na.rm = TRUE)
  observed <- W*(eta+epsilon)
  tss <- sum((observed-mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- round((1 - rss/tss)*100, 4)
}

#' Plot PLC Fit Diagnostics
PLCplotFitDiagnostics <- function(
    dev, epsilon, eta, theta, maxit, outlier, rejected_bx, population,
    theme = PLCtosTheme
) {
  
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
    mar = theme$mar,
    oma = c(0, 0, 2, 0)
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
  PlotMatrix(eta, type = 'c', clip = c(-14, 0),
             main = 'Estimated log m(x,t)', xlab = 't', ylab = 'x',
             concol = theme$color_scale_cols_log_mx,
             N = theme$color_scale_n_log_mx)
  plot(x = 1:N, y = theta$ax, xlab = 'x', ylab = 'a(x)',
       main = 'a(x) estimates')
  plot(x = 1:N, y = theta$bx, xlab = 'x', ylab = 'b(x)',
       main = 'b(x) estimates', pch = ifelse(rejected_bx, 4, 1))
  plot(x = 1:m, y = theta$kt, xlab = 't', ylab = 'k(t)',
       main = 'k(t) estimates', pch = ifelse(outlier, 4, 1))
  mtext(
    paste(
      population, format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z'), sep = ' — '
    ),
    side = 3, outer = TRUE, col = theme$col.main
  )
}

# Outlier detection -----------------------------------------------

#' Detect Outliers in Lee-Carter kt Vector
PLCdetectKtOutliers <- function (kt, k = 5, z = 3, absolute = TRUE) {
  kt_smooth <- runmed(kt, k = k, endrule = 'constant')
  res <- kt-kt_smooth
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

# distribution of magnitudes of mortality crises
PLCforecast <- function (theta, h, nsim,
                         sd_estimation = 'classic',
                         drift_estimation = 'classic',
                         jumpoff_estimation = 'classic',
                         kt_exclude = NULL,
                         p_crisis = NULL,
                         m_crisis = NULL,
                         Ext_forecast = NULL,
                         jumpchoice = c('fit', 'actual'),
                         jumpoff_calibration_vector = NULL,
                         kt_lookback = c(drift = Inf, sd = Inf)) {

  jumpchoice <- match.arg(jumpchoice)
  names(kt_lookback) <- c('drift', 'sd')
  
  kt <- theta$kt
  kt_lookback_effective <- pmin(kt_lookback, length(kt))
  names(kt_lookback_effective) <- names(kt_lookback)

    # exclude kt elements from variance and drift estimation
  kt[kt_exclude] <- NA
  GetdktWindow <- function(n) {
    na.exclude(diff(tail(kt, n)))
  }
  dkt_drift_window <- GetdktWindow(kt_lookback_effective['drift'])
  dkt_sd_window <- GetdktWindow(kt_lookback_effective['sd'])

  # estimate linear drift of kt term
  kt_drift <- switch(
    drift_estimation,
    classic = mean(dkt_drift_window),
    robust = median(dkt_drift_window)
  )
  # estimate standard deviation of kt innovations
  dkt_sd <- switch(
    sd_estimation,
    classic = sd(dkt_sd_window),
    weighted = sqrt(tail(
      EMVar(dkt_sd_window, length(dkt_sd_window)-1), 1
    )),
    robust = 1.48*median(abs(dkt_sd_window-median(dkt_sd_window)))
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
  if (identical(jumpchoice, 'actual')) {
    Eta_sim <- sweep(Eta_sim, 1, jumpoff_calibration_vector, FUN = '+')
  }

  Y_sim <- NULL
  if (!is.null(Ext_forecast)) {
    expected_dimensions <- c(length(theta$ax), h)
    if (!is.matrix(Ext_forecast) ||
        !all(dim(Ext_forecast) == expected_dimensions)) {
      stop(
        sprintf(
          'Ext_forecast must be a %d x %d age-by-horizon matrix.',
          expected_dimensions[1], expected_dimensions[2]
        ),
        call. = FALSE
      )
    }
    exposure_array <- array(Ext_forecast, dim = dim(Eta_sim))
    expected_deaths <- exposure_array * exp(Eta_sim)
    Y_sim <- array(
      rpois(length(expected_deaths), lambda = as.vector(expected_deaths)),
      dim = dim(Eta_sim)
    )
  }
  
  list(
    Eta_forecast_sim = Eta_sim,
    Y_forecast_sim = Y_sim,
    kt_drift = kt_drift,
    dkt_sd = dkt_sd,
    kt_sim = kt_sim,
    kt_lookback = list(
      requested = kt_lookback,
      effective = kt_lookback_effective
    ),
    jumpchoice = jumpchoice,
    jumpoff_calibration_vector = if (identical(jumpchoice, 'actual')) {
      jumpoff_calibration_vector
    } else {
      NULL
    }
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

# Misc ------------------------------------------------------------

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
