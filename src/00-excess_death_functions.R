# Misc ------------------------------------------------------------

# Ensure that the xcod object is sorted
EnsureSorted <- function (X) {
  X[order(xcod_out[['stratum']],
                 xcod_out[['origin_time']],
                 xcod_out[['age']], decreasing = FALSE),]
}

# Return data frame of row-wise quantiles over columns of X
Rowquantiles <- function (X, prob, type = 4, na.rm = TRUE) {
  t(apply(X, 1, quantile, prob = prob, type = type, na.rm = na.rm, names = FALSE))
}
# Return data frame of row-wise p-values: P(X_i>=x_i), were X_i are
# simulations of the test statistic under H0 and x_i is the observed
# test statistic
RowPvalue <- function (X, x, na.rm = TRUE, twotailed = FALSE) {
  X_ <- cbind(x, X)
  if (isTRUE(twotailed)) X_ <- abs(X_)
  apply(X_, 1, function (y) sum(sort(y[-1])>=y[1])/(length(y)-1))
}

GetExcess <- function (
    plc,
    measure = 'absolute',
    quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975),
    cumulative = FALSE, origin_time_start_of_cumulation = 0,
    twotailed = FALSE
) {

  # forgiveness please, this info should really be explicit in xcod_out
  nsim = max(as.integer(
    sub('(^.+_SIM)([[:digit:]]+)(.*$)','\\2',
        grep('[[:digit:]]',names(xcod_out), value = TRUE))
  ))
  nrow = NROW(xcod_out)
  xcod_sorted <- EnsureSorted(xcod_out)

  # data columns
  Y <- xcod_sorted[,grepl('XPC_|OBS_',colnames(xcod_sorted))]
  # label columns
  X <- xcod_sorted[,c('stratum', 'origin_time', 'age', 'cv_flag')]

  # aggregate parts if requested
  if (is.list(name_parts)) {
    amalgamation_names <- names(name_parts)
    amalgamation_parts <- name_parts
    Y <- Map(function (x, y) {
      name <- x
      parts <- unlist(y)
      n_parts <- length(parts)

      # OBS
      OBS_colnames <- grep(paste0('OBS_', parts, collapse = '|'),
                           colnames(Y), value = TRUE)
      OBS_parts <- matrix(unlist(Y[,OBS_colnames]),
                          nrow = nrow, ncol = n_parts)
      OBS_amalgamation <- rowSums(OBS_parts)
      # XPC AVG
      XPC_AVG_colnames <- grep(paste0('XPC_AVG_', parts, collapse = '|'),
                               colnames(Y), value = TRUE)
      XPC_AVG_parts <- matrix(unlist(Y[,XPC_AVG_colnames]),
                              nrow = nrow, ncol = n_parts)
      XPC_AVG_amalgamation <- rowSums(XPC_AVG_parts)
      # XPC SIM
      XPC_SIM_colnames <- grep(paste0('XPC_SIM.+_', parts, collapse = '|'),
                               colnames(Y), value = TRUE)
      XPC_SIM_parts <- array(
        unlist(Y[,XPC_SIM_colnames]),
        dim = c(nrow, nsim, n_parts)
      )
      XPC_SIM_amalgamation <- apply(XPC_SIM_parts, 1:2, sum)
      # XCOD amalgamation
      xcod_amalgamation <- cbind(
        OBS_amalgamation, XPC_AVG_amalgamation, XPC_SIM_amalgamation
      )
      colnames(xcod_amalgamation) <-
        c(paste0(c('OBS_', 'XPC_AVG_'), name),
          paste0('XPC_SIM', 1:nsim, '_', name))

      return(xcod_amalgamation)
    },
    amalgamation_names, amalgamation_parts)

    # reorder columns
    Y <- do.call('cbind', Y)
    Y <- Y[,c(which(grepl('OBS',colnames(Y))),
              which(grepl('XPC_AVG',colnames(Y))),
              which(grepl('XPC_SIM',colnames(Y))))]
  } else {
    amalgamation_names <- name_parts
  }

  # accumulate data columns if requested
  if (isTRUE(cumulative)) {
    # set data to 0 for time points prior to accumulation start
    # so that we only accumulate from the accumulation start
    vec_timeselect <- X[['origin_time']] < origin_time_start_of_cumulation
    Y[vec_timeselect,] <- 0
    # this restarts the accumulation of a data vector whenever the
    # stratum vector changes in value
    vec_stratum <- X[['stratum']]
    Y <- apply(Y, 2, function (x) ave(x, vec_stratum, FUN = cumsum))
    # set values to NA prior to accumulation start so that derived
    # values become NA as well
    Y[vec_timeselect,] <- NA
  }

  # calculate excess measures, prediction intervals, and associated P-values
  for (part in amalgamation_names) {
    OBS <- Y[,paste0('OBS_', part)]
    AVG <- Y[,paste0('XPC_AVG_', part)]
    j <- grepl(paste0('^XPC_SIM[[:digit:]]+_',part,'$'), colnames(Y))
    if (identical(measure, 'observed')) {
      MEASURE <- as.matrix(OBS)
    }
    if (identical(measure, 'expected')) {
      MEASURE <- Y[,j]
    }
    if (identical(measure, 'absolute')) {
      MEASURE <- apply(Y[,j], 2, function (XPC_SIM) {round(OBS-XPC_SIM,0)})
      # for p-values: distribution of the test statistic under the null
      # hypothesis of expected distribution of deaths
      H0DIST <- apply(Y[,j], 2, function (XPC_SIM) {round(XPC_SIM-AVG,0)})
      # test statistic
      TESTSTAT <- OBS-AVG
    }
    if (identical(measure, 'pscore')) {
      MEASURE <- apply(Y[,j], 2, function (XPC_SIM) {(OBS-XPC_SIM)/XPC_SIM*100})
      H0DIST <- apply(Y[,j], 2, function (XPC_SIM) {(XPC_SIM-AVG)/AVG*100})
      TESTSTAT <- (OBS-AVG)/AVG*100
    }
    if (identical(measure, 'ratio')) {
      MEASURE <- apply(Y[,j], 2, function (XPC_SIM) {OBS/XPC_SIM})
      H0DIST <- apply(Y[,j], 2, function (XPC_SIM) {XPC_SIM/AVG})
      TESTSTAT <- OBS/AVG
    }
    # get quantiles
    Q <- Rowquantiles(MEASURE, quantiles, type = 1)
    # get p-values
    if (identical(measure, 'observed') | identical(measure, 'expected')) {
      Pval <- NA
    } else {
      Pval <- RowPvalue(H0DIST, TESTSTAT, twotailed = twotailed)
    }
    QP <- cbind(Q,Pval)

    # set names
    colnames(QP) <-
      c(
        paste0('Q', substr(formatC(quantiles, format = 'f', digits = 3),
                           start = 3, stop = 5), '_', part),
        paste0('pvalue', '_', part)
      )

    # bind results
    X <- cbind(X,QP)
  }
  return(X)
}
