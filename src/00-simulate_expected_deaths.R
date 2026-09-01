# Forecast expected A/S deaths post 2019 by age, sex, and region
#
# Lee-Carter model is used for stochastic forecast of expected deaths
# given pre-pandemic trends.

# Init --------------------------------------------------------------------

library(yaml)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rje)

# Constants ---------------------------------------------------------------

# input and output paths
setwd(".")
paths <- list()
paths$input <- list(
  ASDCounts_AllMar26.csv = "./dat/ASDCounts_AllMar26.csv"
)
paths$output <- list(
  observed_vs_expected.pdf = "./out/observed_vs_expected.pdf"
)

# constants specific to this analysis
cnst <- within(list(), {
  fitting_years = 2000:2019
  forecast_years = 2020:2023
  forecast_horizon = length(forecast_years)
  nsim = 100
  seed = 1987
})

# list containers for analysis artifacts
dat <- list()

set.seed(cnst$seed)

# Functions ---------------------------------------------------------------

source('src/00-penalized_lee_carter_functions.R')

# Import ------------------------------------------------------------------

asd <- read_csv(paths$input$ASDCounts_AllMar26.csv, col_select = -1)

# Ensure uniform dimensions -----------------------------------------------

skeleton <- expand_grid(
  Definition = unique(asd$Definition),
  Country = unique(asd$Country),
  Sex = unique(asd$Sex),
  Year = c(cnst$fitting_years, cnst$forecast_years),
  Ageband = unique(asd$Ageband),
)

asd <- left_join(skeleton, asd)

# Subset to data of interest ----------------------------------------------

asd_sub <-
  asd |>
  filter(
    Definition == 'Full',
  ) |>
  select(-Definition) |>
  # convert to integer to adhere with Poisson likelihood
  # round to closest integer
  mutate(Deaths = as.integer(round(Deaths, 0)))

# Predict expected deaths -------------------------------------------------

# convert to wide format, "matrix-like", df's separate for counts and exposures
# to use with Lee-Carter

# foo <-
#   asd_sub |>
#   nest(.by = c(Sex, Country)) |>
#   filter(Sex == 'Male', Country == 'Argentina') |>
#   unnest()
#
# .y <- foo |> select(Sex, Country) |> unique()
# .x <- foo |> select(Ageband, Year, Deaths, Pop)

expected_vs_observed_sim <-
  asd_sub |>
  group_by(Country, Sex) |>
  group_modify(~{

    cat(format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
        ' Fit stratum ', paste0(.y, collapse = '-'), '\n', sep = '')

    # dimensions
    period_labels <- unique(.x$Year)
    period_labels_numeric <- as.integer(period_labels)
    period_labels_character <- as.character(period_labels)
    period_length <- length(period_labels)
    age_labels <- unique(.x$Ageband)
    age_labels_numeric <- as.integer(substr(age_labels, 1, 2))
    age_labels_character <- as.character(age_labels)
    age_length <- length(age_labels)
    period_labels_forecast_character <-
      period_labels_character[period_labels_numeric %in% cnst$forecast_years]
    H <- cnst$forecast_horizon

    # ensure uniform output
    output_skeleton <-
      expand_grid(
        period = c(period_labels_character),
        age = age_labels_character,
        sim = as.character(1:cnst$nsim)
      )

    # fit models and capture errors
    result <- tryCatch(
      { # normal fit

        # lexis matrices of counts and exposures
        Y <- matrix(
          .x$Deaths,
          nrow = age_length, ncol = period_length,
          dimnames = list(age = age_labels, period = period_labels_character)
        )
        E <- matrix(
          .x$Pop,
          nrow = age_length, ncol = period_length,
          dimnames = list(age = age_labels, period = period_labels_character)
        )

        # subset to fitting years
        Y_fit <- Y[,as.character(cnst$fitting_years)]
        E_fit <- E[,as.character(cnst$fitting_years)]

        # subset to forecast years
        E_fcst <- E[,period_labels_forecast_character]
        Y_fcst <- Y[,period_labels_forecast_character]

        # fit LC model
        LC_fit <- PLCfit(
          Dxt = Y_fit, Ext = E_fit,
          label = paste0(.y, collapse = '-'),
          config = PLCfitConfig(
            init_mode = 'svd',
            bx_sensitivity_threshold = 0.2,
            lambda_ridge = 0
          )
        )

        # make forecast over excess period
        LC_sim <- PLCforecast(
          theta = LC_fit$model_parameters,
          h = H, Ext_forecast = E_fcst,
          jumpchoice = 'actual', sd_estimation = 'classic',
          jumpoff_calibration_vector = LC_fit$jumpoff_calibration_vector,
          kt_lookback = c(10, 10),
          nsim = cnst$nsim
        )

        # derive simulated counts with added Poisson variability over excess period
        expected_counts_sim <- LC_sim$Y_forecast_sim

        dimnames(expected_counts_sim) <- list(
          age = age_labels_character,
          period = period_labels_forecast_character,
          sim = as.character(1:cnst$nsim)
        )

        # average expected counts over excess period
        expected_counts_avg <-
          apply(expected_counts_sim, 1:2, mean)
        dimnames(expected_counts_avg) <- list(
          age = age_labels_character,
          period = period_labels_forecast_character
        )

        # assemble long format output data frame with observed and expected simulations
        observed_df <- bind_rows(
          array2DF(Y_fit, responseName = 'OBS'),
          array2DF(Y_fcst, responseName = 'OBS')
        )
        expected_avg_df <-
          bind_rows(
            array2DF(expected_counts_avg, responseName = 'XPC_AVG')
          )
        expected_sim_df <-
          bind_rows(
            array2DF(expected_counts_sim, responseName = 'XPC_SIM')
          )

        simulation_skeleton <-
          expand_grid(
            period = period_labels_forecast_character,
            age = age_labels_character,
            sim = as.character(1:cnst$nsim)
          )

        result_if_no_error <-
          output_skeleton |>
          left_join(observed_df, by = c('period', 'age')) |>
          left_join(expected_avg_df, by = c('period', 'age')) |>
          left_join(expected_sim_df, by = c('period', 'age', 'sim')) |>
          mutate(
            stratum = paste0(.y, collapse = '-'),
            period = as.numeric(period)
          )

        return(result_if_no_error)
      },

      error = function (e) { # return in case of error

        cat(format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
            ' Error in stratum', paste0(.y, collapse = '-'),
            ' : ', geterrmessage(), '\n')

        result_if_error <-
          output_skeleton |>
          mutate(OBS = NA, XPC_AVG = NA, XPC_SIM = NA, stratum = paste0(.y, collapse = '-')) |>
          mutate(period = as.numeric(period))

        return(result_if_error)
      }
    )

  })

# Plot observed vs. expected by region and sex ----------------------------

pdf(file=paths$output$observed_vs_expected.pdf)
for (i in unique(asd_sub$Country)) {
  cat("Plot", i, "\n")
  the_plot <-
    expected_vs_observed_sim |>
    filter(Country == i) |>
    ggplot() +
    geom_line(
      aes(x = period, y = XPC_AVG, color = Sex),
      data = . %>% filter(sim == 1)
    ) +
    geom_line(
      aes(x = period, y = XPC_SIM, group = interaction(sim,Sex), color = Sex),
      alpha = 0.1
    ) +
    geom_point(
      aes(x = period, y = OBS, color = Sex),
      data = . %>% filter(sim == 1)
    ) +
    geom_vline(xintercept = 2019.5, color = "grey50") +
    facet_wrap(~age, scale = 'free_y') +
    labs(title = i, y = "A/D deaths", x = 'Year')
  print(the_plot)
}
dev.off()

# Demonstration -----------------------------------------------------------

expected_vs_observed_sim |>
  filter(period %in% 2020:2023) |>
  mutate(excess = OBS - XPC_SIM) |>
  group_by(sim) |>
  summarise(
    excess = sum(excess, na.rm = T)
  ) |>
  ungroup() |>
  summarise(
    excess_med = median(excess),
    excess_avg = mean(excess),
    excess_q025 = quantile(excess, 0.025),
    excess_q975 = quantile(excess, 0.975)
  )

expected_vs_observed_sim |>
  filter(
    period %in% 2020:2023
  ) |>
  group_by(Country) |>
  summarise(
    any_na_in_expected_sims = anyNA(XPC_SIM)
  ) |> View()

# Export ------------------------------------------------------------------

# export results of analysis
