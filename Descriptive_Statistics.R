#################################################################################################################
###   Forecasting Tail Risk in the Brazilian Stock Market: An Out-of-Sample Comparison of VaR and ES Methods  ###
###   Authors: Esther Cleveston and Carlos Trucíos
#################################################################################################################
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(moments)

# Import Data
data <- read_excel("precos_b3.xlsx", skip = 3, col_types = c("date", rep("numeric", 475)), na = c("", "-", NA)) |> 
  filter(Data > "2010-01-01" & Data < "2026-10-01") |> 
  filter(!if_all(where(is.numeric), is.na)) |> 
  select_if(~ !any(is.na(.))) |> 
  rename_with(~ str_remove(.x, "^(?s).*orig\n"), -Data) |> 
  mutate(across(where(is.numeric), ~ (log(.x) - lag(log(.x))) * 100)) |> 
  slice(-1)
names_stocks <- colnames(data)[-1]

# Table
descriptive_statistics <- function(x) {
  round(c(min(x), 
    quantile(x, 0.25),
    mean(x),
    median(x),
    quantile(x, 0.75), 
    max(x),
    sd(x),
    skewness(x),
    kurtosis(x)), 3)
}


xtable::xtable(t(apply(data[, -1], 2, descriptive_statistics)))

# Plot

data |>
  pivot_longer(
    cols = -Data,
    names_to = "Stocks",
    values_to = "Returns"
  ) |>
  ggplot(aes(x = Data, y = Returns, color = Stocks)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2020-02-12"), linetype = "dashed") +
  labs(
    x = NULL,
    y = "log-returns (in percentage)",
    color = "Stocks"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
