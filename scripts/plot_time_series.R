# generate time series plots of cases and deaths

# config
source(here::here("depends.R"))
ggplot2::theme_set(ggplot2::theme_bw())
output_path <- here::here("outputs")

# load data ----
molis_deaths <- readr::read_csv(
  here::here("data/molis_deaths_synthetic.csv")
)

sgss_deaths <- readr::read_csv(
  here::here("data/sgss_deaths_synthetic.csv")
)

# Tseries plot ------------

# filter data between these dates
min_date <- as.Date("2022-07-01")
max_date <- as.Date("2025-04-01")

## Test count time series ----

full_dates <- data.frame(
  date = seq(from = min_date,
             to = max_date,
             by = "day"))

molis_test_count <- molis_deaths |>
  dplyr::summarise(count = dplyr::n(), .by = "specimen_date") |>
  dplyr::right_join(full_dates,
                    by = dplyr::join_by(specimen_date == date)) |>
  dplyr::mutate(count = dplyr::coalesce(count, 0L))

sgss_test_count <- sgss_deaths |>
  dplyr::summarise(count = dplyr::n(), .by = "specimen_date") |>
  dplyr::right_join(full_dates,
                    by = dplyr::join_by(specimen_date == date)) |>
  dplyr::mutate(count = dplyr::coalesce(count, 0L))

sgss_molis_tseries <- ggplot2::ggplot() +
  ggplot2::geom_line(data = sgss_test_count,
                     ggplot2::aes(x = specimen_date, y = count, color = "SGSS")) +
  ggplot2::geom_line(data = molis_test_count,
                     ggplot2::aes(x = specimen_date, y = count, color = "MOLIS")) +
  ggplot2::labs(
    tag = "a)",
    x = "Specimen date",
    y = "Number of tests",
    color = "Data source") +
  ggplot2::scale_color_discrete(breaks = c("SGSS", "MOLIS")) +
  ggplot2::scale_x_date(labels = scales::label_date_short(),
                        breaks = seq(min_date, max_date + 1, by = "3 months"),
                        limits = c(min_date, max_date)) +
  ggplot2::theme(legend.position = "none")

## MOLIS proportion time series ----

molis_tseries_counts <- molis_deaths |>
  dplyr::mutate(specimen_week = lubridate::floor_date(specimen_date, unit = "week")) |>
  dplyr::count(specimen_week, genotype)

genogroup_order <- c("GII.4", "Other", "GII.17")

molis_tseries_counts_complete <- expand.grid(
  specimen_week = seq(
    min(molis_tseries_counts$specimen_week, na.rm = TRUE),
    max(molis_tseries_counts$specimen_week, na.rm = TRUE),
    by = "week"),
  genotype = unique(molis_tseries_counts$genotype)
) |>
  dplyr::left_join(molis_tseries_counts,
                   by = c("specimen_week", "genotype")) |>
  dplyr::mutate(n = dplyr::coalesce(n, 0),
                genotype = factor(
                  genotype,
                  levels = genogroup_order
                ))

# Proportion plot
molis_genotype_prop <- molis_tseries_counts_complete |>
  dplyr::filter(genotype != "not_detected",
                specimen_week >= "2022-06-12") |>
  dplyr::mutate(prop = n / sum(n),
                .by = specimen_week) |>
  ggplot2::ggplot(ggplot2::aes(x = specimen_week, y = prop, fill = genotype)) +
  ggplot2::geom_area(position = "stack") +
  ggplot2::scale_fill_brewer(palette = "Set2", name = "Genotype") +
  ggplot2::labs(
    tag = "c)",
    y = "Proportion",
    x = "Specimen week"
  ) +
  ggplot2::scale_x_date(labels = scales::label_date_short(),
                        breaks = seq(min_date, max_date + 1, by = "3 months"),
                        limits = c(min_date, max_date))


## deaths time series ----

sgss_deaths_28days <- sgss_deaths |>
  dplyr::filter(death_date - specimen_date <= 28,
                death_date - specimen_date >= 0) |>
  # deduplicate
  dplyr::slice_max(specimen_date,
                   by = uid)

sgss_deaths_28days_count <- sgss_deaths_28days |>
  dplyr::summarise(count = dplyr::n(), .by = "death_date") |>
  dplyr::right_join(full_dates,
                    by = dplyr::join_by(death_date == date)) |>
  dplyr::mutate(count = dplyr::coalesce(count, 0L))

molis_deaths_28days <- molis_deaths |>
  dplyr::filter(death_date - specimen_date <= 28,
                death_date - specimen_date >= 0) |>
  # deduplicate
  dplyr::slice_max(specimen_date,
                   by = uid)

molis_deaths_28days_count <- molis_deaths_28days |>
  dplyr::summarise(count = dplyr::n(), .by = "death_date") |>
  dplyr::right_join(full_dates,
                    by = dplyr::join_by(death_date == date)) |>
  dplyr::mutate(count = dplyr::coalesce(count, 0L))


deaths_tseries <- ggplot2::ggplot() +
  ggplot2::geom_line(data = sgss_deaths_28days_count,
                     ggplot2::aes(x = death_date, y = count, color = "SGSS")) +
  ggplot2::geom_line(data = molis_deaths_28days_count,
                     ggplot2::aes(x = death_date, y = count, color = "MOLIS")) +
  ggplot2::labs(
    tag = "b)",
    x = "Date of death",
    y = "Number of deaths",
    color = "Data source") +
  ggplot2::scale_color_discrete(breaks = c("SGSS", "MOLIS")) +
  ggplot2::scale_x_date(labels = scales::label_date_short(),
                        breaks = seq(min_date, max_date + 1, by = "3 months"),
                        limits = c(min_date, max_date)) +
  ggplot2::lims(y = c(0, 15)) +
  ggplot2::theme(legend.position = "none")

tseries_plot <- sgss_molis_tseries / deaths_tseries / molis_genotype_prop
tseries_plot
ggplot2::ggsave(filename = "tseries_plot.png",
                plot = tseries_plot,
                path = output_path,
                width = 10,
                height = 10)
