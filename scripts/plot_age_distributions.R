# plot age distributions of norovirus tests and deaths

# config
source(here::here("depends.R"))
ggplot2::theme_set(ggplot2::theme_bw())
output_path <- here::here("outputs")

# load data ----

# filter data between these dates
min_date <- as.Date("2022-07-01")
max_date <- as.Date("2025-04-01")

molis_deaths <- readr::read_csv(
  here::here("data/molis_deaths_synthetic.csv")
) |>
  dplyr::filter(specimen_date >= min_date,
                specimen_date < max_date)

sgss_deaths <- readr::read_csv(
  here::here("data/sgss_deaths_synthetic.csv")
) |>
  dplyr::filter(specimen_date >= min_date,
                specimen_date < max_date)

# plot age distributions ---------

sgss_age_dist <- sgss_deaths |>
  ggplot(aes(x = age)) +
  geom_histogram(aes(fill = "All tests"), binwidth = 1, alpha = 0.7) +
  geom_histogram(data = sgss_deaths |>
                   dplyr::filter(death_date - specimen_date <= 28),
                 aes(fill = "Tests with \nlinked death"), binwidth = 1, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "All tests" = "#F8766D",
      "Tests with \nlinked death" = "#C74336"
    )
  ) +
  labs(tag = "a)",
       fill = "SGSS tests",
       x = "Age at specimen date",
       y = "Count")

molis_age_dist <- molis_deaths |>
  ggplot(aes(x = age)) +
  geom_histogram(aes(fill = "All tests"), binwidth = 1, alpha = 0.7) +
  geom_histogram(data = molis_deaths |>
                   dplyr::filter(death_date - specimen_date <= 28),
                 aes(fill = "Tests with \nlinked death"), binwidth = 1, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "All tests" = "#00BFC4",
      "Tests with \nlinked death" = "#00828A"
    )
  ) +
  labs(tag = "b)",
       fill = "MOLIS tests",
       x = "Age at specimen date",
       y = "Count")

# combine
age_dist_plot <- sgss_age_dist / molis_age_dist
age_dist_plot

ggsave(filename = "age_dist_plot.png",
       plot = age_dist_plot,
       path = output_path,
       width = 8,
       height = 6)
