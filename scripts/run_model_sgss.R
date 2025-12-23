# Fit case fatality model to SGSS data and produce outputs

# setup ---------

source(here::here("depends.R"))

# plot theme
ggplot2::theme_set(ggplot2::theme_bw())

# set options ----------

# associate deaths with positive test if within x days
death_threshold <- 28

# maximum death date in deaths_full
latest_observed_death_date <- as.Date("2025-04-29")

# reference start date
ref_date <- as.Date("2022-07-01")

# set output path
output_path <- here::here("outputs",
                          paste0("death_threshold_", death_threshold, "_results"))
fs::dir_create(output_path, recurse = TRUE)
# load and process data -------------

# load linked data
sgss_deaths_deduped <- readr::read_csv(
  here::here("data/sgss_deaths_synthetic.csv")
)

# only associate deaths with norovirus if within threshold after specimen date
sgss_deaths_filtered <- sgss_deaths_deduped |>
  # excluding period at beginning of dataset with sparse data, and last few days to have complete months
  dplyr::filter(specimen_date >= as.Date("2022-07-01"),
                specimen_date < as.Date("2025-04-01")) |>
  # set death dates to NA if further away than death_threshold days
  dplyr::mutate(
    death_date = ifelse(death_date - specimen_date <= death_threshold,
                        death_date,
                        NA))


requesting_org_order <- c("Primary care", "Secondary care", "Other", "Unknown")

sgss_deaths_for_pamm <- sgss_deaths_filtered |>
  dplyr::filter(!is.na(age)) |>
  dplyr::mutate(
    # factor variables
    sgss_lab_region = as.factor(sgss_lab_region),
    sgss_lab_region = relevel(
      sgss_lab_region, ref = "North West"), # consistent with other models

    requesting_organisation = as.factor(requesting_organisation),
    requesting_organisation = relevel(
      requesting_organisation, ref = "Secondary care"), # largest group

    # specimen date coarse groupings
    specimen_quarter = lubridate::floor_date(specimen_date, unit = "quarter"),
    specimen_quarter = as.factor(specimen_quarter)
  ) |>
  dplyr::arrange(uid, specimen_date) |>
  dplyr::group_by(uid) |>
  dplyr::mutate(
    start_time = 0, # define as time since positive test
    # risk interval ends at the earliest of:
    stop_date = pmin(specimen_date + death_threshold, # death censoring date
                     dplyr::lead(specimen_date), # next test
                     latest_observed_death_date, # administrative censoring date
                     death_date,
                     na.rm = TRUE),
    stop_time = as.numeric(stop_date - specimen_date),
    # fatality = 1 if death occurred and stop date equals death date
    fatality = ifelse(!is.na(death_date) &
                        stop_date == death_date, 1, 0),
    # if stop time <= 0, treat as 1
    stop_time = ifelse(stop_time <= 0, 1, stop_time)
  ) |>
  dplyr::ungroup()



# fit model ---------

# Choose cutpoints for time intervals
cut_vec <- seq(0, death_threshold, by = 2)    # every 2 days for speed

# Build PED (piecewise exponential data)
pamm_data <- as_ped(
  data = sgss_deaths_for_pamm,
  formula = Surv(start_time, stop_time, fatality) ~
    age + requesting_organisation + specimen_quarter + uid + sgss_lab_region,
  cut = cut_vec
)

sgss_pamm_fit <- bam(
  ped_status ~
    specimen_quarter +
    s(tend, by = specimen_quarter, k = 5) +
    sgss_lab_region +
    s(age, k = 4) +
    s(age, by = requesting_organisation, k = 5) +
    requesting_organisation +
    s(uid, bs = "re"),
  data = pamm_data,
  family = poisson(),
  offset = offset,
  method = "fREML", # for large models with random effects
  discrete = TRUE,
  nthreads = 4
)

# get estimates and predictions ----------

time_period_labels <- list(
  "2023-01-01" = "2022/23\n(GII.4 dominant)",
  "2024-01-01" = "2023/24\n(GII.4 dominant)",
  "2024-04-01" = "2023/24\n(GII.17 dominant)",
  "2025-01-01" = "2024/25\n(GII.17 dominant)"
)

# get hazard ratios and survival probabilites
pamm_estimates <- pamm_data |>
  make_newdata(
    tend = unique(pamm_data$tend),
    age = seq(0, 100),
    sgss_lab_region = unique(pamm_data$sgss_lab_region),
    requesting_organisation = c("Primary care", "Secondary care", "Other"),
    specimen_quarter = names(time_period_labels)) |>
  group_by(age, specimen_quarter, requesting_organisation, sgss_lab_region) |>
  add_surv_prob(sgss_pamm_fit) |>
  add_hazard(sgss_pamm_fit, reference = list(requesting_organisation = c("Secondary care"),
                                             sgss_lab_region = c("North West"))) |>
  dplyr::filter(tend == death_threshold) |>
  dplyr::mutate(
    time_period = time_period_labels[[specimen_quarter]],
    time_period = factor(time_period, levels = unlist(time_period_labels)),
    requesting_organisation = factor(requesting_organisation,
                                     levels = requesting_org_order),
    cfr = 1 - surv_prob,
    cfr_lower = 1 - surv_upper,
    cfr_upper = 1 - surv_lower) |>
  dplyr::ungroup()

# plot HRs -------


request_org_hr_plot <- pamm_estimates |>
  dplyr::filter(specimen_quarter == "2023-01-01",
                sgss_lab_region == "North West") |>
  ggplot(aes(x = age, group = requesting_organisation)) +
  geom_line(aes(y = hazard, color = requesting_organisation)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = requesting_organisation), alpha = 0.2) +
  labs(x = "Age (years)",
       y = "Hazard ratio",
       color = "Healthcare level",
       fill = "Healthcare level") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")

region_hr_plot <- pamm_estimates |>
  dplyr::mutate(reference = ifelse(sgss_lab_region == "North West", TRUE, FALSE)) |>
  # fix other variables
  dplyr::filter(specimen_quarter == "2023-01-01",
                requesting_organisation == "Secondary care") |>
  ggplot(aes(x = sgss_lab_region)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(y = hazard, color = reference), size = 3) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper)) +
  labs(x = "Lab region",
       y = "Hazard ratio") +
  ggplot2::scale_color_manual(values = c("black", "salmon")) +
  ggplot2::theme(legend.position = "none")

hr_plot <- request_org_hr_plot / region_hr_plot

ggsave(filename = paste0("pamm_sgss_hr_", death_threshold, ".png"),
       plot = hr_plot,
       path = output_path,
       width = 14,
       height = 8,
       dpi = 300)

# plot CFRs ------
cfr_plot <- pamm_estimates |>
  dplyr::filter(sgss_lab_region == "North West") |>
  ggplot(aes(x = age, group = requesting_organisation)) +
  geom_line(aes(y = cfr, color = requesting_organisation)) +
  geom_ribbon(aes(ymin = cfr_upper, ymax = cfr_lower, fill = requesting_organisation), alpha = 0.2) +
  facet_wrap(~time_period, nrow = 1) +
  labs(x = "Age (years)",
       y = paste0(death_threshold, "-day CFR"),
       color = "Healthcare level",
       fill = "Healthcare level") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")

ggsave(filename = paste0("pamm_sgss_cfr_", death_threshold, ".png"),
       plot = cfr_plot,
       path = output_path,
       width = 8,
       height = 4,
       dpi = 300)


# save table --------

# HR
hr_at_median_age <- pamm_estimates |>
  dplyr::filter(specimen_quarter == "2023-01-01",
                sgss_lab_region == "North West",
                requesting_organisation %in% c("Primary care", "Other"),
                age == median_age) |>
  dplyr::select(requesting_organisation, hazard, ci_lower, ci_upper)
hr_at_median_age

# CFR

median_age <- median(pamm_data$age)
modal_region <- names(sort(-table(pamm_data$sgss_lab_region)))[1]

cfr_table <- pamm_estimates |>
  dplyr::filter(
    age == median_age,
    sgss_lab_region == modal_region,
    requesting_organisation %in% c("Secondary care", "Primary care")) |>
  dplyr::mutate(
    cfr = round(cfr * 100, digits = 2) |>
      format(digits = 3),
    cfr_lower = round(cfr_lower * 100, digits = 2) |>
      format(digits = 3),
    cfr_upper = round(cfr_upper * 100, digits = 2) |>
      format(digits = 3),
    cfr_string = glue::glue("{cfr} ({cfr_lower}, {cfr_upper})")
  ) |>
  dplyr::arrange(specimen_quarter) |>
  dplyr::select(
    `Healthcare level` = requesting_organisation,
    `CFR % (95% CI)` = cfr_string) |>
  gt::gt() |>
  gt::tab_row_group(
    label = "2022/23 (GII.4 dominant)",
    rows = 1:2
  ) |>
  gt::tab_row_group(
    label = "2023/24 (GII.4 dominant)",
    rows = 3:4
  ) |>
  gt::tab_row_group(
    label = "2023/24 (GII.17 dominant)",
    rows = 5:6
  ) |>
  gt::tab_row_group(
    label = "2024/25 (GII.17 dominant)",
    rows = 7:8
  ) |>
  gt::row_group_order(
    groups = c("2022/23 (GII.4 dominant)",
               "2023/24 (GII.4 dominant)",
               "2023/24 (GII.17 dominant)",
               "2024/25 (GII.17 dominant)")) |>
  gt::tab_options(column_labels.font.weight = "bold",
                  row_group.font.weight = "bold")

cfr_table |>
  gt::gtsave(filename = paste0("pamm_sgss_cfr_table_", death_threshold, ".png"),
             path = output_path)
