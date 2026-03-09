# generate plots and tables relating to linked outcomes:
# - reporting delay for hospitalisation/deaths
# - summary tables of linked hospitalisation/deaths by reporting delay, genotype, age

# config
source(here::here("depends.R"))
output_path <- here::here("outputs")
ggplot2::theme_set(ggplot2::theme_bw())

# load data ----

# linked deaths data

molis_deaths <- readr::read_csv(
  here::here("data/molis_deaths_synthetic.csv")
) |>
  dplyr::mutate(specimen_to_death_days = as.numeric(death_date - specimen_date))

sgss_deaths <- readr::read_csv(
  here::here("data/sgss_deaths_synthetic.csv")
) |>
  dplyr::mutate(specimen_to_death_days = as.numeric(death_date - specimen_date))

molis_sgss_deaths <- readr::read_csv(
  here::here("data/molis_sgss_deaths_synthetic.csv")
) |>
  dplyr::mutate(specimen_to_death_days = as.numeric(death_date - specimen_date))

# linked hospitalisations data

molis_hosp_arrival <- readr::read_csv(
  here::here("data/molis_hosp_arrival_synthetic.csv")
)

molis_hosp_noro_episode <- readr::read_csv(
  here::here("data/molis_hosp_noro_episode_synthetic.csv")
)

# time to outcome histograms ------------------

# time to hospital spell admission

# Plot histogram
t_to_hospitalisation_plot <- molis_hosp_arrival |>
  ggplot2::ggplot(ggplot2::aes(x = specimen_to_arrival_days, fill = "MOLIS")) +
  ggplot2::geom_vline(xintercept = 0) +
  ggplot2::geom_histogram(binwidth = 1, color = "black") +
  ggplot2::labs(
    tag = "a)",
    x = "Days from specimen date to hospital spell admission",
    y = "Number of admissions",
  ) +
  ggplot2::coord_cartesian(x = c(-50, 50), y = c(-5, 220), expand = FALSE) +
  ggplot2::theme(legend.position = "none")

# time to death

# Plot histogram
t_to_death_plot <- molis_deaths |>
  ggplot2::ggplot(ggplot2::aes(x = specimen_to_death_days, fill = "MOLIS")) +
  ggplot2::geom_vline(xintercept = 0) +
  ggplot2::geom_histogram(binwidth = 1, color = "black") +
  ggplot2::labs(
    tag = "b)",
    x = "Days from specimen date to death",
    y = "Number of deaths",
  ) +
  ggplot2::coord_cartesian(x = c(-12, 150), y = c(-0.2, 12), expand = FALSE) +
  ggplot2::theme(legend.position = "none")

t_to_outcome_plot <- t_to_hospitalisation_plot / t_to_death_plot
t_to_outcome_plot

ggplot2::ggsave(filename = "t_to_outcome_plot.png", plot = t_to_outcome_plot, path = output_path, width = 8, height = 6)

# linkage to outcome summary table ------------

n_tests <- nrow(molis_deaths)

## overall test count ----------

overall_test_count <- data.frame(
  `Days from specimen to outcome date` = "",
  `Count` = n_tests,
  `Proportion of positive tests` = "",
  check.names = FALSE
)

## hospital spells ------------

hosp_spell_linkage_count <- molis_hosp_arrival |>
  dplyr::mutate(
    # Categorise time to hospital admission
    t_to_outcome_group = cut(
      specimen_to_arrival_days,
      breaks = c(-28, 1, 15, 29),
      labels = c("[-28, 0]", "[1,14]", "[15,28]"),
      right = FALSE
    )
  ) |>
  dplyr::filter(!is.na(t_to_outcome_group)) |>
  dplyr::summarise(
    `Count` = dplyr::n(),
    .by = t_to_outcome_group
  ) |>
  dplyr::mutate(
    `Proportion of positive tests` = round(`Count` / n_tests, 3)
  ) |>
  dplyr::arrange(t_to_outcome_group) |>
  dplyr::rename(`Days from specimen to outcome date` = t_to_outcome_group)

## hospital episodes with norovirus diagnosis code -------------

hosp_noro_episode_linkage_count <- molis_hosp_noro_episode |>
  dplyr::mutate(
    # Categorise time to hospital admission
    t_to_outcome_group = cut(
      specimen_to_arrival_days,
      breaks = c(-28, 1, 15, 29),
      labels = c("[-28, 0]", "[1,14]", "[15,28]"),
      right = FALSE
    ),
    t_to_outcome_group = factor(
      t_to_outcome_group,
      levels = c("[-28, 0]", "[1,14]", "[15,28]")
    )
  ) |>
  dplyr::filter(!is.na(t_to_outcome_group)) |>
  # keep all groups even if count is 0
  dplyr::group_by(t_to_outcome_group, .drop = FALSE) |>
  dplyr::summarise(
    `Count` = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    `Proportion of positive tests` = round(`Count` / n_tests, 3)
  ) |>
  dplyr::arrange(t_to_outcome_group) |>
  dplyr::rename(`Days from specimen to outcome date` = t_to_outcome_group)

## deaths -----------

deaths_linkage_count <- molis_deaths |>
  dplyr::mutate(
    # Categorise time to death
    t_to_outcome_group = cut(
      specimen_to_death_days,
      breaks = c(0, 14, 29, 61, 101),
      labels = c("[0,14]", "[15,28]", "[29,60]", "[61,100]"),
      right = FALSE
    )
  ) |>
  dplyr::filter(!is.na(t_to_outcome_group)) |>
  dplyr::summarise(
    `Count` = dplyr::n(),
    .by = t_to_outcome_group
  ) |>
  dplyr::mutate(
    `Proportion of positive tests` = round(`Count` / n_tests, 3)
  ) |>
  dplyr::arrange(t_to_outcome_group) |>
  dplyr::rename(`Days from specimen to outcome date` = t_to_outcome_group)

## combine subtables -------

# merge tables and format gt
outcome_linkage_count_gt <- rbind(
  overall_test_count,
  hosp_spell_linkage_count,
  hosp_noro_episode_linkage_count,
  deaths_linkage_count
) |>
  gt::gt() |>
  gt::tab_row_group(
    label = "Positive tests",
    rows = 1:1
  ) |>
  gt::tab_row_group(
    label = "Hospital spell admission",
    rows = 2:4
  ) |>
  gt::tab_row_group(
    label = "Hospital episode with norovirus diagnosis",
    rows = 5:7
  ) |>
  gt::tab_row_group(
    label = "Death",
    rows = 8:11
  ) |>
  gt::row_group_order(
    groups = c("Positive tests", "Hospital spell admission", "Hospital episode with norovirus diagnosis", "Death")
  ) |>
  gt::tab_options(column_labels.font.weight = "bold", row_group.font.weight = "bold")

outcome_linkage_count_gt |>
  gt::gtsave(filename = "outcome_linkage_count_table.png", path = output_path)

# MOLIS case fatality ratio summary table -----------

death_threshold <- 28

molis_deaths_filtered <- molis_deaths |>
  # set death dates to NA if further away than death_threshold days
  dplyr::mutate(
    death_date = ifelse(death_date - specimen_date <= death_threshold, death_date, NA)
  ) |>
  dplyr::mutate(
    # code fatality as 0 or 1
    fatality = ifelse(is.na(death_date), 0, 1)
  )

## Overall ---------

# Calculate case counts and deaths overall
molis_overall_cfr <- molis_deaths_filtered |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3)
  ) |>
  dplyr::mutate(Covariate = "") |>
  dplyr::select(Covariate, `Positive tests`, `Deaths`, `Deaths / positive tests`)

# Calculate case counts and deaths by genotype
molis_genotype_cfr <- molis_deaths_filtered |>
  dplyr::mutate(
    genotype = factor(genotype, levels = c("GII.4", "GII.17", "Other", "Unknown"))
  ) |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = genotype
  ) |>
  dplyr::arrange(genotype) |>
  dplyr::rename(Covariate = genotype)

# Calculate case counts and deaths by region
molis_region_cfr <- molis_deaths_filtered |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = region
  ) |>
  dplyr::arrange(region) |>
  dplyr::rename(Covariate = region)

# Calculate case counts and deaths by age
molis_age_cfr <- molis_deaths_filtered |>
  dplyr::filter(!is.na(age)) |>
  dplyr::mutate(
    # Categorise age
    age_group = cut(
      age,
      breaks = c(0, 5, 20, 50, 75, Inf),
      labels = c("0-4", "5-19", "20-49", "50-74", "75+"),
      right = FALSE
    )
  ) |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = age_group
  ) |>
  dplyr::arrange(age_group) |>
  dplyr::rename(Covariate = age_group)

# merge tables and format gt
molis_cfr_gt <- rbind(molis_overall_cfr, molis_genotype_cfr, molis_age_cfr, molis_region_cfr) |>
  gt::gt() |>
  gt::tab_row_group(
    label = "Overall",
    rows = 1:1
  ) |>
  gt::tab_row_group(
    label = "Genotype",
    rows = 2:5
  ) |>
  gt::tab_row_group(
    label = "Age group",
    rows = 6:10
  ) |>
  gt::tab_row_group(
    label = "Region",
    rows = 11:17
  ) |>
  gt::row_group_order(groups = c("Overall", "Genotype", "Age group", "Region")) |>
  gt::tab_options(column_labels.font.weight = "bold", row_group.font.weight = "bold")

molis_cfr_gt |>
  gt::gtsave(filename = "molis_raw_cfr_table.png", path = output_path)


# SGSS case fatality ratio summary table -----------

death_threshold <- 28

sgss_deaths_filtered <- sgss_deaths |>
  # set death dates to NA if further away than death_threshold days
  dplyr::mutate(
    death_date = ifelse(death_date - specimen_date <= death_threshold, death_date, NA)
  ) |>
  dplyr::mutate(
    # code fatality as 0 or 1
    fatality = ifelse(is.na(death_date), 0, 1)
  )

## Overall ---------

# Calculate case counts and deaths overall
sgss_overall_cfr <- sgss_deaths_filtered |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3)
  ) |>
  dplyr::mutate(Covariate = "") |>
  dplyr::select(Covariate, `Positive tests`, `Deaths`, `Deaths / positive tests`)

# Calculate case counts and deaths by requesting_organisation
sgss_requesting_organisation_cfr <- sgss_deaths_filtered |>
  dplyr::mutate(
    requesting_organisation = factor(
      requesting_organisation,
      levels = c("Primary care", "Secondary care", "Other", "Unknown")
    )
  ) |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = requesting_organisation
  ) |>
  dplyr::arrange(requesting_organisation) |>
  dplyr::rename(Covariate = requesting_organisation)

# Calculate case counts and deaths by age
sgss_age_cfr <- sgss_deaths_filtered |>
  dplyr::filter(!is.na(age)) |>
  dplyr::mutate(
    # Categorise age
    age_group = cut(
      age,
      breaks = c(0, 5, 20, 50, 75, Inf),
      labels = c("0-4", "5-19", "20-49", "50-74", "75+"),
      right = FALSE
    )
  ) |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = age_group
  ) |>
  dplyr::arrange(age_group) |>
  dplyr::rename(Covariate = age_group)

# Calculate case counts and deaths by region
sgss_region_cfr <- sgss_deaths_filtered |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = sgss_lab_region
  ) |>
  dplyr::arrange(sgss_lab_region) |>
  dplyr::rename(Covariate = sgss_lab_region)

# merge tables and format gt
sgss_cfr_gt <- rbind(
  sgss_overall_cfr,
  sgss_age_cfr,
  sgss_requesting_organisation_cfr,
  sgss_region_cfr
) |>
  gt::gt() |>
  gt::tab_row_group(
    label = "Overall",
    rows = 1:1
  ) |>
  gt::tab_row_group(
    label = "Age group",
    rows = 2:6
  ) |>
  gt::tab_row_group(
    label = "Healthcare level",
    rows = 7:10
  ) |>
  gt::tab_row_group(
    label = "Lab region",
    rows = 11:19
  ) |>
  gt::row_group_order(
    groups = c("Overall", "Age group", "Healthcare level")
  ) |>
  gt::tab_options(column_labels.font.weight = "bold", row_group.font.weight = "bold")

sgss_cfr_gt |>
  gt::gtsave(filename = "sgss_raw_cfr_table.png", path = output_path)


# MOLIS-SGSS case fatality ratio summary table -----------

death_threshold <- 28

molis_sgss_deaths_filtered <- molis_sgss_deaths |>
  # set death dates to NA if further away than death_threshold days
  dplyr::mutate(
    death_date = ifelse(death_date - specimen_date <= death_threshold, death_date, NA)
  ) |>
  dplyr::mutate(
    # code fatality as 0 or 1
    fatality = ifelse(is.na(death_date), 0, 1)
  )

## Overall ---------

# Calculate case counts and deaths overall
molis_sgss_overall_cfr <- molis_sgss_deaths_filtered |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3)
  ) |>
  dplyr::mutate(Covariate = "") |>
  dplyr::select(Covariate, `Positive tests`, `Deaths`, `Deaths / positive tests`)

# Calculate case counts and deaths by requesting_organisation
molis_sgss_requesting_organisation_cfr <- molis_sgss_deaths_filtered |>
  dplyr::mutate(
    requesting_organisation = factor(
      requesting_organisation,
      levels = c("Primary care", "Secondary care", "Other", "Unknown")
    )
  ) |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = requesting_organisation
  ) |>
  dplyr::arrange(requesting_organisation) |>
  dplyr::rename(Covariate = requesting_organisation)

# Calculate case counts and deaths by age
molis_sgss_age_cfr <- molis_sgss_deaths_filtered |>
  dplyr::filter(!is.na(age)) |>
  dplyr::mutate(
    # Categorise age
    age_group = cut(
      age,
      breaks = c(0, 5, 20, 50, 75, Inf),
      labels = c("0-4", "5-19", "20-49", "50-74", "75+"),
      right = FALSE
    )
  ) |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = age_group
  ) |>
  dplyr::arrange(age_group) |>
  dplyr::rename(Covariate = age_group)

# Calculate case counts and deaths by genotype
molis_sgss_genotype_cfr <- molis_sgss_deaths_filtered |>
  dplyr::mutate(
    genotype = factor(genotype, levels = c("GII.4", "GII.17", "Other", "Unknown"))
  ) |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = genotype
  ) |>
  dplyr::arrange(genotype) |>
  dplyr::rename(Covariate = genotype)

# Calculate case counts and deaths by broad region
molis_sgss_broad_region_cfr <- molis_sgss_deaths_filtered |>
  dplyr::mutate(
    sgss_lab_broad_region = factor(sgss_lab_broad_region, levels = c("North", "Midlands", "South"))
  ) |>
  dplyr::summarise(
    `Positive tests` = dplyr::n(),
    `Deaths` = sum(fatality, na.rm = TRUE),
    `Deaths / positive tests` = round(`Deaths` / `Positive tests`, 3),
    .by = sgss_lab_broad_region
  ) |>
  dplyr::arrange(sgss_lab_broad_region) |>
  dplyr::rename(Covariate = sgss_lab_broad_region)


# merge tables and format gt
molis_sgss_cfr_gt <- rbind(
  molis_sgss_overall_cfr,
  molis_sgss_genotype_cfr,
  molis_sgss_age_cfr,
  molis_sgss_requesting_organisation_cfr,
  molis_sgss_broad_region_cfr
) |>
  gt::gt() |>
  gt::tab_row_group(
    label = "Overall",
    rows = 1:1
  ) |>
  gt::tab_row_group(
    label = "Genotype",
    rows = 2:5
  ) |>
  gt::tab_row_group(
    label = "Age group",
    rows = 6:10
  ) |>
  gt::tab_row_group(
    label = "Healthcare level",
    rows = 11:14
  ) |>
  gt::tab_row_group(
    label = "Lab broad region",
    rows = 15:17
  ) |>
  gt::row_group_order(
    groups = c("Overall", "Age group", "Genotype", "Healthcare level", "Lab broad region")
  ) |>
  gt::tab_options(column_labels.font.weight = "bold", row_group.font.weight = "bold")

molis_sgss_cfr_gt |>
  gt::gtsave(filename = "molis_sgss_raw_cfr_table.png", path = output_path)
