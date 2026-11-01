options(survey.lonely.psu = "adjust")

required_packages <- c("broom", "dplyr", "ggplot2", "openxlsx", "readr", "readxl", "rms", "survey")
invisible(lapply(required_packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) stop("Missing package: ", x)
}))

source("config.R")
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readxl::read_excel(config$nhanes_analysis_file, sheet = config$nhanes_sheet)

required_columns <- c(
  "SEQN", "SDMVSTRA", "SDMVPSU", "WTDR2D", "Cyanidin_mean", "MASLD",
  "Gender", "Age", "Race", "Marital", "EDUcation", "PIR", "Smoke", "ALQ",
  "DRTKCAL", "DRTPROT", "BMI", "CAP", "LSM", "HSI", "FLI"
)
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns)) stop("Missing columns: ", paste(missing_columns, collapse = ", "))

# Complete-case analysis used for the reported population model.
analysis_vars <- c(
  "SEQN", "SDMVSTRA", "SDMVPSU", "WTDR2D", "Cyanidin_mean", "MASLD",
  "Gender", "Age", "Race", "Marital", "EDUcation", "PIR", "Smoke", "ALQ",
  "DRTKCAL", "DRTPROT", "BMI"
)
analysis_dat <- dat |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(analysis_vars), ~ !is.na(.x))) |>
  dplyr::mutate(
    cyanidin_log1p = log1p(Cyanidin_mean),
    cyanidin_quartile = factor(
      dplyr::ntile(Cyanidin_mean, 4), levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4")
    ),
    MASLD = as.integer(MASLD),
    Gender = factor(Gender),
    Race = factor(Race),
    Marital = factor(Marital),
    Education = factor(EDUcation),
    Smoke = factor(Smoke),
    Alcohol = factor(ALQ),
    BMI_category = cut(
      BMI, breaks = c(-Inf, 25, 30, Inf), right = FALSE,
      labels = c("<25", "25 to <30", ">=30")
    )
  )

if (nrow(analysis_dat) != 3882L) {
  warning("Analysis row count is ", nrow(analysis_dat), "; manuscript reports 3,882. Reconcile exclusions.")
}

# NHANES 2017-2018 is a single two-year cycle; use the supplied two-day dietary weight.
# Do not rescale WTDR2D for population estimates.
nhanes_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDR2D,
  nest = TRUE,
  data = analysis_dat
)

primary_formula <- MASLD ~ cyanidin_quartile + Gender + Age + Race + Marital +
  Education + PIR + Smoke + Alcohol + DRTKCAL + DRTPROT
bmi_formula <- update(primary_formula, . ~ . + BMI_category)

primary_model <- survey::svyglm(primary_formula, design = nhanes_design, family = quasibinomial())
bmi_model <- survey::svyglm(bmi_formula, design = nhanes_design, family = quasibinomial())

tidy_or <- function(model, model_name) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::mutate(model = model_name) |>
    dplyr::select(model, term, estimate, conf.low, conf.high, p.value)
}

model_results <- dplyr::bind_rows(
  tidy_or(primary_model, "Primary model (without BMI)"),
  tidy_or(bmi_model, "Sensitivity model (+ BMI category)")
)
readr::write_csv(model_results, file.path(config$results_dir, "logistic_models.csv"))

# Restricted cubic spline. The four knots are survey-weighted quantiles of log1p intake.
knot_probs <- c(0.05, 0.35, 0.65, 0.95)
knots <- as.numeric(survey::svyquantile(~cyanidin_log1p, nhanes_design, knot_probs, ci = FALSE))
spline_primary <- survey::svyglm(
  MASLD ~ rms::rcs(cyanidin_log1p, knots) + Gender + Age + Race + Marital +
    Education + PIR + Smoke + Alcohol + DRTKCAL + DRTPROT,
  design = nhanes_design,
  family = quasibinomial()
)
spline_bmi <- update(spline_primary, . ~ . + BMI_category)

saveRDS(
  list(primary = primary_model, bmi = bmi_model, spline_primary = spline_primary,
       spline_bmi = spline_bmi, knots = knots),
  file.path(config$results_dir, "fitted_models.rds")
)

sample_counts <- data.frame(
  analysis = c("Main analysis", "CAP", "LSM", "HSI", "FLI"),
  n = c(nrow(dat), sum(!is.na(dat$CAP)), sum(!is.na(dat$LSM)),
        sum(!is.na(dat$HSI)), sum(!is.na(dat$FLI)))
)
readr::write_csv(sample_counts, file.path(config$results_dir, "analysis_sample_counts.csv"))

# Basic reproducibility record.
capture.output(sessionInfo(), file = file.path(config$results_dir, "sessionInfo.txt"))
