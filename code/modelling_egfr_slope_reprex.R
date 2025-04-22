#*******************************************************************************
#
# Project:      Modelling eGFR slope
# Last updated: 19-Jan-2025
# Authors:      Robert Fletcher and Niels Jongs
# Purpose:      Model eGFR slope using synthetic data
# Contact:      rfletcher@georgeinstitute.org.au | n.jongs@umcg.nl
#
#*******************************************************************************


# Notes -------------------------------------------------------------------

# This reproducible example illustrates how to calculate estimated glomerular 
# filtration rate (eGFR) slope. It demonstrates the required data format and the
# steps needed to replicate the analyses in your own clinical trial data

# For additional details or support, please consult the repository’s README.md 
# or use the contact information provided above

# Please note, this example assumes your dataset includes repeated eGFR
# measurements over time, in addition to any relevant covariates for the 
# analyses


# Install dependencies (if not currently installed) -----------------------

# Library names
libs <- c("glue", "multcomp", "lme4", "tidyverse")

# Install libraries
install.packages(setdiff(libs, rownames(installed.packages())))


# Load libraries ----------------------------------------------------------

library(glue)
library(nlme)
library(multcomp)
library(tidyverse)


# Define variables --------------------------------------------------------

# Path to cloned repository
# This assumes you're using RStudio. It extracts the current script's directory
# and removes the "/code" suffix. If you're not using RStudio, specify the path 
# manually
path <- sub("/code$", "", dirname(rstudioapi::getSourceEditorContext()$path))

# Set spline knot point
# In the CREDENCE trial, 21 days (3 weeks) corresponds to the first study visit,
# which marks the initial acute drop in eGFR with SGLT2 inhibition. Adjust this
# value according to the available date in your trial. For example, in DAPA-CKD, 
# the first post-randomisation visit is at 14 days (2 weeks), so `k` would be 
# set to 14
k <- 14


# Source functions --------------------------------------------------------

source(glue::glue("{path}/src/generate_synthetic_data.R"))
source(glue::glue("{path}/src/compute_slope.R"))


# Load data ---------------------------------------------------------------

synth_bl <- generate_synthetic_data(.table = "baseline")
synth_fu <- generate_synthetic_data(.table = "follow-up")

# If you encounter issues running the function above, the two datasets generated 
# by this function are also available in CSV format in the `data` sub-directory 
# of this repository. Use the following functions to load those files:

# bl <- readr::read_csv("{path}/data/synthetic_trial_baseline.csv")
# fu <- readr::read_csv("{path}/data/synthetic_trial_egfr_follow_up.csv")

bl <- haven::read_dta("/Volumes/JHund/DAPA-HF/SMART-C/slope_bl.dta")
bl$hfhospn = factor(bl$hfhospn, labels = c("Y", "N"))

fu <- haven::read_dta("/Volumes/JHund/DAPA-HF/SMART-C/slope_fu.dta")
fu

# Prepare data ------------------------------------------------------------

# Treatment Arms
# Here, we filter for on-treatment individuals only. NOTE: This step is 
# optional. Depending on your study design, you might prefer an 
# intention-to-treat population instead
arm <- bl |> 
  dplyr::select(usubjid, trt01pn = trtpn, egfrgr1, t2dblfl, stratan, hfhospn, acearbarni, t2dblfl) |> 
  dplyr::mutate(
    egfrgr1 = factor(
      egfrgr1,
      levels = c("< 30", "30 and < 45", "45 and < 60", ">= 60")
      )
  )

# Repeat eGFR
# Here, we filter for on-treatment measurements and those flagged for analysis. 
# Adjust these criteria according to your chosen analysis strategy
gfr <- fu |> 
  dplyr::distinct(usubjid, aval,ady, base, avisitn)

# Join
gfr_c <- arm |> 
  dplyr::left_join(gfr, by = dplyr::join_by(usubjid), multiple = "all") |> 
  dplyr::mutate(
    # Convert days from baseline to years
    time = ady / 365.25,
    # Generate spline term
    spline = dplyr::if_else(time >= k / 365.25, time - k / 365.25, 0)
  )


gfr_c$stratan <- factor(gfr_c$stratan)

# Formatting the data for the model is likely the most crucial step. The 
# synthetic data used in this reprex should give you an idea of the required 
# format, but consult the repository README if you're unsure about anything or 
# would like a more detailed explanation


# Fit model (whole population) --------------------------------------------

# Fit mixed effects model with unstructured residual variance-covariance matrix
fit <- lme4::lmer(
  aval ~ base + stratan + hfhospn + time * trt01pn + spline * trt01pn - 1 + 
  (time | usubjid), data = gfr_c
)

# The above model is for the overall trial population. Subgroup-specific 
# analyses are detailed further down

# Please see the repository README if you'd like more information on 
# specification of the model parameters


# Extract model results (whole population) --------------------------------

# Define the proportion of the total slope attributed to the chronic slope. In 
# this example, 1095.75 represents the total trial follow-up period (roughly 3 
# years in days). By subtracting 21 days (acute slope), we calculate the
# fraction that remains for the chronic slope
maxfup <- 835 # max(t2dth) in DAPA
prop <- (835 - k) / 835

# Compute eGFR slope for the whole population
all <- compute_slope(
  .model_obj = fit, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline", 
  .prop = prop, 
  .output = "all"
)

all

# Fit model and compute slope by binary subgroup variable -----------------

# Recode baseline GLP-1RA use as factor and set `blglp1` == "Yes" as reference
gfr_c <- gfr_c  |> 
  dplyr::mutate(
    acearbarni = factor(acearbarni, levels = c(1, 0), labels = c("Yes", "No")),
    t2dblfl = factor(t2dblfl, levels = c("Y", "N"), labels = c("Yes", "No"))
  )

# Fit mixed effects model with unstructured residual variance-covariance matrix
# this time with adjustment and interactions with `acearbarni`
fit_binary_rasblock <- lme4::lmer(
  aval ~ base + stratan + hfhospn + time * trt01pn + spline * trt01pn + time * acearbarni +
  spline * acearbarni + trt01pn * acearbarni + time * trt01pn * acearbarni + 
  spline * trt01pn * acearbarni - 1 + (time | usubjid), data = gfr_c
)

fit_binary_diabetes <- lme4::lmer(
  aval ~ base + stratan + hfhospn + time * trt01pn + spline * trt01pn + time * t2dblfl +
    spline * t2dblfl + trt01pn * t2dblfl + time * trt01pn * t2dblfl + 
    spline * trt01pn * t2dblfl - 1 + (time | usubjid), data = gfr_c
)

# Compute eGFR slope by binary subgroups (baseline GLP-1RA use)
total_ras <- compute_slope(
  .model_obj = fit_binary_rasblock, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline", 
  .prop = prop, 
  .by = "acearbarni", # `.by` argument specified with subgroup variable
  .output = "total"
) |> print()

chronic_ras <- compute_slope(
  .model_obj = fit_binary_rasblock, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline", 
  .prop = prop, 
  .by = "acearbarni", # `.by` argument specified with subgroup variable
  .output = "chronic"
) |> print()

total_t2db <- compute_slope(
  .model_obj = fit_binary_diabetes, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline", 
  .prop = prop, 
  .by = "t2dblfl", # `.by` argument specified with subgroup variable
  .output = "total"
) |> print()

total_t2db <- compute_slope(
  .model_obj = fit_binary_diabetes, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline", 
  .prop = prop, 
  .by = "t2dblfl", # `.by` argument specified with subgroup variable
  .output = "chronic"
) |> print()


# Fit model and compute slope by multi-level subgroup variable ------------

# Fit mixed effects model with unstructured residual variance-covariance matrix
# with adjustment and interactions with `gfr_grp`
fit_ordinal <- lme4::lmer(
  aval ~ base + stratan + hfhospn + time * trt01pn + spline * trt01pn + time * egfrgr1 +
  spline * egfrgr1 + trt01pn * egfrgr1 + time * trt01pn * egfrgr1 + 
  spline * trt01pn * egfrgr1 - 1 + (time | usubjid), data = gfr_c
)

# Compute eGFR slope by multiple subgroups (baseline eGFR quartiles)
ordinal_subgroups <- compute_slope(
  .model_obj = fit_ordinal, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline",
  .prop = prop, 
  .by = "egfrgr1", # `.by` argument specified with subgroup variable
  .output = "total"
) |> print()

chronic_gfrgroup <- compute_slope(
  .model_obj = fit_ordinal, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline",
  .prop = prop, 
  .by = "egfrgr1", # `.by` argument specified with subgroup variable
  .output = "chronic"
) |> print()




# egfr by diabetes --------------------------------------------------------
gfr_c_nodiab <- gfr_c |> 
  filter(t2dblfl == "No", egfrgr1 != "< 30")

gfr_c_diab <- gfr_c |> 
  filter(t2dblfl == "Yes", egfrgr1 != "< 30")

db_egfr <- lme4::lmer(
  aval ~ base + stratan + hfhospn + time * trt01pn + spline * trt01pn + time * egfrgr1 +
    spline * egfrgr1 + trt01pn * egfrgr1 + time * trt01pn * egfrgr1 + 
    spline * trt01pn * egfrgr1 - 1 + (time | usubjid), data = gfr_c_diab
)

no_db_egfr <- lme4::lmer(
  aval ~ base + stratan + hfhospn + time * trt01pn + spline * trt01pn + time * egfrgr1 +
    spline * egfrgr1 + trt01pn * egfrgr1 + time * trt01pn * egfrgr1 + 
    spline * trt01pn * egfrgr1 - 1 + (time | usubjid), data = gfr_c_nodiab
)
  
# Compute eGFR slope by multiple subgroups (baseline eGFR quartiles)
compute_slope(
  .model_obj = db_egfr, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline",
  .prop = prop, 
  .by = "egfrgr1", # `.by` argument specified with subgroup variable
  .output = "total"
) |> print()

compute_slope(
  .model_obj = db_egfr, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline",
  .prop = prop, 
  .by = "egfrgr1", # `.by` argument specified with subgroup variable
  .output = "chronic"
) |> print()  


# Compute eGFR slope by multiple subgroups (baseline eGFR quartiles)
compute_slope(
  .model_obj = no_db_egfr, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline",
  .prop = prop, 
  .by = "egfrgr1", # `.by` argument specified with subgroup variable
  .output = "total"
) |> print()

compute_slope(
  .model_obj = no_db_egfr, 
  .time_var = "time", 
  .intervention_var = "trt01pn", 
  .spline_var = "spline",
  .prop = prop, 
  .by = "egfrgr1", # `.by` argument specified with subgroup variable
  .output = "chronic"
) |> print()  
