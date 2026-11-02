
# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tidyverse)


# ------------------------------------------------------------
# 2. Set working directory
# ------------------------------------------------------------

# setwd!!!!!
# setwd()


# ------------------------------------------------------------
# 3. Read ESS data
# ------------------------------------------------------------

ess <- read_csv("Data/ESS11e04_2-subset.csv")

# Check for parsing problems
problems(ess)


# ------------------------------------------------------------
# 4. Select and clean variables
# ------------------------------------------------------------

ess_data <- ess %>%
  select(
    idno,
    cntry,
    gndr,
    agea,
    edctn,
    wkhtot,
    wkhct,
    pdwrk,
    pspwght
  ) %>%
  mutate(
    # Country factor
    country = factor(cntry),

    # Gender (Ref: Male)
    gen = case_when(
      gndr == 1 ~ "Male",
      gndr == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    gen = factor(gen, levels = c("Male", "Female")),

    # Age (18–70 years)
    age = if_else(
      agea >= 18 & agea <= 70,
      as.numeric(agea),
      NA_real_
    ),

    # Current participation in education
    edu = case_when(
      edctn == 1 ~ "in education",
      edctn == 0 ~ "not in education",
      TRUE ~ NA_character_
    ),
    edu = factor(edu, levels = c("not in education", "in education")),

    # DV 1 for H1a/H2a: Employment Status (1 = Paid work, 0 = No paid work)
    emp = case_when(
      pdwrk == 1 ~ 1L,
      pdwrk == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),

    # Flag: total hours replaced by contracted hours (for sensitivity analysis)
    hours_subst = coalesce(
      emp == 1L & wkhtot == 0 & wkhct > 0 & wkhct <= 168,
      FALSE
    ),

    # Weekly working hours for employed respondents
    weekly_hours = case_when(
      wkhtot > 0 & wkhtot <= 168 ~ as.numeric(wkhtot),
      wkhtot == 0 & wkhct > 0 & wkhct <= 168 ~ as.numeric(wkhct),
      TRUE ~ NA_real_
    ),

    # DV 2 for H1b/H2b: Working hours ONLY for employed (NA for non-employed)
    wwh = if_else(emp == 1L, weekly_hours, NA_real_)
  )


# ------------------------------------------------------------
# 5. Define EU-27 Member States
# ------------------------------------------------------------

eu_countries <- c(
  "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "ES", "FI", 
  "FR", "GR", "HR", "HU", "IE", "IT", "LT", "LU", "LV", "MT", 
  "NL", "PL", "PT", "RO", "SE", "SI", "SK"
)


# ------------------------------------------------------------
# 6. Construct final analysis samples
# ------------------------------------------------------------

# Sample A: all respondents with known employment status (part a)
ess_full <- ess_data %>%
  filter(
    cntry %in% eu_countries,
    !is.na(emp), 
    !is.na(gen), 
    !is.na(age), 
    !is.na(edu), 
    !is.na(pspwght)
  ) %>%
  mutate(
    country = droplevels(country),
    age_c   = age - mean(age),
    age_c2  = age_c^2
  )

# Sample B: employed respondents with valid hours (part b)
ess_hours <- ess_full %>%
  filter(emp == 1L, !is.na(wwh)) %>%
  mutate(
    age_c  = age - mean(age),
    age_c2 = age_c^2
  )



# ------------------------------------------------------------
# 8. Export Data
# ------------------------------------------------------------

saveRDS(ess_full, "Data/ess_full.rds")
saveRDS(ess_hours, "Data/ess_hours.rds")

