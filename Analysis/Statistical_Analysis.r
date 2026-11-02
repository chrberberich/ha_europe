# ------------------------------------------------------------
# 1. Packages and settings
# ------------------------------------------------------------

library(tidyverse)
library(lme4)
library(marginaleffects)
library(DHARMa)

# setwd!!!!!
# setwd()

# Optimizer settings (random slopes with ~22 countries can be fragile)
ctrl_glmer <- glmerControl(optimizer = "bobyqa",
                           optCtrl = list(maxfun = 2e5))
ctrl_lmer  <- lmerControl(optimizer = "bobyqa",
                          optCtrl = list(maxfun = 2e5))


# ------------------------------------------------------------
# 2. Read prepared analysis data
# ------------------------------------------------------------

# ess_full : all respondents with known employment status (part a)
# ess_hours: employed respondents with valid hours (part b)
ess_full  <- readRDS("Data/ess_full.rds")
ess_hours <- readRDS("Data/ess_hours.rds")

# Quick sample overview
nrow(ess_full)
nrow(ess_hours)
n_distinct(ess_full$country)
n_distinct(ess_hours$country)


# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------


lrt_mixture <- function(m0, m1) {
  lr <- max(as.numeric(2 * (logLik(m1) - logLik(m0))), 0)
  tibble(
    LR        = lr,
    p_mixture = 0.5 * pchisq(lr, df = 1, lower.tail = FALSE) +
                0.5 * pchisq(lr, df = 2, lower.tail = FALSE)
  )
}

# SD of the random slope for gender
slope_sd <- function(m) {
  vc <- as.data.frame(VarCorr(m))
  vc$sdcor[vc$grp == "country" & vc$var1 == "genFemale" & is.na(vc$var2)]
}

# Convergence / singularity check
check_fit <- function(m, label) {
  cat("\n---", label, "---\n")
  cat("Singular fit:", isSingular(m), "\n")
  msg <- m@optinfo$conv$lme4$messages
  if (is.null(msg)) cat("No optimizer messages\n") else print(msg)
}

# Compact summary of a random-slope model (for robustness tables)
summarise_rs <- function(m, label) {
  tibble(
    model    = label,
    n        = nobs(m),
    beta_gen = fixef(m)[["genFemale"]],
    se_gen   = sqrt(diag(as.matrix(vcov(m))))[["genFemale"]],
    sd_slope = slope_sd(m)
  )
}


# ============================================================
# PART (a): EMPLOYMENT STATUS  --  P(employed)
# ============================================================

# Model 1: random intercept
hurdle_binary_ri <- glmer(
  emp ~ gen + age_c + I(age_c^2) + edu + (1 | country),
  family  = binomial(link = "logit"),
  data    = ess_full,
  control = ctrl_glmer
)

# Model 2: random intercept + random slope for gender
hurdle_binary_rs <- glmer(
  emp ~ gen + age_c + I(age_c^2) + edu + (1 + gen | country),
  family  = binomial(link = "logit"),
  data    = ess_full,  
  control = ctrl_glmer
)

summary(hurdle_binary_ri)
summary(hurdle_binary_rs)

check_fit(hurdle_binary_ri, "Part (a), model 1")
check_fit(hurdle_binary_rs, "Part (a), model 2")

# LR test (boundary-corrected), AIC, BIC
lrt_binary <- lrt_mixture(hurdle_binary_ri, hurdle_binary_rs)
lrt_binary
AIC(hurdle_binary_ri, hurdle_binary_rs)
BIC(hurdle_binary_ri, hurdle_binary_rs)

final_binary <- hurdle_binary_rs


# ============================================================
# PART (b): WEEKLY HOURS OF THE EMPLOYED  --  E(hours | employed)
# ============================================================

# ML (REML = FALSE) so that the LR test is valid
hurdle_hours_ri <- lmer(
  wwh ~ gen + age_c + I(age_c^2) + edu + (1 | country),
  data    = ess_hours,
  REML    = FALSE,
  control = ctrl_lmer
)

hurdle_hours_rs <- lmer(
  wwh ~ gen + age_c + I(age_c^2) + edu + (1 + gen | country),
  data    = ess_hours,
  REML    = FALSE,
  control = ctrl_lmer
)

summary(hurdle_hours_ri)
summary(hurdle_hours_rs)

check_fit(hurdle_hours_ri, "Part (b), model 1")
check_fit(hurdle_hours_rs, "Part (b), model 2")

lrt_hours <- lrt_mixture(hurdle_hours_ri, hurdle_hours_rs)
lrt_hours
AIC(hurdle_hours_ri, hurdle_hours_rs)
BIC(hurdle_hours_ri, hurdle_hours_rs)

final_hours <- hurdle_hours_rs

# Check: N used must equal the sample sizes
stopifnot(nobs(final_binary) == nrow(ess_full))
stopifnot(nobs(final_hours)  == nrow(ess_hours))




# ============================================================
# RANDOM EFFECTS AND VARIANCE OF THE GENDER GAP
# ============================================================

VarCorr(final_binary)
VarCorr(final_hours)

# Fixed effects with Wald confidence intervals
fixef(final_binary)
fixef(final_hours)
confint(final_binary, parm = "beta_", method = "Wald")
confint(final_hours,  parm = "beta_", method = "Wald")


ci_sd_hours <- try(confint(final_hours, parm = "theta_",
                           method = "profile"), silent = TRUE)
ci_sd_hours



# Range containing about 95 % of the country-specific gaps (part b)
b_h <- fixef(final_hours)[["genFemale"]]
s_h <- slope_sd(final_hours)
b_h + c(-1.96, 1.96) * s_h

# Conditional modes of the random effects (with conditional SDs)
ranef_binary <- as.data.frame(ranef(final_binary, condVar = TRUE))
ranef_hours  <- as.data.frame(ranef(final_hours,  condVar = TRUE))


# ============================================================
# EXPLORATORY: COUNTRY-SPECIFIC GENDER GAPS
# ============================================================

gap_employed <- avg_comparisons(
  final_binary, variables = "gen", by = "country", type = "response"
)
gap_hours <- avg_comparisons(
  final_hours, variables = "gen", by = "country"
)

plot_gap <- function(gap, xlab) {
  gap |>
    mutate(country = fct_reorder(country, estimate)) |>
    ggplot(aes(x = estimate, y = country)) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(x = xlab, y = NULL) +
    theme_minimal()
}

plot_gap(gap_employed, "Gender gap in employment probability (women - men)")
plot_gap(gap_hours,    "Gender gap in weekly hours (women - men)")

gap_rank <- gap_employed |>
  select(country, gap_emp = estimate) |>
  left_join(select(gap_hours, country, gap_hrs = estimate),
            by = "country")

cor(gap_rank$gap_emp, gap_rank$gap_hrs, method = "spearman")


# ============================================================
# MODEL DIAGNOSTICS (same formulas as the reported models)
# ============================================================

sim_binary <- simulateResiduals(final_binary, n = 250)
plot(sim_binary)   # dispersion tests are not informative for binary data

sim_hours <- simulateResiduals(final_hours, n = 250)
plot(sim_hours)
testDispersion(sim_hours)
# Expect deviations: hours cluster at 20 and 40 (not Gaussian).






# ============================================================
# AMEs berechnen und exportieren
# ============================================================


# ------------------------------------------------------------
# Logit RI – Population Level
# ------------------------------------------------------------

ame_logit_ri <- avg_comparisons(
  hurdle_binary_ri,
  variables = c("gen", "edu"),
  type = "response",
  re.form = NA
)

ame_logit_ri_age <- avg_slopes(
  hurdle_binary_ri,
  variables = "age_c",
  type = "response",
  re.form = NA
)


# ------------------------------------------------------------
# Logit RI – inklusive Random Effects
# ------------------------------------------------------------

ame_logit_ri_re <- avg_comparisons(
  hurdle_binary_ri,
  variables = c("gen", "edu"),
  type = "response",
  re.form = NULL
)

ame_logit_ri_age_re <- avg_slopes(
  hurdle_binary_ri,
  variables = "age_c",
  type = "response",
  re.form = NULL
)


# ------------------------------------------------------------
# Logit RS – Population Level
# ------------------------------------------------------------

ame_logit_rs <- avg_comparisons(
  hurdle_binary_rs,
  variables = c("gen", "edu"),
  type = "response",
  re.form = NA
)

ame_logit_rs_age <- avg_slopes(
  hurdle_binary_rs,
  variables = "age_c",
  type = "response",
  re.form = NA
)


# ------------------------------------------------------------
# Logit RS – inklusive Random Effects
# ------------------------------------------------------------

ame_logit_rs_re <- avg_comparisons(
  hurdle_binary_rs,
  variables = c("gen", "edu"),
  type = "response",
  re.form = NULL
)

ame_logit_rs_age_re <- avg_slopes(
  hurdle_binary_rs,
  variables = "age_c",
  type = "response",
  re.form = NULL
)


# ------------------------------------------------------------
# Gauss RI – Population Level
# ------------------------------------------------------------

ame_gauss_ri <- avg_comparisons(
  hurdle_hours_ri,
  variables = c("gen", "edu"),
  re.form = NA
)

ame_gauss_ri_age <- avg_slopes(
  hurdle_hours_ri,
  variables = "age_c",
  re.form = NA
)


# ------------------------------------------------------------
# Gauss RI – inklusive Random Effects
# ------------------------------------------------------------

ame_gauss_ri_re <- avg_comparisons(
  hurdle_hours_ri,
  variables = c("gen", "edu"),
  re.form = NULL
)

ame_gauss_ri_age_re <- avg_slopes(
  hurdle_hours_ri,
  variables = "age_c",
  re.form = NULL
)


# ------------------------------------------------------------
# Gauss RS – Population Level
# ------------------------------------------------------------

ame_gauss_rs <- avg_comparisons(
  hurdle_hours_rs,
  variables = c("gen", "edu"),
  re.form = NA
)

ame_gauss_rs_age <- avg_slopes(
  hurdle_hours_rs,
  variables = "age_c",
  re.form = NA
)


# ------------------------------------------------------------
# Gauss RS – inklusive Random Effects
# ------------------------------------------------------------

ame_gauss_rs_re <- avg_comparisons(
  hurdle_hours_rs,
  variables = c("gen", "edu"),
  re.form = NULL
)

ame_gauss_rs_age_re <- avg_slopes(
  hurdle_hours_rs,
  variables = "age_c",
  re.form = NULL
)


# ============================================================
# Get Together
# ============================================================

ames <- list(

  # Population-Level AMEs
  logit_ri = bind_rows(
    ame_logit_ri,
    ame_logit_ri_age
  ),

  logit_rs = bind_rows(
    ame_logit_rs,
    ame_logit_rs_age
  ),

  gauss_ri = bind_rows(
    ame_gauss_ri,
    ame_gauss_ri_age
  ),

  gauss_rs = bind_rows(
    ame_gauss_rs,
    ame_gauss_rs_age
  ),


  # AMEs inklusive Random Effects
  logit_ri_re = bind_rows(
    ame_logit_ri_re,
    ame_logit_ri_age_re
  ),

  logit_rs_re = bind_rows(
    ame_logit_rs_re,
    ame_logit_rs_age_re
  ),

  gauss_ri_re = bind_rows(
    ame_gauss_ri_re,
    ame_gauss_ri_age_re
  ),

  gauss_rs_re = bind_rows(
    ame_gauss_rs_re,
    ame_gauss_rs_age_re
  )

)




# ============================================================
# EXPORT FOR QUARTO
# ============================================================

dir.create("Output", showWarnings = FALSE)

# AMEs

saveRDS(ames, "Output/ames.rds")

# Model comparison
saveRDS(lrt_binary, "Output/lrt_binary.rds")
saveRDS(lrt_hours,  "Output/lrt_hours.rds")

# Final models and all four models for the regression table
saveRDS(final_binary, "Output/final_binary.rds")
saveRDS(final_hours,  "Output/final_hours.rds")
saveRDS(hurdle_binary_ri, "Output/hurdle_binary_ri.rds")
saveRDS(hurdle_binary_rs, "Output/hurdle_binary_rs.rds")
saveRDS(hurdle_hours_ri,  "Output/hurdle_hours_ri.rds")
saveRDS(hurdle_hours_rs,  "Output/hurdle_hours_rs.rds")




# Country-specific results
saveRDS(ranef_binary, "Output/ranef_binary.rds")
saveRDS(ranef_hours,  "Output/ranef_hours.rds")
saveRDS(as.data.frame(gap_employed), "Output/gap_employed.rds")
saveRDS(as.data.frame(gap_hours),    "Output/gap_hours.rds")


# Diagnostics (cached simulations)
saveRDS(sim_binary, "Output/sim_binary.rds")
saveRDS(sim_hours,  "Output/sim_hours.rds")