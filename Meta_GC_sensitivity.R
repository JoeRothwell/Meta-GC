# MetaGC sensitivity analyses (as suggested by Mazda). 
# Run appropriate commands first from MetaGC_hp_subsets.R to get each relevant data subset

# (a)	run sensitivity analyses only on the n=105 case-control pairs with Hppos = 1,
# Subset metabolomics data only using prefix. Log transform data.
# Model with essential covariates
library(tidyverse)
ints <- concordant11 %>% select(starts_with("Untg_Rp")) %>% as.matrix() %>% log2()

# Get feature names binding together mode and m/z from data frame labels
mode <- c(rep("pos", ncol(pos)), rep("neg", 2689-ncol(pos)))
features <- sapply(posneg[ , 1:ncol(ints)], attr, "label") %>% unname
featnames <- paste(mode, features, sep = "_")

clr <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + QE_ENERGY + QE_ALC + 
                                L_School + Fasting_C + strata(Match_Caseset), data = concordant11)

library(survival)
varlist <- c("Smoke_Stat", "Alc_Drinker", "Center", "L_School", "Fasting_C")
concordant11 <- concordant11 %>% mutate(across((varlist), as.factor))
# Apply model across subset
mod.adj2 <- apply(ints, 2, clr)

library(broom)
res.adj2 <- map_df(mod.adj2, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.adj2))

# Format OR and CI
tab.adj <- res.adj2 %>%
  select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, conf.low, hyph, conf.high, B2, sep = "")

# No significant features after p-value adjustment

# (b)	run sensitivity analyses on the n=105 + n=71 (i.e., n=176) case-control pairs where the case is Hppos=1, 
# but the control could be 1 or 0.

# This is concordant11 + discordant10 (168 + 134 obs = actually 151 pairs)
discordant10 <- discordant10 %>% mutate(across((varlist), as.factor)) %>% select(-disc.type)
caseposall <- bind_rows(concordant11, discordant10)
ints <- caseposall %>% select(starts_with("Untg_Rp")) %>% as.matrix() %>% log2()

clr1 <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + 
                                strata(Match_Caseset), data = caseposall)

# Apply model across subset
mod.adj3 <- apply(ints, 2, clr1)

res.adj3 <- map_df(mod.adj3, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.adj3))
# No significant features after p-value adjustment

# Format OR and CI
tab.adj <- res.adj3 %>%
  select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, conf.low, hyph, conf.high, B2, sep = "")


# (c)	Possibly consider an analysis on the n=11 + n=17 where the cases are Hppos negative, but this will be weak.

# This is concordant00 + discordant01 (28 pairs)
casenegall <- bind_rows(concordant00, discordant01)
ints <- casenegall %>% select(starts_with("Untg_Rp")) %>% as.matrix() %>% log2()

clr2 <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + 
                                QE_ENERGY + QE_ALC + L_School + Fasting_C + 
                                strata(Match_Caseset), data = casenegall)

# Apply model across subset
pos.adj <- apply(ints, 2, clr2)

res.adj <- map_df(pos.adj, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(pos.adj))

# Format OR and CI
tab.adj <- res.adj3 %>%
  select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, conf.low, hyph, conf.high, B2, sep = "")


# (d)	Break matching and run unconditional analyses by Hppos, but with adjustment for country instead of by centre.

# Matching factors are study centre, age at blood collection, sex, fasting status, 
# time of blood collection, menopausal status, exogenous hormone use, phase of menstrual cycle
# Model with essential covariates

# Better to leave unknown HPPOS status
#posdat$HPPOS2 <- as_factor(posdat$HPPOS)
#posdat$HPPOS2 <- fct_explicit_na(posdat$HPPOS2)

# Unmatched LR model adjusting additionally for Hppos and matching factors where possible: sex, country,
# hormone use, fasting status
glm.hpp <- function(x) glm(Cncr_Caco_Stom ~ x + Sex + Country + #Use_Horm + 
                             Bmi_C + Smoke_Stat + Pa_Total + 
                                Age_Blood + QE_ENERGY + QE_ALC + L_School + Fasting_C + HPPOS,
                                data = posdat, family = "binomial")

# Reset posdat from MetaGC_hp_subsets.R
ints <- posdat %>% select(starts_with("Untg_Rp")) %>% as.matrix() %>% log2()

# Convert categorical variables to factors
varlist <- c("Smoke_Stat", "Alc_Drinker", "Country", "L_School", "Fasting_C", "Use_Horm", "Sex", "HPPOS")
posdat <- posdat %>% mutate(across((varlist), as.factor))

# Apply model across subset
pos.adj <- apply(ints, 2, glm.hpp)

#res.adj <- map_df(pos.adj, tidy, exponentiate = T, conf.int = T) %>% 
res.adj <- map_df(pos.adj, tidy,  exponentiate = T) %>% 
  filter(term == "x") %>% 
  mutate(OR = exp(estimate), p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(pos.adj))

# Format OR and CI
tab.adj <- res.adj %>%
  select(feature, estimate, everything()) %>%
  #select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:OR), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, hyph, B2, sep = "")


# (e)	In all sub-group analyses a, b, c and d above, also run analyses on the subset of n=230 who have missing data, 
# possibly with heterogeneity analyses?

ints <- allmiss.hp %>% select(starts_with("Untg_Rp")) %>% as.matrix() %>% log2()

clr <- function(x) clogit(Cncr_Caco_Stom ~ x + Bmi_C + Smoke_Stat + Pa_Total + QE_ENERGY + QE_ALC + 
                                L_School + Fasting_C + strata(Match_Caseset), data = allmiss.hp)

# Apply model across subset
library(survival)
mod.adj3 <- apply(ints, 2, clr)

library(broom)
res.adj3 <- map_df(mod.adj3, tidy, exponentiate = T, conf.int = T) %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"), feature = features) %>%
  select(feature, everything(), -term) %>% mutate(feat.no = 1:length(mod.adj3))

# Format OR and CI
tab.adj <- res.adj3 %>%
  select(feature, estimate, conf.low, conf.high, everything()) %>%
  mutate_at(vars(estimate:conf.high), ~ format(round(., digits = 2), nsmall = 2)) %>% 
  mutate(B1 = " (", hyph = "-", B2 = ")") %>%
  unite("CI95", B1, conf.low, hyph, conf.high, B2, sep = "")
  