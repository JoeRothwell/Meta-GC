# Read in participant data and pos and neg data and join
dat <- read_sas("meta_gc_caco.sas7bdat")
pos <- read_sas("meta_gc_untg_rp_pos_metabo.sas7bdat") %>% select(-(Idepic_Bio:Cncr_Caco_Stom))
neg <- read_sas("meta_gc_untg_rp_neg_metabo.sas7bdat") %>%
  select(-(Idepic_Bio:Center)) %>% select(-Match_Caseset, -Cncr_Caco_Stom)

posneg <- bind_cols(pos, neg)

# Both modes ----
# Join participant data to metabolomics data
library(tidyverse)
posdat <- right_join(dat, posneg, by = "Idepic")

# Keep CC pairs with no HPPOS missings
nonmiss.hp <- posdat %>% group_by(Match_Caseset) %>% filter(!anyNA(HPPOS)) %>% ungroup() # 400

# Condordant case-control pairs
concordant11 <- nonmiss.hp %>% group_by(Match_Caseset) %>% filter(sum(HPPOS) == 2) %>% ungroup() # 210/2 = 105
concordant00 <- nonmiss.hp %>% group_by(Match_Caseset) %>% filter(sum(HPPOS) == 0) %>% ungroup() # 22/2 = 11

# Discordant case-control pairs: all
discordant   <- nonmiss.hp %>% group_by(Match_Caseset) %>% filter(sum(HPPOS) == 1) %>% ungroup() # 168/2 = 84
table(discordant$HPPOS, discordant$Cncr_Caco_Stom) # 67 and 17 of each discordance type

# Make new variable to distinguish 2 discordance types
discordant <- discordant %>% mutate(disc.type = ifelse(HPPOS == 1 & Cncr_Caco_Stom == 1, 1, 0))

# Get discordance types
discordant10 <- discordant %>% group_by(Match_Caseset) %>% filter(sum(disc.type) == 1) %>% ungroup() # 134/2 = 67
discordant01 <- discordant %>% group_by(Match_Caseset) %>% filter(sum(disc.type) == 0) %>% ungroup() # 34/2 = 17

# Identify CC pairs with at least one missing
anymiss.hp <- posdat %>% group_by(Match_Caseset) %>% filter(is.na(sum(HPPOS))) %>% ungroup() # 474/2 = 237

# Identify CC pairs with both missings
allmiss.hp <- posdat %>% group_by(Match_Caseset) %>% filter(sum(is.na(HPPOS))==2) %>% ungroup() # 460/2 = 230

# Identify CC pairs with an HPPOS == 1 and a missing
miss.hp1 <- anymiss.hp %>% group_by(Match_Caseset) %>% filter(sum(HPPOS, na.rm = TRUE) == 1) %>% ungroup() # 10/2 = 5
table(HP = miss.hp1$HPPOS, case = miss.hp1$Cncr_Caco_Stom)
# 5 pairs have HPPOS either in the case or the control, of which 4 pairs have HPPOS for the case

# Make new variable to distinguish 2 discordance types
miss.hp1 <- miss.hp1 %>% mutate(disc.type = ifelse(HPPOS == 1 & Cncr_Caco_Stom == 1, 1, 0))

# Get discordance types
miss.disc10 <- miss.hp1 %>% group_by(Match_Caseset) %>% filter(sum(disc.type, na.rm = TRUE) == 1) %>% ungroup() # 8/2 = 4
miss.disc01 <- miss.hp1 %>% group_by(Match_Caseset) %>% filter(sum(disc.type, na.rm = TRUE) == 0) %>% ungroup() # 1/2 = 1

