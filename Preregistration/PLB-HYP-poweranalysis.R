

# based on guide: 
# https://easystats.github.io/bayestestR/articles/bayes_factors.html

library(MASS)
library(tidyverse)
library(BayesFactor)
library(HDInterval)
library(pkgcond)
library(bayestestR)
library(lme4)


# predictors of pain
# age + female + height + expected_CPT_pain + FPQ_total

# potential moderators of the effect
# expectancy_hypnoanalgesia + EHS_total + motivation_to_experience_hypnosis

n = 60
sd_outcome = 70
mean_hyp = 136
mean_plb = 100
mean_baseline = 70

sd_age = 3
mean_age = 24
sd_height = 10
mean_height = 178
sd_expected_CPT_pain = 2
mean_expected_CPT_pain = 4
sd_FPQ_total = 1
mean_FPQ_total = 2


varnames_sim_pre = c(
  "outcome_baseline",
  "outcome_hyp",
  "outcome_plb",
  "age",	
  "female",
  "height",
  "expected_CPT_pain",
  "FPQ_total"
)

data_sim_pre = as.data.frame(
  mvrnorm(n = n,
          mu = rep(0, length(varnames_sim_pre)),
          Sigma = matrix(c(1,0.7,0.7,-0.2,-0.2,-0.2,0.2,0.2,
                           0.7,1,0.7,-0.2,-0.2,-0.2,-0.2,-0.2,
                           0.7,0.7,1,-0.2,-0.2,-0.2,-0.2,-0.2,
                           -0.2,-0.2,-0.2,1,0,0,0,0,
                           -0.2,-0.2,-0.2,0,1,-0.2,0,0,
                           -0.2,-0.2,-0.2,0,-0.2,1,0,0,
                           0.2,-0.2,-0.2,0,0,0,1,0,
                           0.2,-0.2,-0.2,0,0,0,0,1), nrow = length(varnames_sim_pre)))
)

names(data_sim_pre) = varnames_sim_pre

data_pre = as.data.frame(matrix(NA, nrow = n, ncol = length(varnames_sim_pre)+1))
names(data_pre) = c("ID", varnames_sim_pre)

data_pre[,"outcome_baseline"] = data_sim_pre[,"outcome_baseline"]*sd_outcome+mean_baseline
data_pre[,"outcome_hyp"] = data_sim_pre[,"outcome_hyp"]*sd_outcome+mean_hyp
data_pre[,"outcome_plb"] = data_sim_pre[,"outcome_plb"]*sd_outcome+mean_plb
data_pre[,"age"] = data_sim_pre[,"age"]*sd_age+mean_age
data_pre[data_sim_pre[,"female"] < 0,"female"] = 0
data_pre[data_sim_pre[,"female"] > -0.00000000001,"female"] = 1
data_pre[,"height"] = data_sim_pre[,"height"]*sd_height+mean_height
data_pre[,"expected_CPT_pain"] = data_sim_pre[,"expected_CPT_pain"]*sd_expected_CPT_pain+mean_expected_CPT_pain
data_pre[,"FPQ_total"] = data_sim_pre[,"FPQ_total"]*sd_FPQ_total+mean_FPQ_total
data_pre[,"ID"] = 1:n

data_long = data_pre %>% gather(key = condition, value = outcome, outcome_hyp:outcome_plb) %>% arrange(ID)
data_long = data_long %>% 
  mutate(condition = as.factor(condition),
         ID = as.factor(ID))

data_long

mod_simple = lm(outcome ~ outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total, data = data_long)
summary(mod_simple)
mod_complex = lm(outcome ~ condition + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total, data = data_long)
summary(mod_complex)

AIC(mod_simple, mod_complex)[1,2]-AIC(mod_simple, mod_complex)[2,2]

