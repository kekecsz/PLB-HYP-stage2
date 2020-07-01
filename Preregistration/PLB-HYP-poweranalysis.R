

# based on guide: 
# https://easystats.github.io/bayestestR/articles/bayes_factors.html

library(MASS)
library(tidyverse)
library(BayesFactor)
library(HDInterval)
library(pkgcond)
library(bayestestR)
library(lme4)
library(cAIC4)


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


BF_threshold = 3
AIC_threshold = 2
SESOI = 15
decision_criterion = "BF_lmBF" # either "CI_lme4", "AIC_lme4", "BF_lme4", "HDI_lmBF", "BF_lmBF"




simul = function(n,
                 mean_hyp,
                 mean_plb,
                 mean_baseline,
                 sd_outcome,
                 SESOI,
                 
                 sd_age,
                 mean_age,
                 sd_height,
                 mean_height,
                 sd_expected_CPT_pain,
                 mean_expected_CPT_pain,
                 sd_FPQ_total,
                 mean_FPQ_total,
                 
                 BF_threshold,
                 AIC_threshold,
                 
                 decision_criterion
){
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
            Sigma = matrix(c(1,0.5,0.5,-0.2,-0.2,-0.2,0.2,0.2,
                             0.5,1,0.7,-0.2,-0.2,-0.2,-0.2,-0.2,
                             0.5,0.7,1,-0.2,-0.2,-0.2,-0.2,-0.2,
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
           ID = as.factor(ID),
           condition_hyp = recode(condition,
                                  "outcome_hyp" = 1,
                                  "outcome_plb" = 0))
  
  ### models with lme4
  
  if(decision_criterion == "CI_lme4"){
    mod_complex = lmer(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data_long)
    CI_lme4 = confint(mod_complex)["condition_hyp",]
    decision = if(CI_lme4[2] < SESOI){"H0"} else if (CI_lme4[1] > 0){"H1"} else {"inconclusive"} 
  } else if(decision_criterion == "AIC_lme4"){
    mod_simple = lmer(outcome ~ outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data_long)
    mod_complex = lmer(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data_long)
    cAIC_coparison_lme4 = cAIC(mod_simple)$caic - cAIC(mod_complex)$caic
    decision = if(cAIC_coparison_lme4 > AIC_threshold){"H1"} else if(cAIC_coparison_lme4 < -AIC_threshold){"H0"} else {"inconclusive"}
  } else if(decision_criterion == "BF_lme4"){
    mod_simple = lmer(outcome ~ outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data_long)
    mod_complex = lmer(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data_long)
    BF_comparison_lme4 = as.numeric(bayesfactor_models(mod_simple, mod_complex))[2]
    decision = if(BF_comparison_lme4 > BF_threshold){"H1"} else if(BF_comparison_lme4 < 1/BF_threshold){"H0"} else {"inconclusive"}
    
    
    ### models with lmBF    
    
  } else if(decision_criterion == "BF_lmBF"){
    mod_simple_BF = lmBF(outcome ~ outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + ID, whichRandom="ID", data = data_long)
    mod_complex_BF = lmBF(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + ID, whichRandom="ID", data = data_long)
    BF_comparison_lmBF = extractBF(mod_complex_BF/mod_simple_BF)[1,1]
    decision = if(BF_comparison_lmBF > BF_threshold){"H1"} else if(BF_comparison_lmBF < 1/BF_threshold){"H0"} else {"inconclusive"}
  } else if(decision_criterion == "HDI_lmBF"){
    mod_complex_BF = lmBF(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + ID, whichRandom="ID", data = data_long)
    posteriors = posterior(mod_complex_BF, iterations = 1000)
    HDI = summary(posteriors)[[2]]["condition_hyp",c(1,5)]
    decision = if(HDI[1] > SESOI){"H1"} else if(HDI[2] < SESOI){"H0"} else {"inconclusive"}
  }
  
  return(decision)
}




iter = 100

results = replicate(iter, simul(  n = 50,
                                  mean_hyp = 136,
                                  mean_plb = 100,
                                  mean_baseline = 70,
                                  sd_outcome = 70,
                                  SESOI = 15,
                                  
                                  sd_age = 3,
                                  mean_age = 24,
                                  sd_height = 10,
                                  mean_height = 178,
                                  sd_expected_CPT_pain = 2,
                                  mean_expected_CPT_pain = 4,
                                  sd_FPQ_total = 1,
                                  mean_FPQ_total = 2,
                                  
                                  BF_threshold = 3,
                                  AIC_threshold = 2,

                                  decision_criterion = "BF_lmBF" # either "CI_lme4", "AIC_lme4", "BF_lme4", "HDI_lmBF", "BF_lmBF"
                    ))

                                  

H1_detection_rate = mean(results=="H1")
H0_detection_rate = mean(results=="H0")
Inconclusive_rate = mean(results=="inconclusive")

H1_detection_rate
H0_detection_rate
Inconclusive_rate

