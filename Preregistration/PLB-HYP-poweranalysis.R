

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

simulate_data_with_moderators = function(
                                                           n,
                                                           mean_hyp,
                                                           mean_plb,
                                                           mean_baseline,
                                                           sd_outcome,
                                                           sd_age,
                                                           mean_age,
                                                           sd_height,
                                                           mean_height,
                                                           sd_expected_CPT_pain,
                                                           mean_expected_CPT_pain,
                                                           sd_FPQ_total,
                                                           mean_FPQ_total,
                                                           sd_expectancy_pain_reduction,
                                                           mean_expectancy_pain_reduction_hyp,
                                                           mean_expectancy_pain_reduction_plb,
                                                           sd_expectancy_hypnotic_depth,
                                                           mean_expectancy_hypnotic_depth_hyp,
                                                           mean_expectancy_hypnotic_depth_plb,
                                                           sd_hypnotic_depth,
                                                           mean_hypnotic_depth_hyp,
                                                           mean_hypnotic_depth_plb,
                                                           sd_hypnotizability,
                                                           mean_hypnotizability
                                         ){
  varnames_sim_pre = c(
    "outcome_baseline", "outcome_hyp", "outcome_plb", "age", "female", "height", "expected_CPT_pain", "FPQ_total", "expectancy_pain_reduction_hyp", "expectancy_pain_reduction_plb", "expectancy_hypnotic_depth_hyp", "expectancy_hypnotic_depth_plb", "hypnotic_depth_hyp", "hypnotic_depth_plb", "hypnotizability"
  )
  
  data_sim_pre = as.data.frame(
    mvrnorm(n = n,
            mu = rep(0, length(varnames_sim_pre)),
            Sigma = matrix(c(1,0.5,0.5,-0.2,-0.2,-0.2,0.2,0.2,0,0,0,0,0,0,0,
                             0.5,1,0.7,-0.2,-0.2,-0.2,-0.2,-0.2,0.5,0.3,0.4,0.2,0.5,0.3,0.4,
                             0.5,0.7,1,-0.2,-0.2,-0.2,-0.2,-0.2,0.3,0.5,0.2,0.4,0.3,0.5,0.4,
                             -0.2,-0.2,-0.2,1,0,0,0,0,0,0,0,0,0,0,0,
                             -0.2,-0.2,-0.2,0,1,-0.2,0,0,0,0,0,0,0,0,0,
                             -0.2,-0.2,-0.2,0,-0.2,1,0,0,0,0,0,0,0,0,0,
                             0.2,-0.2,-0.2,0,0,0,1,0,0,0,0,0,0,0,0,
                             0.2,-0.2,-0.2,0,0,0,0,1,0,0,0,0,0,0,0,
                             0,0.5,0.3,0,0,0,0,0,1,0.6,0.5,0.4,0.3,0.3,0.4,
                             0,0.3,0.5,0,0,0,0,0,0.6,1,0.4,0.5,0.3,0.3,0.4,
                             0,0.4,0.2,0,0,0,0,0,0.5,0.4,1,0.6,0.4,0.3,0.4,
                             0,0.2,0.4,0,0,0,0,0,0.4,0.5,0.6,1,0.3,0.4,0.4,
                             0,0.5,0.3,0,0,0,0,0,0.3,0.3,0.4,0.3,1,0.6,0.6,
                             0,0.3,0.5,0,0,0,0,0,0.3,0.3,0.3,0.4,0.6,1,0.6,
                             0,0.4,0.4,0,0,0,0,0,0.4,0.4,0.4,0.4,0.6,0.6,1
                             
            ), nrow = length(varnames_sim_pre)))
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
  data_pre[,"expectancy_pain_reduction_hyp"] = data_sim_pre[,"expectancy_pain_reduction_hyp"]*sd_expectancy_pain_reduction+mean_expectancy_pain_reduction_hyp
  data_pre[,"expectancy_pain_reduction_plb"] = data_sim_pre[,"expectancy_pain_reduction_plb"]*sd_expectancy_pain_reduction+mean_expectancy_pain_reduction_plb
  data_pre[,"expectancy_hypnotic_depth_hyp"] = data_sim_pre[,"expectancy_hypnotic_depth_hyp"]*sd_expectancy_hypnotic_depth+mean_expectancy_hypnotic_depth_hyp
  data_pre[,"expectancy_hypnotic_depth_plb"] = data_sim_pre[,"expectancy_hypnotic_depth_plb"]*sd_expectancy_hypnotic_depth+mean_expectancy_hypnotic_depth_plb
  data_pre[,"hypnotic_depth_hyp"] = data_sim_pre[,"hypnotic_depth_hyp"]*sd_hypnotic_depth+mean_hypnotic_depth_hyp
  data_pre[,"hypnotic_depth_plb"] = data_sim_pre[,"hypnotic_depth_plb"]*sd_hypnotic_depth+mean_hypnotic_depth_plb
  data_pre[,"hypnotizability"] = data_sim_pre[,"hypnotizability"]*sd_hypnotizability+mean_hypnotizability
  
  
  data = data_pre %>% 
    gather(key = "condition", value = "outcome", outcome_hyp:outcome_plb) %>%
    gather(key = "condition_exp_painred", value = "expectancy_pain_reduction", expectancy_pain_reduction_hyp:expectancy_pain_reduction_plb) %>%
    filter(substring(condition, 9, 11) == substring(condition_exp_painred, 27, 29)) %>% 
    gather(key = "condition_exp_hypdepth", value = "expectancy_hypnotic_depth", expectancy_hypnotic_depth_hyp:expectancy_hypnotic_depth_plb) %>%
    filter(substring(condition, 9, 11) == substring(condition_exp_hypdepth, 27, 29)) %>%
    gather(key = "condition_hypdepth", value = "hypnotic_depth", hypnotic_depth_hyp:hypnotic_depth_plb) %>%
    filter(substring(condition, 9, 11) == substring(condition_hypdepth, 16, 18)) %>%
    arrange(ID)

  
  
  
  data = data %>% 
    mutate(condition = as.factor(condition),
           ID = as.factor(ID),
           condition_hyp = recode(condition,
                                  "outcome_hyp" = 1,
                                  "outcome_plb" = 0))
  return(data)
}





within_subjects_analysis = function(data,
                                    SESOI,
                                    BF_threshold,
                                    AIC_threshold,
                                    decision_criterion){
  ### models with lme4

  if(decision_criterion == "CI_lme4"){
    mod_complex = lmer(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data)
    CI_lme4 = confint(mod_complex)["condition_hyp",]
    decision = if(CI_lme4[2] < SESOI){"H0"} else if (CI_lme4[1] > 0){"H1"} else {"inconclusive"} 
  } else if(decision_criterion == "AIC_lme4"){
    mod_simple = lmer(outcome ~ outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data)
    mod_complex = lmer(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data)
    cAIC_coparison_lme4 = cAIC(mod_simple)$caic - cAIC(mod_complex)$caic
    decision = if(cAIC_coparison_lme4 > AIC_threshold){"H1"} else if(cAIC_coparison_lme4 < -AIC_threshold){"H0"} else {"inconclusive"}
  } else if(decision_criterion == "BF_lme4"){
    mod_simple = lmer(outcome ~ outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data)
    mod_complex = lmer(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + (1|ID), data = data)
    BF_comparison_lme4 = as.numeric(bayesfactor_models(mod_simple, mod_complex))[2]
    decision = if(BF_comparison_lme4 > BF_threshold){"H1"} else if(BF_comparison_lme4 < 1/BF_threshold){"H0"} else {"inconclusive"}
    
    
    ### models with lmBF    
    
  } else if(decision_criterion == "BF_lmBF"){
    mod_simple_BF = lmBF(outcome ~ outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + ID, whichRandom="ID", data = data)
    mod_complex_BF = lmBF(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + ID, whichRandom="ID", data = data)
    BF_comparison_lmBF = extractBF(mod_complex_BF/mod_simple_BF)[1,1]
    decision = if(BF_comparison_lmBF > BF_threshold){"H1"} else if(BF_comparison_lmBF < 1/BF_threshold){"H0"} else {"inconclusive"}
  } else if(decision_criterion == "HDI_lmBF"){
    mod_complex_BF = lmBF(outcome ~ condition_hyp + outcome_baseline + age + female + height + expected_CPT_pain + FPQ_total + ID, whichRandom="ID", data = data)
    posteriors = posterior(mod_complex_BF, iterations = 1000)
    HDI = summary(posteriors)[[2]]["condition_hyp",c(1,5)]
    decision = if(HDI[1] > SESOI){"H1"} else if(HDI[2] < SESOI){"H0"} else {"inconclusive"}
  }
  
  return(decision)
}


simul = function(
                  n,
                  mean_hyp,
                  mean_plb,
                  mean_baseline,
                  sd_outcome,
                  sd_age,
                  mean_age,
                  sd_height,
                  mean_height,
                  sd_expected_CPT_pain,
                  mean_expected_CPT_pain,
                  sd_FPQ_total,
                  mean_FPQ_total,
                  sd_expectancy_pain_reduction,
                  mean_expectancy_pain_reduction_hyp,
                  mean_expectancy_pain_reduction_plb,
                  sd_expectancy_hypnotic_depth,
                  mean_expectancy_hypnotic_depth_hyp,
                  mean_expectancy_hypnotic_depth_plb,
                  sd_hypnotic_depth,
                  mean_hypnotic_depth_hyp,
                  mean_hypnotic_depth_plb,
                  sd_hypnotizability,
                  mean_hypnotizability,
                 
                  BF_threshold,
                  AIC_threshold,
                  SESOI,
                  decision_criterion
){

  data = simulate_data_with_moderators(
                  n = n,
                  mean_hyp = mean_hyp,
                  mean_plb = mean_plb,
                  mean_baseline = mean_baseline,
                  sd_outcome = sd_outcome,
                  sd_age = sd_age,
                  mean_age = mean_age,
                  sd_height = sd_height,
                  mean_height = mean_height,
                  sd_expected_CPT_pain = sd_expected_CPT_pain,
                  mean_expected_CPT_pain = mean_expected_CPT_pain,
                  sd_FPQ_total = sd_FPQ_total,
                  mean_FPQ_total = mean_FPQ_total,
                  sd_expectancy_pain_reduction = sd_expectancy_pain_reduction,
                  mean_expectancy_pain_reduction_hyp = mean_expectancy_pain_reduction_hyp,
                  mean_expectancy_pain_reduction_plb = mean_expectancy_pain_reduction_plb,
                  sd_expectancy_hypnotic_depth = sd_expectancy_hypnotic_depth,
                  mean_expectancy_hypnotic_depth_hyp = mean_expectancy_hypnotic_depth_hyp,
                  mean_expectancy_hypnotic_depth_plb = mean_expectancy_hypnotic_depth_plb,
                  sd_hypnotic_depth = sd_hypnotic_depth,
                  mean_hypnotic_depth_hyp = mean_hypnotic_depth_hyp,
                  mean_hypnotic_depth_plb = mean_hypnotic_depth_plb,
                  sd_hypnotizability = sd_hypnotizability,
                  mean_hypnotizability = mean_hypnotizability
  )
  
  decision = within_subjects_analysis(data = data,
                                       SESOI = SESOI,
                                       BF_threshold = BF_threshold,
                                       AIC_threshold = AIC_threshold,
                                       decision_criterion = decision_criterion)
  
  return(decision)

}



################
# Hypothesis I #
################

# H1 is simulated to be true

iter = 100

results = replicate(iter, simul(  n = 50,
                                  mean_hyp = 136,
                                  mean_plb = 100,
                                  mean_baseline = 70,
                                  sd_outcome = 70,
                                  sd_age = 3,
                                  mean_age = 24,
                                  sd_height = 10,
                                  mean_height = 178,
                                  sd_expected_CPT_pain = 2,
                                  mean_expected_CPT_pain = 4,
                                  sd_FPQ_total = 1,
                                  mean_FPQ_total = 2,
                                  sd_expectancy_pain_reduction = 2,
                                  mean_expectancy_pain_reduction_hyp = 5,
                                  mean_expectancy_pain_reduction_plb = 5,
                                  sd_expectancy_hypnotic_depth = 2,
                                  mean_expectancy_hypnotic_depth_hyp = 5,
                                  mean_expectancy_hypnotic_depth_plb = 5,
                                  sd_hypnotic_depth = 2,
                                  mean_hypnotic_depth_hyp = 5,
                                  mean_hypnotic_depth_plb = 5,
                                  sd_hypnotizability = 2,
                                  mean_hypnotizability = 6,
                                  
                                  SESOI = 15,
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

