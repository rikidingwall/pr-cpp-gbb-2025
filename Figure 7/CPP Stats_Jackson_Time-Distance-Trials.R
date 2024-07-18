# load in libraries
library(dplyr)
library(lme4)
library(emmeans)

# load in data file 
d <- read.csv("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Conditioned Place Preference/Jackson_CPP.csv")

# filter data file for habituation
cleaned_d <- d %>% 
  filter(Stage == "Habituation")
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "NL3"))

# LM for Total Distance
distance_model <- lm(Total_Distance ~ Genotype, 
             data = cleaned_d)
cis <- confint(distance_model, level = 0.95)
summary(distance_model)
print(cis)

distance_shapiro <- shapiro.test(resid(distance_model))
print(distance_shapiro)

# LM for Time Preference
time_model <- lm(Time_Difference ~ Genotype, 
                 data = cleaned_d)
cis <- confint(time_model, level = 0.95)
summary(time_model)
print(cis)

time_shapiro <- shapiro.test(resid(time_model))
print(time_shapiro)

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_cpp-habit.xlsx")






# repeat for Saline Conditioning, load in data file 
d <- read.csv("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Conditioned Place Preference/Jackson_CPP.csv")

# filter data file for habituation
cleaned_d <- d %>% 
  filter(Stage == "Saline_Conditioning")
cleaned_d$Date <- as.Date(cleaned_d$Date)
cleaned_d$Date <- scale(cleaned_d$Date)
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "NL3"))

# LMM for Total Distance
distance_model <- lmer(Total_Distance ~ Genotype * Date + (1 | Animal_ID), 
                   data = cleaned_d)
cis <- confint(distance_model, level = 0.95)
summary(distance_model)
print(cis)

distance_shapiro <- shapiro.test(resid(distance_model))
print(distance_shapiro)







# repeat for Cocaine Conditioning, load in data file 
d <- read.csv("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Conditioned Place Preference/Jackson_CPP.csv")

# filter data file for habituation
cleaned_d <- d %>% 
  filter(Stage == "Cocaine_Conditioning")
cleaned_d$Date <- as.Date(cleaned_d$Date)
cleaned_d$Date <- scale(cleaned_d$Date)
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "NL3"))

# GLMM for Total Distance
distance_model <- lmer(Total_Distance ~ Genotype * Date + (1 | Animal_ID), 
                       data = cleaned_d)
cis <- confint(distance_model, level = 0.95)
summary(distance_model)
print(cis)

distance_shapiro <- shapiro.test(resid(distance_model))
print(distance_shapiro)

# Posthoc test for interaction effect
cleaned_d$Date <- as.factor(cleaned_d$Date)
posthoc_model <- lmer(Total_Distance ~ Genotype * Date + (1 | Animal_ID), 
                       data = cleaned_d)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)








# repeat for CPP, load in data file 
d <- read.csv("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Conditioned Place Preference/Jackson_CPP.csv")

# filter data file for habituation
cleaned_d <- d %>% 
  filter(Stage == "Conditioned_Place_Preference")
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "NL3"))

# GLM for Total Distance
distance_model <- lm(Total_Distance ~ Genotype, 
                      data = cleaned_d)
cis <- confint(distance_model, level = 0.95)
summary(distance_model)
print(cis)

distance_shapiro <- shapiro.test(resid(distance_model))
print(distance_shapiro)

# GLM for Time Preference
time_model <- lm(Time_Difference ~ Genotype, 
                  data = cleaned_d, 
                  family = gaussian)
cis <- confint(time_model, level = 0.95)
summary(time_model)
print(cis)

time_shapiro <- shapiro.test(resid(time_model))
print(time_shapiro)

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_cpp-cpp-test.xlsx")







# repeat for EPP, load in data file 
d <- read.csv("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Conditioned Place Preference/Jackson_CPP.csv")

# filter data file for extinction test
cleaned_d <- d %>% 
  filter(Stage == "Extinction_Test")
cleaned_d$Date <- as.Date(cleaned_d$Date)
cleaned_d$Date <- scale(cleaned_d$Date)
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "NL3"))

# GLMM for Total Distance
distance_model <- lmer(Total_Distance ~ Genotype * Date + (1 | Animal_ID), 
                       data = cleaned_d)
cis <- confint(distance_model, level = 0.95)
summary(distance_model)
print(cis)

distance_shapiro <- shapiro.test(resid(distance_model))
print(distance_shapiro)

# GLMM for Time Preference 
time_model <- lmer(Time_Difference ~ Genotype * Date + (1 | Animal_ID), 
                       data = cleaned_d)
cis <- confint(time_model, level = 0.95)
summary(time_model)
print(cis)

time_shapiro <- shapiro.test(resid(time_model))
print(time_shapiro)

# filter data file for progression extinction trial
cleaned_d <- cleaned_d %>%
  group_by(Animal_ID) %>%
  mutate(Number_of_Extinction_Trials = n()) %>%
  ungroup()
cleaned_d <- cleaned_d %>%
  arrange(Animal_ID, desc(Date)) %>%
  group_by(Animal_ID) %>%
  slice(1)

# LM for Total Distance
distance_model <- lm(Total_Distance ~ Genotype, 
                      data = cleaned_d)
cis <- confint(distance_model, level = 0.95)
summary(distance_model)
print(cis)

distance_shapiro <- shapiro.test(resid(distance_model))
print(distance_shapiro)

# LM for Time Preference
time_model <- lm(Time_Difference ~ Genotype, 
                  data = cleaned_d, 
                  family = gaussian)
cis <- confint(time_model, level = 0.95)
summary(time_model)
print(cis)

time_shapiro <- shapiro.test(resid(time_model))
print(time_shapiro)


# GLM for Number of Extinction Trials 
extinction_model <- glm(Number_of_Extinction_Trials ~ Genotype, 
                      data = cleaned_d, 
                      family = poisson)
cis <- confint(extinction_model, level = 0.95)
summary(extinction_model)
print(cis)

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_cpp-epp-test.xlsx")






# repeat for Reinstatement, load in data file 
d <- read.csv("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Conditioned Place Preference/Jackson_CPP.csv")

# filter data file for habituation
cleaned_d <- d %>% 
  filter(Stage == "Reinstatement")
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "NL3"))

# GLM for Total Distance
distance_model <- lm(Total_Distance ~ Genotype, 
                      data = cleaned_d)
cis <- confint(distance_model, level = 0.95)
summary(distance_model)
print(cis)

distance_shapiro <- shapiro.test(resid(distance_model))
print(distance_shapiro)

# GLM for Time Preference
time_model <- glm(Time_Difference ~ Genotype, 
                  data = cleaned_d, 
                  family = gaussian)
cis <- confint(time_model, level = 0.95)
summary(time_model)
print(cis)

time_shapiro <- shapiro.test(resid(time_model))
print(time_shapiro)

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_cpp-reinstat.xlsx")

