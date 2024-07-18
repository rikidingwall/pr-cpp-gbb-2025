# Load the dependencies
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(patchwork)
library(emmeans)
library(DHARMa)

# Read in CSV for FR1 SL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR1 SL"
cleaned_d <- d[d$Stage == "FR1 SL", ]

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)

# Fit GLMM with poisson regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))

# Plot active lever presses
ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_fr1-sl.xlsx")





# Read in CSV for FR1 DL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR1 DL" and remove faulty inactive lever values
cleaned_d <- d[d$Stage == "FR1 DL", ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1501] <- NA

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with poisson regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Fit GLMM with poisson regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Convert faulty lever values to -1 (for removal before plotting) and export data file
cleaned_d$Inactive_Lever_Presses[is.na(cleaned_d$Inactive_Lever_Presses)] <- -1
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_fr1-dl-test.xlsx")






# Read in CSV for FR3 DL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR3"
cleaned_d <- d[d$Stage == "FR3", ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1501] <- NA

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with poisson regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Fit GLMM with poisson regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Convert faulty lever values to -1 (for removal before plotting) and export data file
cleaned_d$Inactive_Lever_Presses[is.na(cleaned_d$Inactive_Lever_Presses)] <- -1
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_fr3.xlsx")








# Read in CSV for FR5 DL (all sessions)
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR5"
cleaned_d <- d[d$Stage == "FR5", ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1501] <- NA

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with poisson regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Fit GLMM with negative binomial regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Convert faulty lever values to -1 (for removal before plotting) and export data file
cleaned_d$Inactive_Lever_Presses[is.na(cleaned_d$Inactive_Lever_Presses)] <- -1
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_fr5-all.xlsx")








# Read in CSV for FR5 DL (before PR)
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR5"
cleaned_d <- d[d$Stage == "FR5" & d$Session_Number <= 8, ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1501] <- NA

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with poisson regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Fit GLMM with poisson regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot








# Read in CSV for FR5 DL, baselining FR5
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR5"
cleaned_d <- d[d$Stage == "FR5" & d$Session_Number > 8, ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1501] <- NA

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with negative binomial regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Fit GLMM with poisson regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot










# Read in CSV for PR 2HR
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "PR 2HR"
cleaned_d <- d[d$Stage == "PR 2HR", ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1] <- NA 

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with poisson regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Fit GLMM with poisson regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Session_Number + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Convert faulty lever values to -1 (for removal before plotting) and export data file
cleaned_d$Inactive_Lever_Presses[is.na(cleaned_d$Inactive_Lever_Presses)] <- -1
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_pr.xlsx")









# Read in CSV for PR 2HR + ATO
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "PR 2HR + ATO"
cleaned_d <- d[d$Stage == "PR 2HR + ATO", ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1] <- NA 

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Dose <- as.numeric(cleaned_d$Dose)
cleaned_d$Drug <- as.factor(cleaned_d$Drug)
cleaned_d$Drug <- factor(cleaned_d$Drug, levels= c("Saline", "ATO"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)
cleaned_d$Drug_Dose <- paste(cleaned_d$Drug, cleaned_d$Dose, sep = ".")

# Fit GLMM with negative binomial regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype * Drug + Genotype * Dose + (1 | Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for dose effect
cleaned_d$Dose <- as.factor(cleaned_d$Dose)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Drug + Genotype * Dose + (1 | Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
emm_int <- emmeans(posthoc_model, pairwise ~ Dose)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Fit GLMM with negative binomial regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Drug + Genotype * Dose + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate the mean and SEM of each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Drug_Dose) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Drug_Dose) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create the bar chart
active.plot <- ggplot(active_grouped_data, aes(x = Drug_Dose, y = Active_Mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), 
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Drug and Dose", y = "Mean Active Lever Presses", fill = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Drug_Dose, y = Inactive_Mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), 
                position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Drug and Dose", y = "Mean Inactive Lever Presses", fill = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

active.plot + inactive.plot

# Convert faulty lever values to -1 (for removal before plotting) and export data file
cleaned_d$Inactive_Lever_Presses[is.na(cleaned_d$Inactive_Lever_Presses)] <- -1
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_pr-ato.xlsx")






# Read in CSV for FR10 DL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR10"
cleaned_d <- d[d$Stage == "FR10", ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1501] <- NA 

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with negative binomial regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype, 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Fit GLMM with negative binomial regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype, 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Convert faulty lever values to -1 (for removal before plotting) and export data file
cleaned_d$Inactive_Lever_Presses[is.na(cleaned_d$Inactive_Lever_Presses)] <- -1
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_fr10.xlsx")







# Read in CSV for FR20 DL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR20"
cleaned_d <- d[d$Stage == "FR20", ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1501] <- NA 

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with negative binomial regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype, 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Fit GLMM with negative binomial regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype, 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Convert faulty lever values to -1 (for removal before plotting) and export data file
cleaned_d$Inactive_Lever_Presses[is.na(cleaned_d$Inactive_Lever_Presses)] <- -1
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_fr20.xlsx")







# Read in CSV for FR40 DL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Jackson_2017_PR.csv")

# Only keep rows where Stage is "FR40"
cleaned_d <- d[d$Stage == "FR40", ]
cleaned_d$Inactive_Lever_Presses[cleaned_d$Inactive_Lever_Presses == 1501] <- NA 

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit GLMM with negative binomial regression for active lever presses
model <- glmmTMB(Active_Lever_Presses ~ Genotype, 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Fit GLMM with negative binomial regression for inactive lever presses
model <- glmmTMB(Inactive_Lever_Presses ~ Genotype, 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype),
                                    levels = c("WT", "KI"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype),
                                      levels = c("WT", "KI"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Convert faulty lever values to -1 (for removal before plotting) and export data file
cleaned_d$Inactive_Lever_Presses[is.na(cleaned_d$Inactive_Lever_Presses)] <- -1
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/jackson_2017_fr40.xlsx")
