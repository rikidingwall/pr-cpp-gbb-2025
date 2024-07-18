# Load the dependencies
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(patchwork)
library(emmeans)

# Read in CSV for FR1 SL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Riki_FR_PR.csv")

# Only keep rows where Stage is "FR1 SL"
cleaned_d <- d[d$Stage == "FR1 SL" & d$Progression == "1", ]

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Position <- as.factor(cleaned_d$Position)
cleaned_d$Progression_Session <- scale(cleaned_d$Progression_Session)
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)

# Fit GLMM with poisson regression for active lever presses
active_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(active_model, level = 0.95)
summary(active_model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Progression_Session) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype, active_grouped_data$Position, sep = "."),
                             levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))

# Plot active lever presses
ggplot(active_grouped_data, aes(x = Progression_Session, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Progression Session", y = "Active Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/riki-fr1.xlsx")






# Read in CSV for FR1 DL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Riki_FR_PR.csv")

# Only keep rows where Stage is "FR1 DL"
cleaned_d <- d[d$Stage == "FR1 DL" & d$Progression == "1", ]

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Position <- as.factor(cleaned_d$Position)
cleaned_d$Progression_Session <- scale(cleaned_d$Progression_Session)
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit a poisson GLMM for active lever presses
active_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID),
                 data = cleaned_d,
                 family = poisson)
cis <- confint(active_model, level = 0.95)
summary(active_model)
print(cis)

# Fit a poisson GLMM for inactive lever presses
inactive_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID),
                          data = cleaned_d,
                          family = poisson)
cis <- confint(inactive_model, level = 0.95)
summary(inactive_model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Progression_Session <- as.factor(cleaned_d$Progression_Session)
posthoc_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID),
                         data = cleaned_d,
                         family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Position | Progression_Session)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Progression_Session) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Progression_Session) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype, active_grouped_data$Position, sep = "."),
                             levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype, inactive_grouped_data$Position, sep = "."),
                                    levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Progression_Session, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Progression Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Progression_Session, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Progression Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/riki-fr1-dl.xlsx")







# Read in CSV for FR3 DL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Riki_FR_PR.csv")

# Only keep rows where Stage is "FR3"
cleaned_d <- d[d$Stage == "FR3" & d$Progression == "1", ]

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Position <- as.factor(cleaned_d$Position)
cleaned_d$Progression_Session <- scale(cleaned_d$Progression_Session)
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit a poisson GLMM for active lever presses
active_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID),
                        data = cleaned_d,
                        family = poisson)
cis <- confint(active_model, level = 0.95)
summary(active_model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Progression_Session <- as.factor(cleaned_d$Progression_Session)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID),
                        data = cleaned_d,
                        family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Progression_Session)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Fit a poisson GLMM for inactive lever presses
inactive_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID),
                          data = cleaned_d,
                          family = poisson)
cis <- confint(inactive_model, level = 0.95)
summary(inactive_model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Progression_Session) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Progression_Session) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype, active_grouped_data$Position, sep = "."),
                                    levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype, inactive_grouped_data$Position, sep = "."),
                                      levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Progression_Session, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Progression Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Progression_Session, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Progression Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/riki-fr3.xlsx")








# Read in CSV for FR5 DL
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Riki_FR_PR.csv")

# Only keep rows where Stage is "FR5"
cleaned_d <- d[d$Stage == "FR5" & d$Progression == "1", ]

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Position <- as.factor(cleaned_d$Position)
cleaned_d$Progression_Session <- scale(cleaned_d$Progression_Session)
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit a poisson GLMM for active lever presses
active_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID),
                        data = cleaned_d,
                        family = poisson)
cis <- confint(active_model, level = 0.95)
summary(active_model)
print(cis)

# Fit a poisson GLMM for inactive lever presses
inactive_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position + Position * Progression_Session + Genotype * Progression_Session + (1|Animal.ID),
                          data = cleaned_d,
                          family = poisson)
cis <- confint(inactive_model, level = 0.95)
summary(inactive_model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Progression_Session) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Progression_Session) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype, active_grouped_data$Position, sep = "."),
                                    levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype, inactive_grouped_data$Position, sep = "."),
                                      levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Progression_Session, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Progression Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Progression_Session, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Progression Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/riki-fr5.xlsx")








# Read in CSV for PR 2HR
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Riki_FR_PR.csv")

# Only keep rows where Stage is "PR 2HR"
cleaned_d <- d[d$Stage == "PR 2HR", ]

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Position <- as.factor(cleaned_d$Position)
cleaned_d$Session_Number <- scale(cleaned_d$Session_Number)
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit a poisson GLMM for active lever presses
active_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position * Session_Number + (1|Animal.ID),
                        data = cleaned_d,
                        family = poisson)
cis <- confint(active_model, level = 0.95)
summary(active_model)
print(cis)

# Run post-hoc analysis for genotype by session interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position * Session_Number + (1|Animal.ID),
                        data = cleaned_d,
                        family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Run post-hoc analysis for session by lever location interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position * Session_Number + (1|Animal.ID),
                         data = cleaned_d,
                         family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Position | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Fit a poisson GLMM for active lever presses collapsing sessions
active_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position  + (1|Animal.ID),
                        data = cleaned_d,
                        family = poisson)
cis <- confint(active_model, level = 0.95)
summary(active_model)
print(cis)

# Fit a poisson GLMM for inactive lever presses
inactive_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position * Session_Number + (1|Animal.ID),
                          data = cleaned_d,
                          family = poisson)
cis <- confint(inactive_model, level = 0.95)
summary(inactive_model)
print(cis)

# Run post-hoc analysis for genotype by session interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position * Session_Number + (1|Animal.ID),
                          data = cleaned_d,
                          family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Run post-hoc analysis for session by lever location interaction effect
cleaned_d$Session_Number <- as.factor(cleaned_d$Session_Number)
posthoc_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position * Session_Number + (1|Animal.ID),
                         data = cleaned_d,
                         family = poisson)
emm_int <- emmeans(posthoc_model, pairwise ~ Position | Session_Number)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Fit a poisson GLMM for inactive lever presses collapsing sessions
inactive_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position + (1|Animal.ID),
                          data = cleaned_d,
                          family = poisson)
cis <- confint(inactive_model, level = 0.95)
summary(inactive_model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Session_Number) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position, Session_Number) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype, active_grouped_data$Position, sep = "."),
                                    levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype, inactive_grouped_data$Position, sep = "."),
                                      levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))

# Plot active and inactive lever presses
active.plot <- active <- ggplot(active_grouped_data, aes(x = Session_Number, y = Active_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Active Lever Presses (Mean ± SEM)") +
  theme(legend.position = "none", legend.title = element_blank())

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Session_Number, y = Inactive_Mean, group = Group, color = Group)) +
  geom_line() +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), width = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                                "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Session", y = "Inactive Lever Presses (Mean ± SEM)", 
       color = "Group") +
  theme(legend.title = element_blank())

active.plot + inactive.plot

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/riki-pr-2hr.xlsx")









# Read in CSV for PR 6HR
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Operant/Riki_FR_PR.csv")

# Only keep rows where Stage is "PR 6HR"
cleaned_d <- d[d$Stage == "PR 6HR", ]

# Convert variables
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels= c("WT", "KI"))
cleaned_d$Position <- as.factor(cleaned_d$Position)
cleaned_d$Active_Lever_Presses <- as.numeric(cleaned_d$Active_Lever_Presses)
cleaned_d$Inactive_Lever_Presses <- as.numeric(cleaned_d$Inactive_Lever_Presses)

# Fit a negative binomial GLMM for active lever presses collapsing sessions
active_model <- glmmTMB(Active_Lever_Presses ~ Genotype * Position,
                        data = cleaned_d,
                        family = nbinom2)
cis <- confint(active_model, level = 0.95)
summary(active_model)
print(cis)

# Fit a negative binomial GLMM for inactive lever presses collapsing sessions
inactive_model <- glmmTMB(Inactive_Lever_Presses ~ Genotype * Position,
                          data = cleaned_d,
                          family = nbinom2)
cis <- confint(inactive_model, level = 0.95)
summary(inactive_model)
print(cis)

# Calculate means and SEM for each group
active_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position) %>%
  summarise(
    Active_Mean = mean(Active_Lever_Presses),
    Active_SEM = sd(Active_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

inactive_grouped_data <- cleaned_d %>%
  group_by(Genotype, Position) %>%
  summarise(
    Inactive_Mean = mean(Inactive_Lever_Presses),
    Inactive_SEM = sd(Inactive_Lever_Presses)/sqrt(n())
  ) %>%
  ungroup()

# Create a factor with levels in the specific order for plotting
active_grouped_data$Group <- factor(paste(active_grouped_data$Genotype, active_grouped_data$Position, sep = "."),
                                    levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))
inactive_grouped_data$Group <- factor(paste(inactive_grouped_data$Genotype, inactive_grouped_data$Position, sep = "."),
                                      levels = c("WT.Opposite", "WT.Beside", "KI.Opposite", "KI.Beside"))

# Plot active and inactive lever presses
active.plot <- ggplot(active_grouped_data, aes(x = Group, y = Active_Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Active_Mean - Active_SEM, ymax = Active_Mean + Active_SEM), 
                position = position_dodge(0.9), width = 0.25) +
  scale_fill_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                               "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Group", y = "Active Lever Presses", fill = "Group") +
  theme(legend.position = "none")

inactive.plot <- ggplot(inactive_grouped_data, aes(x = Group, y = Inactive_Mean, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Inactive_Mean - Inactive_SEM, ymax = Inactive_Mean + Inactive_SEM), 
                position = position_dodge(0.9), width = 0.25) +
  scale_fill_manual(values = c("WT.Opposite" = "black", "WT.Beside" = "purple", 
                               "KI.Opposite" = "grey", "KI.Beside" = "blue")) +
  theme_minimal() +
  labs(x = "Group", y = "Inactive Lever Presses", fill = "Group") +
  theme(legend.position = "none")

active.plot + inactive.plot

# Export data file
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/riki-pr-6hr.xlsx")
