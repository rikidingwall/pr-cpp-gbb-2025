# Load the dependencies
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(emmeans)
library(multcomp)
library(multcompView)

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Magazine_Entries != "", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the negative binomial regression 
model <- glmmTMB(Magazine_Entries ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)

grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries),
    sem = sd(Magazine_Entries)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries")









# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up date file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "1", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries")]

# Fit the negative binomial regression 
model <- glmmTMB(Magazine_Entries ~ Genotype, 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)

genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries, na.rm = TRUE),
    sem = sd(Magazine_Entries, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")) +
  labs(title = "Mean Magazine Entries by Genotype",
       x = "Genotype",
       y = "Magazine Entries (#)")







# repeat for FR2, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up date file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "2", ]

# Covnert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries")]

# Fit the negative binomial regression 
model <- glmmTMB(Magazine_Entries ~ Genotype, 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)

genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries, na.rm = TRUE),
    sem = sd(Magazine_Entries, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")) +
  labs(title = "Mean Magazine Entries by Genotype",
       x = "Genotype",
       y = "Magazine Entries (#)")








# repeat for FR3, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "3", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries")]

# Fit the poisson regression 
model <- glmmTMB(Magazine_Entries ~ Genotype, 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)

genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries, na.rm = TRUE),
    sem = sd(Magazine_Entries, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")
  ) +
  labs(title = "Mean Magazine Entries by Genotype",
       x = "Genotype",
       y = "Magazine Entries (#)")






# repeat for FR4, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "4", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the poisson regression 
model <- glmmTMB(Magazine_Entries ~ Genotype, 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)

genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries, na.rm = TRUE),
    sem = sd(Magazine_Entries, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")
  ) +
  labs(title = "Mean Magazine Entries by Genotype",
       x = "Genotype",
       y = "Magazine Entries (#)")








# repeat for FR5 (all sessions), read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "5", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the negative binomial regression 
model <- glmmTMB(Magazine_Entries ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)

grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries),
    sem = sd(Magazine_Entries)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries")






# repeat for FR5 (before PR), read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Define the start and end dates
d$Date.Time <- as.Date(d$Date.Time)
start_date <- as.Date("2016-10-02")
end_date <- as.Date("2016-10-05")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "5" & d$Date.Time >= start_date & d$Date.Time <= end_date, ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the poisson regression 
model <- glmmTMB(Magazine_Entries ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)

grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries),
    sem = sd(Magazine_Entries)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries")






# repeat for FR5 (after PR), read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Define the start and end dates
d$Date.Time <- as.Date(d$Date.Time)
start_date <- as.Date("2016-10-21")
end_date <- as.Date("2016-10-24")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "5" & d$Date.Time >= start_date & d$Date.Time <= end_date, ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the poisson regression 
model <- glmmTMB(Magazine_Entries ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = poisson)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)

grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries),
    sem = sd(Magazine_Entries)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries")







# repeat for FR10, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "10", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the negative binomial regression 
model <- glmmTMB(Magazine_Entries ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)

genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries, na.rm = TRUE),
    sem = sd(Magazine_Entries, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")
  ) +
  labs(title = "Mean Magazine Entries by Genotype",
       x = "Genotype",
       y = "Magazine Entries (#)")






# repeat for FR20, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "20", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the negative binomial regression 
model <- glmmTMB(Magazine_Entries ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmmTMB(Magazine_Entries ~ Genotype * Date.Time + (1|Animal.ID), data = cleaned_d, family = nbinom2)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date.Time)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)

genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries, na.rm = TRUE),
    sem = sd(Magazine_Entries, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")
  ) +
  labs(title = "Mean Magazine Entries by Genotype",
       x = "Genotype",
       y = "Magazine Entries (#)")






# repeat for FR40, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Magazine_Entries) & d$Ratio == "40", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the negative binomial regression 
model <- glmmTMB(Magazine_Entries ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)

genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries, na.rm = TRUE),
    sem = sd(Magazine_Entries, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = mean, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")
  ) +
  labs(title = "Mean Magazine Entries by Genotype",
       x = "Genotype",
       y = "Magazine Entries (#)")
