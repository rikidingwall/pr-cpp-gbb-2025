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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Magazine_Entries_Rate != "" & d$Trial == "None", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Magazine_Entries_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), data = cleaned_d, family = beta_family)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date.Time)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries_Rate),
    sem = sd(Magazine_Entries_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries Rate by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries Rate")








# repeat for PR + ATO, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Magazine_Entries_Rate != "" & d$Trial == "ATO", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Magazine_Entries_Rate", "Drug", "Schedule_Run_ID", "Machine_Name")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Drug <- as.factor(cleaned_d$Drug)
cleaned_d$Drug <- factor(cleaned_d$Drug, levels = c("Saline", "ATO"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Create interaction terms for 'Genotype' and 'Drug'
cleaned_d$Genotype_Drug <- interaction(cleaned_d$Genotype, cleaned_d$Drug)
cleaned_d$Genotype_Drug <- as.factor(cleaned_d$Genotype_Drug)
cleaned_d$Genotype_Drug <- factor(cleaned_d$Genotype_Drug, levels = c('WT.Saline', 'WT.ATO', 'KI.Saline', 'KI.ATO'))

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Drug * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group and each date
group_summary <- cleaned_d %>%
  group_by(Date.Time, Genotype_Drug) %>%
  summarise(Mean = mean(Magazine_Entries_Rate), SEM = sd(Magazine_Entries_Rate)/sqrt(n())) %>%
  ungroup()

# Plot the line graph
ggplot(group_summary, aes(x = Date.Time, y = Mean, group = Genotype_Drug, color = Genotype_Drug)) +
  geom_line(size = 1.5) + # Bold lines
  geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM, fill = Genotype_Drug), alpha = 0.2) + # SEM ribbon
  scale_color_manual(values = c('black', 'purple', 'grey', 'blue')) + # Custom colors for lines
  scale_fill_manual(values = c('black', 'purple', 'grey', 'blue')) + # Custom colors for ribbons
  labs(title = "Line Graph of Mean Magazine Entries Rate by Date", x = "Date", y = "Mean Magazine Entries Rate") +
  theme_minimal() +
  theme(legend.title = element_blank()) # Remove legend title







# repeat for PR + MPH, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Magazine_Entries_Rate != "" & d$Trial == "MPH", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Magazine_Entries_Rate", "Drug", "Schedule_Run_ID", "Machine_Name")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Drug <- as.factor(cleaned_d$Drug)
cleaned_d$Drug <- factor(cleaned_d$Drug, levels = c("Saline", "MPH"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Create interaction terms for 'Genotype' and 'Drug'
cleaned_d$Genotype_Drug <- interaction(cleaned_d$Genotype, cleaned_d$Drug)
cleaned_d$Genotype_Drug <- as.factor(cleaned_d$Genotype_Drug)
cleaned_d$Genotype_Drug <- factor(cleaned_d$Genotype_Drug, levels = c('WT.Saline', 'WT.MPH', 'KI.Saline', 'KI.MPH'))

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Drug * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group and each date
group_summary <- cleaned_d %>%
  group_by(Date.Time, Genotype_Drug) %>%
  summarise(Mean = mean(Magazine_Entries_Rate), SEM = sd(Magazine_Entries_Rate)/sqrt(n())) %>%
  ungroup()

# Plot the line graph
ggplot(group_summary, aes(x = Date.Time, y = Mean, group = Genotype_Drug, color = Genotype_Drug)) +
  geom_line(size = 1.5) + # Bold lines
  geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM, fill = Genotype_Drug), alpha = 0.2) + # SEM ribbon
  scale_color_manual(values = c('black', 'purple', 'grey', 'blue')) + # Custom colors for lines
  scale_fill_manual(values = c('black', 'purple', 'grey', 'blue')) + # Custom colors for ribbons
  labs(title = "Line Graph of Mean Magazine Entries Rate by Date", x = "Date", y = "Mean Magazine Entries Rate") +
  theme_minimal() +
  theme(legend.title = element_blank()) # Remove legend title








# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Ratio == "1", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)
genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries_Rate, na.rm = TRUE),
    sem = sd(Magazine_Entries_Rate, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Magazine Entries Rate by Genotype",
       x = "Genotype",
       y = "Magazine Entries Rate (#)")








# repeat for FR2, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Ratio == "2", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)
genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries_Rate, na.rm = TRUE),
    sem = sd(Magazine_Entries_Rate, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Magazine Entries Rate by Genotype",
       x = "Genotype",
       y = "Magazine Entries Rate (#)")







# repeat for FR3, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Ratio == "3", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype, 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype
# Create the bar chart with error bars (without time factor)
genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    mean = mean(Magazine_Entries_Rate, na.rm = TRUE),
    sem = sd(Magazine_Entries_Rate, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Magazine Entries Rate by Genotype",
       x = "Genotype",
       y = "Magazine Entries Rate (#)")








# repeat for FR5 (all sessions), read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), data = cleaned_d, family = beta_family)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date.Time)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries_Rate),
    sem = sd(Magazine_Entries_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries Rate by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries Rate")








# repeat for FR5 before PR + ATO, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-02") & Date.Time <= as.Date("2018-06-09"))

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), data = cleaned_d, family = beta_family)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date.Time)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)

grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries_Rate),
    sem = sd(Magazine_Entries_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries Rate by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries Rate")








# repeat for FR5 during PR + ATO, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-13") & Date.Time <= as.Date("2018-07-02"))

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries_Rate),
    sem = sd(Magazine_Entries_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries Rate by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries Rate")







# repeat for FR5 during PR + MPH, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-07-16") & Date.Time <= as.Date("2018-08-02"))

# Create Magazine Entry Rate and clean up data file
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]

# Assuming 'Genotype' is a factor with two levels and 'Date.Time' is properly formatted as a Date
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Magazine_Entries_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEMs for each Genotype at each Date
# Create line graph with error ribbon (with time factor)
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    mean = mean(Magazine_Entries_Rate),
    sem = sd(Magazine_Entries_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Magazine Entries Rate by Genotype and Date",
       x = "Date",
       y = "Mean Magazine Entries Rate")





# Export PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate) & d$Magazine_Entries_Rate != "" & d$Trial == "None", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-pr-magaz-entry-rate.xlsx")

# Export FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d$Magazine_Entries_Rate <- d$Magazine_Entries / d$Schedule_Length
cleaned_d <- d[!is.na(d$Magazine_Entries_Rate), ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Magazine_Entries_Rate")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-fr-magaz-entry-rate.xlsx")
