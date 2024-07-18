# Load necessary libraries
library(lqmm)
library(dplyr)
library(ggplot2)
library(car)
library(DHARMa)
library(lme4)
library(statmod)
library(quantreg)
library(lmtest)
library(emmeans)

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "" & d$Trial == "None", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Run inverse gaussian GLMM 
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency across Date.Time separated by Genotype
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
ggplot(cleaned_d, aes(x = Date.Time, y = Reward_Collection_Latency, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA) + # Remove outliers
  theme_minimal() +
  labs(title = "Boxplot of Reward Collection Latency by Date and Genotype",
       x = "Date",
       y = "Reward Collection Latency") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  ylim(0, 2) # Set the limits of the Y axis






# repeat for PR + ATO, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "" & d$Reward_Collection_Latency != "0" & d$Trial == "ATO", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Trial", "Reward_Collection_Latency", "Drug", "Schedule_Run_ID", "Machine_Name")]

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

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Drug * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Run post-hoc analysis for Genotype x Drug interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Drug * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"), control = glmerControl(optimizer = "bobyqa"))
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Drug)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, method = "Wald")
summary(posthoc_results)
print(posthoc_cis)

# Plotting the boxplot
ggplot(cleaned_d, aes(x = Date.Time, y = Reward_Collection_Latency, fill = Genotype_Drug)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  ylim(0,2) +
  labs(title = "Boxplot of Reward Collection Latency by Date and Genotype/Drug",
       x = "Date.Time",
       y = "Reward Collection Latency")








# repeat for PR + MPH, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Clean data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "" & d$Trial == "MPH", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Trial", "Reward_Collection_Latency", "Drug", "Schedule_Run_ID", "Machine_Name")]

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

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Drug * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Run post-hoc analysis for Genotype x Drug interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Drug + (1 + Date.Time | Animal.ID), data = cleaned_d, family=inverse.gaussian(link="log"))
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Drug)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, method = "Wald")
summary(posthoc_results)
print(posthoc_cis)

# Get the estimated marginal means for the interaction of Genotype and Date.Time
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
model <- glmer(Reward_Collection_Latency ~ Genotype * Drug * Date.Time + (1 | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="identity"),
               control = glmerControl(optimizer = "nloptwrap", optCtrl = list(xtol_rel = 1e-4)))
emm <- emmeans(model, ~ Drug | Date.Time)
pairs <- pairs(emm)
summary(pairs)

# Plotting the boxplot
ggplot(cleaned_d, aes(x = Date.Time, y = Reward_Collection_Latency, fill = Genotype_Drug)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplot of Reward Collection Latency by Date and Genotype/Drug",
       x = "Date.Time",
       y = "Reward Collection Latency")






# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-05-28") & Date.Time <= as.Date("2018-05-29"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Plotting bar chart
genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    median = median(Reward_Collection_Latency, na.rm = TRUE),
    sem = sd(Reward_Collection_Latency, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = median, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = median - sem, ymax = median + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")
  ) +
  labs(title = "Median Reward Collection Latency by Genotype",
       x = "Genotype",
       y = "Median Reward Collection Latency (s)")







# repeat for FR2, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-05-30") & Date.Time <= as.Date("2018-05-31"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Summary of the model
summary(model)

genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    median = median(Reward_Collection_Latency, na.rm = TRUE),
    sem = sd(Reward_Collection_Latency, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = median, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = median - sem, ymax = median + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")
  ) +
  labs(title = "Median Reward Collection Latency by Genotype",
       x = "Genotype",
       y = "Median Reward Collection Latency (s)")






# repeat for FR3, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-01") & Date.Time <= as.Date("2018-06-01"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype + (1 | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Plotting bar chart
genotype_stats <- cleaned_d %>%
  group_by(Genotype) %>%
  summarise(
    median = median(Reward_Collection_Latency, na.rm = TRUE),
    sem = sd(Reward_Collection_Latency, na.rm = TRUE) / sqrt(n())
  )

ggplot(genotype_stats, aes(x = Genotype, y = median, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = median - sem, ymax = median + sem), width = 0.2, size = 0.5, linewidth=1, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, color = "black")
  ) +
  labs(title = "Median Reward Collection Latency by Genotype",
       x = "Genotype",
       y = "Median Reward Collection Latency (s)")






# repeat for FR5 (all sessions), read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Clean up data file
specific_dates <- as.Date(c("2018-06-02", "2018-06-04", "2018-06-05", "2018-06-06", "2018-06-07", "2018-06-08", "2018-06-09", "2018-06-13", "2018-06-18", "2018-06-21", "2018-06-26", "2018-07-02", "2018-07-15", "2018-07-16", "2018-07-19", "2018-07-23", "2018-07-26", "2018-08-02")) 
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & 
                 as.Date(d$Date.Time) %in% specific_dates, ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date.Time)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, method = "Wald")
summary(posthoc_results)
print(posthoc_cis)

# Plotting the barchart 
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    median = median(Reward_Collection_Latency),
    sem = sd(Reward_Collection_Latency)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = median, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = median - sem, ymax = median + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Median Reward Collection Latency by Genotype and Date",
       x = "Date",
       y = "Median Reward Collection Latency")







# repeat for FR5 before PR + ATO, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-02") & Date.Time <= as.Date("2018-06-09"))

# Remove rows where Reward_Collection_Latency is blank (NA) and Ratio is not "5"
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Assuming 'Genotype' is a factor with two levels and 'Date.Time' is properly formatted as a Date
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Plotting a barchart
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    median = median(Reward_Collection_Latency),
    sem = sd(Reward_Collection_Latency)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = median, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = median - sem, ymax = median + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Median Reward Collection Latency by Genotype and Date",
       x = "Date",
       y = "Median Reward Collection Latency")






# repeat for FR5 during PR + ATO, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-13") & Date.Time <= as.Date("2018-07-02"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time) 

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
                       data = cleaned_d, family=inverse.gaussian(link="log"))
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date.Time)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, method = "Wald")
summary(posthoc_results)
print(posthoc_cis)

# Plotting a barchart
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    median = median(Reward_Collection_Latency),
    sem = sd(Reward_Collection_Latency)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = median, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = median - sem, ymax = median + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Median Reward Collection Latency by Genotype and Date",
       x = "Date",
       y = "Median Reward Collection Latency")







# repeat for FR5 during PR + MPH, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-07-16") & Date.Time <= as.Date("2018-08-02"))

# Remove rows where Reward_Collection_Latency is blank (NA) and Ratio is not "5"
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Assuming 'Genotype' is a factor with two levels and 'Date.Time' is properly formatted as a Date
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Run inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Plotting a barchart
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    median = median(Reward_Collection_Latency),
    sem = sd(Reward_Collection_Latency)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = median, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = median - sem, ymax = median + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Median Reward Collection Latency by Genotype and Date",
       x = "Date",
       y = "Median Reward Collection Latency")




# Export PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "", ]
cleaned_d <- cleaned_d %>%
  group_by(Animal.ID, Date.Time) %>%
  mutate(Median_Reward_Collection_Latency = median(Reward_Collection_Latency, na.rm = TRUE)) %>%
  ungroup()
cleaned_d <- cleaned_d[!is.na(cleaned_d$Schedule_Length), ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Median_Reward_Collection_Latency")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-pr-reward-col-laten.xlsx")

# Export FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "", ]
cleaned_d <- cleaned_d %>%
  group_by(Animal.ID, Date.Time) %>%
  mutate(Median_Reward_Collection_Latency = median(Reward_Collection_Latency, na.rm = TRUE)) %>%
  ungroup()
cleaned_d <- cleaned_d[!is.na(cleaned_d$Schedule_Length), ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Median_Reward_Collection_Latency")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-fr-reward-col-laten.xlsx")
