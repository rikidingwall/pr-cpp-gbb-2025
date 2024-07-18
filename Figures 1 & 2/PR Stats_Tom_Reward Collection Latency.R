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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Run log-adjusted inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency across Date.Time separated by Genotype and Date.Time
ggplot(cleaned_d, aes(x = Date.Time, y = Reward_Collection_Latency, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA) + # Remove outliers
  theme_minimal() +
  labs(title = "Boxplot of Reward Collection Latency by Date and Genotype",
       x = "Date",
       y = "Reward Collection Latency") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  ylim(0, 2) # Set the limits of the Y axis







# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2016-09-28") & Date.Time <= as.Date("2016-09-28"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the inverse gaussian GLM
model <- glmer(Reward_Collection_Latency ~ Genotype + (1 | Animal.ID), 
             data = cleaned_d, 
             family = inverse.gaussian)
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency separated by Genotype
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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2016-09-29") & Date.Time <= as.Date("2016-09-29"))

# Remove rows where Reward_Collection_Latency is blank (NA) and Ratio is not "2"
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Assuming 'Genotype' is a factor with two levels and 'Date.Time' is properly formatted as a Date
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the inverse gaussian GLM
model <- glmer(Reward_Collection_Latency ~ Genotype + (1|Animal.ID), 
             data = cleaned_d, 
             family = inverse.gaussian)
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency separated by Genotype
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







# repeat for FR3, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2016-09-30") & Date.Time <= as.Date("2016-09-30"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the inverse gaussian GLM
model <- glmer(Reward_Collection_Latency ~ Genotype + (1|Animal.ID), 
             data = cleaned_d, 
             family = inverse.gaussian)
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency separated by Genotype
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







# repeat for FR4, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2016-10-01") & Date.Time <= as.Date("2016-10-01"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the inverse gaussian GLM
model <- glmer(Reward_Collection_Latency ~ Genotype + (1|Animal.ID), 
             data = cleaned_d, 
             family = inverse.gaussian)
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency separated by Genotype
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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
specific_dates <- as.Date(c("2016-10-02", "2016-10-03", "2016-10-04", "2016-10-05", "2016-10-21", "2016-10-24")) 
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & 
                 as.Date(d$Date.Time) %in% specific_dates, ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the log-adjusted inverse gaussian GLMM
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency across Date.Time separated by Genotype
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






# repeat for FR5 (before PR), read in data file
cleaned_d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Define the start and end dates
d$Date.Time <- as.Date(d$Date.Time)
start_date <- as.Date("2016-10-02")
end_date <- as.Date("2016-10-05")

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Date.Time >= start_date & d$Date.Time <= end_date, ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the inverse gaussian GLMM 
model <- glmer(Reward_Collection_Latency ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian)
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency across Date.Time separated by Genotype
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






# repeat for FR5 after PR, read in FR data file
cleaned_d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Define the start and end dates
d$Date.Time <- as.Date(d$Date.Time)
start_date <- as.Date("2016-10-21")
end_date <- as.Date("2016-10-24")

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Date.Time >= start_date & d$Date.Time <= end_date, ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns 
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the inverse gaussian GLMM 
model <- glmer(Reward_Collection_Latency ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian)
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmer(Reward_Collection_Latency ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date.Time)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, method = "Wald")
summary(posthoc_results)
print(posthoc_cis)

# Create a boxplot of Reward_Collection_Latency across Date.Time separated by Genotype
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







# repeat for FR10, read in FR data file
cleaned_d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2016-10-25") & Date.Time <= as.Date("2016-10-27"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the inverse gaussian GLMM 
model <- glmer(Reward_Collection_Latency ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian)
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency across Date.Time separated by Genotype
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








# repeat for FR20, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2016-10-28") & Date.Time <= as.Date("2016-10-30"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "0", ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the inverse gaussian GLMM 
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

# Create a boxplot of Reward_Collection_Latency across Date.Time separated by Genotype
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






# repeat for FR40, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2016-10-31") & Date.Time <= as.Date("2016-11-02"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the inverse gaussian GLMM 
model <- glmer(Reward_Collection_Latency ~ Genotype * Date.Time + (1 + Date.Time | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian)
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Create a boxplot of Reward_Collection_Latency across Date.Time separated by Genotype
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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_PR.csv")
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d %>%
  group_by(Animal.ID, Date.Time) %>%
  mutate(Median_Reward_Collection_Latency = median(Reward_Collection_Latency, na.rm = TRUE)) %>%
  ungroup()
cleaned_d <- cleaned_d[!is.na(cleaned_d$Schedule_Length), ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Median_Reward_Collection_Latency")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/tom-pr-reward-col-laten.xlsx")

# Export FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "", ]
cleaned_d <- cleaned_d %>%
  group_by(Animal.ID, Date.Time) %>%
  mutate(Median_Reward_Collection_Latency = median(Reward_Collection_Latency, na.rm = TRUE)) %>%
  ungroup()
cleaned_d <- cleaned_d[!is.na(cleaned_d$Schedule_Length), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Median_Reward_Collection_Latency")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/tom-fr-reward-col-laten.xlsx")
