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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Remove outliers
Q1 <- quantile(cleaned_d$Reward_Collection_Latency, 0.25)
Q3 <- quantile(cleaned_d$Reward_Collection_Latency, 0.75)
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR
outliers <- cleaned_d$Reward_Collection_Latency < lower_bound | cleaned_d$Reward_Collection_Latency > upper_bound
cleaned_d <- cleaned_d[!outliers, ]

# Run log-adjusted inverse gaussian GLMM 
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
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = median - sem, ymax = median + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Median Reward Collection Latency by Genotype and Date",
       x = "Date",
       y = "Median Reward Collection Latency")






# repeat for PR only using days with 5 min timeout, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")

# Define the start and end dates
start_date <- as.Date("2017-11-30")
end_date <- as.Date("2017-12-01") 

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "" &
                 d$Date.Time >= start_date & d$Date.Time <= end_date, ]
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

# Plotting a barchart
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    median = median(Reward_Collection_Latency),
    sem = sd(Reward_Collection_Latency)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = median, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = median - sem, ymax = median + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Median Reward Collection Latency by Genotype and Date",
       x = "Date",
       y = "Median Reward Collection Latency")







# repeat for PR only using days without 5 min timeout, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")

# Define the start and end dates
start_date <- as.Date("2017-12-02")
end_date <- as.Date("2017-12-05") 

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "" &
                 d$Date.Time >= start_date & d$Date.Time <= end_date, ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Remove outliers
Q1 <- quantile(cleaned_d$Reward_Collection_Latency, 0.25)
Q3 <- quantile(cleaned_d$Reward_Collection_Latency, 0.75)
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR
outliers <- cleaned_d$Reward_Collection_Latency < lower_bound | cleaned_d$Reward_Collection_Latency > upper_bound
cleaned_d <- cleaned_d[!outliers, ]

# Run log-adjusted inverse gaussian GLMM 
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
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = median - sem, ymax = median + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Median Reward Collection Latency by Genotype and Date",
       x = "Date",
       y = "Median Reward Collection Latency")









# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2017-11-17") & Date.Time <= as.Date("2017-11-23"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

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

# Create a line graph of Reward_Collection_Latency across Date.Time separated by Genotype
grouped_stats <- cleaned_d %>%
  group_by(Genotype, Date.Time) %>%
  summarise(
    median = median(Reward_Collection_Latency),
    sem = sd(Reward_Collection_Latency)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = median, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = median - sem, ymax = median + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Median Reward Collection Latency by Genotype and Date",
       x = "Date",
       y = "Median Reward Collection Latency")







# repeat for FR2, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2017-11-24") & Date.Time <= as.Date("2017-11-24"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Run log-adjusted inverse gaussian GLMM 
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype + (1 | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Plotting a barchart
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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2017-11-27") & Date.Time <= as.Date("2017-11-27"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Run log-adjusted inverse gaussian GLMM 
model <- glmer(log(Reward_Collection_Latency+1) ~ Genotype + (1 | Animal.ID), 
               data = cleaned_d, family=inverse.gaussian(link="log"))
cis <- confint(model, method = "Wald")
summary(model)
print(cis)

# Plotting a barchart
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






# repeat for FR5, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2017-11-28") & Date.Time <= as.Date("2017-11-29"))

# Clean up data file
cleaned_d <- d[!is.na(d$Reward_Collection_Latency), ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Reward_Collection_Latency")]

# Assuming 'Genotype' is a factor with two levels and 'Date.Time' is properly formatted as a Date
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

# Plot line graph separated by Genotype & Date.Time
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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "", ]
cleaned_d <- cleaned_d %>%
  group_by(Animal.ID, Date.Time) %>%
  mutate(Median_Reward_Collection_Latency = median(Reward_Collection_Latency, na.rm = TRUE)) %>%
  ungroup()
cleaned_d <- cleaned_d[!is.na(cleaned_d$Schedule_Length), ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Median_Reward_Collection_Latency")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/carlos-pr-reward-col-laten.xlsx")

# Export FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")
cleaned_d <- d[!is.na(d$Reward_Collection_Latency) & d$Reward_Collection_Latency != "", ]
cleaned_d <- cleaned_d %>%
  group_by(Animal.ID, Date.Time) %>%
  mutate(Median_Reward_Collection_Latency = median(Reward_Collection_Latency, na.rm = TRUE)) %>%
  ungroup()
cleaned_d <- cleaned_d[!is.na(cleaned_d$Schedule_Length), ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Median_Reward_Collection_Latency")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/carlos-fr-reward-col-laten.xlsx")
