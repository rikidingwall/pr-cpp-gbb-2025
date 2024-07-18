# Load necessary libraries
library(lme4)
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(emmeans)

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Breakpoint")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit GLMM with negative binomial regression
model <- glmmTMB(Breakpoint ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Breakpoint),
    sem = sd(Breakpoint)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Breakpoint by Genotype and Date",
       x = "Date",
       y = "Mean Breakpoint")





# repeat for PR with 5 min timeout, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")

# Define the start and end dates
start_date <- as.Date("2017-11-30")
end_date <- as.Date("2017-12-01") 

# Clean up data file
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "" &
                 d$Date.Time >= start_date & d$Date.Time <= end_date, ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Breakpoint")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit GLMM with negative binomial regression
model <- glmmTMB(Breakpoint ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Breakpoint),
    sem = sd(Breakpoint)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Breakpoint by Genotype and Date",
       x = "Date",
       y = "Mean Breakpoint")









# repeat for PR only using days without 5 min timeout, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")

# Define the start and end dates
start_date <- as.Date("2017-12-02")
end_date <- as.Date("2017-12-05") 

# Clean up data file
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "" &
                 d$Date.Time >= start_date & d$Date.Time <= end_date, ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Breakpoint")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit GLMM with negative binomial regression
model <- glmmTMB(Breakpoint ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Breakpoint),
    sem = sd(Breakpoint)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Breakpoint by Genotype and Date",
       x = "Date",
       y = "Mean Breakpoint")






# Export PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Breakpoint")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/carlos-pr-breakpoint.xlsx")
