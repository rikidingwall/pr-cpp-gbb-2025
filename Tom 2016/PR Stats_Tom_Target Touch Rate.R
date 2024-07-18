# Load the dependencies
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(emmeans)
library(multcomp)
library(multcompView)
library(gamlss)
library(broom.mixed)
library(openxlsx)

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_PR.csv")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Target_Touches_Rate != "", ]
cleaned_d <- cleaned_d[, c( "Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Target_Touches_Rate),
    sem = sd(Target_Touches_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Target Touches Rate by Genotype and Date",
       x = "Date",
       y = "Mean Target Touches Rate")










# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "1", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype, 
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
    mean = mean(Target_Touches_Rate, na.rm = TRUE),
    sem = sd(Target_Touches_Rate, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Target Touches Rate by Genotype",
       x = "Genotype",
       y = "Target Touches Rate (#)")







# repeat for FR2, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Remove rows where Target_Touches_Rate is blank (NA) and Ratio is not "2"
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "2", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Assuming 'Genotype' is a factor with two levels and 'Date.Time' is properly formatted as a Date
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype, 
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
    mean = mean(Target_Touches_Rate, na.rm = TRUE),
    sem = sd(Target_Touches_Rate, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Target Touches Rate by Genotype",
       x = "Genotype",
       y = "Target Touches Rate (#)")








# repeat for FR3, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Create Target Touch Rate and clean up data file 
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "3", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Animal.ID <- as.factor(cleaned_d$Animal.ID)

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype, 
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
    mean = mean(Target_Touches_Rate, na.rm = TRUE),
    sem = sd(Target_Touches_Rate, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Target Touches Rate by Genotype",
       x = "Genotype",
       y = "Target Touches Rate (#)")






# repeat for FR4, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "4", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Animal.ID <- as.factor(cleaned_d$Animal.ID)

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype, 
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
    mean = mean(Target_Touches_Rate, na.rm = TRUE),
    sem = sd(Target_Touches_Rate, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Target Touches Rate by Genotype",
       x = "Genotype",
       y = "Target Touches Rate (#)")








# repeat for FR5 (all sessions), read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert column
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Target_Touches_Rate),
    sem = sd(Target_Touches_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Target Touches Rate by Genotype and Date",
       x = "Date",
       y = "Mean Target Touches Rate")








# repeat for FR5 (before PR), read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Define the start and end dates
d$Date.Time <- as.Date(d$Date.Time)
start_date <- as.Date("2016-10-02")
end_date <- as.Date("2016-10-05")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "5" & d$Date.Time >= start_date & d$Date.Time <= end_date, ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert column
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Target_Touches_Rate),
    sem = sd(Target_Touches_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Target Touches Rate by Genotype and Date",
       x = "Date",
       y = "Mean Target Touches Rate")






# repeat for FR5 after PR, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Define the start and end dates
d$Date.Time <- as.Date(d$Date.Time)
start_date <- as.Date("2016-10-21")
end_date <- as.Date("2016-10-24")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "5" & d$Date.Time >= start_date & d$Date.Time <= end_date, ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert column
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Target_Touches_Rate),
    sem = sd(Target_Touches_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Target Touches Rate by Genotype and Date",
       x = "Date",
       y = "Mean Target Touches Rate")









# repeat for FR10, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "10", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Target_Touches_Rate),
    sem = sd(Target_Touches_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Target Touches Rate by Genotype and Date",
       x = "Date",
       y = "Mean Target Touches Rate")








# repeat for FR20, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "20", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d <- cleaned_d[cleaned_d$Target_Touches_Rate != 0, ] # Removes unexpected 0 value - likely system error

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Target_Touches_Rate),
    sem = sd(Target_Touches_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Target Touches Rate by Genotype and Date",
       x = "Date",
       y = "Mean Target Touches Rate")









# repeat for FR40, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Create Target Touch Rate and clean up data file
d$Target_Touches_Rate <- d$Target_Touches / d$Schedule_Length
cleaned_d <- d[!is.na(d$Target_Touches_Rate) & d$Ratio == "40", ]
cleaned_d <- cleaned_d[, c("Schedule_Length", "Genotype", "Animal.ID", "Date.Time", "Target_Touches_Rate")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the beta regression 
model <- glmmTMB(Target_Touches_Rate ~ Genotype * Date.Time + (1|Animal.ID), 
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
    mean = mean(Target_Touches_Rate),
    sem = sd(Target_Touches_Rate)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Target Touches Rate by Genotype and Date",
       x = "Date",
       y = "Mean Target Touches Rate")
