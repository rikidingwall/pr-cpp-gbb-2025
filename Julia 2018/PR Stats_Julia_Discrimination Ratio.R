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

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Discrimination_Ratio != "" & d$Trial == "None", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio", "Trial", "Drug")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Animal.ID <- as.factor(cleaned_d$Animal.ID)

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
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
    mean = mean(Discrimination_Ratio),
    sem = sd(Discrimination_Ratio)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Discrimination Ratio by Genotype and Date",
       x = "Date",
       y = "Mean Discrimination Ratio")








# repeat for PR + ATO, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Discrimination_Ratio != "" & d$Trial == "ATO", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio", "Trial", "Drug")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Animal.ID <- as.factor(cleaned_d$Animal.ID)
cleaned_d$Drug <- as.factor(cleaned_d$Drug)
cleaned_d$Drug <- factor(cleaned_d$Drug, levels = c("Saline", "ATO"))

# Create interaction terms for 'Genotype' and 'Drug'
cleaned_d$Genotype_Drug <- interaction(cleaned_d$Genotype, cleaned_d$Drug)
cleaned_d$Genotype_Drug <- as.factor(cleaned_d$Genotype_Drug)
cleaned_d$Genotype_Drug <- factor(cleaned_d$Genotype_Drug, levels = c('WT.Saline', 'WT.ATO', 'KI.Saline', 'KI.ATO'))

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Drug * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group and each date
group_summary <- cleaned_d %>%
  group_by(Date.Time, Genotype_Drug) %>%
  summarise(Mean = mean(Discrimination_Ratio), SEM = sd(Discrimination_Ratio)/sqrt(n())) %>%
  ungroup()

# Plot the line graph
ggplot(group_summary, aes(x = Date.Time, y = Mean, group = Genotype_Drug, color = Genotype_Drug)) +
  geom_line(size = 1.5) + # Bold lines
  geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM, fill = Genotype_Drug), alpha = 0.2) + # SEM ribbon
  scale_color_manual(values = c('black', 'purple', 'grey', 'blue')) + # Custom colors for lines
  scale_fill_manual(values = c('black', 'purple', 'grey', 'blue')) + # Custom colors for ribbons
  labs(title = "Line Graph of Mean Breakpoint by Date", x = "Date", y = "Mean Breakpoint") +
  theme_minimal() +
  theme(legend.title = element_blank()) # Remove legend title








# repeat for PR + MPH, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Discrimination_Ratio != "" & d$Trial == "MPH", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio", "Trial", "Drug")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Animal.ID <- as.factor(cleaned_d$Animal.ID)
cleaned_d$Drug <- as.factor(cleaned_d$Drug)
cleaned_d$Drug <- factor(cleaned_d$Drug, levels = c("Saline", "MPH"))

# Create interaction terms for 'Genotype' and 'Drug'
cleaned_d$Genotype_Drug <- interaction(cleaned_d$Genotype, cleaned_d$Drug)
cleaned_d$Genotype_Drug <- as.factor(cleaned_d$Genotype_Drug)
cleaned_d$Genotype_Drug <- factor(cleaned_d$Genotype_Drug, levels = c('WT.Saline', 'WT.MPH', 'KI.Saline', 'KI.MPH'))

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Drug * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Calculate means and SEM for each group and each date
group_summary <- cleaned_d %>%
  group_by(Date.Time, Genotype_Drug) %>%
  summarise(Mean = mean(Discrimination_Ratio), SEM = sd(Discrimination_Ratio)/sqrt(n())) %>%
  ungroup()

# Plot the line graph
ggplot(group_summary, aes(x = Date.Time, y = Mean, group = Genotype_Drug, color = Genotype_Drug)) +
  geom_line(size = 1.5) + # Bold lines
  geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM, fill = Genotype_Drug), alpha = 0.2) + # SEM ribbon
  scale_color_manual(values = c('black', 'purple', 'grey', 'blue')) + # Custom colors for lines
  scale_fill_manual(values = c('black', 'purple', 'grey', 'blue')) + # Custom colors for ribbons
  labs(title = "Line Graph of Mean Breakpoint by Date", x = "Date", y = "Mean Breakpoint") +
  theme_minimal() +
  theme(legend.title = element_blank()) # Remove legend title







# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Ratio == "1", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
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
    mean = mean(Discrimination_Ratio, na.rm = TRUE),
    sem = sd(Discrimination_Ratio, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Discrimination Ratio by Genotype",
       x = "Genotype",
       y = "Discrimination Ratio (#)")








# repeat for FR2, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Ratio == "2", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
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
    mean = mean(Discrimination_Ratio, na.rm = TRUE),
    sem = sd(Discrimination_Ratio, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Discrimination Ratio by Genotype",
       x = "Genotype",
       y = "Discrimination Ratio (#)")








# repeat for FR3, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Ratio == "3", ]
cleaned_d <- cleaned_d[, c("Ratio", "Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype, 
                 ziformula = ~1,
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
    mean = mean(Discrimination_Ratio, na.rm = TRUE),
    sem = sd(Discrimination_Ratio, na.rm = TRUE) / sqrt(n())
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
  labs(title = "Mean Discrimination Ratio by Genotype",
       x = "Genotype",
       y = "Discrimination Ratio (#)")








# repeat for FR5, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Create Discrimination Ratio and read in data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
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
    mean = mean(Discrimination_Ratio),
    sem = sd(Discrimination_Ratio)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Discrimination Ratio by Genotype and Date",
       x = "Date",
       y = "Mean Discrimination Ratio")








# repeat for FR5 before PR + ATO, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-02") & Date.Time <= as.Date("2018-06-09"))

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
                 data = cleaned_d, 
                 family = beta_family)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmmTMB(Discrimination_Ratio ~ Genotype * Date.Time + (1|Animal.ID), ziformula = ~1, data = cleaned_d, family = beta_family)
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
    mean = mean(Discrimination_Ratio),
    sem = sd(Discrimination_Ratio)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting
ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Discrimination Ratio by Genotype and Date",
       x = "Date",
       y = "Mean Discrimination Ratio")








# repeat for FR5 during PR + ATO, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-13") & Date.Time <= as.Date("2018-07-02"))

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
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
    mean = mean(Discrimination_Ratio),
    sem = sd(Discrimination_Ratio)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Discrimination Ratio by Genotype and Date",
       x = "Date",
       y = "Mean Discrimination Ratio")








# repeat for FR5 during PR + MPH, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-07-16") & Date.Time <= as.Date("2018-08-02"))

# Create Discrimination Ratio and clean up data file
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Ratio == "5", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the zero-inflated beta regression 
model <- glmmTMB(Discrimination_Ratio ~ Genotype * Date.Time + (1|Animal.ID), 
                 ziformula = ~1,
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
    mean = mean(Discrimination_Ratio),
    sem = sd(Discrimination_Ratio)/sqrt(n())
  ) %>%
  ungroup() # Ungroup for plotting

ggplot(grouped_stats, aes(x = Date.Time, y = mean, group = Genotype, color = Genotype)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Genotype), alpha = 0.2) +
  scale_color_manual(values = c("WT" = "black", "KI" = "blue")) +
  scale_fill_manual(values = c("WT" = "black", "KI" = "blue")) +
  theme_minimal() +
  labs(title = "Mean Discrimination Ratio by Genotype and Date",
       x = "Date",
       y = "Mean Discrimination Ratio")





# Export PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio) & d$Discrimination_Ratio != "" & d$Trial == "None", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio", "Drug")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-pr-discrim-ratio.xlsx")

# Export FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d$Discrimination_Ratio <- (d$Blank_Touches / (d$Blank_Touches + d$Target_Touches))
cleaned_d <- d[!is.na(d$Discrimination_Ratio), ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Discrimination_Ratio")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-fr-discrim-ratio.xlsx")
