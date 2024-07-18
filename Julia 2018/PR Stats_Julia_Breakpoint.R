# Load necessary libraries
library(lme4)
library(ggplot2)
library(dplyr)

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "" & d$Trial == "None", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Breakpoint")]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)

# Fit the negative binomial regression 
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








# repeat for PR + ATO, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "" & d$Trial == "ATO", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Breakpoint", "Drug", "Schedule_Run_ID", "Machine_Name")]

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

# Fit the negative binomial regression 
model <- glmmTMB(Breakpoint ~ Genotype * Drug * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
posthoc_model <- glmmTMB(Breakpoint ~ Genotype * Drug * Date.Time + (1|Animal.ID), data = cleaned_d, family = nbinom2)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype * Drug)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Run post-hoc analysis for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- glmmTMB(Breakpoint ~ Genotype * Drug * Date.Time + (1|Animal.ID), data = cleaned_d, family = nbinom2)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Drug)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEM for each group and each date
group_summary <- cleaned_d %>%
  group_by(Date.Time, Genotype_Drug) %>%
  summarise(Mean = mean(Breakpoint), SEM = sd(Breakpoint)/sqrt(n())) %>%
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

# Clean up data file
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "" & d$Trial == "MPH", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Breakpoint", "Drug", "Schedule_Run_ID", "Machine_Name")]

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

# Fit the negative binomial regression 
model <- glmmTMB(Breakpoint ~ Genotype * Drug * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
cis <- confint(model, level = 0.95)
summary(model)
print(cis)

# Run post-hoc analysis for interaction effect
posthoc_model <- glmmTMB(Breakpoint ~ Genotype * Drug * Date.Time + (1|Animal.ID), 
                 data = cleaned_d, 
                 family = nbinom2)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype * Drug)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, level = 0.95)
summary(posthoc_results)
print(posthoc_cis)

# Calculate means and SEM for each group and each date
group_summary <- cleaned_d %>%
  group_by(Date.Time, Genotype_Drug) %>%
  summarise(Mean = mean(Breakpoint), SEM = sd(Breakpoint)/sqrt(n())) %>%
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




# Export PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "" & d$Trial == "None", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Breakpoint")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-pr-breakpoint.xlsx")

# Export PR + ATO data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "" & d$Trial == "ATO", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Breakpoint")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-pr-ato-breakpoint.xlsx")

# Export PR + MPH data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")
cleaned_d <- d[!is.na(d$Breakpoint) & d$Breakpoint != "" & d$Trial == "MPH", ]
cleaned_d <- cleaned_d[, c("Genotype", "Animal.ID", "Date.Time", "Trial", "Breakpoint")]
write_xlsx(cleaned_d, "/Users/rikidingwall/Downloads/julia-pr-mph-breakpoint.xlsx")
