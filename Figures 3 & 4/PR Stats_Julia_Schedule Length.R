# Load the dependencies
library(survival)
library(survminer)
library(ggplot2)
library(emmeans)

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Trial', 'Animal.ID', 'Date.Time', 
           'Schedule_Run_ID', 'Machine_Name', 'Drug')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Schedule_Length != "" & d$Trial == "None", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length < 3600, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Perform posthoc test for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Date.Time)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, method = "Wald")
summary(posthoc_results)
print(posthoc_cis)

# Create Kaplan-Meier plots for each genotype by day
ggsurvplot(survfit(surv_object ~ Genotype + Date.Time, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype and Date",
           facet.by = "Date.Time",
           facet.scales = "free_x")

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x")





# repeat for PR + ATO, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Schedule_Length != "" & d$Trial == "ATO", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Drug <- as.factor(cleaned_d$Drug)
cleaned_d$Drug <- factor(cleaned_d$Drug, levels = c("Saline", "ATO"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length < 3600, 1, 0)

# Create interaction terms for 'Genotype' and 'Drug'
cleaned_d$Genotype_Drug <- interaction(cleaned_d$Genotype, cleaned_d$Drug)
cleaned_d$Genotype_Drug <- as.factor(cleaned_d$Genotype_Drug)
cleaned_d$Genotype_Drug <- factor(cleaned_d$Genotype_Drug, levels = c('WT.Saline', 'WT.ATO', 'KI.Saline', 'KI.ATO'))

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Drug * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Perform posthoc test for interaction effect
cleaned_d$Date.Time <- as.factor(cleaned_d$Date.Time)
posthoc_model <- coxph(surv_object ~ Genotype * Drug * Date.Time + cluster(Animal.ID), data = cleaned_d)
emm_int <- emmeans(posthoc_model, pairwise ~ Genotype | Drug)
posthoc_results <- pairs(emm_int)
posthoc_cis <- confint(posthoc_model, method = "Wald")
summary(posthoc_results)
print(posthoc_cis)

# Create Kaplan-Meier plots for each genotype by day
ggsurvplot(survfit(surv_object ~ Genotype_Drug + Date.Time, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           palette = c("black", "red", "blue", "purple"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype and Date",
           facet.by = "Date.Time",
           facet.scales = "free_x")

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype_Drug, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "red", "blue", "purple"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x",
           legend.labs = c("WT.Saline", "WT.ATO", "KI.Saline", "KI.ATO"))







# repeat for PR + MPH, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_PR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Schedule_Length != "" & d$Trial == "MPH", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Drug <- as.factor(cleaned_d$Drug)
cleaned_d$Drug <- factor(cleaned_d$Drug, levels = c("Saline", "MPH"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length < 3600, 1, 0)

# Create interaction terms for 'Genotype' and 'Drug'
cleaned_d$Genotype_Drug <- interaction(cleaned_d$Genotype, cleaned_d$Drug)
cleaned_d$Genotype_Drug <- as.factor(cleaned_d$Genotype_Drug)
cleaned_d$Genotype_Drug <- factor(cleaned_d$Genotype_Drug, levels = c('WT.Saline', 'WT.MPH', 'KI.Saline', 'KI.MPH'))

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Drug * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Create Kaplan-Meier plots for each genotype by day
ggsurvplot(survfit(surv_object ~ Genotype_Drug + Date.Time, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           palette = c("black", "orange", "blue", "green"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype and Date",
           facet.by = "Date.Time",
           facet.by = "Date.Time",
           facet.scales = "free_x")

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype_Drug, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "orange", "blue", "green"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x",
           legend.labs = c("WT.Saline", "WT.MPH", "KI.Saline", "KI.MPH"))






# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "1", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x")








# repeat for FR2, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "2", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x")







# repeat for FR3, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "3", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x")







# repeat for FR5 (all sessions), read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "5", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x")







# repeat for FR5 before PR + ATO, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-02") & Date.Time <= as.Date("2018-06-09"))

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "5", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x")








# repeat for FR5 during PR + ATO, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-13") & Date.Time <= as.Date("2018-07-02"))

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length), ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x")








# repeat for FR5 during PR + MPH, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Julia 2018/Julia_FR and ERC.csv")
d <- d %>%
  filter(Date.Time >= as.Date("2018-06-13") & Date.Time <= as.Date("2018-07-02"))

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length), ]

# Convert columns 
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           xlim = c(0, 3600),
           break.x.by = 600,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
           facet.scales = "free_x")
