# Load the dependencies
library(survival)
library(survminer)
library(ggplot2)
library(emmeans)

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_PR.csv")

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 
           'Schedule_Run_ID', 'Machine_Name')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Schedule_Length != "", ]

# Create events and convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Animal.ID <- as.factor(cleaned_d$Animal.ID)
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length < 3600, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
cox_model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(cox_model)

cleaned_d <- cleaned_d[, c('Schedule_Length', 'Genotype', 'Event')]
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d %>% arrange(Genotype)
write.csv(cleaned_d, file = "/Users/rikidingwall/Desktop/Schedule Length Graphing/tom_PR.csv", row.names = FALSE)

# Create Kaplan-Meier plots for each genotype by day
ggsurvplot(survfit(surv_object ~ Genotype + Date.Time, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           palette = c("black", "grey"),
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







# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "1", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype + cluster(Animal.ID), data = cleaned_d)
summary(model)

cleaned_d <- cleaned_d[, c('Schedule_Length', 'Genotype', 'Event')]
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d %>% arrange(Genotype)
write.csv(cleaned_d, file = "/Users/rikidingwall/Desktop/Schedule Length Graphing/tom_FR1.csv", row.names = FALSE)

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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "2", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype + cluster(Animal.ID), data = cleaned_d)
summary(model)

cleaned_d <- cleaned_d[, c('Schedule_Length', 'Genotype', 'Event')]
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d %>% arrange(Genotype)
write.csv(cleaned_d, file = "/Users/rikidingwall/Desktop/Schedule Length Graphing/tom_FR2.csv", row.names = FALSE)

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








# repeat for FR3, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

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

cleaned_d <- cleaned_d[, c('Schedule_Length', 'Genotype', 'Event')]
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d %>% arrange(Genotype)
write.csv(cleaned_d, file = "/Users/rikidingwall/Desktop/Schedule Length Graphing/tom_FR3.csv", row.names = FALSE)

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







# repeat for FR4, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "4", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype + cluster(Animal.ID), data = cleaned_d)
summary(model)

cleaned_d <- cleaned_d[, c('Schedule_Length', 'Genotype', 'Event')]
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d %>% arrange(Genotype)
write.csv(cleaned_d, file = "/Users/rikidingwall/Desktop/Schedule Length Graphing/tom_FR4.csv", row.names = FALSE)

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
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

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

cleaned_d <- cleaned_d[, c('Schedule_Length', 'Genotype', 'Event')]
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d <- cleaned_d %>% arrange(Genotype)
write.csv(cleaned_d, file = "/Users/rikidingwall/Desktop/Schedule Length Graphing/tom_FR5.csv", row.names = FALSE)

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








# repeat for FR5 (before PR), read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Define the start and end dates
d$Date.Time <- as.Date(d$Date.Time)
start_date <- as.Date("2016-10-02")
end_date <- as.Date("2016-10-05")

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "5" & d$Date.Time >= start_date & d$Date.Time <= end_date, ]

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







# repeat for FR5 (after PR), read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Define the start and end dates
d$Date.Time <- as.Date(d$Date.Time)
start_date <- as.Date("2016-10-21")
end_date <- as.Date("2016-10-24")

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "5" & d$Date.Time >= start_date & d$Date.Time <= end_date, ]

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






# repeat for FR10, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Cleaned up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "10" & d$Schedule_Length > 50, ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0 & cleaned_d$Schedule_Length != 3600, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

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






# repeat for FR20, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "20" & d$Schedule_Length > 50, ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0 & cleaned_d$Schedule_Length != 3600, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

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






# repeat for FR40, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Tom 2016/Carlos_Tom_FR.csv")

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 'Ratio')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "40", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0 & cleaned_d$Schedule_Length != 3600, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(model)

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
