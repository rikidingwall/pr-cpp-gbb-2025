# Load the dependencies
library(survival)
library(survminer)
library(ggplot2)

# read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 
           'Schedule_Run_ID', 'Machine_Name')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Schedule_Length != "", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length < 3600, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
cox_model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(cox_model)

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






# repeat for PR only using days with 5 min timeout, read in PR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_PR.csv")

# Define the start and end dates
start_date <- as.Date("2017-11-30")
end_date <- as.Date("2017-12-01") 

# Clean up data file
d <- d[, c('Schedule_Length', 'Genotype', 'Animal.ID', 'Date.Time', 
           'Schedule_Run_ID', 'Machine_Name')]
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Schedule_Length != "" &
                 d$Date.Time >= start_date & d$Date.Time <= end_date, ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length < 3600, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
cox_model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(cox_model)

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







# repeat for FR1, read in FR data file 
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")

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
cox_model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(cox_model)

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







# repeat for FR2, read in FR data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")

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
cox_model <- coxph(surv_object ~ Genotype, data = cleaned_d)
summary(cox_model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
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






# repeat for FR3, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")

# Clean up data file
cleaned_d <- d[!is.na(d$Schedule_Length) & d$Ratio == "3", ]

# Convert columns
cleaned_d$Genotype <- as.factor(cleaned_d$Genotype)
cleaned_d$Genotype <- factor(cleaned_d$Genotype, levels = c("WT", "KI"))
cleaned_d$Date.Time <- as.Date(cleaned_d$Date.Time)
cleaned_d$Date.Time <- scale(cleaned_d$Date.Time)
cleaned_d$Event <- ifelse(cleaned_d$Schedule_Length > 0, 1, 0)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = cleaned_d$Schedule_Length, event = cleaned_d$Event)
cox_model <- coxph(surv_object ~ Genotype, data = cleaned_d)
summary(cox_model)

# Create Kaplan-Meier plots for each genotype
ggsurvplot(survfit(surv_object ~ Genotype, data = cleaned_d), 
           data = cleaned_d, 
           pval = FALSE, 
           conf.int = TRUE,
           palette = c("black", "blue"),
           xlab = "Schedule Length (s)",
           ylab = "Still in Task (%)",
           title = "Survival Curves by Genotype",
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







# repeat for FR5, read in data file
d = read.csv ("/Volumes/LaCie/Florey PhD/Original Papers/PR-ATO & MPH/Data/Spreadsheets/Touchscreens/Carlos 2017/Carlos_FR.csv")

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
cox_model <- coxph(surv_object ~ Genotype * Date.Time + cluster(Animal.ID), data = cleaned_d)
summary(cox_model)

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
