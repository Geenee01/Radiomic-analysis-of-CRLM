# mrmr feature selection for all aggregation files for overall survival

#Load necessary libraries
library(mRMRe)

### UWA ###
survdata <- readxl::read_xlsx("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/CRLM Clinical data.xlsx",col_names = TRUE) [,-c(1)]
surv_object <- Surv(event = survdata$vital_status, time = survdata$overall_survival_months)
mydata <- read.table("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/UWA_radiomic_aggregation.csv", header=TRUE, sep=",")[,-c(1)]
newdata <- cbind(surv_object, mydata)
#print(newdata)
dd <- mRMR.data(data = newdata)
#print(dd)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
scores(filter)
matrix_w =mim(filter) # to check selected features
#print( matrix_w)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter) # list of feature indices
solutions_indices <- unlist(solutions_indices) # converting list into vector
selected_features <- features[solutions_indices] # getting corresponding feature names
write.csv(selected_features, "C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/mrmr_overall_survival/UWA_10features_mRMR.csv", row.names = FALSE)


### WA ###
survdata <- readxl::read_xlsx("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/CRLM Clinical data.xlsx",col_names = TRUE)
surv_object <- Surv(event = survdata$vital_status, time = survdata$overall_survival_months)
mydata <- read.table("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/WA_radiomic_aggregation.csv", header=TRUE, sep=",")[,-c(1)]
newdata <- cbind(surv_object, mydata)
dd <- mRMR.data(data = newdata)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter)
solutions_indices <- unlist(solutions_indices)
selected_features <- features[solutions_indices]
write.csv(selected_features, "C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/mrmr_overall_survival/WA_10features_mRMR.csv", row.names = FALSE)

### WA of largest 3 ###
survdata <- readxl::read_xlsx("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/CRLM Clinical data.xlsx",col_names = TRUE)
surv_object <- Surv(event = survdata$vital_status, time = survdata$overall_survival_months)
mydata <- read.table("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/largest3_radiomic_aggregation.csv", header=TRUE, sep=",")[,-c(1)]
newdata <- cbind(surv_object, mydata)
dd <- mRMR.data(data = newdata)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter)
solutions_indices <- unlist(solutions_indices)
selected_features <- features[solutions_indices]
write.csv(selected_features, 'C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/mrmr_overall_survival/largest3_10features_mRMR.csv', row.names = FALSE)

### Largest 1 ###
survdata <- readxl::read_xlsx("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/CRLM Clinical data.xlsx",col_names = TRUE)
surv_object <- Surv(event = survdata$vital_status, time = survdata$overall_survival_months)
mydata <- read.table("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/largest1_radiomic.csv", header=TRUE, sep=",")[,-c(1)]
mydata <- mydata[,-ncol(mydata)]
newdata <- cbind(surv_object, mydata)
dd <- mRMR.data(data = newdata)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter)
solutions_indices <- unlist(solutions_indices)
selected_features <- features[solutions_indices]
write.csv(selected_features, 'C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/mrmr_overall_survival/largest1_10features_mRMR.csv', row.names = FALSE)

### Smallest 1 ###
survdata <- readxl::read_xlsx("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/CRLM Clinical data.xlsx",col_names = TRUE)
surv_object <- Surv(event = survdata$vital_status, time = survdata$overall_survival_months)
mydata <- read.table("C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/smallest1_radiomic.csv", header=TRUE, sep=",")[,-c(1)]
newdata <- cbind(surv_object, mydata)
dd <- mRMR.data(data = newdata)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter)
solutions_indices <- unlist(solutions_indices)
selected_features <- features[solutions_indices]
write.csv(selected_features, 'C:/Users/Hemangini/Desktop/MCRCdata/Analysis0-300/mrmr_overall_survival/smallest1_10features_mRMR.csv', row.names = FALSE)
