# mrmr feature selection 

#Load necessary libraries
library(mRMRe)

### UWA ###
survdata <- readxl::read_xlsx("path/to/clinical.xlsx",col_names = TRUE) [,-c(1)]
surv_object <- Surv(event = survdata$fu status, time = survdata$overall surv months)
mydata <- read.table("path/to/UWA_radiomic_aggregation.csv", header=TRUE, sep=",")[,-c(1)]
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
write.csv(selected_features, "path/to/UWA_10features_mRMR.csv", row.names = FALSE)


### WA ###
survdata <- readxl::read_xlsx("path/to/clinical.xlsx",col_names = TRUE)
surv_object <- Surv(event = survdata$fu status, time = survdata$overall surv months)
mydata <- read.table("path/to/WA_radiomic_aggregation.csv", header=TRUE, sep=",")[,-c(1)]
newdata <- cbind(surv_object, mydata)
dd <- mRMR.data(data = newdata)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter)
solutions_indices <- unlist(solutions_indices)
selected_features <- features[solutions_indices]
write.csv(selected_features, "path/to/WA_10features_mRMR.csv", row.names = FALSE)

### WA of largest 3 ###
survdata <- readxl::read_xlsx("path/to/clinical.xlsx",col_names = TRUE)
surv_object <- Surv(event = survdata$fu status, time = survdata$overall surv months)
mydata <- read.table("path/to/largest3_aggregation.csv", header=TRUE, sep=",")[,-c(1)]
newdata <- cbind(surv_object, mydata)
dd <- mRMR.data(data = newdata)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter)
solutions_indices <- unlist(solutions_indices)
selected_features <- features[solutions_indices]
write.csv(selected_features, 'path/to/largest3_10features_mRMR.csv', row.names = FALSE)

### Largest 1 ###
survdata <- readxl::read_xlsx("path/to/clinical.xlsx",col_names = TRUE)
surv_object <- Surv(event = survdata$fu status, time = survdata$overall surv months)
mydata <- read.table("path/to/largest1_radiomic.csv", header=TRUE, sep=",")[,-c(1)]
mydata <- mydata[,-ncol(mydata)]
newdata <- cbind(surv_object, mydata)
dd <- mRMR.data(data = newdata)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter)
solutions_indices <- unlist(solutions_indices)
selected_features <- features[solutions_indices]
write.csv(selected_features, 'path/to/largest1_10features_mRMR.csv', row.names = FALSE)

### Smallest 1 ###
survdata <- readxl::read_xlsx("path/to/clinical.xlsx",col_names = TRUE)
surv_object <- Surv(event = survdata$fu status, time = survdata$overall surv months)
mydata <- read.table("path/to/smallest1_radiomic.csv", header=TRUE, sep=",")[,-c(1)]
newdata <- cbind(surv_object, mydata)
dd <- mRMR.data(data = newdata)
filter <- mRMR.classic("mRMRe.Filter", data = dd, target_indices = c(1), feature_count = 10)
features <- featureNames(filter) # all feature with indices used in mrmr
solutions_indices <- solutions(filter)
solutions_indices <- unlist(solutions_indices)
selected_features <- features[solutions_indices]
write.csv(selected_features, 'path/to/smallest1_10features_mRMR.csv', row.names = FALSE)
