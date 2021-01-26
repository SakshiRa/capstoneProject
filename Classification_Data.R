getwd()
setwd("/Users/sakshirathi/ColorectalTumourClassification")

library(dplyr)
library(tibble)
library(caret)
library(caretEnsemble)
library(mice)
library(doParallel)
library(car)
library(lava)
library(tree)
library(rpart.plot)
library(pROC)
library(ROCR)
library(AppliedPredictiveModeling)

Data <- Data %>%
  select (-c(chr, pos))

Data$label <- as.factor(Data$label)
str(Data)

dim(Data)

table(Data$label)

options(digits = 2)

prop.table(table(Data$label))

summary <- do.call(cbind, lapply(Data[,1:ncol(Data)], summary))
summary <- as.data.frame((t(summary)))
summary

cooksd <- cooks.distance(glm(label ~ ., 
                             family = "binomial", 
                             data = Data))

plot(cooksd, 
     pch="*", 
     cex=2, 
     main="Influential Obs by Cooks distance")

abline(h = 4*mean(cooksd, na.rm=T), col="red")

outliers <- rownames(Data[cooksd > 4*mean(cooksd, na.rm=T), ])
outliers

descrCor <-  cor(Data[1:6])
descrCor

highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > .999)
highCorr

featurePlot(x = Data[,1:6], 
            y = Data$label, 
            plot = "box",
            scales = list(y = list(relation = "free"),
                          x = list(rot = 90)),
            layout = c(6,1),
            ## Add a key at the top
            auto.key = list(columns = 6))


set.seed(1000)
# Split dataset into 70% training, and 30% test
train_index <- createDataPartition(Data$label, ## outcome
                                   p = 0.7, ## percentage of training samples
                                   list = FALSE, ## show subsamples as matrix, not list
                                  
)

x.train <- Data[train_index,1:6] 
y.train <- Data$label[train_index]
x.test <- Data[-train_index,1:6]
y.test <- Data$label[-train_index]


train_control <- trainControl(
  method = "cv",
  number = 10, 
  summaryFunction=twoClassSummary, # computes area under the ROC curve
  classProbs = TRUE ## required for scoring models using ROC
)

set.seed(1000)
rf_train <- train(x = x.train, y = as.factor(y.train),
                   method='rf',
                   metric="ROC", ## default accuracy
                   trControl = train_control
                   )

rf_train
rf_train$resample
rf_test <- predict(rf_train, x.test) 

# compare predicted outcome and true outcome
conf_matrix <- confusionMatrix(rf_test, y.test, positive = 'true_variant')
conf_matrix$table
conf_matrix$byClass

accuracy <- sum(diag(conf_matrix$table))/sum(conf_matrix$table)
accuracy

## Can also spit out probability instead of predicted class
rf_test <- predict(rf_train, x.test, type = "prob")
rf_test
rf_test <- rf_test[,1]

############ Plot ROC curve ############
## ROC curve
rf <- roc(y.test,rf_test) ## pROC package
auc <- rf$auc
auc

plot(rf, col="blue",legacy.axes = TRUE)
impVars <- varImp(rf_train)
ImpMeasure <- data.frame(impVars$importance)
ImpMeasure <- ImpMeasure %>%
  arrange(desc(ImpMeasure))
round(ImpMeasure,2)

plot(varImp(rf_train))


