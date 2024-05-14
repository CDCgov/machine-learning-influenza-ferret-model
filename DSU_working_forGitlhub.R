########################################################################################
########################################################################################
### CDC Data Science Upskilling Project
### Influenza - Ferret Predictive Analytics
### Machine learning approaches for influenza A risk assessment: 
###   identifying predictive correlates using ferret model in vivo data
### 4 November 2022 - 5 October 2023
### Troy J. Kieran, Jessica Belser, Xiangjie Sun
########################################################################################
########################################################################################
### Load packages
library(tidyverse)
library(tidylog) ## detailed log of dplyr/tidyr functions
library(funModeling) ## loads Hmisc
library(DescTools) ## for AUC
library(caret) ## modeling
library(caretEnsemble) ## ensemble/stacked models
library(patchwork)

########################################################################################
### Import Data

## download data
# Pathogenesis Laboratory Team, Influenza Division, CDC. 
# An aggregated dataset of serially collected influenza A virus morbidity and 
# titer measurements from virus-infected ferrets.  
# https://data.cdc.gov/National-Center-for-Immunization-and-Respiratory-D/An-aggregated-dataset-of-serially-collected-influe/cr56-k9wj/about_data

## replace file.csv with name of file
# fullData <- read.csv("file.csv", header = TRUE)

########################################################################################

## if 'describe()' produces an nchar error, can run this
Sys.setlocale('LC_ALL', 'C')

########################################################################################

## calculate Area Under the Curve (AUC)
AUC_inoc <- 
  fullData %>% 
  gather(d1_inoc:d9_inoc, 
         key = day_inoc, 
         value = inoc_titer, 
         factor_key = TRUE) %>%
  #group_by(units) %>%
  mutate(day = case_when(day_inoc == "d1_inoc" ~ 1,
                         day_inoc == "d2_inoc" ~ 2,
                         day_inoc == "d3_inoc" ~ 3,
                         day_inoc == "d4_inoc" ~ 4,
                         day_inoc == "d5_inoc" ~ 5,
                         day_inoc == "d6_inoc" ~ 6,
                         day_inoc == "d7_inoc" ~ 7,
                         day_inoc == "d8_inoc" ~ 8,
                         day_inoc == "d9_inoc" ~ 9))

## filter AUC_inoc to each day inclusion threshold
AUC_4 <- AUC_inoc %>% filter(day %in% c(1, 2, 3, 4))
AUC_6 <- AUC_inoc %>% filter(day %in% c(1, 2, 3, 4, 5, 6))
AUC_8 <- AUC_inoc %>% filter(day %in% c(1, 2, 3, 4, 5, 6, 7, 8))
AUC_9 <- AUC_inoc %>% filter(day %in% c(1, 2, 3, 4, 5 ,6 ,7 ,8 ,9))

## Create a function to calculate trapezoid AUC for each Ferret
calculate_AUC_tr <- function(df) {
  data.frame(
    Ferret = unique(df$Ferret),
    trap_AUC = sapply(unique(df$Ferret), function(x) {
      auc <- AUC(df[df$Ferret == x, "day"], df[df$Ferret == x, "inoc_titer"], 
                 method = "trapezoid", na.rm = TRUE)
      ifelse(is.na(auc), NA_real_, auc)}))}

## Calculate trapezoid AUC for each data frame and combine with fullData
fullData <- fullData %>%
  left_join(calculate_AUC_tr(AUC_4), by = "Ferret") %>%
  rename(AUC_4 = trap_AUC) %>%
  left_join(calculate_AUC_tr(AUC_6), by = "Ferret") %>%
  rename(AUC_6 = trap_AUC) %>%
  left_join(calculate_AUC_tr(AUC_8), by = "Ferret") %>%
  rename(AUC_8 = trap_AUC) %>%
  left_join(calculate_AUC_tr(AUC_9), by = "Ferret") %>%
  rename(AUC_9 = trap_AUC)

## add slope for days 1-3, and peak_inoc
fullData <- fullData %>%
  mutate(fullData, slope13 = (d3_inoc - d1_inoc)/(24*2)) %>%
  mutate(peak_inoc = pmax(d1_inoc, d2_inoc, d3_inoc, d4_inoc,
                          d5_inoc, d6_inoc, d7_inoc, d8_inoc,
                          d9_inoc, na.rm = TRUE))

## rename Origin and add new column
## Origin to avian or mammal only
fullData <- 
  fullData %>%
  rename(Origin_orig = Origin) %>%
  mutate(Origin = if_else(Origin_orig == "avian" , "avian", "mammal"))

### add weight loss category 

## check initial bins
describe(equal_freq(fullData$wt_loss, n_bins = 3))
## remove wt_loss less than 5 and then
## examine the data into thirds for categories
fullData_wt_loss5 <- 
  fullData %>%
  dplyr::filter(wt_loss >= 5) %>% 
  dplyr::select(wt_loss) 
## use adjusted bins
describe(equal_freq(fullData_wt_loss5$wt_loss, n_bins = 3))

fullData <- 
  fullData %>%
  mutate(wt_loss_cat = case_when(wt_loss < 5 ~ 'none',
                                 wt_loss < 9.5 ~ 'low',
                                 wt_loss < 14.5 ~ 'med',
                                 wt_loss < 27.6 ~ 'high'))

fullData <-
  fullData %>%
  mutate(wt_loss_high = ifelse(wt_loss_cat == 'high', 'yes', 'no')) 

## check
str(fullData)
skimr::skim(fullData)

########################################################################################

## Basic exploratory data summaries - funModeling package
status(fullData)
data_integrity(fullData)
plot_num(fullData)
profiling_num(fullData)

########################################################################################

## Data imputation

## Convert NAs to median by units+Virus for wt_loss, temp
## Fill in the d2,4,6,8 inoc with the mean of day on either side
## Convert NAs to median by units+Virus for dx_inoc
## 
## Convert dx_inoc values back to NA for euth individuals
fullData_imputed <- 
  fullData %>%
  mutate(d2_inoc = ifelse(is.na(d2_inoc), 
                          rowMeans(.[, c("d1_inoc", "d3_inoc")]), 
                          d2_inoc),
         d4_inoc = ifelse(is.na(d4_inoc), 
                          rowMeans(.[, c("d3_inoc", "d5_inoc")]), 
                          d4_inoc),
         d6_inoc = ifelse(is.na(d6_inoc), 
                          rowMeans(.[, c("d5_inoc", "d7_inoc")]), 
                          d6_inoc),
         d8_inoc = ifelse(is.na(d8_inoc), 
                          rowMeans(.[, c("d7_inoc", "d9_inoc")]), 
                          d8_inoc),
         d3_inoc = ifelse(is.na(d3_inoc), 
                          rowMeans(.[, c("d2_inoc", "d4_inoc")]), 
                          d3_inoc),
         d5_inoc = ifelse(is.na(d5_inoc), 
                          rowMeans(.[, c("d4_inoc", "d6_inoc")]), 
                          d5_inoc),
         d7_inoc = ifelse(is.na(d7_inoc), 
                          rowMeans(.[, c("d6_inoc", "d8_inoc")]), 
                          d7_inoc)) %>%
  group_by(units, Virus) %>% 
  mutate(wt_loss = ifelse(is.na(wt_loss),
                          median(wt_loss, na.rm = TRUE),
                          wt_loss),
         temp = ifelse(is.na(temp),
                       median(temp, na.rm = TRUE),
                       temp),
         temp_5 = ifelse(is.na(temp_5),
                       median(temp_5, na.rm = TRUE),
                       temp_5),
         d1_inoc = ifelse(is.na(d1_inoc),
                          median(d1_inoc, na.rm = TRUE),
                          d1_inoc),
         d2_inoc = ifelse(is.na(d2_inoc),
                          median(d2_inoc, na.rm = TRUE),
                          d2_inoc),
         d3_inoc = ifelse(is.na(d3_inoc),
                          median(d3_inoc, na.rm = TRUE),
                          d3_inoc),
         d4_inoc = ifelse(is.na(d4_inoc),
                          median(d4_inoc, na.rm = TRUE),
                          d4_inoc),
         d5_inoc = ifelse(is.na(d5_inoc),
                          median(d5_inoc, na.rm = TRUE),
                          d5_inoc),
         d6_inoc = ifelse(is.na(d6_inoc),
                          median(d6_inoc, na.rm = TRUE),
                          d6_inoc),
         d7_inoc = ifelse(is.na(d7_inoc),
                          median(d7_inoc, na.rm = TRUE),
                          d7_inoc),
         d8_inoc = ifelse(is.na(d8_inoc),
                          median(d8_inoc, na.rm = TRUE),
                          d8_inoc),
         d9_inoc = ifelse(is.na(d9_inoc),
                          median(d9_inoc, na.rm = TRUE),
                          d9_inoc)) %>% 
  ungroup()

## check
status(fullData_imputed)
data_integrity(fullData_imputed)
plot_num(fullData_imputed)
profiling_num(fullData_imputed) #%>% View()
skimr::skim_without_charts(fullData_imputed)

## check numeric correlations
fullData_imputed %>% 
  dplyr::select_if(is.numeric) %>% 
  na.omit() %>% 
  cor() %>%
  findCorrelation(cutoff = 0.9, verbose = FALSE, names = TRUE)

## Summary info
fullData %>% filter(!is.na(expt)) %>% janitor::tabyl(HA)

########################################################################################

## predictive modeling

########################################################################################

## predictive modeling - lethality

## tested Origin_orig v Origin - 
## tested temp_5 v temp -
## Origin_orig and temp_5 provide overall better models
leth_test1 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin_orig,
                  wt_loss, temp_5, peak_inoc, AUC_6, slope13))

leth_test2 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin,
                  wt_loss, temp_5, peak_inoc, AUC_6, slope13))

leth_test3 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin_orig,
                  wt_loss, temp, peak_inoc, AUC_6, slope13))

## use Origin_orig + temp_5, swap peak_inoc w/ HA
leth_test4 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin_orig,
                  wt_loss, temp_5, HA, AUC_6, slope13))

## use Origin_orig + temp_5, remove peak_inoc + HA
leth_test5 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin_orig,
                  wt_loss, temp_5, AUC_6, slope13))

## use peak_inoc and HA
leth_test6 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin, HA,
                  wt_loss, temp_5, peak_inoc, AUC_6, slope13))

## test6 with Origin_orig replacing Origin
leth_test7 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin_orig, HA,
                  wt_loss, temp_5, peak_inoc, AUC_6, slope13))

## test6 removing slop13
leth_test8 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin, HA,
                  wt_loss, temp_5, peak_inoc, AUC_6))

## test8 removing peak_inoc
leth_test9 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin, HA,
                  wt_loss, temp_5, AUC_6))

## test9 - swap HA for Subtype
leth_test10 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, Origin, Subtype,
                  wt_loss, temp_5, AUC_6))

## test9 - remove Origin
leth_test11 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6))

###

## one hot code molecular data
leth_test1_dummy <- fastDummies::dummy_cols(
  leth_test1, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin_orig'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test2_dummy <- fastDummies::dummy_cols(
  leth_test2, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test3_dummy <- fastDummies::dummy_cols(
  leth_test3, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin_orig'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test4_dummy <- fastDummies::dummy_cols(
  leth_test4, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin_orig', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test5_dummy <- fastDummies::dummy_cols(
  leth_test5, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin_orig'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test6_dummy <- fastDummies::dummy_cols(
  leth_test6, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test7_dummy <- fastDummies::dummy_cols(
  leth_test7, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin_orig', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test8_dummy <- fastDummies::dummy_cols(
  leth_test8, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test9_dummy <- fastDummies::dummy_cols(
  leth_test9, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test10_dummy <- fastDummies::dummy_cols(
  leth_test10, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'Subtype'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test11_dummy <- fastDummies::dummy_cols(
  leth_test11, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

## set outcome variable to factor for classification models
leth_test1_dummy$lethal <- as.factor(leth_test1_dummy$lethal)
leth_test2_dummy$lethal <- as.factor(leth_test2_dummy$lethal)
leth_test3_dummy$lethal <- as.factor(leth_test3_dummy$lethal)
leth_test4_dummy$lethal <- as.factor(leth_test4_dummy$lethal)
leth_test5_dummy$lethal <- as.factor(leth_test5_dummy$lethal)
leth_test6_dummy$lethal <- as.factor(leth_test6_dummy$lethal)
leth_test7_dummy$lethal <- as.factor(leth_test7_dummy$lethal)
leth_test8_dummy$lethal <- as.factor(leth_test8_dummy$lethal)
leth_test9_dummy$lethal <- as.factor(leth_test9_dummy$lethal)
leth_test10_dummy$lethal <- as.factor(leth_test10_dummy$lethal)
leth_test11_dummy$lethal <- as.factor(leth_test11_dummy$lethal)

## split data into train/test using the tidymodels/rsample package
set.seed(9595)
## rotate through testing sets
#leth_test_dumSplit <- rsample::initial_split(leth_test1_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test2_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test3_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test4_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test5_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test6_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test7_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test8_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test9_dummy, prop = 0.70)
#leth_test_dumSplit <- rsample::initial_split(leth_test10_dummy, prop = 0.70)
leth_test_dumSplit <- rsample::initial_split(leth_test11_dummy, prop = 0.70)

trainData <- rsample::training(leth_test_dumSplit)
testData <- rsample::testing(leth_test_dumSplit)

###

## check backward feature selection - no need to use
rfeCtrl <- rfeControl(functions = rfFuncs,
                      method = "repeatedcv",
                      number = 10,
                      repeats = 2,
                      saveDetails = TRUE,
                      verbose = FALSE)

trainData2 <- na.omit(trainData)

rfeProfile <- rfe(x = trainData2[,2:17],
                 y = trainData2$lethal,
                 sizes = c(1:16),
                 rfeControl = rfeCtrl)
rfeProfile

###

## Predictive Power Scores
## compare AUC, AUC_6 is the best, supporting models
set.seed(123)
fullData %>%
  dplyr::select(lethal, AUC_4, AUC_6, AUC_8, AUC_9) %>%
  ppsr::score_predictors(df = ., y = 'lethal')

set.seed(123)
fullData %>%
  dplyr::select(lethal, AUC_4, AUC_6, AUC_8, AUC_9) %>%
  ppsr::visualize_pps(df = ., y = 'lethal')

## top model metrics - AUC_6 and wt_loss are the most predictive. 
set.seed(123)
fullData %>%
  dplyr::select(lethal, wt_loss, HPAI_MBAA, AUC_6, temp_5, 
                HA, PBS, RBS) %>%
  ppsr::score_predictors(df = ., y = 'lethal')

########################################################################################

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

## train models with split trainData
set.seed(2626)
m1 <- train(lethal ~ ., data = trainData,
            method = "glm",
            family = "binomial",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m2 <- train(lethal ~ ., data = trainData,
            method = "knn",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m3 <- train(lethal ~ ., data = trainData,
            method = "nnet",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m4 <- train(lethal ~ ., data = trainData,
            method = "glmnet",
            family = "binomial",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m5 <- train(lethal ~ ., data = trainData,
            method = "ranger",
            importance = 'impurity',
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m6 <- train(lethal ~ ., data = trainData,
            method = "svmRadial", 
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m7 <- train(lethal ~ ., data = trainData,
            method = "bayesglm",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m8 <- train(lethal ~ ., data = trainData,
            method = "LogitBoost",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m9 <- train(lethal ~ ., data = trainData,
            method = "rf",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m10 <- train(lethal ~ ., data = trainData,
             method = "rpart",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m11 <- train(lethal ~ ., data = trainData,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

# set.seed(2626)
# m11.1 <- train(lethal ~ ., data = trainData,
#              method = "xgbLinear",
#              na.action = na.exclude,
#              preProcess = c("nzv", "scale", "center"),
#              trControl = fitControl)

## resample metrics
resamps <- resamples(list(glm = m1,
                          knn = m2,
                          nnet = m3,
                          glmnet = m4,
                          ran = m5,
                          svm = m6,
                          bglm = m7,
                          BlogReg = m8,
                          rf = m9,
                          rpart = m10,
                          gbm = m11))#,
                          #xgb = m11.1))
summary(resamps)
bwplot(resamps)
dotplot(resamps)
modelCor(resamps)
splom(resamps)

## difference between models
difValues <- diff(resamps)
difValues
summary(difValues)
bwplot(difValues)
dotplot(difValues)

## check for variable importance
ggplot(varImp(m1))
ggplot(varImp(m3))
ggplot(varImp(m4))
ggplot(varImp(m5))
ggplot(varImp(m9))
ggplot(varImp(m10))
## errors:
# ggplot(varImp(m2))
# ggplot(varImp(m6))
# ggplot(varImp(m7))
# ggplot(varImp(m8))
varImp(m1, scale = TRUE)
varImp(m3, scale = TRUE)
varImp(m4, scale = TRUE)
varImp(m5, scale = TRUE)
varImp(m9, scale = TRUE)
varImp(m10, scale = TRUE)
varImp(m11, scale = TRUE) # library(gbm)

## more control over plotting variable importance
varImp(m3, scale = TRUE)$importance %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(desc(Overall)) %>%
  mutate(rowname = forcats::fct_inorder(rowname)) %>%
  head(40)

###

## check probabilities of call
leth_Models <- list(glm = m1,
                    knn = m2,
                    glm = m3,
                    glm = m4,
                    ranger = m5,
                    bglm = m6,
                    boostLogReg = m7,
                    rf = m8,
                    rpart = m9,
                    gbm = m10,
                    xgb = m11.1)
extractPrediction(leth_Models, testX = testData)
extractProb(leth_Models, testX = testData) %>%
  as.data.frame() %>%
  gather(no:yes, key = lethal, value = prob, factor_key = TRUE) %>%
  ggplot(aes(x = lethal, y = prob, fill = model)) +
  geom_jitter(alpha = 0.1) +
  geom_violin(alpha = 0.7) +
  facet_grid(~ model)

extractProb(leth_Models, testX = testData) %>%
  filter(obs == 'yes', pred == 'no') %>% view()

extractProb(leth_Models, testX = testData) %>%
  filter(obs == 'no', pred == 'yes') %>% view()

thresholder(m5, threshold = seq(.01, 0.99, by = 0.01), 
            final = TRUE, statistics = "all") %>% View()


## check which parameters can be changed
modelLookup(model = 'gbm') ## n.trees, interaction.depth, shrinkage, n.minobsinnode

m11Gridl <- expand.grid(n.trees = seq(from = 50, to = 500, by = 50), 
                        interaction.depth = seq(from = 2, to = 4, by = 1),
                        shrinkage = seq(from = 0.1, to = 0.3, by = 0.1),
                        n.minobsinnode = seq(from = 5, to = 10, by = 1))

m11Grid_final <- expand.grid(n.trees = 250, 
                        interaction.depth = 4,
                        shrinkage = 0.2,
                        n.minobsinnode = 8)

set.seed(2626)
m11_tune <- train(lethal ~ ., data = trainData,
             method = "gbm",
             na.action = na.exclude,
             #tuneLength = 10,
             tuneGrid = m11Grid_final,
             metric = "Balanced_Accuracy",
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

## save the model
saveRDS(m11_tune, "./final_models/m11_tune_lethality_final_model.rds")
## read model back in
m11_tune <- readRDS("./final_models/m11_tune_lethality_final_model.rds")


probTest <- predict(m11_tune, testData, type = 'prob')
#probTest <- probTest$lethal
probTest <- factor(ifelse(probTest$no >= 0.97, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$lethal, mode = "everything")

###

## predict testData
TYPE <- 'raw' # yes/no outcome
TYPE <- 'prob' # probability of yes/no outcome
set.seed(2626)
predictions <- data.frame(
  glm = predict(m1, testData, type = TYPE),
  knn = predict(m2, testData, type = TYPE),
  nnet = predict(m3, testData, type = TYPE),
  glmnet = predict(m4, testData, type = TYPE),
  ranger = predict(m5, testData, type = TYPE),
  svm = predict(m6, testData, type = TYPE),
  bglm = predict(m7, testData, type = TYPE),
  boostLogReg = predict(m8, testData, type = TYPE),
  randomForest = predict(m9, testData, type = TYPE),
  rpart = predict(m10, testData, type = TYPE),
  gbm = predict(m11, testData, type = TYPE),
  xgb = predict(m11.1, testData, type = TYPE)
)
predictions

## compile predictions with lethal truth
truth_predict <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., predictions) 

## short hand objects
tp <- truth_predict
tpl <- truth_predict$lethal
## confusion matrices with model metrics
confusionMatrix(tp$glm, tpl, mode = "everything") 
confusionMatrix(tp$knn, tpl, mode = "everything")
confusionMatrix(tp$nnet, tpl, mode = "everything") 
confusionMatrix(tp$glmnet, tpl, mode = "everything") 
confusionMatrix(tp$ranger, tpl, mode = "everything")
confusionMatrix(tp$svm, tpl, mode = "everything") 
confusionMatrix(tp$bglm, tpl, mode = "everything")
confusionMatrix(tp$boostLogReg, tpl, mode = "everything") 
confusionMatrix(tp$randomForest, tpl, mode = "everything") 
confusionMatrix(tp$rpart, tpl, mode = "everything")
confusionMatrix(tp$gbm, tpl, mode = "everything")
confusionMatrix(tp$xgb, tpl, mode = "everything")

compare_models(m3, m5)
compare_models(m5, m9)
compare_models(m3, m9)

########################################################################################

## hyper-parameter tuning - top 3 performing models

## check which parameters can be changed
modelLookup(model = 'knn') ## k
## expand parameter grid to test
knnGrid <-  expand.grid(k = seq(from = 1, to = 6, by = 0.5))
knnGrid2 <-  expand.grid(k = seq(from = 5, to = 6, by = 0.1))
knnGrid_final <- expand.grid(k = 5) 

set.seed(2626)
m2.1 <- train(lethal ~ ., data = trainData,
            method = "knn",
            na.action = na.exclude,
            metric = 'Balanced_Accuracy', 
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl,
            tuneGrid = knnGrid_final)
m2.1

###

## check which parameters can be changed
modelLookup(model = 'nnet') ## size, decay
## expand paramter grid to test
nnetGrid <-  expand.grid(size = seq(from = 3, to = 7, by = 0.5), 
                        decay = seq(from = 0, to = 0.5, by = 0.05))
nnetGrid2 <-  expand.grid(size = seq(from = 3, to = 5, by = 0.5), 
                         decay = seq(from = 0, to = 0.7, by = 0.05))
nnetGrid3 <-  expand.grid(size = seq(from = 3, to = 4, by = 0.1), 
                          decay = seq(from = 0.35, to = 0.55, by = 0.05))
nnetGrid_final <-  expand.grid(size = 3, decay = 0.4)

set.seed(2626)
m3.1 <- train(lethal ~ ., data = trainData,
            method = "nnet",
            na.action = na.exclude,
            metric = 'Balanced_Accuracy',
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl,
            tuneGrid = nnetGrid_final)
m3.1

###

## check which parameters can be changed
modelLookup(model = 'ranger') ## mtry, splitrule, min.node.size
## expand parameters grid to test
rangerGrid <-  expand.grid(mtry = seq(from = 5, to = 18, by = 1), 
                        splitrule = c('gini', 'extratrees'), 
                        min.node.size = seq(from = 1, to = 1, by = 1))
rangerGrid_final <-  expand.grid(mtry = 7, 
                           splitrule = 'extratrees', 
                           min.node.size = 1)

set.seed(2626)
m5.1 <- train(lethal ~ ., data = trainData,
            method = "ranger",
            na.action = na.exclude,
            metric = 'Balanced_Accuracy',
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)
m5.1

###

## check which parameters can be changed
modelLookup(model = 'rpart') ## cp
## expand paramter grid to test
rpartGrid <-  expand.grid(cp = seq(from = 0, to = 1, by = 0.01))
rpartGrid2 <-  expand.grid(cp = seq(from = 0, to = 0.1, by = 0.001))

set.seed(2626)
m10.1 <- train(lethal ~ ., data = trainData,
             method = "rpart",
             na.action = na.exclude,
             metric = 'Balanced_Accuracy',
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl,
             tuneGrid = rpartGrid2)
m10.1

###

## check which parameters can be changed
modelLookup(model = 'gbm') ## n.trees, interaction.depth, shrinkage, n.minobsinnode
## expand parameters grid to test
gbmGrid <-  expand.grid(n.trees = seq(from = 90, to = 110, by = 1), 
                        interaction.depth = seq(from = 1, to = 5, by = 1),
                        shrinkage = seq(from = 0, to = 0.5, by = 0.1),
                        n.minobsinnode = 10)
gbmGrid2 <-  expand.grid(n.trees = seq(from = 99, to = 110, by = 1), 
                        interaction.depth = seq(from = 3, to = 4, by = 1),
                        shrinkage = seq(from = 0.1, to = 0.5, by = 0.1),
                        n.minobsinnode = seq(from = 8, to = 12, by = 1))
gbmGrid3 <-  expand.grid(n.trees = seq(from = 100, to = 110, by = 1), 
                         interaction.depth = seq(from = 3, to = 4, by = 1),
                         shrinkage = seq(from = 0.1, to = 0.4, by = 0.1),
                         n.minobsinnode = seq(from = 9, to = 10, by = 1))
gbmGrid_final <-  expand.grid(n.trees = 102, 
                         interaction.depth = 3,
                         shrinkage = 0.3,
                         n.minobsinnode = 9)

set.seed(2626)
m11.2 <- train(lethal ~ ., data = trainData,
              method = "gbm",
              na.action = na.exclude,
              metric = 'Balanced_Accuracy',
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl,
              tuneGrid = gbmGrid_final)
m11.2

###

resamps_tune <- resamples(list(
                          knn = m2.1,
                          nnet = m3.1,
                          ran = m5.1,
                          rpart = m10.1,
                          gbm = m11.2))
resamps_tune
summary(resamps_tune)
bwplot(resamps_tune)
dotplot(resamps_tune)
modelCor(resamps_tune)

###

set.seed(999)
knnPred <- predict(m2.1, testData)

truth_predictKNN <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., knnPred) 

## short hand objects
tpKNN <- truth_predictKNN
tplKNN <- truth_predictKNN$lethal

confusionMatrix(tpKNN$knn, tplKNN, mode = "everything")

###

set.seed(999)
nnetPred <- predict(m3.1, testData)

truth_predictNNET <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., nnetPred) 

## short hand objects
tpNNET <- truth_predictNNET
tplNNET <- truth_predictNNET$lethal

confusionMatrix(tpNNET$nnet, tplNNET, mode = "everything")

###

set.seed(999)
rangerPred <- predict(m5.1, testData)

truth_predictRAN <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., rangerPred) 

## short hand objects
tpRAN <- truth_predictRAN
tplRAN <- truth_predictRAN$lethal

confusionMatrix(tpRAN$ranger, tplRAN, mode = "everything")

###

set.seed(999)
gbmPred <- predict(m11.2, testData)

truth_predictGBM <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., gbmPred) 

## short hand objects
tpGBM <- truth_predictGBM
tplGBM <- truth_predictGBM$lethal

confusionMatrix(tpGBM$gbm, tplGBM, mode = "everything")

###

#Make a list of all the models
#all.models <- list(m2.1, m3.1, m5.1, m10.1, m11.2)
#all.models <- list(m2.1, m3.1, m11.2)
#all.models <- list(m2.1, m3.1)
#all.models <- list(m2, m3.1, m11)
all.models <- list(m3.1, m5, m11)
names(all.models) <- sapply(all.models, function(x) x$method)
sort(sapply(all.models, function(x) min(x$results$Balanced_Accuracy)))
sort(sapply(all.models, function(x) min(x$results$Specificity)))
sort(sapply(all.models, function(x) min(x$results$Sensitivity)))
sort(sapply(all.models, function(x) min(x$results$Accuracy)))
class(all.models) <- "caretList"

#Make a greedy ensemble - currently can only use RMSE
greedy <- caretEnsemble(all.models)
greedy$error

#Make a linear regression ensemble
set.seed(2626)
nnet_stack <- caretStack(all.models, 
                         method = 'nnet', 
                         metric = "Balanced_Accuracy",
                         #tuneLength = 10, 
                         tuneGrid = expand.grid(size = 9, decay = 0.0001),
                         trControl = trainControl(method = "repeatedcv",   
                                                  number = 10,
                                                  repeats = 2,
                                                  savePredictions = 'final',
                                                  classProbs = TRUE,
                                                  sampling = 'smote',
                                                  summaryFunction = multiClassSummary))
nnet_stack

nnet_stack$error %>% 
  #filter(size == 15) %>%
  View()

## stacked predictions
stack_nnet_predictions <- data.frame(
  nnet_stack = predict(nnet_stack, testData))

## compile predictions with lethal truth
truth_predict_nnet <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., stack_nnet_predictions)

confusionMatrix(truth_predict_nnet$nnet_stack, 
                truth_predict_nnet$lethal, 
                mode = "everything") 


###
set.seed(2626)
probTest <- predict(nnet_stack, testData, type = 'prob')
#probTest <- probTest$lethal
probTest <- factor(ifelse(probTest >= 0.92, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$lethal, mode = "everything")


########################################################################################



########################################################################################

## Ensemble/Stacked lethal models

## leth_test8_dummy
#algoList <- c('knn', 'nnet', 'ranger') 

## leth_test9_dummy
algoList <- c('nnet', 'ranger', 'gbm') 

## leth_test11_dummy
#algoList <- c('knn', 'ranger')

set.seed(12345)
models1 <- caretList(lethal ~ ., 
                     data = trainData, 
                     metric = "Balanced_Accuracy",
                     preProcess = c("nzv", "scale", "center"),
                     trControl = fitControl, 
                     methodList = algoList)

models1
summary(resamples(models1))
modelCor(resamples(models1))

# stack using random forest as the second order model
set.seed(12345)
rand_forest.stack <- caretStack(models1, 
                                method = "rf", 
                                metric = "Balanced_Accuracy", 
                                #tuneLength = 10,
                                trControl = fitControl) ## no change using bootstraps

## stacked predictions
stack_predictions <- data.frame(
  rf_stack = predict(rand_forest.stack, testData))

## compile predictions with lethal truth
truth_predict <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., stack_predictions)

confusionMatrix(truth_predict$rf_stack, truth_predict$lethal, mode = "everything") 

###

probTest <- predict(rand_forest.stack, testData, type = 'prob')
#probTest <- probTest$lethal
probTest <- factor(ifelse(probTest >= 0.9, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$lethal, mode = "everything")

###

#Inspect elements of model
rand_forest.stack$models
rand_forest.stack$ens_model$finalModel$confusion
rand_forest.stack$ens_model$results$mtry
rand_forest.stack$ens_model$results$Accuracy
rand_forest.stack$ens_model$bestTune$mtry
rand_forest.stack$ens_model$resample
rand_forest.stack$ens_model$times$everything # Performance information re the random forest model utilized
rand_forest.stack$ens_model$pred # Returns top level predictions from the random forest 
rand_forest.stack$ens_model$perfNames
rand_forest.stack$ens_model$control
rand_forest.stack$ens_model$preProcess
rand_forest.stack$ens_model$maximize

###

## tuning parameters

set.seed(12345)
models2 <- caretList(lethal ~ ., 
                     data = trainData, 
                     metric = "Balanced_Accuracy",
                     preProcess = c("nzv", "scale", "center"),
                     trControl = fitControl, 
                     #methodList = algoList,
                     tuneList = list(
                       ran = caretModelSpec(method = "ranger", 
                                            tuneLength = 10),
                       gbm = caretModelSpec(method = "gbm", 
                                            tuneLength = 10),
                       nn = caretModelSpec(method = "nnet", 
                                           tuneLength = 2, 
                                           trace = FALSE)))
models2 
summary(resamples(models2))
modelCor(resamples(models2))

## stack using random forest as the second order model
set.seed(12345)
rand_forest.stack2 <- caretStack(models2, 
                                method = "rf",
                                metric = "Balanced_Accuracy", 
                                tuneLength = 10,
                                trControl = fitControl)

## stacked predictions 
stack_predictions2 <- data.frame(
  rf_stack = predict(rand_forest.stack2, testData))

## compile predictions with lethal truth
truth_predict2 <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., stack_predictions2)

confusionMatrix(truth_predict2$rf_stack, truth_predict2$lethal, mode = "everything") 

########################################################################################




########################################################################################

## lethality - molecular correlates models

leth_mol_test1 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, '17':PBS)) %>%
  dplyr::select(-'329')

leth_mol_test2 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, HA, RBS, PBS))

leth_mol_test3 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, Origin, HA,
                  wt_loss, temp_5, AUC_6, '17':PBS))

leth_mol_test4 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, RBS, PBS))

leth_mol_test5 <- fullData_imputed %>% 
  drop_na(lethal) %>%
  dplyr::select(c(lethal, HPAI_MBAA, HA, Origin, '17':PBS)) %>%
  dplyr::select(-'329')


## one hot code molecular data
leth_mol_test1_dummy <- fastDummies::dummy_cols(
  leth_mol_test1, select_columns =
    c('17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_mol_test2_dummy <- fastDummies::dummy_cols(
  leth_mol_test2, select_columns =
    c('HPAI_MBAA',	'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_mol_test3_dummy <- fastDummies::dummy_cols(
  leth_mol_test3, select_columns =
    c('HPAI_MBAA',	'HA', 'Origin', 
      '17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_mol_test4_dummy <- fastDummies::dummy_cols(
  leth_mol_test4, select_columns =
    c('RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_mol_test5_dummy <- fastDummies::dummy_cols(
  leth_mol_test5, select_columns =
    c('HPAI_MBAA',	'HA', 'Origin',
      '17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_mol_test1_dummy$lethal <- as.factor(leth_mol_test1_dummy$lethal)
leth_mol_test2_dummy$lethal <- as.factor(leth_mol_test2_dummy$lethal)
leth_mol_test3_dummy$lethal <- as.factor(leth_mol_test3_dummy$lethal)
leth_mol_test4_dummy$lethal <- as.factor(leth_mol_test4_dummy$lethal)
leth_mol_test5_dummy$lethal <- as.factor(leth_mol_test5_dummy$lethal)

###

## split data into train/test using the tidymodels/rsample package

set.seed(9595)
## rotate through testing sets
#leth_mol_test_dumSplit <- rsample::initial_split(leth_mol_test1_dummy, prop = 0.70)
#leth_mol_test_dumSplit <- rsample::initial_split(leth_mol_test2_dummy, prop = 0.70)
leth_mol_test_dumSplit <- rsample::initial_split(leth_mol_test3_dummy, prop = 0.70)
#leth_mol_test_dumSplit <- rsample::initial_split(leth_mol_test4_dummy, prop = 0.70)
#leth_mol_test_dumSplit <- rsample::initial_split(leth_mol_test5_dummy, prop = 0.70)

trainData_mol <- rsample::training(leth_mol_test_dumSplit)
testData_mol <- rsample::testing(leth_mol_test_dumSplit)

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

## train models with split trainData
set.seed(2626)
m1m <- train(lethal ~ ., data = trainData_mol,
            method = "glm",
            family = "binomial",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m2m <- train(lethal ~ ., data = trainData_mol,
            method = "knn",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m3m <- train(lethal ~ ., data = trainData_mol,
            method = "nnet",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m4m <- train(lethal ~ ., data = trainData_mol,
            method = "glmnet",
            family = "binomial",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m5m <- train(lethal ~ ., data = trainData_mol,
            method = "ranger",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m6m <- train(lethal ~ ., data = trainData_mol,
            method = "svmRadial",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m7m <- train(lethal ~ ., data = trainData_mol,
            method = "bayesglm",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m8m <- train(lethal ~ ., data = trainData_mol,
            method = "LogitBoost",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m9m <- train(lethal ~ ., data = trainData_mol,
            method = "rf",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

# set.seed(2626)
# m10m <- train(lethal ~ ., data = trainData_mol,
#              method = "rpart",
#              na.action = na.exclude,
#              preProcess = c("nzv", "scale", "center"),
#              trControl = fitControl)

set.seed(2626)
m11m <- train(lethal ~ ., data = trainData_mol,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

## resample metrics
resamps_lm <- resamples(list(glm = m1m,
                          knn = m2m,
                          nnet = m3m,
                          glmnet = m4m,
                          ran = m5m,
                          svm = m6m,
                          bglm = m7m,
                          BlogReg = m8m,
                          rf = m9m,
                          #rpart = m10m,
                          gbm = m11m))
summary(resamps_lm)
bwplot(resamps_lm)
dotplot(resamps_lm)
modelCor(resamps_lm)
splom(resamps_lm)

varImp(m11m, scale = TRUE)

## predict testData
set.seed(2626)
predictions_lm <- data.frame(
  glm = predict(m1m, testData_mol),
  knn = predict(m2m, testData_mol),
  nnet = predict(m3m, testData_mol),
  glmnet = predict(m4m, testData_mol),
  ranger = predict(m5m, testData_mol),
  svm = predict(m6m, testData_mol),
  bglm = predict(m7m, testData_mol),
  boostLogReg = predict(m8m, testData_mol),
  randomForest = predict(m9m, testData_mol),
  #rpart = predict(m10m, testData_mol),
  gbm = predict(m11m, testData_mol)
)

## compile predictions with lethal truth
truth_predict_lm <- testData_mol %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., predictions_lm) 

## short hand objects
tp_lm <- truth_predict_lm
tplm <- truth_predict_lm$lethal
## confusion matrices with model metrics
confusionMatrix(tp_lm$glm, tplm, mode = "everything") 
confusionMatrix(tp_lm$knn, tplm, mode = "everything")
confusionMatrix(tp_lm$nnet, tplm, mode = "everything") 
confusionMatrix(tp_lm$glmnet, tplm, mode = "everything") 
confusionMatrix(tp_lm$ranger, tplm, mode = "everything")
confusionMatrix(tp_lm$svm, tplm, mode = "everything") 
confusionMatrix(tp_lm$bglm, tplm, mode = "everything")
confusionMatrix(tp_lm$boostLogReg, tplm, mode = "everything") 
confusionMatrix(tp_lm$randomForest, tplm, mode = "everything") 
confusionMatrix(tp_lm$gbm, tplm, mode = "everything")

###



m11m_gbmGrid <- expand.grid(n.trees = seq(from = 450, to = 450, by = 50), 
                            interaction.depth = seq(from = 10, to = 10, by = 1),
                            shrinkage = seq(from = 0.1, to = 0.3, by = 0.1),
                            n.minobsinnode = seq(from = 8, to = 10, by = 1))

m11m_gbmGrid_final <-
  expand.grid(n.trees = 450,
              interaction.depth = 9,
              shrinkage = 0.1,
              n.minobsinnode = 10)

set.seed(2626)
m11m_tune <- train(lethal ~ ., data = trainData_mol,
                  method = "gbm",
                  na.action = na.exclude,
                  #tuneLength = 10,
                  tuneGrid = m11m_gbmGrid_final,
                  metric = "Balanced_Accuracy",
                  preProcess = c("nzv", "scale", "center"),
                  trControl = fitControl)

## test3 with molecular + other variables
## save the model
saveRDS(m11m_tune, "./final_models/m11m_tune_lethality_molecular_all_final_model.rds")
## read model back in
m11m_tune <- readRDS("./final_models/m11m_tune_lethality_molecular_all_final_model.rds")

## test1 with molecular only
## save the model
saveRDS(m11m_tune, "./final_models/m11m_tune_lethality_molecular_final_model.rds")
## read model back in
m11m_tune <- readRDS("./final_models/m11m_tune_lethality_molecular_final_model.rds")


## calculate 95% confidence intervals using z score 1.96 and n= number of crossfolds 
m11m_tune$results$Balanced_Accuracy+c(-1,1)*1.96*(m11m_tune$results$Balanced_AccuracySD/sqrt(20))
m11m_tune$results$Sensitivity+c(-1,1)*1.96*(m11m_tune$results$SensitivitySD/sqrt(20))
m11m_tune$results$Specificity+c(-1,1)*1.96*(m11m_tune$results$SpecificitySD/sqrt(20))
m11m_tune$results$Kappa+c(-1,1)*1.96*(m11m_tune$results$KappaSD/sqrt(20))

###


###

probTest <- predict(m11m_tune, testData_mol, type = 'prob')
#probTest <- probTest$lethal
probTest <- factor(ifelse(probTest$no >= 0.8, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData_mol %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$lethal, mode = "everything")

###

## mol_test1 - rf stack
#algoList <- c('glm', 'rf', 'gbm') 
#algoList <- c('glm', 'nnet', 'gbm')
#algoList <- c('glm', 'gbm')
algoList <- c('glm', 'glmnet', 'gbm')

## test 3 - rf stack
#algoList <- c('glm', 'glmnet', 'gbm')  
algoList <- c('glm', 'rf', 'gbm') 
#algoList <- c('glm', 'rf', 'gbm', 'bayesglm') 
## glm stack - rf wouldn't work
#algoList <- c('glm', 'rf', 'gbm', 'nnet')

set.seed(12345)
models_lm1 <- caretList(lethal ~., 
                     data = trainData_mol, 
                     metric = "Balanced_Accuracy",
                     preProcess = c("nzv", "scale", "center"),
                     trControl = fitControl,
                     tuneLength = 10,
                     methodList = algoList)
models_lm1
summary(resamples(models_lm1))
modelCor(resamples(models_lm1))

# stack using random forest as the second order model
set.seed(12345)
rand_forest.stack_lm <- caretStack(models_lm1, 
                                method = "glm", 
                                metric = "Balanced_Accuracy", 
                                tuneLength = 10,
                                trControl = fitControl) ## no change using bootstraps

## stacked predictions
stack_predictions_lm <- data.frame(
  rf_stack = predict(rand_forest.stack_lm, testData_mol, type = 'prob'))

## compile predictions with lethal truth
truth_predict_lm <- testData_mol %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., stack_predictions_lm)

confusionMatrix(truth_predict_lm$rf_stack, truth_predict_lm$lethal, mode = "everything") 

########################################################################################





########################################################################################

## predictive modeling - weight loss

wt_test1 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, lethal, HPAI_MBAA, RBS, PBS, 
                  Origin, HA, temp_5, AUC_6))

## remove lethal as non-prior info
wt_test2 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, 
                  Origin, HA, temp_5, AUC_6))

## swap AUC_6 with peak_inoc - AUC_6 performs better
wt_test3 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, 
                  Origin, HA, temp_5, peak_inoc))

## test2 w/ d3_inoc added - performs worse then test2
wt_test4 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, 
                  Origin, HA, temp_5, AUC_6, d3_inoc))

## swap AUC_6 for 9 - AUC_6 performs better
wt_test5 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, 
                  Origin, HA, temp_5, AUC_9))

## use only top three important vars across multiple models
wt_test6 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, temp_5, AUC_6))

## test2 without Origin - improves slightly
wt_test7 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, 
                  HA, temp_5, AUC_6))

## test7 without RBS
wt_test8 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, PBS, 
                  HA, temp_5, AUC_6))

## test7 without PBS
wt_test9 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, 
                  HA, temp_5, AUC_6))

## test7 without RBS & PBS
wt_test10 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, 
                  HA, temp_5, AUC_6))

###

## one hot code molecular data
wt_test1_dummy <- fastDummies::dummy_cols(
  wt_test1, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'HA', 'lethal'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test2_dummy <- fastDummies::dummy_cols(
  wt_test2, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test3_dummy <- fastDummies::dummy_cols(
  wt_test3, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test4_dummy <- fastDummies::dummy_cols(
  wt_test4, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test5_dummy <- fastDummies::dummy_cols(
  wt_test5, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'Origin', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test6_dummy <- fastDummies::dummy_cols(
  wt_test6, select_columns =
    c('HPAI_MBAA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test7_dummy <- fastDummies::dummy_cols(
  wt_test7, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test8_dummy <- fastDummies::dummy_cols(
  wt_test8, select_columns =
    c('HPAI_MBAA', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test9_dummy <- fastDummies::dummy_cols(
  wt_test9, select_columns =
    c('HPAI_MBAA', 'RBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test10_dummy <- fastDummies::dummy_cols(
  wt_test10, select_columns =
    c('HPAI_MBAA', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test1_dummy$wt_loss_high <- as.factor(wt_test1_dummy$wt_loss_high)
wt_test2_dummy$wt_loss_high <- as.factor(wt_test2_dummy$wt_loss_high)
wt_test3_dummy$wt_loss_high <- as.factor(wt_test3_dummy$wt_loss_high)
wt_test4_dummy$wt_loss_high <- as.factor(wt_test4_dummy$wt_loss_high)
wt_test5_dummy$wt_loss_high <- as.factor(wt_test5_dummy$wt_loss_high)
wt_test6_dummy$wt_loss_high <- as.factor(wt_test6_dummy$wt_loss_high)
wt_test7_dummy$wt_loss_high <- as.factor(wt_test7_dummy$wt_loss_high)
wt_test8_dummy$wt_loss_high <- as.factor(wt_test8_dummy$wt_loss_high)
wt_test9_dummy$wt_loss_high <- as.factor(wt_test9_dummy$wt_loss_high)
wt_test10_dummy$wt_loss_high <- as.factor(wt_test10_dummy$wt_loss_high)

###
## Predictive Power Scores

set.seed(123)
fullData %>%
  dplyr::select(wt_loss_high, AUC_4, AUC_6, AUC_8, AUC_9) %>%
  ppsr::score_predictors(df = ., y = 'wt_loss_high')

set.seed(123)
fullData %>%
  dplyr::select(wt_loss_high, AUC_4, AUC_6, AUC_8, AUC_9) %>%
  ppsr::visualize_pps(df = ., y = 'wt_loss_high')


set.seed(123)
fullData %>%
  dplyr::select(wt_loss_high, HPAI_MBAA, RBS, PBS, 
                HA, temp_5, AUC_6) %>%
  ppsr::score_predictors(df = ., y = 'wt_loss_high')



###

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

## split data into train/test using the tidymodels/rsample package
set.seed(9595)
## rotate through testing sets
#wt_test_dumSplit <- rsample::initial_split(wt_test1_dummy, prop = 0.70)
#wt_test_dumSplit <- rsample::initial_split(wt_test2_dummy, prop = 0.70)
#wt_test_dumSplit <- rsample::initial_split(wt_test3_dummy, prop = 0.70)
#wt_test_dumSplit <- rsample::initial_split(wt_test4_dummy, prop = 0.70)
#wt_test_dumSplit <- rsample::initial_split(wt_test5_dummy, prop = 0.70)
#wt_test_dumSplit <- rsample::initial_split(wt_test6_dummy, prop = 0.70)
wt_test_dumSplit <- rsample::initial_split(wt_test7_dummy, prop = 0.70)
#wt_test_dumSplit <- rsample::initial_split(wt_test8_dummy, prop = 0.70)
#wt_test_dumSplit <- rsample::initial_split(wt_test9_dummy, prop = 0.70)
#wt_test_dumSplit <- rsample::initial_split(wt_test10_dummy, prop = 0.70)

trainData_wt <- rsample::training(wt_test_dumSplit)
testData_wt <- rsample::testing(wt_test_dumSplit)


## train models with split trainData
set.seed(2626)
m1w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "glm",
             family = "binomial",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m2w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "knn",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m3w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "nnet",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m4w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "glmnet",
             family = "binomial",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m5w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "ranger",
             importance = 'impurity',
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m6w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "svmRadial",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m7w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "bayesglm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m8w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "LogitBoost",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m9w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "rf",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m10w <- train(wt_loss_high ~ ., data = trainData_wt,
             method = "rpart",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m11w <- train(wt_loss_high ~ ., data = trainData_wt,
              method = "gbm",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

## resample metrics
resamps_w <- resamples(list(glm = m1w,
                             knn = m2w,
                             nnet = m3w,
                             glmnet = m4w,
                             ran = m5w,
                             svm = m6w,
                             bglm = m7w,
                             BlogReg = m8w,
                             rf = m9w,
                             rpart = m10w,
                             gbm = m11w))
resamps_w
summary(resamps_w)
bwplot(resamps_w)
dotplot(resamps_w)
modelCor(resamps_w)
splom(resamps_w)

ggplot(varImp(m1w))
ggplot(varImp(m3w))
ggplot(varImp(m4w))
ggplot(varImp(m5w))
ggplot(varImp(m9w))
ggplot(varImp(m10w))
ggplot(varImp(m11w)) ## need library(gbm)

varImp(m11w, scale = TRUE)

probTest <- predict(m11w, testData_wt, type = 'prob')
probTest <- factor(ifelse(probTest$no >= 0.73, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData_wt %>%
  na.omit() %>%
  dplyr::select(wt_loss_high) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$wt_loss_high, mode = "everything")

###

#algoList <- c('nnet', 'ranger', 'knn') 
algoList <- c('nnet', 'ranger', 'gbm') 
#algoList <- c('nnet', 'gbm') 

nnGrid <- expand.grid(size = seq(from = 15, to = 15, by = 1), 
                     decay = 0)

rGrid <- expand.grid(mtry = seq(from = 10, to = 10, by = 1), 
                     splitrule = 'extratrees',
                     min.node.size = seq(from = 1, to = 1, by = 1))

gGrid <- expand.grid(n.trees = seq(from = 400, to = 400, by = 50), 
            interaction.depth = seq(from = 5, to = 5, by = 1),
            shrinkage = seq(from = 0.1, to = 0.1, by = 0.1),
            n.minobsinnode = seq(from = 10, to = 10, by = 1))

set.seed(12345)
models_w1 <- caretList(wt_loss_high ~ ., 
                     data = trainData_wt, 
                     metric = "Balanced_Accuracy",
                     preProcess = c("nzv", "scale", "center"),
                     trControl = fitControl, 
                     #tuneLength = 10,
                     methodList = algoList,
                     tuneList = list(
                       nnet = caretModelSpec(method = "nnet",
                                             tuneGrid = nnGrid),
                       ranger = caretModelSpec(method = "ranger",
                                            tuneGrid = rGrid),
                       gbm = caretModelSpec(method = "gbm",
                                            tuneGrid = gGrid)))
summary(resamples(models_w1))
modelCor(resamples(models_w1))

## stack using random forest as the second order model
set.seed(12345)
rand_forest.stack_w1 <- caretStack(models_w1, 
                                method = "rf",
                                metric = "Balanced_Accuracy", 
                                tuneLength = 10, # mtry = 2
                                trControl = fitControl) 

## save the model
saveRDS(models_w1, "./final_models/models_w1_wt_loss_final_model.rds")
saveRDS(rand_forest.stack_w1, "./final_models/rand_forest.stack_w1_wt_loss_final_model.rds")
## read model back in
models_w1 <- readRDS("./final_models/models_w1_wt_loss_final_model.rds")
rand_forest.stack_w1 <- readRDS("./final_models/rand_forest.stack_w1_wt_loss_final_model.rds")


## stacked predictions
stack_predictions_w <- data.frame(
  rf_stack = predict(rand_forest.stack_w1, testData_wt))

## compile predictions with lethal truth
truth_predict_w <- testData_wt %>%
  na.omit() %>%
  dplyr::select(wt_loss_high) %>%
  cbind(., stack_predictions_w)

confusionMatrix(truth_predict_w$rf_stack, truth_predict_w$wt_loss_high, mode = "everything") 

###

probTest <- predict(rand_forest.stack_w1, testData_wt, type = 'prob')
probTest <- factor(ifelse(probTest >= 0.74, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData_wt %>%
  na.omit() %>%
  dplyr::select(wt_loss_high) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$wt_loss_high, mode = "everything")

########################################################################################





########################################################################################

## predictive modeling - weight loss - molecular

wt_mol_test1 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, '17':PBS)) %>%
  dplyr::select(-'329')

wt_mol_test2 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, HA, RBS, PBS))

wt_mol_test3 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, Origin, HA,
                  temp_5, AUC_6, '17':PBS))

wt_mol_test4 <- fullData_imputed %>% 
  drop_na(wt_loss_high) %>%
  dplyr::select(c(wt_loss_high, RBS, PBS))


## one hot code molecular data
wt_mol_test1_dummy <- fastDummies::dummy_cols(
  wt_mol_test1, select_columns =
    c('17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_mol_test2_dummy <- fastDummies::dummy_cols(
  wt_mol_test2, select_columns =
    c('HPAI_MBAA',	'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_mol_test3_dummy <- fastDummies::dummy_cols(
  wt_mol_test3, select_columns =
    c('HPAI_MBAA',	'HA', 'Origin', 
      '17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_mol_test4_dummy <- fastDummies::dummy_cols(
  wt_mol_test4, select_columns =
    c('RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_mol_test1_dummy$wt_loss_high <- as.factor(wt_mol_test1_dummy$wt_loss_high)
wt_mol_test2_dummy$wt_loss_high <- as.factor(wt_mol_test2_dummy$wt_loss_high)
wt_mol_test3_dummy$wt_loss_high <- as.factor(wt_mol_test3_dummy$wt_loss_high)
wt_mol_test4_dummy$wt_loss_high <- as.factor(wt_mol_test4_dummy$wt_loss_high)

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

## split data into train/test using the tidymodels/rsample package
set.seed(9595)
## rotate through testing sets
#wt_mol_test_dumSplit <- rsample::initial_split(wt_mol_test1_dummy, prop = 0.70)
#wt_mol_test_dumSplit <- rsample::initial_split(wt_mol_test2_dummy, prop = 0.70)
wt_mol_test_dumSplit <- rsample::initial_split(wt_mol_test3_dummy, prop = 0.70)
#wt_mol_test_dumSplit <- rsample::initial_split(wt_mol_test4_dummy, prop = 0.70)

trainData_wmol <- rsample::training(wt_mol_test_dumSplit)
testData_wmol <- rsample::testing(wt_mol_test_dumSplit)


## train models with split trainData
set.seed(2626)
m1wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "glm",
             family = "binomial",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m2wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "knn",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m3wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "nnet",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m4wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "glmnet",
             family = "binomial",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m5wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "ranger",
             importance = 'impurity',
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m6wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "svmRadial",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m7wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "bayesglm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m8wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "LogitBoost",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m9wm <- train(wt_loss_high ~ ., data = trainData_wmol,
             method = "rf",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

# set.seed(2626)
# m10wm <- train(wt_loss_high ~ ., data = trainData_wmol,
#               method = "rpart",
#               na.action = na.exclude,
#               preProcess = c("nzv", "scale", "center"),
#               trControl = fitControl)

set.seed(2626)
m11wm <- train(wt_loss_high ~ ., data = trainData_wmol,
              method = "gbm",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

## resample metrics
resamps_wm <- resamples(list(glm = m1wm,
                            knn = m2wm,
                            nnet = m3wm,
                            glmnet = m4wm,
                            ran = m5wm,
                            svm = m6wm,
                            bglm = m7wm,
                            BlogReg = m8wm,
                            rf = m9wm,
                            #rpart = m10wm,
                            gbm = m11wm))
summary(resamps_wm)
bwplot(resamps_wm)
dotplot(resamps_wm)
modelCor(resamps_wm)
splom(resamps_wm)

ggplot(varImp(m1wm))
ggplot(varImp(m3wm))
ggplot(varImp(m4wm))
ggplot(varImp(m5wm))
ggplot(varImp(m9wm))
ggplot(varImp(m10wm))
ggplot(varImp(m11wm)) ## need library(gbm)

varImp(m11wm, scale = TRUE)

## check which parameters can be changed
modelLookup(model = 'gbm') ## n.trees, interaction.depth, shrinkage, n.minobsinnode

## wt_mol_test3
m11wmGrid1 <- expand.grid(n.trees = seq(from = 50, to = 500, by = 50), 
                        interaction.depth = seq(from = 1, to = 4, by = 1),
                        shrinkage = seq(from = 0.1, to = 0.3, by = 0.1),
                        n.minobsinnode = seq(from = 8, to = 11, by = 1))

m11wmGrid2 <- expand.grid(n.trees = seq(from = 400, to = 500, by = 25), 
                          interaction.depth = seq(from = 4, to = 4, by = 1),
                          shrinkage = seq(from = 0.2, to = 0.2, by = 0.1),
                          n.minobsinnode = seq(from = 10, to = 11, by = 1))

m11wmGrid3 <- expand.grid(n.trees = c(400, 500), 
                          interaction.depth = seq(from = 4, to = 4, by = 1),
                          shrinkage = seq(from = 0.2, to = 0.6, by = 0.1),
                          n.minobsinnode = seq(from = 10, to = 11, by = 1))

m11wmGrid_final <- expand.grid(n.trees = 500, 
                          interaction.depth = 4,
                          shrinkage = 0.2,
                          n.minobsinnode = seq(from = 11, to = 11, by = 1))
## wt_mol_test1
m11wmGrid4 <- expand.grid(n.trees = seq(from = 250, to = 350, by = 50), 
                          interaction.depth = seq(from = 3, to = 4, by = 1),
                          shrinkage = seq(from = 0.3, to = 0.3, by = 0.1),
                          n.minobsinnode = seq(from = 10, to = 11, by = 1))

m11wmGrid5 <- expand.grid(n.trees = seq(from = 300, to = 300, by = 50), 
                          interaction.depth = seq(from = 4, to = 4, by = 1),
                          shrinkage = seq(from = 0.3, to = 0.3, by = 0.1),
                          n.minobsinnode = seq(from = 10, to = 10, by = 1))

set.seed(2626)
m11wm_tune <- train(wt_loss_high ~ ., data = trainData_wmol,
               method = "gbm",
               na.action = na.exclude,
               metric = "Balanced_Accuracy",
               preProcess = c("nzv", "scale", "center"),
               #tuneGrid = m11wmGrid_final, # test3
               tuneGrid = m11wmGrid5, # test1
               trControl = fitControl)

## test3 molecular + other variables
## save the model
saveRDS(m11wm_tune, "./final_models/m11wm_tune_wt_loss_molecular_all_final_model.rds")
## read model back in
m11wm_tune <- readRDS("./final_models/m11wm_tune_wt_loss_molecular_all_final_model.rds")

## test1 molecular only
## save the model
saveRDS(m11wm_tune, "./final_models/m11wm_tune_wt_loss_molecular_final_model.rds")
## read model back in
m11wm_tune <- readRDS("./final_models/m11wm_tune_wt_loss_molecular_final_model.rds")


probTest <- predict(m11wm_tune, testData_wmol, type = 'prob')
probTest <- factor(ifelse(probTest$no >= 0.88, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData_wmol %>%
  na.omit() %>%
  dplyr::select(wt_loss_high) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$wt_loss_high, mode = "everything")

########################################################################################





########################################################################################

## ROC plot

## get ROC curves for each cross-validation fold
## lethality models
model <- m11_tune
model <- m11m_tune
## wt_loss models
model <- m11wm_tune
  
for_lift <- data.frame(Class = model$pred$obs, 
                       gbm = model$pred$no, 
                       resample = model$pred$Resample)
lift_df <-  data.frame()
for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(Class ~ gbm, data = fold_df, class = "no")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}
lift_obj <- lift(Class ~ gbm, data = for_lift, class = "no")

## get the final model ROC 
for_lift_final <- data.frame(Class = model$pred$obs,  gbm = model$pred$no)
lift_obj_final <-  lift(Class ~ gbm, data = for_lift, class = 'no')

## plot both together
ROC_plot <- lift_df %>% 
  ggplot() +
  geom_line(aes(1-Sp , Sn, color = fold), size = 0.5, alpha = 0.6) +
  geom_line(data = lift_obj_final$data,
            aes(1-Sp , Sn), color = 'blue', size = 2) +
            #aes(1-Sp , Sn, color = liftModelVar), size = 2) +
  xlab('1 - Specificity') +
  ylab('Sensitivity') +
  scale_color_viridis_d(end = 0.85) +
  geom_abline(intercept = 0, slope = 1, 
              size = 1, linetype = 'dashed') +
  theme_bw() +
  theme(legend.position = 'none')

## m11_tune
ROC_plot + ggtitle('lethal GBM')

## m11m_tune
ROC_plot + ggtitle('lethal molecular GBM')

## m11wm_tune
ROC_plot + ggtitle('Wt_loss molecular GBM')

###

## wt_loss models
model <- rand_forest.stack_w1

for_lift <- data.frame(Class = model$ens_model$pred$obs, 
                       ensemble = model$ens_model$pred$no, 
                       resample = model$ens_model$pred$Resample)
lift_df <-  data.frame()
for (fold in unique(for_lift$resample)) {
  fold_df <- dplyr::filter(for_lift, resample == fold)
  lift_obj_data <- lift(Class ~ ensemble, data = fold_df, class = "no")$data
  lift_obj_data$fold = fold
  lift_df = rbind(lift_df, lift_obj_data)
}
lift_obj <- lift(Class ~ ensemble, data = for_lift, class = "no")

## get the final model ROC 
for_lift_final <- data.frame(Class = model$ens_model$pred$obs,  ensemble = model$ens_model$pred$no)
lift_obj_final <-  lift(Class ~ ensemble, data = for_lift, class = 'no')

ROC_plot <- lift_df %>% 
  ggplot() +
  geom_line(aes(1-Sp , Sn, color = fold), size = 0.5, alpha = 0.6) +
  geom_line(data = lift_obj_final$data,
            aes(1-Sp , Sn), color = 'blue', size = 2) +
  #aes(1-Sp , Sn, color = liftModelVar), size = 2) +
  xlab('1 - Specificity') +
  ylab('Sensitivity') +
  scale_color_viridis_d(end = 0.85) +
  geom_abline(intercept = 0, slope = 1, 
              size = 1, linetype = 'dashed') +
  theme_bw() +
  theme(legend.position = 'none')

## rand_forest.stack_w1 - nnet-ranger-gbm stacked w/ rf
ROC_plot + ggtitle('wt_loss ensemble model')

########################################################################################





########################################################################################

## predictive models - RD transmission

## add AUC
RD_trans_dataFilt <- 
  fullData_imputed %>% 
  filter(expt != "path") #%>%

RD_TRANSMISSION <- 'RD_trans'

###

RD_test1 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, HPAI_MBAA, HA, RBS, PBS)) 

## HPAI_MBAA is not important, swap for Origin
RD_test2 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, HA, RBS, PBS)) 

## add d1-d11 AUC
RD_test3 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, HA, RBS, PBS, AUC_RD)) 

## drop AUC and swap HA for Subtype
## Subtype performs better than HA
RD_test4 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, Subtype, RBS, PBS)) 

RD_test5 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, Subtype, RBS, PBS, 
                  peak_inoc, slope13)) 

## best so far
RD_test6 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, Subtype, RBS, PBS, 
                  AUC_6, slope13)) 

RD_test7 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, Subtype, RBS, PBS, 
                  AUC_4, slope13)) 

## test6 wtih pH fusion added - performs better, but not as many data points. Need more pH data. 
RD_test8 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, Subtype, RBS, PBS, 
                  AUC_6, slope13, pH.for.fusion)) %>%
  drop_na(pH.for.fusion)

## test8 adding the top two molecular markers - minor improvement
RD_test9 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, Subtype, RBS, PBS, 
                  AUC_6, slope13, pH.for.fusion, PB2_627, '155')) %>%
  drop_na(pH.for.fusion)

## almost just as good as tests 6/8. Slope13 is unimportant in this model
RD_test10 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, AUC_6, slope13, pH.for.fusion)) %>%
  drop_na(pH.for.fusion)

## test10 drop slope13 - pH fusion becomes unimportant
RD_test11 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, AUC_6, pH.for.fusion)) %>%
  drop_na(pH.for.fusion)

## test11 drop pH fusion - AUC_6 only performs fairly well verifying importance
RD_test12 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, AUC_6)) 

## test12 - swap for AUC_4
RD_test13 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, AUC_4)) 

## test6 - swap subtype for HA
RD_test14 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, HA, RBS, PBS, 
                  AUC_6, slope13)) 

## test7 - swap AUC_6 for AUC_RD
RD_test15 <- RD_trans_dataFilt %>% 
  drop_na(RD_TRANSMISSION) %>%
  dplyr::select(c(RD_TRANSMISSION, Origin, HA, RBS, PBS, 
                  AUC_RD, slope13))

###

RD_test1_dummy <- fastDummies::dummy_cols(
  RD_test1, select_columns =
    c('HPAI_MBAA', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test2_dummy <- fastDummies::dummy_cols(
  RD_test2, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test3_dummy <- fastDummies::dummy_cols(
  RD_test3, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test4_dummy <- fastDummies::dummy_cols(
  RD_test4, select_columns =
    c('Origin', 'Subtype', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test5_dummy <- fastDummies::dummy_cols(
  RD_test5, select_columns =
    c('Origin', 'Subtype', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test6_dummy <- fastDummies::dummy_cols(
  RD_test6, select_columns =
    c('Origin', 'Subtype', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test7_dummy <- fastDummies::dummy_cols(
  RD_test7, select_columns =
    c('Origin', 'Subtype', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test8_dummy <- fastDummies::dummy_cols(
  RD_test8, select_columns =
    c('Origin', 'Subtype', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test9_dummy <- fastDummies::dummy_cols(
  RD_test9, select_columns =
    c('Origin', 'Subtype', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test14_dummy <- fastDummies::dummy_cols(
  RD_test14, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test15_dummy <- fastDummies::dummy_cols(
  RD_test15, select_columns =
    c('Origin', 'HA', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_test1_dummy$RD_trans_ext <- as.factor(RD_test1_dummy$RD_trans_ext)
RD_test2_dummy$RD_trans_ext <- as.factor(RD_test2_dummy$RD_trans_ext)
RD_test3_dummy$RD_trans_ext <- as.factor(RD_test3_dummy$RD_trans_ext)
RD_test4_dummy$RD_trans_ext <- as.factor(RD_test4_dummy$RD_trans_ext)
RD_test5_dummy$RD_trans_ext <- as.factor(RD_test5_dummy$RD_trans_ext)
RD_test6_dummy$RD_trans_ext <- as.factor(RD_test6_dummy$RD_trans_ext)
RD_test7_dummy$RD_trans_ext <- as.factor(RD_test7_dummy$RD_trans_ext)
RD_test8_dummy$RD_trans_ext <- as.factor(RD_test8_dummy$RD_trans_ext)
RD_test9_dummy$RD_trans_ext <- as.factor(RD_test9_dummy$RD_trans_ext)
RD_test10$RD_trans_ext <- as.factor(RD_test10$RD_trans_ext)
RD_test11$RD_trans_ext <- as.factor(RD_test11$RD_trans_ext)
RD_test12$RD_trans_ext <- as.factor(RD_test12$RD_trans_ext)
RD_test13$RD_trans_ext <- as.factor(RD_test13$RD_trans_ext)
RD_test14_dummy$RD_trans_ext <- as.factor(RD_test14_dummy$RD_trans_ext)
RD_test15_dummy$RD_trans_ext <- as.factor(RD_test15_dummy$RD_trans_ext)

###

## Predictive Power Scores
set.seed(123)
RD_trans_dataFilt %>%
  dplyr::select(RD_trans_ext, AUC_4, AUC_6, AUC_8, AUC_9, AUC_RD) %>%
  ppsr::score_predictors(df = ., y = 'RD_trans_ext')

###

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

## split data into train/test using the tidymodels/rsample package
set.seed(9595)
## rotate through testing sets
#RD_test_dumSplit <- rsample::initial_split(RD_test1_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test2_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test3_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test4_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test5_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test6_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test7_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test8_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test9_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test10, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test11, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test12, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test13, prop = 0.70)
RD_test_dumSplit <- rsample::initial_split(RD_test14_dummy, prop = 0.70)
#RD_test_dumSplit <- rsample::initial_split(RD_test15_dummy, prop = 0.70)

trainData_RD <- rsample::training(RD_test_dumSplit)
testData_RD <- rsample::testing(RD_test_dumSplit)

## train models with split trainData
FORMULA <- RD_trans ~ .

set.seed(2626)
m1RD <- train(FORMULA, data = trainData_RD,
               method = "glm",
               family = "binomial",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m2RD <- train(FORMULA, data = trainData_RD,
               method = "knn",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m3RD <- train(FORMULA, data = trainData_RD,
               method = "nnet",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m4RD <- train(FORMULA, data = trainData_RD,
               method = "glmnet",
               family = "binomial",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m5RD <- train(FORMULA, data = trainData_RD,
               method = "ranger",
               importance = 'impurity',
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m6RD <- train(FORMULA, data = trainData_RD,
               method = "svmRadial",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m7RD <- train(FORMULA, data = trainData_RD,
               method = "bayesglm",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m8RD <- train(FORMULA, data = trainData_RD,
               method = "LogitBoost",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m9RD <- train(FORMULA, data = trainData_RD,
               method = "rf",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

# set.seed(2626)
# m10RD <- train(FORMULA, data = trainData_RD,
#               method = "rpart",
#               na.action = na.exclude,
#               preProcess = c("nzv", "scale", "center"),
#               trControl = fitControl)

set.seed(2626)
m11RD <- train(FORMULA, data = trainData_RD,
                method = "gbm",
                na.action = na.exclude,
                preProcess = c("nzv", "scale", "center"),
                trControl = fitControl)

## resample metrics
resamps_RD <- resamples(list(glm = m1RD,
                              knn = m2RD,
                              nnet = m3RD,
                              glmnet = m4RD,
                              ran = m5RD,
                              svm = m6RD,
                              bglm = m7RD,
                              BlogReg = m8RD,
                              rf = m9RD,
                              #rpart = m10RD,
                              gbm = m11RD))
summary(resamps_RD)
bwplot(resamps_RD)
dotplot(resamps_RD)
modelCor(resamps_RD)
splom(resamps_RD)

ggplot(varImp(m1RD))
ggplot(varImp(m3RD))
ggplot(varImp(m4RD))
ggplot(varImp(m5RD))
ggplot(varImp(m9RD))
ggplot(varImp(m11RD)) ## need library(gbm)

probTest <- predict(m9RD, testData_RD, type = 'prob')
probTest <- factor(ifelse(probTest$no >= 0.7, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData_RD %>%
  na.omit() %>%
  dplyr::select(RD_TRANSMISSION) %>%
  cbind(., probTest)
probTruth$RD_trans_ext <- as.factor(probTruth$RD_trans_ext)
confusionMatrix(probTruth$., probTruth$RD_trans_ext, mode = "everything")


###

## best model tuning

modelLookup(model = 'rf') ## mtry

## RD_trans_test12
m9RDGrid1 <- expand.grid(mtry = seq(from = 1, to = 10, by = 1))
m9RDGrid2 <- expand.grid(mtry = seq(from = 8, to = 20, by = 1))
m9RDGrid_final <- expand.grid(mtry = seq(from = 16, to = 16, by = 1))


set.seed(2626)
m9RD_tune <- train(FORMULA, data = trainData_RD,
                    method = "rf",
                    na.action = na.exclude,
                    metric = "Balanced_Accuracy",
                    preProcess = c("nzv", "scale", "center"),
                    tuneGrid = m9RDGrid_final,
                    trControl = fitControl)

saveRDS(m9RD_tune, "./final_models/m9RD_tune_RD_final_model.rds")
m9RD_tune <- readRDS("./final_models/m9RD_tune_RD_final_model.rds")

probTest <- predict(m9RD_tune, testData_RD, type = 'prob')
probTest <- factor(ifelse(probTest$no >= 0.78, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData_RD %>%
  na.omit() %>%
  dplyr::select(RD_TRANSMISSION) %>%
  cbind(., probTest)
probTruth$RD_trans_ext <- as.factor(probTruth$RD_trans_ext)
confusionMatrix(probTruth$., probTruth$RD_trans_ext, mode = "everything")


########################################################################################

## predictive models - RD transmission molecular

RD_mol_test1 <- RD_trans_dataFilt %>% 
  drop_na(RD_trans_ext) %>%
  dplyr::select(c(RD_trans_ext, '17':PBS)) %>%
  dplyr::select(-'329')

RD_mol_test2 <- RD_trans_dataFilt %>% 
  drop_na(RD_trans_ext) %>%
  dplyr::select(c(RD_trans_ext, RBS:PBS)) 

RD_mol_test3 <- RD_trans_dataFilt %>% 
  drop_na(RD_trans_ext) %>%
  dplyr::select(c(RD_trans_ext, Subtype, Origin, '17':PBS)) %>%
  dplyr::select(-'329')

RD_mol_test4 <- RD_trans_dataFilt %>% 
  drop_na(RD_trans_ext) %>%
  dplyr::select(c(RD_trans_ext, Subtype, Origin, '21', '83', 
                  '138', '155', '186', '214', '225', '226',
                  'PB2_627', RBS, PBS)) 

RD_mol_test5 <- RD_trans_dataFilt %>% 
  drop_na(RD_trans_ext) %>%
  dplyr::select(c(RD_trans_ext, Subtype, 'PB2_627', RBS)) 

RD_mol_test6 <- RD_trans_dataFilt %>% 
  drop_na(RD_trans_ext) %>%
  dplyr::select(c(RD_trans_ext, HA, 'PB2_627', RBS)) 

## test3 - swap subtype for HA
RD_mol_test7 <- RD_trans_dataFilt %>% 
  drop_na(RD_trans_ext) %>%
  dplyr::select(c(RD_trans_ext, HA, Origin, '17':PBS)) %>%
  dplyr::select(-'329')

RD_mol_test8 <- RD_trans_dataFilt %>% 
  drop_na(RD_trans_ext) %>%
  dplyr::select(c(RD_trans_ext, Origin, HA, 
                  AUC_6, slope13, '17':PBS)) 

###

RD_mol_test1_dummy <- fastDummies::dummy_cols(
  RD_mol_test1, select_columns =
    c('17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_mol_test2_dummy <- fastDummies::dummy_cols(
  RD_mol_test2, select_columns =
    c('RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_mol_test3_dummy <- fastDummies::dummy_cols(
  RD_mol_test3, select_columns =
    c('Subtype', 'Origin',
      '17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_mol_test4_dummy <- fastDummies::dummy_cols(
  RD_mol_test4, select_columns =
    c('Subtype', 'Origin', '21', '83', 
      '138', '155', '186', '214', '225', 
      '226', 'PB2_627', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_mol_test5_dummy <- fastDummies::dummy_cols(
  RD_mol_test5, select_columns =
    c('Subtype', 'PB2_627', 'RBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_mol_test6_dummy <- fastDummies::dummy_cols(
  RD_mol_test6, select_columns =
    c('HA', 'PB2_627', 'RBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_mol_test7_dummy <- fastDummies::dummy_cols(
  RD_mol_test7, select_columns =
    c('HA', 'Origin',
      '17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_mol_test8_dummy <- fastDummies::dummy_cols(
  RD_mol_test8, select_columns =
    c('HA', 'Origin',
      '17', '18', '21', '83', '110',	'126', '128',
      '137', '138', '143',	'155', '158',	'160', '186',
      '187',	'190', '192', '193',	'196',	'197',
      '214', '225', '226', '227',	'228',	'239',	'255',
      '318', '387',	'443',	'446',	'496', 'PB2_271',
      'PB2_590_591', 'PB2_627',	'PB2_701', 'RBS', 'PBS'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

RD_mol_test1_dummy$RD_trans_ext <- as.factor(RD_mol_test1_dummy$RD_trans_ext)
RD_mol_test2_dummy$RD_trans_ext <- as.factor(RD_mol_test2_dummy$RD_trans_ext)
RD_mol_test3_dummy$RD_trans_ext <- as.factor(RD_mol_test3_dummy$RD_trans_ext)
RD_mol_test4_dummy$RD_trans_ext <- as.factor(RD_mol_test4_dummy$RD_trans_ext)
RD_mol_test5_dummy$RD_trans_ext <- as.factor(RD_mol_test5_dummy$RD_trans_ext)
RD_mol_test6_dummy$RD_trans_ext <- as.factor(RD_mol_test6_dummy$RD_trans_ext)
RD_mol_test7_dummy$RD_trans_ext <- as.factor(RD_mol_test7_dummy$RD_trans_ext)
RD_mol_test8_dummy$RD_trans_ext <- as.factor(RD_mol_test8_dummy$RD_trans_ext)

## split data into train/test using the tidymodels/rsample package
set.seed(9595)
## rotate through testing sets
RD_mol_test_dumSplit <- rsample::initial_split(RD_mol_test1_dummy, prop = 0.70)
#RD_mol_test_dumSplit <- rsample::initial_split(RD_mol_test2_dummy, prop = 0.70)
#RD_mol_test_dumSplit <- rsample::initial_split(RD_mol_test3_dummy, prop = 0.70)
#RD_mol_test_dumSplit <- rsample::initial_split(RD_mol_test4_dummy, prop = 0.70)
#RD_mol_test_dumSplit <- rsample::initial_split(RD_mol_test5_dummy, prop = 0.70)
#RD_mol_test_dumSplit <- rsample::initial_split(RD_mol_test6_dummy, prop = 0.70)
#RD_mol_test_dumSplit <- rsample::initial_split(RD_mol_test7_dummy, prop = 0.70)
#RD_mol_test_dumSplit <- rsample::initial_split(RD_mol_test8_dummy, prop = 0.70)

trainData_RDmol <- rsample::training(RD_mol_test_dumSplit)
testData_RDmol <- rsample::testing(RD_mol_test_dumSplit)


## train models with split trainData
set.seed(2626)
m1RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "glm",
              family = "binomial",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

set.seed(2626)
m2RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "knn",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

set.seed(2626)
m3RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "nnet",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

set.seed(2626)
m4RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "glmnet",
              family = "binomial",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

set.seed(2626)
m5RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "ranger",
              importance = 'impurity',
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

set.seed(2626)
m6RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "svmRadial",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

set.seed(2626)
m7RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "bayesglm",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

set.seed(2626)
m8RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "LogitBoost",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

set.seed(2626)
m9RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
              method = "rf",
              na.action = na.exclude,
              preProcess = c("nzv", "scale", "center"),
              trControl = fitControl)

# set.seed(2626)
# m10RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
#               method = "rpart",
#               na.action = na.exclude,
#               preProcess = c("nzv", "scale", "center"),
#               trControl = fitControl)

set.seed(2626)
m11RDm <- train(RD_trans_ext ~ ., data = trainData_RDmol,
               method = "gbm",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

## resample metrics
resamps_RDm <- resamples(list(glm = m1RDm,
                             knn = m2RDm,
                             nnet = m3RDm,
                             glmnet = m4RDm,
                             ran = m5RDm,
                             svm = m6RDm,
                             bglm = m7RDm,
                             BlogReg = m8RDm,
                             rf = m9RDm,
                             #rpart = m10RDm,
                             gbm = m11RDm))
summary(resamps_RDm)
bwplot(resamps_RDm)
dotplot(resamps_RDm)
modelCor(resamps_RDm)
splom(resamps_RDm)

ggplot(varImp(m1RDm))
ggplot(varImp(m3RDm))
ggplot(varImp(m4RDm))
ggplot(varImp(m5RDm))
ggplot(varImp(m9RDm))
ggplot(varImp(m11RDm)) ## need library(gbm)

probTest <- predict(m9RDm, testData_RDmol, type = 'prob')
probTest <- factor(ifelse(probTest$no >= 0.73, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData_RDmol %>%
  na.omit() %>%
  dplyr::select(RD_trans_ext) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$RD_trans_ext, mode = "everything")


###


modelLookup(model = 'rf') ## mtry

## RD_trans_mol_test7
m9RDmGrid1 <- expand.grid(mtry = seq(from = 100, to = 300, by = 10))
m9RDmGrid2 <- expand.grid(mtry = seq(from = 185, to = 195, by = 1))
m9RDmGrid_final <- expand.grid(mtry = seq(from = 192, to = 192, by = 1))

## RD_trans_mol_test1
m9RDmGrid1 <- expand.grid(mtry = seq(from = 100, to = 300, by = 10))
m9RDmGrid2 <- expand.grid(mtry = seq(from = 265, to = 275, by = 1))
m9RDmGrid_final <- expand.grid(mtry = seq(from = 270, to = 270, by = 1))

set.seed(2626)
m9RDm_tune <- train(FORMULA, data = trainData_RDmol,
                   method = "rf",
                   na.action = na.exclude,
                   metric = "Balanced_Accuracy",
                   preProcess = c("nzv", "scale", "center"),
                   tuneGrid = m9RDmGrid_final,
                   trControl = fitControl)

## RD_trans_mol_test1
saveRDS(m9RDm_tune, "./final_models/m9RDm_tune_RD_molecular_final_model.rds")
m9RDm_tune <- readRDS("./final_models/m9RDm_tune_RD_molecular_final_model.rds")
## RD_trans_mol_test7
saveRDS(m9RDm_tune, "./final_models/m9RDm_tune_RD_molecular_all_final_model.rds")
m9RDm_tune <- readRDS("./final_models/m9RDm_tune_RD_molecular_all_final_model.rds")

###

probTest <- predict(m9RDm_tune, testData_RDmol, type = 'prob')
probTest <- factor(ifelse(probTest$no >= 0.85, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData_RDmol %>%
  na.omit() %>%
  dplyr::select(RD_trans_ext) %>%
  cbind(., probTest)
probTruth$RD_trans_ext <- as.factor(probTruth$RD_trans_ext)
confusionMatrix(probTruth$., probTruth$RD_trans_ext, mode = "everything")


########################################################################################

## Matthews Correlation Coefficient

## manual calculation/function
## T = truth, F = false, P = positive, N = negative
## use data from confusion matrix
mcc_func <- function(TP, FP, FN, TN){
  ((TP*TN) - (FP*FN))/
    sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))}

## lethal
mcc_func(169, 2, 14, 31)
# 0.7644239
## lethal molecular
mcc_func(176, 4, 7, 29)
# 0.811431
## molecular only
mcc_func(174, 6, 9, 27)
# 0.742373

## lethal simData
mcc_func(1718, 235, 178, 624)
# 0.6450128
## lethal standardizedData
mcc_func(79, 1, 2, 2)
# 0.5594309


## weight loss
mcc_func(122, 13, 41, 39)
# 0.4415901
## weight loss molecular
mcc_func(120, 12, 43, 40)
# 0.4445575
## molecular only
mcc_func(112, 8, 51, 44)
# 0.4598228

## weight loss simData
mcc_func(1375, 382, 566, 432)
# 0.2269719
## weight loss standardizedData
mcc_func(37, 9, 28, 10)
# 0.08031151


## RD transmission models
mcc_func(68, 0, 3, 63)
# 0.9561446
mcc_func(71, 0, 0, 63)
# 1
mcc_func(72, 0, 1, 70)
# 0.9861084
mcc_func(71, 0, 2, 70)
# 0.9724125
mcc_func(67, 0, 6, 70)
# 0.919429


########################################################################################

## final models metric heatmap
modelMetrics <- read.csv("inputs/ModelMetrics4Heatmap.csv", header = TRUE)

ggplot(modelMetrics, aes(x = Model, y = Value, fill = Metric)) + 
  geom_bar(stat = "identity", position = "dodge", show.legend = TRUE) + 
  facet_grid(~ Group) +
  scale_fill_viridis_d(end = 0.8) +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.title.x = element_text (size = 20, face = "bold"), 
        axis.title.y = element_text (size = 20, face = "bold"),
        axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1, 
                                   size = 18, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 18)) +
  guides(fill = guide_legend(nrow = 1))

###

modelMetrics <- modelMetrics %>% 
  gather(L1:TM, key = Model, value = Value, factor_key = TRUE) %>%
  dplyr::filter(Metric != 'test') %>%
  drop_na()

modelMetrics$Metric <- factor(modelMetrics$Metric, 
                              levels=c('Accuracy', 'BalancedAccuracy', 
                                            'Sensitivity', 'Specificity', 'F1',
                                            'Precision', 'Recall', 'Kappa',
                                            'MCC', 'Probability'))

modelMetrics_In <- modelMetrics %>% 
  filter(Model %in% c('L1', 'L1M', 'LM',
                      'M1', 'M1M', 'MM',
                      'T1', 'T1M', 'TM'))

modelMetrics_Ex <- modelMetrics %>% 
  filter(Model %in% c('L1.H1N1', 'L1.sim', 'LM.pub',
                      'M1.H1N1', 'M1.sim', 'TM.pub'))

modelMetrics_In$Model <- factor(modelMetrics_In$Model, 
                             levels=c('L1', 'L1M', 'LM',
                                      'M1', 'M1M', 'MM',
                                      'T1', 'T1M', 'TM'))

modelMetrics_Ex$Model <- factor(modelMetrics_Ex$Model, 
                                levels=c('L1.H1N1', 'L1.sim', 'LM.pub',
                                         'M1.H1N1', 'M1.sim', 'TM.pub'))

modelMetrics_In$Value <- as.numeric(modelMetrics_In$Value)
modelMetrics_Ex$Value <- as.numeric(modelMetrics_Ex$Value)

ggplot(modelMetrics_In, aes(Model, Metric, fill = Value)) + 
  geom_tile() +
  geom_text(aes(label = Value), color = "white", size = 3) +
  scale_fill_viridis_c(direction = -1, end = 0.85) +
  #scale_x_discrete(sec.axis = dup_axis()) +
  theme_minimal() +
  theme(legend.position = 'top') 

ggplot(modelMetrics_Ex, aes(Model, Metric, fill = Value)) + 
  geom_tile() +
  geom_text(aes(label = Value), color = "white", size = 3) +
  scale_fill_viridis_c(direction = -1, end = 0.85) +
  theme_minimal() +
  theme(legend.position = 'top') 

########################################################################################

## example of lethality balanced accuracy mean model metric for each set of features for all 11 models
lethalBA_metrics <- read.csv("inputs/FullModelMetricExample.csv", header = TRUE) 

ggplot(lethalBA_metrics, aes(Test, Model, fill = Mean)) + 
  geom_tile() +
  geom_text(aes(label = Mean), color = "white", size = 3) +
  scale_fill_viridis_c(direction = -1, end = 0.85) +
  #scale_x_discrete(sec.axis = dup_axis()) +
  theme_minimal() +
  theme(legend.position = 'none',
        axis.text.x = element_blank()) 

########################################################################################

# features <- read.csv("inputs/Features4Heatmap.csv", header = TRUE)
# 
# ## Select relevant columns for creating the presence-absence matrix
# presence_matrix <- features[, 5:16] != ""
# ## Convert to binary matrix format (True for presence, False for absence)
# binary_matrix <- as.matrix(presence_matrix)
# ## drop text Feature columns
# features <- features %>% dplyr::select(-starts_with('Feature'))
# ## combine binary to features
# features_binary <- cbind(features, binary_matrix)   
# ## convert TRUE/FALSE to 1/0
# features_binary <- features_binary %>%
#   mutate(across(where(is.logical), ~ if_else(.x, 1, 0)))
# 
# ggplot(data = features_binary) +
#   geom_tile(aes(y = Test, x = Classification, fill = Feature1 + Feature2 + Feature3 + Feature4 +
#                   Feature5 + Feature6 + Feature7 + Feature8 +
#                   Feature9 + Feature10 + Feature11 + Feature12)) +
#   labs(x = "Feature", y = "Test") +                           
#   theme_minimal()
  
###
  
features2 <- read.csv("inputs/Features4Heatmap2.csv", header = TRUE)

## Select relevant columns for creating the presence-absence matrix
presence_matrix <- features2[, 5:24] != ""
## Convert to binary matrix format (True for presence, False for absence)
binary_matrix <- as.matrix(presence_matrix)
## drop text Feature columns
features2 <- features2[, 1:4]
## combine binary to features
features_binary <- cbind(features2, binary_matrix)   
## convert TRUE/FALSE to 1/0
features_binary <- features_binary %>%
  mutate(across(where(is.logical), ~ if_else(.x, 1, 0)))

gathered_data <- features_binary %>%
  gather(key = "Feature", value = "Presence", 5:24)

ggplot(data = gathered_data) +
  geom_tile(aes(y = Feature, x = Test, fill = Presence)) +
  labs(x = "Feature", y = "Test") +  
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  facet_grid(Type ~ Classification, scales = 'free', space = 'free') +
  theme(legend.position = 'none',
        axis.text.x = element_blank()) 

## should separate out the Types and then patchwork together rather than facet.

gathered_data$Test <- as.factor(gathered_data$Test)

## subset data
Standard_data <- gathered_data %>% filter(Type == 'Standard')
Molecular_data <- gathered_data %>% filter(Type == 'Molecular')
Combined_data <- gathered_data %>% filter(Type == 'Combined')

## remove Features that were not used at all per group
Standard_data <- Standard_data %>% 
  subset(Feature != 'All_molecular')

Molecular_data <- Molecular_data %>% 
  group_by(Classification, Feature) %>% 
  mutate(SUM = sum(Presence)) %>%
  filter(SUM != 0)

Combined_data <- Combined_data %>% 
filter(!Feature %in% c('temp', 'peak_inoc', 'pH.for.fusion',
                       'd3_inoc', 'Origin_orig', 'AUC_RD', 
                       'AUC_9', 'AUC_4'))

## make subplots
Standard_plot <- 
  ggplot(data = Standard_data) +
  geom_tile(aes(y = Feature, x = Test, fill = Presence)) +
  labs(x = "Test", y = "Standard") +   
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  facet_grid(~ Classification, scales = 'free', space = 'free') +
  theme(legend.position = 'none') 

Molecular_plot <- 
  ggplot(data = Molecular_data) +
  geom_tile(aes(y = Feature, x = Test, fill = Presence)) +
  labs(x = "Test", y = "Molecular") +   
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  facet_grid(~ Classification, scales = 'free', space = 'free') +
  theme(legend.position = 'none') 

Combined_plot <- 
  ggplot(data = Combined_data) +
  geom_tile(aes(y = Feature, x = Test, fill = Presence)) +
  labs(x = "Test", y = "Combined") +  
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  facet_grid(~ Classification, scales = 'free', space = 'free') +
  theme(legend.position = 'none') 

## combine subplots
Standard_plot / Molecular_plot / Combined_plot + 
  patchwork::plot_layout(heights = c(1.8, 0.6, 1.2))

###

## reconfigure figure to combine example balanced accuracy values
## and features used for L1 model

# Standard_data2 <- Standard_data %>%
#   filter(Feature %in% c('HPAI_MBAA', 'RBS', 'PBS', 'HA',
#                         'wt_loss', 'temp', 'AUC_6')) %>%
#   filter(Classification == 'Lethality')

Standard_data2 <- Standard_data %>%
  filter(!(Feature %in% c('pH.for.fusion', 'd3_inoc', 'AUC_RD', 
                          'AUC_9', 'AUC_4', 'Other_molecular',
                          'RBS', 'PA', 'MBAA', 'AUC_6', 'wt_loss'))) %>%
  filter(Classification == 'Lethality')

lethalBA_metrics2 <- lethalBA_metrics %>% 
  filter(Model %in% c('gbm', 'rf', 'svm', 'rpart'))
  
Standard_plot2 <- 
  ggplot(data = Standard_data2, 
         aes(y = Feature, x = Test, fill = Presence)) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  labs(x = "Model Test", y = "Feature") +   
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  theme(legend.position = 'none') 

lethalBA_metrics_plot <- 
  ggplot(lethalBA_metrics2, 
         aes(x = Test, y = Model, fill = Mean)) + 
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  ylab("Algorithm") + 
  geom_text(aes(label = Mean), color = "white", size = 3) +
  scale_fill_viridis_c(direction = -1, end = 0.85) +
  theme_minimal() +
  theme(legend.position = 'top',
        #axis.text.x = element_blank(),
        axis.title.x = element_blank()) 

lethalBA_metrics_plot / Standard_plot2 + 
  plot_annotation(tag_levels = 'A')

########################################################################################

## plot with Standard models (L1, S1, T1) features heatmap and 
## barplot with selected importance variables

Standard_data3 <- Standard_data %>%
  filter((Classification == 'Lethality'& Test == 11) |
           (Classification == 'Morbidity'& Test == 7) |
           (Classification == 'Transmission'& Test == 14)) %>%
  filter(!(Feature %in% c('pH.for.fusion', 'd3_inoc', 'AUC_RD', 
                          'AUC_9', 'AUC_4', 'peak_inoc', 'temp',
                          'Other_molecular','Subtype', 'Origin_orig'))) %>%
  mutate(Test = case_when(Test == 11 ~ 'L1',
                          Test == 7 ~ 'M1',
                          Test == 14 ~ 'T1'))
  
Standard_plot3 <- 
  ggplot(data = Standard_data3, 
         aes(y = Feature, x = Test, fill = Presence)) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  labs(x = "Model Test", y = "Feature") +   
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  theme(legend.position = 'none') +
  coord_flip()


SelectedImportance <- read.csv("inputs/SelectedImportance.csv", header = TRUE)

SelectedImportance_num <- SelectedImportance %>%
  filter(Catergory == 'Numeric') %>%
  mutate(Feature = fct_reorder(factor(Feature), Importance, .desc = FALSE))

SelectedImportance_fac <- SelectedImportance %>%
  filter(Catergory != 'Numeric')

SelectedImportance_num_plot <- 
  ggplot(SelectedImportance_num, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model, scales = "free_y") +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'B', begin = 0.1, end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

SelectedImportance_fac_plot <- 
  ggplot(SelectedImportance_fac, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model) +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'G', end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

Standard_plot3 / 
  SelectedImportance_num_plot / 
  SelectedImportance_fac_plot + 
  plot_annotation(tag_levels = 'A')


########################################################################################
########################################################################################
### End of Code
########################################################################################
########################################################################################
