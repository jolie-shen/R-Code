library(dplyr)
library(exact2x2)
library(fmsb)
library(mice)
library(PropCIs)
library(questionr)
library(regclass)
library(ROCit)
library(stats)
library(tidyverse)
library(MASS)
library(hash)

# Boolean for deleting the 6 rows that are 90% "Unknown"
listwise_deletion <- TRUE

# Should p-values in univariate use imputed data
USE_IMPUTATIONS_FOR_UNIVARIATE <- TRUE

# Decide whether the NAs in the race category are "other" or are missing,
# and should therefore be imputed
TREAT_RACE_NAS_AS_OTHER <- TRUE

# Load the CSV
file_path <- "~/Downloads/clade1.csv"
data <- read.csv(
    file = file_path, 
    header = TRUE, 
    na.strings = c("Unknown"), 
    strip.white = TRUE
)

# Remove clades of zero--we just want a binary comparison. Only 2 rows
data <- data %>% filter(clade != 0)

if (TREAT_RACE_NAS_AS_OTHER) {
    levels(data$Race) <- c(levels(data$Race), "Other")
    data$Race[is.na(data$Race)] <- "Other"
}

data <- data %>% mutate(IsClade1 = ifelse(clade == 1, TRUE, FALSE))

# Set factor variables
data <- data %>%
    mutate(
        Gender = as.factor(Gender),
        Race = as.factor(Race),
        SNF = as.factor(SNF),
        Hospitalized.for.COVID = as.factor(Hospitalized.for.COVID),
        Diabetes = as.factor(Diabetes),
        HTN = as.factor(HTN),
        COPD = as.factor(COPD),
        Asthma = as.factor(Asthma),
        Hx.of.DVT = as.factor(Hx.of.DVT),
        CKD = as.factor(CKD),
        Cancer = as.factor(Cancer),
        Hx.of.MI = as.factor(Hx.of.MI),
        CVD = as.factor(CVD),
        CHF = as.factor(CHF),
        Hypothyroid = as.factor(Hypothyroid),
        Steroids.or.IMT = as.factor(Steroids.or.IMT),
        Baseline.Plaquenil = as.factor(Baseline.Plaquenil),
        ACEI.ARB = as.factor(ACEI.ARB),
        Smoking.History. = as.factor(Smoking.History.),
        Anticoagulation. = as.factor(Anticoagulation.),
        Death = as.factor(Death)
    )

data <- data[, !names(data) %in% c("clade")] 

# Relevel to reference groups
data$Race <- relevel(data$Race, ref = "White")
data$Smoking.History. <- relevel(data$Smoking.History., ref = "Never")

# Set random seed
random_seed_num <- 3249
set.seed(random_seed_num)

# This number was originally set by Rubin, and 5 was believed to be enough. 
# Since then, Bodner (2008), White et al. (2011) and Graham, Olchowski, 
# and Gilreath (2007) have all suggested this can and should be higher. 
# Graham et. al suggests that "researchers using MI should perform many 
# more imputations than previously considered sufficient", and White 
# suggested a lower bound to be 100 * the percent of cases and then to 
# go slightly higher, which here is 28. Graham suggests 20 imputations 
# for 10% to 30% of missing data. The main conclusion of the recent literature
# is, "the number of imputations should be similar to the percentage of 
# cases that are incomplete." Given the computational expense and the above
# literature, plus the small amount of missing data, a value of 40 seems valid
num_imputations <- 20

# Royston and White (2011) and Van Buuren et al. (1999) have all suggested
# that more than 10 cycles are needed for the convergence of the sampling
# distribution of imputed values, but it has also been found that it can be
# satisfactory with just 5-10 (Brand 1999; vanBuuren et al. 2006b). However,
# they also note that while slower, added extra iterations is not a bad thing.
# Van Buuren 2018 says 5-20 iterations is enough to reach convergence. However,
# we ran the well-known method described in "MICE in R" from the Journal of 
# Statistical Software (2011), and found good convergence using just 10 
# iterations. As a precaution, we've upped this to 25.
iterations <- 10

# Simply just set up the methods and predictor matrices, as suggested in Heymans and Eekhout's "Applied Missing Data Analysis"
init <- mice(data, maxit = 0) 
methods <- init$method
predM <- init$predictorMatrix

# For dichotomous variables, use logistic regression predictors, and for
# categorical variables, use polynomial regression
methods[c("Race", "Smoking.History.")] = "polyreg"
methods[c("SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", 
    "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", 
    "Baseline.Plaquenil", "ACEI.ARB", "Anticoagulation.", "IsClade1", 
    "Hospitalized.for.COVID", "Death")] = "logreg" 

# Set all variables to 0 to begin with
predM <- ifelse(predM < 0, 1, 0)

# Variables which will be used for prediction
predictor_vars <- c(
    "Hospitalized.for.COVID", "Age", "Gender", "Race", "SNF", "Diabetes", 
    "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", 
    "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", 
    "ACEI.ARB", "Smoking.History.", "Anticoagulation.", "IsClade1", "Death"
)

# Pick which factors should be involved in imputation. This is a well-known
# issue in multiple imputation. Meng (1994), Rubin (1996), 
# Taylor et al. (2002), and White, Royston, and Wood (2011) advocate 
# including all variables associated with the probability of missingness, 
# along with the variables contained in the dataset, and van Buuren (1999) 
# found that, "Asa a general rule, using all available information yields 
# multiple imputations that have minimal bias and maximal certainty. This 
# principle implies that the number of predictors should be as large as 
# possible."  Enders, Dietz, Montague, and Dixon (2006), Graham (2009), and 
# Jolani, Van Buuren, and Frank (2011), the imputation model should be more
#  general than the analysis model in order to capture more associations 
# between the variables. Finally, it is summed up by Hardt (2012): "the 
# imputation model should include all variables of the analysis, plus those 
# highly correlated with responses or explanatory variables". For this reason,
# we've included all variables
for (predictor_var in predictor_vars) {
    predM[predictor_var, predictor_vars] <- 1
    predM[predictor_var, predictor_var] <- 0
}

# We use multiple imputation using MICE. This is a set of multiple imputations for 
# data that is MAR. For Race and for Smoking history, their probability of being 
# missing is heavily correlated with Hospitalized.for.COVID & HTN (Fisher test has 
# p = 0.007 and 0.01943, respectively). It is difficult to rule out MNAR, though
# MICE still works for MNAR
imputed <- mice(
    data = data, 
    method = methods, 
    predictorMatrix = predM, 
    m = num_imputations, 
    maxit = iterations, 
    seed = random_seed_num
)

### Do univariate analysis
# test_matrices is a list that combines each variable for which we want to do
# univariate analysis. The first term is the name of the category, the second
# term is the default value, and the third is the reference group, or the
# value that we will be comparing against
test_matrices <- list(
    c("Hospitalized.for.COVID", "Yes", "No"),
    c("Gender", "M", "F"),
    c("SNF", "Yes", "No"),
    c("Race", "White", "White"),
    c("Race", "Asian", "White"),  # For example, we use "White" as a ref. group
    c("Race", "African American", "White"),
    c("Race", "Native Hawaiian or Pacific Islander", "White"),
    c("Race", "American Indian or Alaskan Native", "White"),
    c("Race", "Other", "White"),
    c("Diabetes", "Yes", "No"),
    c("HTN", "Yes", "No"),
    c("COPD", "Yes", "No"),
    c("Asthma", "Yes", "No"),
    c("CVD", "Yes", "No"),
    c("CHF", "Yes", "No"),
    c("CKD", "Yes", "No"),
    c("Cancer", "Yes", "No"),
    c("Hx.of.DVT", "Yes", "No"),
    c("Hypothyroid", "Yes", "No"),
    c("Hx.of.MI", "Yes", "No"),
    c("Smoking.History.", "Former", "Never"),
    c("Smoking.History.", "Never", "Never"),
    c("Smoking.History.", "Active", "Never"),
    c("Smoking.History.", "Unknown", "Never"),
    c("Steroids.or.IMT", "Yes", "No"),
    c("Baseline.Plaquenil", "Yes", "No"),
    c("ACEI.ARB", "Yes", "No"),
    c("Anticoagulation.", "Yes", "No")
)

# Function that returns a 2x2 contingency matrix given a data set and a test
# matrix. Computes the contingency table given by the T/F matrix from the
# boolean values of clade and test_matrix[2]
get_matrix <- function(test_matrix, data) {
    if (test_matrix[2] == test_matrix[3]) {
        # In the case when the default and the reference group are the same, 
        # we compare the ref. group against everyone else. This is only to
        # attain a useful p-value and odds ratio. If preferred, we could simply
        # delete the results attained here.
        curr <- matrix(c(
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & IsClade1)),
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & !IsClade1)),
            nrow(subset(data, eval(parse(text = test_matrix[1])) != test_matrix[3] & IsClade1)),
            nrow(subset(data, eval(parse(text = test_matrix[1])) != test_matrix[3] & !IsClade1))), nrow = 2, ncol = 2 )
    } else {
        curr <- matrix(c(
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & IsClade1)),
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & !IsClade1)),
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[3] & IsClade1)),
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[3] & !IsClade1))), nrow = 2, ncol = 2 )
    }
    return(curr)
}

format_p_val <- function(p_val) {
  if (is.na(p_val) || is.nan(p_val) || !is.numeric(p_val)) {
    return("-")
  } else if (p_val < 0.0001) {
    return(paste0(format(round(p_val, 3), nsmall = 3), "****"))
  } else if (p_val < 0.001) {
    return(paste0(format(round(p_val, 3), nsmall = 3), "***"))
  } else if (p_val < 0.01) {
    return(paste0(format(round(p_val, 3), nsmall = 3), "**"))
  } else if (p_val < 0.05) {
    return(paste0(format(round(p_val, 3), nsmall = 3), "*"))
  } else {
    return(paste0(format(round(p_val, 3), nsmall = 3), ""))
  }
}

pct_cols <- c(4, 5, 6, 9, 10, 11, 14, 15, 16)

# Computes all univariate analysis, using Fisher's and Chi-SQ. Also computes 
# N and percentages, and Clopper-Pearson intervals from them. This could also
# use the Barnard's test, as it is been shown to have a higher power, but we
# won't include it since most of the p-values are conclusive without it
print("value,num,total-num,percentage,95-CI-pct-low,95-CI-pct-high,num-clade-1,total-clade-1,pct-clade-1,95-CI-pct-clade-1-low,95-CI-pct-clade-1-high,num-clade-2,total-clade-2,pct-clade-2,95-CI-pct-clade-2-low,95-CI-pct-clade-2-high,odds-ratio,95-CI-odds-ratio-low,95-CI-odds-ratio-high,fisher-p-val,chisq-p-val,num-missing-data,pct-missing-data")
for (test_matrix in test_matrices) {
    # Get the contingency matrix for this term
    curr <- get_matrix(test_matrix, data)

    # Adding all the parts of the matrix together, this is basically what the bottom of the excel spreadsheet is doing when finding percentages of the three groups "all" "clade 1" and "clade 2"
    all_total <- (curr[1][1] + curr[2][1] + curr[3][1] + curr[4][1])
    all_N <- (curr[1][1] + curr[2][1])
    all_ci <- exactci(all_N, all_total, conf.level = 0.95)
    all <- c(all_N, all_total, all_N / all_total, all_ci$conf.int[1][1], all_ci$conf.int[2][1])

    # Creating N and totals for clade 1
    clade_1_total <- curr[1][1] + curr[3][1]
    clade_1_N <- curr[1][1]
    clade_1_CI <- exactci(clade_1_N, clade_1_total, conf.level = 0.95)
    clade_1 <- c(clade_1_N, clade_1_total, clade_1_N / clade_1_total, clade_1_CI$conf.int[1][1], clade_1_CI$conf.int[2][1])

    # Creating N and totals for clade 2
    clade_2_total <- curr[2][1] + curr[4][1]
    clade_2_N <- curr[2][1]
    clade_2_CI <- exactci(clade_2_N, clade_2_total, conf.level = 0.95)
    clade_2 <- c(clade_2_N, clade_2_total, clade_2_N / clade_2_total, clade_2_CI$conf.int[1][1], clade_2_CI$conf.int[2][1])

    # Don't print the standard output because we don't care about this until later
    invisible(capture.output(odds <- oddsratio(curr)))

    # Calculating fisher's and chi-square. The for-loop loops through each 
    # imputated data set and runs the univariate regression. This uses the 
    # median-P approach for univariate analysis on imputed data sets, here:
    # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-017-0404-7
    fisher_p_vals <- c()
    chisq_p_vals <- c()
    if (USE_IMPUTATIONS_FOR_UNIVARIATE) {
        for (i in 1:num_imputations) {
            # Gets the contingency matrix based on the ith imputed data set
            curr <- get_matrix(test_matrix, complete(imputed, i))

            # Uses fisher's mid-P correction, since it tends to be too 
            # conservative. Can find research on this
            fisher <- uncondExact2x2(curr[1][1], 
                curr[1][1]+curr[3][1], 
                curr[2][1], 
                curr[2][1] + curr[4][1], 
                parmtype = "oddsratio", 
                conf.int = FALSE, 
                midp = TRUE
            )
            chisq <- chisq.test(curr)
            fisher_p_vals <- c(fisher_p_vals, fisher$p.value)
            chisq_p_vals <- c(chisq_p_vals, chisq$p.value)
        }
    } else {
        fisher <- uncondExact2x2(curr[1][1], 
            curr[1][1]+curr[3][1], 
            curr[2][1], 
            curr[2][1] + curr[4][1], 
            parmtype = "oddsratio", 
            conf.int = FALSE, 
            midp = TRUE
        )
        chisq <- chisq.test(curr)
        fisher_p_vals <- c(fisher_p_vals, fisher$p.value)
        chisq_p_vals <- c(chisq_p_vals, chisq$p.value)
    }


    row_name <- paste(test_matrix[1], test_matrix[2], sep="-")
    row <- c(row_name, all, clade_1, clade_2, 
        odds$estimate, odds$conf.int[1], odds$conf.int[2], 
        median(fisher_p_vals), median(chisq_p_vals)
    )

    for (i in pct_cols) {
        if (!is.na(row[i])) {
            row[i] = paste0(format(round(100 * as.numeric(row[i]), 1), nsmall = 1), "%")
        }
    }

    for (i in 17:19) {
        if (!is.na(row[i])) {
            row[i] = format(round(as.numeric(row[i]), 3), nsmall = 3)
        }
    }
    
    row[20] = format_p_val(as.numeric(row[20]))
    row[21] = format_p_val(as.numeric(row[21]))

    col_name <- test_matrix[1]
    num_missing <- length(which(is.na(data[col_name])))
    pct_missing <- paste0(format(round(100 * num_missing / nrow(data), 1), nsmall = 1), "%")
    row <- c(row, num_missing, pct_missing)

    print(paste(row, collapse=","))
}

# Perform t-test on age
# In this case, age is not actually normally distrubuted (as can be seen by Shapiro's), but
# the p-value is extremely small, and a large body of literature supports the use of a 
# t-test with enough N, as opposed to the Mann-Whitney U test. However, we have included it
# here, just in case p values differ by a lot. Care must be taken when interpreting this 
# p-value, as the null hypothesis is not the same as a t-test
age_clade_1 <- (data %>% filter(!IsClade1))$Age
age_clade_1 <- (data %>% filter(IsClade1))$Age
t_test <- t.test(age_clade_1, age_clade_1)
wilcox <- wilcox.test(age_clade_1, age_clade_1)
print("Age stats:")
print("total-num,mean,total-clade-1,pct,mean,total-clade-2,pct,mean,t-test-p-value,95-CI-low,95-CI-high,wilcox-p-value")
print(paste(
    nrow(data),
    mean(data$Age),
    length(age_clade_1),
    length(age_clade_1) / nrow(data),
    mean(age_clade_1),
    length(age_clade_1),
    length(age_clade_1) / nrow(data),
    mean(age_clade_1),
    format_p_val(t_test$p.value),
    t_test$conf.int[1],
    t_test$conf.int[2],
    format_p_val(wilcox$p.value),
    sep = ","
))


# Remove the 4 rows that have missing data here. Each of these rows has exactly
# 4 missing values, and they are all from the same 4 patients. These patients have
# only age, gender, and race completed, and therefore probably shouldn't be included
# The paper found here (https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-017-0442-1)
# remarks that studies should use listwise deletion if the % of deleted rows is 
# below 5%. In our case, it is 2.1%, and encompasses data for which imputation makes
# little sense. Furthermore, imputation doesn't really make sense here since there is no
# pattern to the missingness. It is essentially MCAR, for which MI should not be done
if (listwise_deletion) {
    data <- data %>% filter(
        !is.na(Diabetes) 
        & !is.na(HTN) 
        & !is.na(COPD) 
        & !is.na(Asthma) 
        & !is.na(Hx.of.DVT) 
        & !is.na(CKD) 
        & !is.na(Cancer) 
        & !is.na(Hx.of.MI) 
        & !is.na(CVD) 
        & !is.na(CHF) 
        & !is.na(Hypothyroid)
        & !is.na(Baseline.Plaquenil)
    )

    imputed <- mice(
        data = data, 
        method = methods, 
        predictorMatrix = predM, 
        m = num_imputations, 
        maxit = iterations, 
        seed = random_seed_num
    )
}

### Perform LOO calibration
subset <- complete(imputed, sample(1:imputed$m, 1))
subset <- subset[, !names(subset) %in% c("Accession", "Death")]

# The models which we will test LOOCV for
# model_1: age, gender, race,
model_1 <- IsClade1 ~ Age + Gender + Race

# All univariate significance of 0.1 <
model_2 <- IsClade1 ~ Hospitalized.for.COVID + COPD + CVD + Cancer + Hx.of.DVT + Smoking.History. + Steroids.or.IMT + Anticoagulation.

# Fit a model using backwards stepwise AIC. Note that stepwise AIC has been found to be problematic for many
# reasons, and should not be used as anything other than a diagnostic
result <- stepAIC(glm(IsClade1 ~ ., family = binomial(link = logit), data = subset), direction = "backward", trace = FALSE)
terms <- attr(terms(result), "term.labels")
stepwise_model <- as.formula(paste0("IsClade1 ~ ", paste(terms, collapse = " + ")))

# LASSO regression selected. Same for MSE, AUC, and class, exception of Steroids
model_4 <- IsClade1 ~ Cancer + CVD + Steroids.or.IMT + Smoking.History. + Anticoagulation.

# Stepwise with optimizing for ROC_AUC by MDG
model_5 <- IsClade1 ~ Age + Race + Cancer + Smoking.History. + Hospitalized.for.COVID

# Stepwise AUC-RF by MDA
model_6 <- IsClade1 ~ Race + Age + Hospitalized.for.COVID + Hypothyroid + Gender + Cancer + CVD + Steroids.or.IMT + Smoking.History. + Anticoagulation.

# Model removing worst elements
combined_model <- IsClade1 ~ Age + Race + CVD + CKD + Cancer + Smoking.History. + Steroids.or.IMT + Anticoagulation.

models <- c(model_1, combined_model)

# Variables with significance in univariate, stepwise AIC regression, and in LASSO regression:
# Cancer & Smoking

patients <- unique(data$Accession)
predicteds <- matrix(nrow = length(patients), ncol = length(models) + 1)

# For each patient, we will test each model on each imputed data set. For each model, we
# will average the outputs across imputed data sets, and report that number. That will
# also be used for the MSE at the end
for (n in 1:length(patients)) {
    patient <- patients[n]
    if ((data %>% filter(Accession == patient))$IsClade1[1]) {
        is_clade_1 <- 1
    } else {
        is_clade_1 <- 0
    }

    averages <- vector(length = length(models))
    index_num <- 1
    for (model in models) {
        sum <- 0
        for (i in 1:num_imputations) {
            imputed_data <- complete(imputed, i)
            train_data <- imputed_data %>% filter(Accession != patient)
            test_data <- imputed_data %>% filter(Accession == patient)
            output <- glm(model, family = binomial(link=logit), data = train_data)
            response <- predict(output, newdata = test_data, type="response")
            sum <- sum + response
        }

        avg <- sum / num_imputations
        averages[index_num] <- avg
        index_num <- index_num + 1
    }

    predicteds[n,1] <- is_clade_1
    for (i in 1:length(models)) {
        predicteds[n,i+1] <- averages[i]
    }

    print(paste(paste(patient, is_clade_1, sep=","), paste(averages, collapse = ","), sep=","))
}

# Create graphs for ROC and PR data
par(mfrow=c(2,3))
method <- "empirical"
negref <- "0"
roc_1 <- rocit(score = predicteds[,2], class = predicteds[,1], method = method, negref= negref)
roc_2 <- rocit(score = predicteds[,3], class = predicteds[,1], method = method, negref= negref)
measure_1 <- measureit(score = predicteds[,2], class = predicteds[,1], measure = c("ACC", "SENS", "SPEC", "PREC"), negref= negref)
measure_2 <- measureit(score = predicteds[,3], class = predicteds[,1], measure = c("ACC", "SENS", "SPEC", "PREC"), negref= negref)

plot(roc_1, col = c(1,"black"), legend = FALSE, YIndex = FALSE)
lines(roc_2$TPR ~ roc_2$FPR, col = c(2,"red"), lwd = 2)
legend("bottomright", col = c(1,2), c("Model 1", "Model 2"), lwd = 2)

plot(measure_1$ACC ~ measure_1$Cutoff, col = c(1,"black"), legend = FALSE, YIndex = FALSE, type = "l", ylim = c(0.35, 0.65), xlab = "Cutoff", ylab = "Accuracy")
lines(measure_2$ACC ~ measure_2$Cutoff, type = "l", col = c(2,"red"), lwd = 2)
legend("bottomright", title = "Accuracy by Cutoff", col = c(1,2), c("Model 1", "Model 2"), lwd = 2)

plot(measure_1$SENS ~ measure_1$Cutoff, col = c(1,"black"), legend = FALSE, YIndex = FALSE, type = "l", xlab = "Cutoff", ylab = "Sensitivity")
lines(measure_2$SENS ~ measure_2$Cutoff, type = "l", col = c(2,"red"), lwd = 2)
legend("bottomright", title = "Sensitivity by Cutoff", col = c(1,2), c("Model 1", "Model 2"), lwd = 2)

plot(measure_1$SPEC ~ measure_1$Cutoff, col = c(1,"black"), legend = FALSE, YIndex = FALSE, type = "l", xlab = "Cutoff", ylab = "Specificity")
lines(measure_2$SPEC ~ measure_2$Cutoff, type = "l", col = c(2,"red"), lwd = 2)
legend("bottomright", title = "Specification by Cutoff", col = c(1,2), c("Model 1", "Model 2"), lwd = 2)

plot(measure_1$PREC ~ measure_1$Cutoff, col = c(1,"black"), legend = FALSE, YIndex = FALSE, type = "l", ylim = c(0, 1), xlab = "Cutoff", ylab = "Precision")
lines(measure_2$PREC ~ measure_2$Cutoff, type = "l", col = c(2,"red"), lwd = 2)
legend("bottomright", title = "Precision by Cutoff", col = c(1,2), c("Model 1", "Model 2"), lwd = 2)

plot(measure_1$PREC ~ measure_1$SENS, col = c(1,"black"), legend = FALSE, type = "l", YIndex = FALSE, ylim = c(0, 1), xlab = "Recall", ylab = "Precision")
lines(measure_2$PREC ~ measure_2$SENS, type = "l", col = c(2,"red"), lwd = 2)
legend("bottomright", title = "Precision/Recall", col = c(1,2), c("Model 1", "Model 2"), lwd = 2)

# Pool uses Rubin's Rules to pool models built on a matrix of imputed data sets 
# (basically builds a model for every single imputed dataset, i.e. all 25, and 
# then uses this set of rules to pool to coefficents into an average)
# glm is the logistic regression function. adjusted risk ratio is the e^coefficient value provided

# Function that returns the data associated with a coefficient term
get_data <- function(pooled, term) {
    summ <- summary(pooled)
    for (i in 1:length(summ$term)) {
        if ((summ$term)[i] == term) {
            or <- format(round(exp((summ$estimate)[i]), 2), nsmall = 2)
            # error <- format(round((summ$std.error)[i], 3), nsmall = 3)
            p_val <- format_p_val((summ$p.value)[i])
            v <- c(or, p_val)
            return(v)
        }
    }
    return(c("",""))
}

# All variables included
all_adjusted_RR <- pool(with(
    imputed, 
    glm(IsClade1 ~ Hospitalized.for.COVID + Age + Gender + SNF + Race + Diabetes + HTN 
        + COPD + Asthma + CVD + CHF + CKD + Cancer + Hypothyroid + Hx.of.DVT 
        + Hx.of.MI + Smoking.History. + Steroids.or.IMT + Baseline.Plaquenil 
        + ACEI.ARB + Anticoagulation., family = binomial(link=logit))
    ))

# Only using the model from model_1
model_1_adjusted_RR <- pool(with(
    imputed, 
    glm(IsClade1 ~ Age + Gender + Race, family = binomial(link=logit))
    ))

# Only using the model from model_2
model_2_adjusted_RR <- pool(with(
    imputed, 
    glm(IsClade1 ~ Hospitalized.for.COVID + COPD + CVD + Cancer + 
    Hx.of.DVT + Smoking.History. + Steroids.or.IMT + Anticoagulation., family = binomial(link=logit))
    ))

# Stepwise AIC
model_3_adjusted_RR <- pool(with(
    imputed, 
    glm(IsClade1 ~ Diabetes + CKD + Cancer + CVD + Smoking.History. + Steroids.or.IMT + Anticoagulation., family = binomial(link=logit))
    ))

# Using LASSO regression
model_4_adjusted_RR <- pool(with(
    imputed, 
    glm(IsClade1 ~ Cancer + CVD + Steroids.or.IMT + Smoking.History. + Anticoagulation., family = binomial(link=logit))
    ))

# Combined model
model_5_adjusted_RR <- pool(with(
    imputed,
    glm(IsClade1 ~ Age + Race + CVD + CKD + Cancer + Smoking.History. + Steroids.or.IMT + Anticoagulation., family = binomial(link=logit))
))

# Create a list of all the terms
terms <- c("(Intercept)", "Age")
for (test_matrix in test_matrices) {
    terms <- c(terms, paste(test_matrix[1], test_matrix[2], sep=""))
}

# Print out the adjusted risk ratios, standard error, and p-values of the models
print("term,all-var-adjusted-ratio,all-var-std-error,all-var-p-value,model-1-adjusted-ratio,model-1-std-error,model-1-p-value,model-2-adjusted-ratio,model-2-std-error,model-2-p-value,model-3-adjusted-ratio,model-3-std-error,model-3-p-value,model-4-adjusted-ratio,model-4-std-error,model-4-p-value,model-5-adjusted-ratio,model-5-std-error,model-5-p-value")
for (term in terms) {
    aaRR <- paste(get_data(all_adjusted_RR, term), collapse=",")
    m1RR <- paste(get_data(model_1_adjusted_RR, term), collapse=",")
    m2RR <- paste(get_data(model_2_adjusted_RR, term), collapse=",")
    m3RR <- paste(get_data(model_3_adjusted_RR, term), collapse=",")
    m4RR <- paste(get_data(model_4_adjusted_RR, term), collapse=",")
    m5RR <- paste(get_data(model_5_adjusted_RR, term), collapse=",")
    print(paste(term, aaRR, m1RR, m2RR, m3RR, m4RR, m5RR, sep=","))
}

# Important to note that the VIF for all variables in all 3 models were under 2.
# As a result, we can conclude that there was not a high degree of collinearity


## Bibliogrpahy
# Allison, PD. (2002). Missing data. Thousand Oaks, CA: Sage.
# Brand JPL (1999). Development, Implementation and Evaluation of Multiple Imputation Strategies for the Statistical Analysis of Incomplete Data Sets. Ph.D. thesis, Erasmus University, Rotterdam.
# Bodner, Todd E. (2008) “What improves with increased missing data imputations?” Structural Equation Modeling: A Multidisciplinary Journal 15: 651-675.
# Moons, KG, Donders, RA, Stijnen, T, & Harrell, FE, Jr. (2006). Using the outcome for imputation of missing predictor values was preferred. Journal of Clinical Epidemiology, 59, 1092–1101.
# Royston, P, & White, IR. (2011). Multiple imputation by chained equations (MICE): implementation in Stata. Journal of Statistical Software, 45(4), 1–20.
# White, IR, Royston, P, & Wood, AM. (2011). Multiple imputation using chained equations: Issues and guidance for practice. Statistics in Medicine, 30, 377–399
# Graham, JW, Olchowski, AE, & Gilreath, TD. (2007). How many imputations are really needed? Some practical clarifications of multiple imputation theory. Prevention Science, 8, 206–213.
# Van Buuren, S, Boshuizen, HC, & Knook, DL. (1999). Multiple imputation of missing blood pressure covariates in survival analysis. Statistics in Medicine, 18, 681–694.
# van Buuren S, Brand JPL, Groothuis-Oudshoorn CGM, Rubin DB (2006b). “Fully Conditional Specification in Multivariate Imputation.” Journal of Statistical Computation and Simulation, 76(12), 1049–1064.
# Van Buuren, S. 2018. Flexible Imputation of Missing Data. Second Edition. Boca Raton, FL: Chapman & Hall/CRC.
