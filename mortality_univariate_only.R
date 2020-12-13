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

# Boolean for deleting the 4 rows that are 90% "Unknown"
listwise_deletion <- FALSE

# Load the CSV
file_path <- "~/Downloads/COVIDCSV - main.csv"
data <- read.csv(
    file = file_path, 
    header = TRUE, 
    na.strings = c("Unknown"), 
    strip.white = TRUE
)

# Set factor variables
data <- data %>%
    mutate(
        Accession = as.factor(Accession),
        Gender = as.factor(Gender),
        Race = as.factor(Race),
        SNF = as.factor(SNF),
        Death = as.factor(Death),
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
        Anticoagulation. = as.factor(Anticoagulation.)
    )

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
}

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
num_imputations <- 40

# Royston and White (2011) and Van Buuren et al. (1999) have all suggested
# that more than 10 cycles are needed for the convergence of the sampling
# distribution of imputed values, but it has also been found that it can be
# satisfactory with just 5-10 (Brand 1999; vanBuuren et al. 2006b). However,
# they also note that while slower, added extra iterations is not a bad thing.
# Van Buuren 2018 says 5-20 iterations is enough to reach convergence. However,
# we ran the well-known method described in "MICE in R" from the Journal of 
# Statistical Software (2011), and found good convergence using just 10 
# iterations. As a precaution, we've upped this to 25.
iterations <- 25

# Simply just set up the methods and predictor matrices, as suggested in Heymans and Eekhout's "Applied Missing Data Analysis"
init <- mice(data, maxit = 0) 
methods <- init$method
predM <- init$predictorMatrix

# For dichotomous variables, use logistic regression predictors, and for
# categorical variables, use polynomial regression
methods[c("Race", "Smoking.History.")] = "polyreg"
methods[c("SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", 
    "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", 
    "Baseline.Plaquenil", "ACEI.ARB", "Anticoagulation.")] = "logreg" 

# Set all variables to 0 to begin with
predM <- ifelse(predM < 0, 1, 0)

# Variables which will be used for prediction
predictor_vars <- c(
    "Hospitalized.for.COVID", "Age", "Gender", "Race", "SNF", "Diabetes", 
    "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", 
    "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", 
    "ACEI.ARB", "Smoking.History.", "Anticoagulation."
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
# missing is heavily correlated with hospitalized.for.COVID & HTN (Fisher test has 
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
    c("Smoking.History.", "Unclear", "Never"),
    c("Steroids.or.IMT", "Yes", "No"),
    c("Baseline.Plaquenil", "Yes", "No"),
    c("ACEI.ARB", "Yes", "No"),
    c("Anticoagulation.", "Yes", "No")
)

# Function that returns a 2x2 contingency matrix given a data set and a test
# matrix. Computes the contingency table given by the T/F matrix from the
# boolean values of Death and test_matrix[2]
get_matrix <- function(test_matrix, data) {
    if (test_matrix[2] == test_matrix[3]) {
        # In the case when the default and the reference group are the same, 
        # we compare the ref. group against everyone else. This is only to
        # attain a useful p-value and odds ratio. If preferred, we could simply
        # delete the results attained here.
        curr <- matrix(c(
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & Death == 'Yes')),
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & Death == 'No')),
            nrow(subset(data, eval(parse(text = test_matrix[1])) != test_matrix[3] & Death == 'Yes')),
            nrow(subset(data, eval(parse(text = test_matrix[1])) != test_matrix[3] & Death == 'No'))), nrow = 2, ncol = 2 )
    } else {
        curr <- matrix(c(
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & Death == 'Yes')),
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & Death == 'No')),
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[3] & Death == 'Yes')),
            nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[3] & Death == 'No'))), nrow = 2, ncol = 2 )
    }
    return(curr)
}

# Computes all univariate analysis, using Fisher's and Chi-SQ. Also computes 
# N and percentages, and Clopper-Pearson intervals from them. This could also
# use the Barnard's test, as it is been shown to have a higher power, but we
# won't include it since most of the p-values are conclusive without it
print("value,num,total-num,percentage,95-CI-pct-low,95-CI-pct-high,num-dead,total-dead,pct-dead,95-CI-pct-dead-low,95-CI-pct-dead-high,num-alive,total-alive,pct-alive,95-CI-pct-alive-low,95-CI-pct-alive-high,odds-ratio,95-CI-odds-ratio-low,95-CI-odds-ratio-high,fisher-p-val,chisq-p-val")
for (test_matrix in test_matrices) {
    # Get the contingency matrix for this term
    curr <- get_matrix(test_matrix, data)

    # Adding all the parts of the matrix together, this is basically what the bottom of the excel spreadsheet is doing when finding percentages of the three groups "all" "dead" and "alive"
    all_total <- (curr[1][1] + curr[2][1] + curr[3][1] + curr[4][1])
    all_N <- (curr[1][1] + curr[2][1])
    all_ci <- exactci(all_N, all_total, conf.level = 0.95)
    all <- c(all_N, all_total, all_N / all_total, all_ci$conf.int[1][1], all_ci$conf.int[2][1])

    # Creating N and totals for mortality
    dead_total <- curr[1][1] + curr[3][1]
    dead_N <- curr[1][1]
    dead_ci <- exactci(dead_N, dead_total, conf.level = 0.95)
    dead <- c(dead_N, dead_total, dead_N / dead_total, dead_ci$conf.int[1][1], dead_ci$conf.int[2][1])

    # Creating N and totals for alive
    alive_total <- curr[2][1] + curr[4][1]
    alive_N <- curr[2][1]
    alive_ci <- exactci(alive_N, alive_total, conf.level = 0.95)
    alive <- c(alive_N, alive_total, alive_N / alive_total, alive_ci$conf.int[1][1], alive_ci$conf.int[2][1])

    # Don't print the standard output because we don't care about this until later
    invisible(capture.output(odds <- oddsratio(curr)))

    # Calculating fisher's and chi-square. The for-loop loops through each 
    # imputated data set and runs the univariate regression. This uses the 
    # median-P approach for univariate analysis on imputed data sets, here:
    # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-017-0404-7
    fisher_p_vals <- c()
    chisq_p_vals <- c()
    for (i in 1:num_imputations) {
        # Gets the contingency matrix based on the ith imputed data set
        curr <- get_matrix(test_matrix, complete(imputed, i))

        # Uses fisher's mid-P correction, since it tends to be too 
        # conservative. Can find research on this
        fisher <- uncondExact2x2(curr[1][1], 
            curr[1][1]+curr[2][1], 
            curr[3][1], 
            curr[3][1] + curr[4][1], 
            parmtype = "oddsratio", 
            conf.int = FALSE, 
            midp = TRUE
        )
        chisq <- chisq.test(curr)
        fisher_p_vals <- c(fisher_p_vals, fisher$p.value)
        chisq_p_vals <- c(chisq_p_vals, chisq$p.value)
    }

    row_name <- paste(test_matrix[1], test_matrix[2], sep="-")
    row <- c(row_name, all, dead, alive, 
        odds$estimate, odds$conf.int[1], odds$conf.int[2], 
        median(fisher_p_vals), median(chisq_p_vals)
    )
    print(paste(row, collapse=","))
}

# Perform t-test on age
# In this case, age is not actually normally distrubuted (as can be seen by Shapiro's), but
# the p-value is extremely small, and a large body of literature supports the use of a 
# t-test with enough N, as opposed to the Mann-Whitney U test. However, we have included it
# here, just in case p values differ by a lot. Care must be taken when interpreting this 
# p-value, as the null hypothesis is not the same as a t-test
age_alive <- (data %>% filter(Death == "No"))$Age
age_dead <- (data %>% filter(Death == "Yes"))$Age
t_test <- t.test(age_alive, age_dead)
wilcox <- wilcox.test(age_alive, age_dead)
print("Age stats:")
print("total-num,mean,total-dead,pct,mean,total-alive,pct,mean,t-test-p-value,95-CI-low,95-CI-high,wilcox-p-value")
print(paste(
    nrow(data),
    mean(data$Age),
    length(age_dead),
    length(age_dead) / nrow(data),
    mean(age_dead),
    length(age_alive),
    length(age_alive) / nrow(data),
    mean(age_alive),
    t_test$p.value,
    t_test$conf.int[1],
    t_test$conf.int[2],
    wilcox$p.value,
    sep = ","
))
