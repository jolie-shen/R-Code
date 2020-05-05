library(dplyr)
library(mice)
library(PropCIs)
library(questionr)
library(exact2x2)
library(regclass)
library(fmsb)

covid <- read.csv(file = "~/Downloads/COVIDCSV - main.csv", header = TRUE, na.strings=c("Unknown"))

covid <- covid %>%
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
        Anticoagulation. = as.factor(Anticoagulation.)
    )

### Perform multiple imputation
set.seed(3249)
num_imputations <- 25
init <- mice(covid, maxit=0) 
meth <- init$method
predM <- init$predictorMatrix

meth[c("SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Anticoagulation.")]="logreg" 
meth[c("Race", "Smoking.History.")]="polyreg"

predM <- ifelse(predM < 0, 1, 0)

#Pick which factors should be involved in imputation, this could be changed based on what we think best predicts each (i.e. does having a smoking history reallly predict anticoagulation status?). 
#However, conditionoing on alll other data is often reasonable to small data sets contain up to 20-30 variabels. as a genral rule, usign eveyr bit of availabel information yields multiple imputatuons that have minimal bias and maximal efficeincy
#Collins, L. M., J. L. Schafer, and C. M. Kam. 2001. “A Comparison of Inclusive and Restrictive Strategies in Modern Missing Data Procedures.” Psychological Methods 6 (3): 330–51.
predM["Race", c("Age", "Gender", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["Diabetes", c("Age", "Gender", "Race", "SNF", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["HTN", c("Age", "Gender", "Race", "SNF", "Diabetes", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["COPD", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["Asthma", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["Hx.of.DVT", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["CKD", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["Cancer", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["Hx.of.MI", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["CVD", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["CHF", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["Hypothyroid", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["Steroids.or.IMT", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["Baseline.Plaquenil", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "ACEI.ARB", "Smoking.History.", "Anticoagulation.")] = 1
predM["ACEI.ARB", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "Smoking.History.", "Anticoagulation.")] = 1
predM["Smoking.History.", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Anticoagulation.")] = 1
predM["Anticoagulation.", c("Age", "Gender", "Race", "SNF", "Diabetes", "HTN", "COPD", "Asthma", "Hx.of.DVT", "CKD", "Cancer", "Hx.of.MI", "CVD", "CHF", "Hypothyroid", "Steroids.or.IMT", "Baseline.Plaquenil", "ACEI.ARB", "Smoking.History.")] = 1


imputed <- mice(covid, method=meth, predictorMatrix=predM, m=num_imputations, maxit=50)



### Do univariate analysis

# Using median-P from https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-017-0404-7
# c(name of category, default value, reference value) where reference value is what we are comparing against (i.e. white and never for race/smoker)
test_matrices <- list(
c("Gender", "M", "F"),
c("SNF", "Yes", "No"),
c("Race", "White", "White"),
c("Race", "Asian", "White"),
c("Race", "African American", "White"),
c("Race", "American Indian or Alaskan Native", "White"),
c("Race", "Native Hawaiian or Pacific Islander", "White"),
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
c("Steroids.or.IMT", "Yes", "No"),
c("Baseline.Plaquenil", "Yes", "No"),
c("ACEI.ARB", "Yes", "No"),
c("Anticoagulation.", "Yes", "No"))


##prints out results for all univariate analysis, then constructing the matrix for univariate analysis (i.e. fishers, chi-square, etc.)

print(paste("value,num,total-num,percentage,95-CI-pct-low,95-CI-pct-high,num-hospitalized,total-hospitalized,pct-hospitalized,95-CI-pct-hospitalized-low,95-CI-pct-hospitalized-high,num-non-hospitalized,total-non-hospitalized,pct-non-hospitalized,95-CI-pct-non-hospitalized-low,95-CI-pct-non-hospitalized-high,odds-ratio,95-CI-odds-ratio-low,95-CI-odds-ratio-high,fisher-p-val,chisq-p-val", sep=","))
for (test_matrix in test_matrices) {
    #the if statement refers to when value of [2] = [3], which means that 'white' = 'white' and 'never' = 'never'. these are unique cases beause they have no reference group (i.e. are the reference group). Instead, we deided to compare white vs. everybody who is not white. This portion could be deleted and the OR/CI would just be reported as blanks
    if (test_matrix[2] == test_matrix[3]) {
        curr <- matrix(c(
            nrow(subset(covid, eval(parse(text = test_matrix[1])) == test_matrix[2] & Hospitalized.for.COVID == 'Yes')),
            nrow(subset(covid, eval(parse(text = test_matrix[1])) == test_matrix[2] & Hospitalized.for.COVID == 'No')),
            nrow(subset(covid, eval(parse(text = test_matrix[1])) != test_matrix[3] & Hospitalized.for.COVID == 'Yes')),
            nrow(subset(covid, eval(parse(text = test_matrix[1])) != test_matrix[3] & Hospitalized.for.COVID == 'No'))), nrow = 2, ncol = 2 )
    #this is the normal case
    } else {
        curr <- matrix(c(
            nrow(subset(covid, eval(parse(text = test_matrix[1])) == test_matrix[2] & Hospitalized.for.COVID == 'Yes')),
            nrow(subset(covid, eval(parse(text = test_matrix[1])) == test_matrix[2] & Hospitalized.for.COVID == 'No')),
            nrow(subset(covid, eval(parse(text = test_matrix[1])) == test_matrix[3] & Hospitalized.for.COVID == 'Yes')),
            nrow(subset(covid, eval(parse(text = test_matrix[1])) == test_matrix[3] & Hospitalized.for.COVID == 'No'))), nrow = 2, ncol = 2 )
    }

    #adding all the parts of the matrix together, this is basically what the bottom of the excel spreadsheet is doing when finding percentages of the three groups "all" "hospitalized" and "non-hospitalized"
    all_total <- (curr[1][1] + curr[2][1] + curr[3][1] + curr[4][1])
    all <- (curr[1][1] + curr[2][1])
    all_ci <- exactci(all, all_total, conf.level = 0.95)
    all_str <- paste(all, all_total, all / all_total, all_ci$conf.int[1][1], all_ci$conf.int[2][1], sep=",")

    hosp_total <- curr[1][1] + curr[3][1]
    hosp <- curr[1][1]
    hosp_ci <- exactci(hosp, hosp_total, conf.level = 0.95)
    hosp_str <- paste(hosp, hosp_total, hosp / hosp_total, hosp_ci$conf.int[1][1], hosp_ci$conf.int[2][1], sep=",")

    non_total <- curr[2][1] + curr[4][1]
    non <- curr[2][1]
    non_ci <- exactci(non, non_total, conf.level = 0.95)
    non_str <- paste(non, non_total, non / non_total, non_ci$conf.int[1][1], non_ci$conf.int[2][1], sep=",")

    #don't print the standard output because we don't care about this until later
    invisible(capture.output(odds <- oddsratio(curr)))

    #calculating fisher's and chi-square. The for-loop loops through each imputated data set and runs the univariate regression
    fisher_p_vals <- c()
    chisq_p_vals <- c()
    for (i in 1:num_imputations) {
        data <- complete(imputed, i)
        if (test_matrix[2] == test_matrix[3]) {
            curr <- matrix(c(
                nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & Hospitalized.for.COVID == 'Yes')),
                nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & Hospitalized.for.COVID == 'No')),
                nrow(subset(data, eval(parse(text = test_matrix[1])) != test_matrix[3] & Hospitalized.for.COVID == 'Yes')),
                nrow(subset(data, eval(parse(text = test_matrix[1])) != test_matrix[3] & Hospitalized.for.COVID == 'No'))), nrow = 2, ncol = 2 )
        } else {
            curr <- matrix(c(
                nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & Hospitalized.for.COVID == 'Yes')),
                nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[2] & Hospitalized.for.COVID == 'No')),
                nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[3] & Hospitalized.for.COVID == 'Yes')),
                nrow(subset(data, eval(parse(text = test_matrix[1])) == test_matrix[3] & Hospitalized.for.COVID == 'No'))), nrow = 2, ncol = 2 )
        }

        fisher <- uncondExact2x2(curr[1][1], curr[1][1]+curr[2][1], curr[3][1], curr[3][1] + curr[4][1], parmtype = "oddsratio", conf.int = FALSE, midp = TRUE)
        chisq <- chisq.test(curr)
        fisher_p_vals <- c(fisher_p_vals, fisher$p.value)
        chisq_p_vals <- c(chisq_p_vals, chisq$p.value)
    }

    print(paste(paste(test_matrix[1], test_matrix[2], sep="-"), all_str, hosp_str, non_str, odds$estimate, odds$conf.int[1], odds$conf.int[2], median(fisher_p_vals), median(chisq_p_vals), sep=","))
}



### Perform LOO calibration
#model_1: age, gender, race,
#model_2: all variables with <=0.1 significance in the univariate analysis (Age, SNF, HTN, DVT, MI, Cancer, CVD, CKD, CHF, Steroids/IMT, ACEI/ARB, Anticoagulation,  )

model_1 <- as.formula("Hospitalized.for.COVID ~ Age + Gender + Race")
model_2 <- as.formula("Hospitalized.for.COVID ~ Age + SNF + HTN + Hx.of.DVT + Hx.of.MI + CKD + Cancer + CVD + CHF + Steroids.or.IMT + ACEI.ARB + Anticoagulation.")

#This goes through each patient and creates 2 datasets with each imputed dataset made previously: one with only that patient (n=1) and one without that patient (n=189). 
patients <- unique(covid$Accession)
mse_1 = 0
mse_2 = 0
for (patient in patients) {
    curr_row <- covid %>% filter(Accession == patient)
    sum_1 = 0
    sum_2 = 0
    for (i in 1:num_imputations) {
        data_set <- complete(imputed, i)
        train_data <- data_set %>% filter(Accession != patient)
        test_data <- data_set %>% filter(Accession == patient)
    #two outputs that fit the data to model_1 and model_2, then predicts the outcome based on the one patient. sum_1 are the predicted outcomes
        output_1 <- glm(model_1, family = binomial, data = train_data)
        response_1 <- predict(output_1, newdata = test_data, type="response")
        sum_1 = sum_1 + response_1

        output_2 <- glm(model_2, family = binomial, data = train_data)
        response_2 <- predict(output_2, newdata = test_data, type="response")
        sum_2 = sum_2 + response_2
    }
    #then we averaged sum_1 over all the imputations and printed it out
    avg_1 <- sum_1 / num_imputations
    avg_2 <- sum_2 / num_imputations
    if (curr_row$Hospitalized.for.COVID[1] == "Yes") {
        print(paste(patient, 1, avg_1, avg_2, sep=","))
        mse_1 = mse_1 + ((1 - avg_1) * (1 - avg_1))
        mse_2 = mse_2 + ((1 - avg_2) * (1 - avg_2))
    } else {
        print(paste(patient, 0, avg_1, avg_2, sep=","))
        mse_1 = mse_1 + (avg_1 * avg_1)
        mse_2 = mse_2 + (avg_2 * avg_2)
    }
}

#gives the Mean Square Error for model_1 and model_2
print(paste("MSE_1: ", mse_1 / length(patients)))
print(paste("MSE_2: ", mse_2 / length(patients)))


#pool uses Rubin's Rules to pool models built on a matrix of imputed data sets (basically builds a model for every single imputed dataset, i.e. all 25, and then uses this set of rules to pool to coefficents into an average)
#glm is the logistic regression function. adjusted risk ratio is the e^coefficient value provided

adjusted_risk_ratios <- pool(with(imputed, glm(Hospitalized.for.COVID ~ Age + SNF + HTN + Hx.of.DVT + Hx.of.MI + CKD + Cancer + CVD + CHF + Steroids.or.IMT + ACEI.ARB + Anticoagulation., family = binomial)))
print(adjusted_risk_ratios)
