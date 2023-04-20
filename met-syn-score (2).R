library(tidyverse)
library(readr)
library(psych)
library(moments)
library(reshape2)
library("readr")
library("janitor")
library(car)


#-----------Metabolic syndrome score----------------------#

#------------------start here------------------####


setwd("C:/Users/c3371138/Dropbox/activate-food")
act_food = read.csv("act-food-pccv.csv", stringsAsFactors = TRUE, encoding = "UTF-8")



# ------------- NOTE ----------------------------- #

# They have only used medication free adults in their calculation of M


#1) Censor: set all biomarker values for each individual case that are below the clinical threshold to the clinical threshold
#2) Center: subtract relevant clinical threshold from each biomarker value (optional: winsorize to remove any outliers)

# I have winsorized the data, should be relatively easy to do posthoc though. 

#input the relevent clinical thresholds 

waist_f = 80
waist_m = 94
trig = 1.7
hdl_m = 1.3
hdl_f = 1
sBP = 130
dBP = 85
BlGl = 5.6

#work out what each of them are in the data 
colnames(act_food)

new_cols = c("visit_age", 
             "sex",
             "waist_average" ,
             "TRIGLY",
             "sbp_average.x",
             "dbp_average.x",
             "blood_glucose",
             "HDL_mmoll", 
             "record_id")

met_data <- act_food[, new_cols] %>% 
  as_tibble() %>% 
  na.omit() # NAs bugger everything

#we need to work out whether 1 or 2 is female

met_data %>% na.omit(waist_average) %>% group_by(sex) %>% summarise(mean = mean(waist_average))

#looks like it's male = 1, female = 2

summary(met_data)

#Censor  the data 

met_data = met_data %>% 
  mutate(waist_average = if_else(sex == 1, waist_average - waist_m,waist_average),
         waist_average = if_else(sex == 2, waist_average - waist_f, waist_average),
         TRIGLY = TRIGLY-trig, 
         HDL_mmoll = if_else(sex == 1,HDL_mmoll - hdl_m,HDL_mmoll),
         HDL_mmoll = if_else(sex == 2,HDL_mmoll - hdl_f, HDL_mmoll),
         sbp_average.x = sbp_average.x - sBP,
         dbp_average.x = dbp_average.x - dBP,
         blood_glucose = blood_glucose - BlGl)
summary(met_data)

#set all less than threshold data to zero and inverse hdl

met_data = met_data %>% mutate_at(vars(-HDL_mmoll), funs(ifelse(. < 0, as.numeric(0), .))) %>% 
  mutate(HDL_mmoll = abs(ifelse(HDL_mmoll > 0, as.numeric(0), HDL_mmoll)))



summary(met_data)


#3) Standardize: divide each censored and centered biomarker value by
#the matching standard deviation shown in Table A.1

met_data[,3:8] = scale(met_data[,3:8], center = FALSE, scale = apply(met_data[,3:8], 2, sd, na.rm = TRUE))

#check sd = 1
apply(met_data[,3:8], 2, sd)



# 4) Orthogonalize: post-multiply the normalized data matrix by the component loading matrix shown in Table A.2
# post-multiply A by B = AB
# pre-multiply A by B = BA


head(met_data)

pca_met <- principal(met_data[,3:8], rotate = "none")

#check the eigenvalues
pca_met$values

#use a 3 factor

pca_met <- principal(met_data[,3:8], rotate = "varimax", nfactors = 3)
pca_met$loadings


#loadings = scale(pca_met$loadings[1:6,]) 
#loadings
loadings = pca_met$loadings[1:6,]

#check the loadings
apply(loadings,2,sd)

#remove the attributes
attr(loadings, "scaled:scale") <- NULL
loadings

#calculate scores using the scaled loadings
scores = scale(met_data[,3:8]) %*% loadings 
head(scores)

# check the correlation with the psych score output 
cor(scores, pca_met$scores)

ncol(met_data[,3:8])
nrow(loadings)


# In the methods section the explanation is much simpler. 
# 4) Principal Component Analysis was used to extract six 
# uncorrelated component scores which were standardized to a standard deviation of one;

#standardise the scores
scores[,1:3] = scale(scores[,1:3], center = FALSE, scale = apply(scores[,1:3], 2, sd, na.rm = TRUE))

#check that they have an SD of 0
apply(scores,2,sd)

#6) Square: square each individual value

square = function(x){x^2}

square_scores = apply(scores, c(1,2), square)

head(square_scores)

#7) Sum: take the sum of the squared  values for each individual row

row_sum = rowSums(square_scores)

#8) Square root: take the square root of the sum of the squared values

MetSSS = sqrt(row_sum)

met_data = cbind(met_data, MetSSS)

head(met_data)

#---------- join met_data to act_food -=-------- ####

gaucci = met_data %>% 
  select(record_id,MetSSS) %>% 
  left_join(.,act_food)

#---------- MetSSS split samples -------------------- ####

# ----------- women -----------# 

met_f = met_data %>% filter(sex == 2)

#use a 3 factor

pca_met <- principal(met_f[,3:8], rotate = "varimax", nfactors = 3)
pca_met$loadings

rm(met_f)


#loadings = scale(pca_met$loadings[1:6,]) 
#loadings
loadings = pca_met$loadings[1:6,]

#check the loadings
apply(loadings,2,sd)

#remove the attributes
attr(loadings, "scaled:scale") <- NULL
loadings

#calculate scores using the scaled loadings
scores = scale(met_data[,3:8]) %*% loadings 
head(scores)


#standardise the scores
scores[,1:3] = scale(scores[,1:3], center = FALSE, scale = apply(scores[,1:3], 2, sd, na.rm = TRUE))

#check that they have an SD of 0
apply(scores,2,sd)

#6) Square: square each individual value

square_scores = apply(scores, c(1,2), square)

head(square_scores)

#7) Sum: take the sum of the squared  values for each individual row

row_sum = rowSums(square_scores)

#8) Square root: take the square root of the sum of the squared values

MetSSS_f = sqrt(row_sum)

met_data = cbind(met_data, MetSSS_f)

# ------- and now for men -------------#

met_m = met_data %>% filter(sex == 1)

#use a 3 factor

pca_met <- principal(met_m[,3:8], rotate = "varimax", nfactors = 3)
pca_met$loadings

rm(met_m)


#loadings = scale(pca_met$loadings[1:6,]) 
#loadings
loadings = pca_met$loadings[1:6,]

#check the loadings
apply(loadings,2,sd)

#remove the attributes
attr(loadings, "scaled:scale") <- NULL
loadings

#calculate scores using the scaled loadings
scores = scale(met_data[,3:8]) %*% loadings 
head(scores)


#standardise the scores
scores[,1:3] = scale(scores[,1:3], center = FALSE, scale = apply(scores[,1:3], 2, sd, na.rm = TRUE))

#check that they have an SD of 0
apply(scores,2,sd)

#6) Square: square each individual value

square_scores = apply(scores, c(1,2), square)

head(square_scores)

#7) Sum: take the sum of the squared  values for each individual row

row_sum = rowSums(square_scores)

#8) Square root: take the square root of the sum of the squared values

MetSSS_m = sqrt(row_sum)

met_data = cbind(met_data, MetSSS_m)

# ------------- old ---------------- #

avg_age = mean(met_data$visit_age)


met_old = met_data %>% filter(visit_age > avg_age)

#use a 3 factor

pca_met <- principal(met_old[,3:8], rotate = "varimax", nfactors = 3)
pca_met$loadings

rm(met_old)


#loadings = scale(pca_met$loadings[1:6,]) 
#loadings
loadings = pca_met$loadings[1:6,]

#check the loadings
apply(loadings,2,sd)

#remove the attributes
attr(loadings, "scaled:scale") <- NULL
loadings

#calculate scores using the scaled loadings
scores = scale(met_data[,3:8]) %*% loadings 
head(scores)


#standardise the scores
scores[,1:3] = scale(scores[,1:3], center = FALSE, scale = apply(scores[,1:3], 2, sd, na.rm = TRUE))

#check that they have an SD of 0
apply(scores,2,sd)

#6) Square: square each individual value

square_scores = apply(scores, c(1,2), square)

head(square_scores)

#7) Sum: take the sum of the squared  values for each individual row

row_sum = rowSums(square_scores)

#8) Square root: take the square root of the sum of the squared values

MetSSS_old = sqrt(row_sum)

met_data = cbind(met_data, MetSSS_old)

rm(MetSSS_old)


# ------------- old ---------------- #


met_young = met_data %>% filter(visit_age < avg_age)

#use a 3 factor

pca_met <- principal(met_young[,3:8], rotate = "varimax", nfactors = 3)
pca_met$loadings

rm(met_young)


#loadings = scale(pca_met$loadings[1:6,]) 
#loadings
loadings = pca_met$loadings[1:6,]

#check the loadings
apply(loadings,2,sd)

#remove the attributes
attr(loadings, "scaled:scale") <- NULL
loadings

#calculate scores using the scaled loadings
scores = scale(met_data[,3:8]) %*% loadings 
head(scores)


#standardise the scores
scores[,1:3] = scale(scores[,1:3], center = FALSE, scale = apply(scores[,1:3], 2, sd, na.rm = TRUE))

#check that they have an SD of 0
apply(scores,2,sd)

#6) Square: square each individual value

square_scores = apply(scores, c(1,2), square)

head(square_scores)

#7) Sum: take the sum of the squared  values for each individual row

row_sum = rowSums(square_scores)

#8) Square root: take the square root of the sum of the squared values

MetSSS_young = sqrt(row_sum)

met_data = cbind(met_data, MetSSS_young)

rm(MetSSS_young)

# -------- test the correlation ------- #

# we want to look at the top row, the only one that's below their cutoff is old

cor = cor(met_data[10:14])[,1] %>% print()

summary(met_data$MetSSS)

# sex differences 
t.test(MetSSS ~ sex, data = met_data)

ggplot(data = met_data, aes(y = MetSSS, x = visit_age, color = as.factor(sex), group = sex)) + 
  geom_point() + 
  geom_smooth() + 
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) + 
  theme_classic()


# merge it back into the main dataset ##

#act_food = merge(act_food, met_data) --> Seems to lead to 0 observations in act_food

act_food_mss <- merge(act_food, met_data, by = "record_id")

write.csv(act_food_mss, "act_food_mss.csv", row.names = FALSE)


#---------------- GAUCCI ANALYSIS------------------- 

# check what you would like to exclude based on #

#Simple reaction time
# choice reaction time 
#Immediate/delayed recognition 
#Stroop colour-word
#Spatial Working memory 
#Contextual Memory 
# Cognitive domain calculations 

# 2)Dietary Assessment 
#Med Diet
# Dash Diet
# Mind diet. 

#3)Cardiovasular Health Measures, 

#4)Blood bIomarkers - lipids, glucose. 

#5)Anthropomorphic - Waist and Hip 

#6) Metabolic Syndrome Severity Score. 


colnames(gaucci)

gaucci_cols = c("visit_age",
                "bmi.x",
                "sex",
                "educationtotal",
                "waist_average",
                "hips_average",
                "LongTermMem",
                "ShortTermMem",
                "ExecFunc",
                "ProcSpeed",
                "plant.diet",
                "meat.diet",
                "western.diet",
                "sbp_average.x",
                "dbp_average.x",
                "chol_mmoll",
                "HDL_mmoll",
                "TRIGLY",
               "blood_glucose",
               "MetSSS"
                )




gaucci <- gaucci[,gaucci_cols] %>% 
#  as_tibble() %>%  
  na.omit()

p.corr = function(x,y){}


summary(gaucci)

cor = cor(gaucci) 



no_cor <- abs(cor) < .025 #select a pretty low cut-off for correlations. 
cor[no_cor] <- 0
cor


cog_per = c("LongTermMem",
            "ShortTermMem",
            "ExecFunc",
            "ProcSpeed")

attr = c("visit_age",
         "bmi.x",
         "sex",
         "educationtotal",
         "waist_average",
         "hips_average", 
         "MetSSS")

diets = c("plant.diet",
          "meat.diet",
          "western.diet")

blood_measures = c("sbp_average.x",
                   "dbp_average.x",
                   "chol_mmoll",
                   "HDL_mmoll",
                   "TRIGLY",
                   "blood_glucose")


ct.p = gaucci %>% 
  as.matrix() %>% 
  corr.test(method = "pearson") %>% 
  .$p %>% 
  round(4) %>% 
  print()

ct.r = gaucci %>% 
  as.matrix() %>% 
  corr.test(method = "pearson") %>% 
  .$r %>% 
  round(4) %>% 
  print()

#remove all non significant correlations

ct.p[ct.p >.05] <- NA
ct.r[ct.p >.05] <- NA

#filter to look at attributes with cognitive performance
ct.p[attr,cog_per] 
ct.r[attr,cog_per] 

cor.test(gaucci$visit_age, gaucci$ExecFunc,alternative = "two.sided")
cor.test(gaucci$ProcSpeed, gaucci$MetSSS, alternative = "two.sided")
t.test(gaucci$ExecFunc ~ gaucci$sex) # do a t-test to check it sex vs exec function

#can extract a specific correlation or p value 

ct.p["western.diet","ProcSpeed"]
ct.r["western.diet","ProcSpeed"]


#diets with cognitive performance
ct.p[diets,cog_per]
ct.r[diets,cog_per]

cor.test(gaucci$ProcSpeed, gaucci$meat.diet, alternative = "two.sided")

#blood measures with cognitive performance.
ct.p[blood_measures,cog_per]
ct.r[blood_measures,cog_per]

cor.test(gaucci$ProcSpeed, gaucci$TRIGLY, alternative = "two.sided")

#attributes with diets

ct.p[attr,diets]
ct.r[attr,diets]

#blood measures with diets 
ct.p[blood_measures,diets]
ct.r[blood_measures,diets]




## ------- RUNNING THE REGRESSION ---------------


# "Following the correlation analysis, separate hierarchical regressions were conducted, controlling for age, gender,
# education, and energy intake#

#lets see if the proc-speed and meat diet
ct.p[,"ProcSpeed"]

#first fit a model using processing speed, age, sex and education 
PS1 = lm(data = gaucci, ProcSpeed ~ visit_age + sex +educationtotal)
summary(PS1)

#only age was significant so lose the extra variables
PS2 = lm(data = gaucci, ProcSpeed ~ visit_age )
summary(PS2)
#test to see whether removing variables has significantly decreased the fit of the model
anova(PS1,PS2)
# no significant difference

#add in TRIGLYcerides 

PS3 = update(PS2, .~. + TRIGLY)
summary(PS3)

#test to see if the model is improved
anova(PS2,PS3)

#now add in meat_diet

PS4 = update(PS3, .~. + meat.diet)
summary(PS4)
anova(PS3,PS4)
#we can see that the inclusion of meat_diet improved the basic model just containing age. 
#we can say that the meat diet has a negative effect on processing speed, outside of the general
#decline in processing speed associated with age. 

#test to see if there is an interaction 

PS5 = update(PS4, .~. + meat.diet:visit_age)
anova(PS5,PS4)
#no improvement, so there is no interaction term

#diagnostics tests about the model

models = c("PS4, PS3")

for(m in models){
  paste0(m,"_diagn") = as_tibble(m$fitted.values)
  
}

PS4_diag = as_tibble(PS4$fitted.values) %>% 
  mutate(residuals =  resid(PS4),
          sta.resid = rstandard(PS4),
          stu.resid = rstudent(PS4),
          cooks = cooks.distance(PS4),
    #      df.beta = dfbeta(PS4),
          df.fit = dffits(PS4),
          leverage = hatvalues(PS4),
          cov.ratio = covratio(PS4)) %>% 
  round(3) %>% 
  as_tibble(rownames = "case")

head(PS4_diag)


#look at which are outliers using the standardised residuals
PS4_diag$outlier <- abs(PS4_diag$sta.resid) > 2

#count them up
table(PS4_diag$outlier)

#10/302 = 3%

PS4_diag[PS4_diag$outlier,] %>% arrange(residuals)

# check how many are outside 2.5 sd (1 case, this isn't cause for concern)
# check the cook's distance is all less than .5 (high influence) or 1 (extreme influence)
# df.fit should not be > 2*sqrt(k/n) = .20 (these are all high, cause for concern but not critical)
# average leverage = (predictors + 1)/n = 4/300 = .0133 
# we want to look for cases 2 or 3 times as large as this (.027 or .04)
# covariance ration 
# CVR > 1 + [3(k + 1)/n]  = 1.04, CVR < 1 – [3(k + 1)/n]  = 0.96
# so we look for covariance values > 1.04 and  < .96

#df beta looks at how much a coefficient would change if we removed that one case
PS4_dfbeta = dfbeta(PS4)
PS4_dfbeta[abs(PS4_dfbeta[1:4]) > 0.15,] %>% round(3) 

# Test for the assumption of independence of errors using Durbin Watson 
# We want a D-W statistic close to 2 and a non-significant p value

durbinWatsonTest(PS4)

#Assumption of no-collinearity

vif(PS4)
mean(vif(PS4))

# tollerance = 1/VIF
# tollerace needs to be above .2

1/vif(PS4)
mean(1/vif(PS4))

#------------------ residual diagnostics ------------------#

hist(PS4_diag$stu.resid, breaks = seq(-3.5,3.5, by = .5))   # check for normality 
plot(PS4_diag$value, PS4_diag$sta.resid) # check for independence 

qqnorm(PS4_diag$stu.resid)
qqline(PS4_diag$stu.resid, color = "red", lwd = 2)

qqPlot(PS4_diag$stu.resid)







EF1 = lm(data = gaucci, ExecFunc ~ visit_age+sex+educationtotal)
summary(EF1)
#maybe we take out education but for now we'll leave it in. 
EF2 = lm(data = gaucci, ExecFunc ~ visit_age+sex+ meat.diet + plant.diet + western.diet)
summary(EF2)

anova(EF1,EF2)

LM1 = lm(data = gaucci, LongTermMem ~ visit_age + sex + educationtotal)
summary(LM1)
LM2 = update(LM1, .~. + meat.diet + plant.diet + western.diet)
summary(LM2)






#------------- CODE GRAVEYARD -------------- 

# Sex-specific clinical thresholds
# female score is listed first. 

#Waist circumference, cm ≥80 ≥94
#TRIGLYcerides, mmol/L ≥1.7
#HDL cholesterol, mmol/L ≤1.3 ≤1.0
#Systolic BP, mm Hg ≥130
#Diastolic BP, mm Hg ≥85
#Blood glucose, mmol/L ≥5.6

head(met_data)


#SBP_sd = 13.28
#DBP_sd = 3.14
#TRIGLY_sd = 0.30
#HDL_sd = 0.17
#Waist_sd = 8.60
#BG_sd = 0.13

#met_data = met_data %>% 
#  mutate(waist_average = waist_average/Waist_sd,
#         TRIGLY = TRIGLY/TRIGLY_sd,
#         sbp_average_x = sbp_average_x/SBP_sd,
#        dbp_average_x = dbp_average_x/DBP_sd, 
blood_glucose = blood_glucose/BG_sd,
HDL_mmoll = HDL_mmoll/HDL_sd)

summary(met_data)

apply(met_data, 2, sd)

# SKIP THIS PART 
# It's mad complicated and doesn't work

#5) Standardize: post-multiply by a diagonal matrix comprised of the inverse of the standard deviations of the principal components shown
#in Table A.3

#take the sd of the loadings

#sd = matrix(apply(loadings, 2, sd))
sd

scores = scale(met_data[,3:8]) %*% pca_met$loadings 

cor(scores, pca_met$scores )

#and the inverse of it

inv_sd = matrix(pinv(sd), ncol = 3)

inv_sd

#times the scores by the inverse of the sd 

new_scores = scores
new_scores[,1:3] = scores[,1:3] * inv_sd[,1:3]

head(scores)
head(new_scores)

apply(new_scores,2,sd)

#IGNORE THAT LAST NONSENSE 