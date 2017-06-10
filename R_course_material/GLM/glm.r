
#the linear model (lm) is a specialized case of the generalized linear model (glm)
#The dependent variable of a GLM can take the form of any exponential distribution.  COmmonly GLM is used for logistic and poisson regression

#Logistic regression
#In logistic regression the dependent variable is binary.
#The logistic model uses log of odds (Logit) as a default link function, relating the response variable to the predictors by inverse logit 
#logit(Y)=Beta1+Beta2*X
#where logit(y)=y/(1-y)

#Here we will use the data in mtcars.  The following command loads this data set into memory
data(mtcars)

#the vs data describes whether this car has a V-style or straight (inline) cylinder design
#have a look at the vs data to see that it is binomial
mtcars$vs

#by attaching the data we can omit the name of the dataframe hereafter
attach(mtcars)
vs

vsmodel1=glm(vs~wt+disp,family="binomial"(link=logit))
#since logit is the default link for a binomial error distribution it is not necessary to include it
#prove this by querying the link value of the binomial family
binomial()$link
vsmodel1=glm(vs~wt+disp,family="binomial")

#examine the model using summary. 
summary(vsmodel1)
#Look at the coefficients and notice that displacement is negatively associated with use of V-style engines and wt is postively associated 
#the coefficients tell us how an increment by one unit will affect the log-odds (logit)
#also returned are the standard errors of the coefficients and the p value against the null hypothesis that each predictor has no effect
#the null deviance tells us the deviance of the data based on the inclusion of neither of the variables (i.e. based only on the intercept)

#we can predict the log-odds of a v-style engine through a linear combination of predictors.  The first coefficent is the intercept
#the next two are coefficients of wt and disp
newdata <- data.frame(wt=2.1, disp=180)

#the following two commands are equivalent and return the expected log-odds using these predictors
a=coef(vsmodel1)[1]+coef(vsmodel1)[2]*newdata$wt+coef(vsmodel1)[3]*newdata$disp
b=predict(vsmodel1,newdata)
#confirm that 
a==b

#we can translate this to the probability of a v-style engine using the inverse of the link function
#inverse logit is exp(x)/(1+exp(x))
exp(a)/(1+exp(a))

#we can simplify this using the response argument, which automatically converts log-odds to probability
predict(vsmodel1,newdata,type="response")



#POISSON MODEL
docfile=paste(getwd(),"/doctor_data.csv",sep="")

#Macs may need to use method "curl"
download.file("https://vincentarelbundock.github.io/Rdatasets/csv/Ecdat/DoctorAUS.csv", 
              destfile = docfile, method = "libcurl") 

#medical data from Australia 1977
#sex sex
#age age (0 male, 1 female)
#income annual income in tens of thousands of dollars
#insurance (medlevy, levyplus: private, freepoor: low income govt, freerepa old age govt)
#illness number of illness in past 2 weeks
#actdays number of days of reduced activity in past 2 weeks due to illness or injury
#hscore general health score using Goldberg's method (from 0 to 12)
#chcond chronic condition (np : no problem, la : limiting activity, nla : not limiting activity)
#doctorco number of consultations with a doctor or specialist in the past 2 weeks
#nondocco number of consultations with non-doctor health professionals
#hospadmi number of admissions to a hospital #in the past 12 months (up to 5+)
#hospdays number of nights in a hospital during most recent admission
#medecine total number of prescribed and nonprescribed medications used in past 2 days
#prescrib total number of prescribed medications used in past 2 days
#nonpresc total number of nonprescribed medications used in past 2 days

#let's look at the number of medications used by these people in australia
docdata=read.csv(docfile)
attach(docdata)
age=age*100

#first examine the response variable using a histogram
hist(medecine)

#a histogram shows that the response variable is highly skewed to the left
#as it clearly does not follow a normal distribution, the lm() is a poor choice.  Let's model
#this count data instead using a poisson model (used for count or rate data)

model1 <- glm(medecine~age+sex,family="poisson")
summary(model1)
#females take 48% more medication than males when correcting for age
#every year users take 1.9% more medications when correcting for sex

#however, maybe men get sicker as they age?  lets include a multiplicative term
model2 <- glm(medecine~age+sex+age:sex,family="poisson")
summary(model2)
#yes multiplicative term is significant and slope is adjusted down for older females 
#note that with the multiplicative term the number of medications taken by females are now twice as much

#does age:sex improve prediction of the data?
anova(model1, model2, test="Chisq")
#yes, a significant chi squared test indicates that extra parameter has some descriptive ability 
#but compare AIC score in summary to determine if Pmodel2 is more parsimoneous
AIC(model1)
AIC(model2)
#model1 AIC is 15439, model2 AIC is 15379 so we will keep the age:sex term

#Check for overdispersion
#a condition of the poisson is that conditional mean = conditional variance.  
#If variance > mean then the model is overdispersed, if variance < mean then underdispersed
#Check if this is true
mean(medecine)
var(medecine) 
#this isn't a test but its likely these data are overdispersed since variance>>mean

#perform a test for overdisperson using dispersiontest from AER package 
#if p value is very small we should be using quasi poisson or neg. binomial model which include 
#additional parameter for dispersion

#install.packages("AER") #UNCOMMENT ON FIRST RUN
library(AER)

dispersiontest(model2)
#as dispersion > 1 the data are indeed overdispersed and so violate the poisson assumptions.

#lets use a quasi-poisson model
model3 <- glm(medecine~age+sex+age:sex,family="quasipoisson")

summary(model2)$coefficients
summary(model3)$coefficients
#notice that the coefficiences are identical from the poisson and quasipoisson but the standard error
#is greater from the quasipoission, which is a more robust method and better represents the error

#use this model to predict the number of medications taken by a female age 50
predict(model3,data.frame(age=50,sex=1))


#negative binomial is suitable for overdispersion (but not underdispersion)
library(MASS)
negbin <- glm.nb(medecine~age+sex+age:sex,link=log)

#compare the AIC with the previous model2 (15379)
AIC(negbin) #14810

#the quasipoisson and negative binomial produce the similar coefficients but st error is different
summary(model3)$coefficients
summary(negbin)$coefficients

#Note: for zero-inflated data options include trucated poisson, truncated negative binomial, quasi-poisson,
#zero inflated poisson and zero inflated neg binomial (see hurdle() of the AER package)

#Mixed effects
#mixed effects models describe the relationship between a response variable and independent variables with 
#coefficients that can vary with respect to one or more grouping variables.  They contain fixed and random effects
#Fixed effects are the conventional linear regression variables.  Random effects are associated with
#individual experimental units.  
#in this example we will regress oat yield against fertilizer concentration, accounting for crop strain as a
#random effect.  It is a 'random' effect because we have sampled only 3 strains out of a large universe of
#potential strains.  

library(nlme) #contains the lme() function
attach(Oats) #comes with R
head(Oats)

oatsmodel1=lm(yield~nitro) #usual linear model
oatsmodel2=lme(yield~nitro,random=~1|Variety) #accounts for variety as random effect (and nitro as fixed)

coef(oatsmodel1) #provides slope and (global) intercept assuming all samples are independent
coef(oatsmodel2) #provides random intercept for each strain

