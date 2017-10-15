## Save the data file to a folder on your computer
## Now we will write the actual path to the exact file we want
## remember to change the slashes

setwd("C:\\Users\\tim\\Documents\\R")
cheese <- read.csv("C:/Users/tim/Documents/R/Cheese.csv", header=T)

head(cheese)
##this data shows a subjetive rating of the tatse of cheese (taste) and the corresponding levels of 
##Acetic acid (acetic), Hydrogen Sulfide (H2S), Lactic acid (Lactic), the location 
##of where it came from (farm) and the observer id (observer)

plot(cheese) ## there is possible multicollinearity between H2S and Lactic Acid

## from last week
## t.test
t.test1 <- t.test(cheese$taste ~ cheese$farm)
t.test1
## ANOVA
amodel1 <- aov(taste~farm+observer+farm:observer, data=cheese)
summary(amodel1)

amodel2 <- aov(taste~farm*observer, data=cheese)
summary(amodel2)

amodel3 <- lm(taste~farm*observer, data=cheese)
summary(amodel3) ## can compare among catergorical groups

cheese$obs.b.ref <- relevel(cheese$observer, ref="B") ## change the reference group

amodel4 <- lm(taste~farm*obs.b.ref, data=cheese) ## rerun the model, notice you have to change the observer variable
summary(amodel4)

##ANCOVA with one categorical variable and one continuous variable
amodel5 <- lm(taste~Acetic*observer, data=cheese)
summary(amodel5)

amodel6 <- lm(taste~Acetic+observer, data=cheese)
summary(amodel6)

## compare the two models using F-test
anova(amodel5,amodel6)

## Ok on to model selection of multiple linear regression
## Show three ways to do it
## 1. by hand
## 2. using the update function which allows you to take away or add variables
## 3. by machine (uses the step function)

mlmodel1 <- lm(taste~Acetic*H2S*Lactic, data=cheese)
summary(mlmodel1)

mlmodel2 <- lm(taste~Acetic+H2S+Lactic+Acetic:H2S+Acetic:Lactic+H2S:Lactic, data=cheese)
summary(mlmodel2)

mlmodel3 <- lm(taste~Acetic+H2S+Lactic+Acetic:H2S+Acetic:Lactic, data=cheese)
summary(mlmodel3)

mlmodel4 <- lm(taste~Acetic+H2S+Lactic+Acetic:H2S, data=cheese)
summary(mlmodel4)

mlmodel5 <- lm(taste~Acetic+H2S+Lactic, data=cheese)
summary(mlmodel5)

mlmodel6 <- lm(taste~H2S+Lactic, data=cheese)
summary(mlmodel6)

mlmodel7 <- lm(taste~H2S, data=cheese)
summary(mlmodel7)

mlmodel8 <- lm(taste~Lactic, data=cheese)
summary(mlmodel8)

## obtain the AIC (Akaike's information criterion) 
## Use the AICcmodavg package to get the AICc value
AIC(mlmodel1)##235.1119
AIC(mlmodel2)##233.2895
AIC(mlmodel3)##232.1373
AIC(mlmodel4)##230.4944
AIC(mlmodel5)##229.7775
AIC(mlmodel6)##227.7838!!!!
AIC(mlmodel7)##232.0245
AIC(mlmodel8)##236.8724

## use anova to compare the models
anova(mlmodel5, mlmodel6)
anova(mlmodel6, mlmodel7)

anova(mlmodel3,mlmodel4,mlmodel5,mlmodel6,mlmodel7)

#### now that we have our model we can look at some diagnostics
plot(mlmodel6)## shows some diagnostics to evaluate whether you met the assumptions of the test

## With each of these models you have created a object that will creates a variety of output
?lm
AIC(45) ## the AIC function from above can only work on a class of "lm" value
coefficients(mlmodel6)
confint(mlmodel6)
residuals(mlmodel6)
fitted.values(mlmodel6)
plot(mlmodel6$residuals~mlmodel6$fitted.values)

### use the update function 
mlmodel1
mlmodel.update.2=update(mlmodel1, .~. - Acetic:H2S:Lactic)
summary(mlmodel.update.2)

mlmodel.update.3=update(mlmodel2, .~. - H2S:Lactic)
summary(mlmodel.update.3)

mlmodel.update.4=update(mlmodel3, .~. - Acetic:Lactic)
summary(mlmodel.update.4)

mlmodel.update.5=update(mlmodel4, .~. - Acetic:H2S)
summary(mlmodel.update.5)

AIC(mlmodel.update.2)##233.2895
AIC(mlmodel.update.3)##232.1373
AIC(mlmodel.update.4)##230.4944
AIC(mlmodel.update.5)##229.7775


### now to use the step function

null=lm(taste~1, data=cheese)
summary(null)
full=lm(taste~Acetic*H2S*Lactic, data=cheese)
summary(full)

step(null, scope=list(lower=null, upper=full), direction="forward") ##forwards

step(full, data=Housing, direction="backward") ## backwards

step(null, scope = list(upper=full), data=Housing, direction="both") ##both

### now lets practice
##Species on British Islands
brit <- read.csv("C:/Users/TJ/Documents/R/Data/r class/brit.species.csv", header=T)
head(brit)
plot(brit)

brit.2 <- brit[-c(6),]
plot(brit.2)

##now you try