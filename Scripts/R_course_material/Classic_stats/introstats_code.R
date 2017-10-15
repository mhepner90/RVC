### R Coding Clinic - Intro to t-tools and extensions ####
### Jan 30, 2017; Chris Stallings stallings@usf.edu ###

## t-tools are used to determine whether the means of different groups differ
## note that this is not a stats course; we assume you understand the premise and assumptions of these tests
## CMS offers several stats courses for those that need formal training (Data Analysis, Biometry, Multivariate Stats)

## set your working directory and save code and space files accordingly

## read in the csv file
dodata=read.csv("do_data.csv")
names(dodata)

## metadata (note that these are fake data used for this lecture)
### estuary (estuary ID, A-D); site (site ID within each estuary, 1-10);
### do (dissolved oxygen); nutrients (scaled nutrient load); 
### eutrophic (1=no, 2=yes, based on scaled nutrient load w/cutoff at 80);
### hpd (human population density within 5 km^2 of DO measurement)

## lets suppose we are interest in whether DO differs between eutrophic and non-eutrophic (or mesotrohic) estuaries

## we need to tell R that the "eutrophic" column contains categorical info, not numeric data
## do this with "as.factor"
dodata$eutrophic <- as.factor(dodata$eutrophic)

## always start by plotting your data
## since we are interested in the DO data, let's first focus on those

par(mfrow=c(1,3)) ## sets stage for combining plots; par=graphical parameters; mfrow=ncol, nrows
plot(dodata$do) ## scatterplot in order data appear in dataframe
boxplot(dodata$do) ## box and whisker plot (25th, 50th, 75th + whiskers 1.5x1QR)
hist(dodata$do, main="") ## frequency plot/histogram
## we don't see any spurious values, outliers, etc with the DO data

## we can also summarize our data
summary(dodata$do)

## let's look at the DO data between the eutrophic and mesotrophic estuaries
par(mfrow=c(1,1)) ## note this is the default, but we established a 1-row, 3-column format above, so have to reset here
boxplot(dodata$do~dodata$eutrophic)
## notice the single high value above the upper whisker in the mesotrophic group (group 1)

## parametric tests require the data to meet assumptions of having a normal distribution

## one of the simplest ways is a QQ plot (plots our ranked data against that of a norm dist)
qqnorm(dodata$do)
qqline(dodata$do, lty=2)  ## note that lty is the line type (2=medium dashes; try different values b/t 1-6)
## we don't see any major non-linearities (departures from normality) such as S- or banana- shapes

## can test whether data are NOT normally distributed using Shapiro-Wilk Normaility Test
shapiro.test(dodata$do)
## but be aware that this tests the null that data come from normal distribution; thus,
## we have failed to reject the null (but think about what that means regarding the burden or proof)
## there is also a package dedicated to tests for normality; library(nortest)

## skewness is a measure of the extent to which a distribution has long tails on either side
## a normal distribution is symmetrical and a skewness equal to zero
## to check the skewness of your data, first write this function (you can copy and paste this for future use)
skew<-function(x){
  m3<-sum((x-mean(x))^3)/length(x)
  s3<-sqrt(var(x))^3
  m3/s3}
## note that the length expression automatically determines the sample size of your data

## then run the skew function on your data
skew(dodata$do)

## now we need to check whether the calculated skew value is significantly different than zero (i.e., not normal)
## we can do this using a t-test; to calulate the t-value, divide the skew by its standard error:
skew(dodata$do)/sqrt(6/length(dodata$do))
## note that a rough measure of the standard error of the skewness is "sqrt{6/n}" where n is the sample size

## now check the probability of this t-value by chance alone when the skew is zero using the following:
1-pt(-0.02247742, 39) ##pt=cummulative density function of t-dist; first number=calculated t-value; second=df

## we do not find evidence of skewness different from zero, so transformations are not needed on the DO data
## note that possitively skewed data are quite common and can be normalized with both sqrt and log trans
## sqrt trans: "skew(sqrt(values))/sqrt(6/length(values))" -- then check the t-value as in above
## log trans: "skew(log(values))/sqrt(6/length(values))" -- then check the t-value as in above

## parametric tests must also meet assumption of homoscedasticity (equal variance)
## can test this using a few different approaches
var.test(do~eutrophic, data=dodata) ##Fisher test; can only test data with 2 levels, sensitive to outliers
fligner.test(do~eutrophic, data=dodata) ##Barlett test; can test data with >2 levels, sensitive to outliers
bartlett.test(do~eutrophic, data=dodata) ##Fligner-Killeen test; ; can test data with >2 levels, not sensitive to outliers
## we have failed to reject null of homoscedasticity

## we have sufficient evidence to proceed with parametric test(s) without the need for transforming the data
t.test(do~eutrophic, data=dodata)
## we reject the null that dissolved oxygen does not differ between eutrophic/non-eutrophic estuaries

## we could have also run a t-test on data without a factor, if they were in separate columns
## using "t.test(y1,y2)" where y1 and y2 are numeric

## or a paired t-test (assuming a paired design) using "t.test(y1,y2,paired=TRUE)" where y1 & y2 are numeric
## note the default is paired=FALSE

## the non-parametric extension analog of the 2-group t-test is the Mann-Whitney U Test (aka Wilcoxon Rank Sum Test)
wilcox.test(do~eutrophic, data=dodata)
## note the warning (not an error)

## t-tests and Mann-Whitney are good for testing whether two groups differ
## to test whether more than two groups differ, we use ANOVA (if data meet parametric assumptions)

## tell R that estuary is a factor (the A-D structure could be intepretted as ordinal data)
dodata$estuary <- as.factor(dodata$estuary)

## plot the data
boxplot(do~estuary, data=dodata)

## test for homoscedasticity
fligner.test(do~estuary, data=dodata)
bartlett.test(do~estuary, data=dodata)

## run an ANOVA using the aov function
## we can also do this with the lm function (next class)
do1.aov<-aov(do~estuary, data=dodata)

## note that this runs and stores the model we named "do1.aov" in the environment pane
## we use the summary function to call the results
summary(do1.aov)
## we reject the null that DO does not differ among the estuaries
## however, we do not know which one(s) differed based on this omnibus test

## Tukey's HSD allows us to compare among each pair of estuaries
TukeyHSD(do1.aov)

## as with t-tests, there are non-parametric analogs to the ANOVA
## such as the Kruskal Wallis test; "kruskal.test(y~A)" where y1 is numeric and A is a factor

## note that we are conducting simple one-way ANOVAs here
## next week we will build more complicated models with multiple main predictors and interactions
## then procede with model selection

## let's look at two sets of continuous data (rather than continuous response and categorical predictors as above)
## plot the data
plot(dodata$do~dodata$nutrients)

## determine their correlation and whether it is statistically significant
cor.test(dodata$do,dodata$nutrients, method="spearman")

## other simple/classic stats not covered:
## chi-squared test (do observed data differ from expected?); "chisq.test" and "fisher.test"
## correlation matrices (correlation between all vectors in a dataframe); "cor(dataframe)"
## kolmogorov-smirnov test (are two data distributions different?); "ks.test(A,B)"

## Two activities for you to conduct on your own

## Activity 1
## Create Anscombe plots to illustarte the importance of plotting your data prior to analysis
# Code from http://blog.ouseful.info/2011/08/30/the-visual-difference-%e2%80%93-r-and-anscombe%e2%80%99s-quartet/
library(ggplot2)
mydata=with(anscombe,data.frame(xVal=c(x1,x2,x3,x4), yVal=c(y1,y2,y3,y4),group=gl(4,nrow(anscombe))))
aggregate(.~group,data=mydata,mean)
aggregate(.~group,data=mydata,sd)

ggplot(mydata,aes(x=xVal, y=yVal)) + geom_point() + facet_wrap(~group)

## Activity 2
## Determine whether the estimated biomass of Jaguar Sharks observed by Captain Steve Zissou differed among four islands
## data: jagshark_data.xlsx
## metadata: island=four islands visited by Zissou (A-D); jagshark=estimated biomass in kg of observed Jaguar Sharks

