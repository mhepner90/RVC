#Analysis for Glaspie and Seitz (2017) Role of habitat and predators in maintaining functional 
#diversity of estuarine bivalves, Marine Ecology Progress Series 570:113-125

#R version 3.2.5 (2016-04-14)

setwd("....")

mya=read.csv("Bivalve_Survey.csv",header=T)

#Remove metadata

mya=mya[1:216,]

#Remove fall 2011 sampling period due to lack of predator data

mya=mya[mya$Year!="2011",]

#Sort by substrate type and abundance for rarefaction.

mya=mya[order(mya$Sb,mya$Abundance ), ]


#Read in columns
#Abundance converted to density, sample area 0.11 m^2

river=mya$River
site=mya$Site
season=mya$Season
year=as.factor(mya$Year)
sb=mya$Sb
volume=mya$Volume/1000 #T model in L rather than mL
mya$Rays=mya$Rays/100 #To correct for transect area
rays=mya$Rays
mya$Crabs=mya$Crabs/40 #To correct for transect area
crabs=mya$Crabs
mya$Fish=mya$Fish/40 #To correct for transect area
fish=mya$Fish
diversity=mya$Diversity
richness=mya$Richness
abundance=mya$Abundance/0.11 #To correct for sample area
hard=mya$Hard.shelled/0.11 #To correct for sample area
deposit=mya$Deposit.feeders/0.11 #To correct for sample area
DBSF=mya$DBSF/0.11 #To correct for sample area
TSSB=mya$TSSD/0.11 #To correct for sample area
func.div=mya$Func.div
func.rich=mya$Func.rich

#install.packages("iNEXT")

library(iNEXT)

# Rarefaction -------------------------------------------------------------

#Calculate incidence frequency for rarefaction

DT.mya=mya[sb=="DT",]
n.start=which(colnames(mya)=="T.No" )
n.end=which( colnames(mya)=="M.modiolus.No" )
sub.inc.freq=vector()
for(i in seq(n.start,n.end,1))
{
sub.inc.freq[i]=(sum(DT.mya[,i]>=1))
}
inc.freq.DT=na.omit(sub.inc.freq)

GR.mya=mya[sb=="GR",]

sub.inc.freq=vector()
for(i in seq(n.start,n.end,1))
{
  sub.inc.freq[i]=(sum(GR.mya[,i]>=1))
}
inc.freq.GR=na.omit(sub.inc.freq)

OS.mya=mya[sb=="OS",]

sub.inc.freq=vector()
for(i in seq(n.start,n.end,1))
{
  sub.inc.freq[i]=(sum(OS.mya[,i]>=1))
}
inc.freq.OS=na.omit(sub.inc.freq)

SG.mya=mya[sb=="SG",]

sub.inc.freq=vector()
for(i in seq(n.start,n.end,1))
{
  sub.inc.freq[i]=(sum(SG.mya[,i]>=1))
}
inc.freq.SG=na.omit(sub.inc.freq)

SH.mya=mya[sb=="SH",]

sub.inc.freq=vector()
for(i in seq(n.start,n.end,1))
{
  sub.inc.freq[i]=(sum(SH.mya[,i]>=1))
}
inc.freq.SH=na.omit(sub.inc.freq)

# Estimating rarefied richness

obj=iNEXT(c(nrow(DT.mya),inc.freq.DT), q=0, datatype="incidence_freq",size=DT.mya$Abundance)

qD=obj$iNextEst[,4]
t=obj$iNextEst[,1]

rare.rich.DT=vector()

for(i in 1:length(DT.mya$Abundance))
{
  rare.rich.DT[i]=qD[t==DT.mya$Abundance[i]]
  
}

obj=iNEXT(c(nrow(GR.mya),inc.freq.GR), q=0, datatype="incidence_freq",size=GR.mya$Abundance)

qD=obj$iNextEst[,4]
t=obj$iNextEst[,1]

rare.rich.GR=vector()

for(i in 1:length(GR.mya$Abundance))
{
  rare.rich.GR[i]=qD[t==GR.mya$Abundance[i]]
  
}

obj=iNEXT(c(nrow(OS.mya),inc.freq.OS), q=0, datatype="incidence_freq",size=OS.mya$Abundance)

qD=obj$iNextEst[,4]
t=obj$iNextEst[,1]

rare.rich.OS=vector()

for(i in 1:length(OS.mya$Abundance))
{
  rare.rich.OS[i]=qD[t==OS.mya$Abundance[i]]
  
}

obj=iNEXT(c(nrow(SG.mya),inc.freq.SG), q=0, datatype="incidence_freq",size=SG.mya$Abundance)

qD=obj$iNextEst[,4]
t=obj$iNextEst[,1]

rare.rich.SG=vector()

for(i in 1:length(SG.mya$Abundance))
{
  rare.rich.SG[i]=qD[t==SG.mya$Abundance[i]]
  
}

obj=iNEXT(c(nrow(SH.mya),inc.freq.SH), q=0, datatype="incidence_freq",size=SH.mya$Abundance)

qD=obj$iNextEst[,4]
t=obj$iNextEst[,1]

rare.rich.SH=vector()

for(i in 1:length(SH.mya$Abundance))
{
  rare.rich.SH[i]=qD[t==SH.mya$Abundance[i]]
  
}

rare.rich=c(rare.rich.DT,rare.rich.GR,rare.rich.OS,rare.rich.SG,rare.rich.SH)

###############################################################################
#Looking at correlation between variables

library(lattice)


predictor.mat=cbind(diversity,richness,func.rich,func.div,year,season,river,sb,volume,fish,crabs,rays)
head(predictor.mat)

splom(predictor.mat)
###############################
#Diversity and functional diversity analysis

# Average overall sample Gini-Simpson diversity index and rarefied species 
#richness were 0.37 and 8.44, respectively.
mean(diversity)
mean(rare.rich)

#Code for substitution to change which category is included in the intercept
#sb=sub("SG","A",sb,ignore.case=F)

#Bootstrapping function for multiple comparisons
#Will produce slightly different results each time.

boots=function(response,lower="DT",higher="SG")
{
  diversity1=na.omit(response[sb==lower|sb==higher])
  sb1=na.omit(sb[sb==lower|sb==higher])
  calc.s=mean(diversity1[sb1==lower])-mean(diversity1[sb1==higher])
  
  l.Spring=length(diversity1[sb1==lower])
  l.total=length(diversity1)
  
  rand.s=NULL
  #Another way to make a blank data set to put bootstrapped means
  n=9999
  
  for(i in 1:n)
  {
    samp=sample(diversity1, size = l.total, replace = T)
    rand.s[i]=mean(samp[1:l.Spring])-mean(samp[(l.Spring+1):l.total])
  }
  
  sim.diffs=length(rand.s[rand.s<=calc.s])
  sim.diffs
  
  pvalue=(sim.diffs+1)/(n+1)
  pvalue
}

#Function to calculate Cohen's d for effect size
ES=function(a,b)
{
  (mean(a)-mean(b))/sqrt(((length(a)-1)*var(a))+((length(b)-1)*var(b))/(length(a)+length(b)-2))
}


lm11.5=glm(rare.rich~year+season+river+site:river+sb+volume+fish+crabs+rays,family="gaussian"(link="identity"))
par(mfrow=c(2,2))
plot(lm11.5)

summary(lm11.5)
confint(lm11.5)

# (Intercept)                               3.3153     1.6619   1.995 0.047790 *  
# year2013                                  0.4039     0.6549   0.617 0.538232    
# seasonSpring                              1.5215     0.8891   1.711 0.088999 .  
# seasonSummer                              0.5857     0.9396   0.623 0.533969    
# riverMobjack                             -2.2128     1.6474  -1.343 0.181135    
# riverYork                                -2.1372     1.7134  -1.247 0.214141    
# sbGR                                      2.9208     1.5115   1.932 0.055106 .  
# sbOS                                      4.2884     1.3063   3.283 0.001266 ** 
# sbSG                                      6.8827     1.4134   4.870  2.7e-06 ***
# sbSH                                      3.5702     1.1799   3.026 0.002899 ** 
# volume                                    0.7653     0.3749   2.042 0.042867 *  
# fish                                      0.9620     1.0466   0.919 0.359450    
# crabs                                    -1.8701     1.2645  -1.479 0.141169    
# rays                                     -3.9440     4.5431  -0.868 0.386652    
# riverYork:siteKing's Creek                1.6948     1.7935   0.945 0.346127    
# riverLynnhaven:siteLinkhorn Bay          -0.4355     1.5554  -0.280 0.779822    
# riverLynnhaven:siteLynnhaven Bay         -1.7247     1.5115  -1.141 0.255590    
# riverMobjack:siteMobjack 1                0.3617     1.4053   0.257 0.797199    
# riverMobjack:siteMobjack 2               -0.1156     1.4037  -0.082 0.934485    
# riverMobjack:siteMobjack 3                0.3132     1.3611   0.230 0.818310    
# riverLynnhaven:sitePleasure House Creek   0.5287     1.4424   0.367 0.714477    
# riverYork:sitePropotank                   6.1479     1.6573   3.710 0.000288 ***
# riverYork:sitePurtan                      5.4255     1.7119   3.169 0.001837 ** 
 

#McFadden's R squared
1-(lm11.5$deviance/lm11.5$null.deviance)

# 0.2444964

#In comparison with detrital mud (the habitat with lowest diversity), seagrass 
#supported 4.11 - 9.65 more species


# 2.5 %    97.5 %
# (Intercept)                               0.05797225 6.5726864
# year2013                                 -0.87956282 1.6874587
# seasonSpring                             -0.22108109 3.2641150
# seasonSummer                             -1.25594181 2.4273557
# riverMobjack                             -5.44160390 1.0159897
# riverYork                                -5.49546076 1.2210864
# sbGR                                     -0.04160961 5.8831970
# sbOS                                      1.72818983 6.8486429
# sbSG                                      4.11256907 9.6528681
# sbSH                                      1.25754611 5.8828347
# volume                                    0.03059502 1.5000262
# fish                                     -1.08937760 3.0132784
# crabs                                    -4.34849141 0.6083055
# rays                                    -12.84825843 4.9603386
# riverYork:siteKing's Creek               -1.82039766 5.2100045
# riverLynnhaven:siteLinkhorn Bay          -3.48398290 2.6128867
# riverLynnhaven:siteLynnhaven Bay         -4.68726955 1.2378268
# riverMobjack:siteMobjack 1               -2.39265095 3.1161447
# riverMobjack:siteMobjack 2               -2.86677761 2.6356280
# riverMobjack:siteMobjack 3               -2.35444184 2.9808149
# riverLynnhaven:sitePleasure House Creek  -2.29843420 3.3557657
# riverYork:sitePropotank                   2.89963295 9.3961767
# riverYork:sitePurtan                      2.07019390 8.7807339


# Rarefied species richness was significantly greater in seagrass than in 
# detritus (p = 0.04, d = 0.08; Fig. 3a)

boots(response=rare.rich,higher="SG",lower="DT") #0.0125
p.adjust(0.0125, method = "bonferroni", n = 3)  ###0.0375
ES(rare.rich[sb=="SG"],rare.rich[sb=="DT"])  ###0.07577

boots(response=rare.rich,higher="SG",lower="SH") #0.0412
p.adjust(0.0412, method = "bonferroni", n = 3)  ###0.1236
ES(rare.rich[sb=="SG"],rare.rich[sb=="SH"])  ###0.0496

boots(response=rare.rich,higher="OS",lower="DT") #0.0638
p.adjust(0.0638, method = "bonferroni", n = 3)  ###0.1914
ES(rare.rich[sb=="OS"],rare.rich[sb=="SH"])  ###0.0382


lm30=glm(func.rich~year+season+river+site:river+sb+volume+fish+crabs+rays,family="gaussian"(link="identity"))
par(mfrow=c(2,2))
plot(lm30)

summary(lm30)
confint(lm30)

# (Intercept)                              1.16422    0.43822   2.657 0.008705 ** 
# year2013                                 0.08315    0.17267   0.482 0.630804    
# seasonSpring                             0.39793    0.23443   1.697 0.091600 .  
# seasonSummer                             0.05973    0.24776   0.241 0.809800    
# riverMobjack                            -0.36738    0.43437  -0.846 0.398968    
# riverYork                               -0.53956    0.45179  -1.194 0.234174    
# sbGR                                     0.28470    0.39854   0.714 0.476061    
# sbOS                                     0.22759    0.34443   0.661 0.509731    
# sbSG                                     1.46930    0.37267   3.943 0.000121 ***
# sbSH                                     0.38631    0.31112   1.242 0.216209    
# volume                                   0.21257    0.09884   2.151 0.033036 *  
# fish                                     0.10542    0.27597   0.382 0.702981    
# crabs                                   -0.19943    0.33342  -0.598 0.550616    
# rays                                     1.45519    1.19791   1.215 0.226276    
# riverYork:siteKing's Creek              -0.07571    0.47290  -0.160 0.873019    
# riverLynnhaven:siteLinkhorn Bay          0.23803    0.41011   0.580 0.562478    
# riverLynnhaven:siteLynnhaven Bay         0.09576    0.39856   0.240 0.810444    
# riverMobjack:siteMobjack 1              -0.01630    0.37055  -0.044 0.964971    
# riverMobjack:siteMobjack 2              -0.28895    0.37012  -0.781 0.436166    
# riverMobjack:siteMobjack 3              -0.18660    0.35888  -0.520 0.603826    
# riverLynnhaven:sitePleasure House Creek  0.01404    0.38033   0.037 0.970592    
# riverYork:sitePropotank                  0.79033    0.43699   1.809 0.072431 .  
# riverYork:sitePurtan                     1.23789    0.45139   2.742 0.006808 ** 

1-(lm30$deviance/lm30$null.deviance)
#0.213533

# 2.5 %    97.5 %
# (Intercept)                              0.30533708 2.0231128
# year2013                                -0.25528285 0.4215799
# seasonSpring                            -0.06155095 0.8574127
# seasonSummer                            -0.42586704 0.5453312
# riverMobjack                            -1.21873790 0.4839765
# riverYork                               -1.42505874 0.3459356
# sbGR                                    -0.49641491 1.0658162
# sbOS                                    -0.44748359 0.9026586
# sbSG                                     0.73887754 2.1997231
# sbSH                                    -0.22347639 0.9961026
# volume                                   0.01884519 0.4062994
# fish                                    -0.43546810 0.6463051
# crabs                                   -0.85292412 0.4540658
# rays                                    -0.89266670 3.8030384
# riverYork:siteKing's Creek              -1.00258037 0.8511701
# riverLynnhaven:siteLinkhorn Bay         -0.56577312 1.0418269
# riverLynnhaven:siteLynnhaven Bay        -0.68539774 0.8769098
# riverMobjack:siteMobjack 1              -0.74256890 0.7099700
# riverMobjack:siteMobjack 2              -1.01437435 0.4364796
# riverMobjack:siteMobjack 3              -0.88999205 0.5167887
# riverLynnhaven:sitePleasure House Creek -0.73139584 0.7594827
# riverYork:sitePropotank                 -0.06615766 1.6468270
# riverYork:sitePurtan                     0.35317990 2.1225903


#On average, seagrass supported 0.74 - 2.20 more functional groups than detrital 
#mud, 0.41 - 1.96 more than coarse sand, 0.56 - 1.93 more than oyster shell, and 
#0.53 - 1.64 more than shell hash (95% CIs).

# sbDT                                    -2.19972312 -0.7388775
# sbGR                                    -1.96110133 -0.4080980
# sbOS                                    -1.92595139 -0.5574743
# sbSH                                    -1.63534634 -0.5306281

#Functional richness was greater in seagrass than in oyster shell (p = 0.04, d = 0.07), 
#shell hash (p = 0.007, d = 0.08) and detrital mud (p = 0.05, d = 0.08; Fig. 3b). 

boots(response=func.rich,lower="OS") #0.0119
p.adjust(0.0119, method = "bonferroni", n = 3)  ###0.0357
ES(func.rich[sb=="SG"],func.rich[sb=="OS"])  ###0.074

boots(response=func.rich,lower="SH") #0.0022
p.adjust(0.0022, method = "bonferroni", n = 3)  ###0.0066
ES(func.rich[sb=="SG"],func.rich[sb=="SH"])  ###0.078

boots(response=func.rich) #0.0152
p.adjust(0.0152, method = "bonferroni", n = 3)  ###0.0456
ES(func.rich[sb=="SG"],func.rich[sb=="DT"])  ###0.078


lm10=glm(diversity~year+season+river+site:river+sb+volume+fish+crabs+rays,family="gaussian"(link="identity"))
par(mfrow=c(2,2))
plot(lm10)

summary(lm10)
confint(lm10)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                              0.162820   0.107439   1.515  0.13167    
# year2013                                -0.078763   0.042335  -1.860  0.06469 .  
# seasonSpring                             0.075924   0.057477   1.321  0.18844    
# seasonSummer                            -0.015587   0.060744  -0.257  0.79783    
# riverMobjack                            -0.045750   0.106497  -0.430  0.66808    
# riverYork                               -0.117011   0.110768  -1.056  0.29243    
# sbGR                                     0.046051   0.097711   0.471  0.63808    
# sbOS                                     0.094198   0.084446   1.115  0.26634    
# sbSG                                     0.372115   0.091370   4.073 7.35e-05 ***
# sbSH                                     0.079328   0.076279   1.040  0.29995    
# volume                                   0.032685   0.024234   1.349  0.17936    
# fish                                     0.089721   0.067660   1.326  0.18675    
# crabs                                   -0.005353   0.081747  -0.065  0.94788    
# rays                                     0.767913   0.293696   2.615  0.00980 ** 
# riverYork:siteKing's Creek              -0.103954   0.115944  -0.897  0.37131    
# riverLynnhaven:siteLinkhorn Bay          0.037147   0.100548   0.369  0.71230    
# riverLynnhaven:siteLynnhaven Bay         0.104572   0.097716   1.070  0.28619    
# riverMobjack:siteMobjack 1              -0.006096   0.090850  -0.067  0.94659    
# riverMobjack:siteMobjack 2              -0.060715   0.090745  -0.669  0.50443    
# riverMobjack:siteMobjack 3              -0.076278   0.087988  -0.867  0.38731    
# riverLynnhaven:sitePleasure House Creek -0.064810   0.093248  -0.695  0.48807    
# riverYork:sitePropotank                  0.186084   0.107140   1.737  0.08438 .  
# riverYork:sitePurtan                     0.294372   0.110669   2.660  0.00863 ** 

#McFadden's R squared
1-(lm10$deviance/lm10$null.deviance)
#0.2820634

#In comparison with detrital mud (the habitat with lowest diversity), seagrass had 
#0.19 - 0.55 units higher Gini-Simpson diversity index

# 2.5 %      97.5 %
# (Intercept)                             -0.04775713 0.373397559
# year2013                                -0.16173756 0.004211882
# seasonSpring                            -0.03672871 0.188577672
# seasonSummer                            -0.13464309 0.103469891
# riverMobjack                            -0.25448110 0.162980941
# riverYork                               -0.33411245 0.100090079
# sbGR                                    -0.14545902 0.237560106
# sbOS                                    -0.07131170 0.259708631
# sbSG                                     0.19303397 0.551195955
# sbSH                                    -0.07017684 0.228832744
# volume                                  -0.01481188 0.080181973
# fish                                    -0.04289086 0.222332258
# crabs                                   -0.16557284 0.154867662
# rays                                     0.19227976 1.343546573
# riverYork:siteKing's Creek              -0.33119968 0.123292532
# riverLynnhaven:siteLinkhorn Bay         -0.15992466 0.234217762
# riverLynnhaven:siteLynnhaven Bay        -0.08694668 0.296091174
# riverMobjack:siteMobjack 1              -0.18415849 0.171966905
# riverMobjack:siteMobjack 2              -0.23857105 0.117141245
# riverMobjack:siteMobjack 3              -0.24873181 0.096174867
# riverLynnhaven:sitePleasure House Creek -0.24757224 0.117953057
# riverYork:sitePropotank                 -0.02390622 0.396073812
# riverYork:sitePurtan                     0.07746457 0.511278753

 
boots(response=diversity) #0.0024
p.adjust(0.0024, method = "bonferroni", n = 4)  ###0.0096
ES(diversity[sb=="SG"],diversity[sb=="DT"])  ###0.112

boots(response=diversity,lower="SH") #0.0002
p.adjust(0.0002, method = "bonferroni", n = 4)  ###0.0008
ES(diversity[sb=="SG"],diversity[sb=="SH"])  ###0.104

boots(response=diversity,lower="OS") #0.011
p.adjust(0.011, method = "bonferroni", n = 4)  ###0.044
ES(diversity[sb=="SG"],diversity[sb=="OS"])  ###0.074

boots(response=diversity,lower="DT",higher="GR") #0.1855
p.adjust(0.1855, method = "bonferroni", n = 4)  ###0.742
ES(diversity[sb=="OS"],diversity[sb=="DT"])  ###0.0454


lm20=glm(func.div~year+season+river+site:river+sb+volume+fish+crabs+rays,family="gaussian"(link="identity"))
par(mfrow=c(2,2))
plot(lm20)

summary(lm20)
confint(lm20)

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                              0.15347    0.09843   1.559  0.12097    
# year2013                                -0.01818    0.03878  -0.469  0.63992    
# seasonSpring                             0.01356    0.05266   0.258  0.79711    
# seasonSummer                            -0.08551    0.05565  -1.537  0.12639    
# riverMobjack                            -0.02888    0.09757  -0.296  0.76763    
# riverYork                               -0.09564    0.10148  -0.942  0.34739    
# sbGR                                     0.04606    0.08952   0.515  0.60762    
# sbOS                                     0.01902    0.07736   0.246  0.80613    
# sbSG                                     0.34076    0.08371   4.071  7.4e-05 ***
# sbSH                                     0.05486    0.06988   0.785  0.43362    
# volume                                   0.01652    0.02220   0.744  0.45785    
# fish                                     0.08510    0.06199   1.373  0.17176    
# crabs                                    0.01452    0.07489   0.194  0.84647    
# rays                                     0.84207    0.26907   3.130  0.00209 ** 
# riverYork:siteKing's Creek              -0.03751    0.10622  -0.353  0.72449    
# riverLynnhaven:siteLinkhorn Bay          0.08893    0.09212   0.965  0.33580    
# riverLynnhaven:siteLynnhaven Bay         0.10831    0.08952   1.210  0.22812    
# riverMobjack:siteMobjack 1              -0.01250    0.08323  -0.150  0.88078    
# riverMobjack:siteMobjack 2              -0.08062    0.08313  -0.970  0.33365    
# riverMobjack:siteMobjack 3              -0.08063    0.08061  -1.000  0.31871    
# riverLynnhaven:sitePleasure House Creek -0.03221    0.08543  -0.377  0.70662    
# riverYork:sitePropotank                  0.09519    0.09815   0.970  0.33365    
# riverYork:sitePurtan                     0.16508    0.10139   1.628  0.10549    
#  
1-(lm20$deviance/lm20$null.deviance)
#[1] 0.2678447

#In comparison with detrital mud (the habitat with lowest diversity), seagrass had 
#0.18 - 0.51 units higher functional diversity index

# 2.5 %     97.5 %
# (Intercept)                             -0.03944920 0.34638574
# year2013                                -0.09419517 0.05783707
# seasonSpring                            -0.08964533 0.11676594
# seasonSummer                            -0.19458599 0.02355786
# riverMobjack                            -0.22010420 0.16234776
# riverYork                               -0.29453602 0.10325251
# sbGR                                    -0.12939129 0.22150629
# sbOS                                    -0.13261105 0.17064856
# sbSG                                     0.17669508 0.50482013
# sbSH                                    -0.08210694 0.19182648
# volume                                  -0.02699096 0.06003633
# fish                                    -0.03639392 0.20658651
# crabs                                   -0.13225917 0.16130788
# rays                                     0.31471584 1.36943271
# riverYork:siteKing's Creek              -0.24569504 0.17068160
# riverLynnhaven:siteLinkhorn Bay         -0.09160929 0.26947874
# riverLynnhaven:siteLynnhaven Bay        -0.06714248 0.28377226
# riverMobjack:siteMobjack 1              -0.17563252 0.15062674
# riverMobjack:siteMobjack 2              -0.24356128 0.08231953
# riverMobjack:siteMobjack 3              -0.23862274 0.07735865
# riverLynnhaven:sitePleasure House Creek -0.19964946 0.13522139
# riverYork:sitePropotank                 -0.09719203 0.28756676
# riverYork:sitePurtan                    -0.03363786 0.36379489



boots(response=func.div) #0.0001
p.adjust(0.0001, method = "bonferroni", n = 5)  ###0.0005
ES(func.div[sb=="SG"],func.div[sb=="DT"])  ###0.155

boots(response=func.div,lower="SH") #0.0001
p.adjust(0.0001, method = "bonferroni", n = 5)  ###0.0005
ES(func.div[sb=="SG"],func.div[sb=="SH"])  ###0.108

boots(response=func.div,lower="OS") #0.0012
p.adjust(0.0012, method = "bonferroni", n = 5)  ###0.006
ES(func.div[sb=="SG"],func.div[sb=="OS"])  ###0.104

boots(response=func.div,lower="GR") #0.0662
p.adjust(0.0662, method = "bonferroni", n = 5)  ###0.331
ES(func.div[sb=="SG"],func.div[sb=="GR"])  ###0.067

boots(response=func.div,higher="GR",lower="DT") #0.0402
p.adjust(0.0402, method = "bonferroni", n = 5)  ###0.201
ES(func.div[sb=="GR"],func.div[sb=="DT"])  ###0.157


#############################################################
png("diversity and sb.png",height=9,width=6.25,units="in",res=300)

par(mfrow=c(2,1))

means.diversity=c(mean(diversity[sb=="DT"]),mean(diversity[sb=="SH"]),mean(diversity[sb=="OS"]),mean(diversity[sb=="GR"]),mean(diversity[sb=="SG"]))
ses.diversity=c(sd(diversity[sb=="DT"])/sqrt(length(diversity[sb=="DT"])),sd(diversity[sb=="SH"])/sqrt(length(diversity[sb=="SH"])),sd(diversity[sb=="OS"])/sqrt(length(diversity[sb=="OS"])),sd(diversity[sb=="GR"])/sqrt(length(diversity[sb=="GR"])),sd(diversity[sb=="SG"])/sqrt(length(diversity[sb=="SG"])))
txt1=expression(paste("Mean diversity \u00b1 SE"))
txt2=expression(paste("a.                                                                                                                     "))
par(mar=c(4.0,6.0,3.0,1),mgp=c(3,1.5,0),cex=1)
onewayplot2=barplot(means.diversity,names.arg=c("Detrital\nmud","Shell\nhash","Oyster\nshell","Coarse\nsand","Seagrass"),space=c(.05,.05,.05,.05,.05),ylim=c(0,0.6),col=c("gray28","gray94","black","darkgrey","white"),border=c("black","black","black","black","black"),xlab="Habitat",ylab=txt1,cex.lab=1.4,cex.axis=1.2,cex.names=1.2,main=txt2,cex.main=1.2)
arrows(onewayplot2,means.diversity+ses.diversity,onewayplot2,means.diversity-ses.diversity,code=3,angle=90,length=0.05,lwd=1.5)
text(0.55, means.diversity[1]+ses.diversity[1]+0.03, "a",cex=1.2)
text(1.61, means.diversity[2]+ses.diversity[2]+0.03, "a",cex=1.2)
text(2.655, means.diversity[3]+ses.diversity[3]+0.03, "a",cex=1.2)
text(3.69, means.diversity[4]+ses.diversity[4]+0.03, "ab",cex=1.2)
text(4.755, means.diversity[5]+ses.diversity[5]+0.03, "b",cex=1.2)

means.func.div=c(mean(func.div[sb=="DT"]),mean(func.div[sb=="SH"]),mean(func.div[sb=="OS"]),mean(func.div[sb=="GR"]),mean(func.div[sb=="SG"]))
ses.func.div=c(sd(func.div[sb=="DT"])/sqrt(length(func.div[sb=="DT"])),sd(func.div[sb=="SH"])/sqrt(length(func.div[sb=="SH"])),sd(func.div[sb=="OS"])/sqrt(length(func.div[sb=="OS"])),sd(func.div[sb=="GR"])/sqrt(length(func.div[sb=="GR"])),sd(func.div[sb=="SG"])/sqrt(length(func.div[sb=="SG"])))
txt1=expression(paste("Mean functional \n  diversity \u00b1 SE"))
txt2=expression(paste("b.                                                                                                                     "))
par(mar=c(4.0,6.0,3.0,1),mgp=c(3,1.5,0),cex=1)
onewayplot2=barplot(means.func.div,names.arg=c("Detrital\nmud","Shell\nhash","Oyster\nshell","Coarse\nsand","Seagrass"),space=c(.05,.05,.05,.05,.05),ylim=c(0,0.60),col=c("gray28","gray94","black","darkgrey","white"),border=c("black","black","black","black","black"),xlab="Habitat",ylab=txt1,cex.lab=1.4,cex.axis=1.2,cex.names=1.2,main=txt2,cex.main=1.2)
arrows(onewayplot2,means.func.div+ses.func.div,onewayplot2,means.func.div-ses.func.div,code=3,angle=90,length=0.05,lwd=1)
text(0.55, means.func.div[1]+ses.func.div[1]+0.03, "a",cex=1.2)
text(1.61, means.func.div[2]+ses.func.div[2]+0.03, "a",cex=1.2)
text(2.655, means.func.div[3]+ses.func.div[3]+0.03, "a",cex=1.2)
text(3.69, means.func.div[4]+ses.func.div[4]+0.03, "ab",cex=1.2)
text(4.755, means.func.div[5]+ses.func.div[5]+0.03, "b",cex=1.2)

dev.off()

png("richness and sb.png",height=9,width=6.25,units="in",res=300)

par(mfrow=c(2,1))

means.rare.rich=c(mean(rare.rich[sb=="DT"]),mean(rare.rich[sb=="SH"]),mean(rare.rich[sb=="OS"]),mean(rare.rich[sb=="GR"]),mean(rare.rich[sb=="SG"]))
ses.rare.rich=c(sd(rare.rich[sb=="DT"])/sqrt(length(rare.rich[sb=="DT"])),sd(rare.rich[sb=="SH"])/sqrt(length(rare.rich[sb=="SH"])),sd(rare.rich[sb=="OS"])/sqrt(length(rare.rich[sb=="OS"])),sd(rare.rich[sb=="GR"])/sqrt(length(rare.rich[sb=="GR"])),sd(rare.rich[sb=="SG"])/sqrt(length(rare.rich[sb=="SG"])))
txt1=expression(paste("Mean rarefied \nrichness \u00b1 SE"))
txt2=expression(paste("a.                                                                                                                     "))
par(mar=c(4.0,6.0,3.0,1),mgp=c(3,1.5,0),cex=1)
onewayplot2=barplot(means.rare.rich,names.arg=c("Detrital\nmud","Shell\nhash","Oyster\nshell","Coarse\nsand","Seagrass"),space=c(.05,.05,.05,.05,.05),ylim=c(0,12),col=c("gray28","gray94","black","darkgrey","white"),border=c("black","black","black","black","black"),xlab="Habitat",ylab=txt1,cex.lab=1.4,cex.axis=1.2,cex.names=1.2,main=txt2,cex.main=1.2)
arrows(onewayplot2,means.rare.rich+ses.rare.rich,onewayplot2,means.rare.rich-ses.rare.rich,code=3,angle=90,length=0.05,lwd=1.5)
text(0.55, means.rare.rich[1]+ses.rare.rich[1]+0.8, "a",cex=1.2)
text(1.61, means.rare.rich[2]+ses.rare.rich[2]+0.8, "ab",cex=1.2)
text(2.655, means.rare.rich[3]+ses.rare.rich[3]+0.8, "ab",cex=1.2)
text(3.69, means.rare.rich[4]+ses.rare.rich[4]+0.8, "ab",cex=1.2)
text(4.755, means.rare.rich[5]+ses.rare.rich[5]+0.8, "b",cex=1.2)

means.func.rich=c(mean(func.rich[sb=="DT"]),mean(func.rich[sb=="SH"]),mean(func.rich[sb=="OS"]),mean(func.rich[sb=="GR"]),mean(func.rich[sb=="SG"]))
ses.func.rich=c(sd(func.rich[sb=="DT"])/sqrt(length(func.rich[sb=="DT"])),sd(func.rich[sb=="SH"])/sqrt(length(func.rich[sb=="SH"])),sd(func.rich[sb=="OS"])/sqrt(length(func.rich[sb=="OS"])),sd(func.rich[sb=="GR"])/sqrt(length(func.rich[sb=="GR"])),sd(func.rich[sb=="SG"])/sqrt(length(func.rich[sb=="SG"])))
txt1=expression(paste("Mean functional \n  richness \u00b1 SE"))
txt2=expression(paste("b.                                                                                                                        "))
par(mar=c(4.0,6.0,3.0,1),mgp=c(3,1.5,0),cex=1)
onewayplot2=barplot(means.func.rich,names.arg=c("Detrital\nmud","Shell\nhash","Oyster\nshell","Coarse\nsand","Seagrass"),space=c(.05,.05,.05,.05,.05),ylim=c(0,3.5),col=c("gray28","gray94","black","darkgrey","white"),border=c("black","black","black","black","black"),xlab="Habitat",ylab=txt1,cex.lab=1.4,cex.axis=1.2,cex.names=1.2,main=txt2,cex.main=1.2)
arrows(onewayplot2,means.func.rich+ses.func.rich,onewayplot2,means.func.rich-ses.func.rich,code=3,angle=90,length=0.05,lwd=1)
text(0.55, means.func.rich[1]+ses.func.rich[1]+0.2, "a",cex=1.2)
text(1.61, means.func.rich[2]+ses.func.rich[2]+0.2, "a",cex=1.2)
text(2.655, means.func.rich[3]+ses.func.rich[3]+0.2, "a",cex=1.2)
text(3.69, means.func.rich[4]+ses.func.rich[4]+0.2, "ab",cex=1.2)
text(4.755, means.func.rich[5]+ses.func.rich[5]+0.2, "b",cex=1.2)

dev.off()

##############################################################

#Function to graph linear model output

mod.lm=function(ind,dep,col1,col2,col3,x.lab,y.lab,ti)
{
  x=match(ind,names(mya))
  y=match(dep,names(mya))
  mya.cut=mya[,-x] 
  names(mya.cut)
  dependent=mya.cut[,y]
  var1=as.factor(mya.cut[,6])
  var2=mya.cut[,5]
  var3=mya.cut[,2]
  var5=mya.cut[,20]
  var6=mya.cut[,26]
  var7=mya.cut[,27]
  var8=mya.cut[,28]
  independent=mya[,x]
  
  lm.dep=lm(dependent~var1+var2+var3+var5+
              var6+var7+var8+independent)
  
  new.seq=seq(min(na.omit(independent)), max(na.omit(independent)), length.out=1000)
  newdat = expand.grid(independent = new.seq, var1=c("2012","2013"),var2=c("Fall","Spring","Summer"),
                       var3=c("Mobjack","York","Lynnhaven"),var5=c("DT","GR","OS","SG","SH"),
                       var6=mean(na.omit(var6)),var7=mean(na.omit(var7)),var8=mean(na.omit(var8)))
  pred.mat=data.frame()
  lower.mat=data.frame()
  upper.mat=data.frame()
  ys=vector(length=1000)
  
  for(i in 1:2){
    if (i==1) {
      year.string="2012"
    } else {
      year.string="2013"
    }
    for(j in 1:3){
      if(j==1) {
        season.string="Spring"
      } else {
        if(j==2) {
          season.string="Summer"
        } else {
          season.string="Fall"
        }
      }
      for(k in 1:3){
        if(k==1) {
          river.string="Lynnhaven"
        } else {
          if(k==2) {
            river.string="Mobjack"
          } else {
            river.string="York"
          }
        }
        for(l in 1:5){
          if(l==1) {
            sb.string="DT"
          } else { 
          if(l==2) {
            sb.string="GR"
          } else {
            if(l==3) {
              sb.string="OS"
            } else {
              if(l==4) {
                sb.string="SG"
              } else {
                sb.string="SH"
              }
            }
           }
          }
          new.row.y=predict(lm.dep, newdata=newdat[newdat$var1==year.string&newdat$var2==season.string&newdat$var3==river.string&newdat$var5==sb.string,],interval="confidence",level=0.95)[,1]
          new.row.upper=predict(lm.dep, newdata=newdat[newdat$var1==year.string&newdat$var2==season.string&newdat$var3==river.string&newdat$var5==sb.string,],interval="confidence",level=0.95)[,3]
          new.row.lower=predict(lm.dep, newdata=newdat[newdat$var1==year.string&newdat$var2==season.string&newdat$var3==river.string&newdat$var5==sb.string,],interval="confidence",level=0.95)[,2]
          
          pred.mat=rbind(pred.mat,unname(new.row.y))
          lower.mat=rbind(lower.mat,unname(new.row.lower))
          upper.mat=rbind(upper.mat,unname(new.row.upper))
          
          
        }
      }
    }
  }
  ys=apply(pred.mat,2,mean)
  lower=apply(lower.mat,2,mean)
  upper=apply(upper.mat,2,mean)
  # png(filename=paste(dep,ind,".png",sep=""),height=5,width=6,unit="in",res=300)
  # plot(mya[,x],mya[,y],col="white",xlab=x.lab,ylab=y.lab,main=ti)
  # polygon(c(new.seq , rev(new.seq)), c(upper, rev(lower)), col = col3, border = NA)
  # lines(new.seq,ys,lwd=2,col=col2)
  # points(jitter(mya[,x]),jitter(mya[,y]),pch=1,col=col1)
  # dev.off()
  plot(mya[,x],mya[,y],col="white",xlab=x.lab,ylab=y.lab,main=ti)
  polygon(c(new.seq , rev(new.seq)), c(upper, rev(lower)), col = col3, border = NA)
  lines(new.seq,ys,lwd=2,col=col2)
  points(jitter(mya[,x]),jitter(mya[,y]),pch=1,col=col1)
}

#Function to graph glm output

mod.glm=function(ind,dep,linkfunc="quasipoisson",correct=1,col1,col2,col3,x.lab,y.lab,ti)
{
  x=match(ind,names(mya))
  y=match(dep,names(mya))
  mya.cut=mya[,-x] 
  names(mya.cut)
  dependent=mya.cut[,y]
  var1=as.factor(mya.cut[,6])
  var2=mya.cut[,5]
  var3=mya.cut[,2]
  var5=mya.cut[,20]
  var6=mya.cut[,26]
  var7=mya.cut[,27]
  var8=mya.cut[,28]
  independent=mya[,x]
  
  glm.dep=glm(dependent~var1+var2+var3+var5+
              var6+var7+var8+independent,family=linkfunc)
  
  new.seq=seq(min(na.omit(independent)), max(na.omit(independent)), length.out=1000)
  newdat = expand.grid(independent = new.seq, var1=c("2012","2013"),var2=c("Fall","Spring","Summer"),
                       var3=c("Mobjack","York","Lynnhaven"),var5=c("DT","GR","OS","SG","SH"),
                       var6=mean(na.omit(var6)),var7=mean(na.omit(var7)),var8=mean(na.omit(var8)))
  pred.mat=data.frame()
  lower.mat=data.frame()
  upper.mat=data.frame()
  ys=vector(length=1000)
  
  for(i in 1:2){
    if (i==1) {
      year.string="2012"
    } else {
      year.string="2013"
    }
    for(j in 1:3){
      if(j==1) {
        season.string="Spring"
      } else {
        if(j==2) {
          season.string="Summer"
        } else {
          season.string="Fall"
        }
      }
      for(k in 1:3){
        if(k==1) {
          river.string="Lynnhaven"
        } else {
          if(k==2) {
            river.string="Mobjack"
          } else {
            river.string="York"
          }
        }
        for(l in 1:5){
          if(l==1) {
            sb.string="DT"
          } else { 
          if(l==2) {
            sb.string="GR"
          } else {
            if(l==3) {
              sb.string="OS"
            } else {
              if(l==4) {
                sb.string="SG"
              } else {
                sb.string="SH"
              }
            }
          }
          }
          new.row.y=predict(glm.dep, newdata=newdat[newdat$var1==year.string&newdat$var2==season.string&newdat$var3==river.string&newdat$var5==sb.string,],se.fit=T)$fit
          new.row.upper=new.row.y+(qnorm(0.975)*predict(glm.dep, newdata=newdat[newdat$var1==year.string&newdat$var2==season.string&newdat$var3==river.string&newdat$var5==sb.string,],se.fit=T)$se.fit)
          new.row.lower=new.row.y-(qnorm(0.975)*predict(glm.dep, newdata=newdat[newdat$var1==year.string&newdat$var2==season.string&newdat$var3==river.string&newdat$var5==sb.string,],se.fit=T)$se.fit)
          
          pred.mat=rbind(pred.mat,exp(unname(new.row.y)))
          lower.mat=rbind(lower.mat,exp(unname(new.row.lower)))
          upper.mat=rbind(upper.mat,exp(unname(new.row.upper)))
        }
      }
    }
  }
  
  
  ys=apply(pred.mat,2,mean) 
  lower=apply(lower.mat,2,mean)
  upper=apply(upper.mat,2,mean)
  # png(filename=paste(dep,ind,".png",sep=""),height=5,width=6,unit="in",res=300)
  # plot(mya[,x],mya[,y]/correct,col="white",xlab=x.lab,ylab=y.lab,main=ti)
  # polygon(c(new.seq , rev(new.seq)), c(upper/correct, rev(lower)/correct), col = col3, border = NA)
  # lines(new.seq,ys/correct,lwd=2,col=col2)
  # points(jitter(mya[,x]),jitter(mya[,y])/correct,pch=1,col=col1)
  # dev.off()
  plot(mya[,x],mya[,y]/correct,col="white",xlab=x.lab,ylab=y.lab,main=ti)
  polygon(c(new.seq , rev(new.seq)), c(upper/correct, rev(lower)/correct), col = col3, border = NA)
  lines(new.seq,ys/correct,lwd=2,col=col2)
  points(jitter(mya[,x]),jitter(mya[,y])/correct,pch=1,col=col1)
}

png(filename="scatter.plots.png",height=4*3.5,width=4*4,unit="in",res=300)
par(mfrow=c(4,4),mar=c(3,4,2,1),mgp=c(2,0.6,0),cex=1.3)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
mtext(text="Diversity Indices",side=1,line=3,cex=par()$cex+0.3)
mod.lm("Rays","Diversity",col1="black",col2="black",col3="darkgrey",x.lab=expression(paste("Number of Ray Pits m"^"-2")),y.lab="Gini-Simpson Diversity",ti=expression(paste("a.                                       ")))
mod.lm("Rays","Func.div",col1="black",col2="black",col3="darkgrey",x.lab=expression(paste("Number of Ray Pits m"^"-2")),y.lab="Functional Diversity Index",ti=expression(paste("b.                                       ")))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
mod.lm("Volume","Rare.rich",col1="black",col2="black",col3="darkgrey",x.lab=expression(paste("Habitat Volume (mL)")),y.lab="Rarefied Species Richness",ti=expression(paste("c.                                       ")))
mod.lm("Volume","Func.rich",col1="black",col2="black",col3="darkgrey",x.lab="Habitat Volume (mL)",y.lab="Functional Richness",ti=expression(paste("d.                                       ")))
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
mtext("Functional Group\nDensities",side=1,line=3,cex=par()$cex+0.3)
mod.glm("Crabs","Deposit.feeders",col1="black",col2="black",col3="darkgrey",correct=0.11,x.lab=expression(paste("Number of Blue Crabs m"^"-2")),y.lab="Total DF Bivalve Density",ti=expression(paste("e.                                       ")))
mod.glm("Rays","Deposit.feeders",col1="black",col2="black",col3="darkgrey",correct=0.11,x.lab=expression(paste("Number of Ray Pits m"^"-2")),y.lab="Total DF Bivalve Density",ti=expression(paste("f.                                       ")))
mod.glm("Volume","DBSF",col1="black",col2="black",col3="darkgrey",correct=0.11,x.lab="Habitat Volume (mL)",y.lab="Total DBSF Bivalve Density",ti=expression(paste("g.                                       ")))

plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
mod.glm("Volume","TSSD",col1="black",col2="black",col3="darkgrey",correct=0.11,x.lab="Habitat Volume (mL)",y.lab="Total TSSD Bivalve Density",ti=expression(paste("h.                                       ")))
mod.glm("Crabs","TSSD",col1="black",col2="black",col3="darkgrey",correct=0.11,x.lab=expression(paste("Number of Blue Crabs m"^"-2")),y.lab="Total TSSD Bivalve Density",ti=expression(paste("i.                                       ")))
mod.glm("Volume","Hard.shelled",col1="black",col2="black",col3="darkgrey",correct=0.11,x.lab="Habitat Volume (mL)",y.lab="Total ARM Bivalve Density",ti=expression(paste("j.                                       ")))

dev.off()


###############################################################################################################

#Analysis of functional groups

lm18=glm(deposit*0.11~year+season+river+site:river+sb+volume+fish+crabs+rays,family="quasipoisson")
par(mfrow=c(2,2))
plot(lm18)

#Examining if a zero inflated model is a better fit:
library(pscl)

fm1=glm(deposit*0.11~year+season+river+site:river+sb+volume+fish+crabs+rays,family="quasipoisson")
## set up version of data with non-identified regressors omitted 
quine3 <- as.data.frame(model.matrix(fm1)) ## all regressors 
quine3 <- quine3[, !is.na(coef(fm1))]      ## only identified 
quine3 <- quine3[, -1]                     ## omit intercept 
quine3$deposit <- deposit*0.11             ## add response 
# tail(quine3)
# head(quine3)
fm1a <- glm(deposit~as.factor(year2013)+as.factor(seasonSpring)+as.factor(seasonSummer)+as.factor(riverMobjack)+as.factor(riverYork)+as.factor(sbGR)+as.factor(sbOS)+as.factor(sbSG)+as.factor(sbSH)+as.factor(quine3[,14])+as.factor(quine3[,21])+as.factor(quine3[,22])+volume+fish+crabs+rays, data = quine3, family="quasipoisson") 
fm2a <- glm(deposit~as.factor(year2013)+as.factor(seasonSpring)+as.factor(seasonSummer)+as.factor(riverMobjack)+as.factor(riverYork)+as.factor(sbGR)+as.factor(sbOS)+as.factor(sbSG)+as.factor(sbSH)+as.factor(quine3[,14])+as.factor(quine3[,21])+as.factor(quine3[,22])+volume+fish+crabs+rays, data = quine3, family="poisson") 
fm3a <- zeroinfl(deposit~as.factor(year2013)+as.factor(seasonSpring)+as.factor(seasonSummer)+as.factor(riverMobjack)+as.factor(riverYork)+as.factor(sbGR)+as.factor(sbOS)+as.factor(sbSG)+as.factor(sbSH)+as.factor(quine3[,14])+as.factor(quine3[,22])+volume+fish+crabs+rays, data = quine3, dist = "poisson") 
#Comparing quasipoisson (reduced) to zeroinflated
plot(fitted(fm1a), fitted(fm3a))
boxplot(resid(fm1a) - resid(fm3a))

boxplot(resid(fm1a))
boxplot(resid(fm3a))

#comparing quasipoisson (with all nested terms) to zeroinflated

plot(fitted(fm3a), fitted(fm1))
boxplot(resid(fm3a) - resid(fm1))

boxplot(resid(fm3a))
boxplot(resid(fm1))

#comparing quasipoisson (reduced) to quasipoisson (with all nested terms)

plot(fitted(fm1a), fitted(fm1))
boxplot(resid(fm1a) - resid(fm1))

boxplot(resid(fm1a))
boxplot(resid(fm1))

#Quasipoisson with all nested terms reduces residuals and has similar fitted values to the rest of the possible models

summary(lm18)
exp(confint(lm18))

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                              -0.5131     1.1071  -0.464 0.643636    
# year2013                                  0.7413     0.2408   3.079 0.002454 ** 
# seasonSpring                              1.2226     0.4437   2.755 0.006554 ** 
# seasonSummer                              0.2469     0.5102   0.484 0.629202    
# riverMobjack                             -0.7311     1.4968  -0.488 0.625925    
# riverYork                                -1.1120     1.5586  -0.713 0.476614    
# sbGR                                     -0.7593     0.4362  -1.741 0.083716 .  
# sbOS                                     -1.0317     0.3990  -2.586 0.010620 *  
# sbSG                                      1.0236     1.0112   1.012 0.312967    
# sbSH                                     -0.3056     0.3092  -0.989 0.324422    
# volume                                   -0.1993     0.1462  -1.364 0.174671    
# fish                                      0.4275     0.3062   1.396 0.164645    
# crabs                                    -2.0900     0.6036  -3.463 0.000689 ***
# rays                                     -5.0191     2.0078  -2.500 0.013454 *  
# riverYork:siteKing's Creek                5.2987     1.2114   4.374 2.22e-05 ***
# riverLynnhaven:siteLinkhorn Bay           1.6619     1.1916   1.395 0.165092    
# riverLynnhaven:siteLynnhaven Bay          0.3919     1.2334   0.318 0.751141    
# riverMobjack:siteMobjack 1                0.2922     0.9242   0.316 0.752246    
# riverMobjack:siteMobjack 2                0.6443     0.8104   0.795 0.427786    
# riverMobjack:siteMobjack 3                0.2483     0.8897   0.279 0.780575    
# riverLynnhaven:sitePleasure House Creek   1.0828     1.1295   0.959 0.339182    
# riverYork:sitePropotank                   4.5924     1.2054   3.810 0.000199 ***
# riverYork:sitePurtan                      3.9711     1.2141   3.271 0.001318 ** 


1-(lm18$deviance/lm18$null.deviance)
#0.7902575

# 2.5 %       97.5 %
# (Intercept)                             3.115996e-02    3.6073613
# year2013                                1.320668e+00    3.4059089
# seasonSpring                            1.482737e+00    8.6311870
# seasonSummer                            4.802599e-01    3.6251619
# riverMobjack                            2.009657e-02   13.5045932
# riverYork                               1.217628e-02    9.9794169
# sbGR                                    1.935881e-01    1.0820626
# sbOS                                    1.621429e-01    0.7771115
# sbSG                                    5.106657e-01   37.7359224
# sbSH                                    4.046357e-01    1.3645621
# volume                                  6.064009e-01    1.0799508
# fish                                    8.258696e-01    2.7596169
# crabs                                   3.177632e-02    0.3578556
# rays                                    6.464747e-05    0.2046574
# riverYork:siteKing's Creek              2.571440e+01 3697.7277612
# riverLynnhaven:siteLinkhorn Bay         6.454104e-01  109.5549278
# riverLynnhaven:siteLynnhaven Bay        1.361263e-01   31.5460925
# riverMobjack:siteMobjack 1              1.854879e-01    8.7707740
# riverMobjack:siteMobjack 2              3.864129e-01   10.9213383
# riverMobjack:siteMobjack 3              1.993093e-01    8.0612730
# riverLynnhaven:sitePleasure House Creek 4.191452e-01   57.3786801
# riverYork:sitePropotank                 1.290097e+01 1810.9400319
# riverYork:sitePurtan                    6.768260e+00  983.5580061

#Detrital mud supported higher densities of the facultative deposit feeders (DF)
#Limecola balthica and Ameritella mitchelli than seagrass (p = 0.0006, d = 0.15; Fig. 4a). 

boots(response=deposit, higher="DT",lower="SG") #0.0002
p.adjust(0.0002, method = "bonferroni", n = 3)  ###0.0006
ES(deposit[sb=="DT"],deposit[sb=="SG"])  ###0.15

boots(response=deposit, higher="DT",lower="GR") #0.0883
p.adjust(0.0883, method = "bonferroni", n = 3)  ###0.2649
ES(deposit[sb=="DT"],deposit[sb=="OS"])  ###0.084

boots(response=deposit, higher="SH",lower="SG") #0.0341
p.adjust(0.0341, method = "bonferroni", n = 3)  ###0.1023
ES(deposit[sb=="SH"],deposit[sb=="SG"])  ###0.033

lm14=glm(DBSF*0.11~year+season+river+site:river+sb+volume+fish+crabs+rays,family="quasipoisson")
par(mfrow=c(2,2))
plot(lm14)
summary(lm14) 
exp(confint(lm14))

# Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                              0.58761    0.62175   0.945   0.3461  
# year2013                                -0.18391    0.24247  -0.758   0.4493  
# seasonSpring                             0.22706    0.31489   0.721   0.4719  
# seasonSummer                            -0.09561    0.35830  -0.267   0.7899  
# riverMobjack                            -0.41023    0.57870  -0.709   0.4794  
# riverYork                                0.28450    0.53725   0.530   0.5972  
# sbGR                                     0.18986    0.61348   0.309   0.7574  
# sbOS                                     0.34083    0.53973   0.631   0.5286  
# sbSG                                     0.83431    0.57042   1.463   0.1456  
# sbSH                                     0.93132    0.47318   1.968   0.0508 .
# volume                                   0.28157    0.12460   2.260   0.0252 *
# fish                                    -0.08305    0.38651  -0.215   0.8301  
# crabs                                   -0.58262    0.54480  -1.069   0.2865  
# rays                                    -0.96224    1.86304  -0.516   0.6062  
# riverYork:siteKing's Creek              -1.55939    0.72768  -2.143   0.0337 *
# riverLynnhaven:siteLinkhorn Bay         -0.16996    0.54130  -0.314   0.7540  
# riverLynnhaven:siteLynnhaven Bay         0.14126    0.47055   0.300   0.7644  
# riverMobjack:siteMobjack 1               0.15346    0.50949   0.301   0.7637  
# riverMobjack:siteMobjack 2              -0.26532    0.56214  -0.472   0.6376  
# riverMobjack:siteMobjack 3              -0.29748    0.54836  -0.542   0.5882  
# riverLynnhaven:sitePleasure House Creek -0.43660    0.53614  -0.814   0.4167  
# riverYork:sitePropotank                 -0.07171    0.50973  -0.141   0.8883  
# riverYork:sitePurtan                    -1.49092    0.75187  -1.983   0.0491 *

1-(lm14$deviance/lm14$null.deviance)
#0.1926388

# 2.5 %     97.5 %
# (Intercept)                             0.498946759  5.7714396
# year2013                                0.515311053  1.3376389
# seasonSpring                            0.681262695  2.3538590
# seasonSummer                            0.451114204  1.8470465
# riverMobjack                            0.205469233  2.0249612
# riverYork                               0.456526275  3.8080855
# sbGR                                    0.354475523  4.0781144
# sbOS                                    0.504514515  4.2738640
# sbSG                                    0.794676653  7.5410871
# sbSH                                    1.055825304  6.8869211
# volume                                  1.029574097  1.6836083
# fish                                    0.404908485  1.8649141
# crabs                                   0.170507621  1.4891750
# rays                                    0.007644229 12.7591541
# riverYork:siteKing's Creek              0.044629082  0.8135639
# riverLynnhaven:siteLinkhorn Bay         0.285063183  2.4215809
# riverLynnhaven:siteLynnhaven Bay        0.454425199  2.9203158
# riverMobjack:siteMobjack 1              0.422210305  3.2221977
# riverMobjack:siteMobjack 2              0.240032934  2.2861508
# riverMobjack:siteMobjack 3              0.239415518  2.1606668
# riverLynnhaven:sitePleasure House Creek 0.210378184  1.7891291
# riverYork:sitePropotank                 0.342534596  2.5681377
# riverYork:sitePurtan                    0.042183652  0.8894158


#Deep-burrowing suspension-feeding (DBSF) bivalves, such as Tagelus plebeius, 
#Ensis directus, Mya arenaria, Petricolaria pholadiformis, and Tagelus divisus 
#(Table 1), had similar densities in all habitats (Fig. 4b).

boots(response=DBSF,higher="SH",lower="DT") #0.0848
ES(DBSF[sb=="SH"],DBSF[sb=="DT"])  ###0.04


lm17=glm(TSSB*0.11~year+season+river+site:river+sb+volume+fish+crabs+rays,family="quasipoisson")
par(mfrow=c(2,2))
plot(lm17)

#Examining if a zero inflated model is a better fit:

fm1=glm(TSSB*0.11~year+season+river+site:river+sb+volume+fish+crabs+rays,family="quasipoisson")
## set up version of data with non-identified regressors omitted 
quine3 <- as.data.frame(model.matrix(fm1)) ## all regressors 
quine3 <- quine3[, !is.na(coef(fm1))]      ## only identified 
quine3 <- quine3[, -1]                     ## omit intercept 
quine3$TSSB <- TSSB*0.11             ## add response 
# tail(quine3)
# head(quine3)
fm1a <- glm(TSSB~as.factor(year2013)+as.factor(seasonSpring)+as.factor(seasonSummer)+as.factor(riverMobjack)+as.factor(riverYork)+as.factor(sbGR)+as.factor(sbOS)+as.factor(sbSG)+as.factor(sbSH)+as.factor(quine3[,14])+as.factor(quine3[,21])+as.factor(quine3[,22])+volume+fish+crabs+rays, data = quine3, family="quasipoisson") 
fm2a <- glm(TSSB~as.factor(year2013)+as.factor(seasonSpring)+as.factor(seasonSummer)+as.factor(riverMobjack)+as.factor(riverYork)+as.factor(sbGR)+as.factor(sbOS)+as.factor(sbSG)+as.factor(sbSH)+as.factor(quine3[,14])+as.factor(quine3[,21])+as.factor(quine3[,22])+volume+fish+crabs+rays, data = quine3, family="poisson") 
fm3a <- zeroinfl(TSSB~as.factor(year2013)+as.factor(seasonSpring)+as.factor(seasonSummer)+as.factor(riverMobjack)+as.factor(riverYork)+as.factor(sbGR)+as.factor(sbOS)+as.factor(sbSG)+as.factor(sbSH)+as.factor(quine3[,14])+as.factor(quine3[,22])+volume+fish+crabs+rays, data = quine3, dist = "poisson") 

#Comparing quasipoisson (reduced) to zeroinflated
plot(fitted(fm1a), fitted(fm3a))
boxplot(resid(fm1a) - resid(fm3a))

boxplot(resid(fm1a))
boxplot(resid(fm3a))

#comparing quasipoisson (with all nested terms) to zeroinflated

plot(fitted(fm3a), fitted(fm1))
boxplot(resid(fm3a) - resid(fm1))

boxplot(resid(fm3a))
boxplot(resid(fm1))

#comparing quasipoisson (reduced) to quasipoisson (with all nested terms)

plot(fitted(fm1a), fitted(fm1))
boxplot(resid(fm1a) - resid(fm1))

boxplot(resid(fm1a))
boxplot(resid(fm1))

#Quasipoisson with all nested terms reduces residuals and has similar fitted values to the rest of the possible models


summary(lm17)
exp(confint(lm17))

# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                              1.31905    0.97199   1.357  0.17671   
# year2013                                 0.13012    0.33882   0.384  0.70148   
# seasonSpring                             0.02326    0.49442   0.047  0.96254   
# seasonSummer                             1.19232    0.58285   2.046  0.04246 * 
# riverMobjack                            -1.53616    0.82714  -1.857  0.06515 . 
# riverYork                               -2.01355    0.98130  -2.052  0.04184 * 
# sbGR                                     0.34759    0.94599   0.367  0.71379   
# sbOS                                     0.65287    0.85268   0.766  0.44503   
# sbSG                                     0.67815    0.82236   0.825  0.41083   
# sbSH                                     0.18265    0.79397   0.230  0.81835   
# volume                                   0.43480    0.22120   1.966  0.05110 . 
# fish                                    -1.02069    0.69599  -1.467  0.14450   
# crabs                                   -1.89825    1.02222  -1.857  0.06518 . 
# rays                                    -5.00153    3.09805  -1.614  0.10845   
# riverYork:siteKing's Creek              -1.19221    1.59673  -0.747  0.45639   
# riverLynnhaven:siteLinkhorn Bay         -2.56134    0.98433  -2.602  0.01015 * 
# riverLynnhaven:siteLynnhaven Bay        -2.73264    0.98616  -2.771  0.00626 **
# riverMobjack:siteMobjack 1               1.18845    0.64583   1.840  0.06763 . 
# riverMobjack:siteMobjack 2               1.19429    0.61348   1.947  0.05335 . 
# riverMobjack:siteMobjack 3               1.45889    0.54364   2.684  0.00807 **
# riverLynnhaven:sitePleasure House Creek -0.22878    0.62929  -0.364  0.71669   
# riverYork:sitePropotank                 -2.94078    3.13870  -0.937  0.35023   
# riverYork:sitePurtan                    -0.41879    1.32018  -0.317  0.75150   

 
1-(lm17$deviance/lm17$null.deviance)
#0.4731316

# 2.5 %     97.5 %
# (Intercept)                             4.661217e-01 22.3527103
# year2013                                5.796150e-01  2.2080802
# seasonSpring                            3.964661e-01  2.8223967
# seasonSummer                            1.093420e+00 10.8875644
# riverMobjack                            3.938130e-02  1.0525181
# riverYork                               1.603283e-02  0.8275789
# sbGR                                    2.308888e-01 10.3953635
# sbOS                                    4.129221e-01 12.6852577
# sbSG                                    4.676965e-01 12.6567840
# sbSH                                    2.956229e-01  7.2023763
# volume                                  9.970597e-01  2.4209277
# fish                                    7.824324e-02  1.2006792
# crabs                                   1.339514e-02  0.8024106
# rays                                    1.068875e-05  2.1932813
# riverYork:siteKing's Creek              6.019547e-03  5.3994668
# riverLynnhaven:siteLinkhorn Bay         8.582028e-03  0.4360440
# riverLynnhaven:siteLynnhaven Bay        6.350352e-03  0.3727603
# riverMobjack:siteMobjack 1              9.573665e-01 12.7744270
# riverMobjack:siteMobjack 2              1.039924e+00 12.2292910
# riverMobjack:siteMobjack 3              1.618452e+00 14.3976923
# riverLynnhaven:sitePleasure House Creek 2.253002e-01  2.7256339
# riverYork:sitePropotank                           NA  2.6738967
# riverYork:sitePurtan                    3.441701e-02  8.4170259

#Thin-shelled and surface-dwelling (TSSD) bivalves, such as Kelliopsis 
#elevata and Arcuatula papyria (Table 1), had higher densities in seagrass 
#habitat than detrital mud (p = 0.05, d = 0.06) or shell hash 
#(p = 0.04, d = 0.04, Fig. 4c). 

boots(response=TSSB) #0.0159
p.adjust(0.0159, method = "bonferroni", n = 3)  ###0.0477
ES(TSSB[sb=="SG"],TSSB[sb=="DT"])  ###0.057

boots(response=TSSB,lower="SH") #0.0129
p.adjust(0.0129, method = "bonferroni", n = 3)  ###0.0387
ES(TSSB[sb=="SG"],TSSB[sb=="SH"])  ###0.044

boots(response=TSSB,higher="OS") #0.1945
p.adjust(0.1945, method = "bonferroni", n = 3)  ###0.5835
ES(TSSB[sb=="SG"],TSSB[sb=="GR"])  ###0.0259


lm16=glm(hard*0.11~year+season+river+site:river+sb+volume+fish+crabs+rays,family="quasipoisson")
par(mfrow=c(2,2))
plot(lm16)
summary(lm16)
exp(confint(lm16))

# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)                              -2.3454     0.9508  -2.467  0.01471 * 
# year2013                                  0.8604     0.3075   2.798  0.00579 **
# seasonSpring                             -0.5607     0.3849  -1.457  0.14713   
# seasonSummer                             -0.7722     0.4236  -1.823  0.07019 . 
# riverMobjack                              0.7662     0.8451   0.907  0.36596   
# riverYork                                 1.0048     0.8550   1.175  0.24170   
# sbGR                                      0.6979     0.8354   0.835  0.40472   
# sbOS                                      1.9880     0.6499   3.059  0.00261 **
# sbSG                                      1.5694     0.7526   2.085  0.03866 * 
# sbSH                                      1.0338     0.6442   1.605  0.11057   
# volume                                    0.3999     0.1530   2.613  0.00985 **
# fish                                      0.5927     0.4127   1.436  0.15297   
# crabs                                     0.1022     0.5040   0.203  0.83962   
# rays                                      3.4242     2.1695   1.578  0.11650   
# riverYork:siteKing's Creek               -2.9863     1.2022  -2.484  0.01404 * 
# riverLynnhaven:siteLinkhorn Bay           0.1695     0.7939   0.213  0.83126   
# riverLynnhaven:siteLynnhaven Bay         -0.4013     0.9196  -0.436  0.66319   
# riverMobjack:siteMobjack 1               -0.3415     0.5979  -0.571  0.56874   
# riverMobjack:siteMobjack 2               -1.5131     0.8699  -1.739  0.08392 . 
# riverMobjack:siteMobjack 3                0.3781     0.5211   0.726  0.46916   
# riverLynnhaven:sitePleasure House Creek   0.4381     0.7611   0.576  0.56566   
# riverYork:sitePropotank                   0.6406     0.6282   1.020  0.30939   
# riverYork:sitePurtan                     -0.3107     0.8053  -0.386  0.70012   

1-(lm16$deviance/lm16$null.deviance)
#0.3584024

# 2.5 %       97.5 %
# (Intercept)                             0.012459178    0.5374612
# year2013                                1.309455648    4.3972589
# seasonSpring                            0.266543974    1.2157598
# seasonSummer                            0.199772363    1.0590811
# riverMobjack                            0.424077898   12.6687315
# riverYork                               0.535383810   16.4767135
# sbGR                                    0.367910877   10.6589820
# sbOS                                    2.243370580   29.8910621
# sbSG                                    1.214345049   23.8964686
# sbSH                                    0.864417239   11.3249998
# volume                                  1.095815220    2.0051013
# fish                                    0.764320821    3.9241453
# crabs                                   0.384197559    2.8331136
# rays                                    0.479254344 2799.3084062
# riverYork:siteKing's Creek              0.003399163    0.4177504
# riverLynnhaven:siteLinkhorn Bay         0.265249989    6.4376795
# riverLynnhaven:siteLynnhaven Bay        0.099272098    4.2572112
# riverMobjack:siteMobjack 1              0.210679944    2.3159135
# riverMobjack:siteMobjack 2              0.027277499    1.0203960
# riverMobjack:siteMobjack 3              0.527327058    4.2449220
# riverLynnhaven:sitePleasure House Creek 0.366365778    7.9406778
# riverYork:sitePropotank                 0.573184677    6.8696551
# riverYork:sitePurtan                    0.143841686    3.5250855

#Armored (ARM) bivalves, such as Geukensia demissa, Mercenaria mercenaria, 
#Mulinia lateralis, Noetia ponderosa, and Modiolus modiolus (Table 1), 
#had the highest densities in oyster shell, and densities of ARM bivalves 
#were significantly higher in oyster shell than in detrital mud 
#(p = 0.04, d = 0.09) or shell hash (p = 0.03, d = 0.07; Fig. 4d). 



boots(response=hard, higher="OS",lower="SH") #0.0095
p.adjust(0.0095, method = "bonferroni", n = 3)  ###0.0285
ES(hard[sb=="OS"],hard[sb=="SH"])  ###0.070

boots(response=hard, higher="OS") #0.0138
p.adjust(0.0138, method = "bonferroni", n = 3)  ###0.0414
ES(hard[sb=="OS"],hard[sb=="DT"])  ###0.088

boots(response=hard) #0.0623
p.adjust(0.0623, method = "bonferroni", n = 3)  ###0.1869
ES(hard[sb=="SG"],hard[sb=="DT"])  ###0.044


#################################################Func div######################################################3
png("func groups and sb.png",height=9,width=12.5,units="in",res=300)
par(mfrow=c(2,2))

means.deposit=c(mean(deposit[sb=="DT"]),mean(deposit[sb=="SH"]),mean(deposit[sb=="OS"]),mean(deposit[sb=="GR"]),mean(deposit[sb=="SG"]))
ses.deposit=c(sd(deposit[sb=="DT"])/sqrt(length(deposit[sb=="DT"])),sd(deposit[sb=="SH"])/sqrt(length(deposit[sb=="SH"])),sd(deposit[sb=="OS"])/sqrt(length(deposit[sb=="OS"])),sd(deposit[sb=="GR"])/sqrt(length(deposit[sb=="GR"])),sd(deposit[sb=="SG"])/sqrt(length(deposit[sb=="SG"])))
txt1=expression(paste("Mean density DF\n  bivalves \u00b1 SE"))
txt2=expression(paste("a.                                                                                                                     "))
par(mar=c(4.0,6.0,3.0,1),mgp=c(3,1.5,0),cex=1)
onewayplot2=barplot(means.deposit,names.arg=c("Detrital\nmud","Shell\nhash","Oyster\nshell","Coarse\nsand","Seagrass"),space=c(.05,.05,.05,.05,.05),ylim=c(0,250),col=c("gray28","gray94","black","darkgrey","white"),border=c("black","black","black","black","black"),xlab="Habitat",ylab=txt1,cex.lab=1.4,cex.axis=1.2,cex.names=1.2,main=txt2,cex.main=1.2)
arrows(onewayplot2,means.deposit+ses.deposit,onewayplot2,means.deposit-ses.deposit,code=3,angle=90,length=0.05,lwd=1)
text(0.55, means.deposit[1]+ses.deposit[1]+15, "a",cex=1.2)
text(1.61, means.deposit[2]+ses.deposit[2]+15, "ab",cex=1.2)
text(2.655, means.deposit[3]+ses.deposit[3]+15, "ab",cex=1.2)
text(3.69, means.deposit[4]+ses.deposit[4]+15, "ab",cex=1.2)
text(4.755, means.deposit[5]+ses.deposit[5]+15, "b",cex=1.2)

txt1=expression(paste("Mean density DBSF \n     bivalves \u00b1 SE"))
txt2=expression(paste("b.                                                                                                                     "))
par(mar=c(4.0,6.0,3.0,1),mgp=c(3,1.5,0),cex=1)
means.DBSF=c(mean(DBSF[sb=="DT"]),mean(DBSF[sb=="SH"]),mean(DBSF[sb=="OS"]),mean(DBSF[sb=="GR"]),mean(DBSF[sb=="SG"]))
ses.DBSF=c(sd(DBSF[sb=="DT"])/sqrt(length(DBSF[sb=="DT"])),sd(DBSF[sb=="SH"])/sqrt(length(DBSF[sb=="SH"])),sd(DBSF[sb=="OS"])/sqrt(length(DBSF[sb=="OS"])),sd(DBSF[sb=="GR"])/sqrt(length(DBSF[sb=="GR"])),sd(DBSF[sb=="SG"])/sqrt(length(DBSF[sb=="SG"])))
onewayplot2=barplot(means.DBSF,names.arg=c("Detrital\nmud","Shell\nhash","Oyster\nshell","Coarse\nsand","Seagrass"),space=c(.05,.05,.05,.05,.05),ylim=c(0,40),col=c("gray28","gray94","black","darkgrey","white"),border=c("black","black","black","black","black"),xlab="Habitat",ylab=txt1,cex.lab=1.4,cex.axis=1.2,cex.names=1.2,main=txt2,cex.main=1.2)
arrows(onewayplot2,means.DBSF+ses.DBSF,onewayplot2,means.DBSF-ses.DBSF,code=3,angle=90,length=0.05,lwd=1)
lengths.DBSF=c(length(DBSF[sb=="DT"]),length(DBSF[sb=="SH"]),length(DBSF[sb=="GR"]),length(DBSF[sb=="OS"]),length(DBSF[sb=="SG"]))
# text(0.55, 26, "a",cex=1.2)
# text(1.61, 38, "ab",cex=1.2)
# text(2.655, 34, "ab",cex=1.2)
# text(3.69, 35, "ab",cex=1.2)
# text(4.755, 36, "b",cex=1.2)

means.TSSB=c(mean(TSSB[sb=="DT"]),mean(TSSB[sb=="SH"]),mean(TSSB[sb=="OS"]),mean(TSSB[sb=="GR"]),mean(TSSB[sb=="SG"]))
ses.TSSB=c(sd(TSSB[sb=="DT"])/sqrt(length(TSSB[sb=="DT"])),sd(TSSB[sb=="SH"])/sqrt(length(TSSB[sb=="SH"])),sd(TSSB[sb=="OS"])/sqrt(length(TSSB[sb=="OS"])),sd(TSSB[sb=="GR"])/sqrt(length(TSSB[sb=="GR"])),sd(TSSB[sb=="SG"])/sqrt(length(TSSB[sb=="SG"])))
txt1=expression(paste("Mean density TSSD \n     bivalves \u00b1 SE"))
txt2=expression(paste("c.                                                                                                                     "))
par(mar=c(4.0,6.0,3.0,1),mgp=c(3,1.5,0),cex=1)
onewayplot2=barplot(means.TSSB,names.arg=c("Detrital\nmud","Shell\nhash","Oyster\nshell","Coarse\nsand","Seagrass"),space=c(.05,.05,.05,.05,.05),ylim=c(0,80),col=c("gray28","gray94","black","darkgrey","white"),border=c("black","black","black","black","black"),xlab="Habitat",ylab=txt1,cex.lab=1.4,cex.axis=1.2,cex.names=1.2,main=txt2,cex.main=1.2)
arrows(onewayplot2,means.TSSB+ses.TSSB,onewayplot2,means.TSSB-ses.TSSB,code=3,angle=90,length=0.05,lwd=1)
text(0.55, means.TSSB[1]+ses.TSSB[1]+5, "a",cex=1.2)
text(1.61, means.TSSB[2]+ses.TSSB[2]+5, "a",cex=1.2)
text(2.655, means.TSSB[3]+ses.TSSB[3]+5, "ab",cex=1.2)
text(3.69, means.TSSB[4]+ses.TSSB[4]+5, "ab",cex=1.2)
text(4.755, means.TSSB[5]+ses.TSSB[5]+5, "b",cex=1.2)

means.hard=c(mean(hard[sb=="DT"]),mean(hard[sb=="SH"]),mean(hard[sb=="OS"]),mean(hard[sb=="GR"]),mean(hard[sb=="SG"]))
ses.hard=c(sd(hard[sb=="DT"])/sqrt(length(hard[sb=="DT"])),sd(hard[sb=="SH"])/sqrt(length(hard[sb=="SH"])),sd(hard[sb=="OS"])/sqrt(length(hard[sb=="OS"])),sd(hard[sb=="GR"])/sqrt(length(hard[sb=="GR"])),sd(hard[sb=="SG"])/sqrt(length(hard[sb=="SG"])))
txt1=expression(paste("Mean density ARM\n   bivalves \u00b1 SE"))
txt2=expression(paste("d.                                                                                                                     "))
par(mar=c(4.0,6.0,3.0,1),mgp=c(3,1.5,0),cex=1)
onewayplot2=barplot(means.hard,names.arg=c("Detrital\nmud","Shell\nhash","Oyster\nshell","Coarse\nsand","Seagrass"),space=c(.05,.05,.05,.05,.05),ylim=c(0,30),col=c("gray28","gray94","black","darkgrey","white"),border=c("black","black","black","black","black"),xlab="Habitat",ylab=txt1,cex.lab=1.4,cex.axis=1.2,cex.names=1.2,main=txt2,cex.main=1.2)
arrows(onewayplot2,means.hard+ses.hard,onewayplot2,means.hard-ses.hard,code=3,angle=90,length=0.05,lwd=1)
text(0.55, means.hard[1]+ses.hard[1]+1.5, "a",cex=1.2)
text(1.61, means.hard[2]+ses.hard[2]+1.5, "a",cex=1.2)
text(2.655, means.hard[3]+ses.hard[3]+1.5, "b",cex=1.2)
text(3.69, means.hard[4]+ses.hard[4]+1.5, "ab",cex=1.2)
text(4.755, means.hard[5]+ses.hard[5]+1.5, "ab",cex=1.2)

dev.off()



#####################################################
#Calculating % difference in diversity between seagrass and other habitats.
div.means=function(subs="SG")
{
dependent=mya[,which(colnames(mya)=="Diversity")]
var1=as.factor(mya[,which(colnames(mya)=="Year")])
var2=mya[,which(colnames(mya)=="Season")]
var3=mya[,which(colnames(mya)=="River")]
var5=mya[,which(colnames(mya)=="Sb")]
var6=mya[,which(colnames(mya)=="Volume")]
var7=mya[,which(colnames(mya)=="Fish")]
var8=mya[,which(colnames(mya)=="Crabs")]
var9=mya[,which(colnames(mya)=="Rays")]

lm.dep=lm(dependent~var1+var2+var3+var5+var6+var7+var8+var9)

newdat = expand.grid(var1=c("2012","2013"),var2=c("Fall","Spring","Summer"),
                     var3=c("Mobjack","York","Lynnhaven"),var5=c(subs),
                     var6=mean(na.omit(var6)),var7=mean(na.omit(var7)),
                     var8=mean(na.omit(var8)),var9=mean(na.omit(var9)))
pred.mat=data.frame()
ys=vector(length=1000)

for(i in 1:2){
  if (i==1) {
    year.string="2012"
  } else {
    year.string="2013"
  }
  for(j in 1:2){
    if(j==1) {
      season.string="Spring"
    } else {
      if(j==2) {
        season.string="Summer"
      } else {
        season.string="Fall"
      }
    }
    for(k in 1:3){
      if(k==1) {
        river.string="Lynnhaven"
      } else {
        if(k==2) {
          river.string="Mobjack"
        } else {
          river.string="York"
        }
      }
      
          sb.string=subs
        
        new.row.y=predict(lm.dep, newdata=newdat[newdat$var1==year.string&newdat$var2==season.string&newdat$var3==river.string&newdat$var5==sb.string,],interval="confidence",level=0.95)[,1]
         
        pred.mat=rbind(pred.mat,unname(new.row.y))

      }
    }
  }

ys=apply(pred.mat,2,mean)
ys
}

SG.mean=as.numeric(div.means())
# SG 0.5334274
OS.mean=as.numeric(div.means("OS"))
# OS 0.3176492
SH.mean=as.numeric(div.means("SH"))
# SH 0.2756183
GR.mean=as.numeric(div.means("GR"))
# GR 0.2853604
DT.mean=as.numeric(div.means("DT"))
# DT 0.3029859

# Seagrass increased bivalve diversity by 68%, 76%, 87%, and 94% when compared 
# to oyster shell, detrital mud, coarse sand, and shell hash habitats, 
# respectively, in models. 

SG.mean/OS.mean-1
# 0.6788923
SG.mean/SH.mean-1
# 0.935185
SG.mean/GR.mean-1
# 0.8692186
SG.mean/DT.mean-1
# 0.7603195


# Samples were collected in detrital mud (n = 21), shell hash (n = 61), 
# oyster shell (n = 30), coarse sand (n = 13), and seagrass (n = 55). 
# Letters denote significant differences at ?? = 0.05.

nrow(mya[mya$Sb=="DT",])
nrow(mya[mya$Sb=="SH",])
nrow(mya[mya$Sb=="OS",])
nrow(mya[mya$Sb=="GR",])
nrow(mya[mya$Sb=="SG",])
