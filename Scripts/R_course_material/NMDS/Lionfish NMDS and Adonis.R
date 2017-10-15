##Permutation test for differences in groups

## start with initial survey (IS)
IS <- read.table ("C:/Users/TJ/Documents/R/Data/PTVOEPST/NMSplot/ReefFish_Initial_RVert_Log.txt", header = TRUE)
## load final survey (FS)
FS <- read.table ("C:/Users/TJ/Documents/R/Data/PTVOEPST/NMSplot/ReefFish_Final_RVert_Log.txt", header = TRUE)
##Load environmental matrix
Env <- read.table ("C:/Users/TJ/Documents/R/Data/PTVOEPST/RFFinalEnv.txt", header = TRUE, colClasses=c("character", "numeric", "character", "character"))
attach(Env)
## keys, first column with the site/sample/community names cannot have a title
## I have already converted my data in excel and then run it in here, in this data I have removed some as well as
## log transformed the data too.
## I also find importing it as text makes things run smoother, not sure why.
## leave the community matrix alone, but attached the Env matrix

##get adonis from the vegan package
library(vegan)

## running adonis
## adonis is a perumational multivariate analysis of variance using distance matrices
Initial.1 <- adonis (IS ~ PoolTreat + ReefType + PTVOabund,
	data = IS,
	permutations = 999, 
	method = "bray")
Initial.1
Initial.2<- adonis (IS ~ PoolTreat + ReefType,
	data = IS,
	permutations = 999, 
	method = "bray")
Initial.2

Final.1 <- adonis (FS ~ PoolTreat + ReefType + PTVOabund,
	data = FS,
	permutations = 999, 
	method = "bray")
Final.1
Final.2<- adonis (FS ~ PoolTreat + ReefType,
	data = FS,
	permutations = 999, 
	method = "bray")
Final.2

##now to plot the NMS ordination
library(lattice)
modIS.n<-metaMDS(IS, Env$poolTreat, distance = "bray", k=3, try=10, trymax=50)
modIS.n

modIS <- MDSrotate(modIS, Env$PoolTreat) ## rotates the ordination based on some environmental factors

nf<-layout(matrix(c(1,2), nrow=1,ncol=2, byrow=TRUE))
ordiplot(modIS.n)
ordihull(modIS.n,groups=PoolTreat,draw="polygon",col="grey40",
         label=FALSE)
ordiplot(modIS)
ordihull(modIS,groups=PoolTreat,draw="polygon",col="grey80",
         label=FALSE)

modFS<-metaMDS(FS, distance = "bray", k=3, try=10, trymax=50)
modFS
MDSrotate(modFS, Env$PoolTreat)



windows(6.5,6.5) 
nf<-layout(matrix(c(1,2,3,4), nrow=2,ncol=2, byrow=TRUE))
layout.show(nf)

#############
############# initial reef type
par (lwd=2, mar=c(2,2,1,1), xpd=NA)
plot(modIS, type="n", xlim=c(-0.8,1.2), ylim=c(-0.5,.5),xaxt="n", yaxt="n",ylab="", xlab="")
points(modIS, display = "sites", cex=1.25, pch=c(15,15,15,2,15,2,2,2,2,15,15,15,15,2,2,2,2,15,15,15,15,15,2,15,15,15,2,2), 
                                          col= c("darkmagenta","darkmagenta","darkmagenta","darkorange1","darkmagenta","chocolate2","chocolate2","chocolate2","chocolate2","darkmagenta","darkmagenta","darkmagenta","darkmagenta",
                                                                                  "chocolate2","chocolate2","chocolate2","chocolate2","darkmagenta","darkmagenta","darkmagenta","darkmagenta","darkmagenta",
                                                                                                            "chocolate2","darkmagenta","darkmagenta","darkmagenta","chocolate2","chocolate2"))

ordihull(modIS, ReefType=="1", draw="lines", label=F, lwd=3, lty=3, col="chocolate2", show.groups=T)
ordihull(modIS, ReefType=="2", draw="lines", label=F, lwd=3, lty=2, col="darkmagenta", show.groups=T)

legend(.2,1, legend=c("artificial reef", "transplant reef"), pch=c(15,2), lty=c(2,3), cex=1, bty="n", col=c("darkmagenta","chocolate2"), pt.cex=1, lwd=2, merge=FALSE)


#############
############# final reef type
par (lwd=2, mar=c(2,1,1,2), xpd=NA)
plot(modFS, type="n", xlim=c(-1.3,.5), ylim=c(-.5,1),xaxt="n", yaxt="n",ylab="", xlab="")
points(modFS, display = "sites", cex=1.25, pch=c(15,15,15,2,15,2,2,2,2,15,15,15,15,2,2,2,2,15,15,15,15,15,2,15,15,15,2,2),
                                            col= c("darkmagenta","darkmagenta","darkmagenta","darkorange1","darkmagenta","chocolate2","chocolate2","chocolate2","chocolate2","darkmagenta","darkmagenta","darkmagenta","darkmagenta",
                                                                                  "chocolate2","chocolate2","chocolate2","chocolate2","darkmagenta","darkmagenta","darkmagenta","darkmagenta","darkmagenta",
                                                                                                            "chocolate2","darkmagenta","darkmagenta","darkmagenta","chocolate2","chocolate2"))

ordihull(modFS, ReefType=="1", scaling = 3, draw="lines", label=F, lwd=3, lty=3, col="chocolate2", show.groups=T)
ordihull(modFS, ReefType=="2", scaling = 3, draw="lines", label=F, lwd=3, lty=2, col="darkmagenta", show.groups=T)

#############
############# initial pool treat
par (lwd=2, mar=c(2,2,1,1), xpd=NA)
plot(modIS, type="n", xlim=c(-0.7,1.3), ylim=c(-0.7,1),xaxt="n", yaxt="n",ylab="", xlab="")
points(modIS, display = "sites", cex=1.25, pch=c(18,18,18,18,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,16,16,16,16,16,16), 
                                          col= c("red","red","red","red",
                                                              "blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue",
                                                                                        "green","green","green","green","green","green",
                                                                                                    "black","black","black","black","black","black"))

ordihull(modIS, PoolTreat=="0", scaling = 3, draw="lines", label=F, lwd=4, lty=2, col="red", show.groups=T)
ordihull(modIS, PoolTreat=="1", scaling = 3, draw="lines", label=F, lwd=4, lty=9, col="blue", show.groups=T)
ordihull(modIS, PoolTreat=="2", scaling = 3, draw="lines", label=F, lwd=4, lty=6, col="green", show.groups=T)
ordihull(modIS, PoolTreat=="3", scaling = 3, draw="lines", label=F, lwd=4, lty=1, col="black", show.groups=T)

legend(-.2,1.3, legend=c("lionfish only", "lionfish & low grouper", "lionfish & high grouper", "control"), pch=c(18,0,2,16), lty=c(2,9,6,1),  
      cex=1, bty="n", col=c("red","blue","green","black"), pt.cex=2, lwd=3, merge=FALSE)
text (-.8, .25, labels=c("NMDS 2"), cex=2.15, srt=90)

#############
############# final pool treat
par (lwd=2, mar=c(2,3,2,2), xpd=NA)
plot(modFS, type="n", ylab="", xlab="", xlim=c(-1.3,.5), ylim=c(-.5,1), xaxt="n", yaxt="n")

points(modFS, display = "sites", cex=1.25, pch=c(18,18,18,18,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,16,16,16,16,16,16),
                                            col= c("white","white","white","white",
                                                              "white","white","white","white","white","white","white","white","white","white","white","white","white",
                                                                                        "white","white","white","white","white","white",
                                                                                                    "black","black","black","black","black","black"))
ordihull(modFS, PoolTreat=="3", scaling = 3, draw="lines", label=F, lwd=4, lty=1, col="black", show.groups=T)


points(modFS, display = "sites", cex=1.25, pch=c(18,18,18,18,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,16,16,16,16,16,16),
                                            col= c("red","red","red","red",
                                                              "white","white","white","white","white","white","white","white","white","white","white","white","white",
                                                                                        "white","white","white","white","white","white",
                                                                                                    "black","black","black","black","black","black"))
ordihull(modFS, PoolTreat=="0", scaling = 3, draw="lines", label=F, lwd=4, lty=2, col="red", show.groups=T)

points(modFS, display = "sites", cex=1.25, pch=c(18,18,18,18,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,16,16,16,16,16,16),
                                            col= c("red","red","red","red",
                                                              "blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue",
                                                                                        "white","white","white","white","white","white",
                                                                                                    "black","black","black","black","black","black"))
ordihull(modFS, PoolTreat=="1", scaling = 3, draw="lines", label=F, lwd=4, lty=9, col="blue", show.groups=T)

points(modFS, display = "sites", cex=1.25, pch=c(18,18,18,18,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,16,16,16,16,16,16),
                                            col= c("red","red","red","red",
                                                              "blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue",
                                                                                        "green","green","green","green","green","green",
                                                                                                    "black","black","black","black","black","black"))
ordihull(modFS, PoolTreat=="2", scaling = 3, draw="lines", label=F, lwd=4, lty=6, col="green", show.groups=T)


text(-.5, -.8, labels=c("NMDS 1"), cex=2.15)




### looks at the correlation of each factor with the nmds axis 

ord.fit.IS <-envfit(modIS ~ ReefType + EPSTTreat, data=Env, perm =1000)
ord.fit.IS
ord.fit.FS <-envfit(modFS ~ ReefType + EPSTTreat, data=Env, perm =1000)
ord.fit.FS
