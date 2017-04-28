#######################################################################################################
#                                                                                                     #
#                   ANALYSIS OF TAXONOMIC, FUNCTIONAL, AND PHYLOGENETIC DIVERSITY                     #
#                              IN CHESAPEAKE BAY DEMERSAL FISHES                                      #
#                                                                                                     #
#######################################################################################################

#Author: Jon Lefcheck, jslefche@vims.edu
#Please contact at address above if any of the below code is to be used in any other published analyses
#Last updated: 2013-10-31

#######################################################################################################
#                                        TABLE OF CONTENTS                                            #
#   Line  26: Required libraries                                                                      #
#   Line  56: Importing and formatting the data                                                       #
#   Line 125: Taxonomic distinctness (taxonomic tree)                                                 #
#   Line 147: Functional dendrogram                                                                   #
#   Line 220: Phylogenetic tree                                                                       #
#   Line 286: Calculating local diversity                                                             #
#   Line 376: Characterizing bivariate relationships (SPLOMs, Mantel tests)                           #
#   Line 398: Mapping local diversity                                                                 #
#   Line 408: Fitting GAMs and examining spatial signal in residuals                                  #
#   Line 630: Fitting GAMs and partitioning deviance                                                  #
#                                                                                                     #
#######################################################################################################

library(ade4) #Calls: mantel.rtest, dudi.mix
library(ape) #Calls: read.tree, read.nexus, chronopl, as.phylo
library(automap) #Calls: autoKrige
library(car) #Calls: recode
library(clue) #Calls: cl_consensus, cl_ultrametric
library(FD) #Calls: gowdis
library(fpc) #Calls: pamk
library(ggplot2) #Calls: ggplot
library(gridExtra) #Calls: grid.arrange
library(mgcv) #Calls: gam
library(MuMIn) #Calls: dredge
library(parallel) #Calls: clusterEvalQ,clusterExport,detectCores,makeCluster,parLapply
library(picante) #Calls: phylosignal
library(plotrix) #Calls: std.error
library(plyr) #Calls: ddply
library(psych) #Calls: pairs.panels
library(reshape2) #Calls: dcast, melt
library(rgdal) #Calls: proj4string, readOGR, spTransform
library(sp) #Calls: coordinates, spsample
library(stringr) #Calls: str_split_fixed
library(tree) #Calls: tree
library(vegan) #Calls: metaMDS, ordiplot, taxa2dist

#Create socket cluster of multiple cores 
cl=makeCluster(detectCores())
#Load automap on all cores (to access function `autoKrige`)
clusterEvalQ(cl,c(library(automap),library(mgcv),library(plotrix)) )

#######################################################################################################
#                                IMPORTING AND FORMATTING THE DATA                                    #
#######################################################################################################

#Import catch (abundance) from file: ChesMMAP_CatchData_3-20-12
catch_abund=read.csv("chesmmap_catch_abund.csv")
#Remove column "s9999" which is just a placeholder
catch_abund=catch_abund[,-which(colnames(catch_abund)=="s9999")]
#Import catch (biomass) from file: ChesMMAP_CatchData_3-20-12
catch_biomass=read.csv("chesmmap_catch_biomass.csv")
catch_biomass=catch_biomass[,-which(colnames(catch_biomass)=="s9999")]
#Import species list from file: ChesMMAP Species for Diversity class_3-20-12
species_list=read.csv("chesmmap_species_list.csv")
#Import functional trait matrix from file: Fish Functional Traits MASTER V6.3
traits=read.csv("chesmmap_species_traits.csv")

#Replace species codes with species names in the catch datasets
#First, recode species codes so they match columns in catch_abund/catch_biomass
species_list$newcode=sprintf("s%04d",species_list$VIMSCODE)
#Extract vector of species codes from the column "variable"
sp_codes=colnames(catch_abund)[18:ncol(catch_abund)]
#Pull out of a vector of names appearing in sp_codes from the species_list
sp_names=species_list[match(sp_codes,species_list$newcode),3]
#Apply the vector of names as the levels of the factor "variable" containing the codes
colnames(catch_abund)[18:ncol(catch_abund)]=as.vector(sp_names)
#Repeat for biomass data
colnames(catch_biomass)[18:ncol(catch_biomass)]=as.vector(sp_names)

#Set species names as row.names and remove extra columns
rownames(traits)=traits$LATIN_NAME; traits=traits[,-c(1:3)]
#Arrange trait matrix alphabetically by species
traits=traits[order(rownames(traits)),]
#Arrange species columns of catch_abund and catch_biomass alphabetically
catch_abund=cbind(catch_abund[,1:17],
 catch_abund[,18:ncol(catch_abund)][,order(names(catch_abund[,18:ncol(catch_abund)]))])
catch_biomass=cbind(catch_biomass[,1:17],
 catch_biomass[,18:ncol(catch_biomass)][,order(names(catch_biomass[,18:ncol(catch_biomass)]))])

#Convert months from numbers to names
catch_abund=transform(catch_abund,Month=recode(Month,"'3'='Mar';'5'='May';'7'='July';'9'='Sept';'11'='Nov'"))
catch_abund$Month=factor(catch_abund$Month,levels=c("Mar","May","July","Sept","Nov"))
catch_biomass=transform(catch_biomass,Month=recode(Month,"'3'='Mar';'5'='May';'7'='July';'9'='Sept';'11'='Nov'"))
catch_biomass$Month=factor(catch_biomass$Month,levels=c("Mar","May","July","Sept","Nov"))

#Calculate species richness at each station
catch_abund$richness=rowSums(catch_abund[,18:ncol(catch_abund)]>0)
catch_biomass$richness=rowSums(catch_biomass[,18:ncol(catch_biomass)]>0)

#Remove all rows where S = 0 or 1 (does not inform about biodiversity or evenness)
catch_abund=catch_abund[catch_abund$richness>=2,]
catch_biomass=catch_biomass[catch_biomass$richness>=2,]

#Convert observations into presence/absence
catch_pres.abs=catch_abund
catch_pres.abs[,18:67][catch_pres.abs[,18:67]>0]=1

#Remove SOURCE columns from the trait matrix
remove=grep("SOURCE",names(traits),value=F)
traits=traits[,-remove]
#Change misclassified categorical traits to factors
traits[,c(16,20)]=lapply(traits[,c(16,20)],function(x) as.factor(x))
#Change some traits to ordinal
traits[,10]=ordered(traits[,10],levels=c("short","medium","long","verylong"))
traits[,16]=ordered(traits[,16],levels=c("3","5","7","9","11"))
traits[,17]=ordered(traits[,17],levels=c("2","3","4","5"))
traits[,18]=ordered(traits[,18],levels=c("River","Estuary","Estuary, Ocean","Ocean"))
traits[,19]=ordered(traits[,19],levels=c("Spring","Summer","Spring, Fall","Fall","Winter","Winter, Summer"))
#str(traits)

#######################################################################################################
#                                        TAXONOMIC TREE                                               #
#######################################################################################################

#Read in classification table
taxo=read.csv("taxonomic_table.csv")
taxo=taxo[order(taxo[,1]),]; rownames(taxo)=taxo[,1]
#Generate cladogram
taxo.tree=as.phylo(~Class/Subclass/Infraclass/Superorder/Order/Family/Genus/DBDGS.genus.species,data=taxo)
#Plot cladogram
#plot(taxo.tree)

#Create ultrametric species x species distance matrix (from Clarke and Warwick 1998 J Appl Ecol)
taxo.dist=cl_ultrametric(hclust(taxa2dist(taxo,check=TRUE,varstep=TRUE)/100))
#Plot tree from distance matrix
# par(mar=c(3,1,0,16))
# plot(taxo.dist,horiz=TRUE)
#Save plot: 10" x 15"

#Write newick tree
#write.tree(taxo.tree,"Taxonomic tree.new")

#######################################################################################################
#                                     FUNCTIONAL DENDROGRAM                                           #
#######################################################################################################

#Remove traits that are highly correlated (>0.8)
traits[,c("MEAN_WEIGHT","MEAN_LENGTH","MAX_LENGTH","DORSAL_SPINE_LENGTH")]=list(NULL)

#Downweight diet information (since they are not independent) 
gower_weights=rep(1,ncol(traits))
gower_weights[grep("DIET",names(traits),value=F)]=1/length(grep("DIET",names(traits),value=F))
#Calculate Gower distances
traits.dist=gowdis(traits,w=gower_weights,ord="podani")

#Build each tree for the different clusting methods (Mouchet et al 2008 Oikos)
#Methods "ward","median", and "centroid" are not recommended for nonmetric distances (such as Gower)
tree_methods=c("single","complete","average","mcquitty","ward") #,"median","centroid")  
trees=lapply(tree_methods,function(i) hclust(traits.dist, method=i))
par(mfrow=c(3,2))
for(i in 1:length(trees)) {plot(trees[[i]])}

#Convert trees to ultrametric 
trees.ultra=lapply(trees,function(i) cl_ultrametric(as.hclust(i)))
#Plot each tree
par(mfrow=c(3,2))
for (i in 1:length(trees.ultra)) {plot(trees.ultra[[i]])}

#Build the consensus tree (Mouchet et al 2008 Oikos) use Clue package 
ensemble.trees=cl_ensemble(list=trees)
class(ensemble.trees)
consensus.tree=cl_consensus(ensemble.trees) 
#plot(consensus.tree)

#Calculate dissimilarity values for each tree (using 2-norm, Merigot et al 2010 Ecology)
all.trees=c(trees.ultra,consensus.tree[1])
names(all.trees)=c(tree_methods,"consensus")
(trees.dissim=lapply(all.trees,function(i) cl_dissimilarity(i,traits.dist,method="spectral")))
#Identify best tree and isolate
trees.dissim2=do.call(rbind,trees.dissim)
min.tree=which.min(trees.dissim2)
names(all.trees)[min.tree]
func.dist=all.trees[names(all.trees)==names(all.trees)[min.tree]][[1]]
#Confirm lowest 2-norm value, spectral norm (2-norm) of the differences if the ultrametrics 
cl_dissimilarity(func.dist,traits.dist,method="spectral")
#Scale between 0-1
func.dist=func.dist/max(func.dist)
#Plot the best tree
# par(mfrow=c(1,1))
# par(mar=c(3,1,0,16))
# plot(func.dist,horiz=TRUE)
#Save plot: 10" x 15"

#Write newick tree
write.tree(as.phylo(as.hclust(func.dist)),"Functional dendrogram.new")

#Visualize species' differences in multivariate trait space
#Perform k-means clustering with no a priori specification for k
traits.kclus=pamk(traits.dist,krange=2:10)
# #Perform multidimensional scaling on functional dendrogram
traits_nmds=metaMDS(traits.dist,k=traits.kclus$nc,trymax=500)
# #Plot in two dimensions
par(mar=c(4,4,1,1))
ordiplot(traits_nmds,type="n")
#Assign colors to different groups
groups=levels(factor(traits.kclus$pamobject$clustering))
points.symbols=15:16
points.colors=c("firebrick3","cornflowerblue")
for(i in seq_along(groups)) {
   points(traits_nmds$points[traits.kclus$pamobject$clustering==groups[i],],
          pch=points.symbols[i],col=points.colors[i],cex=1.4) }
 ordispider(traits_nmds,factor(traits_kclus$pamobject$clustering),label=F)
 ordihull(traits_nmds,factor(traits.kclus$pamobject$clustering),lty="dotted")
 orditorp(traits_nmds,dis="sites",pcex=0,air=0.5,col="grey10",cex=0.8)

#######################################################################################################
#                                       PHYLOGENETIC TREE                                             #
#######################################################################################################

#Read in the Newick tree from PHYLIP
phylo.tree=read.tree("1105_RAxML_bipartitionsBS_BestTree_FULLNAME.phylip"); phylo.tree
#par(mar=c(0,0,0,0))
#Check to see if rooted
is.rooted(phylo.tree)
phylo.tree=root(phylo.tree,"outgrpLamp",resolve.root=T)
#Plot tree
#plot(phylo.tree,show.node.label=T); add.scale.bar(2,0,length=0.3,lwd=3)
#Check to see if ultrametric
is.ultrametric(phylo.tree) #False

#Estimate divergence times using penalized-likelihood (PL)

#Determine lambda (smoothing parameter) through cross-validation
#Modified from: http://www.r-phylo.org/wiki/HowTo/Divergence_Time_Estimation
# l=10^(-2:10); cv=numeric(length(l))
# for(i in 1:length(l))
#   cv[i]=sum(attr(chronopl(phylo.tree,lambda=l[i],CV=TRUE),"D2")) #S = the number of substitution sites
# l[which.min(cv)] #lambda = 10

#Convert the tree to ultrametric using the smoothing parameter lambda that results in the lowest CV
phylo.ultratree=chronopl(phylo.tree,lambda=100) 
#Root ultrametric tree
is.ultrametric(phylo.ultratree) #TRUE
is.rooted(phylo.ultratree) #FALSE -> will root when outgroups are dropped
is.binary.tree(phylo.ultratree) #TRUE

#Prune outgroups
phylo.ultratree=drop.tip(phylo.ultratree,grep("outgrp",phylo.ultratree$tip.label,value=T)); phylo.ultratree
#Replace species names with those from the trait matrix
phylo.ultratree$tip.label=rownames(traits)[match(phylo.ultratree$tip.label,make.names(rownames(traits)))]
#Plot tree
# par(mar=c(4,4,1,1))
# plot(rotate(phylo.ultratree,53),show.node.label=F); axisPhylo(side=1)

#Write newick tree
#write.tree(phylo.ultratree,"Phylogenetic tree.new")

#Extract cophenetic distance matrix and alphabetize labels
phylo.dist=as.matrix(cophenetic(phylo.ultratree))
ordering=sort(rownames(phylo.dist))
phylo.dist=as.dist(as.matrix(phylo.dist)[ordering,ordering])
#Scale between 0-1
phylo.dist=phylo.dist/max(phylo.dist) 

# #Examine phylogenetic signal of traits (Blomberg's K)
# #Can only use continuous values, so subset only continuous traits
# traits_cont=traits[,sapply(traits, is.numeric)]
# traits_cont[is.na(traits_cont)]=0
# #Loop over each continuous trait value and calculate K and p-value (if different from 0)
# traits_K=data.frame()
# for(i in 1:ncol(traits_cont)) {
#   trait=traits_cont[,i]
#   names(trait)=rownames(traits_cont)
#   K=phylosignal(trait[phylo.ultratree$tip.label],phylo.ultratree)
#   traits_K[i,1]=K$K; traits_K[i,2]=K$PIC.variance.P }
# #Value close to 0 indicate no phylogenetic signal, values close to 1 indicate phylogenetic signal
# names(traits_K)=c("K","p")
# rownames(traits_K)=colnames(traits[,sapply(traits, is.numeric)])
# traits_K[order(traits_K$K),]
# sum(traits_K$p<0.05)/nrow(traits_K)

#######################################################################################################
#                                   CALCULATING LOCAL DIVERSITY                                       #
#######################################################################################################

#Calculate alpha diversity 
alphadiv.list=lapply(list(catch_abund,catch_biomass,catch_pres.abs),
 function(i) {
   #Extract community matrix from object in list
   mat=i[,c(1,18:(ncol(i)-1))] #MH: calling all rows and column 1 and >18
   rownames(mat)=mat[,1]; mat=mat[,-1]
   #Calculate relative values for community matrix
   rel.mat=mat/apply(mat,1,sum) #MH: counts divided by sum of rows to give relative abundance values  
   #Compute species diversity
   species.dist=matrix(1,ncol(mat),ncol(mat))-diag(rep(1,ncol(mat))) #MH: creates a matrix with number of columns of matrix for rows and columns  
   species.div=1/(1-apply(rel.mat,1, function(x) t(x) %*% species.dist %*% x)) #MH apply to all rows in rel.mat the function x, transpose, multiply matrix
   #Compute evenness (as in Jost 2010 Diversity)
   evenness=log(species.div)/log(rowSums(mat>0))
   #Ensure that all evenness calculations are not >1 (bug)
   evenness[evenness>1]=1
   #Compute functional diversity 
   func.div=1/(1-apply(rel.mat,1, function(x) t(x) %*% as.matrix(func.dist) %*% x))
   #Compute phylogenetic diversity
   phylo.div=1/(1-apply(rel.mat,1, function(x) t(x) %*% as.matrix(phylo.dist) %*% x))
   #Compute taxonomic diversity
   taxo.div=1/(1-apply(rel.mat,1,function(x) t(x) %*% as.matrix(taxo.dist) %*% x))
   #Bind all to the original dataframe
   cbind(i,evenness,species.div,func.div,phylo.div,taxo.div) } ) 
names(alphadiv.list)=c("abund","biomass","pres.abs")

#Generate figure legend based on CB map
#Read in shapefile and fortify
cb.shp=readOGR("./Chesapeake Bay Shapefiles/ChesMMAP buffer",layer="CMapbuffer25k")
cb.shp.fortify=fortify(cb.shp)
#Create data frame for legend points and titles
data=data.frame(labels=c("Lower","Lower-\nMiddle","Middle","Upper-\nMiddle","Upper"),
                shape=c(15:18,6),color=c("grey15","grey30","grey45","grey60","black"),
                long=c(-76.15,-76.135,-76.21,-76.41,-76.34),lat=c(37.2,37.7,38.15,38.67,39.16),
                lat2=c(37.08,37.54,38.03,38.51,39.04))
data$labels=factor(data$labels,levels=c("Lower","Lower-\nMiddle","Middle","Upper-\nMiddle","Upper"))
#Generate map legend for inset
inset.map=
  ggplot()+
  #Plot outline of Chesapeake Bay
  geom_polygon(data=subset(cb.shp.fortify,group==0.1),aes(x=long,y=lat,group=group),col="grey50",fill="white",lwd=0.4)+
  geom_segment(aes(x=-75.98281,y=37.39372,xend=-76.24769,yend=37.39372),col="grey50",lwd=1)+
  geom_segment(aes(x=-76.24572,y=37.87455,xend=-76.01332,yend=37.87455),col="grey50",lwd=1)+
  geom_segment(aes(x=-76.41753,y=38.35538,xend=-76.26395,yend=38.35538),col="grey50",lwd=1)+
  geom_segment(aes(x=-76.49554,y=38.83621,xend=-76.35090,yend=38.83621),col="grey50",lwd=1)+
  #Add legend points
  geom_point(data=data,aes(x=long,y=lat,shape=labels,color=labels),size=8)+
  scale_color_manual(values=c("black","grey60","grey45","grey30","grey15"))+
  scale_shape_manual(values=c(6,18:15))+
  #Add legend labels
  geom_text(data=data,aes(x=long,y=lat2,label=labels),size=4.6,fontface="bold")+
  theme(
    plot.margin=unit(c(0,0,0,0),"cm"),
    panel.background=element_blank(),
    legend.position="none",panel.border=element_blank(),
    panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
    axis.title.x=element_blank(),axis.title.y=element_blank(),
    axis.text.x=element_blank(),axis.text.y=element_blank(),
    axis.line=element_blank(),axis.ticks=element_blank())    

#Plot mean alpha diversity by region and month for each diversity index
lapply(alphadiv.list,function(i) {
  #Melt dataframe so diversity indices are in a single column
  x=melt(i,id.vars=c("Month","Year","Reg"),measure.vars=c("richness","evenness","species.div","func.div","phylo.div","taxo.div"))
  #Summarize means and SEs by month, region, and variable
  x=ddply(x,c("Month","Reg","variable"),summarize,value.mean=mean(value),value.se=std.error(value))
  #Rename levels for figure plotting
  levels(x$variable)=c("A) Richness","B) Evenness","C) Gini-Simpson","D) Functional","E) Phylogenetic","F) Taxonomic")
  #Rename regions for figure plotting
  x$Reg=factor(x$Reg); levels(x$Reg)=c("Upper","Upper-Middle","Middle","Lower-Middle","Lower")
  #Plot line graph
  p=ggplot(x,aes(x=as.factor(Month),y=value.mean,col=as.factor(Reg),shape=as.factor(Reg)))+
    #annotation_custom(grob=inset.map,xmin=19,xmax=22.5,ymin=0,ymax=1.5)+
    geom_line(aes(group=as.factor(Reg)),lwd=1)+geom_point(size=5)+
    scale_color_manual(values=c("grey15","grey30","grey45","grey60","black"),name="Region")+
    scale_shape_manual(values=c(15:18,6),name="Region")+
    geom_errorbar(aes(ymax=value.mean+value.se,ymin=value.mean-value.se),width=0.1)+facet_wrap(~variable,scales="free")+
    labs(x="",y="Diversity")+theme_bw(base_size=16)+
    theme(
      plot.margin=unit(c(0,1,0,0),"cm"),legend.position="none",
      panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.margin=unit(0.8,"lines"),
      strip.background=element_blank(),strip.text.x=element_text(size=18,hjust=0))
  #Add legend map
  grid.arrange(p,inset.map,nrow=1,widths=c(1,0.2))
} )
#Save PDF (13" x 8")
  
#######################################################################################################
#                               CHARACTERIZING BIVARIATE RELATIONSHIPS                                #
#######################################################################################################

#Create scatterplot matrix of different diversity measures against one another (9" x 9")
lapply(seq_along(alphadiv.list),function(i) {
  pairs.panels(alphadiv.list[[i]][,68:73],#main=names(alphadiv.list)[i],
               labels=c("Richness","Evenness","Gini-Simpson","Functional","Phylogenetic","Taxonomic"),
               hist.col="grey80",rug=T,smooth=F,ellipses=F,lm=F,method="spearman") } )
   
#Test significance of Spearman rank correlations (H0 that rho = 0)
do.call(rbind,lapply(list("richness","evenness","species.div","func.div","phylo.div","taxo.div"),function(j) {
  do.call(rbind,lapply(rev(list("richness","evenness","species.div","func.div","phylo.div","taxo.div")),function(k) {
    data.frame(title=paste(j,k,sep="~"),
               correlation=cor.test(alphadiv.list[[2]][,j],alphadiv.list[[2]][,k])$estimate,
               p.value=cor.test(alphadiv.list[[2]][,j],alphadiv.list[[2]][,k])$p.value) } ) ) } ) ) 

#Examine correlations between matrices using Mantel's test
mantel.rtest(taxo.dist,func.dist) #R^2 = 0.74
mantel.rtest(taxo.dist,phylo.dist) #R^2 = 0.90
mantel.rtest(func.dist,phylo.dist) #R^2 = 0.78

#######################################################################################################
#                                    MAPPING LOCAL DIVERSITY                                          #
#######################################################################################################

#Create new list of data frames to convert to SpatialPointsDataFrame
alphasp.list=alphadiv.list
names(alphasp.list)=names(alphadiv.list)
#Convert class("dataframe") into class("SpatialPointsDataFrame")
for(i in seq_along(alphasp.list)) { coordinates(alphasp.list[[i]])=~Long+Lat } 
lapply(alphasp.list,class)
#Project latlong coordinates onto an ellipse
for(i in seq_along(alphasp.list)) { 
  proj4string(alphasp.list[[i]])="+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs" }
#Transform the projection into cartesian coordinates
alphasp.list=lapply(seq_along(alphasp.list),function(i) {
  spTransform(alphasp.list[[i]],CRS("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")) } )
#Create a spatial grid onto which interpolated values can be plotted
#Read in 2.5 km buffered polygon based on position of ChesMMAP stations from ESRI shapefile
cb.shp=readOGR("./Chesapeake Bay Shapefiles/ChesMMAP buffer",layer="CMapbuffer25k")
#Generate interpolation grid
cb.shp=spTransform(cb.shp,CRS("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0 ")) 
grd=spsample(cb.shp,type="regular",n=2500,nsig=2) 
#par(mfrow=c(1,1))
#plot(grd)
#Send interpolation grid to each core
clusterExport(cl,varlist=c("grd"))

#Interpolate across all months and years by passing each matrix to a different core  
# #Perform kriging interpolation for each matrix (abundance-,biomass-, and presence/absence-weighted)
# krigeall.list=parLapply(cl,alphasp.list,function(i) {
#   lapply(names(i@data)[66:71],function(j) {
#     as.data.frame(autoKrige(formula(paste("log10(",j,")~1")),i,grd)$krige_output) } ) } )
# 
# #Apply names for titling during plotting
# names(krigeall.list)=c("abund","biomass","pres.abs")
# for(i in seq_along(krigeall.list)) {
#   names(krigeall.list[[i]])=c("richness","species","taxonomic","functional","phylogenetic","evenness") }  
# 
# #Plot all components across all years and months
# lapply(seq_along(krigeall.list),function(i) {
#   do.call(grid.arrange,lapply(seq_along(krigeall.list[[i]]),function(j) {
#     ggplot()+
#       geom_raster(data=krigeall.list[[i]][[j]],aes(x=x1,y=x2,fill=var1.pred))+
#       scale_fill_gradientn(colours=rev(rainbow(3)))+
#       coord_equal()+labs(x="",y="",title=paste(names(krigeall.list[[i]])[j]))+
#       theme_bw()+theme(plot.margin=unit(c(0,0,0,0),"cm"),
#                        axis.text.x=element_blank(),axis.text.y=element_blank(),
#                        axis.ticks=element_blank(),legend.title=element_blank()) } ) ) } )

#Now repeat but instead interpolate by month across years
#Perform kriging interpolation for each matrix (abundance-,biomass-, and presence/absence-weighted)
krigemonth.list=parLapply(cl,alphasp.list,function(i) {
  lapply(names(i@data)[66:71],function(j) {
    lapply(unique((i@data)$Month),function(k) {
      if(j=="evenness") {
        as.data.frame(autoKrige(formula(paste("asin(sqrt(",j,"))~1")),subset(i,i@data$Month==k),grd)$krige_output)
      } else {
        as.data.frame(autoKrige(formula(paste("log(",j,")~1")),subset(i,i@data$Month==k),grd)$krige_output)  }
} ) } ) } )

#Apply names for titling during plotting
names(krigemonth.list)=c("abund","biomass","pres.abs")
for(i in seq_along(krigemonth.list)) {
  names(krigemonth.list[[i]])=c("A) Richness","B) Evenness","C) Gini-Simpson","D) Functional","E) Phylogenetic","F) Taxonomic") }
for(i in seq_along(krigemonth.list)) {
  for(j in seq_along(krigemonth.list[[i]])) {
    names(krigemonth.list[[i]][[j]])=c("Mar","May","July","Sept","Nov") } }
#Collapse lists into data frame
krigemonth.list=lapply(krigemonth.list,function(i) lapply(i,function(j) do.call(rbind,j)) ) 
#Split rownames into month value and bind to dataframe
krigemonth.list=lapply(krigemonth.list,function(i) {
  lapply(i,function(j) {
    cbind(data.frame(str_split_fixed(rownames(j),"\\.",2)[,1]),j) } ) } )
for(i in seq_along(krigemonth.list)) {
  for(j in seq_along(krigemonth.list[[i]])) { names(krigemonth.list[[i]][[j]])[1]=c("month") } }
#Order levels of month for facetting
for(i in seq_along(krigemonth.list)) {
  for(j in seq_along(krigemonth.list[[i]])) { 
    krigemonth.list[[i]][[j]]$month=factor(krigemonth.list[[i]][[j]]$month,levels=c("Mar","May","July","Sept","Nov")) } }

#Plot and facet by month (12" x 18")
lapply(seq_along(krigemonth.list),function(i) {
  do.call(grid.arrange,c(lapply(seq_along(krigemonth.list[[i]]),function(j) {
    if(j!=2) {
      ggplot()+
        geom_raster(data=krigemonth.list[[i]][[j]],aes(x=x1,y=x2,fill=exp(var1.pred)))+
        scale_fill_gradientn(colours=rev(rainbow(3)))+
        facet_wrap(~month,nrow=1)+
        coord_equal()+labs(x="",y="",title=names(krigemonth.list[[i]])[j])+
        theme_bw()+theme(plot.margin=unit(c(0,0,0,0),"cm"),
                   plot.title=element_text(size=18,hjust=0,vjust=1),
                   strip.text=element_text(size=12),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                   axis.text.x=element_blank(),axis.text.y=element_blank(),
                   axis.ticks=element_blank(),legend.title=element_blank()) 
    } else {
      ggplot()+
        geom_raster(data=krigemonth.list[[i]][[j]],aes(x=x1,y=x2,fill=sin(var1.pred)^2))+
        scale_fill_gradientn(colours=rev(rainbow(3)))+
        facet_wrap(~month,nrow=1)+
        coord_equal()+labs(x="",y="",title=names(krigemonth.list[[i]])[j])+
        theme_bw()+theme(plot.margin=unit(c(0,0,0,0),"cm"),
                     plot.title=element_text(size=18,hjust=0,vjust=1),
                     strip.text=element_text(size=12),
                     panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                     axis.text.x=element_blank(),axis.text.y=element_blank(),
                     axis.ticks=element_blank(),legend.title=element_blank()) }
  } ), list(ncol=2)) )
} )

#######################################################################################################
#                                 EXAMINING SPATIAL SIGNAL IN RESIDUALS                               #
#######################################################################################################

#Extract residuals for each metric vs richness, evenness, and species diversity after fitting GAM
resids.list=lapply(alphadiv.list,function(i) { cbind(i,do.call(cbind,
  lapply(list(
    "log(species.div)~log(richness)","log(taxo.div)~log(richness)","log(func.div)~log(richness)","log(phylo.div)~log(richness)",
    "log(species.div)~log(taxo.div)","log(species.div)~log(func.div)","log(species.div)~log(phylo.div)",
    "log(func.div)~log(species.div)","log(func.div)~log(taxo.div)","log(func.div)~log(phylo.div)",
    "log(phylo.div)~log(species.div)","log(phylo.div)~log(taxo.div)","log(phylo.div)~log(func.div)",
    "log(taxo.div)~log(species.div)","log(taxo.div)~log(func.div)","log(taxo.div)~log(phylo.div)"),
    function(j) { res=as.data.frame(residuals(gam(formula(paste(j)),data=i))); return(res) } ) ) ) } )

#Replace column names for residuals
for(i in seq_along(resids.list)) { colnames(resids.list[[i]])[74:89]=c(
  "species.div.richness","taxo.div.richness","func.div.richness","phylo.div.richness",
  "species.div.taxo.div","species.div.func.div","species.div.phylo.div",
  "func.div.species.div","func.div.taxo.div","func.div.phylo.div",
  "phylo.div.species.div","phylo.div.taxo.div","phylo.div.func.div",
  "taxo.div.species.div","taxo.div.func.div","taxo.div.phylo.div") }
                   
#Save residuals as another object for mapping
residssp.list=resids.list
#Convert class("dataframe") into class("SpatialPointsDataFrame")
for(i in seq_along(residssp.list)) { coordinates(residssp.list[[i]])=~Long+Lat } 
lapply(residssp.list,class)
#Project latlong coordinates onto an ellipse
for(i in seq_along(residssp.list)) { 
  proj4string(residssp.list[[i]])="+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs" }
#Transform the projection into cartesian coordinates
residssp.list=lapply(seq_along(residssp.list),function(i) {
  spTransform(residssp.list[[i]],CRS("+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")) } )
# 
# #Interpolate across all months and years by passing each matrix to a different core 
# #Perform kriging interpolation for each matrix in the list (abundance- and biomass-weighted))
# krigeall.resids.list=parLapply(cl,resids.list,function(i) {
#   lapply(names(i@data)[16:19],function(j) {
#     as.data.frame(autoKrige(formula(paste(j,"~1")),i,grd)$krige_output) } ) } )
# 
# #Apply names for titling during plotting
# names(krigeall.resids.list)=c("abund","biomass")
# for(i in seq_along(krigeall.resids.list)) {
#   names(krigeall.resids.list[[i]])=c("species.div.resids","taxo.div.resids","func.div.resids","phylo.div.resids") }  
# 
# #Plot all components across all years and months
# lapply(seq_along(krigeall.resids.list),function(i) {
#   do.call(grid.arrange,lapply(seq_along(krigeall.resids.list[[i]]),function(j) {
#     ggplot()+
#       geom_raster(data=krigeall.resids.list[[i]][[j]],aes(x=x1,y=x2,fill=var1.pred))+
#       scale_fill_gradientn(colours=rev(rainbow(3)))+
#       coord_equal()+labs(x="",y="",title=paste(names(krigeall.resids.list[[i]])[j]))+
#       theme_bw()+theme(plot.margin=unit(c(0,0,0,0),"cm"),
#                        axis.text.x=element_blank(),axis.text.y=element_blank(),
#                        axis.ticks=element_blank(),legend.title=element_blank()) } ) ) } )

#Interpolate across by month across all years by passing each matrix to a different core 
#Perform kriging interpolation for each matrix in the list (abundance- and biomass-weighted))
krigemonth.resids.list=parLapply(cl,residssp.list,function(i) {
  lapply(names(i@data)[72:87],function(j) {
    lapply(unique((i@data)$Month),function(k) {
      as.data.frame(autoKrige(formula(paste(j,"~1")),subset(i,i@data$Month==k),grd)$krige_output) 
} ) } ) } )

#Apply names for titling during plotting
names(krigemonth.resids.list)=c("abund","biomass","pres.abs")

for(i in seq_along(krigemonth.resids.list)) {
  names(krigemonth.resids.list[[i]])=c(
    "Gini-Simpson~Richness","Taxonomic~Richness","Functional~Richness","Phylogenetic~Richness",
    "Gini-Simpson~Taxonomic","Gini-Simpson~Functional","Gini-Simpson~Phylogenetic",
    "Functional~Gini-Simpson","Functional~Taxonomic","Functional~Phylogenetic",
    "Phylogenetic~Gini-Simpson","Phylogenetic~Taxonomic","Phylogenetic~Functional",
    "Taxonomic~Gini-Simpson","Taxonomic~Functional","Taxonomic~Phylogenetic") }

for(i in seq_along(krigemonth.resids.list)) {
  for(j in seq_along(krigemonth.resids.list[[i]])) {
    names(krigemonth.resids.list[[i]][[j]])=c("Mar","May","July","Sept","Nov") } }
#Collapse lists into data frame
krigemonth.resids.list=lapply(krigemonth.resids.list,function(i) lapply(i,function(j) do.call(rbind,j)) ) 
#Split rownames into month value and bind to dataframe
krigemonth.resids.list=lapply(krigemonth.resids.list,function(i) {
  lapply(i,function(j) {
    cbind(data.frame(str_split_fixed(rownames(j),"\\.",2)[,1]),j) } ) } )
for(i in seq_along(krigemonth.resids.list)) {
  for(j in seq_along(krigemonth.resids.list[[i]])) { names(krigemonth.resids.list[[i]][[j]])[1]=c("month") } }
#Order levels of month for facetting
for(i in seq_along(krigemonth.resids.list)) {
  for(j in seq_along(krigemonth.resids.list[[i]])) { 
    krigemonth.resids.list[[i]][[j]]$month=factor(krigemonth.resids.list[[i]][[j]]$month,
                                                  levels=c("Mar","May","July","Sept","Nov")) } }

#Plot and facet by month (12.5" x 20")
lapply(seq_along(krigemonth.resids.list),function(i) {
  do.call(grid.arrange,c(lapply(seq_along(krigemonth.resids.list[[i]]),function(j) {
    ggplot()+
      geom_raster(data=krigemonth.resids.list[[i]][[j]],aes(x=x1,y=x2,fill=exp(var1.pred)))+
      scale_fill_gradientn(colours=rev(rainbow(3)))+
      facet_wrap(~month,nrow=1)+
      coord_equal()+labs(x="",y="",title=paste(names(krigemonth.resids.list[[i]])[j]))+
      theme_bw()+theme(plot.margin=unit(c(0,0,0,0),"cm"),
                       panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                       axis.text.x=element_blank(),axis.text.y=element_blank(),
                       axis.ticks=element_blank(),legend.title=element_blank()) } ), list(ncol=3)) ) } )

#Repeat but only for subset for figure in main text
figure3.subset=krigemonth.resids.list[[2]][c("Gini-Simpson~Richness","Functional~Phylogenetic","Functional~Taxonomic","Phylogenetic~Taxonomic")]
names(figure3.subset)=c("A) Gini-Simpson~Richness","B) Functional~Phylogenetic","C) Functional~Taxonomic","D) Phylogenetic~Taxonomic")
#Generate figure ("12 x 12")
do.call(grid.arrange,c(lapply(seq_along(figure3.subset),function(j) {
  ggplot()+
    geom_raster(data=figure3.subset[[j]],aes(x=x1,y=x2,fill=exp(var1.pred)))+
    scale_fill_gradientn(colours=rev(rainbow(3)))+
    facet_wrap(~month,nrow=1)+
    coord_equal()+labs(x="",y="",title=paste(names(figure3.subset)[j]))+
    theme_bw()+theme(plot.margin=unit(c(0.2,0,0,0),"cm"),
                     plot.title=element_text(size=18,hjust=0,vjust=1),
                     strip.text=element_text(size=12),
                     panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                     axis.text.x=element_blank(),axis.text.y=element_blank(),
                     axis.ticks=element_blank(),legend.title=element_blank()) } ), list(ncol=2))) 

#######################################################################################################
#                                  DRIVERS OF PATTERNS IN DIVERSITY                                   #
#######################################################################################################

#Run generalized additive models to identify the contribution of each predictor to % deviance explained
#Code modified from: https://stat.ethz.ch/pipermail/r-help/2011-November/295324.html

#Send alphadiv.list to each core
clusterExport(cl,varlist=c("alphadiv.list"))

#Compute models and store summary table of partial deviances in a list
partial.deviance.list=parLapply(cl,seq_along(alphadiv.list),function(i) {
  do.call(rbind,lapply(list("richness","evenness","species.div","func.div","phylo.div","taxo.div"),function(j) {

  #Subset only columns that are necessary to run the models and omit any NAs so all models are run on the same data
  data=na.omit(alphadiv.list[[i]][,c("Month","Year","Long","Lat","SA","WT","DO","TDEPTH","richness","evenness",
                                     "species.div","func.div","phylo.div","taxo.div")])
  
  #Build full model and null model based on (transformed) response variable
  if(j=="richness") {
    fullmod=gam(formula(paste(j,"~Month+Year+s(Long,Lat,by=Month,bs='ts')+s(SA,bs='ts')+s(WT,bs='ts')+
                  s(DO,bs='ts')+s(TDEPTH,bs='ts')")),family=poisson,data=data)
    nullmod=gam(formula(paste(j,"~1")),family=poisson,data=data)
  } else if(j=="evenness") {
    fullmod=gam(formula(paste("asin(sqrt(",j,"))~Month+Year+s(Long,Lat,by=Month,bs='ts')+s(SA,bs='ts')+s(WT,bs='ts')+
                  s(DO,bs='ts')+s(TDEPTH,bs='ts')")),data=data)
    nullmod=gam(formula(paste("asin(sqrt(",j,"))~1")),data=data)
  } else {
    fullmod=gam(formula(paste("log(",j,")~Month+Year+s(Long,Lat,by=Month,bs='ts')+s(SA,bs='ts')+s(WT,bs='ts')+
                  s(DO,bs='ts')+s(TDEPTH,bs='ts')")),data=data)
    nullmod=gam(formula(paste("log(",j,")~1")),data=data) }
  
  #Set up empty vectors for deviances
  time.dev=c(); space.dev=c(); env.dev=c()
  
  #Sequence 1: time->space->environment
  time.dev[1]=
    (deviance(update(fullmod,.~.-Month-Year,sp=fullmod$sp))-
       deviance(fullmod))/deviance(nullmod)
  space.dev[1]=
    (deviance(update(fullmod,.~.-Month-Year-s(Long,Lat,by=Month,bs='ts'),sp=fullmod$sp[6:9]))-
      deviance(update(fullmod,.~.-Month-Year,sp=fullmod$sp)))/deviance(nullmod)
  env.dev[1]=
    (deviance(nullmod)-
       deviance(update(fullmod,.~.-Month-Year-s(Long,Lat,by=Month,bs='ts'),sp=fullmod$sp[6:9])))/deviance(nullmod)

  #Sequence 2: time->environment->space
  time.dev[2]=
    (deviance(update(fullmod,.~.-Month-Year,sp=fullmod$sp))-
       deviance(fullmod))/deviance(nullmod)
  env.dev[2]=
    (deviance(update(fullmod,.~.-Month-Year-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts'),sp=fullmod$sp[1:5]))-
       deviance(update(fullmod,.~.-Month-Year,sp=fullmod$sp)))/deviance(nullmod)
  space.dev[2]=
    (deviance(nullmod)-
       deviance(update(fullmod,.~.-Month-Year-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts'),sp=fullmod$sp[1:5])))/deviance(nullmod)
       
  #Sequence 3: space->time->environment
  space.dev[3]=
    (deviance(update(fullmod,.~.-s(Long,Lat,by=Month,bs='ts'),sp=fullmod$sp[6:9]))-
       deviance(fullmod))/deviance(nullmod)
  time.dev[3]=
    (deviance(update(fullmod,.~.-s(Long,Lat,by=Month,bs='ts')-Month-Year,sp=fullmod$sp[6:9]))-
       deviance(update(fullmod,.~.-s(Long,Lat,by=Month,bs='ts'),sp=fullmod$sp[6:9])))/deviance(nullmod)
  env.dev[3]=
    (deviance(nullmod)-
       deviance(update(fullmod,.~.-s(Long,Lat,by=Month,bs='ts')-Month-Year,sp=fullmod$sp[6:9])))/deviance(nullmod)
  
  #Sequence 4: space->environment->time
  space.dev[4]=
    (deviance(update(fullmod,.~.-s(Long,Lat,by=Month,bs='ts'),sp=fullmod$sp[6:9]))-
       deviance(fullmod))/deviance(nullmod)
  env.dev[4]=
    (deviance(update(fullmod,.~.-s(Long,Lat,by=Month,bs='ts')-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts')))-
       deviance(update(fullmod,.~.-s(Long,Lat,by=Month,bs='ts'),sp=fullmod$sp[6:9])))/deviance(nullmod)
  time.dev[4]=
    (deviance(nullmod)-
       deviance(update(fullmod,.~.-s(Long,Lat,by=Month,bs='ts')-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts'))))/deviance(nullmod)
  
  #Sequence 5: environment->time->space
  env.dev[5]=
    (deviance(update(fullmod,.~.-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts'),sp=fullmod$sp[1:5]))-
       deviance(fullmod))/deviance(nullmod)
  time.dev[5]=
    (deviance(update(fullmod,.~.-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts')-Month-Year,sp=fullmod$sp[1:5]))-
       deviance(update(fullmod,.~.-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts'),sp=fullmod$sp[1:5])))/deviance(nullmod)
  space.dev[5]=
    (deviance(nullmod)-
       deviance(update(fullmod,.~.-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts')-Month-Year,sp=fullmod$sp[1:5])))/deviance(nullmod)
       
  #Sequence 6: environment->space->time
  env.dev[6]=
    (deviance(update(fullmod,.~.-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts'),sp=fullmod$sp[1:5]))-
       deviance(fullmod))/deviance(nullmod)
  space.dev[6]=
    (deviance(update(fullmod,.~.-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts')-s(Long,Lat,by=Month,bs='ts')))-
       deviance(update(fullmod,.~.-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts'),sp=fullmod$sp[1:5])))/deviance(nullmod)
  time.dev[6]=
    (deviance(nullmod)-
       deviance(update(fullmod,.~.-s(SA,bs='ts')-s(WT,bs='ts')-s(DO,bs='ts')-s(TDEPTH,bs='ts')-s(Long,Lat,by=Month,bs='ts'))))/deviance(nullmod)
  
  #Bind results into data.frame
  prop.deviance=data.frame(
    Diversity=paste(j),
    Predictor=c("All","Space","Time","Environment"),
    Partial.Deviance=c(summary(fullmod)$dev.expl,mean(space.dev),mean(time.dev),mean(env.dev)),
    Partial.Deviance.SE=c(0,std.error(space.dev),std.error(time.dev),std.error(env.dev)) )
} ) ) } )
#Close the cluster
stopCluster(cl)

#Name objects in list
names(partial.deviance.list)=c("abund","biomass","pres.abs")
#Rename diversity levels for better plotting
for(i in seq_along(partial.deviance.list)) {
  partial.deviance.list[[i]]=transform(partial.deviance.list[[i]],Diversity=recode(Diversity,
    "'richness'='Richness';'evenness'='Evenness';'species.div'='Gini-Simpson';
     'func.div'='Functional';'phylo.div'='Phylogenetic';'taxo.div'='Taxonomic'")) }
#Reorder diversity levels for better plotting
for(i in seq_along(partial.deviance.list)) {
  partial.deviance.list[[i]][,"Diversity"]=factor(partial.deviance.list[[i]][,"Diversity"],
    levels=c("Richness","Evenness","Gini-Simpson","Functional","Phylogenetic","Taxonomic")) }

#Create barplot (8" x 6")
lapply(seq_along(partial.deviance.list),function(i) {
  data=subset(partial.deviance.list[[i]],Predictor!="All")
  ggplot(data,aes(x=Diversity,y=Partial.Deviance*100,fill=Predictor))+
    geom_bar(position="dodge",color="black")+
    geom_errorbar(aes(max=(Partial.Deviance+Partial.Deviance.SE)*100,
                      min=(Partial.Deviance-Partial.Deviance.SE)*100),
                  width=0.3,size=0.6,position=position_dodge(width=0.9))+
    labs(x="",y="Total % Deviance Explained\n")+
    scale_fill_manual(values=c("grey50","grey80","white"))+
    theme_bw(base_size=18)+theme(
      panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
      legend.position="bottom",legend.title=element_blank()) } )