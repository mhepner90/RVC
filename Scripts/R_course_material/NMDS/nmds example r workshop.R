
# Let's lay some conceptual groundwork first to make sure everyone is on the same page
# Consider a single axis of abundance representing a single species:
plot(0:10,0:10,type="n",axes=F,xlab="Abundance of Species 1",ylab="")
axis(1)
# We can plot each community on that axis depending on the abundance of 
# species 1 within it
points(5,0); text(5.5,0.5,labels="community A")
points(3,0); text(3.2,0.5,labels="community B")
points(0,0); text(0.8,0.5,labels="community C")

# Now consider a second axis of abundance representing a different species
# Communities can be plotted along both axes depending on the abundance of
# species within it
plot(0:10,0:10,type="n",xlab="Abundance of Species 1",
     ylab="Abundance of Species 2")
points(5,5); text(5,4.5,labels="community A")
points(3,3); text(3,3.5,labels="community B")
points(0,5); text(0.8,5.5,labels="community C")

# Now consider a THIRD axis of abundance representing yet another species
# (For this we're going to need to load another package)
#install.packages("scatterplot3d")
library(scatterplot3d)
d=scatterplot3d(0:10,0:10,0:10,type="n",xlab="Abundance of Species 1",
                ylab="Abundance of Species 2",
                zlab="Abundance of Species 3"); d
d$points3d(5,5,0); text(d$xyz.convert(5,5,0.5),labels="community A")
d$points3d(3,3,3); text(d$xyz.convert(3,3,3.5),labels="community B")
d$points3d(0,5,5); text(d$xyz.convert(0,5,5.5),labels="community C")


# There are lots of ordination such as PCA, CCA, ect.  each on of these
# has its own assumptions.  For example PCA assumes linearity, which is rare
# for community data.  NMDS, non metric multi dimmeniosnal scaling uses the
# rank order, which has few assumptions and thus a powerful way to ordinate
# communities is species space
# WAIT! species space?!?!?! What does that mean?  We plot each site/sample/point
# in terms of the relative abundances of each species.  In the example above
# we had a three species space.  You can have as many species (which become
# the axis) as you could want.  Eventhough out minds can really only imagine 
# three dimensions (for you really smart ones four dimesnion), but beyond that
# it is difficult for our minds to think about. However mathatically we can just
# keep adding species (aka axes) to figure out the ordination of the communities
# in relation to one another based the the relative abundnaces of species.
# Once we have this n-dimensional space that we can't imagine, we take a 2-D or
# 3-D projection of it so we can see it

# The use of ranks omits some of the issues associated with using absolute
# distance (e.g., sensitivity to transformation), and as a result is much 
# more flexible technique that accepts a variety of types of data
# (It is also where the "non-metric" part of the name comes from)

# The NMDS procedure is iterative and takes place over several steps:
# (1) Define the original positions of communities in multidimensional space 
# (2) Specify the number m of reduced dimensions (typically 2)
# (3) Construct an initial configuration of the samples in 2-dimensions
# (4) Regress distances in this initial configuration against the observed (measured) distances
# (5) Determine the stress (disagreement between 2-D configuration and 
# predicted values from regression)
# If the 2-D configuration perfectly preserves the original rank 
# orders, then a plot ofone against the other must be monotonically 
# increasing. The extent to which the points on the 2-D configuration
# differ from this monotonically increasing line determines the
# degree of stress (see Shepard plot)
# (6) If stress is high, reposition the points in m dimensions in the 
#direction of decreasing stress, and repeat until stress is below 
#some threshold 

# Generally, stress < 0.05 provides an excellent represention in reduced 
# dimensions, < 0.1 is great, < 0.2 is good, and stress > 0.3 provides a 
# poor representation

# NOTE: The final configuration may differ depending on the initial 
# configuration (which is often random) and the number of iterations, so
# it is advisable to run the NMDS multiple times and compare the 
# interpretation from the lowest stress solutions

# To begin, NMDS requires a distance matrix, or a matrix of dissimilarities
# Raw Euclidean distances are not ideal for this purpose: they are 
# sensitive to totalabundances, so may treat sites with a similar number
# of species as more similar, even though the identities of the species
# are different
# They are also sensitive to species absences, so may treat sites with
# the same number of absent species as more similar

# Consequently, ecologists use the Bray-Curtis dissimilarity calculation, 
# which has many ideal properties:
# It is invariant to changes in units
# It is unaffected by additions/removals of species that are not 
# present in two communities
# It is unaffected by the addition of a new community
# It can recognize differences in total abudnances when relative 
# abundances are the same

# To run the NMDS, we will use the function `metaMDS` from the vegan 
# package
#install.packages("vegan")
library(vegan)
# `metaMDS` requires a community-by-species matrix
# Let's create that matrix with some randomly sampled data
set.seed(2)
community_matrix=matrix(sample(1:100,300,replace=T),nrow=10,
                        dimnames=list(paste("community",1:10,sep=""),
                                      paste("sp",1:30,sep="")))


# The function `metaMDS` will take care of most of the distance 
# calculations, iterative fitting, etc. We need simply to supply:
example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                     k=2) # The number of reduced dimensions
# You should see each iteration of the NMDS until a solution is reached
# (i.e., stress was minimized after some number of reconfigurations of 
# the points in 2 dimensions). You can increase the number of default 
# iterations using the argument "trymax=##"
example_NMDS=metaMDS(community_matrix,k=2,trymax=100)

# And we can look at the NMDS object
example_NMDS # metaMDS has automatically applied a square root 
# transformation and calculated the Bray-Curtis distances for our 
# community-by-site matrix
# Let's examine a Shepard plot, which shows scatter around the regression
# between the interpoint distances in the final configuration (distances 
# between each pair of communities) against their original dissimilarities
stressplot(example_NMDS)
# Large scatter around the line suggests that original dissimilarities are
# not well preserved in the reduced number of dimensions

#Now we can plot the NMDS
plot(example_NMDS)
# It shows us both the communities ("sites", open circles) and species 
# (red crosses), but we  don't know which are which!

# We can use the functions `ordiplot` and `orditorp` to add text to the 
# plot in place of points
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01) ##air is the amount of empty space between text labels, <1 allows overlapping
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)

# There are some additional functions that might of interest
# Let's suppose that communities 1-5 had some treatment applied, and 
# communities 6-10 a different treatment
# We can draw convex hulls connecting the vertices of the points made by
# these communities on the plot
# First, let's create a vector of treatment values:
treat=c(rep("Treatment1",5),rep("Treatment2",5))
ordiplot(example_NMDS,type="n")
ordihull(example_NMDS,groups=treat,draw="polygon",col="grey80",
         label=FALSE)
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)
# I find this an intuitive way to understand how communities and species 
# cluster based on treatments
# One can also plot ellipses and "spider graphs" using the functions 
# `ordiellipse` and `orderspider` which emphasize the centroid of the 
# communities in each treatment


## create a plot with ellipses of various sizes
ordiplot(example_NMDS,type="n")
ordiellipse(example_NMDS,groups=treat,kind="ehull",draw="polygon",col="dodgerblue",
            label=FALSE)
ordiellipse(example_NMDS,groups=treat,kind="sd",draw="polygon",col="darkgreen",
            label=FALSE)
ordiellipse(example_NMDS,groups=treat,kind="se",draw="polygon",col="darkorchid3",
            label=FALSE)

legend(-0.25,.125,legend=c("ellipse","sd","se"),fill=c("dodgerblue","darkgreen","darkorchid3"), bty="n", cex=1)
orditorp(example_NMDS,display="species",col="1",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("12",5),rep("14",5)),
         air=0.01,cex=1.25)
# Another alternative is to plot a minimum spanning tree (from the 
# function `hclust`), which clusters communities based on their original 
# dissimilarities and projects the dendrogram onto the 2-D plot
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)
ordicluster(example_NMDS,hclust(vegdist(community_matrix,"bray"))) 
# Note that clustering is based on Bray-Curtis distances
# This is one method suggested to check the 2-D plot for accuracy


