## setting the working directory
setwd("~/Honours_2024/Garthia_DATA")

library(dartRverse)

#convert to genlight object
gl <- gl.read.dart(filename="Report_DGart24_9342_SNP_2.csv",
                   ind.metafile="Garthia_metadata_UPDATED.csv")

#save genlight object for future reference
saveRDS(gl, file="R_Garthia_data_populations")

gl <- gl.load("R_Garthia_data_populations")

#### Explore Data ####
gl@other$loc.metrics$RepAvg[1:10]

names(gl@other$loc.metrics)

gl@pop

#Only first 10 entries shown
gl@other$ind.metrics$lat[1:10]

#number of individuals
nInd(gl)

#number of loci
nLoc(gl)

#number of populations
nPop(gl)

#names of populations
popNames(gl)

#names of individuals
indNames(gl)

#names of loci
locNames(gl)

#table of populations
table(pop(gl))

#set verbosity
gl.set.verbosity(3)

#visualise gaps in data 
gl.smearplot(gl)

#gl.report.callrate(gl, method='ind', ind.to.list=44)
#checking to see call rate for particular individuals


#### Filtering ####

gl1 <- gl

gl.report.reproducibility(gl1)
gl2 <- gl.filter.reproducibility(gl1, threshold = 0.99)
# Get rid of unreliable loci

gl.report.rdepth(gl2)
gl3 <- gl.filter.rdepth(gl2, lower = 3, upper = 100)
# Get rid of low and super high read depth loci
## DAMIEN COMMENT: Once you have the analyses figured, have a little play to see if increasing the minimum read deoth a bit, to like 5 or 7, changes anything.
## In theory wanna be stringent with this to ensure reliable SNPs, but if you look at the report, a filtering a minimum rwad depth of 3 of only lose abut 3000 SNPs, but if you increase to 5 you lose about 12.000! So thats why its improtant to test the impact of these filters on your results


gl.report.callrate(gl3) 
gl4 <- gl.filter.callrate(gl3, method = "loc", threshold = 0.4)
# Get rid of really poorly sequenced loci, but don’t cut hard

gl.report.maf(gl4)
gl5 <- gl.filter.maf(gl4,  threshold = 1/(2*nInd(gl4)))
#formula for the number of individuals

gl.report.secondaries(gl5)
gl6 <- gl.filter.secondaries(gl5)
# Always do this as the last loci filter so that you’ve cut for quality before you cut because there are two SNPs

gl.report.callrate(gl6,method="ind") 
gl7 <- gl.filter.callrate(gl6, method = "ind", threshold = 0.4)
# low filter – this is only to get rid of really bad individuals
# this does not remove any individuals currently  

gl8 <- gl.filter.monomorphs(gl7)
# Always run this after removing individuals – removes loci that are no longer variable

gl.report.callrate(gl8) 
gl9 <- gl.filter.callrate(gl8, method = "loc", threshold = 0.6)
#more stringent since we still have 90,000 loci

#show smearplot
gl.smearplot(gl9)

#finalised gl
gl <- gl9

saveRDS(gl9, file = "Garthia_SNPs_FILTERED.rds")

gl <- gl.load("Garthia_SNPs_FILTERED.rds")

#### Visualise on a map ####

## quicker way to do that
install.packages("leaflet")
install.packages("leaflet.minicharts")
library(leaflet)
library(leaflet.minicharts)
install.packages("leaflet.providers")
library(leaflet.providers)

latlong_df <- data.frame(lat = gl$other$ind.metrics$lat,

                         long = gl$other$ind.metrics$long)

colnames(latlong_df) <- c("lat", "lon")




gl@other$latlon <- data.frame(lat=gl@other$ind.metrics$lat, lon=gl@other$ind.metrics$lon)

map <- gl.map.interactive(gl, provider = "Esri.WorldTopoMap")

map


# Load the necessary libraries
library(leaflet)
library(leaflet.minicharts)
library(leaflet.providers)

# Assuming gl@other$ind.metrics$lat and gl@other$ind.metrics$long exist and contain latitude and longitude
latlong_df <- data.frame(lat = gl$other$ind.metrics$lat,
                         lon = gl$other$ind.metrics$long)

# Add these coordinates to the gl object (if required)
gl@other$latlon <- latlong_df

# Create a leaflet map using the "Esri.WorldTopoMap" provider
map <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addCircleMarkers(data = latlong_df, lat = ~lat, lng = ~lon, radius = 5, color = "blue", popup = ~paste("Lat:", lat, "<br>Lon:", lon))

# Display the map
map










?gl.map.interactive

#### Subsampling ####

#just penai clade
south <- gl.keep.pop(gl, pop.list = c("penai_combarbala", "penai_illapel", "gaudichaudii_losmolles", mono.rm = T))

table(pop(south))

#just gaudi clade 
north <- gl.drop.pop(gl, pop.list = c("penai_combarbala", "penai_illapel", "gaudichaudii_losmolles", mono.rm = T))

table(pop(north))

#### Pop Structure analyses ####
library("gplots")
library("ggplot2")
library(HardyWeinberg)
library(devtools)
library("ggtern")


# Calculate an Euclidean Distance Matrix on individuals
D <- gl.dist.ind(gl)

# Represent visually as a heat map
gl.plot.heatmap(D)

#report heterozygosity
glHe <- gl.report.heterozygosity(gl)

#hardy weinberg 
hwe <- gl.report.hwe(gl, multi_comp = TRUE)

#pairwise FST
#100 bootstraps
FST_boot <- gl.fst.pop(gl, nboots = 100, percent = 95, nclusters = 1, verbose = NULL)

#explore results
FST <- as.matrix(as.dist(FST_boot$Fsts))
P_values <- as.matrix(as.dist(FST_boot$Pvalues))
boot <- as.matrix(FST_boot$Bootstraps)

#export results
write.csv(FST, "FST.csv")
write.csv(P_values, "FST_Pvalues.csv")
write.csv(boot, "FST_bootstraps.csv")


##isolation by distance analysis (mantel test)###

latlong_df <- data.frame(lat = north$other$ind.metrics$lat,
                         long = north$other$ind.metrics$long)
colnames(latlong_df) <- c("lat", "lon")

mantel <- gl.ibd(north, distance = "euclidean", coordinates = latlong_df, permutations = 1000, paircols = 'pop')

?gl.ibd

#### PCA ####
library("plotly")
library("ggplot2")
library("ggbiplot")

#PCA analysis
pc <- gl.pcoa(gl)

#PCA plot
gl.pcoa.plot(pc,gl,interactive=F)

# Plot the first and third dimensions of the PCA
gl.pcoa.plot(pc,gl,xaxis=1,yaxis=3)

# or you can do a 3D plot
gl.pcoa.plot.3d(pc,gl)




#### EXPORT DATA FOR PHYLOGENETICS #####

#TO CREATE PARTITIONS FILE FOR SVDQ (SENT BY CARLOS)
options(scipen = 999)
trim.seqs <- gl$other$loc.metrics$TrimmedSequence
trim.seqs <- as.vector(trim.seqs)
loc.l <- nchar(trim.seqs,type="chars")
par.t <- data.frame(loc.l)
par.t$beg[[1]] <- 1
par.t$end[[1]] <- loc.l[[1]]
par.t$part[[1]] <- paste("DNA, L",1," = ",par.t$beg[[1]],"-",par.t$end[[1]],sep="")
for (i in 2:length(loc.l)){
  par.t$beg[[i]] <- par.t$end[[i-1]]+1
  par.t$end[[i]] <- par.t$beg[[i]]-1+par.t$loc.l[[i]]
  par.t$part[[i]] <- paste("DNA, L",i," = ",par.t$beg[[i]],"-",par.t$end[[i]],sep="")
}
write(unlist(par.t$part),file="partitions.txt")



#for IQ-Tree
tmp <- gl2fasta(gl, outfile = "GarthiaAll.fas", method=1, outpath = getwd()) #use method=1 to get the full comncatenated sequence tags and heterozygous sites coded as ambiguity codes rather than choosing one random allele.

tmp <- gl2fasta(gl, outfile = "GarthiaAllSNPs.fas", method=3, outpath = getwd()) #use method=3 to get the full comncatenated SNPs and heterozygous sites coded as ambiguity codes rather than choosing one random allele. Using these data with ascertainment bias correction on iqtree is much faster than using the full sequences

#For SVDQuartets
tmp <- gl2svdquartets(gl, outfile = "Garthia_svd.nex", outpath = getwd(), method=2, v=3)

#tmp <- gl2svdquartets(gl, outfile = "Garthia_svd_TEST.nex", outpath = getwd(), method=1, v=3)


#For SNAPP/snapper

library(adegenet)
library(ape)
library(phangorn)

gl2snapp(gl, outfile = "Garthia_snapp.nex", outpath = getwd(), verbose = NULL)


####### EXPORT FOR BPP #########

library(adegenet)

gl_1000 <- gl.subsample.loc(gl, 1000, replace = TRUE, error.check = TRUE, verbose = NULL)

tmp <- gl2bpp(gl_1000, method=1, outfile = "output_bpp.txt", imap = "Imap_bpp.txt", outpath = getwd())




#### sNMF analysis ####
#takes quite a while to run

#drop outgroup
gl <- gl.drop.pop(gl, "homonota")
gl <- gl.filter.monomorphs(gl)
popNames(gl)

#install LEA
BiocManager::install("LEA")
library("LEA")

# LEA requires the genotype style file
gl2faststructure(gl, outfile = "gl_structure.fstr", outpath = "C:\\Users\\Kamryn\\Documents\\Honours_2024\\Garthia_DATA")

struct2geno("gl_structure.fstr", ploidy = 2, FORMAT = 2)
###this hates any loci with all heterozygotes


#run snmf across a range of tolerance and alpha values to determine best parameter set.
#this uses a 5% masked genotypes in each iteration, but could also test a range of masked values
#to do this you have to copy and rename the gl9_fstr.geno file for the different runs
#if you don't do this is just overwrites everything for each run. Has to end in ".geno"
all_gl.snmf_1 = snmf("gl_structure.fstr1.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.00001, alpha = 10,  percentage = 0.05) 
all_gl.snmf_2 = snmf("gl_structure.fstr2.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.00001, alpha = 0,  percentage = 0.05) 
all_gl.snmf_3 = snmf("gl_structure.fstr3.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.00001, alpha = 1,  percentage = 0.05) 
all_gl.snmf_4 = snmf("gl_structure.fstr4.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.00001, alpha = 100,  percentage = 0.05) 
all_gl.snmf_5 = snmf("gl_structure.fstr5.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.00001, alpha = 1000,  percentage = 0.05) 
all_gl.snmf_6 = snmf("gl_structure.fstr6.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.00001, alpha = 10000,  percentage = 0.05) 

all_gl.snmf_7 = snmf("gl_structure.fstr.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.0001, alpha = 10,  percentage = 0.05) 
all_gl.snmf_8 = snmf("gl_structure.fstr8.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.0001, alpha = 0,  percentage = 0.05) 
all_gl.snmf_9 = snmf("gl_structure.fstr9.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.0001, alpha = 1,  percentage = 0.05) 
all_gl.snmf_10 = snmf("gl_structure.fstr10.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.0001, alpha = 100,  percentage = 0.05) 
all_gl.snmf_11 = snmf("gl_structure.fstr11.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.0001, alpha = 1000,  percentage = 0.05) 
all_gl.snmf_12 = snmf("gl_structure.fstr12.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.0001, alpha = 10000,  percentage = 0.05) 

all_gl.snmf_13 = snmf("gl_structure.fstr13.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.000001, alpha = 10,  percentage = 0.05) 
all_gl.snmf_14 = snmf("gl_structure.fstr14.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.000001, alpha = 0,  percentage = 0.05) 
all_gl.snmf_15 = snmf("gl_structure.fstr15.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.000001, alpha = 1,  percentage = 0.05) 
all_gl.snmf_16 = snmf("gl_structure.fstr16.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.000001, alpha = 100,  percentage = 0.05) 
all_gl.snmf_17 = snmf("gl_structure.fstr17.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.000001, alpha = 1000,  percentage = 0.05) 
all_gl.snmf_18 = snmf("gl_structure.fstr18.geno", K = 1:12, entropy = T, ploidy = 2, project="new", repetitions = 10, tolerance = 0.000001, alpha = 10000,  percentage = 0.05) 



#this is to get the lowest cross-entropy value from each set to figure out which values are best
#have to manually look at which K is best for each iteration
ce_estimation = all_gl.snmf_18
summary(ce_estimation)


#run 1 0.2408441
#run 2 0.2409724   
#run 3 0.2415422   
#run 4 0.2443424   
#run 5 0.2446891   
#run 6 0.2463506  
#run 7 0.2454262 
#run 8 0.2447721  
#run 9 0.2422795   
#run 10 0.2404922    
#run 11 0.2417642  
#run 12 0.2453382   
#run 13 0.2412626   
#run 14 0.2434475   
#run 15 0.2409264   
#run 16 0.2423933 
#run 17 0.2436080   
#run 18 0.2441303   



#run 10 0.2404922    
#tolerance = 0.0001, alpha = 100
#lowest cross entropy 

# sNMF <-load.snmfProject("gl_structure.fstr1.snmfproject")

ce_estimation = all_gl.snmf_10
summary(ce_estimation)
plot(ce_estimation)

#K=7


##### Run full run with the best cross entropy (ce) value

plot(all_gl.snmf_10)
summary(all_gl.snmf_10) #k = 7 has the lowest mean
ce_7 = cross.entropy(all_gl.snmf_10, K = 7)
best_7 <- which.min(ce_7) ## run 14 has the lowest CE

barplot(t(Q(all_gl.snmf_10, K = 7, run = best_7)), col = 1:7)  ### plot admixture


k7 <- as.data.frame(Q(all_gl.snmf_10, K = 7, run = best_7))
k7$ind <- gl$ind.names
k7$pop <- gl$pop    
k7$lat <- gl$other$ind.metrics$lat
k7$long <- gl$other$ind.metrics$long

k7$lon <- gl@other$latlong$lon
k7 <- k7[,c(8:11, 1:7)]

colnames(k7) <- c("ind", "pop", "lat", "long", "P1", "P2", "P3", "P4", "P5", "P6", "P7")

write.csv(k7, "Garthia_k7.csv", row.names = F)


##### AFTER WRITING THESE RESULTS AND ADDING COLUMNS OF CLADES WE READ THEM AGAIN

snmf7 <- read.csv("Garthia_k7.csv")

snmf7ord <- snmf7[order(snmf7$pop),]

XX <- t(as.matrix(snmf7ord[,5:11])) # transpose the matric for the barplot

library(RColorBrewer)
display.brewer.all()
colours_7 = c("#E69F00", "#F71945", "#F7F752", "#0072B2", "#F980F9", "#009E73", "#73EEF5")


#barplot(XX, col=brewer.pal(n = 7, name = "Set1"), names.arg = snmf7ord$ind, las=2, cex.names = 0.5, legend = T)

# Adjust the margins to fit long labels
par(mar = c(12.5, 3, 3, 0))  # Increase bottom margin to 8 lines


barplot(XX, col= colours_7, names.arg = snmf7ord$ind, las=2, cex.names = 0.5)

###### plot on map:

library(ggplot2)
library(scatterpie)
library(raster)
library(ggnewscale)


world <- map_data('world')

?map_data

## first make sure no two samples have identical coordinates, otherwise they will be combined into one pie. If two have same coordiantes just add a small decimal to one of them

k7_all = read.csv("Garthia_K7.csv")

k7_all <- k7
k7 <- k7_all




k7_all_plot <- ggplot(world, aes(long, lat)) +
  lims(x = c(-75.5, -66.5), y = c(-33, -24))+
  geom_map(map=world, aes(map_id=region), fill="grey", color="black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_quickmap()+ 
  geom_scatterpie(aes(x=long, y=lat),
                  data=k7, cols=c("P1", "P2", "P3", "P4", "P5", "P6", "P7"), alpha=1, sorted_by_radius = F, pie_scale = 8)+
  scale_fill_manual(values=colours_7 )+
  theme(legend.position='none')

k7_all_plot 



########## FIXED DIFFERENCE ANALYSIS ###########

# Let's revisit the full dataset and look for structure, using fixed difference analysis
# to do this, need to remove all individuals identified as hybrids in the snmf to be 
# testing specifically for whether, outside of hybrids, there is a measure of admixture
# back into parental groups

fd <- gl.fixed.diff(gl,test=T)
fd


D <- as.dist(fd$pcfd)
gl.plot.heatmap(D, values=FALSE)
gl.plot.heatmap(D, values=FALSE, rank=TRUE)



#fixdif <- fd$fd
#class(fixdif) <- "data.frame"

write.csv(fd$pval, "FixedDifferencesGarthiaPValues.csv")
write.csv(fd$expfpos, "FixedDifferencesGarthiaExpectedFalsePositives.csv")
write.csv(as.matrix(fd$fd), "FixedDifferencesGarthia.csv")




#### REORGANISING CSV FOR UPDATED STRUCTURE PLOT ####

test_data <- read.csv("Garthia_k7.csv")

test_data <- test_data[, -c(3, 4)]

library("reshape2")

melted_df <- melt(test_data, id.vars = c("ind", "pop"), measure.vars = c("P1", "P2", "P3", "P4", "P5", "P6", "P7"))

sorted_df <- melted_df[order(melted_df$ind), ]

reordered_df <- sorted_df[, c(1, 3, 4, 2)]

reordered_df$ind <- paste(reordered_df$ind, reordered_df$pop, sep = "_")

write.csv(reordered_df, file = "FIXED_data.csv", row.names = F)




############# MAKE TREE FIGURE WITH ANCESTRY COEFFICIENTS ##########

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#
BiocManager::install("ggtree")

library(ape)
library(ggtree)
library(picante)
library(ggstance)
library(RColorBrewer)
library(ggplot2)



#### UPDATED STRUCTURE PLOT - IQTREE ####

#read IQTREE file
tree <- ladderize(read.tree("GarthiaAll.fas.TREEFILE"))

plot(tree)


#root tree (IQTREE)
rooted_tree <- root(tree, outgroup = c("Garthia_gaudichaudii_PlayalosMolles_SVN_1262_D710_gaudichaudii_losmolles", "Garthia_gaudichaudii_PlayalosMolles_SVN_1263_D711_gaudichaudii_losmolles", "Garthia_gaudichaudii_PlayalosMolles_SVN_1264_D712_gaudichaudii_losmolles", "Garthia_penai_MontePatria_SVN_1256_D704_penai1", "Garthia_penai_MontePatria_SVN_1257_D705_penai1", "Garthia_penai_Combarbala_SVN_1258_D706_penai1", "Garthia_penai_Combarbala_SVN_1261_D709_penai1", "Garthia_penai_Combarbala_SVN_1260_D708_penai1", "Garthia_penai_Illapel__D866_penai2", "Garthia_penai_Illapel__D867_penai2", "Garthia_penai_Illapel__D870_penai2", "Garthia_penai_Illapel__D868_penai2", "Garthia_penai_Illapel__D869_penai2"))

plot(rooted_tree, cex=0.3)

#drop outgroup (IQTREE)
#tree <- drop.tip(rooted_tree, c("Homonota_cf_darwinii_Heucu_PACH575_D1381_homonota", "Homonota_cf_darwinii_Heucu_PACH576_D1382_homonota", "Homonota_cf_darwinii_Heucu_PACH577_D1325_homonota"))

plot(tree)

#read in ancestry coefficients (IQTREE)
data <- read.csv("FIXED_data.csv")


data1 <- data[,1:3]


# NOT WORKING CURRENTLY
display.brewer.all()
a<- brewer.pal(n = 7, name = "Set1")

p <- ggtree(tree, branch.length = "none", right=T) + theme(legend.position='none') + geom_tiplab()
p2 <- facet_plot(p, panel = 'bar', data = data1, geom = geom_barh,
                 mapping = aes(x = value, fill = as.factor(variable), color=a),
                 width = 0.8, stat='identity')
facet_widths(p2, widths = c(1, 2))
p2


# CHATGPT VERSION WITHOUT DEFINING COLOUR

# Create the tree plot with tip labels
p <- ggtree(rooted_tree, branch.length = "none", right = TRUE) +
  theme(legend.position = 'none', plot.margin = margin(1, 1, 1, 1)) +
  geom_tiplab(size = 2, align = TRUE) +  coord_cartesian(clip = 'off')

# Optionally, you can rescale branch lengths for better visibility
p <- p + scale_x_continuous(expand = expansion(mult = c(0.05, 0.55)))  # Adjust expansion to add space

#p <- p + scale_y_continuous(expand = expansion(mult = c(-0.1, -0.1)))


plot(p)


# Create the faceted plot with a stacked bar chart
p2 <- facet_plot(p, panel = 'bar', data = data1, geom = geom_barh,
                 mapping = aes(x = value, fill = as.factor(variable), color = as.factor(variable)),
                 width = 0.8, stat = 'identity') +
  scale_color_manual(values = rep("black", length(unique(data1$variable))))

# Adjust facet widths if needed
facet_widths(p2, widths = c(1, 1))



plot(p2)

p2



#### UPDATED STRUCTURE PLOT - SVDQ ####

#read SVDQ file
tree <- ladderize(read.nexus("Garthia_svd.tre"))

plot(tree)

#root tree (SVDQ)
rooted_tree <- root(tree, outgroup = c("Garthia_penai_Illapel__D866", "Garthia_penai_Illapel__D867", "Garthia_penai_Illapel__D868", "Garthia_penai_Illapel__D869", "Garthia_penai_Illapel__D870", "Garthia_penai_MontePatria_SVN_1256_D704", "Garthia_penai_MontePatria_SVN_1257_D705", "Garthia_penai_Combarbala_SVN_1258_D706", "Garthia_penai_Combarbala_SVN_1261_D709", "Garthia_penai_Combarbala_SVN_1260_D708", "Garthia_gaudichaudii_PlayalosMolles_SVN_1262_D710", "Garthia_gaudichaudii_PlayalosMolles_SVN_1263_D711", "Garthia_gaudichaudii_PlayalosMolles_SVN_1264_D712"))

plot(rooted_tree, cex=0.3)

#drop outgroup (SVDQ)
#tree <- drop.tip(rooted_tree, c("LpapuanaABTC68240", "LpapuanaABTC68241", "LpapuanaABTC68242"))

tree <- rooted_tree

plot(tree)

#read in ancestry coefficients (SVDQ)
data <- read.csv("FIXED_data.csv") 

data1 <- data[,1:3]


# CHATGPT VERSION WITHOUT DEFINING COLOUR

# Create the tree plot with tip labels
p <- ggtree(tree, branch.length = "none", right = TRUE) +
  theme(legend.position = 'none', plot.margin = margin(1, 1, 1, 1)) +
  geom_tiplab(size = 2, align = TRUE) +  coord_cartesian(clip = 'off')

# Optionally, you can rescale branch lengths for better visibility
p <- p + scale_x_continuous(expand = expansion(mult = c(0.05, 0.55)))  # Adjust expansion to add space

#p <- p + scale_y_continuous(expand = expansion(mult = c(-0.1, -0.1)))


plot(p)


# Create the faceted plot with a stacked bar chart
p2 <- facet_plot(p, panel = 'bar', data = data1, geom = geom_barh,
                 mapping = aes(x = value, fill = as.factor(variable), color = as.factor(variable)),
                 width = 0.8, stat = 'identity') +
  scale_color_manual(values = rep("black", length(unique(data1$variable))))



# Adjust facet widths if needed
facet_widths(p2, widths = c(1.5, 1))



plot(p2)

p2




#### compare Phylo ####
library(ape)
library(ggtree)
library(treeio)

x <- read.nexus("Garthia16S.TRE")
y <- read.nexus("GarthiaND2.TRE")


comparePhylo(x, y, plot = TRUE, force.rooted = FALSE,
             use.edge.length = FALSE, commons = TRUE,
             location = "bottomleft")

## S3 method for class 'comparePhylo'
print(x, ...)


merge_tree(x, y)

#### elevation map ####
# Load necessary libraries
library(ggplot2)
library(raster)
library(rnaturalearth)
library(elevatr)
library(rnaturalearthdata)
library(ggplot2)
library(scatterpie)
library(raster)
library(ggnewscale)


# Read your CSV file
k7_all <- read.csv("Garthia_K7.csv")
k7 <- k7_all

# Overlay scatter pie charts on the elevation map
elev_map_with_pies <- elev_map +
  geom_scatterpie(aes(x=long, y=lat),
                  data=k7, cols=c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
                  alpha=1, sorted_by_radius = FALSE, pie_scale = 8) +
  scale_fill_manual(values=colours_7) +
  theme(legend.position='none')

# Display the final map
elev_map_with_pies


colours_7 = c("#E69F00", "#F71945", "#F7F752", "#0072B2", "#F980F9", "#009E73", "#73EEF5")




# Get a world map and focus on a specific region
chile <- ne_countries(scale = "medium", type = "countries", continent = "South America", returnclass = "sf")


# Download elevation data for the region
elev_data <- get_elev_raster(locations = chile, z = 6, clip = "locations")

# Convert the elevation raster data to a data frame
elev_df <- as.data.frame(elev_data, xy = TRUE)

# Rename the columns to something more intuitive
colnames(elev_df) <- c("long", "lat", "elevation")


# Convert raster to data frame
elev_map <- ggplot() +
  geom_raster(data = elev_df, aes(x = long, y = lat, fill = elevation)) +
  scale_fill_gradient(low = "black", high = "white") +  # Change to black and white scale
  geom_sf(data = chile, fill = NA, color = "black") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#D9F1F9", color = NA)) +  # Change ocean color to blue
  coord_sf(xlim = c(-75.5, -66.5), ylim = c(-33, -24))

# Display the map
elev_map
