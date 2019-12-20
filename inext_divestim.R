# Import Required Packages
#library("png")
#library("grid")
library("tidyverse")   
library("vegan")
library("viridis")
#library("cowplot")
#library("ggrepel")
library("iNEXT")
library("broom")
#library("ggpmisc")
#library("pander")
library("lubridate")
#library("betapart")
#library("VennDiagram")

source("bin/mothur_tools.R")
se <- function(x, ...){sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))}


#Next, we'll set the aesthetics of the figures we will produce. 

my.cols <- RColorBrewer::brewer.pal(n = 4, name = "Greys")[3:4]

# Set theme for figures in the paper
theme_set(theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 15, b = 15)),
        axis.title.y = element_text(margin = margin(l = 15, r = 15)),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5)),
        #axis.line.x = element_line(size = 1),
        #axis.line.y = element_line(size = 1),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        axis.ticks.length = unit(.1, "in"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        strip.background = element_blank()
        ))

## Import Data
#Here, we read in the processed sequence files from mothur (shared and taxonomy) and a design of the sampling. We also load in the environmental data. We then remove the mock community from the dataset and ensure the the design and OTU table are aligned by row.

# Define Inputs
# Design = general design file for experiment
# shared = OTU table from mothur with sequence similarity clustering
# Taxonomy = Taxonomic information for each OTU
design <- "data/UL.design.txt"
shared <- "data/ul_resgrad.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared"
taxon  <- "data/ul_resgrad.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy"

# Import Design
design <- read.delim(design, header=T, row.names=1)

# Import Shared Files
OTUs <- read.otu(shared = shared, cutoff = "0.03")    # 97% Similarity

# Import Taxonomy
OTU.tax <- read.tax(taxonomy = taxon, format = "rdp")

# Load environmental data
env.dat <- read.csv("data/ResGrad_EnvDat.csv", header = TRUE)
env.dat <- env.dat[-c(16,17,18),]

# Subset to just the reservoir gradient sites
OTUs <- OTUs[str_which(rownames(OTUs), "RG"),]
OTUs <- OTUs[-which(rownames(OTUs) == "RGMockComm"),]

# make sure OTU table matches up with design order
design <- design[-c(34:39),]
OTUs <- OTUs[match(rownames(design), rownames(OTUs)),]
design$distance <- max(na.omit(design$distance)) - design$distance
env.dat$distance <- max(na.omit(env.dat$dist.dam)) - env.dat$dist.dam

## Clean and transform OTU table
#Here, we remove OTUs with low incidence across sites, we remove any samples with low coverage, and we standardize the OTU table by log-transforming the abundances and relativizing by site. 

# Remove OTUs with less than two occurences across all sites
#OTUs <- OTUs[, which(colSums(OTUs) >= 2)]

# Sequencing Coverage
coverage <- rowSums(OTUs)

# Remove Low Coverage Samples (This code removes two sites: Site 5DNA, Site 6cDNA)
lows <- which(coverage < 10000)
OTUs <- OTUs[-which(coverage < 10000), ]
design <- design[-which(coverage < 10000), ]
otus.for.inext <-  t(OTUs)
# Remove OTUs with < 2 occurences across all sites
OTUs <- OTUs[, which(colSums(OTUs) >= 2)]
coverage <- rowSums(OTUs)

# Rarify the community, nest RNA in DNA, and reorganize OTU table
set.seed(47405)
OTUs <- rrarefy(OTUs, min(coverage))
OTUs.w.dna <- OTUs[which(design$type == "water" & design$molecule == "DNA"),]
rowSums((OTUs.w.dna > 1))
OTUs.w.rna <- OTUs[which(design$type == "water" & design$molecule == "RNA"),]
rowSums((OTUs.w.rna > 1))

OTUs.w.dna <- OTUs.w.dna + as.matrix(decostand(OTUs.w.rna, method = "pa"))
rowSums((OTUs.w.dna > 1))
OTUs <- rbind(OTUs[1:3,],
              OTUs.w.dna,
              OTUs.w.rna)
OTUs <- OTUs[match(rownames(design), rownames(OTUs)),]

# Make Relative Abundance Matrices
OTUsREL <- decostand(OTUs, method = "total")

# Log Transform Relative Abundances
OTUsREL.log <- decostand(OTUs, method = "log")


# Reservoir environmental gradients
#Just to see if there are any strong underlying resource or nutrient gradients in the reservoir, we'll plot them along the distance of the reservoir.

#So, there are some weak gradients, but nothing too prevailing. 

# Analyze Diversity 
# Now, we will analyze the bacterial diversity in the reservoir and nearby soils to figure out how well they support different mechanisms of community assembly. 

## How does \(\alpha\)-diversity vary along the reservoir?

# First, we use the method of rarefaction and extrapolation developed by Chao et al. in the iNEXT package.

# Observed Richness
S.obs <- rowSums((OTUs > 0) * 1)

# Simpson's Evenness
SimpE <- function(x = ""){
  x <- as.data.frame(x)
  D <- diversity(x, "inv")
  S <- sum((x > 0) * 1) 
  E <- (D)/S 
  return(E)
}
simpsE <- round(apply(OTUs, 1, SimpE), 3)
shan <- diversity(OTUs, index = "shannon")
exp.shan <- exp(shan)
alpha.div <- cbind(design, S.obs, simpsE, shan, exp.shan)

# define singleton estimator from Chiu and Chao 2016 PeerJ
source("bin/Chao_functions.R")

# # estimate richness
singleton.apply <- function(x){
  singleton.Est(x, "abundance")$corrected.data
}

otus.for.inext <- apply(otus.for.inext, MARGIN = 2, singleton.apply)
divestim <- estimateD(otus.for.inext, datatype = "abundance",
           base = "size", conf = 0.95)
saveRDS(divestim, file = "intermediate-data/inext-output.rda")

