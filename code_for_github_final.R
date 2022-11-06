##########################################
#BASIC FILTER/data cleaning and factor declaration

library(BiocManager)
library(Biobase)
library(BiocGenerics)
library(dplyr)
library(data.table)
library(edgeR)
library(limma)

#Read the published data
reinfectionSR <- as.matrix(fread("GSE206397_countsAll.txt"),
                           rownames=1)

#samples to keep - those who were reinfected only !@!
idx_keep <- c(4:7, 9:12, 16:19, 20:23, 24:27, 28:30, 34:37)

reinfectionSR <- reinfectionSR[,idx_keep]
#Function for basic filtering 
dgeWay <- function(count_mat){
  library(edgeR)
  
  dgee <- DGEList(count_mat)
  #dgee$samples$lib.size <- colSums(dgee$counts)
  #dgee <- calcNormFactors(dgee)
  keep <- rowSums(cpm(dgee) > 1) >= 4
  dgee <- dgee[keep, ,keep.lib.sizes = FALSE]
  
  dgee$samples$lib.size <- colSums(dgee$counts)
  dgee <- calcNormFactors(dgee)
  
  return(dgee)
}
##

SR_dge <- dgeWay(reinfectionSR %>% as.matrix)

#FACTOR DECLERATION
####################
pID <- c(rep("SR1",4), rep("SR2",4), rep("SR3",4), rep("SR5",4), rep("SR6",4),
         rep("SR7",3), rep("SR",4))

casetime <- c(rep(c("NI","EA","LA","FU"),5), "NI","EA","LA",
              "NI","EA","LA","FU") %>% as.factor()

viremia <- c("min","plus","min","min",
             "min","plus","min","min",
             "min","plus","plus","min",
             "min","plus","min","min",
             "min","plus","min","min",
             "min","plus","plus",
             "min","plus","min","min") 

viremia <- factor(viremia, levels = c("min", "plus"))
casetime <- factor(casetime, levels = c("NI", "EA", "LA", "FU"))

treat <- paste(casetime, viremia, sep = "_") 
treat <- factor(treat, levels = c("NI_min","EA_plus","LA_min","LA_plus","FU_min"))

#Put it into the dge object for design matrix and downstream analysis
SR_dge$samples$patientID <- as.factor(pID)
SR_dge$samples$treat <- as.factor(treat)

#Make design matrix
designo <- model.matrix(~patientID + treat, data = SR_dge$samples)

#######################
#VOOM VOOM VOOM 
#######################
library(limma)
library(Matrix)
library(lme4)
library(zoo)
library(multcomp)
library(MASS)
library(broom)
library(sva)

patientID <- SR_dge$samples$patientID 
treat <- SR_dge$samples$treat

#use weights from voom for building the mixed weighted linear model::
vSR <- voom(SR_dge, design = designo, plot = TRUE)

#################
#PCA PCA PCA
library(ggplot2)
library(factoextra)
#################
# Remove patientID batch effect
expr.set <- removeBatchEffect(vSR$E, 
                              batch = pID)

# Select top x most variable genes
num.genes <- 500

rv <- rowVars(expr.set)

select <- order(rv,decreasing=TRUE)[seq_len(num.genes)]

# Plot top x most variable genes
pca <- prcomp(t(expr.set[select,]), scale.=TRUE)

p.pca <- ggbiplot(pca,
                  groups=treat,
                  choices=c(1,2),
                  circle=FALSE,
                  ellipse=TRUE,
                  var.axes=FALSE
)

# Edit ggbiplot object to make aesthetics work
p.pca$layers[[1]] <- NULL

p.pca + 
  geom_point(aes(color= factor(treat), shape=factor(viremia)), size = 4, alpha = 0.8) +
  ggtitle(label="Top 500 most variable genes (expressed)") +
  theme_bw() +
  coord_fixed()
