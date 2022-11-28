# Final R Code for the Re-infection article 

# Load Libraries ####################################
suppressPackageStartupMessages(library(dplyr) )
suppressPackageStartupMessages(library(assertthat)) 
suppressPackageStartupMessages(library(Biobase) ) 
suppressPackageStartupMessages(library(limma) ) 
suppressPackageStartupMessages(library(edgeR) )
suppressPackageStartupMessages(library(openxlsx) ) 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(qusage))
suppressPackageStartupMessages(library(BiocGenerics))
suppressPackageStartupMessages(library(matrixStats)) 
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(circlize))



# END ##############################################################################################################################################################################################

# Data Normalisation ##################################################################################################

#Enter here directory with your data::
load("C:/Users/User/Downloads/eset.Rdata")

dgeList <- DGEList(counts = exprs(esetRaw))
# of Ge that zero counts across all samples
table(rowSums(dgeList$counts==0)==ncol(dgeList$counts))

# Remove any ge that are expressed at levels less than 1 CPM in less than 3 samples
cpm.dgeList <- cpm(dgeList)
keep.exprs  <- rowSums(cpm.dgeList>1)>=3
dgeList     <- dgeList[keep.exprs,, keep.lib.sizes=FALSE]
rm(cpm.dgeList, keep.exprs)

# Norm 
dgeList <- edgeR::calcNormFactors(dgeList)
FC=factor(pData(esetRaw)$CClass,levels=unique(pData(esetRaw)$CClass))
blo=factor(pData(esetRaw)$"Patient.ID")
design=model.matrix(~-1+blo+FC)
vAll=voom(dgeList,design=design,plot=FALSE) 
rownames(vAll$weights)=rownames(vAll$E)
colnames(vAll$weights)=colnames(vAll$E)
fdataN=fData(esetRaw)[rownames(vAll$E),]
fdataN=fData(esetRaw)[rownames(vAll$E),]
phenoData   <- new("AnnotatedDataFrame", data = pData(esetRaw))
featureData <- new("AnnotatedDataFrame", data =fdataN)

#remove batch effect of different sequencing times
batch=pData(esetRaw)$Sequencing.Batch

mm=sva::ComBat(vAll$E,batch=as.character(batch))
esetNorm <- ExpressionSet(assayData = mm,
                          phenoData = phenoData,featureData=featureData)

# END ############################################################################################################################################################################################ 

# Figure 1B | PCA ##################

exprs.set<-exprs(esetNorm)
pheno.all<-pData(esetNorm)

View(pheno.all)
# Samples of interests

# Resolvers Secondary Infection
resII <- c("SR1.4", "SR1.5", "SR1.6", "SR1.7", 
           "SR2.2","SR2.3","SR2.4","SR2.5", 
           "SR3.4", "SR3.5", "SR3.6", "SR3.7",
           "SR5.1", "SR5.2", "SR5.3", "SR5.4",  
           "SR6.1", "SR6.2", "SR6.3", "SR6.4", 
           "SR7.1", "SR7.2","SR7.3",
           "CI1.4", "CI1.5", "CI1.6", "CI1.7")

expr.set <- exprs.set[,resII]
pheno <- pheno.all[resII,]

pheno$Class <- c("Pre-reinfection", "Early acute (+)", "Late acute (-)", "Follow up",
                 "Pre-reinfection", "Early acute (-)", "Late acute (-)", "Follow up",
                 "Pre-reinfection", "Early acute (+)", "Late acute (+)", "Follow up",
                 "Pre-reinfection", "Early acute (+)", "Late acute (-)", "Follow up",
                 "Pre-reinfection", "Early acute (+)", "Late acute (-)", "Follow up",
                 "Pre-reinfection", "Early acute (+)", "Late acute (+)",
                 "Pre-reinfection", "Early acute (+)", "Late acute (-)", "Follow up")


# Remove atientID batch effect
expr.set <- removeBatchEffect(expr.set, batch = pheno$Patient.ID)

# Select top x most variable genes

num.genes <- 500

rv <- rowVars(expr.set)

select <- order(rv,decreasing=TRUE)[seq_len(num.genes)]

# Plot top 500 most variable genes

pca <- prcomp(t(expr.set[select,]), scale.=TRUE)


p.pca <- ggbiplot(pca,
                  groups=pheno$Class,
                  choices=c(1,2),
                  circle=FALSE,
                  ellipse=TRUE,
                  var.axes=FALSE
)



# Edit ggbiplot object to make aesthetics work
p.pca$layers[[1]] <- NULL

p.pca + 
  geom_point(aes(color= factor(pheno$Class), 
                 shape=factor(pheno$Viral.load)), 
             size = 4, alpha = 0.8) +
  ggtitle(label="Top 500 most variable genes (expressed)") +
  theme_bw() +
  coord_fixed() 
