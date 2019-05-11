#' ---
#' title: " R script for analyzing RNA-seq data with edgeR "
#' output:
#'  pdf_document:
#'   keep_tex: true
#' ---

#---------------------------------------
# Remove all variables from the base environment  

rm (list = ls ())

#  Loading Packages
#install.packages("pacman")  # Packages manager (Recommended)
pacman::p_load(edgeR, limma) # DGE analysis
pacman::p_load(Glimma)       # DGE report 
pacman::p_load(openxlsx)     # Save results as excel file,
                             # any error: 'sudo apt-get install r-cran-rjava'

#---------------------------------------
# 1.- Creation of DGEList-object 

#  Matrix count
matrix_counts <- read.delim("matrix_count_geneID.txt",
                            row.names = "Geneid",
                            comment.char = "#")
matrix_counts <- matrix_counts[ ,1:8] # Keeping only root samples
head(matrix_counts)

#---
#  Matrix design
SraRunTable <- read.csv("SraRunTable.csv")
SraRunTable <- SraRunTable[1:8, ] # Keeping only root samples

group <- rep(c("ControlRoots", "RootInfect48h", "RootInfect07d", "RootInfect15d"),
             each =2)
group <- factor(group, levels= unique(group)) 

#  Annotations
annot <- read.delim("olea_annotation.csv", sep = ",")
genes <- annot[,c("EntrezID", "GeneName", "Protein", "GO_terms")]
genes <- genes[!duplicated(genes$EntrezID),]

urlNCBI <- ifelse(genes$GeneName=="NA", NA, paste(
                      "<a href='http://www.ncbi.nlm.nih.gov/gene/?term=",
                       genes$GeneName,
                       "'>",
                       genes$GeneName,
                       "</a>",
                       sep=""
                      )
                 )
urlKEGG <- ifelse(genes$EntrezID =="NA", NA, paste(
                       "<a href='https://www.genome.jp/dbget-bin/www_bget?oeu:",
                       genes$EntrezID,
                       "'>",
                       "oeu:",
                       genes$EntrezID,
                       "</a>",
                       sep=""
                       )
                 )

genes$GeneName <- urlNCBI
genes$KEGG <- urlKEGG
rownames(genes) <- genes$EntrezID

#  DGEList-object
dge <- DGEList(counts = matrix_counts[rownames(genes), ],
               group = group)
dge$genes <- genes[rownames(dge), 1:5]

dge

#---------------------------------------
# 2.- Data pre-processing

# Eliminate genes with low level of expression
dge <- dge[rowSums(cpm(dge)>1) >= 2, , keep.lib.sizes=FALSE]
dim(dge)

#---
# Collapse technical replicates
dge <- sumTechReps(dge, ID = dge$samples$group)

# Normalization of raw data by TMM method
dge <- calcNormFactors(dge, method = "TMM")
dge$samples$norm.factors

#---------------------------------------
# 3.- Differential expression analysis

# Estimation of the common dispersion from the
# housekeeping genes and all the libraries as one group:
# Source: 'Carmona, R., et.al (2017).Identification of reference genes based on RNA-seq data.
# Biomedical engineering online, 16, 65. doi:10.1186/s12938-017-0356-5'

housekeeping <- c(111372741, 111386620, 111405060, 111393675, 111397761, 111396244,
                  111405236, 111384356, 111408661, 111394484, 111396638, 111383659,
                  111412389, 111399585, 111409889, 111397323, 111379122, 111369920,
                  111385724, 111367905)

housekeeping <- as.character(housekeeping)

# checked housekeeping constant expression
logCPM <- cpm(dge, prior.count = 2, log = TRUE)
logCPM[housekeeping, ]

#---
# Estimation of common dispersion
y1 <- dge
y1$samples$group <- 1
y0 <- estimateDisp(y1[housekeeping,], trend = "none", tagwise = FALSE)

dge$common.dispersion <- y0$common.dispersion
dge$common.dispersion

#---
# Creating a design matrix and contrasts
group0 <- unique(group)
design <- model.matrix(~0 + group0)
colnames(design) <- gsub("group0", "", colnames(design))
design

contr.matrix <- makeContrasts(
  Root_Early = RootInfect48h - ControlRoots,
  Root_Late = (RootInfect15d + RootInfect07d)/2 - RootInfect48h,
  Root_InfectvsCtrl = (RootInfect48h + RootInfect07d + RootInfect15d)/3 - ControlRoots, 
  levels = colnames(design)
)

contr.matrix

#---
#  Comparison
fit <- glmFit(dge, design)

res_R_early  <- glmLRT(fit, contrast = contr.matrix[,"Root_Early"])
res_R_late  <- glmLRT(fit, contrast = contr.matrix[,"Root_Late"])
res_R_InfvsCtrl <- glmLRT(fit, contrast = contr.matrix[,"Root_InfectvsCtrl"])

# Summary
DGE_res <- topTags(res_R_late, p.value = 0.5, n = 50)
head(DGE_res$table)

#---
# Glimma report

glMDPlot(res_R_late,
         status = decideTests.DGELRT(res_R_late , p.value=0.05),
         counts=dge$counts[,c(2:4)], groups=group0[2:4],
         side.xlab="Infecction stages",
         side.main="Protein",
         transform=TRUE,
         display.columns=c("GeneName", "Protein", "KEGG"),
         cols =  c("black","grey","red")
         )

#  Saving DGE results in a excel file

list_of_datasets <- list("early-stage infection" = topTags(res_R_early, n = Inf),
                         "late-stage infection" = topTags(res_R_late, n = Inf),
                         "Common response" = topTags(res_R_InfvsCtrl, n = Inf)
                         )
write.xlsx(list_of_datasets, file = "DGE_res.xlsx", row.names=TRUE)



#---------------------------------------
# 3.- Pathways enrichment analysis

keg_R_early <- topKEGG(kegga(res_R_early, species.KEGG = "oeu", plot = F, trend = F),
                       n = 25,
                       truncate = 60)
keg_R_late <- topKEGG(kegga(res_R_late, species.KEGG="oeu",plot = F, trend = F),
                      n = 25,
                      truncate = 60)
keg_R_InfVsCtrl <- topKEGG(kegga(res_R_InfvsCtrl, species.KEGG="oeu", plot = F, trend = F),
                           n = 25,
                           truncate = 60)

list_of_datasets <- list("early-stage infection" = keg_R_early,
                         "late-stage infection" = keg_R_late,
                         "Common response" = keg_R_InfVsCtrl
                         )
write.xlsx(list_of_datasets, file = "KEGG_res.xlsx", row.names = TRUE)

#---------------------------------------
sessionInfo()
