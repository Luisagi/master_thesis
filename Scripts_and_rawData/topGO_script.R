#' ---
#' title: "  "
#' output:
#'  pdf_document:
#'   keep_tex: true
#' ---


# Packages
pacman::p_load(topGO)
pacman::p_load(readxl)

# Load file data with genes expresion analysis
# Sheets: Common response, early-stage infection, late-stage infection.
data_raw <- read_excel("DGE_res.xlsx", sheet = "Common response")

# Remove rows without GO terms
data <- data_raw[!is.na(data_raw$GO_terms),]

# topGOdata object 

mygene2GO <- as.list(strsplit(as.character( data$GO_terms),","))
mygene2GO <- setNames(mygene2GO, data$EntrezID)

# Name numeric vector, list of interesting genes where the identifiers
# are stored in the names attribute of the vector.

geneList <- as.vector(data$FDR)
names(geneList) <- data$EntrezID

# Building the topGo object

# Function to select the significant genes provided
topDiffGenes <- function(allScore){
    return(allScore < 0.05)
    }
sum(topDiffGenes(geneList))

# ontology: MF, BP and CC
GOdata <- new("topGOdata",
              description = "Olea europaea time course RNA-seq analysis",
              ontology = "BP",
              allGenes = geneList,
              geneSel = topDiffGenes,
              annot = annFUN.gene2GO,
              gene2GO = mygene2GO,
              nodeSize = 10)
GOdata

# Running the enrichment test based on gene scores (p-value)
# using the Kolmogorov-Smirnov test. We will use the both the classic and the weight01 method.

resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "ks")

# Analysis of results

allRes <- GenTable(GOdata,
                      classicKS = resultKS,
                      weight01 = resultKS.weight01,
                      orderBy = "weight01",
                      topNodes = 20)
write.table(allRes, file = "./res_topGO_BP.csv")

sessionInfo()

      

