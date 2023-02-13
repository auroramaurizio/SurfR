

# BiocManager::install("SPsimSeq")
library(SPsimSeq)
library(SurfR)
?DGE
data("zhang.data.sub")
zhang.counts <- zhang.data.sub$counts
MYCN.status  <- zhang.data.sub$MYCN.status
sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts,
                          group = MYCN.status, n.genes = 1000, batch.config = 1,
                          group.config = c(0.5, 0.5), tot.samples = ncol(zhang.counts),
                          pDE = 0.2, lfc.thrld = 0.5, result.format = "list", return.details = TRUE)

sim.data.bulk1 <- sim.data.bulk$sim.data.list[[1]]

count = sim.data.bulk1$counts
row.names(count) <- row.names(zhang.counts)

meta = sim.data.bulk1$colData
meta$Group = as.factor(meta$Group)
head(count)
length(unique(row.names(count)))
table(meta$Group)
head(meta)

?Annotate_SPID

df <- DGE(expression = count,
          metadata = meta,
          Nreplica = 50,
          design = "~Group",
          condition = "Group",
          alpha = 0.05,
          TEST = "1", CTRL =  "0")
head(df)


fdr_GeneID = df[df$padj < 0.05,"GeneID"]
fdrUP_GeneID = df[df$padj < 0.05 & df$log2FoldChange > 0,"GeneID"]
fdrDW_GeneID = df[df$padj < 0.05 & df$log2FoldChange < 0,"GeneID"]

SP = Gene2SProtein(genes = fdr_GeneID, input_type ="gene_name")

SPup = Gene2SProtein(genes = fdrUP_GeneID, input_type ="gene_name")
SPdw = Gene2SProtein(genes = fdrDW_GeneID, input_type ="gene_name")


dfhead(sim.data.bulk1$counts[, seq_len(5)])


# -------- GEO ---------
#GSE121810

library(SurfR)
library(stringr)

mGSE121810 = GEOmetadata(GSE = "GSE121810")
fx = function(x) {str_split(x, " ")[[1]][1]}
mGSE121810$condition = sapply(mGSE121810$therapy, fx)
mGSE121810$condition = as.factor(mGSE121810$condition)
table(mGSE121810$condition)
cGSE121810 = DownloadArchS4(mGSE121810$GSM, species = "human", print_tsv = FALSE, filename = NULL)
#cGSE121810_2 = as.data.frame(cGSE121810)
df <- DGE(expression = cGSE121810,
          metadata = mGSE121810,
          Nreplica = 14,
          design = "~condition",
          condition = "condition",
          alpha = 0.05,
          TEST = "neoadjuvant", CTRL =  "adjuvant",
          output_tsv = F)

fdr_GeneID = df[df$padj < 0.05,"GeneID"]
fdrUP_GeneID = df[df$padj < 0.05 & df$log2FoldChange > 0,"GeneID"]
fdrDW_GeneID = df[df$padj < 0.05 & df$log2FoldChange < 0,"GeneID"]

SP = Gene2SProtein(genes = fdr_GeneID, input_type ="gene_name")
SPup = Gene2SProtein(genes = fdrUP_GeneID, input_type ="gene_name")
SPdw = Gene2SProtein(genes = fdrDW_GeneID, input_type ="gene_name")

SPdw
df[SPup$UniProt.gene,]
as.data.frame(df[SPdw$UniProt.gene,])

# --------- TCGA ------------
?TCGA_download

TCGA_UVM = TCGA_download(project = "TCGA-UVM")

TCGA_UVM = TCGA_download(project = "TCGA-UVM")
cTCGA_UVM <- TCGA_UVM[[1]]
mTCGA_UVM <- TCGA_UVM[[2]]




TCGA.PRAD = TCGA_download(project = "TCGA-PRAD")
cTCGA.PRAD <- TCGA.PRAD[[1]]
mTCGA.PRAD <- TCGA.PRAD[[2]]
table(mTCGA.PRAD$shortLetterCode)

TCGA.LUAD = TCGA_download(project = "TCGA-LUAD")
cTCGA.LUAD <- TCGA.LUAD[[1]]
mTCGA.LUAD <- TCGA.LUAD[[2]]
table(mTCGA.LUAD$shortLetterCode)

TCGA.LUSC = TCGA_download(project = "TCGA-LUSC")
cTCGA.LUSC <- TCGA.LUSC[[1]]
mTCGA.LUSC <- TCGA.LUSC[[2]]
table(mTCGA.LUSC$shortLetterCode)


TCGA.THYM = TCGA_download(project = "TCGA-THYM")
cTCGA.THYM <- TCGA.THYM[[1]]
mTCGA.THYM <- TCGA.THYM[[2]]
table(mTCGA.THYM$shortLetterCode)

mTCGA.THYM$shortLetterCode = as.factor(mTCGA.THYM$shortLetterCode)

df <- DGE(expression = cTCGA.THYM,
          metadata = mTCGA.THYM,
          Nreplica = 2,
          design = "~shortLetterCode",
          condition = "shortLetterCode",
          alpha = 0.05,
          TEST = "TP", CTRL =  "NT",
          output_tsv = F)

head(df)

fdr_GeneID = df[df$padj < 0.05,"GeneID"]
SP = Gene2SProtein(genes = fdr_GeneID, input_type ="gene_name")

fdrUP_GeneID = df[df$padj < 0.05 & df$log2FoldChange > 0,"GeneID"]
SPup = Gene2SProtein(genes = fdrUP_GeneID, input_type ="gene_name")

fdrDW_GeneID = df[df$padj < 0.05 & df$log2FoldChange < 0,"GeneID"]
SPdw = Gene2SProtein(genes = fdrDW_GeneID, input_type ="gene_name")

library(SurfR)
mGSE177522 = GEOmetadata("GSE177522")
head(mGSE177522 )
cGSE177522 = DownloadArchS4(mGSE177522$GSM, species = "human")
cGSE177522



# --------- Breast cancer meta-analalysis ---------------
library(SurfR)

TCGA.BRCA = TCGA_download(project = "TCGA-BRCA")
cTCGA.BRCA <- TCGA.BRCA[[1]]
mTCGA.BRCA <- TCGA.BRCA[[2]]

table(mTCGA.BRCA$shortLetterCode)

mGSE58135 = GEOmetadata("GSE58135")
mGSE58135 = mGSE58135[mGSE58135$tissue != "Breast Cancer Cell Line",]
mGSE58135$condition = "NT"
mGSE58135$condition[mGSE58135$tissue %in% c("ER+ Breast Cancer Primary Tumor", "Triple Negative Breast Cancer Primary Tumor")] <- "TP"
cGSE58135 <- DownloadArchS4(mGSE58135$GSM, species = "human")
cGSE58135_bk = cGSE58135
table(mGSE58135$condition)

df.TCGA <- DGE(expression = cTCGA.BRCA,
               metadata = mTCGA.BRCA,
               Nreplica = 2,
               design = "~shortLetterCode",
               condition = "shortLetterCode",
               alpha = 0.05,
               TEST = "TP", CTRL =  "NT",
               output_tsv = F)

head(df.TCGA)

# GSE58135 DGE
dim(cGSE58135)
dim(mGSE58135)
cGSE58135 = cGSE58135[,row.names(mGSE58135)]
table(mGSE58135$condition)

mGSE58135$condition = as.factor(mGSE58135$condition)
df.GSE58135 <- DGE(expression = cGSE58135,
                   metadata = mGSE58135,
                   Nreplica = 56,
                   design = "~condition",
                   condition = "condition",
                   alpha = 0.05,
                   TEST = "TP", CTRL =  "NT",
                   output_tsv = F)

mTCGA.BRCA$shortLetterCode = as.factor(mTCGA.BRCA$shortLetterCode)
df.TCGA <- DGE(expression = cTCGA.BRCA,
               metadata = mTCGA.BRCA,
               Nreplica = 113,
               design = "~shortLetterCode",
               condition = "shortLetterCode",
               alpha = 0.05,
               TEST = "TP", CTRL =  "NT",
               output_tsv = F)
head(df.TCGA)

libr
L_fishercomb = metaRNAseq(ind_deg = list(TCGA.BRCA =  df.TCGA, GEO.GSE58135 = df.GSE58135),
                          test_statistic = "fishercomb",
                          BHth = 0.05,
                          adjpval.t = 0.05)

L_invnorm = metaRNAseq(ind_deg = list(TCGA.BRCA =  df.TCGA, GEO.GSE58135 = df.GSE58135),
                       test_statistic = "invnorm",
                       BHth = 0.05,
                       adjpval.t = 0.05,
                       nrep = c(102, 56))

?combine_fisher_invnorm

metacomb <- combine_fisher_invnorm(ind_deg = list(TCGA.BRCA =  df.TCGA, GEO.GSE58135 = df.GSE58135),
                                   invnorm = L_invnorm,
                                   fishercomb = L_fishercomb,
                                   adjpval = 0.05)

metacomb_GeneID = metacomb[metacomb$signFC != 0, "GeneID"]
SP = Gene2SProtein(genes = metacomb_GeneID, input_type ="gene_name")

metacombUP_GeneID = metacomb[metacomb$signFC == 1, "GeneID"]
SPup = Gene2SProtein(genes = metacombUP_GeneID, input_type ="gene_name")

metacombDW_GeneID = metacomb[metacomb$signFC == -1, "GeneID"]
SPdw = Gene2SProtein(genes = metacombDW_GeneID, input_type ="gene_name")


## plotting and enrichment

?Enrichment

library(SPsimSeq)

data("zhang.data.sub")
zhang.counts <- zhang.data.sub$counts
MYCN.status  <- zhang.data.sub$MYCN.status
# Simulation of bulk RNA data with
sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = zhang.counts,
                          group = MYCN.status, n.genes = 1000, batch.config = 1,
                          group.config = c(0.5, 0.5), tot.samples = ncol(zhang.counts),
                          pDE = 0.2, lfc.thrld = 0.5, result.format = "list", return.details = TRUE)

sim.data.bulk1 <- sim.data.bulk$sim.data.list[[1]]

countMatrix = sim.data.bulk1$counts; row.names(countMatrix) <- row.names(zhang.counts)
metadata = sim.data.bulk1$colData; metadata$Group = as.factor(metadata$Group)

df <- DGE(expression = countMatrix,
          metadata = metadata,
          Nreplica = 50,
          design = "~Group",
          condition = "Group",
          alpha = 0.05,
          TEST = "1", CTRL =  "0",
          output_tsv = FALSE)

fdr_GeneID = df[df$padj < 0.05,"GeneID"]
SP = Gene2SProtein(genes = fdr_GeneID, input_type ="gene_name")

# upregulated fdr
fdrUP_GeneID = df[df$padj < 0.05 & df$log2FoldChange > 0,"GeneID"]
SPup = Gene2SProtein(genes = fdrUP_GeneID, input_type ="gene_name")

# dowregulated fdr
fdrDW_GeneID = df[df$padj < 0.05 & df$log2FoldChange < 0,"GeneID"]
SPdw = Gene2SProtein(genes = fdrDW_GeneID, input_type ="gene_name")
SPdw
Splot(SPdw,main = "prova")

library(enrichR)
dfList = list(exp1 = df)
test = SurfR::Enrichment(dfList)

Enrichment_barplot(test$exp1,
                   enrich.databases = c("GO_Biological_Process_2021","GO_Cellular_Component_2021"),
                   p_adj = 0.5,
                   num_term = 5,
                   cond = "UP")


DGE = df
annotated = Annotate_SPID(DGE, c("GO_Biological_Process_2021"))

SurfR::plotPCA(matrix = cpm(countMatrix), metadata = metadata, dims = c(1,2),
               shape.by = "Group", color.by = "Group", nTOP = 500)
