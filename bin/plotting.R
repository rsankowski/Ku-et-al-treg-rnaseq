library(DESeq2)
library(tidyverse)
library(pheatmap)
library(viridis)
library(clusterProfiler)
library(tidytext)
library(org.Mm.eg.db)

source("/home/roman/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R")

counts <- read.delim("data/Galaxy94-[Column_Join_on_data_92,_data_87,_and_others].tabular", sep="\t")

genes <- bitr(counts$X.KEY, fromType = "ENTREZID", toType = "SYMBOL",
              OrgDb = org.Mm.eg.db)
colnames(counts)[1] <- "ENTREZID"
counts$ENTREZID <- as.character(counts$ENTREZID)
counts <- inner_join(counts, genes) %>% as.data.frame()
rownames(counts) <- counts$SYMBOL

counts <- as.matrix(counts[,-c(1,10)])

#gene diffgenes
metadata <- data.frame(row.names=colnames(counts), condition=c(rep('T138',2), rep('T143',3), rep('WT',3)))
all(rownames(metadata) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

#normalized counts
normalized_counts <- counts(dds, normalized=TRUE)

#log transformation
vsd_wt <- vst(dds, blind=TRUE)

# Extract the vst matrix from the object
vsd_cor_wt <- vsd_wt %>% 
  assay() %>%
  cor()

# plot heatmap
pheatmap(vsd_cor_wt, annotation = select(metadata, condition))

pdf('plots/sample-cor-heatmap.pdf')
pheatmap(vsd_cor_wt, annotation = select(metadata, condition))
dev.off()


# Plot PCA
plotPCA(vsd_wt, intgroup="condition")

pdf('plots/pca.pdf')
plotPCA(vsd_wt, intgroup="condition")
dev.off()


## Run analysis
dds_wt <- DESeq(dds)

## Plot dispersion estimates
plotDispEsts(dds_wt)

#extract highly variable genes
library("genefilter")
topVarGenes <- head(order(-rowVars(assay(dds_wt))),1500)

genes <- rownames(counts)[topVarGenes]
write_csv(data.frame(ID=genes), "data/highly_variable_genes.csv")

#
res <- results(dds_wt,
               contrast = c("condition", "WT", "T138"),
               alpha = 0.1,lfcThreshold = 0.32) #
#MA plot
plotMA(res, ylim=c(-8,8))

#shrink fc
res <- lfcShrink(dds_wt,
                 contrast=c("condition", "WT", "T138"),
                 res=res)

#summarize res
summary(res)

#genes <- read_csv('data/gencode.vM18.annotation.tx2gene.csv')
#colnames(genes) <- c("Ensembl_Trancr_ID","Ensembl_Gene_ID")
#genes$Ensembl_Gene_ID <- gsub("\\..*", '', genes$Ensembl_Gene_ID) 
#genes2 <- read.delim('data/to_deliver/gene_annotation_info.txt', stringsAsFactors = F)
#genes2$Ensembl_Gene_ID <- gsub("\\..*", '', genes2$Ensembl_Gene_ID) 

#genes <- genes %>% left_join(genes2)

res_all <- data.frame(res) %>%
  rownames_to_column(var = "Ensembl_Gene_ID") 
#res_all <- res_all %>% left_join(x = res_all, y = genes[,c("Ensembl_Trancr_ID","Ensembl_Gene_ID", "Gene_name")], by = "Ensembl_Trancr_ID")

#View(res_all)

#extract significant genes 
res_sig <- subset(res_all, pvalue < 0.05)
res_sig <- res_sig %>%
  arrange(padj)

#plot sig genes
# Subset normalized counts to significant genes
sig_norm_counts <- normalized_counts[res_sig$Ensembl_Gene_ID, ]
sig_norm_counts <- sig_norm_counts[grepl("\\:", res_sig$Ensembl_Gene_ID),]
#rownames(sig_norm_counts) <- res_sig$Ensembl_Gene_ID
# Define color palette
#heat_colors <- brewer.pal(6, "YlOrRd")
# Run pheatmap
annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6))),
  AdditionalAnnotation = c(rep("random1", 10), rep("random2", 10))
)

pheat <- pheatmap(sig_norm_counts,
                  color = viridis(100),
                  cluster_rows = T,
                  show_rownames = T,
                  annotation = dplyr::select(metadata, condition),
                  scale = "row")

pdf('plots/final-erv-diff-gene-heatmap.pdf', height = 15,width=8)
pheat
dev.off()

res_sig[res_sig$log2FoldChange>0,]

write_csv(res_sig, 'data/final-res_sig-te.csv')
write_csv(data.frame("gene"=rownames(sig_norm_counts), sig_norm_counts), 'data/final-norm_sig_counts-te.csv')
