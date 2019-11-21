library(DESeq2)
library(tidyverse)
library(pheatmap)
library(viridis)
library(clusterProfiler)
library(tidytext)

source("/home/roman/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R")

counts <- read.delim("/home/roman/Documents/bulk-rna-seq/stz-wt/data/counts.txt", skip = 1)
#exclude stz
counts <- counts[,-c(18,10,25,19,21,11,24,8)]
counts <- counts[,!grepl("S353", colnames(counts))]

rownames(counts) <- counts$Geneid
counts <- counts[,c(12, 8,7,14,13,10,11,9)]
colnames(counts )<- c("Veh1", "Veh2", "Veh3", "Veh4","STZ1", "STZ2", "STZ3", "STZ4")
counts <- as.matrix(counts[,-8])

write.csv(counts, "data/4veh-vs-3stz-mavsko-genomic-counts.csv")

#gene diffgenes
metadata <- data.frame(row.names=colnames(counts), condition=c(rep('Veh',4), rep('STZ',3)))
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

#
res <- results(dds_wt,
               contrast = c("condition", "STZ",
                            "Veh"),
               alpha = 0.05) #,lfcThreshold = 0.32
#MA plot
plotMA(res, ylim=c(-8,8))

#shrink fc
res <- lfcShrink(dds_wt,
                 contrast=c("condition", "STZ", "Veh"),
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
res_sig <- subset(res_all, padj < 0.05)
res_sig <- res_sig %>%
  arrange(log2FoldChange)

#plot sig genes
# Subset normalized counts to significant genes
sig_norm_counts <- normalized_counts[res_sig$Ensembl_Gene_ID, ]
#rownames(sig_norm_counts) <- res_sig$Ensembl_Gene_ID
# Define color palette
#heat_colors <- brewer.pal(6, "YlOrRd")
# Run pheatmap
pheat <- pheatmap(sig_norm_counts,
                  color = viridis(100, option = 4),
                  cluster_rows = T,
                  show_rownames = T,
                  annotation = dplyr::select(metadata, condition),
                  scale = "row")

pdf('plots/non-erv-diff-gene-heatmap.pdf', height = 15,width=8)
pheat
dev.off()

res_sig[res_sig$log2FoldChange>0,]

write_csv(res_sig, 'data/final-res_sig.csv')
write_csv(data.frame("gene"=rownames(sig_norm_counts), sig_norm_counts), 'data/final-norm_sig_counts.csv')

#Coposition of the dysregulation of ERVs
genes <- read_csv('data/final-res_sig.csv') %>%
  filter(grepl("\\:",Ensembl_Gene_ID)) %>%
  tidyr::separate(Ensembl_Gene_ID, c("geneid", "family", "class"), sep=':') %>%
  mutate(condition = ifelse(log2FoldChange > 0, "STZ", "Veh")) 
#mutate(condition = factor(ifelse(log2FoldChange > 0, "STZ", "Veh"), levels = c("Veh", "STZ"))) 


ggplot(genes, aes(x=condition, fill=family)) +
  geom_bar(position='fill') +
  theme_minimal() +
  scale_fill_manual(values=colors_many)

genes %>%
  group_by(condition, family) %>%
  summarise(freq = length(family)) %>%
  ggplot(aes(x=1, y=freq,fill=family)) +
    geom_bar(position = 'fill', stat = 'identity', width = 1, color='black', lwd=0.1) +
    coord_polar(theta='y') +
    theme_void() +
    scale_fill_manual(values=colors_many) +
  facet_wrap(~condition)

#bar plots
genes %>%
  group_by(condition, family) %>%
  summarise(freq = length(family)) %>%
  ungroup() %>%
  mutate(condition = factor(condition, levels = c("Veh", "STZ")),
                            family = reorder_within(family, freq, condition)) %>%
  ggplot(aes(x=family, y=freq,fill=condition)) +
  geom_col(show.legend = F) +
  theme_minimal() +
  scale_fill_manual(values=colors_many) +
  facet_wrap(~condition, scales = "free") +
  coord_flip() +
  scale_x_reordered()

bar_plot <- genes %>%
  group_by(condition, class) %>%
  summarise(freq = length(class)) %>%
  ungroup() %>%
  mutate(condition = factor(condition, levels = c("Veh", "STZ")),
         class = reorder_within(class, freq, condition)) %>%
  ggplot(aes(x=class, y=freq,fill=condition)) +
  geom_col(show.legend = F, color = 'black', lwd= 0.25) +
  theme_minimal() +
  facet_wrap(~condition, scales = "free") +
  coord_flip() +
  scale_x_reordered() +
  scale_fill_manual(values=colors_many) 
  


genes <- read_csv2("data/TEtools/mm10_rmsk_TE.gtf", col_names = F)[,c(2:4)]
colnames(genes) <- c("transcript_id", "family_id", "class_id")
genes1[,1] <- gsub('transcript_id \""', "", genes)
library(refGenome)
ens <- ensemblGenome()
basedir(ens) <- system.file("extdata", package="refGenome")
read.gtf(ens, "data/TEtools/mm10_rmsk_TE.gtf")
tableSeqids(ens)
genes <- getGenePositions(ens)           
