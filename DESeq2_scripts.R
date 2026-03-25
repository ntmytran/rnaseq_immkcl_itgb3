

library(tidyverse)
library('DESeq2')
library(readxl)
library(ggrepel)
library(biomaRt)
library(dplyr)
library("RColorBrewer")
library(pheatmap)
library("tximport")
library(EnhancedVolcano)
library(ggplot2)
library(ggtext)

set.seed(265)

#####------------ import ---------------

#import sample data
samples <- read.csv("samples.csv")
samples <- samples %>% column_to_rownames("sample")


#import count data (from *.isoform.results)
#create tx2gene table using ensembldb
library("EnsDb.Hsapiens.v110")
library(ensembldb)
txdb <- transcripts(EnsDb.Hsapiens.v110, return.type = "DataFrame")
tx2gene <- as.data.frame(txdb$tx_id)
tx2gene$gene_id <- txdb$gene_id
colnames(tx2gene) <- c("transcript_id","gene_id")
#import files
files <- list.files("./count")
names(files) <- c("ITGB3KO4","ITGB3KO8","ITGB3KO9",
                  "ITGB3SPL1_1","ITGB3SPL1_2","ITGB3SPL6_1",
                  "ITGB3SPL6_2",
                  "WTKR_1","WTKR_2",
                  "WTKW_1","WTKW_2")
print(files)
setwd("/count")
txi_rsem_isoforms <- tximport(files, type = "rsem", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)
setwd("/merge")
head(txi_rsem_isoforms)


#####----------- filtering ----------------

#preparing for DESeq
samples$condition <- factor(samples$condition, levels = c("wt","mutant"))

#create dds accounting for batch effect
dds <- DESeqDataSetFromTximport(txi_rsem_isoforms, samples, ~batch+condition)

#remove non-expressed genes
head(counts(dds), n = 30)
keep <- rowSums(counts(dds) >=5) >= 3 #more than 5 reads in at least 3 samples (smallest group sizes is KO: 3)
dds <- dds[keep,]
head(assay(dds), n = 30)

#####-------- vst and pca-------------------

#transform data
vsd <- vst(dds, blind = F)
vsd$ID <- rownames(samples)
plotPCA(vsd, intgroup = "mutation", returnData = F) +
  geom_text_repel(aes(label = vsd$ID))

#remove batch effect within PCA with limma
mat <- assay(vsd)
mm <- model.matrix(~condition, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch = vsd$batch, design = mm)
assay(vsd) <- mat

vsd$mutation <- factor(vsd$mutation,
                       levels = c("wt","ko","spl"))

plotPCA(vsd, intgroup = "mutation", returnData = F) +
  geom_text_repel(aes(label = vsd$ID))



# Use ggplot to plot PCA
# Extract PCA data
pcaData <- plotPCA(vsd, intgroup = "mutation", returnData = TRUE)

# Define labels with expressions

mutation_labels <- c(
  "wt" = "WT",
  "ko" = expression(paste(italic("ITGB"), italic("3"))^{"-/-"}),
  "spl" = expression(paste(italic("ITGB"), italic("3"))^{"WT/D673_E713del"})
)

# Define fill and outline colors
mutation_fills <- c(
  "wt"  = "#D3D3D3",
  "ko"  = "#F3C7C7",
  "spl" = "#CEE8F1"
)

mutation_outlines <- c(
  "wt"  = "#434343",
  "ko"  = "#E3615F",
  "spl" = "#499CB8"
)


samples$clone <- c("ko4","ko8","ko9",
                   rep("spl1",2),
                   rep("spl6",2),
                   rep("wt",4))


clone_fills <- c("#FFDABE","#FFD2D1","#B9A8CC",
                 rep("#8FC4D7",2),
                 rep("#C4E3EE",2),
                 rep("#DCDCDC",4))

clone_outlines <- c("#FE810D","#FD8586","#714FA6",
                 rep("#447C90",2),
                 rep("#3B8AA9",2),
                 rep("#5D5D5D",4))







# Build ggplot with shape=21 (filled circle with border)
ggplot(pcaData, aes(x = PC1, y = PC2, fill = mutation, color = mutation, label = name)) +
  geom_point(shape = 21, size = 5, stroke = 1.2, alpha = 0.9) +
  scale_fill_manual(values = mutation_fills, labels = mutation_labels) +
  scale_color_manual(values = mutation_outlines, labels = mutation_labels) +
  labs(
    x = paste0("PC1: ", round(attr(pcaData, "percentVar")[1] * 100), "% variance"),
    y = paste0("PC2: ", round(attr(pcaData, "percentVar")[2] * 100), "% variance"),
    fill = "Mutation",
    color = "Mutation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 12),
    legend.position = "bottom",
  ) +
  theme(axis.text.x = element_markdown())




#####--------------- deseq2 ------------


#differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)

#plot dispersion
#plotDispEsts(dds)

#extract results
res <- results(dds, name = "condition_mutant_vs_wt")
summary(res, alpha = 0.05)

res_df <- as.data.frame(res)
res_df <- res_df %>% rownames_to_column("ids")



#####-------------- convert gene names ----------------
###Convert Ensembl IDs to gene symbols with biomaRt
#build a mart
ensembl.con <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#list attributes and filters to get gene names

res_genename <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                      filters = 'ensembl_gene_id',
                      values = res_df$ids,
                      mart = ensembl.con)
colnames(res_genename) <- c("ids","gene_name")
res_df <- full_join(res_df,res_genename, by = "ids")



#####-------------- volcano ---------------
# Volcano plot
sub_0.7_up <- subset(res_df, padj < 0.05 & log2FoldChange > 0.7)
sub_0.7_down <- subset(res_df, padj < 0.05 & log2FoldChange < -0.7)

sub_label_down <- subset(res_df, gene_name == "PHGDH"| gene_name == "PSAT1"| gene_name == "SHMT2"| gene_name == "ASNS"|gene_name == "SLC1A4"|
                           gene_name == "MTHFD2"|gene_name=="CBS"|gene_name == "SLC3A2")

nudge_values <- data.frame(gene_name = c("PHGDH",'PSAT1','SHMT2','ASNS',
                                         'SLC1A4','MTHFD2','CBS','SLC3A2'),
                           nudge_y = c(2.6, 2.5, 4, 3.6, 5, -2.2, -2, -3),
                           nudge_x = c(-3.5, -4, -0.5, 1.3, -2, 0, -2, 6.5))

sub_label_down <- sub_label_down %>%
  left_join(nudge_values, by = "gene_name")

#volcano for |log2FC| > 0.7
png("ITGB3-merged_volcano_log2FC0.7.png",unit="in", res=1200, height=3, width=5)
ggplot(res_df) +
  geom_point(aes(log2FoldChange, -log10(padj)), color = "gray50", shape = 21, size = 3) +
  geom_point(data=sub_0.7_up, aes(log2FoldChange, -log10(padj)), colour= "red3", shape = 21, size = 3) +
  geom_point(data=sub_0.7_down, aes(log2FoldChange, -log10(padj)), colour= "blue4", shape = 21, size = 3) +
  theme(text= element_text(size=12)) +
  theme(axis.text.y = element_text(size=12, face=1, colour="black")) +
  theme(axis.text.x = element_text(size=12, face=1, colour="black")) +
  theme(panel.background = element_rect(fill="white")) +
  theme(axis.line = element_line(colour="black", size=0.3)) +
  xlab("\n-log2(FoldChange)") +
  ylab("\n-log10(adjusted p-value)") +
  theme(axis.text.x = ggtext::element_markdown()) +
  theme(axis.title.y = ggtext::element_markdown()) +
  geom_text_repel(data=sub_label_down, aes(log2FoldChange, -log10(padj), label=gene_name),
                  size = 3.5,
                  nudge_y = sub_label_down$nudge_y,
                  nudge_x = sub_label_down$nudge_x,
                  box.padding = 0.5,
                  point.padding = 0.2) +
  theme(plot.margin = margin(t = 5, r = 5, b = 5,l = 5, unit = "pt"))
dev.off()





#####---------------normalized count plots-----------------------
library(edgeR)

countsmatrix <- counts(dds)
group <- samples$condition
dge <- DGEList(counts=countsmatrix,group=group)
dge <- calcNormFactors(dge, method = "TMM")
tmm_cpm <- cpm(dge)
tmm_df <- as.data.frame(tmm_cpm)
colnames(tmm_df) <- c("KO4","KO8","KO9","SPL1_1","SPL1_2","SPL6_1","SPL6_2","WTKR1","WTKR2","WTKDW1","WTKDW2")
tmm_df <- tmm_df[,c((1:7),10,11,8,9)]



plot_genes_tmm <- function(tmm_df, gene_ids, gene_names = NULL, group_vector, clone_vector, gene_grouping = NULL) {
  # Subset TMM for selected genes
  plot_df <- tmm_df[rownames(tmm_df) %in% gene_ids, ]
  
  # Convert to long format
  plot_df <- as.data.frame(t(plot_df))
  plot_df$sample <- rownames(plot_df)
  
  # Add group info
  plot_df$group <- group_vector
  
  # Add clones info
  plot_df$clone <- clone_vector
  
  # Pivot longer for multiple genes
  plot_df_long <- plot_df %>%
    pivot_longer(cols = all_of(gene_ids), names_to = "gene_id", values_to = "tmm")
  
  # If gene_names provided, map them
  if(!is.null(gene_names)) {
    gene_map <- setNames(gene_names, gene_ids)
    plot_df_long$gene <- gene_map[plot_df_long$gene_id]
  } else {
    plot_df_long$gene <- plot_df_long$gene_id
  }
  
  plot_df_long$group <- factor(plot_df_long$group, levels = c("WT","ITGB3KO","ITGB3SPL"))
  plot_df_long$clone <- factor(plot_df_long$clone,
                               levels = c("WTKDW1","WTKDW2",
                                          "WTKR1","WTKR2",
                                          "KO4","KO8","KO9",
                                          "SPL1_1","SPL1_2",
                                          "SPL6_1","SPL6_2"))
  plot_df_long$gene_label <- paste0("italic('", plot_df_long$gene, "')")
  
  # Add gene grouping if provided
  if(!is.null(gene_grouping)) {
    # Map gene grouping to the data
    grouping_map <- setNames(gene_grouping, gene_ids)
    plot_df_long$gene_group <- grouping_map[plot_df_long$gene_id]
    
    # Use facet_grid for gene groups and genes
    plot <- ggplot(plot_df_long, aes(x = group, y = log(tmm+1), color = clone, fill = clone)) +
      # Increased stroke thickness for geom_point (removed width parameter)
      geom_point(position = position_jitter(width = 0.3, height = 0),
                 size = 4,shape=21, stroke = 1.5) +
      scale_color_manual(values = mutation_outlines) +
      scale_fill_manual(values = mutation_fills) +
      ggnewscale::new_scale_color()+
      stat_summary(aes(group = group, colour = group), fun = mean, geom = "crossbar", width = 0.5, size = 0.3, alpha = 0.6) +
      scale_color_manual(values = clone_outlines) +
      facet_grid(gene_group ~ gene_label, scales = "free_y", labeller = labeller(gene_label = label_parsed, gene_group = label_value),nrow = 1) +
      scale_x_discrete(labels = c(
        expression(WT^phantom("/")),
        expression(paste(italic("ITGB"), italic("3"))^{"-/-"}),
        expression(paste(italic("ITGB"), italic("3"))^{"WT/D673_E713del"}))) +
      theme_minimal() +
      labs(x = NULL, y = "log (normalized count + 1)") +
      theme(legend.position = "none",
            legend.title = element_blank(),
            text = element_text(size = 12),
            axis.text.x = element_text(size = 12, face = 1, colour = "black"))
  } else {
    # Use only gene_label for faceting
    plot <- ggplot(plot_df_long, aes(x = group, y = log(tmm+1), color = clone, fill = clone)) +
      # Increased stroke thickness for geom_point (removed width parameter)
      geom_point(position = position_jitter(width = 0.3, height = 0),
                 size = 4,shape=21, stroke = 1.5) +
      scale_color_manual(values = mutation_outlines) +
      scale_fill_manual(values = mutation_fills) +
      ggnewscale::new_scale_color()+
      stat_summary(aes(group = group, colour = group), fun = mean, geom = "crossbar", width = 0.5, size = 0.3, alpha = 0.6) +
      scale_color_manual(values = clone_outlines) +
      facet_wrap(~ gene_label, scales = "free_y", labeller = label_parsed, nrow = 1) +
      scale_x_discrete(labels = c(
        expression(WT^phantom("/")),
        expression(paste(italic("ITGB"), italic("3"))^{"-/-"}),
        expression(paste(italic("ITGB"), italic("3"))^{"WT/D673_E713del"}))) +
      theme_minimal() +
      labs(x = NULL, y = "log (normalized count + 1)") +
      theme(legend.position = "none",
            legend.title = element_blank(),
            text = element_text(size = 12),
            axis.text.x = element_text(size = 12, face = 1, colour = "black"))
  }
  
  return(plot)
}



# Define labels for samples
mutation_labels <- c(
  "WTKDW1"  = "WT",
  "WTKDW2"  = "WT",
  "WTKR1"   = "WT",
  "WTKR2"   = "WT",
  "KO4" = expression(paste(italic("ITGB"), italic("3"))^{"-/-"} ~ "1"),
  "KO8" = expression(paste(italic("ITGB"), italic("3"))^{"-/-"} ~ "2"),
  "KO9" = expression(paste(italic("ITGB"), italic("3"))^{"-/-"} ~ "3"),
  "SPL1_1" = expression(paste(italic("ITGB"), italic("3"))^{"WT/D673_E713del"} ~ "1"),
  "SPL1_2" = expression(paste(italic("ITGB"), italic("3"))^{"WT/D673_E713del"} ~ "1"),
  "SPL6_1" = expression(paste(italic("ITGB"), italic("3"))^{"WT/D673_E713del"} ~ "2"),
  "SPL6_2" = expression(paste(italic("ITGB"), italic("3"))^{"WT/D673_E713del"} ~ "2")
)

clone_outlines <- c(
  "WT"  = "#5D5D5D",
  "ITGB3KO"  = "#FE810D",
  "ITGB3SPL" = "#C4E3EE"
)

# Define fill and outline colors
mutation_outlines <- c("WTKDW1" = "#5D5D5D",
                       "WTKDW2" = "#5D5D5D",
                       "WTKR1" = "#5D5D5D",
                       "WTKR2" = "#5D5D5D",
                       "KO4" = "#FE810D",
                       "KO8" = "#FF6666",
                       "KO9" = "#714FA6",
                       "SPL1_1" = "#447C90",
                       "SPL1_2" = "#447C90",
                       "SPL6_1" = "#3B8AA9",
                       "SPL6_2" = "#3B8AA9")

mutation_fills <- c("WTKDW1" = "#DCDCDC",
                    "WTKDW2" = "#DCDCDC",
                    "WTKR1" = "#DCDCDC",
                    "WTKR2" = "#DCDCDC",
                    "KO4" = "#FFDABE",
                    "KO8" = "#FFD2D1",
                    "KO9" = "#B9A8CC",
                    "SPL1_1" = "#8FC4D7",
                    "SPL1_2" = "#8FC4D7",
                    "SPL6_1" = "#C4E3EE",
                    "SPL6_2" = "#C4E3EE")

mutation_fills2 <- c("WTKDW1" = "#5D5D5D",
                     "WTKDW2" = "#5D5D5D",
                     "WTKR1" = "#5D5D5D",
                     "WTKR2" = "#5D5D5D",
                     "KO4" = "#FE810D",
                     "KO8" = "#FF6666",
                     "KO9" = "#714FA6",
                     "SPL1_1" = "#447C90",
                     "SPL1_2" = "#447C90",
                     "SPL6_1" = "#C4E3EE",
                     "SPL6_2" = "#C4E3EE")

library(cowplot)
group_vector <- c(rep("ITGB3KO",3), rep("ITGB3SPL",4), rep("WT",4))
clone_vector <- c("KO4","KO8","KO9","SPL1_1","SPL1_2","SPL6_1","SPL6_2","WTKDW1","WTKDW2","WTKR1","WTKR2")

# Downstream metabolism
gene_ids <- c("ENSG00000182199", "ENSG00000065911","ENSG00000160200","ENSG00000070669")
gene_names <- c("SHMT2","MTHFD2","CBS","ASNS")
downstream <- plot_genes_tmm(tmm_df, gene_ids, gene_names, group_vector,clone_vector) #remember to add nrow=1 in facet wrap
plot(downstream) 
title_plot <- ggplot() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "gray90", color = "black", linewidth = 1) +
  annotate("text", x = 0, y = 0, label = "Downstream metabolism", 
           size = 5, fontface = "bold") +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

final_plot <- plot_grid(title_plot, downstream, 
                        ncol = 1, rel_heights = c(0.1, 0.9))

png("ITGB3-merged_tmm-downstream-metabolism.png",unit="in", res=1200, height=5, width=16.8)
plot(final_plot)
dev.off()


 # Serine synthesis
gene_ids <- c("ENSG00000092621", "ENSG00000135069")
gene_names <- c("PHGDH", "PSAT1")
serinesynthesis <- plot_genes_tmm(tmm_df, gene_ids, gene_names, group_vector,clone_vector)
plot(serinesynthesis)

title_plot <- ggplot() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "gray90", color = "black", linewidth = 1) +
  annotate("text", x = 0, y = 0, label = "Serine biosynthesis", 
           size = 5, fontface = "bold") +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

final_plot <- plot_grid(title_plot, serinesynthesis, 
                        ncol = 1, rel_heights = c(0.1, 0.9))
png("ITGB3-merged_tmm-serine-synthesis.png",unit="in", res=1200, height=5, width=8.4)
plot(final_plot)
dev.off()

# Serine uptake
gene_ids <- c("ENSG00000115902")
gene_names <- c("SLC1A4")
slc <- plot_genes_tmm(tmm_df, gene_ids, gene_names, group_vector, clone_vector)
plot(slc)

title_plot <- ggplot() + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           fill = "gray90", color = "black", linewidth = 1) +
  annotate("text", x = 0, y = 0, label = "Serine uptake", 
           size = 5, fontface = "bold") +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

final_plot <- plot_grid(title_plot, slc, 
                        ncol = 1, rel_heights = c(0.1, 0.9))
png("ITGB3-merged_tmm-serine-uptake.png",unit="in", res=1200, height=5, width=4.3)
plot(final_plot)
dev.off()


# Just SLC3A2
gene_ids <- "ENSG00000168003"
gene_names <- "SLC3A2"


slc3a2 <- plot_genes_tmm(tmm_df, gene_ids, gene_names, group_vector, clone_vector)

png("ITGB3-merged_tmm-SLC3A2.png",unit="in", res=1200, height=3.34, width=4.3)
plot(slc3a2)
dev.off()

slc3a2_df <- tmm_df[rownames(tmm_df) %in% gene_ids, ]






#####--------------- TPM (PCA then plots)
ko <- read.table("/tpm/rsem.merged.gene_tpm_ko.tsv", sep = "\t",
                 header = T)
ko <- ko[,-2]
colnames(ko) <- c("gene_id","KO4","KO8","KO9","WTKR1")


spl_ms <- read.table("/tpm/rsem.merged.gene_tpm_spl-ms.tsv", sep = "\t",
                     header = T)
spl <- spl_ms[,-(2:6)]
colnames(spl) <- c("gene_id","SPL1_1","SPL1_2","SPL6_1","SPL6_2","WTKR2")

wt <- read.table("/tpm/rsem.merged.gene_tpm_wt.tsv", sep = "\t",
                 header = T)
wt <- wt[,-2]
colnames(wt) <- c("gene_id","WTKDW1","WTKDW2")

tpm <- inner_join(ko, spl, by = "gene_id")
tpm <- inner_join(tpm, wt, by = "gene_id")
tpm <- tpm %>% column_to_rownames("gene_id")
tpm <- tpm[, order(colnames(tpm))]


#log transform (log2) the TPM matrix
tpm <- as.matrix(tpm)
log_tpm <- log2(tpm+1)
#filter out non-expressed genes
keep <- rowMeans(log_tpm) >= 0
log_tpm <- log_tpm[keep, ]

# adjusting for batch effect with combat seq
library(sva)
#batch <- c(rep(1, 3), rep(2, 4), rep (3,2), rep(1,1), rep(2,1))
batch2 <- c(rep(1,7), rep (2,2), rep (1,2))

adjusted_tpm <- ComBat_seq(tpm, batch = batch2, group = NULL)
adjusted_log <- log2(adjusted_tpm+1)

keep <- rowMeans(adjusted_log) >= 0
adjusted_log <- adjusted_log[keep, ]


#PCA on the scaled tpm
pcs <- prcomp(t(adjusted_log), center = TRUE, scale. = FALSE)
percent_variance <- round((pcs$sdev^2) / sum(pcs$sdev^2) * 100, 2)
pcs_df <- as.data.frame(pcs$x)



# Build ggplot with shape=21 (filled circle with border)
png("all-mutants_4WT_PCA.png",unit="in", res=1200, height=3, width=5)
ggplot(pcs_df, aes(x = PC1, y = PC2, fill = rownames(pcs_df), color = rownames(pcs_df))) +
  geom_point(shape = 21, size = 5, stroke = 1.2, alpha = 0.9) +
  scale_fill_manual(values = mutation_fills2, labels = mutation_labels) + #same filling colors
  scale_color_manual(values = mutation_outlines, labels = mutation_labels) +
  labs(
    x = paste0("PC1: ", round(percent_variance[1]), "% variance"),
    y = paste0("PC2: ", round(percent_variance[2]), "% variance"),
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 12),
    legend.position = "none"
  )
dev.off()









