## Load count matrix and metadata file ##
covid.counts <- fc$counts
covid.anno <- read_xlsx("~/Documents/Covid/RNA_FullDepth/Covid_SampleSheet.xlsx")
covid.anno <- as.data.frame(covid.anno)

colnames(covid.counts) <- gsub("_R1_trimmed.fastq.gzAligned.sortedByCoord.out.bam", "", colnames(covid.counts))
#colnames(covid.counts) <- gsub("_S*", "", colnames(covid.counts))


## match the colnames of the gene expression matrix to the annotation matrix
well.id <- str_extract(colnames(covid.counts), "_(.*)_")
well.id <- gsub("_", "", well.id)
covid.anno <- covid.anno[match(well.id, covid.anno$ID),] ## should be the same order now
##


# FIlter gene expression low exp genes
filt.covid = rowSums(covid.counts) >= 10
f.covid.counts = covid.counts[filt.covid, ]

# PCA plot
pca <- prcomp(t(f.covid.counts[, -c(1,14,4,13)])) ## removing outlier samples for downstream analysis

aux <- as.data.frame(pca$x)

# PCA Plot
ggplot(aux, aes(PC1, PC2, label = rownames(aux), color = covid.anno$Group[-c(1,14,4,13)])) + geom_point(size =3) +
  geom_label_repel(aes(label = rownames(pca$x), vjust= "inward", hjust ="inward")) +
  ylab(paste0("PC2: ", prettyNum(summary(pca)$importance[2,2]*100), " %")) +
  xlab(paste0("PC1: ", prettyNum(summary(pca)$importance[2,1]*100), " %")) +
  theme_bw() + labs(color = "Sample type") + scale_color_manual(values = brewer.pal(3, "Set2"))
#



## DEG Prep and start here ##

library(DESeq2)

coldata <- data.frame(row.names = colnames(f.covid.counts)[-c(1,14,4,13)], condition = covid.anno$Group[-c(1,14,4,13)])
condition <- data.frame(condition = coldata$condition, row.names = rownames(coldata))
dds <- DESeqDataSetFromMatrix(countData = f.covid.counts[,-c(1,14,4,13)],
                              colData = condition,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Cal.20c", "Non.cal.20c"))


res <- as.data.frame(na.omit(res))
res$STATUS <- "Not sig"
res[res$log2FoldChange <= 0 & res$padj <= 0.05, ]$STATUS <- "Down"
res[res$log2FoldChange >= 0 & res$padj <= 0.05, ]$STATUS <- "Up"
# ADd gene name to result table
aux <- gtf[gtf$gene_id %in% rownames(res), ]
res$GeneSymbol <- aux[match(rownames(res), aux$gene_id), ]$gene_name

## to plot top DEGs - volcano ##
up <- res[res$STATUS == "Up", ]
up <- arrange(up, -log2FoldChange, padj)
#up <- arrange(up,  padj)
up$label <- NA
up$label[1:10] <- up$GeneSymbol[1:10]

dn <- res[res$STATUS == "Down", ]
dn <- arrange(dn, log2FoldChange, padj)
#dn <- arrange(dn, padj)
dn$label <- NA
dn$label[1:10] <- dn$GeneSymbol[1:10]
##

# transfer label to the plotting obj
res$label <- NA
res$label[which(res$GeneSymbol %in% up$label[1:10])] <- up$label[1:10]
res$label[which(res$GeneSymbol %in% dn$label[1:10])] <- dn$label[1:10]


ggplot(res, aes(log2FoldChange, -log10(padj), color = STATUS)) + geom_point() +
  theme_bw() + scale_color_manual(values = c("darkgreen", "grey", "darkred")) +
  #geom_label_repel(aes(label = res$label), vjust = "inward", hjust = "inward") +
  xlab("log2 Fold Change") + ylab("-log10 adjusted Pvalue") + ggtitle("NonCal20C vs Cal20C ")

hm <- f.covid.counts[rownames(f.covid.counts) %in% rownames(res[res$STATUS != "Not sig", ]), -c(1,14,4,13)]  
#hm <- hm[, colnames(hm) %in% covid.anno$sample[s.anno$group %in% c("Cal20C", "NonCal20C")]]
#hm <- hm[, 2:7]
hm.anno <- data.frame("Group" = covid.anno$Group[-c(1,14,4,13)], row.names = colnames(f.covid.counts)[-c(1,14,4,13)])
# creat colours for each group
cols <- colorRampPalette(brewer.pal(2, "Set2"))
mycolors <- cols(length(unique(covid.anno$Group[-c(1,4,14,13)])))
names(mycolors) <- unique(covid.anno$Group[-c(1,4,14,13)])
mycolors <- list(Group = mycolors)
mycolors$Group[2] <- "#FC8E62"

pheatmap(hm %>% t %>% scale %>% t %>% na.omit, border_color = NA, color = cividis(50), show_rownames = F,
         main = "NonCal20C vs Cal20C - DEGs (n=109)", annotation_col = hm.anno, annotation_colors = mycolors)  

