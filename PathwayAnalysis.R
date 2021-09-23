### Take results from DEG script and run ORA analysis and specific pathway analysis here ###


library(clusterProfiler)  

ego.down <- enrichGO(gene         = substr(rownames(res[res$STATUS == "Down", ]),1,15),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.01, readable = T)  
dotplot(ego.down, showCategory=50)  


ego.up <- enrichGO(gene         = substr(rownames(res[res$STATUS == "Up", ]),1,15),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)  
dotplot(ego.up, showCategory=50)  


gene.df <- bitr(substr(rownames(res[res$STATUS != "Not sig",]),1,15), fromType = "ENSEMBL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)


kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
dotplot(kk, showCategory=20)

kk2 <- enrichKEGG(gene         = gene.df2$ENTREZID,
                  organism     = 'hsa',
                  pvalueCutoff = 0.05)
dotplot(kk2, showCategory=20)
goplot(ego.down)


mydf <- data.frame(Entrez=rownames(res), FC=res$log2FoldChange)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"

aux <- bitr(substr(mydf$Entrez,1,15), fromType = "ENSEMBL",
            toType = c("SYMBOL", "ENTREZID"),
            OrgDb = org.Hs.eg.db)

mydf <- mydf[substr(mydf$Entrez,1,15) %in% aux$ENSEMBL,]
mydf$Entrez <- aux$ENTREZID[match(substr(mydf$Entrez,1,15), aux$ENSEMBL)]

formula_res <- compareCluster(Entrez~group, data=mydf, fun="enrichKEGG")

head(as.data.frame(formula_res), 20)
dotplot(formula_res)
###



## nfkb pathway ##
nfkb <- read.delim("~/Documents/Covid/nfkb_genes.txt", skip = 2, header = F)
tnf <- read.delim("~/Documents/Covid/tnf_nfkb_genes.txt", skip = 2, header = F)
cd40 <- read.delim("~/Documents/Covid/cd40_nfkb_genes.txt", skip = 2, header = F)
immune <- read_excel("~/Documents/Covid/innatedb_curated_genes.xls") ##### from http://www.innatedb.com/
immune <- immune[immune$Species %in% '9606', ]
immune <- immune$`Gene Symbol` %>% unique
tcell <- read.delim("~/Documents/Covid/tcell_pathway_genes.txt", skip =2 , header = F)


## Annotate gene names ##
gtf = readGFF("~/reference_files/gencode.v37.annotation.gtf")
gtf = gtf[gtf$type == 'gene',]
gtf[gtf$gene_name %in% nfkb$V1,]$gene_id

nfkb.all <- c(nfkb$V1, tnf$V1, cd40$V1)
covid.nfkb <- covid.counts[rownames(covid.counts) %in% gtf[gtf$gene_name %in% nfkb.all,]$gene_id, ]
rownames(covid.nfkb) <- gtf[gtf$gene_name %in% nfkb.all,]$gene_name

immune.covid <- covid.counts[rownames(covid.counts) %in% gtf[gtf$gene_name %in% immune, ]$gene_id, ]

covid.tcell <- covid.counts[rownames(covid.counts) %in% gtf[gtf$gene_name %in% tcell$V1, ]$gene_id, ]

pheatmap(covid.nfkb[, -c(1,4,14,13)] %>% t %>% scale %>% t, , border_color = NA, color = viridis(50), show_rownames = T,
         main = "NFKB Pathway/Network genes", annotation_col = hm.anno, annotation_colors = mycolors)

pheatmap(immune.covid[, -c(1,4,14,13)] %>% t %>% scale %>% t %>% na.omit(), , border_color = NA, color = viridis(50), show_rownames = F ,
         main = "Innate immune  genes", annotation_col = hm.anno, annotation_colors = mycolors)


pheatmap(covid.tcell[, -c(1,4,14,13)] %>% t %>% scale %>% t %>% na.omit(), , border_color = NA, color = viridis(50), show_rownames = F ,
         main = "T cell signaling genes", annotation_col = hm.anno, annotation_colors = mycolors)


## Boxplot of genes 
box <- as.data.frame(t(covid.nfkb[, -c(1,4,14,13)]))
box$Group <- covid.anno$Group[-c(1,4,14,13)]

m.box <- melt(box )
ggplot(m.box, aes(Group, log1p(value), fill = Group)) + geom_boxplot() + 
  scale_fill_manual(values = brewer.pal(2, "Set2")) + ggtitle("NFKB Pathway/Network expression")


box2 <- as.data.frame(t(immune.covid[, -c(1,4,14,13)]))
box2$Group <- covid.anno$Group[-c(1,4,14,13)]

m.box2 <- melt(box2 )
ggplot(m.box2, aes(Group, log1p(value), fill = Group)) + geom_boxplot() + 
  scale_fill_manual(values = brewer.pal(2, "Set2")) + ggtitle("Innate immune genes expression")

wilcox.test(m.box2[m.box2$Group == "Cal.20c", ]$value, m.box2[m.box2$Group == "Non cal.20c", ]$value, correct = T)

box3 <- as.data.frame(t(covid.tcell[, -c(1,4,14,13)]))
box3$Group <- covid.anno$Group[-c(1,4,14,13)]

m.box3 <- melt(box3)
ggplot(m.box3, aes(Group, log1p(value), fill = Group)) + geom_boxplot() + 
  scale_fill_manual(values = brewer.pal(2, "Set2")) + ggtitle("T cell signaling genes expression")

wilcox.test(m.box3[m.box3$Group == "Cal.20c", ]$value, m.box3[m.box3$Group == "Non cal.20c", ]$value, correct = T)
