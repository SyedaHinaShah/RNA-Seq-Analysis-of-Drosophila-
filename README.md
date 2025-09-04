# RNA-Seq-Analysis-of-Drosophila-
I performed RNA seq analysis and I have countmatrix file plus sample information in which two group (Treatment and Non-treatment), I made PCA plot, Heatmap, Box plot and valcano
# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Ensure that counts_Drosophila is a data frame with rows as genes and columns as samples
str(counts_Drosophila)
# You should see something like: 
# 'data.frame': 14869 obs. of 7 variables: 
#  $ SRR031714.bam: int  157 0 14 1666 13 732 5 252 488 289 ...
#  $ SRR031716.bam: int  142 3 18 1948 19 755 4 325 584 295 ...
#  and so on...

# Ensure that SampleInfo_Drosophila is a data frame with 7 observations (samples)
str(SampleInfo_Drosophila)
# You should see something like:
# 'data.frame': 7 obs. of 3 variables: 
#  $ SampleName: chr  "SRR031714" "SRR031716" "SRR031724" ...
#  $ Group: chr  "Untreated" "Untreated" "Treated" "Treated" ...
#  $ Library: chr  "PE" "PE" "PE" "PE" ...

# Ensure the SampleName column matches the column names of counts_Drosophila
# If needed, set the column names of counts_Drosophila to match the SampleName column from SampleInfo_Drosophila
colnames(counts_Drosophila) <- SampleInfo_Drosophila$SampleName

# Create a DESeqDataSet
colData <- DataFrame(SampleInfo_Drosophila)
rownames(colData) <- colData$SampleName

# Create DESeq2 dataset from counts and sample info
dds <- DESeqDataSetFromMatrix(countData = counts_Drosophila,
                              colData = colData,
                              design = ~ Group)
### Ensure that the 'Group' variable is a factor
colData$Group <- as.factor(colData$Group)

# Now, create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_Drosophila,
                              colData = colData,
                              design = ~ Group)

# Quality Control: Perform PCA
rld <- vst(dds)  # Variance stabilizing transformation
pcaData <- plotPCA(rld, intgroup = "Group", returnData = TRUE)
ggplot(pcaData, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) + 
  theme_minimal()

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)
head(res)

# Volcano plot of results
res$padj[is.na(res$padj)] <- 1  # Fix NA values in adjusted p-value
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.6, color = "red") +
  theme_minimal() +
  xlim(c(-5, 5)) +
  ylim(c(0, 10))
##
# Volcano plot for DEGs
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(p-adj)") +
  theme(legend.position = "none")

# Heatmap for significant genes (adjusted p-value < 0.05)
vsd <- vst(dds)
select <- rownames(res)[which(res$padj < 0.05)]  # Get significant genes
mat <- assay(vsd)[select, ]
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE)

# Subset results to only include significant genes
resSig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(resSig)  # View significant genes
##boxplot
# Variance Stabilizing Transformation (VST) for better visualization
vsd <- vst(dds)

# Boxplot of log-transformed counts
boxplot(assay(vsd), notch = TRUE, col = c("blue", "red", "blue", "red", "blue", "red", "blue"),
        main = "Expression Distribution by Sample", 
        ylab = "Log2 Count", xlab = "Samples")
###
# Use variance-stabilizing transformation to normalize expression data
vsd <- vst(dds)

# Select genes that are significantly differentially expressed (padj < 0.05)
select <- rownames(res)[which(res$padj < 0.05)]

# Extract the matrix for these genes
mat <- assay(vsd)[select, ]

# Plot a heatmap
pheatmap(mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = FALSE, 
         show_colnames = TRUE,
         annotation_col = data.frame(Group = colData(dds)$Group))
###
# Check the length and contents of Group in colData
length(colData(dds)$Group)  # Should match the number of columns in 'mat'
table(colData(dds)$Group)   # Check how many unique groups there are
# Make sure 'Group' is a factor
colData(dds)$Group <- factor(colData(dds)$Group, levels = c("Untreated", "Treated"))

# Create annotation data frame
annotation_col <- data.frame(Group = colData(dds)$Group)
rownames(annotation_col) <- colnames(mat)  # Match row names with sample names

# Plot the heatmap
pheatmap(mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         show_rownames = FALSE, 
         show_colnames = TRUE,
         annotation_col = annotation_col)
         
<img width="286" height="245" alt="PCA" src="https://github.com/user-attachments/assets/822e0e74-ab93-484f-b059-3d44b63f3b3c" />
<img width="286" height="245" alt="valcano2" src="https://github.com/user-attachments/assets/0b4c41b5-dfe4-426e-8130-046754bf1104" />
<img width="286" height="245" alt="Valcano" src="https://github.com/user-attachments/assets/14c8e457-bb1c-4934-ad6e-381b5e66ea1b" />
<img width="286" height="245" alt="heatmap ok" src="https://github.com/user-attachments/assets/3c846392-6f01-4123-923b-7fb5a7e94f13" />
<img width="286" height="245" alt="boxplot" src="https://github.com/user-attachments/assets/325c3e04-519b-4475-bbaf-e90e552caa94" />




<img width="286" height="245" alt="PCA" src="https://github.com/user-attachments/assets/8a75cb1d-2ee3-4845-84bf-a183979df3c5" />

