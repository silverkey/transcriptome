# Load edgeR library
library("edgeR")

# Load the sample file
samples = read.table(file='sample_file', sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F)

# Load the raw counts file
counts = read.table('transcripts_counts_filtered.txt', row.names=1, head=T, stringsAsFactors=F)

# Put the counts file column in the same order as the sample file
counts = counts[,samples$name]

# Create the edgeR object
dge = DGEList(counts, group=samples$condition)

# Normalize
dge = calcNormFactors(dge)
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)

# Clustering of the experiments
exprs = dge$pseudo.counts
ecTr = dist(t(exprs), method = "euclidean")
hecTr = hclust(ecTr, method = "average")
plot(hecTr, main = "Hierarchical clustering dendrogram for cpm", xlab = "", sub = "Average linkage, Euclidean distance for cpm") 
hecTr = hclust(ecTr, method = "complete")
plot(hecTr, main = "Hierarchical clustering dendrogram for cpm", xlab = "", sub = "Average linkage, Euclidean distance for cpm")

# Espressione differenziale M vs F
# Execute the statistical test
diff.mf = exactTest(dge,pair=c('M','F'))
# Extract the results
tt.mf = topTags(diff.mf,n=nrow(diff.mf))
# Put the results into a table
res.mf = tt.mf$table
# Select the significant ones according to the cutoffs
sig.mf = res.mf[res.mf$FDR<=0.1 & abs(res.mf$logFC)>=log2(2),]
# Add the normalized counts to the statistical results
sig.mf = merge(sig.mf, exprs, by=0, all.x=T, sort=F)
# Write a table with all the info for the significantly differentially expressed genes
write.table(sig.mf, file='M_F_significant.txt', sep='\t', quote=F, row.names=F)

# Espressione differenziale M vs B
# Execute the statistical test
diff.mb = exactTest(dge,pair=c('M','B'))
# Extract the results
tt.mb = topTags(diff.mb,n=nrow(diff.mb))
# Put the results into a table
res.mb = tt.mb$table
# Select the significant ones according to the cutoffs
sig.mb = res.mb[res.mb$FDR<=0.1 & abs(res.mb$logFC)>=log2(2),]
# Add the normalized counts to the statistical results
sig.mb = merge(sig.mb, exprs, by=0, all.x=T, sort=F)
# Write a table with all the info for the significantly differentially expressed genes
write.table(sig.mb, file='M_B_significant.txt', sep='\t', quote=F, row.names=F)

# Espressione differenziale F vs B
# Execute the statistical test
diff.fb = exactTest(dge,pair=c('F','B'))
# Extract the results
tt.fb = topTags(diff.fb,n=nrow(diff.fb))
# Put the results into a table
res.fb = tt.fb$table
# Select the significant ones according to the cutoffs
sig.fb = res.fb[res.fb$FDR<=0.1 & abs(res.fb$logFC)>=log2(2),]
# Add the normalized counts to the statistical results
sig.fb = merge(sig.fb, exprs, by=0, all.x=T, sort=F)
# Write a table with all the info for the significantly differentially expressed genes
write.table(sig.fb, file='F_B_significant.txt', sep='\t', quote=F, row.names=F)
