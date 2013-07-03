# Load libraries to use
library("edgeR")

# Read the table with the counts
counts = read.table(file='transcripts.counts.matrix', head=T, quote='', comment.char='', sep='\t')
rownames(counts) = counts[,1]
counts = counts[,-1]

# Read the table containing samples metadata
samples = read.table(file='samples.txt',sep='\t',quote='',comment.char='',head=T)
samples$largepVSsmallp = c(2,0,0,2,0,1)

# Create the DGE object
dge = DGEList(counts, group=samples$largepVSsmallp)

# Calculate the normalization factors
dge = calcNormFactors(dge)

# Filter to consider only genes with at least 1 read per million in at least 2 experiments
m = 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx = rowSums(m > 1) >= 2
table(ridx)
dge = dge[ridx,]

# Estimate dispersions
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)

# Calculate the differential expression based on mating types
diff = exactTest(dge,pair=c(1,2))
tt = topTags(diff,n=nrow(diff))

# Prepare table with selected
res = tt$table
sig = res[res$FDR<=0.01,]

# Merge together results and counts
counts.sig = merge(counts,sig,by='row.names')

# Order based on FDR
counts.sig = counts.sig[order(counts.sig$FDR),]

# Write the table with the significantly differentially expressed transcripts
write.table(counts.sig,file='large_plus_vs_small_plus_differentially_expressed_transcripts.txt',row.names=F,sep='\t',quote=F)



