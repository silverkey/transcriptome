# Load libraries to use
library("edgeR")

# Read the table with the counts
counts = read.table(file='../transcripts.counts.matrix', head=T, quote='', comment.char='', sep='\t')
rownames(counts) = counts[,1]
counts = counts[,-1]

# Read the table containing samples metadata
samples = read.table(file='samples.txt',sep='\t',quote='',comment.char='',head=T)

# -------------------- #
# MATING TYPE ANALYSIS #
# -------------------- #

# Create the DGE object
dge = DGEList(counts, group=samples$mating)

# Calculate the normalization factors
dge = calcNormFactors(dge)

# Filter to consider only genes with at least 1 read per million in at least 2 experiments
m = 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx = rowSums(m > 1) >= 2
table(ridx)
dge = dge[ridx,]

# Write the counts table only with the transcripts passing the filter
write.table(dge$counts,file='transcripts_counts_filtered.txt',sep='\t',quote=F)

# Estimate dispersions
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)

# Work on normalized
norm.counts = dge$pseudo.counts
trexprs = t(norm.counts)
dist = dist(trexprs)
plot(hclust(dist))

# Write the normalized counts table only with the transcripts passing the filter
write.table(dge$pseudo.counts,file='transcripts_normalized_counts_filtered.txt',sep='\t',quote=F)

# Calculate the differential expression based on mating types
diff = exactTest(dge,pair=c('+','-'))
tt = topTags(diff,n=nrow(diff))

# Prepare table with selected
res = tt$table
sig = res[res$FDR<=0.01,]

# Merge together results and counts
counts.sig = merge(counts,sig,by='row.names')

# Order based on FDR
counts.sig = counts.sig[order(counts.sig$FDR),]

# Write the table with the significantly differentially expressed transcripts
write.table(counts.sig,file='mating_differentially_expressed_transcripts.txt',row.names=F,sep='\t',quote=F)


# ------------- #
# SIZE ANALYSIS #
# ------------- #

# Create the DGE object
dge = DGEList(counts, group=samples$size)

# Calculate the normalization factors
dge = calcNormFactors(dge)

# Filter to consider only genes with at least 1 read per million in at least 2 experiments
m = 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx = rowSums(m > 1) >= 2
table(ridx)
dge = dge[ridx,]

# Write the counts table only with the transcripts passing the filter
write.table(dge$counts,file='transcripts_counts_filtered.txt',sep='\t',quote=F)

# Estimate dispersions
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)

# Calculate the differential expression based on mating types
diff = exactTest(dge,pair=c('L','S'))
tt = topTags(diff,n=nrow(diff))

# Prepare table with selected
res = tt$table
sig = res[res$FDR<=0.1,]

# Merge together results and counts
counts.sig = merge(counts,sig,by='row.names')

# Order based on FDR
counts.sig = counts.sig[order(counts.sig$FDR),]

# Write the table with the significantly differentially expressed transcripts
write.table(counts.sig,file='size_differentially_expressed_transcripts.txt',row.names=F,sep='\t',quote=F)

# Check for common trnscripts in the two significant lists
one = read.table(file='mating_differentially_expressed_transcripts.txt',sep='\t',head=T)
two = read.table(file='size_differentially_expressed_transcripts.txt',sep='\t',head=T)
sum(one$Row.names %in% two$Row.names)











pdf(file='clustering.pdf',paper='a4r',width=8.3,height=11.7,pointsize=8)
