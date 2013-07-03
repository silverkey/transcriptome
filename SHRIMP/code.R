library("edgeR")

pair = c('test','control')
fdr.cutoff = 0.05
logfc.cutoff = 1

counts = read.table(file='genes.counts.matrix', head=T, quote='', comment.char='', sep='\t')
rownames(counts) = counts[,1]
counts = counts[,2:ncol(counts)]

samples = read.table(file='sample_file',sep='\t',quote='',comment.char='',head=T)

cn = colnames(counts)
cn = sub( ".genes.results", "", cn)
colnames(counts) = cn

counts = counts[,samples$sample.name]

dge = DGEList(counts, group=samples$condition)
dge = calcNormFactors(dge)

m = 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx = rowSums(m > 1) >= 2
table(ridx)
dge = dge[ridx,]
write.table(dge$counts,file='genes_counts_filtered.txt',sep='\t',quote=F)

dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)

diff = exactTest(dge,pair=pair)
tt = topTags(diff,n=nrow(diff))
res = tt$table

sig = subset(res, FDR <= fdr.cutoff & abs(logFC) >= logfc.cutoff)
nrow(sig)

counts.sig = merge(sig[,c(1,4)],dge$pseudo.counts,by='row.names',sort=F)

write.table(counts.sig,file='differentially_expressed_genes.txt',row.names=F,sep='\t',quote=F)

sessionInfo()

#trexprs = t(dge$pseudo.counts)
#dist = dist(trexprs)
#plot(hclust(dist))

#trexprs = t(dge$pseudo.counts)
#dist = dist(trexprs)
#plot(hclust(dist))
