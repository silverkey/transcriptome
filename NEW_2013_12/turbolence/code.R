library("edgeR")

samples = read.table(file='samples.txt', sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F)
counts = read.table('COUNTS.txt', row.names=1, head=T, stringsAsFactors=F)
counts = counts[,samples$sample]
counts = counts[1:(nrow(counts)-5),]

dge = DGEList(counts, group=samples$condition)
cpm = cpm(dge, normalized.lib.sizes=F, log=FALSE, prior.count=0)
ridx = rowSums(cpm >= 1) >= 2
table(ridx)

dge = dge[ridx,]
cpm = cpm[ridx,]
write.table(dge$counts, file='transcripts_counts_filtered.txt',sep='\t', quote=F)
write.table(cpm, file='transcripts_cpm_filtered.txt', sep='\t', quote=F)

pdf(file='tree.pdf',paper='a4r',width=8.3,height=11.7,pointsize=8)
exprs = cpm
ecTr = dist(t(exprs), method = "euclidean")
hecTr = hclust(ecTr, method = "average")
plot(hecTr, main = "Hierarchical clustering dendrogram for cpm", xlab = "", sub = "Average linkage, Euclidean distance for cpm")
dev.off()

dge$samples$lib.size = colSums(dge$counts)
dge = calcNormFactors(dge)
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)
diff = exactTest(dge,pair=c('t','s'))
tt = topTags(diff,n=nrow(diff))

res = tt$table
sig = res[res$FDR<=0.1 & abs(res$logFC)>=log2(2),]
sig2 = merge(sig, cpm, by=0, all.x=T, sort=F)
write.table(sig2, file='significant.txt', sep='\t', quote=F, row.names=F)

