library("edgeR")

samples = read.table(file='samples.txt', sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F)
counts = read.table('COUNTS.txt', row.names=1, head=T, stringsAsFactors=F)
counts = counts[,samples$sample]
anno = read.table(file='pse_mul_oct_13_uniref_filt_ann_out.txt', sep='\t', head=T, comment.char='', quote='', stringsAsFactors=F)
map = read.table(file='mapping_locations.txt', sep='\t', head=T, comment.char='', quote='', stringsAsFactors=F)

dge = DGEList(counts, group=samples$mating)
cpm = cpm(dge, normalized.lib.sizes=F, log=FALSE, prior.count=0)
ridx = rowSums(cpm >= 0.5) >= 2
table(ridx)

dge2 = dge[ridx,]
cpm2 = cpm[ridx,]
write.table(dge2$counts, file='transcripts_counts_filtered.txt',sep='\t', quote=F)
write.table(cpm2, file='transcripts_cpm_filtered.txt', sep='\t', quote=F)

dge2$samples$lib.size = colSums(dge2$counts)
dge2 = calcNormFactors(dge2)
dge2 = estimateCommonDisp(dge2)
dge2 = estimateTagwiseDisp(dge2)
diff = exactTest(dge2,pair=c('+','-'))
tt = topTags(diff,n=nrow(diff))

res = tt$table
sig = res[res$FDR<=0.1 & abs(res$logFC)>=log2(2),]
sig2 = merge(sig, cpm, by=0, all.x=T, sort=F)
write.table(sig2, file='significant.txt', sep='\t', quote=F, row.names=F)

res2 = merge(res, cpm, by=0, all.x=T, sort=F)
write.table(res2, file='all_counts_and_statistics.txt', sep='\t', quote=F, row.names=F)

res3 = merge(res2, anno, by.x='Row.names', by.y='CompName', sort=F)
write.table(res3, file='all_statistics_annotated.txt', sep='\t', quote=F, row.names=F)

sig3 = merge(sig2, anno, by.x='Row.names', by.y='CompName', sort=F)
write.table(sig3, file='significant_statistics_annotated.txt', sep='\t', quote=F, row.names=F)

sig4 = merge(sig3, map, by.x='Row.names', by.y='transcript', sort=F)
write.table(sig4, file='significant_statistics_annotated_mapped.txt', sep='\t', quote=F, row.names=F)
