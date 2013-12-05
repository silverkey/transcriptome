library("edgeR")

samples = read.table(file='samples.txt', sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F)
counts = read.table('COUNTS.txt', row.names=1, head=T, stringsAsFactors=F)
counts = counts[,samples$sample]
anno = read.table(file='sem_rob_oct_13_uniref_filt_ann_out.txt', sep='\t', head=T, comment.char='', quote='', stringsAsFactors=F)

dge = DGEList(counts, group=samples$mating_type)
cpm = cpm(dge,normalized.lib.sizes=F,log=FALSE,prior.count=0)
ridx = rowSums(cpm >= 0.5) >= 2
table(ridx)

dge2 = dge[ridx,]
cpm2 = cpm[ridx,]

write.table(dge2$counts, file='transcripts_counts_filtered.txt', sep='\t', quote=F)
write.table(cpm2, file='transcripts_cpm_filtered.txt', sep='\t', quote=F)

# DONE UNTIL HERE JUST TO FILTER THE COUNTS TABLE!!!!!!!

dge2$samples$lib.size = colSums(dge2$counts)
dge2 = calcNormFactors(dge2)
dge2 = estimateCommonDisp(dge2)
dge2 = estimateTagwiseDisp(dge2)

diff = exactTest(dge2, pair=c('+','-'))
tt = topTags(diff, n=nrow(diff))

res = tt$table
sig = res[res$FDR<=0.01 & abs(res$logFC)>=log2(5),]

sig2 = merge(sig, cpm, by=0, all.x=T, sort=F)
write.table(sig2, file='significant.txt', sep='\t', quote=F, row.names=F)

res2 = merge(res, cpm, by=0, all.x=T, sort=F)
write.table(res2, file='all_counts_and_statistics.txt', sep='\t', quote=F, row.names=F)

res3 = merge(res2, anno, by.x='Row.names', by.y='CompName', sort=F)
write.table(res3, file='all_statistics_annotated.txt', sep='\t', quote=F, row.names=F)

sig3 = merge(sig2, anno, by.x='Row.names', by.y='CompName', sort=F)
write.table(sig3, file='significant_statistics_annotated.txt', sep='\t', quote=F, row.names=F)

gtab = read.table(file='transcripts_on_genome_part', sep='\t', stringsAsFactors=F)
gtab = na.omit(res[gtab$V1,])
gtab = merge(gtab, cpm, by=0, all.x=T, sort=F)
gtab = merge(gtab, anno, by.x='Row.names', by.y='CompName', sort=F)
write.table(gtab, file='transcripts_scaf1897_stat_anno.txt', sep='\t', quote=F, row.names=F)

#hsigp = res[res$FDR<=0.001 & res$logFC>=3,]
#hsigp = merge(hsigp,cpm,by=0,all.x=T,sort=F)
#hsigm = res[res$FDR<=0.001 & res$logFC<=-3,]
#hsigm = merge(hsigm,cpm,by=0,all.x=T,sort=F)
#nrow(hsigp)
#nrow(hsigm)
#hsigp = merge(hsigp,anno,by.x='Row.names',by.y='CompName',sort=F)
#hsigm = merge(hsigm,anno,by.x='Row.names',by.y='CompName',sort=F)
#write.table(hsigp,file='higly_significant_plus_statistics_annotated.txt',sep='\t',quote=F,row.names=F)
#write.table(hsigm,file='higly_significant_minus_statistics_annotated.txt',sep='\t',quote=F,row.names=F)
