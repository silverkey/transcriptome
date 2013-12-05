library(edgeR)
library(gplots)
library(RColorBrewer)

# Load the counts
t = read.table(file='tetdev.combined.length.count.txt',sep='\t',head=T,stringsAsFactors=F)
counts = t[,c(3:5)]
rownames(counts) = t$name
# Load the samples file
samples = read.table(file='samples.txt',sep='\t',quote='',comment.char='',head=T)
# Load the annotations
anno = read.table(file='tetraodon_annocript/output/tetraodon_uniref_filt_ann_out.txt',sep='\t',head=T,stringsAsFactors=F,quote='',comment='',row.names=1)
# Filter out all the transcripts shorter than 200 bp
ok = rownames(anno[anno$QueryLength>=200,])
counts = counts[ok,]
anno = anno[ok,]

# Select transcripts to be used to calculate the dispersion

# Get the number of reads totally mapping on each transcript
s = rowSums(counts)
# Calculate the reads per million
c = s * 1000000 / sum(s)
# Select the transcript recognized with at least 1 per million mapped reads
sel = s[c>=1]
v = counts[names(sel),]
# Calculate the standard deviation for the selected transcripts using counts per million (cpm)
dge = DGEList(v, group=samples$type)
cpm = cpm(dge,normalized.lib.sizes=F,log=FALSE,prior.count=0)
v = apply(cpm,1,function(x)sd(as.numeric(x)))
# Order by decreasing standard deviation
v = v[order(v)]
# Give a look to the selected with the smallest standard deviations ideally the housekeeping
counts[names(head(v)),]

write.table(cpm, file='counts_per_million.txt', sep='\t', quote=F)

# Take only info related to the hsp id in order to get a grouping variable for alternative transcripts
hsp = anno[,c(3,9)]
# Merge annotations with counts
ac = merge(counts,hsp,by=0)
rownames(ac) = ac$Row.names
# Take the first 100 different annotation id from the ordered selected
uhsp = unique(ac[names(v),]$HSPName)[1:101]
# Take out the ones tha d not get hits in uniref(annotated as '-')
uhsp = uhsp[grep('UniRef',uhsp)]
# Create the variable that will contain the transcripts to be used as housekeeping
check = rownames(ac[ac$HSPName %in% uhsp,])

# DE analysis using edgeR, calculating the common dispersion on the list of small variance transcripts (check) and multipling 
# its value by a factor of 10, in order to be sure to give a high dispersion to the value to limit the false positives
# Following code is taken by the vedgeR vignette in case of analysis withou replicates
d = DGEList(counts, group=samples$type)
d1 = d
d1$samples$group = 1
d0 = estimateCommonDisp(d1[check,])
d$common.dispersion = d0$common.dispersion * 10

# Filter lowly expressed genes
cpm = cpm(d,normalized.lib.sizes=F,log=FALSE,prior.count=0)
ridx = rowSums(cpm >= 0.5) >= 1
table(ridx)
d = d[ridx,]
cpm = cpm[ridx,]

# Comparison eggs vs 30epi
diff.ee = exactTest(d,pair=c('eggs','X30epi'))
tt = topTags(diff.ee,n=nrow(diff.ee))
res = tt$table
sig.ee = res[res$FDR<=0.05 & abs(res$logFC)>=log2(2),]
write.table(merge(sig.ee, cpm, by=0, all.x=T, sort=F), file='significant_eggs_vs_30epi.txt', sep='\t', quote=F, row.names=F)

# Comparison eggs vs 24hpf
diff.eh = exactTest(d,pair=c('eggs','X24hpf'))
tt = topTags(diff.eh,n=nrow(diff.eh))
res = tt$table
sig.eh = res[res$FDR<=0.05 & abs(res$logFC)>=log2(2),]
write.table(merge(sig.eh, cpm, by=0, all.x=T, sort=F), file='significant_eggs_vs_24hpf.txt', sep='\t', quote=F, row.names=F)

# Comparison 30epi ve 24hpf
diff.eh2 = exactTest(d,pair=c('X30epi','X24hpf'))
tt = topTags(diff.eh2,n=nrow(diff.eh2))
res = tt$table
sig.eh2 = res[res$FDR<=0.05 & abs(res$logFC)>=log2(2),]
write.table(merge(sig.eh2, cpm, by=0, all.x=T, sort=F), file='significant_30epi_vs_24hpf.txt', sep='\t', quote=F, row.names=F)

# Classify groups of differentially expressed according to expression dinamycs
maternal.ee = rownames(sig.ee[sig.ee$logFC<0,])
maternal.eh = rownames(sig.eh[sig.eh$logFC<0,])
embrionic.ee = rownames(sig.ee[sig.ee$logFC>0,])
embrionic.eh = rownames(sig.eh[sig.eh$logFC>0,])
decreasing.eh2 = rownames(sig.eh2[sig.eh2$logFC<0,])
increasing.eh2 = rownames(sig.eh2[sig.eh2$logFC>0,])
maternal = intersect(maternal.ee,maternal.eh)
embrionic = intersect(embrionic.ee,embrionic.eh)

# Write tables with groups of genes showing dinamyc expression patterns
cpm2 = as.data.frame(cpm)
write.table(cpm2[maternal,], file='maternal.txt', sep='\t', quote=F, row.names=T)
write.table(cpm2[embrionic,], file='embrionic.txt', sep='\t', quote=F, row.names=T)

# Drawings
pdf(file='experimental_charts.pdf',paper='a4r',width=8.3,height=11.7,pointsize=8)
# Venn diagram
egg_vs_epi = ifelse(rownames(cpm) %in% rownames(sig.ee),1,0)
egg_vs_hpf = ifelse(rownames(cpm) %in% rownames(sig.eh),1,0)
epi_vs_hpf = ifelse(rownames(cpm) %in% rownames(sig.eh2),1,0)
vennDiagram(data.frame(egg_vs_epi,egg_vs_hpf,epi_vs_hpf),circle.col=c(2,3,4),lwd=3,cex=c(2,2.5,3))
# CLUSTERING
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
tocl = cpm[unique(c(maternal,embrionic)),]
d = dist(tocl)
c = hclust(d,method='s')
heatmap.2(tocl[c$labels,],dendrogram='none',Rowv=F,Colv=F,scale='row',trace="none",col = rev(hmcol),keysize=1,density.info='none',labRow=NA,margins=c(10,1))
# Histogram for the distribution of cpm
hist(cpm,breaks=500000,xlim=c(0,10))
abline(v=0.5,col='green')
abline(v=0.05,col='red')
dev.off()
