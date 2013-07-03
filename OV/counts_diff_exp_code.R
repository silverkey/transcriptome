library("edgeR")

counts = read.table(file='transcripts.counts.matrix', head=T, quote='', comment.char='', sep='\t')
rownames(counts) = counts[,1]
counts = counts[,2:ncol(counts)]

samples = read.table(file='sample_file.txt',sep='\t',quote='',comment.char='',head=T)

samples = edit(samples)
sample.name = gsub('D1WF1ACXX_','',samples$sample.name)
sample.name = gsub('_13s0024.+$','',sample.name)
samples$name = sample.name

target = unique(samples[,-1])
write.table(target,file='sample_file',quote=F,row.names=F,sep='\t')



----------------------------------

library("edgeR")

counts = read.table(file='transcripts.counts.matrix', head=T, quote='', comment.char='', sep='\t')
rownames(counts) = counts[,1]
counts = counts[,2:ncol(counts)]

target = read.table(file='sample_file',sep='\t',quote='',comment.char='',head=T)

for(i in 1:length(target$name)) {
  counts[,as.character(target$name[i])] = apply(counts[,grep(as.character(target$name[i]),colnames(counts))],1,sum)
  print(as.character(target$name[i]))
}

counts = counts[,c(25:36)]
counts = counts[,target$name]

dge = DGEList(counts, group=target$condition)
dge = calcNormFactors(dge)

m = 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx = rowSums(m > 1) >= 2
table(ridx)
dge = dge[ridx,]
write.table(dge$counts,file='transcripts_counts_filtered2.txt',sep='\t',quote=F)

dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)

supra.sub.diff = exactTest(dge,pair=c('supra','sub'))
supra.sub.tt = topTags(supra.sub.diff,n=nrow(supra.sub.diff))
supra.sub.res = supra.sub.tt$table

supra.sub.sig = supra.sub.res[supra.sub.res$FDR<=0.05,]
nrow(supra.sub.sig)

counts.sig = merge(dge$pseudo.counts,sig,by='row.names',sort=F)

write.table(counts.sig,file='differentially_expressed_transcripts.txt',row.names=F,sep='\t',quote=F)

data = dge$pseudo.counts[,target$name]
names = target$condition

trexprs = t(data)
dist = dist(trexprs)
plot(hclust(dist),labels=names)

data.pca = prcomp(data,scale=F,center=F)

plot3d(data.pca$rotation[,1:3],xlab="PC1",ylab="PC2",zlab="PC3",pch=1,type='s',radius=0.02,col=rainbow(19),size=1)
