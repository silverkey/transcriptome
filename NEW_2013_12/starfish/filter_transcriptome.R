library(edgeR)

targetfile = 'sample_file'
datafile = 'COUNTS.txt'
milfilt = 1
minlibfilt = 1

samples = read.table(file=targetfile,sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
rownames(samples) = samples$name
data = read.table(datafile,header=T,quote='')

rownames(data) = data$tid
data = data[,samples$name]

tot = colSums(data)
lib.size=tot[samples$name]

dge = DGEList(data, group=samples$condition,lib.size=lib.size[samples$name])

cpm = cpm(dge,normalized.lib.sizes=F,log=FALSE,prior.count=0)
ridx = rowSums(cpm >= milfilt) >= minlibfilt
print(table(ridx))

dge = dge[ridx,]
cpm = cpm[ridx,]
write.table(dge$counts,file='transcripts_counts_filtered.txt',sep='\t',quote=F)
write.table(cpm,file='transcripts_cpm_filtered.txt',sep='\t',quote=F)

