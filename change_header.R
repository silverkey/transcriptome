t = read.table(file='transcripts.counts.matrix',head=T)
colnames(t) = gsub('.isoforms.results','',colnames(t))
write.table(t,file='transcripts.counts.matrix',sep="\t",quote=F)

g = read.table(file='genes.counts.matrix',head=T)
colnames(g) = gsub('.genes.results','',colnames(g))
write.table(g,file='genes.counts.matrix',sep="\t",quote=F)

