# Prepare files based on counts and GO classifications 
file = 'sm7_annot_GOs_20130531_1504.txt'
t = read.table(file=file,head=F,comment.char='',quote='',sep='\t')
colnames(t) = c('tid','desc','goid','godef','group')

goidcount = table(t$goid)
goidmap = unique(t[,c('goid','godef','group')])
goidsort = goidcount[order(as.vector(goidcount),decreasing=T)]
sortcount = data.frame(goid=names(goidsort),counts=as.character(goidsort))

merged = merge(sortcount,goidmap,by.x='goid',by.y='goid',sort=F)
write.table(merged,file='Sm_Fe7_GO_counts.txt',sep='\t',quote=F,row.names=F)

file = 'sm7_contigs.txt'
c = read.table(file=file,sep='\t',quote='',comment.char='',head=T)
tc = merge(t,c[,c('contig_id','paired_aligned')],by.x='tid',by.y='contig_id',all.x=T)

write.table(tc,file='Sm_Fe7_annot_readscount.txt',sep='\t',quote=F,row.names=F)
write.table(unique(tc[,c(1,2,6)]),file='Sm_Fe7_desc_readscount.txt',sep='\t',quote=F,row.names=F)


file = 'sm60_annot_GOs_20131108_1536.txt'
t = read.table(file=file,head=F,comment.char='',quote='',sep='\t')
colnames(t) = c('tid','desc','goid','godef','group')

goidcount = table(t$goid)
goidmap = unique(t[,c('goid','godef','group')])
goidsort = goidcount[order(as.vector(goidcount),decreasing=T)]
sortcount = data.frame(goid=names(goidsort),counts=as.character(goidsort))

merged = merge(sortcount,goidmap,by.x='goid',by.y='goid',sort=F)
write.table(merged,file='Sm_Fe60_GO_counts.txt',sep='\t',quote=F,row.names=F)

file = 'sm60_contigs.txt'
c = read.table(file=file,sep='\t',quote='',comment.char='',head=T)
tc = merge(t,c[,c('contig_id','paired_aligned')],by.x='tid',by.y='contig_id',all.x=T)

write.table(tc,file='Sm_Fe60_annot_readscount.txt',sep='\t',quote=F,row.names=F)
write.table(unique(tc[,c(1,2,6)]),file='Sm_Fe60_desc_readscount.txt',sep='\t',quote=F,row.names=F)



# Fisher for functional enrichment analyses
t7 = read.table(file='Sm_Fe7_GO_counts.txt',head=T,comment.char='',quote='',stringsAsFactors=F,sep='\t')
t60 = read.table(file='Sm_Fe60_GO_counts.txt',head=T,comment.char='',quote='',stringsAsFactors=F,sep='\t')

tot7 = 6606
tot60 = 6351

b = merge(t7[,1:2],t60[,1:2],by='goid')
colnames(b) = c('goid','t7','t60')
b = merge(b,t60[,c(1,3)],by='goid',all.x=T)

pval = apply(b,1,function(x) prop.test(c(as.numeric(x[2]),as.numeric(x[3])),c(tot7,tot60))$p.value)
b$pval = pval
b = b[order(b$pval,decreasing=F),]
b$perc7 = b$t7/tot7*100
b$perc60 = b$t60/tot60*100
write.table(b,file='GO_comparison.txt',sep='\t',quote=F,row.names=F)



# GO t-test based on counts
t7 = read.table(file='Sm_Fe7_annot_readscount.txt',head=T,comment.char='',quote='',stringsAsFactors=F,sep='\t')
t60 = read.table(file='Sm_Fe60_annot_readscount.txt',head=T,comment.char='',quote='',stringsAsFactors=F,sep='\t')
go = unique(c(t7$goid,t60$goid))
gotab = unique(rbind(t7[,c('goid','godef')],t60[,c('goid','godef')]))

pval = c()
for(i in 1:length(go)){
  go7 = t7[t7$goid == go[i],]$paired_aligned
  go60 = t60[t60$goid == go[i],]$paired_aligned
  if(length(go7) > 2 & length(go60) > 2) {
    p = t.test(go7,go60)$p.val
  }
  else p = 10
  pval = c(pval,p)
}

df = data.frame(go,pval)
df = df[order(df$pval,decreasing=F),]
df$adj.p = p.adjust(df$pval,method='fdr')
df = merge(df,gotab,by.x='go',by.y='goid')
df = df[order(df$pval,decreasing=F),]

pdf(paper='a4',width=8.3,height=11.7,pointsize=8)
par(mfrow=c(3,2))

avg7 = c()
avg60 = c()
for(i in 1:nrow(df)){
  go7 = t7[t7$goid == df$go[i],]$paired_aligned
  go60 = t60[t60$goid == df$go[i],]$paired_aligned
  boxplot(go7,go60,main=paste(t7[t7$goid == df$go[i],]$godef[1],df[df$go==df$go[i],]$pval),col=8)
  avg7 = c(avg7,mean(go7))
  avg60 = c(avg60,mean(go60))
}
dev.off()

df$avg7 = avg7
df$avg60 = avg60

write.table(df,file='GO_counts_ttest.txt',sep='\t',quote=F,row.names=F)
