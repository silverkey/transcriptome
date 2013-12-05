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


