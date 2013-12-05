ip = read.table(file='mart_export_interpro.txt',sep='\t',head=T,quote='',comment.char='',stringsAsFactors=F)
colnames(ip) = c('id','ipdesc')
a = tapply(ip$ipdesc,ip$id,paste,collapse=', ')
b = as.data.frame(a)
colnames(b) = 'ipdesc'
b$id = rownames(b)
b = b[,c(2,1)]
anno = read.table(file='mart_export_geneinfo.txt',sep='\t',head=T,quote='',comment.char='',stringsAsFactors=F)
colnames(anno) = c('id','desc','chr','start','end','name')
anno = merge(anno,b,by.x='id',by.y='id')
write.table(anno,file='mart_export_anno.txt',sep='\t',quote=F,row.names=F)
annosig = merge(sig,anno,by.x='Row.names',by.y='id',sort=F)
write.table(annosig,file='significant_annotated.txt',sep='\t',quote=F,row.names=F)

