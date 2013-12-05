library(edgeR)

targetfile = 'sample_file'
datafile = 'COUNTS.txt'
milfilt = 0.5
minlibfilt = 2

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

dge2 = dge[ridx,]
cpm2 = cpm[ridx,]
write.table(dge2$counts,file='transcripts_counts_filtered.txt',sep='\t',quote=F)
write.table(cpm2,file='transcripts_cpm_filtered.txt',sep='\t',quote=F)













dge2$samples$lib.size = colSums(dge2$counts)

dge2 = calcNormFactors(dge2)
dge2 = estimateCommonDisp(dge2)
dge2 = estimateTagwiseDisp(dge2)

diff = exactTest(dge2,pair=c('supra','sub'))
tt = topTags(diff,n=nrow(diff))

res = tt$table
res = merge(res,cpm,by=0,all.x=T,sort=F)






































#------------------------------------------------------------------------------
# EDGER ANALYSIS
#------------------------------------------------------------------------------

library(edgeR)
library(org.Hs.eg.db)

targetfile = 'targets.txt'
datafile = 'COUNTS.txt'

dbidname = 'REFSEQ'
idfiltchar = 'N'

milfilt = 1
minlibfilt = 2
filtname = 'filt102'

fdrfilt = 0.1
fdrname = 'fdr01'

fcfilt = 1.5
fcname = 'fc15.xls'

comparisons = list(c('FC','FV'),
                   c('FC','FP'),
                   c('FC','FPV'),
                   c('HC','HV'),
                   c('HC','HP'),
                   c('HC','HPV'))

analyze.comparison = function(comparison) {

  colData = samples[samples$condition %in% comparison,]
  countData = data[,colData$name]

  dge = DGEList(countData, group=colData$condition,lib.size=lib.size[colData$name])

  cpm = cpm(dge,normalized.lib.sizes=F,log=FALSE,prior.count=0)
  ridx = rowSums(cpm >= milfilt) >= minlibfilt
  print(comparison)
  print(table(ridx))
  dge2 = dge[ridx,]

  dge2$samples$lib.size = colSums(dge2$counts)

  dge2 = calcNormFactors(dge2)
  dge2 = estimateCommonDisp(dge2)
  dge2 = estimateTagwiseDisp(dge2)

  diff = exactTest(dge2,pair=comparison)
  tt = topTags(diff,n=nrow(diff))

  res = tt$table
  res = merge(res,cpm,by=0,all.x=T,sort=F)

  res = merge(res, annots,by.x='Row.names',by.y=dbidname,sort=F)
  res$FC = ifelse(res$logFC>0,2^abs(res$logFC),-2^abs(res$logFC))
  res.out = paste(comparison[1],comparison[2],'all',filtname,'expressed.xls',sep='_')

  sig = res[res$FDR<=fdrfilt & abs(res$logFC)>=log2(fcfilt),]
  sig = sig[,c('Row.names','SYMBOL','FC','FDR','GENENAME',colnames(cpm))]
  colnames(sig) = c(dbidname,'symbol','fold_change','fdr','gene_name',colnames(cpm))

  sig.out = paste(comparison[1],comparison[2],filtname,fdrname,fcname,sep='_')

  write.table(res,file=res.out,sep="\t",quote=F,row.names=F)
  write.table(sig,file=sig.out,sep="\t",quote=F,row.names=F)
  print(' ')
  res
}

samples = read.table(file=targetfile,sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
rownames(samples) = samples$name

data = read.table(datafile,header=T,quote='')
colnames(data) = gsub('.','-',colnames(data),fixed=T)
rownames(data) = data$geneid
data = data[,samples$name]

tot = colSums(data)
lib.size=tot[samples$name]

data = data[grep(idfiltchar,rownames(data)),]
notass = data[-grep(idfiltchar,rownames(data)),]

id = rownames(data)
annots = select(org.Hs.eg.db, keys=id, cols=c('SYMBOL','GENENAME'), keytype=dbidname)

comres = list()
for(i in 1:length(comparisons)) {
  comres[[i]] = analyze.comparison(comparisons[[i]])
}

dge = DGEList(data)
cpm = cpm(dge,normalized.lib.sizes=F,log=FALSE,prior.count=0)
tot = merge(cpm,annots,by.x=0,by.y=dbidname,sort=F)
tot.out = paste('ALL_EXPRESSION',datafile,sep='_')
write.table(tot,file=tot.out,sep="\t",quote=F,row.names=F)
