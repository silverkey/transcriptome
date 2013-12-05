#-------------------------------------------------
# PARAMETERS TO SET UP BEFORE TO RUN THE ANALYSIS
#-------------------------------------------------

# Indicate name and path of the ANNOCRIPT TABLE containing the annotations of the whole transcriptome by ANNOCRIPT
anno = '/home/remo/ANALYSIS/tetraodon/tetraodon_annocript/output/tetraodon_uniref_filt_ann_out.txt'

# Indicate name and path of the file containing a column with the id of a set of transcripts to test for GO enrichments
# They may be differentially expressed transcripts or stage specific transcripts or whatever
# The important thing is that they must derive from the same id present into the ANNOCRIPT TABLE
diff = '/home/remo/ANALYSIS/tetraodon/tetraodon_DE_analysis_20131101/significant_30epi_vs_24hpf.txt'

# Indicate the title of the analysis
# It must be consistent with the "diff" file!!!!
# NO MORE THAN 2 WORDS!!!!
title = '30EPI_vs_24HPF Transcripts'

# Indicate the column containing the id of the transcripts into the "diff" file
# Put 0 if it corresponds to rownames
col.diff = 1

# Indicate name and path of the file containing the mapping between the GO id and the definition
gomap = 'go_map'

# Indicate the minimum number of transcripts associated to a GO class to take it into consideration for the analysis
# It is relative to the number of transcript for each class into the "diff" variable 
min = 5

# Indicate the fold the proportion of "diff" for a class shuold be higher than the universe to analyze it
mult = 1

# Indicate what kind of enrichment you are looking for in the "diff" table
# To look for enrichments use 'g', for impoverishment use 'l', for both use 't'
prop.alt = 'g'

# Adjusted p-value cutoff to consider a class as significant
p.filt = 0.1

# Indicate the maximum number of significant classes to display into the plot
topn.toplot = 15

#-------------------------------------------------------
# END OF PARAMETERS TO SET UP BEFORE TO RUN THE ANALYSIS
#-------------------------------------------------------


#-------------------------------------------------------
# ANALYSIS
#-------------------------------------------------------

map.go.bp = function(table) {
  bpid = strsplit(table$BPId,']---[',fixed=T)
}

map.go.mf = function(table) {
  mfid = strsplit(table$MFId,']---[',fixed=T)
}

map.go.cc = function(table) {
  ccid = strsplit(table$CCId,']---[',fixed=T)
}

get.counts = function(counts) {
  counts = as.data.frame(table(unlist(counts)))
  colnames(counts) = c('goid','count')
  counts = counts[grep('GO:',counts$goid),]
  counts
}

calculate.enrichments = function(div.sel,div.uni,n.sel,n.uni,go) {
  div = merge(div.sel,div.uni,by.x='goid',by.y='goid',all.y=T)
  div$count.x[is.na(div$count.x)] = 0
  div$pval = apply(div,1,function(x) prop.test(c(as.numeric(x[2]),as.numeric(x[3])),c(n.sel,n.uni),alternative=prop.alt)$p.value)
  div = merge(div,go,by.x='goid',by.y='go_id')
  if(prop.alt == 'g') div = div[(div$count.x/n.sel >= mult*(div$count.y/n.uni)) & div$count.x>=min,]
  if(prop.alt == 'l') div = div[(div$count.x/n.sel <= mult*(div$count.y/n.uni)) & div$count.x>=min,]
  if(prop.alt == 't') div = div[div$count.x>=min,]
  div$padj = p.adjust(div$pval)
  div[order(div$padj,decreasing=T),]
}

go = read.table(file=gomap,sep='\t',head=T,quote='',comment.char='',stringsAsFactors=F)
uni.t = read.table(file=anno,sep='\t',head=T,quote='',comment.char='',stringsAsFactors=F)
d = read.table(file=diff,sep='\t',head=T,quote='',comment.char='',stringsAsFactors=F)
if(col.diff==0) {
  sel = unique(rownames(d))
} else sel = unique(d[,col.diff])
sel.t = uni.t[uni.t$CompName %in% sel,]

n.uni = length(unique(uni.t$CompName))
n.sel = length(sel)

# Biological Process
map.uni.bp = map.go.bp(uni.t)
uni.bp = get.counts(map.uni.bp)
map.sel.bp = map.go.bp(sel.t)
sel.bp = get.counts(map.sel.bp)
res.bp = calculate.enrichments(sel.bp,uni.bp,n.sel,n.uni,go)
sig.bp = subset(res.bp,padj<=p.filt)

# Molecular Function
map.uni.mf = map.go.mf(uni.t)
uni.mf = get.counts(map.uni.mf)
map.sel.mf = map.go.mf(sel.t)
sel.mf = get.counts(map.sel.mf)
res.mf = calculate.enrichments(sel.mf,uni.mf,n.sel,n.uni,go)
sig.mf = subset(res.mf,padj<=p.filt)

# Cellular Component
map.uni.cc = map.go.cc(uni.t)
uni.cc = get.counts(map.uni.cc)
map.sel.cc = map.go.cc(sel.t)
sel.cc = get.counts(map.sel.cc)
res.cc = calculate.enrichments(sel.cc,uni.cc,n.sel,n.uni,go)
sig.cc = subset(res.cc,padj<=p.filt)

restab = rbind(sig.bp,sig.mf,sig.cc)
write.table(restab,file=paste(gsub(' ','_',title),'GO_significant.txt',sep='_'),row.names=F,sep='\t',quote=F)

pdf(file=paste(gsub(' ','_',title),'GO_significant.pdf',sep='_'),width=15,height=10)
par(las=2,mar=c(5,25,5,5))

sig.bp = tail(sig.bp,topn.toplot)
barplot(
  t(data.frame(selected=sig.bp$count.x/n.sel*100,universe=sig.bp$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.bp$definition,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'GO Biological Process Enriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

sig.mf = tail(sig.mf,topn.toplot)
barplot(
  t(data.frame(selected=sig.mf$count.x/n.sel*100,universe=sig.mf$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.mf$definition,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'GO Molecular Function Enriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

sig.cc = tail(sig.cc,topn.toplot)
barplot(
  t(data.frame(selected=sig.cc$count.x/n.sel*100,universe=sig.cc$count.y/n.uni*100)),
  beside=T,horiz=T,names=unlist(lapply(sig.cc$definition,substring,1,50)),xlab='percentage of transcripts',
  main=paste('Top',topn.toplot,'GO Cellular Component Anriched Among',title,sep=' '),
  cex.axis=1.2,cex.names=1.2,cex.main=1.2,legend=c('Selected','Transcriptome'))

dev.off()
