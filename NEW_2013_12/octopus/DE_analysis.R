library(edgeR)

targetfile = 'sample_file'
datafile = 'COUNTS.txt'
annotfile = 'octopus_rnaseq2_uniref_filt_ann_out.txt'

milfilt = 0.5
minlibfilt = 2
filtname = 'filt0502'

fdrfilt = 0.001
fdrname = 'fdr001'

fcfilt = 2
fcname = 'fc2.xls'

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
  
  res = merge(res, annots,by.x='Row.names',by.y='CompName',sort=F)
  res$FC = ifelse(res$logFC>0,2^abs(res$logFC),-2^abs(res$logFC))
  res.out = paste(comparison[1],comparison[2],'all',filtname,'expressed.xls',sep='_')
  
  sig = res[res$FDR<=fdrfilt & abs(res$logFC)>=log2(fcfilt),]

  sig.out = paste(comparison[1],comparison[2],filtname,fdrname,fcname,sep='_')
  
  write.table(res,file=res.out,sep="\t",quote=F,row.names=F)
  write.table(sig,file=sig.out,sep="\t",quote=F,row.names=F)
  print(' ')
  res
}

# Make pairwise comparisons
samples = read.table(file=targetfile,sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
rownames(samples) = samples$name
data = read.table(datafile,header=T,quote='')
rownames(data) = data[,1]
data = data[,samples$name]
colsums = colSums(data)
lib.size = colsums[samples$name]
annots = read.table(file=annotfile,sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
comparisons = combn(unique(samples$condition),2)
comres = list()
for(i in 1:ncol(comparisons)) {
  comres[[i]] = analyze.comparison(comparisons[,i])
}

# Collect and save general data for all the transcriptome
dge = DGEList(data)
cpm = cpm(dge,normalized.lib.sizes=F,log=FALSE,prior.count=0)
ridx = rowSums(cpm >= milfilt) >= minlibfilt
print(table(ridx))
cpm = cpm[ridx,]
tot = merge(cpm,annots,by.x=0,by.y='CompName',sort=F)
row.names(tot) = tot$Row.names
tot.out = paste('ALL_EXPRESSION',datafile,sep='_')
write.table(tot,file=tot.out,sep="\t",quote=F,row.names=F)

# Counting the whole transcriptome for different features
# 1) Percentages of noncoding for each tissue on the total of the expressed genes for that tissue
sup.ex = sum(tot$IZ10>=0.5 & tot$IZ15>=0.5 & tot$IZ20>=0.5)
sup.nc = sum(tot$IZ10>=0.5 & tot$IZ15>=0.5 & tot$IZ20>=0.5 & tot$lncRNA4Annocript==1)
sup.nc.p = sup.nc/sup.ex*100

sub.ex = sum(tot$IZ11>=0.5 & tot$IZ16>=0.5 & tot$IZ21>=0.5)
sub.nc = sum(tot$IZ11>=0.5 & tot$IZ16>=0.5 & tot$IZ21>=0.5 & tot$lncRNA4Annocript==1)
sub.nc.p = sub.nc/sub.ex*100

opt.ex = sum(tot$IZ12>=0.5 & tot$IZ17>=0.5 & tot$IZ22>=0.5)
opt.nc = sum(tot$IZ12>=0.5 & tot$IZ17>=0.5 & tot$IZ22>=0.5 & tot$lncRNA4Annocript==1)
opt.nc.p = opt.nc/opt.ex*100

arm.ex = sum(tot$IZ13>=0.5 & tot$IZ18>=0.5 & tot$IZ23>=0.5)
arm.nc = sum(tot$IZ13>=0.5 & tot$IZ18>=0.5 & tot$IZ23>=0.5 & tot$lncRNA4Annocript==1)
arm.nc.p = arm.nc/arm.ex*100

barplot(c(sup.nc.p,sub.nc.p,opt.nc.p,arm.nc.p))

# Isolate tissue specific genes
names(comres) = apply(comparisons,2,function(x) paste(x[1],x[2],sep=''))

supsub = comres$suprasub
supopt = comres$supraoptic
suparm = comres$supraarm
subopt = comres$suboptic
subarm = comres$subarm
optarm = comres$opticarm

sup1 = supsub[supsub$FC<=-10 & supsub$FDR<=0.05,'Row.names']
sup2 = supopt[supopt$FC<=-10 & supopt$FDR<=0.05,'Row.names']
sup3 = suparm[suparm$FC<=-10 & suparm$FDR<=0.05,'Row.names']
sup.en = Reduce(intersect,list(sup1,sup2,sup3))
write.table(tot[sup.en,],file='supra_enriched.xls',sep="\t",quote=F,row.names=F)

sub1 = supsub[supsub$FC>= 10 & supsub$FDR<=0.05,'Row.names']
sub2 = subopt[subopt$FC<=-10 & subopt$FDR<=0.05,'Row.names']
sub3 = subarm[subarm$FC<=-10 & subarm$FDR<=0.05,'Row.names']
sub.en = Reduce(intersect,list(sub1,sub2,sub3))
write.table(tot[sub.en,],file='sub_enriched.xls',sep="\t",quote=F,row.names=F)

opt1 = supopt[supopt$FC>= 10 & supopt$FDR<=0.05,'Row.names']
opt2 = subopt[subopt$FC>= 10 & subopt$FDR<=0.05,'Row.names']
opt3 = optarm[optarm$FC<=-10 & optarm$FDR<=0.05,'Row.names']
opt.en = Reduce(intersect,list(opt1,opt2,opt3))
write.table(tot[opt.en,],file='optic_enriched.xls',sep="\t",quote=F,row.names=F)

arm1 = suparm[suparm$FC>= 10 & suparm$FDR<=0.05,'Row.names']
arm2 = subarm[subarm$FC>= 10 & subarm$FDR<=0.05,'Row.names']
arm3 = optarm[optarm$FC>= 10 & optarm$FDR<=0.05,'Row.names']
arm.en = Reduce(intersect,list(arm1,arm2,arm3))
write.table(tot[arm.en,],file='arm_enriched.xls',sep="\t",quote=F,row.names=F)

# > samples
#       animal sex weight                 tissue condition name
# IZ10  13/01   F    296 supra-oesophageal mass     supra IZ10
# IZ11  13/01   F    296   sub-oesophageal mass       sub IZ11
# IZ12  13/01   F    296             Optic lobe     optic IZ12
# IZ13  13/01   F    296          anterior arm        arm IZ13
# IZ15  13/02   M    336 supra-oesophageal mass     supra IZ15
# IZ16  13/02   M    336   sub-oesophageal mass       sub IZ16
# IZ17  13/02   M    336             Optic lobe     optic IZ17
# IZ18  13/02   M    336          anterior arm        arm IZ18
# IZ20  13/03   M    288 supra-oesophageal mass     supra IZ20
# IZ21  13/03   M    288   sub-oesophageal mass       sub IZ21
# IZ22  13/03   M    288             Optic lobe     optic IZ22
# IZ23  13/03   M    288          anterior arm        arm IZ23

# > comparisons
#      [,1]    [,2]    [,3]    [,4]    [,5]  [,6]   
# [1,] "supra" "supra" "supra" "sub"   "sub" "optic"
# [2,] "sub"   "optic" "arm"   "optic" "arm" "arm"  