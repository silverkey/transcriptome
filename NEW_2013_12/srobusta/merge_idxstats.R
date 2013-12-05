get.name = function(name) {
  name = gsub('.counts','',name)
  name
}

get.df = function(filename) {
  t = read.table(file=filename)[,c(1,3)]
  name = get.name(filename)
  colnames(t) = c('tid',name)
  t
}

cfiles = dir()[grep('.counts',dir())]

counts = get.df(cfiles[1])

for(i in 2:length(cfiles)) {
  t = get.df(cfiles[i])
  counts = merge(counts,t,by='tid',all.x=T,all.y=T,sort=F)
}

write.table(counts,file='COUNTS.txt',sep="\t",quote=F,row.names=F)


#reference sequence name, sequence length, # mapped reads and # unmapped reads
