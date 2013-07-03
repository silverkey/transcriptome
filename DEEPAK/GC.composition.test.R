# READ THE TABLES AND EXTRACT THE GC PROPORTION FOR EACH TRANSCRIPT
apo = read.table(file='aporus.txt',sep='\t',comment.char='',quote='',stringsAsFactors=F,head=T)$GC
dan = read.table(file='danicus.txt',sep='\t',comment.char='',quote='',stringsAsFactors=F,head=T)$GC

# SORT THE VALUES IN RESPECT TO THE GC PROPORTION
apo = apo[order(apo)]
dan = dan[order(dan)]

# CREATE A GROUPING FACTOR TO SPLIT EACH GROUP OF SORTED VALUES IN 20 EQUALLY POPULATED INTERVALS
apo.gr = cut(c(1:length(apo)),20)
dan.gr = cut(c(1:length(dan)),20)

# CALCULATE THE MEAN GC PROPORTION IN EACH INTERVAL
apo.m = tapply(apo,apo.gr,mean)
dan.m = tapply(dan,dan.gr,mean)

# EXECUTE THE Mann-Whitney U-test
# IN R THE TEST IS DONE USING THE wilcox.test FUNCTION ON 2 VECTORS OF DATA
# SEE ?wilcox.test FOR EXPLAINATION
wilcox.test(apo.m,dan.m)


# p-value = 0.0009334

