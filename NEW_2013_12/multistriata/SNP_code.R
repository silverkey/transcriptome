t = read.table(file='transcriptome_mapping_SNP_calls.snptab',sep='\t',header=T)
s = read.table(file='../significant_statistics_annotated.txt',sep='\t',head=T,comment.char='', quote='')
m = read.table(file='../mapping_locations.txt',sep='\t',head=T)

t$h = ifelse((t$CIIO_dp>=100 & t$CIIP_dp>=100 & t$HATT_dp>=100 & t$HCUH_dp>=100 & t$HCUN_dp>=100 & t$HCUO_dp>=100),1,0)
t$hdiv = ifelse((t$h==1 & t$CIIP_gt==t$HATT_gt & t$HATT_gt==t$HCUN_gt & t$CIIO_gt==t$HCUO_gt & t$HCUO_gt==t$HCUH_gt & t$CIIP_gt!=t$CIIO_gt),1,0)

cand = setdiff(t[t$h==1 & t$hdiv==1,]$seqid,t[t$h==1 & t$hdiv==0,]$seqid)
loc = table(m[m$transcript %in% cand,]$chr)

aaa = as.data.frame(table(t[t$seqid %in% cand & t$h==1 & t$hdiv==1,]$seqid))
aaa[aaa$Freq>7,]

# END

> aaa[aaa$Freq>7,]
                   Var1 Freq
10290 comp27708_c0_seq1    8
10291 comp27708_c0_seq2    9
11784 comp28728_c0_seq1    9
13000 comp29486_c0_seq1   11
14393 comp30273_c0_seq1   12
14799 comp30480_c0_seq1    8
15074 comp30647_c0_seq3    8


h = t[t$CIIO_dp>=100 & t$CIIP_dp>=100 & t$HATT_dp>=100 & t$HCUH_dp>=100 & t$HCUN_dp>=100 & t$HCUO_dp>=100,]
l = t[t$CIIO_dp>=50 & t$CIIP_dp>=50 & t$HATT_dp>=50 & t$HCUH_dp>=50 & t$HCUN_dp>=50 & t$HCUO_dp>=50,]

ldiv = l[l$CIIP_gt==l$HATT_gt & l$HATT_gt==l$HCUN_gt & l$CIIO_gt==l$HCUO_gt & l$HCUO_gt==l$HCUH_gt & l$CIIP_gt!=l$CIIO_gt,]
hdiv = h[h$CIIP_gt==h$HATT_gt & h$HATT_gt==h$HCUN_gt & h$CIIO_gt==h$HCUO_gt & h$HCUO_gt==h$HCUH_gt & h$CIIP_gt!=h$CIIO_gt,]

lequ = l[l$CIIP_gt==l$HATT_gt & l$HATT_gt==l$HCUN_gt & l$CIIO_gt==l$HCUO_gt & l$HCUO_gt==l$HCUH_gt & l$CIIP_gt==l$CIIO_gt,]
hequ = h[h$CIIP_gt==h$HATT_gt & h$HATT_gt==h$HCUN_gt & h$CIIO_gt==h$HCUO_gt & h$HCUO_gt==h$HCUH_gt & h$CIIP_gt==h$CIIO_gt,]

s = read.table(file='../significant_statistics_annotated.txt',sep='\t',head=T,comment.char='', quote='')
m = read.table(file='../mapping_locations.txt',sep='\t',head=T)
c = table(m[m$transcript %in% hdiv$seqid,]$chr)

hdiv.m = merge(cand,m,by.x='seqid',by.y='transcript')
