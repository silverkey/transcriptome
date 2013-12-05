# Print just some numbers
nrow(sig.ee)
nrow(sig.eh)
nrow(sig.eh2)

sum(sig.ee$logFC>0)
sum(sig.ee$logFC<0)
sum(sig.eh$logFC>0)
sum(sig.eh$logFC<0)
sum(sig.eh2$logFC>0)
sum(sig.eh2$logFC<0)

length(unique(anno$HSPName))

nrow(cpm2)
nrow(cpm2[cpm2$X24hpf<=0.05 & cpm2$X30epi<=0.05 & cpm2$egg>=0.5,])
nrow(cpm2[cpm2$X24hpf<=0.05 & cpm2$X30epi>=0.5 & cpm2$egg>=0.5,])
nrow(cpm2[cpm2$X24hpf>=0.5 & cpm2$X30epi>=0.5 & cpm2$egg>=0.5,])
nrow(cpm2[cpm2$X24hpf<=0.05 & cpm2$X30epi>=0.5 & cpm2$egg<=0.05,])
nrow(cpm2[cpm2$X24hpf>=0.5 & cpm2$X30epi>=0.5 & cpm2$egg<=0.05,])
nrow(cpm2[cpm2$X24hpf>=0.5 & cpm2$X30epi<=0.05 & cpm2$egg<=0.05,])
nrow(cpm2[cpm2$X24hpf>=0.5 & cpm2$X30epi<=0.05 & cpm2$egg>=0.5,])

hist(cpm,breaks=500000,xlim=c(0,10))
abline(v=0.5,col='red')
abline(v=0.05,col='red')

pdf(file='clustering.pdf',paper='a4r',width=8.3,height=11.7,pointsize=8)

