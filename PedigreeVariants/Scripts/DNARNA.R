################

library(vcfR)
library(ggpubr)
library(VennDiagram)
library(ggplot2)
library(dplyr)
library(ggpmisc)
library(tidyverse)
library(ggsignif)
library(broom)
################################################################
################################################################

RNA <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/FILES/SRR1153470.vcf", verbose = FALSE)
DNA <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ERR194147.vcf", verbose = FALSE)

dna <- create.chromR(name="Supercontig", vcf=DNA, verbose=FALSE)
rna <- create.chromR(name="Supercontig", vcf=RNA, verbose=FALSE)

################################################################
#One function for all ohahahaha
################################################################
VAFcalculation <- function(Cr, path) {
    dp  <-  extract.gt(Cr, element="DP", as.numeric = TRUE)
    ad  <-  extract.gt(Cr, element = 'AD')
    RA <- masplit(ad, record = 1, sort = 0)
    RAF <- RA/dp
    VAF <- 1-RAF
    dframe=data.frame(cbind(RAF, VAF))
    #write.table(dframe, file=path,sep="\t",quote=T, na = "NA")
    return(dframe)
}
#AAAAAMMMMIIINNNNNAAAA do not run it will delet all the files already modified
dframedna = VAFcalculation(dna,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/results/dnaVAF.txt")
dframerna = VAFcalculation(rna,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/results/rnaVAF.txt")
Tdna<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/results/dnaVAF.txt", sep="\t", header=T)
Trna<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/results/rnaVAF.txt", sep="\t", header=T)
y=Tdna$pos
z=Trna$pos
################################################################
#different between  rDNA and rRNA variants
################################################################
#### Variants present only in the daughter
getdiff <- function( cc,x,y,dataf) {
    dd=c()
    ddall=c()
    cc
    print (length(cc))
    for (i in 1:length (cc)){
        ccc=paste("KY962518.1_",cc[i],sep = "",collapse = ' ')
        print (ccc)
        vv=cbind(ccc,dataf[ccc,"X20.1"] )
        print ( vv)
        ddall=rbind(ddall, vv)
        print(dataf[ccc,"X20.1"])
        dd=c(dd,dataf[ccc,"X20.1"])
    }
    print (dd)
    return(list(val1=ddall, val2=dd))
}
#how to run
#__________________________________________________________
#__________________________________________________________
cc=setdiff(y,z)
differences=getdiff(cc,y,z, dframedna) #(choisir le dataframe qui contiens les differents ;)
ddall=differences$val1
dF=differences$val2
#__________________________________________________________
diffDFrnadna=data.frame(ddall)
#write.table(diffDFrnadna, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/Figures/diffDNA_RNA.txt",sep="\t",quote=T, na = "NA")
#__________________________________________________________
#__________________________________________________________
cc=setdiff(z,y)
differences=getdiff(cc,y,z, dframerna) #(choisir le dataframe qui contiens les differents ;)
drall=differences$val1
drF=differences$val2
#__________________________________________________________
diffDFdnarna=data.frame(drall)
#write.table(diffDFdnarna, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/Figures/diffRNA_DNA.txt",sep="\t",quote=T, na = "NA")
diff=c(rep("RNA", nrow(drall)), rep("DNA", nrow(ddall)))
all=c( drF,dF)
differenceonly=data.frame(diff, all)

################################################################
#intersect
################################################################
getintersect <- function(CC,dataS,dataF) {
    intD=c()
    intF=c()
    dd=c()
    vv=c()
    for (i in 1:length (CC)){
        ccc=paste("KY962518.1_",CC[i],sep = "",collapse = ' ')
        vv=c(vv,ccc)
        intD=c(intD,dataS[ccc,"X20.1"])
        intF=c(intF,dataF[ccc,"X20.1"])
        
    }
    return(list(val0=vv, val1=intD, val2=intF))
}


y=Tdna$pos
z=Trna$pos
CC=intersect(z,y)
intersect=getintersect(CC, dframedna, dframerna) #(choisir le dataframe qui contiens les differents ;)
vvv=intersect$val0
dddna=intersect$val1
ddrna=intersect$val2
intDFrnadna=data.frame(cbind(vvv, dddna, ddrna))
#write.table(intDFrnadna, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/Figures/intDNA_RNA.txt",sep="\t",quote=T, na = "NA")
#__________________________________________________________
#Creat intersect table for plotting
#__________________________________________________________
intersect<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/Figures/intDNA_RNA.txt", sep="\t", header=T)
NA12878=c(rep("RNA", length(intersect$dddna.1)), rep("DNA", length(intersect$ddrna)))
all=c(intersect$dddna.1,intersect$ddrna)
intersectonly=data.frame(NA12878, all)
#_________________________
fligner.test(all ~ NA12878, data = intersectonly)
bartlett.test(all ~ NA12878, data = intersectonly)

#__________________________________________________________
# plot intersecton
#__________________________________________________________

my_comparisons <- list( c("RNA", "DNA"))
p<-ggboxplot(intersectonly, x = "NA12878", y ="all", fill="NA12878",palette = c("#EFC000FF","#7AA6DCFF","#003C67FF"), add = "jitter",ylab="Variable Alelle Frequency (VAF) ", xlab="")+
theme_light()+
theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold"),legend.position="none")
#+
#stat_compare_means(paired = TRUE)
stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test", aes(label=..p.adj..) )+
stat_compare_means(label.y = 1.2)
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/Figures/DNA_RNA_boxplot.pdf",width=8,height=4)
plot(p)
dev.off()

#__________________________________________________________
#Plot difference only
#__________________________________________________________

my_comparisons <- list( c("RNA", "DNA"))
p<-ggboxplot(differenceonly, x = "diff", y ="all", fill="diff",palette = c("#EFC000FF","#7AA6DCFF","#003C67FF"), add = "jitter",ylab="Variable Alelle Frequency (VAF) ", xlab="")+
theme_light()+
theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold"),legend.position="none")+
stat_compare_means(paired = TRUE)
#+
#stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test", aes(label=..p.adj..) )+
#stat_compare_means(label.y = 1.2)
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/Figures/diff_DNA_RNA_boxplot.pdf",width=8,height=4)
plot(p)
dev.off()


#__________________________________________________________
#______________________Correlation ____________________________________
#__________________________________________________________
intersect<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/Figures/intDNA_RNA.txt", sep="\t", header=T)
fit <- lm(intersect$"dddna.1" ~ intersect$"ddrna")
formula <- y ~ x
p2<-ggplot(intersect, aes(x=intersect$"dddna.1", y=intersect$"ddrna")) + geom_point()+theme_light()+
#geom_smooth(method=lm )+
theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"), axis.title.y = element_text(color="black", size=10, face="bold"),legend.position="top")+ #ggtitle("Luad WXS") +
geom_smooth(method = "lm", formula = formula, se = F) +
stat_poly_eq(aes(label = paste(..rr.label..)),label.x.npc = "right", label.y.npc = 0.25,formula = formula, parse = TRUE, size = 4)+
xlab("rDNA") +
ylab("rRNA") +
stat_fit_glance(method = 'lm',method.args = list(formula = formula),geom = 'text',aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.35, size = 4)
plot(p2)
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/RNASeq/Figures/corr_DNA_RNA_.pdf",width=4,height=3)
plot(p2)
dev.off()




####extrac

NA12878=c(rep("DNA", 1), rep("RNA", 1))
x= c(0.2,1,0.54, 0.71 )
y=c(0.46,0.4,0.24,0.54)
all=c(x,y)
intersectonly=data.frame(NA12878, all)

group_by(intersectonly, NA12878) %>%
summarise(
count = n(),
median = median(all, na.rm = TRUE),
IQR = IQR(all, na.rm = TRUE)
)


res <- t.test(all~NA12878, data=intersectonly)#, paired = TRUE)
res

fisher.test(matrix(c(0.2,1,0.54, 0.71,0.46,0.4,0.24,0.54), nrow = 2))
