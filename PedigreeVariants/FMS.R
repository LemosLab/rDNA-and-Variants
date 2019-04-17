################
library(vcfR)
library(ggpubr)
library(VennDiagram)
library(ggplot2)

vcf <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/PPSmerged.vcf", verbose = FALSE)
chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
chrom <- masker(chrom, min_DP = 300, max_DP = 700)
chrom <- proc.chromR(chrom, verbose = FALSE)
################
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
head(dp)
heatmap.bp(dp[10:50,])
is.na(dp[na.omit(dp == 0)]) <- TRUE
heatmap.bp(dp[1001:1500,])
ad <- extract.gt(chrom, element="AD", as.numeric=TRUE)
ad <- extract.gt(vcf, element = 'AD')
################
par(mar=c(8,4,4,2))
barplot(apply(dp, MARGIN=2, mean, na.rm=TRUE), las=3)
################
par(mar=c(8,4,1,1))
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/PPSdepth.pdf",width=8,height=4)
boxplot(dp, col=c("#C0C0C0", "#808080"), ylab="Depth",log='y', las=2)
#abline(h=seq(0,1e4, by=100), col="#C0C0C088")
dev.off()

################################################################
################################################################
################################################################
#Already calculated and the files are created
M <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/NA12890.vcf", verbose = FALSE)
F <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/NA12889.vcf", verbose = FALSE)
D <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/NA12877.vcf", verbose = FALSE)

cM <- create.chromR(name="Supercontig", vcf=M, verbose=FALSE)
cF <- create.chromR(name="Supercontig", vcf=F, verbose=FALSE)
cD <- create.chromR(name="Supercontig", vcf=D, verbose=FALSE)

#Mother
dpM  <-  extract.gt(cM, element="DP", as.numeric = TRUE)
adM  <-  extract.gt(cM, element = 'AD')
RAM <- masplit(adM, record = 1, sort = 0)
RAFM <- RAM/dpM
VAFM <- 1-RAFM
dframeM=data.frame(cbind(RAFM, VAFM))

#Dauther
dpD  <-  extract.gt(cD, element="DP", as.numeric = TRUE)
adD <-  extract.gt(cD, element = 'AD')
RAD <- masplit(adD, record = 1, sort = 0)
RAFD <- RAD/dpD
VAFD <- 1-RAFD
dframeD=data.frame(cbind(RAFD, VAFD))

#Father
dpF  <-  extract.gt(cF, element="DP", as.numeric = TRUE)
adF  <-  extract.gt(cF, element = 'AD')
RAF <- masplit(adF, record = 1, sort = 0)
RAFF <- RAF/dpF
VAFF <- 1-RAFF
dframeF=data.frame(cbind(RAFF, VAFF))

#write.table(dframeF, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addf.txt",sep="\t",quote=T, na = "NA")
#write.table(dframeM, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addM.txt",sep="\t",quote=T, na = "NA")
#write.table(dframeD, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addS.txt",sep="\t",quote=T, na = "NA")

################################################################
#Avant de lire les fichier, j'ai modifier manuellement les ficier pour ajouter une colonne pos
TF<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addf.txt", sep="\t", header=T)
TD<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addS.txt", sep="\t", header=T)
TM<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addM.txt", sep="\t", header=T)

#names=c(rep("PF", nrow(TF)), rep("P", nrow(TM)), rep("S", nrow(TD)))
names=c(rep("NA12889", nrow(TF)), rep("NA12890", nrow(TM)), rep("NA12877", nrow(TD)))
TTALL=data.frame(names, rbind(TF,TM,TD))

#__________________________________________________________
##Veenplot

pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/PPS-veen.pdf")
venn.plot <- venn.diagram(list(TF$pos, TD$pos, TM$pos),NULL,col = "transparent", fill=c("#CD534CFF",  "#0073C2FF","#868686FF"), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface = "bold", category.names=c("NA12890", "NA12877","NA12889" ) , na.rm = FALSE)
grid.draw(venn.plot)
dev.off()

################################################################
#different between  M, D, F # variants
################################################################
####Variants present only in the daughter
dd=c()
ddall=c()
cc=setdiff(setdiff(TD$pos, TF$pos), TM$pos)
for (i in 1:length (cc)){
    ccc=paste("KY962518.1_",cc[i],sep = "",collapse = ' ')
    vv=cbind(ccc,dframeD[ccc,"X20.1"] )
    ddall=rbind(ddall, vv)
    print(dframeD[ccc,"X20.1"])
    dd=c(dd,dframeD[ccc,"X20.1"])
}
dd
ddall
####Variants present only in the Mother #375
dM=c()
dmall=c()
cc=setdiff(setdiff(TM$pos, TF$pos), TD$pos)
for (i in 1:length (cc)){
    ccc=paste("KY962518.1_",cc[i],sep = "",collapse = ' ')
    dM=c(dM,dframeM[ccc,"X20.1"])
    vv=cbind(ccc,dframeM[ccc,"X20.1"] )
    dmall=rbind(dmall, vv)
}
dM
dmall
####Variants present only in the Father #12
dF=c()
dfall=c()
cc=setdiff(setdiff(TF$pos, TM$pos), TD$pos)
for (i in 1:length (cc)){
    ccc=paste("KY962518.1_",cc[i],sep = "",collapse = ' ')
    dF=c(dF,dframeF[ccc,"X20.1"])
    vv=cbind(ccc,dframeF[ccc,"X20.1"] )
    dfall=rbind(dfall, vv)
}
dF
dfall
#__________________________________________________________
names=c(rep("NA12889", length(dF)), rep("NA12890", length(dM)), rep("NA12877", length(dd)))
alldfmd=c( dF, dM, dd)
diffonly=data.frame(names, c( dF, dM, dd))
#__________________________________________________________
#creat the fale that contain the locationand the VAF of each difference.
names=c(rep("NA12889", length(dfall)), rep("NA12890", length(dmall)), rep("NA12877", length(ddall)))
diffonltall=data.frame(cbind(names, c(dfall,dmall,ddall)))
write.table(diffonltall, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/alldifference.txt",sep="\t",quote=T, na = "NA")

################################################################
#intersect only present in all M, D, F #217 variants
################################################################
intD=c()
intF=c()
intM=c()
cc=intersect(intersect(TD$pos, TF$pos), intersect(TD$pos, TM$pos))

for (i in 1:length (cc)){
    ccc=paste("KY962518.1_",cc[i],sep = "",collapse = ' ')
    vv=c(vv,ccc)
    intD=c(intD,dframeD[ccc,"X20.1"])
    intF=c(intF,dframeF[ccc,"X20.1"])
    intM=c(intM,dframeM[ccc,"X20.1"])
    
}
#__________________________________________________________
#creat the fale that contain the locationand the VAF of each difference.
ddd=cbind (vv, intD,intF,intM)
write.table(ddd, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/MFDintersectall.txt",sep="\t",quote=T, na = "NA")
#__________________________________________________________
length(intD)
length(intF)
length(intM)

names=c(rep("NA12889", length(intF)), rep("NA12890", length(intM)), rep("NA12877", length(intD)))
intersectonly=data.frame(names, c( intF, intM, intD))
library(ggplot2)

################################################################
################### Ploter avec des couleurs ###################
################################################################

#colors jco : c("NA12877"="#0073C2FF","NA12878"="#EFC000FF","NA12890"="#868686FF","NA12889"="#CD534CFF","N112891"="#7AA6DCFF","N112892"="#003C67FF"))
#__________________________________________________________

names=c(rep("N112889", nrow(TF)), rep("N112890", nrow(TM)), rep("NA12877", nrow(TD)))
TTALL=data.frame(names, rbind(TF,TM,TD))
TTALL$names=factor(TTALL$names, levels=c("NA12877","N112890","N112889"))
#__________________________________________________________
#all
p<-ggboxplot(TTALL, x = "names", y ="X20.1", shape="names", color="names",palette = c("#0073C2FF","#CD534CFF","#868686FF"), add = "jitter",ylab="Variable Alelle Frequency (VAF) ", xlab="")+theme_light()+theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"),axis.title.y = element_text(color="black", size=10, face="bold"))

pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/PPS-colors-VAF.pdf",width=7,height=3)
plot(p)
dev.off()

#__________________________________________________________
#difference only
diffonly$names=factor(diffonly$names, levels=c("NA12877","NA12890","NA12889"))
p<-ggboxplot(diffonly, x = "names", y ="c.dF..dM..dd.", shape="names", color="names",palette = c("#0073C2FF","#CD534CFF","#868686FF"), add = "jitter",ylab="Variable Alelle Frequency (VAF) ", ylim=c(0,1),xlab="")+
theme_light()+theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"),axis.title.y = element_text(color="black", size=10, face="bold"))

pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/PPSdiff-colored.pdf",width=7,height=3)
plot(p)
dev.off()

#__________________________________________________________
#intersect
intersectonly$names=factor(intersectonly$names, levels=c("NA12877","NA12890","NA12889"))
p<-ggboxplot(intersectonly, x = "names", y ="c.intF..intM..intD.", shape="names", color="names",palette = c("#0073C2FF","#CD534CFF","#868686FF"), add = "jitter",ylab="Variable Alelle Frequency (VAF) ", ylim=c(0,1),xlab="")+
theme_light()+theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"),axis.title.y = element_text(color="black", size=10, face="bold"))

pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/PPSintersect-color.pdf",width=7,height=3)
plot(p)
dev.off()

################################################################
################################################################
