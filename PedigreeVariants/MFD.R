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
# general analyis that might be interesting
################################################################
vcf <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/ALLVCFfiles/PPSMFDmerged.vcf", verbose = FALSE)
chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
chrom <- masker(chrom, min_DP = 300, max_DP = 700)
chrom <- proc.chromR(chrom, verbose = FALSE)
################
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
head(dp)
heatmap.bp(dp[1:179,])
is.na(dp[na.omit(dp == 0)]) <- TRUE
heatmap.bp(dp[1:179,])
ad <- extract.gt(chrom, element="AD", as.numeric=TRUE)
ad <- extract.gt(vcf, element = 'AD')

################
par(mar=c(8,4,4,2))
barplot(apply(dp, MARGIN=2, mean, na.rm=TRUE), las=3)
################
par(mar=c(8,4,1,1))
pdf("/Users/MoiMeamina/Documents/MDF/181011-Figures/MFDdepth.pdf",width=8,height=4)
boxplot(dp, col=c("#C0C0C0", "#808080"), ylab="Depth",log='y', las=2)
#abline(h=seq(0,1e4, by=100), col="#C0C0C088")
dev.off()


################################################################
################################################################
################################################################
#Already calculated and the files are created

M <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ERR194161.vcf", verbose = FALSE)
F <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ERR194160.vcf", verbose = FALSE)
D <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ERR194147.vcf", verbose = FALSE)

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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#write the datafrime into files
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
write.table(dframeF, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addF.txt",sep="\t",quote=T, na = "NA")
write.table(dframeM, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addM.txt",sep="\t",quote=T, na = "NA")
write.table(dframeD, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addD.txt",sep="\t",quote=T, na = "NA")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Avant de lire les fichier, j'ai modifier manuellement les ficier pour ajouter une colonne pos
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TF<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addF.txt", sep="\t", header=T)
TD<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addD.txt", sep="\t", header=T)
TM<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addM.txt", sep="\t", header=T)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Some statistical tests
wilcox.test(TM$X20.1, TF$X20.1)
names=c(rep("NA12891", nrow(TF)), rep("NA12892", nrow(TM)), rep("NA12878", nrow(TD)))
TTALL=data.frame(names, rbind(TF,TM,TD))
res.ftest <- bartlett.test(X20.1 ~ names, data = TTALL)

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
names=c(rep("NA12891", length(dF)), rep("NA12892", length(dM)), rep("NA12878", length(dd)))
alldfmd=c( dF, dM, dd)
diffonly=data.frame(names, c( dF, dM, dd))
#__________________________________________________________
#creat the file that contain the locationand the VAF of each difference.
names=c(rep("NA12891", length(dfall)), rep("NA12892", length(dmall)), rep("NA12878", length(ddall)))
alldf=cbind(rep("NA12891", length(dfall)), dfall)
diffonltall=data.frame(cbind(names, c(dfall,dmall,ddall)))
#write.table(diffonltall, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/alldifference.txt",sep="\t",quote=T, na = "NA")
#__________________________________________________________________________________________________________

################################################################
#intersect only present in all M, D, F #107 variants
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
#write.table(ddd, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/MFDintersectall.txt",sep="\t",quote=T, na = "NA")
#__________________________________________________________
length(intD)
length(intF)
length(intM)
names=c(rep("NA12891", length(intF)), rep("NA12892", length(intM)), rep("NA12878", length(intD)))
intersectonly=data.frame(names, c( intF, intM, intD))
#____________________________________________________________________________________
################################################################
################### Ploter avec des couleurs ###################
################################################################
#Avant de lire les fichier, j'ai modifier manuellement les ficier pour ajouter une colonne pos
TF<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addF.txt", sep="\t", header=T)
TD<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addD.txt", sep="\t", header=T)
TM<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addM.txt", sep="\t", header=T)
#colors jco : c("NA12877"="#0073C2FF","NA12878"="#EFC000FF","NA12889"="#868686FF","NA12890"="#CD534CFF","NA12891"="#7AA6DCFF","NA12892"="#003C67FF"))

#__________________________________________________________
#__________________________________________________________

names=c(rep("NA12891", nrow(TF)), rep("NA12892", nrow(TM)), rep("NA12878", nrow(TD)))
TTALL=data.frame(names, rbind(TF,TM,TD))
TTALL$names=factor(TTALL$names, levels=c("NA12878","NA12891","NA12892"))

p<-ggboxplot(TTALL, x = "names", y ="X20.1", shape="names", color="names",palette = c("#EFC000FF","#7AA6DCFF","#003C67FF"), add = "jitter",ylab="Variable Alelle Frequency (VAF) ", xlab="")+
theme_light()+theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"),axis.title.y = element_text(color="black", size=10, face="bold"))

pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/MFD-colors-VAF.pdf",width=7,height=3)
plot(p)
dev.off()
#__________________________________________________________
#__________________________________________________________
diffonly$names=factor(diffonly$names, levels=c("NA12878","NA12891","NA12892"))

p<-ggboxplot(diffonly, x = "names", y ="c.dF..dM..dd.", shape="names", color="names",palette = c("#EFC000FF","#7AA6DCFF","#003C67FF"), add = "jitter",ylab="Variable Alelle Frequency (VAF) ", ylim=c(0,1),xlab="")+
theme_light()+theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"),axis.title.y = element_text(color="black", size=10, face="bold"))

pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/MFDdiff-colored.pdf",width=7,height=3)
plot(p)
dev.off()
#__________________________________________________________
#__________________________________________________________

intersectonly$names=factor(intersectonly$names, levels=c("NA12878","NA12891","NA12892"))
p<-ggboxplot(intersectonly, x = "names", y ="c.intF..intM..intD.", shape="names", color="names",palette = c("#EFC000FF","#7AA6DCFF","#003C67FF"), add = "jitter",ylab="Variable Alelle Frequency (VAF) ", ylim=c(0,1),xlab="")+theme_light()+theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"),axis.title.y = element_text(color="black", size=10, face="bold"))
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/MFDintersect-color.pdf",width=7,height=3)
plot(p)
dev.off()

################################################################################
################################################################################

intersectable<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/MFDintersectall.txt", sep="\t", header=T)
#+++++++++++++++++++++
fit <- lm(intersectable$intF ~ intersectable$intM)
formula <- y ~ x
p4<-ggplot(intersectable, aes(x=intF, y=intM)) + geom_point()+geom_smooth(method=lm )+theme_bw()+
theme(axis.text.x=element_text(size=rel(1.25)),axis.text.y=element_text(size=rel(1.25)), text = element_text(size=14))+
labs( x="intersection Father",y = "intersection Mother")+theme(legend.position="top")+
geom_smooth(method = "lm", formula = formula, se = F) +
stat_poly_eq(aes(label = paste(..rr.label..)),label.x.npc = "right", label.y.npc = 0.9, formula = formula, parse = TRUE, size = 4)+
stat_fit_glance(method = 'lm',method.args = list(formula = formula),geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.85, size = 4)

################################################################################
################################################################################
#comparaison old and new reference
#plot in colors
#__________________________________________________________
##Veenplot
library(VennDiagram)
library(ggpubr)

DO<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/0-OLDREFERENCE/Results-FMD/addD.txt", sep="\t", header=T)
DN<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addD.txt", sep="\t", header=T)

pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/MFD-veen.pdf")
venn.plot <- venn.diagram(list(DO$pos, DN$pos),NULL, fill=c("gray23",  "gray43"), alpha=c(1,1), cex = 2, cat.fontface=4, category.names=c( "Daughter O","Daughter N" ) , na.rm = FALSE)
grid.draw(venn.plot)
dev.off()

