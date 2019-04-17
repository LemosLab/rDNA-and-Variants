################
library(vcfR)
library(ggpubr)
library(dplyr)
library(ggpmisc)
library(ggplot2)
library(VennDiagram)
library(venn)

################################################################
################################################################

vcf <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/ALLVCFfiles/ALL14samp.vcf", verbose = FALSE)
chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
chrom <- masker(chrom, min_QUAL = 30)
chrom <- proc.chromR(chrom, verbose = FALSE)

################
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE)
rownames(dp) <- 1:nrow(dp)
head(dp)
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/Figures/14samplesheatmap-VAF.pdf",width=12,height=20)
is.na(dp[na.omit(dp == 0)]) <- TRUE
heatmap.bp(dp[1:233,])
dev.off()
################
par(mar=c(8,4,4,2))
barplot(apply(dp, MARGIN=2, mean, na.rm=TRUE), las=3)
################
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/Figures/14sampesprogenydepth.pdf",width=8,height=5)
boxplot(dp, col=c("#C0C0C0", "#808080"), ylab="Depth",log='y', las=2)
dev.off()
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
    write.table(dframe, file=path,sep="\t",quote=T, na = "NA")
    return(dframe)
}
#AAAAAMMMMIIINNNNNAAAA do not run it will delet all the files already modified
dframe14S = VAFcalculation(chrom,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/14samples-test.txt")

################################################################
################################################################

vcf <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/PROGmerged.vcf", verbose = FALSE)
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
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/Figures/progenydepth.pdf",width=8,height=5)
boxplot(dp, col=c("#C0C0C0", "#808080"), ylab="Depth",log='y', las=2)
#abline(h=seq(0,1e4, by=100), col="#C0C0C088")
dev.off()
#ALL
dpvcf  <-  extract.gt(vcf, element="DP", as.numeric = TRUE)
advcf  <-  extract.gt(vcf, element = 'AD')
RAMvcf <- masplit(advcf, record = 1, sort = 0)
RAFvcf <- RAMvcf/dpvcf
VAFvcf <- 1-RAFvcf
dframevcf=data.frame(VAFvcf))
#write.table(dframevcf, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/ALLVAF.txt",sep="\t",quote=T, na = "NA")
TABLEVCF<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/ALLVAF.txt", sep="\t", header=T)

vaf=subset(TABLEVCF, TABLEVCF$Mean >0.2)
ggplot(TABLEVCF, aes(x=SRR4435254, y=SRR4435252)) + geom_point()+theme_bw()+
geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.6,position=position_dodge(0.05))

################################################################
################################################################
#the fuction works fin tested April 17, 2019
#Already calculated and the files are created

M <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/SRR4435252.vcf", verbose = FALSE)
F <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/SRR4435253.vcf", verbose = FALSE)
D <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/SRR4435254.vcf", verbose = FALSE)
PR   <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/SRR4435255.vcf", verbose = FALSE)
PRR  <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/SRR4435257.vcf", verbose = FALSE)
PPR  <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/SRR4435258.vcf", verbose = FALSE)
PPRR <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/SRR4435260.vcf", verbose = FALSE)
SS   <- read.vcfR("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/SRR4435266.vcf", verbose = FALSE)

cM <- create.chromR(name="Supercontig", vcf=M, verbose=FALSE)
cF <- create.chromR(name="Supercontig", vcf=F, verbose=FALSE)
cD <- create.chromR(name="Supercontig", vcf=D, verbose=FALSE)
C1 <- create.chromR(name="Supercontig", vcf=PR, verbose=FALSE)
C2 <- create.chromR(name="Supercontig", vcf=PRR, verbose=FALSE)
C3 <- create.chromR(name="Supercontig", vcf=PPR, verbose=FALSE)
C4 <- create.chromR(name="Supercontig", vcf=PPRR, verbose=FALSE)
C5 <- create.chromR(name="Supercontig", vcf=SS, verbose=FALSE)

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
dframeM = VAFcalculation(cM,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addMSRR4435252.txt")
dframeD = VAFcalculation(cD,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addDSRR4435254.txt")
dframeF = VAFcalculation(cF,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addFSRR4435253.txt")
dframeC1= VAFcalculation(C1,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC1SRR4435255.txt")
dframeC2= VAFcalculation(C2,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC1SRR4435257.txt")
dframeC3= VAFcalculation(C3,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC1SRR4435258.txt")
dframeC4= VAFcalculation(C4,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC1SRR4435260.txt")
dframeC5= VAFcalculation(C5,"/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC1SRR4435266.txt")
################################################################
#Avant de lire les fichier, j'ai modifier manuellement les ficier pour ajouter une colonne pos
################################################################
#Les 4 grandsparents and les 2 parents
TPF<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addf.txt", sep="\t", header=T)
TS<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addS.txt", sep="\t", header=T)
TPM<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/PPSVCFfiles/Results-PPS/addM.txt", sep="\t", header=T)
TF<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addF.txt", sep="\t", header=T)
TD<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addD.txt", sep="\t", header=T)
TM<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/MFDVCFfiles/ResultsMFD/addM.txt", sep="\t", header=T)

#la progeny
TTF<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addFSRR4435253.txt", sep="\t", header=T)
TTD<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addDSRR4435254.txt", sep="\t", header=T)
TTM<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addMSRR4435252.txt", sep="\t", header=T)
c1<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC1SRR4435255.txt", sep="\t", header=T)
c2<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC2SRR4435257.txt", sep="\t", header=T)
c3<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC3SRR4435258.txt", sep="\t", header=T)
c4<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC4SRR4435260.txt", sep="\t", header=T)
c5<- read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/progenyfiles/addC5SRR4435266.txt", sep="\t", header=T)

#read all and creat one table
names=c(rep("NA12891", nrow(TF)), rep("NA12892", nrow(TM)), rep("NA12878", nrow(TD)),rep("NA12889", nrow(TPF)), rep("NA12890", nrow(TPM)), rep("NA12877", nrow(TS)),rep("NA12880", nrow(TTF)), rep("NA12879", nrow(TTM)), rep("NA12881", nrow(TTD)),rep("NA12882", nrow(c1)),rep("NA12884", nrow(c2)),rep("NA12885", nrow(c3)),rep("NA12887", nrow(c4)),rep("NA12893", nrow(c5)))
TTALL=data.frame(names, rbind(TF,TM,TD, TPF, TPM,TS,TTF,TTM,TTD,c1,c2,c3,c4,c5))
#test
pairwise.wilcox.test(TTALL$"X20.1", TTALL$"names",p.adjust.method = "BH")

#__________________________________________________________
#__________________________________________________________
#Correlation Age_association_coefficient [top to 0.25]
fit <- lm(TTALL$pos ~ TTALL$X20.1)
formula <- y ~ x
p1<-ggplot(TTALL, aes(x=pos, y=X20.1, fill=names, color=names)) + geom_point()+geom_smooth(method=lm )+theme_bw()+
theme(axis.text.x=element_text(size=rel(1.25)),axis.text.y=element_text(size=rel(1.25)), text = element_text(size=14))+
labs(title="promoter")+theme(legend.position="top")+geom_smooth(method = "lm", formula = formula, se = F) +
stat_poly_eq(aes(label = paste(..rr.label..)),label.x.npc = "right", label.y.npc = 0.9, formula = formula, parse = TRUE, size = 4)+
stat_fit_glance(method = 'lm',method.args = list(formula = formula),geom = 'text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 0.85, size = 4)

#__________________________________________________________
#__________________________________________________________

my_comparisons <- list( c("NA12878", "NA12891"), c("NA12891", "NA12892"), c("NA12878", "NA12892"))
p<-ggboxplot(TTALL, x = "names", y ="X20.1", shape="names", fill="names",add = "jitter",ylab="Variable Alelle Frequency (VAF) ",xlab="")+scale_shape_manual(values=1:nlevels(TTALL$names)) +
stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test", aes(label=..p.adj..))+stat_compare_means(label.y = 1.4)+
theme_light()+
theme(axis.text.x = element_text(face="bold", color="black",size=9),axis.text.y = element_text(face="bold", color="black",size=9),axis.title.x = element_text(color="black", size=10, face="bold"),axis.title.y = element_text(color="black", size=10, face="bold"))

pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/Figures/progeny-VAF.pdf",width=12,height=4)
plot(p)
dev.off()
#__________________________________________________________
##Veenplot
#__________________________________________________________
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/Figures/all-veen.pdf")
ALL=list(TTF$pos, TTD$pos, TTM$pos,c1$pos, c2$pos, c3$pos,c4$pos,c5$pos)
names(ALL) <- c("NA12880","NA12881","NA12879","NA12882","NA12884","NA12885","NA12887","NA12893")
cc<-venn(ALL,ilab=TRUE, zcolor = "style",size = 19, cexil = 1.3, cexsn = 1.25, borders = F)
dev.off()

################################################################
#different
################################################################
x=c1$pos
y=TD$pos
z=TS$pos
dataf=dframeC1

getdiff <- function( x,y,z,dataf) {
    dd=c()
    ddall=c()
    cc=setdiff(setdiff(x, y), z)
    for (i in 1:length (cc)){
        ccc=paste("KY962518.1_",cc[i],sep = "",collapse = ' ')
        vv=cbind(ccc,dataf[ccc,"X20.1"] )
        ddall=rbind(ddall, vv)
        print(dataf[ccc,"X20.1"])
        dd=c(dd,dframeC1[ccc,"X20.1"])
    }
    return(list(val1=ddall, val2=dd))
}
#how to run
differences=getdiff(c1$pos,TD$pos,TS$pos, dframeC1)
ddall=differences$val1
dF=differences$val2

#il faut creer les 3 list dF, dM, dd avant de poursuivre
#__________________________________________________________
names=c(rep("NA12891", length(dF)), rep("NA12892", length(dM)), rep("NA12878", length(dd)))
alldfmd=c( dF, dM, dd)
diffonly=data.frame(names, c( dF, dM, dd))
#__________________________________________________________
#creat the fale that contain the locationand the VAF of each difference.
names=c(rep("NA12891", length(dfall)), rep("NA12892", length(dmall)), rep("NA12878", length(ddall)))
alldf=cbind(rep("NA12891", length(dfall)), dfall)
diffonltall=data.frame(cbind(names, c(dfall,dmall,ddall)))
#write.table(diffonltall, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/MFDVCFfiles/ResultsMFD/alldifference.txt",sep="\t",quote=T, na = "NA")

#__________________________________________________________________________________________________________
#______________________plot difference with the comparison p value_________________________________________

library(ggplot2)
my_comparisons <- list( c("NA12878", "NA12891"), c("NA12891", "NA12892"), c("NA12878", "NA12892"))
p<-ggboxplot(diffonly, x = "names", y ="c.dF..dM..dd.", add = "jitter",shape = "names", label.y = 2, ylab="Variable Alelle Frequency (VAF) ")+
#ggtitle("difference between M.(31), F.(10) adn D.(17)")+
stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test", aes(label=..p.adj..) )+
stat_compare_means(label.y = 0.5)
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/Figures/progeny-diff.pdf",width=8,height=4)
plot(p)
dev.off()

#__________________________________________________________
#__________________________________________________________

pairwise.wilcox.test(diffonly$"c.dF..dM..dd.", diffonly$"names",p.adjust.method = "BH")

################################################################
#intersect only present in all M, D, F #107 variants
################################################################

getintersect <- function( x,y,z,dataS,dataM, dataF) {
    intD=c()
    intF=c()
    intM=c()
    cc=intersect(intersect(x,y), intersect(x,z))
    for (i in 1:length (cc)){
        ccc=paste("KY962518.1_",cc[i],sep = "",collapse = ' ')
        vv=c(vv,ccc)
        intD=c(intD,dataS[ccc,"X20.1"])
        intF=c(intF,dataF[ccc,"X20.1"])
        intM=c(intM,dataM[ccc,"X20.1"])
        
    }
    return(list(val0=vv, val1=intD, val2=intF, val3=intM))
}

#__________________________________________________________
#creat the fale that contain the locationand the VAF of each difference.
ddd=cbind (vv, intD,intF,intM)
#write.table(ddd, file="/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/MFDVCFfiles/ResultsMFD/MFDintersectall.txt",sep="\t",quote=T, na = "NA")
#__________________________________________________________
length(intD)
length(intF)
length(intM)

names=c(rep("NA12891", length(intF)), rep("NA12892", length(intM)), rep("NA12878", length(intD)))
intersectonly=data.frame(names, c( intF, intM, intD))
library(ggplot2)

#____________________________________________________________________________________
#__________________plot intersection with the comparison p value_____________________

my_comparisons <- list( c("NA12878", "NA12891"), c("NA12891", "NA12892"), c("NA12878", "NA12892"))
p<-ggboxplot(intersectonly, x = "names", y ="c.intF..intM..intD.", shape = "names",add = "jitter",ylab="Variable Alelle Frequency (VAF) ")+
#ggtitle("Intersect between M.(218), F.(218) adn D.(218)")+
stat_compare_means(comparisons = my_comparisons,  method = "wilcox.test", aes(label=..p.adj..) )+
stat_compare_means(label.y = 1.6)
pdf("/Volumes/AMINABE/0-2018-RNASEQ/0-MFDPPS/1-NEWREFERENCE/Progeny-vcffiles/Figures/progeny-intersect.pdf",width=8,height=4)
plot(p)
dev.off()

#__________________________________________________________
pairwise.wilcox.test(intersectonly$"c.intF..intM..intD.", intersectonly$"names",p.adjust.method = "BH")


