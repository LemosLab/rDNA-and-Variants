#!/usr/bin/env Rscript
##########################################################################################
#GNU GENERAL PUBLIC LICENSE
#Version 3, 29 June 2007

#Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
#Everyone is permitted to copy and distribute verbatim copies
#of this license document, but changing it is not allowed.

        #AUTHORS: Bedrat amina & Bernardo Lemos (2018)
                            #HARVARD.
#The GNU General Public License is a free, copyleft license for
#software and other kinds of works.
###########################################################################################

source("http://www.Bioconductor.org/biocLite.R")
#biocLite( c("ShortRead","DESeq", "edgeR") )

library( "DESeq2" )
library("RColorBrewer")
library("ggplot2") #Best plots
library("ggrepel") #Avoid overlapping labels
library("ggpubr")

#/////////////////
#Read input file
inputcondition="PATH/TO/TCGACSV4.csv"
inputMcount="PATH/TO/CGTAMergedCounts4.csv"
inputlengthgene="PATH/TO/genelengthMeng.txt"
gf = "PATH/TO/Homo_sapiens.GRCh37.75.gtf"
#/////////////////
outputcpm="PATH/TO/DEG_CPM.txt"
outputrpkm="PATH/TO/DEG_rpkm.txt"
PATHTOFIG="PATH/TO/Results/"
#/////////////////

#//////////////////////////////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////  Analysis with EdgrR            ////////////////////////////////////
#//////////////////////////////////////////////////////////////////////////////////////////////////////
sri = read.csv(inputcondition, stringsAsFactors=FALSE, sep=";")
sri = sri[order(sri$patients),]
MergedCount=read.table(inputMcount, header=T,sep=";",row.names = 1)
gene=read.table(inputlengthgene, header=T,sep="\t",row.names = 1)
#==================
dds = DESeqDataSetFromMatrix(countData = MergedCount , colData = sri , design =~groups+patients)
#filter
idx <- rowSums(fpm(dds) >= 2 ) >= 8
dds <- dds[idx,]
#_________________________________________________
##Differential expression analysis
#_________________________________________________
dds <- DESeq(dds)
#_________________________________________________
#ploting PCA
#_________________________________________________
rld = rlogTransformation(dds)
rld <- rlogTransformation(dds, blind=TRUE)
#change path to the results path
setwd(PATHTOFIG)
png("PCA.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
p <- plotPCA(rld,intgroup=c('groups'),ntop = 6000, returnData = F)+theme_light()
p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")+theme_light()
dev.off()
#_________________________________________________
#Results choose one of the contrast
#_________________________________________________
#res <- results(dds) #print DEGs for donors
##or
res <- results(dds,  contrast=c("groups","PTumor","STNormal"))
summary(res)
#print results
mcols(res, use.names=TRUE)
up=sum(res$padj < 0.05 & res$log2FoldChange >= 0.5 , na.rm=TRUE)
down=sum(res$padj < 0.05 & res$log2FoldChange <= -0.5, na.rm=TRUE)
up
down

respadjfc <- res[ order(res$padj, decreasing = FALSE), ]
sum(res$padj < 0.05, na.rm=TRUE)
head(respadjfc)

################################################################################
################################################################################
#write results and add the name of the genes
################################################################################
################################################################################

respadjfc <- res[order(res$padj, decreasing = FALSE), ]
data2=fpm(dds)
countsfile=respadjfc
Length=read.table("/Volumes/AMINABE/0-2018-RNASEQ/0-2018_Mike/outputs/edgR/counts/R-gene_name_length.txt", header=T,sep="\t",row.names = 1)

#////////////////////////////all this to just convert names ////////////////////////////////////////
outputfile3="BAmRNAcountsResults.txt"
outputfile2="BAMergedCountsmRNAgenename.csv"
outputfile1="BACountsmRNAgenename.csv"
length=data.frame(rownames(Length),Length$name, Length$length)

#table of results:
idx=match(length$rownames.Length,rownames(countsfile))
Lengthmatch = as.data.frame(length[!is.na(idx),])
tedgeRmatch = countsfile[idx[!is.na(idx)],]
tt=cbind(tedgeRmatch, Lengthmatch)
#write.table(tt[order(tt$padj),],file=outputfile2,sep=";",quote=F)

#table of counts
idx=match(length$rownames.Length,rownames(data2))
Lengthmatch = as.data.frame(length[!is.na(idx),])
tedgeRmatch = data2[idx[!is.na(idx)],]
tt2=cbind(tedgeRmatch, Lengthmatch)
#write.table(tt2[order(tt2$padj),],file=outputfile1,sep=";",quote=F)

#THE RESULT TABLE
idx=match(rownames(tt),rownames(tt2))
Lengthmatch = as.data.frame(tt[!is.na(idx),])
tedgeRmatch = tt2[idx[!is.na(idx)],]
tt3=cbind(tedgeRmatch, Lengthmatch)
write.table(tt3[order(tt3$padj),],file=outputfile3,sep="\t",quote=F)

################################################################################
################################################################################
#plot
################################################################################
################################################################################
#Ploting dispertionALL
#_________________________________________________
png("Dispersions.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
plotDispEsts(dds)
dev.off()
#_________________________________________________
#plot pval
png("Pval_Distribution.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
hist(res$pvalue,breaks=100,col="grey", xlab="p-value",main="p-value distribution")
dev.off()

########################################
name=read.table("BAmRNAcountsResults.txt", header=T,sep="\t",row.names = 1)
#MAplot
png("MAplot.png", width = 7*300, height = 5*300, res = 400, pointsize = 12)
ggmaplot(name, main = expression(""), #for mRNA
fdr = 0.05, fc = 0.75, size = 1,
palette = c("#B31B21", "#1465AC", "darkgray"),
genenames = as.vector(name[,"Length.name"]),#for mRNA
legend = "top", top = 30,
font.label = c("bold", 5),
font.legend = "bold",
font.main = "bold",
ggtheme = ggplot2::theme_light())
dev.off()
########################################
########################################
png("Volcanoplot.png", width = 7*400, height = 5*300, res = 400, pointsize = 12)
mydata<-name%>%mutate(threshold = ifelse(log2FoldChange >= 1 & padj<0.05 ,"UP", ifelse(log2FoldChange<=-1 & padj<0.05, "Down", "NA")))
mydata2=subset(mydata, mydata$threshold=="UP" | mydata$threshold=="Down")
ggplot(mydata, aes(x=log2FoldChange, y=-log10(padj)),fdr = 0.05, fc = 1, size = 100,) + xlim(-2, 2)+#ylim(0,10)+
geom_point(aes(colour = threshold), size=1.5) +theme_light()+
theme(axis.text.x=element_text(size=rel(1.5)),axis.text.y=element_text(size=rel(1.5)), text = element_text(size=14))+
scale_colour_manual(values = c("UP"= "#B31B21", "Down"="#1465AC",  "NA"= "darkgray"))+
geom_text_repel(data=head(mydata2, 40), aes(label=Length.name) ,size = 2)+ #adding text for the top 20 genes
geom_vline(xintercept=0, linetype="dashed",color = "black", size=.5)

dev.off()

#END
