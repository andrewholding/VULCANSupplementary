
######################
#
# Functions
#
######################



write.bed<-function(dba, filename){
    score <- -10*(log10(dba$FDR)) 
    write.table(cbind(dba[,1:3],rownames(dba),score),
                filename, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE)
}


getTFs <-function(vobj,threshold=NA){
    if (is.na(threshold)){
        TFs<-names(vobj$mrs[,"NES"]) 
    }else if (threshold>0) {
        TFs<-names(vobj$mrs[,"NES"][(vobj$mrs[,"NES"]>p2z(threshold))]) 
    }else if (threshold<0) {
        TFs<-names(vobj$mrs[,"NES"][(vobj$mrs[,"NES"]<sign(threshold)*p2z(-threshold))])
    }
    if (length(TFs) == 0)
    {   #plot Function doesn't like no labels.
        TFs<-names(vobj$mrs[,"NES"]) 
    }
    return(TFs)
}

vulcan.import.dba<- function (dbaobj, samples,intervals = NULL)  
{
    dbcounts <- dbaobj
    listcounts <- dbcounts$peaks
    names(listcounts) <- dbcounts$samples[, 1]
    first <- listcounts[[1]]
    rawmat <- matrix(NA, nrow = nrow(first), ncol = length(listcounts) + 
                         3)
    colnames(rawmat) <- c("Chr", "Start", "End", names(listcounts))
    rownames(rawmat) <- 1:nrow(rawmat)
    rawmat <- as.data.frame(rawmat)
    rawmat[, 1] <- as.character(first[, 1])
    rawmat[, 2] <- as.integer(first[, 2])
    rawmat[, 3] <- as.integer(first[, 3])
    for (i in 1:length(listcounts)) {
        rawmat[, names(listcounts)[i]] <- as.numeric(listcounts[[i]]$RPKM)
    }
    peakrpkms <- rawmat
    rm(rawmat)
    first <- listcounts[[1]]
    rawmat <- matrix(NA, nrow = nrow(first), ncol = length(listcounts) + 
                         3)
    colnames(rawmat) <- c("Chr", "Start", "End", names(listcounts))
    rownames(rawmat) <- 1:nrow(rawmat)
    rawmat <- as.data.frame(rawmat)
    rawmat[, 1] <- as.character(first[, 1])
    rawmat[, 2] <- as.integer(first[, 2])
    rawmat[, 3] <- as.integer(first[, 3])
    for (i in 1:length(listcounts)) {
        rawmat[, names(listcounts)[i]] <- as.integer(listcounts[[i]]$Reads)
    }
    peakcounts <- rawmat
    rm(rawmat)
    vobj <- list(peakcounts = peakcounts, samples = samples, 
                 peakrpkms = peakrpkms)
    return(vobj)
}

prependSampleNames <- function(vobj,prependString){
    colnames(vobj$peakcounts)[0:-3]<-paste0(prependString,colnames(vobj$peakcounts)[0:-3])
    vobj$samples[[2]]<-paste0(prependString,vobj$samples[[2]])
    vobj$samples[[1]]<-paste0(prependString,vobj$samples[[1]])
    colnames(vobj$peakrpkms)[0:-3]<-paste0(prependString,colnames(vobj$peakrpkms)[0:-3])
    return(vobj)
}


loadVulcanNetworks<-function(){
    regulons<-list()
    load("networks/laml-tf-regulon.rda")
    regulons$laml<-regul
    rm(regul)
    load("networks/brca-tf-regulon.rda")
    regulons$tcga<-regul
    rm(regul)
    load("networks/metabric-regulon-tfs.rda")
    regulons$metabric<-regulon
    rm(regulon)
    return(regulons)
}


plotVulcan <-function(vobj,threshold,network_title,title,plotTF,xlim,ylim) {
    threshold<-sign(threshold)*p2z(abs(threshold))
    network=vobj$mrs[,"NES"]
    tfs<-names(network)
    networkmat<-cbind(rep(0,length(network)),network[tfs])
    colnames(networkmat)<-c("0h","45min")
    matplot(t(networkmat),type="l",col="grey",ylab="VULCAN NES",xaxt="n",lty=1,main=title,xlim=xlim,ylim=ylim)
    axis(1,at=c(1:2),labels=colnames(networkmat))
    abline(h=c(0,threshold,-threshold),lty=2)
    text(2,networkmat[plotTF,2],label=names(networkmat[,2][plotTF]),pos=4,cex=0.6,font=2,col="red3")
    mtext(network_title)    
}    

######################
#
# Code
#
######################

library(DiffBind)
setwd("/Users/holdin01/Dropbox (Cambridge University)/Work/CRUK/Manuscripts/Vulcan/grhl2")


if(!file.exists("rdata/003_diffbind.rda")) {
    GRHL2 <- dba(sampleSheet="samplesheet/samplesheet.csv")
    plot(GRHL2)
    GRHL2 <- dba.count(GRHL2, summits=250)
    GRHL2 <- dba.contrast(GRHL2)
    GRHL2<- dba.analyze(GRHL2)
    save(GRHL2,file="rdata/003_diffbind.rda")
} else {
    load("rdata/003_diffbind.rda")  
}
#}

plot(GRHL2, contrast=1)
dba.plotMA(GRHL2, bFlip=1)
dba.plotPCA(GRHL2,components = 2:3)
dba.report(GRHL2)

r<-dba.report(GRHL2,th=1)
length(r[r$Fold>0])
length(r[r$Fold<0])

png("plots/venn.png", width=1000, height=500)
par(mfrow=c(1,2))
dba.plotVenn(GRHL2,GRHL2$masks$ER, label1="Rep1", label2="Rep2", main="GRHL2 +E2")
dba.plotVenn(GRHL2,GRHL2$masks$none, label1="Rep1", label2="Rep2", main="GRHL2 -E2")
dev.off()

called_none<-rowSums(GRHL2$called[,c(1:3)])
called_ER<-rowSums(GRHL2$called[,c(4:6)])
length(called_none[called_none>0])
length(called_ER[called_ER>0])
called_both<-called_none*called_none
length(called_both[called_both>0])


png("plots/dba.png")
dba.plotMA(GRHL2,bFlip=1)
dev.off()

#Code currently does not annotate correctly.
library(ChIPpeakAnno)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

rep=dba.report(GRHL2)
ucsc.hg38.knownGene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
chip.anno <- annotatePeakInBatch(rep, 
                                 AnnotationData=ucsc.hg38.knownGene)
chip.anno <- addGeneIDs(annotatedPeak=chip.anno, 
                       orgAnn="org.Hs.eg.db", 
                       feature_id_type="entrez_id",
                       IDs2Add="symbol")
write.table(as.data.frame(chip.anno), file="txt/report.txt")
library(vulcan)


dist_calc<-function(method,dfanno,genematrix,genesmore,allsamples){
    # This function structure was strongly suggested
    # by the Bioconductor reviewer
    supportedMethods<-c(
        "closest",
        "strongest",
        "sum",
        "topvar",
        "farthest",
        "lowvar"
    )
    if(!method%in%supportedMethods){
        stop("unsupported method ", method)
    }
    
    # Method closest: when multiple peaks are
    # found, keep only the closest to the TSS
    # as the representative one
    
    # Method farthest: when multiple peaks
    # are found, keep only the closest to the
    # TSS as the representative one
    
    # Method sum: when multiple peaks are
    # found, sum their contributions
    
    # Method strongest: when multiple peaks
    # are found, keep the strongest as the
    # representative one
    
    # Method topvar: when multiple peaks are
    # found, keep the most varying as the
    # representative one
    
    # Method lowvar: when multiple peaks are
    # found, keep the least varying as the
    # representative one
    
    
    for (gene in genesmore) {
        subanno <- dfanno[dfanno$feature == gene, ]
        
        if (method == "closest") {
            closest <- which.min(subanno$distanceToStart)
            genematrix[gene, allsamples] <- as.numeric(subanno[closest,
                                                               allsamples])
        }
        
        if (method == "farthest") {
            farthest <- which.max(subanno$distanceToStart)
            genematrix[gene, allsamples] <- as.numeric(subanno[farthest,
                                                               allsamples])
        }
        
        if (method == "sum") {
            sums <- apply(subanno[, allsamples], 2, sum)
            genematrix[gene, allsamples] <- as.numeric(sums)
        }
        
        if (method == "strongest") {
            sums <- apply(subanno[, allsamples], 1, sum)
            top <- which.max(sums)
            genematrix[gene, allsamples] <- as.numeric(subanno[top,
                                                               allsamples])
        }
        
        if (method == "topvar") {
            vars <- apply(subanno[, allsamples], 1, var)
            top <- which.max(vars)
            genematrix[gene, allsamples] <- as.numeric(subanno[top,
                                                               allsamples])
        }
        
        if (method == "lowvar") {
            vars <- apply(subanno[, allsamples], 1, var)
            top <- which.min(vars)
            genematrix[gene, allsamples] <- as.numeric(subanno[top,
                                                               allsamples])
        }
        
    }
    return(genematrix)
}


#SLOW vobj<-vulcan.import("samplesheet/samplesheet.csv")
#load(file="003_vobj.Rda")
samples <- list()
samples[['ER']]<-c('1a','2a','3a')
samples[['none']]<-c('1b','2b','3b')

vobj <-vulcan.import.dba(GRHL2,samples)

#vobj<-vulcan.annotate(vobj,lborder=-10000,rborder=10000,method='sum')

vobj <-prependSampleNames(vobj,"X")

lborder=-10000
rborder=10000
method='sum'
#DEBUG
#source("https://bioconductor.org/biocLite.R")
#biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

annotation <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene, 
                        feature = "gene")
gr <- GRanges(vobj$peakcounts)
seqlevels(annotation)<- sub('chr','',seqlevels(annotation))
anno <- annotatePeakInBatch(gr, AnnotationData = annotation, 
                            output = "overlapping", FeatureLocForDistance = "TSS", 
                            bindingRegion = c(lborder, rborder))
dfanno <- anno
names(dfanno) <- seq_len(length(dfanno))
dfanno <- as.data.frame(dfanno)
allsamples <- unique(unlist(vobj$samples))
genes <- unique(dfanno$feature)
peakspergene <- table(dfanno$feature)
rawcounts <- matrix(NA, nrow = length(genes), ncol = length(allsamples))
colnames(rawcounts) <- allsamples
rownames(rawcounts) <- genes
genesone <- names(peakspergene)[peakspergene == 1]
for (gene in genesone) {
    rawcounts[gene, allsamples] <- as.numeric(dfanno[dfanno$feature == 
                                                         gene, allsamples])
}

genesmore <- names(peakspergene)[peakspergene > 1]
rawcounts <- dist_calc(method, dfanno, rawcounts, genesmore, 
                       allsamples)

gr <- GRanges(vobj$peakrpkms)
anno <- annotatePeakInBatch(gr, AnnotationData = annotation, 
                            output = "overlapping", FeatureLocForDistance = "TSS", 
                            bindingRegion = c(lborder, rborder))
dfanno <- anno
names(dfanno) <- seq_len(length(dfanno))
dfanno <- as.data.frame(dfanno)
allsamples <- unique(unlist(vobj$samples))
genes <- unique(dfanno$feature)
peakspergene <- table(dfanno$feature)
rpkms <- matrix(NA, nrow = length(genes), ncol = length(allsamples))
rownames(rpkms) <- genes
genesone <- names(peakspergene)[peakspergene == 1]
colnames(rpkms)<-allsamples
for (gene in genesone) {
    rpkms[gene, allsamples] <- as.numeric(dfanno[dfanno$feature == 
                                                     gene, allsamples])
}
genesmore <- names(peakspergene)[peakspergene > 1]

rpkms <- dist_calc(method, dfanno, rpkms, genesmore, allsamples)
rawcounts <- matrix(as.numeric(rawcounts), nrow = nrow(rawcounts), 
                    dimnames = dimnames(rawcounts))
rpkms <- matrix(as.numeric(rpkms), nrow = nrow(rpkms), dimnames = dimnames(rpkms))
vobj$rawcounts <- rawcounts
colnames(rpkms)<-allsamples
vobj$rpkms <- rpkms

#DEBUG ENDS

vobj<-vulcan.normalize(vobj)

load(file="networks/metabric-regulon-tfs.rda")

regulons <- loadVulcanNetworks()
library("org.Hs.eg.db")
list_eg2symbol<-as.list(org.Hs.egSYMBOL[mappedkeys(org.Hs.egSYMBOL)]) 

vobj_results<-list()

test<-vulcan(vobj,       network=regulon, contrast=c("none","ER")    )

networks<-c("tcga","metabric") #,"laml")

for (network in networks) {
    vobj_results[[network]]<-vulcan(vobj,
                                       network=regulons[[network]],
                                       contrast=c("none","ER"),
                                       annotation=list_eg2symbol)
}


vobj_objects<-list(vobj_results)
names(vobj_objects)<-c("GRHL2")
networks<-c("tcga","metabric")

png("plots/vulcan_timelines.png", width=2000, height=1000, pointsize=20)
par(mfrow=c(2,3))
for(network in networks){
    for (vobj_name in names(vobj_objects))
    {

        vobj<-vobj_objects[[vobj_name]]
        vobj<-vobj[[network]]
        #TFs<-getTFs(vobj)
        TFs=c("ESR1","PGR","FOXA1","GRHL2") #Overide TFS
        plotVulcan(vobj,0.05, paste0(network," Network"),vobj_name,TFs,xlim=c(1,2.3), ylim=c(min(vobj$mrs[,"NES"]),max(vobj$mrs[,"NES"])))
        
        #TFs<-getTFs(vobj,0.05)
        plotVulcan(vobj,0.05, paste0(network," Network"),vobj_name,TFs,xlim=c(1,2.3), ylim=c(0,max(vobj$mrs[,"NES"])))
        
        #TFs<-getTFs(vobj,-0.05)
        plotVulcan(vobj,-0.05, paste0(network," Network"),vobj_name,TFs,xlim=c(1,2.3), ylim=c(min(vobj$mrs[,"NES"]),0))
    }
}
dev.off()


png("plots/vulcan_scatter.png", width=800, height=800, pointsize=15)
interect<-intersect(rownames(vobj_results$metabric$mrs), 
                       rownames(vobj_results$tcga$mrs))
plot(-log10(vobj_results$metabric$mrs[interect,"pvalue"]),
     -log10(vobj_results$tcga$mrs[interect,"pvalue"]),
     pch=20, xlab="Metabric Network Enrichment Score", ylab="TCGA Network Enrichment Score",main="VULCAN analysis of GRHL2 ChIP-Seq")

filtered<-(log10(vobj_results$metabric$mrs[interect,"pvalue"])^2+log10(vobj_results$tcga$mrs[interect,"pvalue"])^2)>2

text(-log10(vobj_results$metabric$mrs[interect,"pvalue"][filtered]),
     -log10(vobj_results$tcga$mrs[interect,"pvalue"][filtered])-0.1,
     labels=interect[filtered])

points(-log10(vobj_results$metabric$mrs[interect,"pvalue"][c("GATA3","PGR","FOXA1","ESR1","GRHL2")]),
       -log10(vobj_results$tcga$mrs[interect,"pvalue"][c("GATA3","PGR","FOXA1","ESR1","GRHL2")]),pch=20,col="red")
text(-log10(vobj_results$metabric$mrs[interect,"pvalue"][c("GATA3","PGR","FOXA1","ESR1","GRHL2")]),
     -log10(vobj_results$tcga$mrs[interect,"pvalue"][c("GATA3","PGR","FOXA1","ESR1","GRHL2")])-0.1,col="red",
     labels=c("GATA3","PGR","FOXA1","ESR1","GRHL2"))
dev.off()

#Note homer needs this as Hg format so you need to change from "1" to "chr1" etc.
df<-as.data.frame(dba.report(GRHL2))
df$seqnames<-paste0("chr",df$seqnames)
write.bed(df,"bed/up.bed")

df_all<-as.data.frame(dba.report(GRHL2, th=1))
df_all$seqnames<-paste0("chr",df_all$seqnames)
write.bed(df_all,"bed/all.bed")

