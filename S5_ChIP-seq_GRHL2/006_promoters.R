require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
p<-promoters(genes(txdb), upstream = 1500, downstream = 500)

write.bed<-function(df, filename){
    write.table(as.data.frame(df)[,1:3],
                filename, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE)
}

write.bed(p,"bed/promoter.bed")
