library("ggpubr")
eRNA<-read.csv("txt/eRNA_values.csv")

oe<-eRNA[eRNA$Experiment=="overexpression",]

oe_comparisons <- list( c("Control","GRHL2") )
oe_p <- ggboxplot(oe, x = "Treatment", y = "eRNA",
               add = "jitter", color="Treatment",
               facet.by="Gene", ylab="Relative eRNA levels")
oe_p <- oe_p + stat_compare_means(comparisons = oe_comparisons,
                       method = "wilcox.test" ,
                       paired=TRUE)
oe_p

si<-eRNA[eRNA$Experiment=="knockdown"& eRNA$Tissue=="MCF7",]
si_comparisons <- list( c("siCtrl","siGRHL2") )
MCF7p <- ggboxplot(si, x = "Treatment", y = "eRNA",
               add = "jitter", color="Treatment",
                facet.by="Gene",ylab="Relative eRNA levels")
MCF7p<-MCF7p + stat_compare_means(comparisons = si_comparisons,
                       method = "t.test" ,label.y=0.0014,
                       paired = TRUE)+theme(axis.title.x=element_blank())

si<-eRNA[eRNA$Experiment=="knockdown"& eRNA$Tissue=="T47D",]
si_comparisons <- list( c("siCtrl","siGRHL2") )
T47Dp <- ggboxplot(si, x = "Treatment", y = "eRNA",
               add = "jitter", color="Treatment",
               facet.by="Gene",ylab="Relative eRNA levels")
T47Dp<- T47Dp + stat_compare_means(comparisons = si_comparisons,
                          method = "t.test" ,
                          paired = TRUE,label.y=0.0005)+ guides(fill=FALSE, color=FALSE)+theme(axis.title.x=element_blank())

si<-eRNA[eRNA$Experiment=="knockdown"& eRNA$Tissue=="ZR75",]
si_comparisons <- list( c("siCtrl","siGRHL2") )
ZR75p <- ggboxplot(si, x = "Treatment", y = "eRNA",
                   add = "jitter", color="Treatment",
                   facet.by="Gene",ylab="Relative eRNA levels")
ZR75p<-ZR75p + stat_compare_means(comparisons = si_comparisons,
                              method = "t.test" ,
                              paired = TRUE,label.y=0.00025)+ guides(fill=FALSE, color=FALSE)
si_p<-ggarrange(MCF7p,T47Dp,ZR75p,nrow=3,labels=c("MCF7","T47D","ZR75"))
si_p

#Figure for manuscript
png("plot/eRNA.png",point=15)
si_p
dev.off()


wilcox.test(si$eRNA[si$Treatment=='siCtrl'],si$eRNA[si$Treatment=='siGRHL2'],paired=TRUE, alternative="less")
