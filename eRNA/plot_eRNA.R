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

si<-eRNA[eRNA$Experiment=="knockdown",]

si_comparisons <- list( c("siCtrl","siGRHL2") )
p <- ggboxplot(si, x = "Treatment", y = "eRNA",
               add = "jitter", color="Treatment",
               facet.by="Gene", ylab="Relative eRNA levels")
p + stat_compare_means(comparisons = si_comparisons,
                       method = "wilcox.test" ,
                       paired = TRUE)

#Figure for manuscript
png("plot/eRNA.png",point=15)
oe_p
dev.off()


wilcox.test(si$eRNA[si$Treatment=='siCtrl'],si$eRNA[si$Treatment=='siGRHL2'],paired=TRUE, alternative="less")
