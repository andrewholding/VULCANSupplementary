
library("ggpubr")
eRNA<-read.csv("txt/oe.csv")

eRNA<-eRNA[1:135,]

p <- ggboxplot(eRNA, x = "Condition", y = "Fold",
                   add = "jitter", color="Condition",
                   facet.by=c("Gene","CellTye"),ylab="Relative eRNA levels to empty vector")
comp.ev<-list(c("delta","OE"),c("delta","EV"),c("EV","OE"))
p2<-p + stat_compare_means(comparisons=comp.ev,
                           method = "t.test",paired=TRUE,alternative="greater",label.y=c(1.25,1.35,1.45))
p2

pdf("plot/oe_delta.pdf")
p2
dev.off()
