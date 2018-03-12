qpcr<-read.csv("txt/QPCR GREB1 XBP1 TFF1.csv", header=TRUE)
df<-data.frame(qpcr)
library("ggpubr")
p <- ggboxplot(df, x = "Sample", y = "Expression", ylab="Relative Expression",
               xlab="Condition",
               color = "Sample", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = "Sample", facet.by="Target",
               outlier=FALSE)
p

my_comparisons <- list( c("Ctrl", "siGRHL2"), c("Ctrl","GRHL2-DDK"))
p + 
    stat_compare_means(comparisons = my_comparisons,method="t.test", paired=TRUE)  +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

