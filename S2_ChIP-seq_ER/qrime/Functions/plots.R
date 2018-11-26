library(ggplot2)
library(dplyr)

# intensity distribution plot

# intensities is a data frame containing columns for each sample
# sampleColours is a named vector that maps samples to colours

intensityDistributionPlot <- function(
  intensities,
  samples,
  sampleColours,
  title = "",
  xlab = "intensity",
  minIntensity = NA,
  maxIntensity = NA,
  showLegend = TRUE
)
{
  if (is.null(samples) || length(samples) == 0) return(ggplot() + theme_bw())

  # subset of samples within the given intensity data frame
  samples <- samples[which(samples %in% colnames(intensities))]
  if (length(samples) == 0) return(ggplot() + theme_bw())

  colours <- sampleColours[samples] %>% as.character

  intensities <- intensities %>%
    select(one_of(samples)) %>%
    gather(sample, intensity) %>%
    filter(complete.cases(.))

  intensities$sample <- factor(intensities$sample, levels = samples)

  plot <- ggplot(intensities, aes(x = intensity, colour = sample))
  plot <- plot + stat_density(geom = "line", position = "identity")
  plot <- plot + scale_colour_manual(values = colours, drop = FALSE)
  plot <- plot + ggtitle(title)
  plot <- plot + xlab(xlab)
  plot <- plot + theme_bw()
  plot <- plot + theme(
    plot.title = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.text = element_text(size = 6),
    legend.justification = c(1,1),
    legend.position = c(1,1),
    legend.background = element_rect(fill = "transparent")
  )

  if (!is.na(minIntensity) && !is.na(maxIntensity))
    plot <- plot + xlim(minIntensity, maxIntensity)

  if (!showLegend) plot <- plot + theme(legend.position = "none")

  return(plot)
}


# PCA plot

# groupColours is a named vector that maps groups to colours

pcaPlot <- function(intensities, samples, groups, groupColours, title = "", labels = NULL, legend = FALSE)
{
  if (is.null(samples) || length(samples) == 0) return(ggplot() + theme_bw())

  # subset of samples within the given intensity data frame
  indexes <- which(samples %in% colnames(intensities))
  if (length(indexes) == 0) return(ggplot() + theme_bw())
  samples <- samples[indexes]
  groups <- factor(groups[indexes])
  colours <- groupColours[levels(groups)] %>% as.character
  labels <- labels[indexes]

  intensities <- intensities %>%
    select(one_of(samples)) %>%
    filter(complete.cases(.))
  if (nrow(intensities) == 0) return(ggplot() + theme_bw())

  pca <- intensities %>% t %>% prcomp

  pcaVariance <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100)
  pca <- as.data.frame(pca$x)

  pca$Group <- groups

  plot <- ggplot(pca)
  plot <- plot + ggtitle(title)
  
  if (!is.null(labels))
    plot <- plot + geom_text(aes(x = PC1, y = PC2, label = paste("  ", labels)), hjust = 0, size = 2.5)
  plot <- plot + geom_point(aes(x = PC1, y = PC2, colour = Group), size = 3)
  plot <- plot + scale_colour_manual(values = colours)
  plot <- plot + scale_x_continuous(expand = c(0.2, 0))
  plot <- plot + scale_y_continuous(expand = c(0.2, 0))
  plot <- plot + xlab(paste("PC1, ", pcaVariance[1], "% variance", sep = ""))
  plot <- plot + ylab(paste("PC2, ", pcaVariance[2], "% variance", sep = ""))
  plot <- plot + theme_bw()
  plot <- plot + theme(text = element_text(size = 9))

  if (legend)
    plot <- plot + theme(
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.position = "right"
    )
  else
    plot <- plot + theme(legend.position = "none")

  return(plot)
}


# MA plot

maPlot <- function(differentialExpressionResults, selectedGenes = NULL,
                   xlab = expression(average~log[2](expression)),
                   ylab = expression(log[2](fold~change)),
                   significanceLevel = 0.05,
                   minLogFoldChangeForLabelling = 1.5,
                   maxNumberLabelledProteins = 50,
                   controlLogFoldChangeThreshold = -Inf,
                   pointSize = 2.5)
{
  if (!"logFCcontrol" %in% colnames(differentialExpressionResults))
    differentialExpressionResults$logFCcontrol <- Inf

  differentialExpressionResults$group <-
    ifelse(differentialExpressionResults$Gene %in% selectedGenes, 0,
    ifelse(!is.na(differentialExpressionResults$adj.P.Val) & differentialExpressionResults$adj.P.Val <= significanceLevel & differentialExpressionResults$logFCcontrol >= controlLogFoldChangeThreshold, 1,
    ifelse(!is.na(differentialExpressionResults$adj.P.Val) & differentialExpressionResults$adj.P.Val <= significanceLevel, 2,
    ifelse(abs(differentialExpressionResults$logFC) >= minLogFoldChangeForLabelling & differentialExpressionResults$logFCcontrol >= controlLogFoldChangeThreshold, 3,
    ifelse(abs(differentialExpressionResults$logFC) >= minLogFoldChangeForLabelling, 4,
    ifelse(differentialExpressionResults$logFCcontrol >= controlLogFoldChangeThreshold, 5, 6))))))

  # 0 - selected
  # 1 - significant, specific
  # 2 - significant, non-specific
  # 3 - non-significant, large logfc, specific
  # 4 - non-significant, large logfc, non-specific
  # 5 - specific
  # 6 - non-specific

  differentialExpressionResults <- differentialExpressionResults %>% arrange(desc(group))

  differentialExpressionResults$group <- factor(differentialExpressionResults$group, 0:6)

  differentialExpressionResults <- differentialExpressionResults %>%
    mutate(label = paste(" ", Gene)) %>%
    mutate(logFCcontrolAlpha = pmin(logFCcontrol, 2.0))

  labelSubset <- bind_rows(
    differentialExpressionResults %>%
      filter(group != 0 & adj.P.Val <= significanceLevel) %>%
      arrange(desc(abs(logFC))) %>%
      slice(0:(maxNumberLabelledProteins / 2)),
    differentialExpressionResults %>%
      filter(group != 0 & abs(logFC) >= minLogFoldChangeForLabelling) %>%
      arrange(desc(abs(logFC))) %>%
      slice(0:maxNumberLabelledProteins)
  ) %>%
    unique %>%
    slice(0:maxNumberLabelledProteins)

  plot <- ggplot(differentialExpressionResults, aes(x = AveExpr, y = logFC, colour = group, size = group, shape = group, alpha = logFCcontrolAlpha))
  plot <- plot + geom_hline(yintercept = 0, color="gray50", size = 0.5)
  plot <- plot + geom_point()
  plot <- plot + geom_text(data = labelSubset, aes(label = label), hjust = 0, size = 3.5)
  plot <- plot + geom_point(data = subset(differentialExpressionResults, group == 0), alpha = 1.0)
  plot <- plot + geom_text(data = subset(differentialExpressionResults, group == 0), aes(label = label), hjust = 0, size = 4.0, alpha = 1.0)
  plot <- plot + scale_colour_manual(values = c("blue", "deeppink3", "deeppink3", "gray50", "gray50", "gray50", "gray50"), drop = FALSE)
  plot <- plot + scale_size_manual(values = c(1.2 * pointSize, pointSize, pointSize, 0.8 * pointSize, 0.8 * pointSize, 0.6 * pointSize, 0.6 * pointSize), drop = FALSE)
  plot <- plot + scale_shape_manual(values = c(16, 16, 1, 16, 1, 16, 1), drop = FALSE)
  plot <- plot + scale_alpha(range = c(0.2, 1.0))
  plot <- plot + xlab(xlab)
  plot <- plot + ylab(ylab)
  plot <- plot + theme_bw()
  plot <- plot + theme(
    text = element_text(size = 14),
    legend.position = "none"
  )

  return(plot)
}


# Volcano plot

volcanoPlot <- function(differentialExpressionResults, selectedGenes = NULL,
                   xlab = expression(log[2](fold~change)),
                   ylab = expression(-log[10](adjusted~p-value)),
                   significanceLevel = 0.05,
                   minLogFoldChange = NULL,
                   maxLogFoldChange = NULL,
                   minLogFoldChangeForLabelling = 1.5,
                   maxNumberLabelledProteins = 50,
                   controlLogFoldChangeThreshold = -Inf,
                   pointSize = 2.5)
{
  if (!"logFCcontrol" %in% colnames(differentialExpressionResults))
    differentialExpressionResults$logFCcontrol <- Inf

  differentialExpressionResults$group <-
    ifelse(differentialExpressionResults$Gene %in% selectedGenes, 0,
    ifelse(!is.na(differentialExpressionResults$adj.P.Val) & differentialExpressionResults$adj.P.Val <= significanceLevel & differentialExpressionResults$logFCcontrol >= controlLogFoldChangeThreshold, 1,
    ifelse(!is.na(differentialExpressionResults$adj.P.Val) & differentialExpressionResults$adj.P.Val <= significanceLevel, 2,
    ifelse(abs(differentialExpressionResults$logFC) >= minLogFoldChangeForLabelling & differentialExpressionResults$logFCcontrol >= controlLogFoldChangeThreshold, 3,
    ifelse(abs(differentialExpressionResults$logFC) >= minLogFoldChangeForLabelling, 4,
    ifelse(differentialExpressionResults$logFCcontrol >= controlLogFoldChangeThreshold, 5, 6))))))

  # 0 - selected
  # 1 - significant, specific
  # 2 - significant, non-specific
  # 3 - non-significant, large logfc, specific
  # 4 - non-significant, large logfc, non-specific
  # 5 - specific
  # 6 - non-specific

  differentialExpressionResults <- differentialExpressionResults %>% arrange(desc(group))

  differentialExpressionResults$group <- factor(differentialExpressionResults$group, 0:6)

  differentialExpressionResults <- differentialExpressionResults %>%
    mutate(label = paste(" ", Gene)) %>%
    mutate(logFCcontrolAlpha = pmin(logFCcontrol, 2.0))

  differentialExpressionResults <- differentialExpressionResults %>%
    filter(!is.na(adj.P.Val)) %>%
    mutate(minusLog10AdjustedPValue = -log10(adj.P.Val))

  labelSubset <- bind_rows(
    differentialExpressionResults %>%
      filter(group != 0 & adj.P.Val <= significanceLevel) %>%
      arrange(desc(abs(logFC))) %>%
      slice(1:(maxNumberLabelledProteins / 2)),
    differentialExpressionResults %>%
      filter(group != 0 & abs(logFC) >= minLogFoldChangeForLabelling) %>%
      arrange(desc(abs(logFC))) %>%
      slice(0:maxNumberLabelledProteins)
  ) %>%
    unique %>%
    slice(0:maxNumberLabelledProteins)

  plot <- ggplot(differentialExpressionResults, aes(x = logFC, y = minusLog10AdjustedPValue, colour = group, size = group, shape = group, alpha = logFCcontrolAlpha))
  # plot <- plot + geom_hline(yintercept = 0, color="gray50", size = 0.5)
  plot <- plot + geom_point()
  plot <- plot + geom_text(data = labelSubset, aes(label = label), hjust = 0, size = 3.5)
  plot <- plot + geom_point(data = subset(differentialExpressionResults, group == 0), alpha = 1.0)
  plot <- plot + geom_text(data = subset(differentialExpressionResults, group == 0), aes(label = label), hjust = 0, size = 4.0, alpha = 1.0)
  plot <- plot + scale_colour_manual(values = c("blue", "deeppink3", "deeppink2", "gray50", "gray50", "gray50", "gray50"), drop = FALSE)
  plot <- plot + scale_size_manual(values = c(1.2 * pointSize, pointSize, pointSize, 0.8 * pointSize, 0.8 * pointSize, 0.6 * pointSize, 0.6 * pointSize), drop = FALSE)
  plot <- plot + scale_shape_manual(values = c(16, 16, 1, 16, 1, 16, 1), drop = FALSE)
  plot <- plot + scale_alpha(range = c(0.2, 1.0))
  if (!is.null(maxLogFoldChange))
  {
    if (is.null(minLogFoldChange)) minLogFoldChange <- -abs(maxLogFoldChange)
    plot <- plot + xlim(minLogFoldChange, maxLogFoldChange)
  }
  plot <- plot + xlab(xlab)
  plot <- plot + ylab(ylab)
  plot <- plot + theme_bw()
  plot <- plot + theme(
    text = element_text(size = 14),
    legend.position = "none"
  )

  return(plot)
}


# Generates a histogram for the specified column within a data frame

histogram <- function(data, column,
                      xlab = "x",
                      ylab = "count")
{
  plot <- ggplot(data, aes_string(x = column))
  plot <- plot + geom_histogram(colour = "white", fill = "gray50")
  plot <- plot + xlab(xlab)
  plot <- plot + ylab(ylab)
  plot <- plot + theme_bw()
  return(plot)
}


# Generates a scatter plot for the specified columns within a data frame

scatterPlot <- function(data, xcolumn, ycolumn,
                        xlab = xcolumn,
                        ylab = ycolumn)
{
  plot <- ggplot(data, aes_string(x = xcolumn, y = ycolumn))
  plot <- plot + geom_point()
  plot <- plot + xlab(xlab)
  plot <- plot + ylab(ylab)
  plot <- plot + theme_bw()
  plot <- plot + theme(text = element_text(size = 8))
  return(plot)
} 
