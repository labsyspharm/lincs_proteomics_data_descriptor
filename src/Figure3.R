library(UpSetR)
library(corrplot)
library(ggplotify)

# -- Figure 3a ---------
dfs <- read.csv(file='../data/upsetplot_input.csv')

f3a <- upset(dfs, sets=c('batch1', 'batch2', 'batch3', 'batch4', 'batch5'), 
          mainbar.y.label = "proteins identified in common", sets.x.label = "batch size", 
          text.scale = c(1.3, 1.3, 1, 1, 2, 0.75),
          order.by="freq", keep.order=TRUE)

# Figure 3b ------
data <- read.csv('../data/non_malignants_samples.csv')
x <- cor(data)
colnames(x) = c("MCF10A rep1", "MCF10A rep2", 'MCF10A rep3', 'hME1')
rownames(x) = c("MCF10A rep1", "MCF10A rep2", 'MCF10A rep3', 'hME1')
#f3b <- ggcorrplot(x, method = "circle", type='lower', lab=TRUE)
f3b <- corrplot.mixed(x, lower.col = "black", number.cex = .7, tl.cex=0.8)


#f3ab <- ggarrange(f3a, f3b, nrow=1, ncol=2, widths=c(0.7, 0.3))
