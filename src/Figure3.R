library(UpSetR)
library(corrplot)

# -- Figure 3a ---------
dfs <- read.csv(file='../data/upset_plot_input_8batches.csv')
cols <- c('Batch 1', 'Batch 2', 'Batch 3', 'Batch 4', 'Batch 5', 'Batch 6', 'Batch 7', 'Batch 8')
colnames(dfs) <- cols

f3a <- upset(dfs, sets=cols, sets.bar.color = "#56B4E9",
          mainbar.y.label = "Proteins identified in common", sets.x.label = "Batch size", 
          text.scale = c(2, 1.5, 2, 1.5, 1.5, 1.25),
          order.by="freq", keep.order=TRUE)

pdf(file = 'upset_plot_allbatches.pdf', width=12, height=7)
f3a
dev.off()

 # Fig 3B
data <- read.csv('../data/input_qc_fig3b.csv')

s1 <- ggbarplot(data, x='batch', y='avg_intensities',
                 xlab='', ylab='Mean protein intensities', width=0.4, 
                 title="") +
  theme(axis.title.y = element_text( size = 12, face = "bold" ))

s2 <- ggbarplot(data, x='batch', y='avg_num_peptides',
                xlab='', ylab='Mean peptides per protein', width=0.4, 
                title="") +
  theme(axis.title.y = element_text( size = 12, face = "bold" ))

s12 <- ggarrange(s1, s2, nrow=2)
ggsave('Figure3B.pdf', dpi=300)

# --- Figure 3C --------
data <- read.csv('../data/unscaled_non_malignants_samples.csv')

x <- cor(data)
colnames(x) = c("MCF10A rep1", "MCF10A rep2", 'MCF10A rep3', 
                'MCF10A rep4', 'MCF10A rep5', 
                'hME1')
rownames(x) = c("MCF10A rep1", "MCF10A rep2", 'MCF10A rep3', 
                'MCF10A rep4', 'MCF10A rep5', 
                'hME1')

pdf(file = 'Figure3B.pdf', width=8, height=8)
f3b <- corrplot.mixed(x, lower.col = "black", number.cex = .7, tl.cex=0.8)
dev.off()


