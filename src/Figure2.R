library(ggpubr)
library(ggrepel)
                                     
# ---Figure 2a ----
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

ms_colors <- add.alpha(c("#ca7832", "#6d73b0", #"#868686B2",
               "#7da245", "#c85f61",  "#868686B2"), 0.7)
names(ms_colors) <- c('Basal A', 'Basal B',# 'Bridge',
                      'Luminal', 'Non malignant, Basal', 'Basal')
dm <- read.csv('../data/mol_subtype_breakdown.csv')
dimnames(dm) <- list(as.character(dm$molecular_subtype),
                     colnames(dm))
dm$color <- ms_colors[as.character(dm$molecular_subtype)]


f2a <- ggbarplot(dm, x='molecular_subtype', y='num_cells', palette=c(dm$color),
                 fill='molecular_subtype', order=as.character(dm$molecular_subtype),
                 xlab='', ylab='number of cell lines', width=0.5, legend="none",
                 title="Transcriptional subtypes") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(angle=45, hjust=1))

#---Figure 2b -----

rs_colors <- add.alpha(c('#c85f61', '#4dab9b', '#ca7832', '#7da245'), 0.7)
names(rs_colors) <- c('NM', 'HER2amp', 'TNBC', 'HR+')

dr <- read.csv('../data/receptor_status_breakdown.csv')
dimnames(dr) <- list(as.character(dr$receptor_status),
                     colnames(dr))
dr$color <- rs_colors[as.character(dr$receptor_status)]

f2b <- ggbarplot(dr, x='receptor_status', y='num_cells', palette=c(dr$color),
                 fill='receptor_status', order=as.character(dr$receptor_status),
                 xlab='', ylab='number of cell lines', width=0.4, legend="none",
                 title="Receptor status") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(angle=45, hjust=1))


#---Figure 2c------------------------------------------------------------------
ms_colors <- c("#ca7832", "#6d73b0", #"#868686B2",
               "#7da245", "#c85f61",  "#868686B2",  "#868686B2")  #"#7AA6DCB2
names(ms_colors) <- c('Basal A', 'Basal B',# 'Bridge',
                      'Luminal', 'Non malignant, Basal', 'Basal', 'unconfirmed')

dfr <- read.csv(file='../data/mass_spec_pca.csv',
                stringsAsFactors = TRUE,header=TRUE)
dfr <- dfr[dfr$molecular_subtype != 'Bridge', ]

f2c <- ggplot(dfr, aes(x=pc_1, y=pc_2, color=molecular_subtype)) + 
  geom_point(size=3, alpha=0.7) + theme_bw() + 
  scale_color_manual(values=ms_colors) +
  #scale_color_jco(name='transcriptional subtype') + 
  xlab('PC 1 (15.7%)') + ylab('PC 2 (8.6%)') +
  ggtitle("Proteome") +
  geom_label_repel(aes(label = sample),
                   size=5,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_hline(yintercept=0, linetype="dashed", color = "grey") +
  geom_vline(xintercept=0, linetype="dashed", color = "grey") + 
  theme( axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.x = element_text( size = 12, face = "bold" ),
         axis.title.y = element_text( size = 12, face = "bold" ),
         plot.title = element_text( size=15, hjust = 0.5, face = "bold" ),
         panel.spacing = unit(c(1,1,1,1),"cm"),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.position = 'none' ) 


f2ab <- ggarrange(f2a, f2b, nrow=2, ncol=1)
f2abc <- ggarrange(f2ab, f2c, ncol=2, nrow=1, widths=c(0.3, 0.7))

ggsave('Figure2.pdf')
