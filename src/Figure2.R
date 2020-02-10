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


dm <- read.csv('../data/metadata.csv')
dm <- dm[!dm$molecular_subtype %in% c('Bridge', 'ovarian'), ]

dm1 <- data.frame(table(dm$molecular_subtype))
colnames(dm1) <- c('molecular_subtype', 'num_cell_lines')
dm1 <- dm1[dm1$num_cell_lines > 0, ]
ms_colors <- add.alpha(c("#ca7832", "#6d73b0",
                         '#4dab9b',
                         #'#868686B2',
                         'lightblue',
                         #'red',
                         "#7da245", "#c85f61"  ), 0.7)
names(ms_colors) <- c('Basal A', 'Basal B',
                      'Basal',
                      #'Bridge',
                      'Basal A/Luminal',
                      #'ovarian',
                      'Luminal', 'Non malignant, Basal' )
dm1$color <- ms_colors[as.character(dm1$molecular_subtype)]

f2a <- ggbarplot(dm1, x='molecular_subtype', y='num_cell_lines', palette=c(dm1$color),
                 fill='molecular_subtype', order=as.character(dm1$molecular_subtype),
                 xlab='', ylab='number of cell lines', width=0.5, legend="none",
                 title="Transcriptional subtypes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_text( size = 12, face = "bold" ))


#---Figure 2b -----

dm2 <- data.frame(table(dm$receptor_status))
colnames(dm2) <- c('receptor_status', 'num_cell_lines')
dm2 <- dm2[dm2$num_cell_lines > 0, ]
rs_colors <- add.alpha(c('#c85f61', '#4dab9b', '#ca7832', '#7da245'
                         #'#868686B2', 'red'
                         ), 0.7)
names(rs_colors) <- c('NM', 'Her2amp', 'TNBC', 'HR+')#, 'Bridge', 'ovarian')


dm2$color <- rs_colors[as.character(dm2$receptor_status)]

f2b <- ggbarplot(dm2, x='receptor_status', y='num_cell_lines', palette=c(dm2$color),
                 fill='receptor_status', order=as.character(dm2$receptor_status),
                 xlab='', ylab='number of cell lines', width=0.4, legend="none",
                 title="Receptor status") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          axis.title.y = element_text( size = 12, face = "bold" ))


#---Figure 2c------------------------------------------------------------------

dfr <- read.csv(file='../data/mass_spec_pca_all_batches.csv',        
                stringsAsFactors = TRUE,header=TRUE)
dfr <- dfr[!dfr$molecular_subtype %in% c('Bridge', 'ovarian'), ]

f2c <- ggplot(dfr, aes(x=pc_1, y=pc_2, color=molecular_subtype)) + 
  geom_point(size=3, alpha=0.7) + theme_bw() + 
  scale_color_manual(values=ms_colors) +
  #scale_color_jco(name='transcriptional subtype') + 
  xlab('PC 1 (14.1%)') + ylab('PC 2 (7.4%)') +
  ggtitle("Proteome") +
  geom_label_repel(aes(label = cell_line),
                   size=2.5,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.size=0.25) +
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

pdf(file = 'Figure2.pdf', width=12, height=7)
f2abc
dev.off()

