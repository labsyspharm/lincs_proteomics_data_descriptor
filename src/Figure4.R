library(ggplot2)

data <- read.csv('../data/Olaparib_AZD2281_imp.csv')
data$color = 'red'
data[data$pearman_rho  > 0, 'color'] = 'blue'
data$alpha <- mapply(function(x, y){add.alpha(x, y)}, data$color, abs(data$spearman_rho))
data <- data[1:10, ]

s1 <- ggbarplot(data, x='features', y='Olaparib_AZD2281',
                order=as.character(data$features),
                legend='none',
                #palette=c(data$alpha), fill='features',
                xlab='', ylab='feature importance', width=0.4, 
                title="Proteome based predictors of resistance to Olaparib") +
  theme(plot.title = element_text(hjust = 0.5, face='bold')) +
  theme(axis.title.y = element_text( size = 12, face = "bold" ))
  
  
# iBaq figure
data <- read.csv('../data/Bridge_UPS2_AbsoluteQuant_Results.csv')
dr <- read.csv('../data/UPS2_report.csv')


s1 <- ggplot(data, aes_string(x='log10_ibaq', y='Bridge_fmol_per_ug')) + 
  geom_point(color='grey', alpha=0.5) + 
  theme_bw() + xlab('log10 (iBAQ') + ylab('fmol / ug')

d2 <- data[data$uniprot_id %in% dr$UniprotID, ]
s1 <- s1 + geom_point(data=ds, aes_string(x='log10_ibaq', y='t'), color='red', alpha=.7 ) +
  geom_label_repel(data=ds, aes(label = uniprot_id),
                     size=2.5,
                     box.padding   = 0.35, 
                     poinds <t.padding = 0.5,
                     segment.color = 'grey50',
                     segment.size=0.25)

ds <- merge(d2, dr, by.x='uniprot_id', by.y='UniprotID')

#--------------------

pw <- data.table(x = c("EGFR", "GRB2", "SOS1", 'KRAS', 'BRAF', 'MAP2K1', 'MAPK1', 'DUSP4', 'SOS1', 'DUSP4', 'SOS1', 'SOS1', 'KRAS', 'CRAF', 'NRAS', 'NRAS' ), 
                 y = c('GRB2', 'SOS1', 'KRAS', 'BRAF', 'MAP2K1', 'MAPK1', 'DUSP4', 'MAPK1', 'MAPK1', 'MAPK1', 'NRAS', 'KRAS', 'CRAF', 'MAP2K1', 'CRAF', 'BRAF'))
net <- network(pw)

n <- ggnet2(net, color = "phono", palette = "Set1", size = 0, edge.size=0.8) +
  geom_point(aes(color = color), size = 1.7**s, color = "white") +
  geom_point(aes(color = color), size = 1.7**s, alpha = 0.5) +
  geom_point(aes(color = color), size = 1.7**s-5) +
  geom_text(aes(label = label), color = "black", fontface = "bold", size=3) +
  geom_text(aes(label = label), color = "white", fontface = "bold", size=3) +
  guides(color = FALSE) 

dt = read.csv('../../cancerbrowser_rshiny/tnbc_log10molecules_cell.csv')

s = c(5.2, 0, 3.9, 4.6, 6.42, 4.25, 5.4, 5.5, 6.03, 4.42 )

dt3$color <- rs_colors[as.character(dt3$receptor_status)]
dt3$ERBB2 <- as.numeric(as.character(dt3$ERBB2))

f5 <- ggbarplot(dt3, x='cell_line', y='ERBB2', palette=c(dt3$color),
                fill='cell_line', order=as.character(dt3$cell_line),
                xlab='', ylab='log10 (molecules per cell)', width=0.4, legend="none",
                title="ERBB2 protein abundance",
                ylim=c(4, 7.5)) +
  #ylim(c(2,7)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.y = element_text( size = 12, face = "bold" ))