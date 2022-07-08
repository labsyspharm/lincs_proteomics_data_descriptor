library(data.table)
library(ggplot2)
library(UpSetR)
library(tidyr)
library(stringr)
library(dplyr)


#Load batch 1 and batch2 data and create individual dataframes for each set
df1 <- fread('../data/Sets1-4_protein_quant_with_peptides_24348.tsv', data.table=FALSE)
df1a <- df1[df1$Class == 'Set1', ]
df1b <- df1[df1$Class == 'Set2', ]
df1c <- df1[df1$Class == 'Set3', ]
df1d <- df1[df1$Class == 'Set4', ]

df2 <- fread('../data/Sets5-8_protein_quant_with_peptides_24344.tsv', data.table=FALSE)
df2a <- df2[df2$Class == 'Set5', ]
df2b <- df2[df2$Class == 'Set6', ]
df2c <- df2[df2$Class == 'Set7', ]
df2d <- df2[df2$Class == 'Set8', ]

# GEnerate a column of unique identifiers by merging uniprot Id and peptide sequence
generate_identifier_column <- function(data){
  data$identifier <- mapply(function(x, y){sprintf("%s__%s", x, y)}, data$ProteinId, data$PeptideSequence)
  return(data)
}

df1a <- generate_identifier_column(df1a)
df1b <- generate_identifier_column(df1b)
df1c <- generate_identifier_column(df1c)
df1d <- generate_identifier_column(df1d)

df2a <- generate_identifier_column(df2a)
df2b <- generate_identifier_column(df2b)
df2c <- generate_identifier_column(df2c)
df2d <- generate_identifier_column(df2d)

# Load metadata file
dfm <- fread('../data/metadata.csv', data.table=FALSE)
dfm$tmt_tag <- sapply(dfm$tmt_label, function(x){unlist(str_split(x, '~'))[2]})
dfm$tmt_tag <- sapply(dfm$tmt_tag, function(x){unlist(str_split(x, '_sum'))[1]})
dfm$tmt_tag <- sapply(dfm$tmt_tag, function(x){str_replace(x, '131C', '131c')})


# Function to remove reverse hits and contaminants
remove_revs_cont <- function(data){
  data <- data[!data$ProteinId %like% '_contaminant', ]
  data <- data[!data$ProteinId %like% '##', ]
  return(data)
}

# Function to rename columns by tmt tag to cell line name
rename_columns <- function(df, meta, run){
  meta <- meta[meta$run == run, ]
  tags <- meta$tmt_tag
  cs <- meta$cell_line
  
  colnames(df) <- dplyr::recode(
    colnames(df), 
    !!!setNames(as.character(meta$cell_line), meta$tmt_tag)
  )
  return(df)}


# Fuctnin that aggregates duplicate protein-peptide  rows
magg <- function(data, meta, run){
  groupvars <- c("ProteinId", "GeneSymbol", "Description", "PeptideSequence", "Class")
  sub_meta <- meta[meta$run == run,]
  tags <- sub_meta$tmt_tag
  dbg <-data %>% 
    group_by(ProteinId, PeptideSequence, GeneSymbol, Description, Class) %>% 
    dplyr::summarise(across(tags, mean))
  return(data.frame(dbg))
}

# For any given set, normalize each peptide in each sample by the multiplying with ratio of total 
# sample intensity and total bridge intensity
within_set_norm <- function(data, data_bridge, meta, run){
  data <- magg(data, meta, run)
  #ref_data <- magg(ref_data, meta, 'run4')
  data$identifier <- mapply(function(x, y){sprintf("%s__%s", x, y)}, 
                            data$ProteinId, data$PeptideSequence)
  
  
  sub_meta <- meta[meta$run == run,]
  tags <- sub_meta$tmt_tag
  nrm_fac <- colSums(data[, tags]) / sum(data[, data_bridge])
  data[tags] <- data[tags]*nrm_fac
  return(data)
}

dw1 <- within_set_norm(df1a, 'rq_131_sn', dfm, 'run1')
dw2 <- within_set_norm(df1b, 'rq_131_sn', dfm, 'run2')
dw3 <- within_set_norm(df1c, 'rq_131_sn', dfm, 'run3')
dw4 <- within_set_norm(df1d, 'rq_131_sn', dfm, 'run4')
dw5 <- within_set_norm(df2a, 'rq_131_sn', dfm, 'run5')
dw6 <- within_set_norm(df2b, 'rq_131_sn', dfm, 'run6')
dw7 <- within_set_norm(df2c, 'rq_131_sn', dfm, 'run7')
dw8 <- within_set_norm(df2d, 'rq_131_sn', dfm, 'run8')


# For any given set that has been normalized within batch, normalize peptide intensity by 
# multiplying with ratio of total intensity of bridge for that set and total intensity of
# bridge for reference set
across_set_norm <- function(data, ref_data, data_bridge, ref_bridge, meta, run){
  #data <- agg(data, dfm, run)
  #ref_data <- agg(ref_data, dfm , 'run4')
  norm_factor <- sum(data[[data_bridge]], na.rm=TRUE)/sum(ref_data[[ref_bridge]], na.rm=TRUE)
  sub_meta <- meta[meta$run == run,]
  tags <- sub_meta$tmt_tag
  data[tags] <- data[tags]/norm_factor
  #data$identifier <- mapply(function(x, y){sprintf("%s__%s", x, y)}, 
  #                          data$ProteinId, data$PeptideSequence)
  data <- rename_columns(data, meta, run)
  return(data)
}


db1 <- across_set_norm(dw1, dw4, 'rq_131_sn', 'rq_131_sn', dfm, 'run1')
db2 <- across_set_norm(dw2, dw4, 'rq_131_sn', 'rq_131_sn', dfm, 'run2')
db3 <- across_set_norm(dw3, dw4, 'rq_131_sn', 'rq_131_sn', dfm, 'run3')
db4 <- across_set_norm(dw4, dw4, 'rq_131_sn', 'rq_131_sn', dfm, 'run4')
db5 <- across_set_norm(dw5, dw4, 'rq_131c_sn', 'rq_131_sn', dfm, 'run5')
db6 <- across_set_norm(dw6, dw4, 'rq_131c_sn', 'rq_131_sn', dfm, 'run6')
db7 <- across_set_norm(dw7, dw4, 'rq_131c_sn', 'rq_131_sn', dfm, 'run7')
db8 <- across_set_norm(dw8, dw4, 'rq_131c_sn', 'rq_131_sn', dfm, 'run8')

# Scale peptide intenisty such that total intensity for any given peptide across
# all samples in a set is 100
scale <- function(data, dfm){
  csm <- dfm$cell_line
  cs <- colnames(data)
  cs <- cs[cs %in% csm]
  rs <- rowSums(data[cs])
  data[cs] <- 100 * data[cs]/rs
  return(data)
}

ds1 <- scale(db1, dfm)
ds2 <- scale(db2, dfm)
ds3 <- scale(db3, dfm)
ds4 <- scale(db4, dfm)
ds5 <- scale(db5, dfm)
ds6 <- scale(db6, dfm)
ds7 <- scale(db7, dfm)
ds8 <- scale(db8, dfm)



dfl <- list(ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8)
merge_dfs <- function(df_list, groupvars){
  dfa <- Reduce(
    function(x, y, ...) merge(x, y, by=groupvars, all = TRUE, ...), 
    df_list
  )
  return(dfa)
}

groupvars < c("ProteinId", "GeneSymbol" , "Description", "PeptideSequence", "Class")
dfr  <- merge_dfs(dfl, groupvars[1:4])
dfr <- remove_revs_cont(dfr)
pca_cols <- dfm$cell_line#[dfms$sample %in% colnames(dfr)]
pca_cols  <- append(c("ProteinId", "GeneSymbol", "Description", "PeptideSequence"), pca_cols)
dfr <- dfr[, pca_cols]


write.csv(dfr, "../data/merged_peptide_data_normalized.csv", row.names=F)

dfr2 <- dfr %>% drop_na()

run_pca <- function(data, meta , color_by){
  #data <- remove_revs_cont(data)
  
  cs <- meta$cell_line
  xi <- na.omit(data[, cs])
  pca <- prcomp(xi, center = TRUE)#,scale. = TRUE, na.action=na.omit)
  pca_var <- data.frame(summary(pca)$importance)
  pc1_var <- pca_var['Proportion of Variance', 'PC1'] * 100
  pc2_var <- pca_var['Proportion of Variance', 'PC2'] * 100
  xlab <- sprintf('PC1 (%s%%)', pc1_var)
  ylab <- sprintf('PC2 (%s%%)', pc2_var)
  df_pca <- data.frame(pca$rotation)
  df_pca$cell_line <- rownames(df_pca)
  df_pca <- merge(df_pca, meta, by.x='cell_line', by.y='cell_line')
  pl <- ggplot(df_pca) + geom_point(aes_string(x = 'PC1', y='PC2', color=color_by)) + theme_bw() +
    scale_color_manual(values=ms_colors) + 
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    geom_vline(xintercept=0, linetype="dashed", color = "grey") + 
    xlab(xlab) + ylab(ylab) + 
    theme( axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks.x = element_blank(),
           axis.ticks.y = element_blank(),
           axis.title.x = element_text( size = 12, face = "bold" ),
           axis.title.y = element_text( size = 12, face = "bold" ),
           plot.title = element_text( size=15, hjust = 0.5, face = "bold" ),
           panel.spacing = unit(c(1,1,1,1),"cm"),
           legend.position = 'none',
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank()
           ) 
  pl  <- pl + geom_label_repel(data=df_pca,
                               aes_string(x = 'PC1', y='PC2', label='cell_line'), size=2)
  return(pl)
}

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


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
run_pca(dfr2, dfms, 'molecular_subtype')
ggsave('PCA_peptide_level.pdf', 
       height = 6, width=8,
       dpi=600)


# PLot peptide level correlations
mcf10a_samples <- pca_cols[pca_cols %like% 'MCF10A']
nm_samples <- append(mcf10a_samples, 'hME1')
xp <- dfr2[, nm_samples]

xc <- cor(xp)
colnames(xc) = c("MCF10A rep1", "MCF10A rep2", 'MCF10A rep3', 
                 'MCF10A rep4', 'MCF10A rep5', 
                 'hME1')
rownames(xc) = c("MCF10A rep1", "MCF10A rep2", 'MCF10A rep3', 
                 'MCF10A rep4', 'MCF10A rep5', 
                 'hME1')
pdf(file = 'corrplot_peptide_level.pdf', 
    width=8, height=8)
corrplot::corrplot.mixed(xc, lower.col = "black", number.cex = .7, tl.cex=0.8)
dev.off()




