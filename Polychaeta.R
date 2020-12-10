#' ---
#' title: "Polychaeta molecular evolution"
#' author: Jason Moggridge
#' output:
#'   pdf_document
#' ---
#'
#'
#+ echo=FALSE, message=FALSE, warning=FALSE

# Polychaete COI-5p barcode notebook

# 
# rmarkdown::render(
#   "./R code/Polychaeta.R",
#   output_dir = './',
#   output_file = "Polychaeta report"
# )

rm(list=ls())
# Libraries
library(tidyverse)
library(rentrez)
library(bold)
library(Biostrings)
library(DECIPHER)
library(dendextend)
library(ape)
library(beepr)
library(rcartocolor)
library(ggthemes)
library(patchwork)




filter_CO1_seqs <- function(bold.df, markers, minlen, maxlen, n_threshold){
  ## Filter data and tidy up COI sequences
  
  bold.df <- bold.df %>%
    # replace blanks with na's
    mutate(across(where(is.character), ~na_if(.x, ''))) %>%
    # make strings factors, except for sequences
    mutate(across(where(is.character), as.factor)) %>%
    mutate(nucleotides = as.character(nucleotides))  %>%
    # Keep only specimens with desired marker gene sequences
    filter(!is.na(nucleotides) & markercode %in% markers) %>%
    # Remove gaps and outer N's
    mutate(nucleotides.trim = str_remove_all(nucleotides, '\\s|-|^N+|N+$')) %>%
    # Filter to parameters shortest, longest, with proportion n < n_threshold 
    mutate(slen = nchar(nucleotides.trim)) %>%
    filter(slen > minlen & slen < maxlen & 
             !str_count(nucleotides.trim, 'N')/slen > n_threshold) %>%
    # keep only specimens with taxonomy info to species-rank
    filter(!is.na(genus_name) & !is.na(species_name)) %>%
    filter(!str_detect(species_name, ' sp.| cf.'))

  return(bold.df)  
}

distVector <- function(model, dnabin){
  # return numeric vector of the upper triangle of the distance matrix
  dist <- dist.dna(x = dnabin,
                   model = model,
                   as.matrix = FALSE,
                   pairwise.deletion = TRUE)
  return(as.numeric(dist))
}


## Download data ----
#' 
#' ```
#' # BOLD database results for taxon 'Polychaeta'; downloaded on 2020/10/31.
#' bold::bold_stats(taxon='Polychaeta')
#' polychaeta.bold <- bold::bold_seqspec(taxon = 'Polychaeta')
#' write_rds(polychaeta.bold, './data/polychaeta.bold.rds', compress = 'gz')
#' ```
#'
#+ echo=FALSE, message=FALSE, warning=FALSE

# Prepare BOLD data for analysis ----
polychaeta.raw <- read_rds('./data/polychaeta.bold.rds')
dim(polychaeta.raw)

polychaeta.CO1 <- filter_CO1_seqs(polychaeta.raw, 'COI-5P', 600, 700, 0.01)
dim(polychaeta.CO1)
rm(polychaeta.raw)

# some sequences that did not align well, could be mislabelled
poorly_aligned <- c(
  'GBAN1575-07', 'GBAN1576-07', 'GBAN2779-10', 'GBAN6823-15', 'GBAN13543-19', 'GBAN13546-19', 'GBAN13547-19', 'GBAN13548-19', 'GBAN13550-19', 'GBAN13551-19', 'GBAN13552-19', 'GBAN13553-19', 'GBMNB34178-20', 'GBMNB34179-20', 'GBMNB34180-20', 'GBMNB34182-20', 'GBMNB34183-20', 'CCFZP144-19', 'HZPLY376-13'
)
# check and see if all from a single family or what not
polychaeata.notaligned <- polychaeta.CO1 %>%
  filter(processid %in% poorly_aligned) %>%
  mutate(family_name = fct_drop(family_name))
summary(polychaeata.notaligned$family_name)


# Prepare down-sampled dataframe for alignment
set.seed(1)
polychaeta.CO1.sample <- polychaeta.CO1 %>%
  filter(family_name == 'Polynoidae') %>%
  # take 1 sequence per genus
  group_by(species_name) %>%
  sample_n(1) %>%
  # exclude sequences known to align poorly
  filter(!processid %in% poorly_aligned) %>%
  data.frame() %>%
  # create labels for the alignment
  mutate(name = paste(word(species_name, 1)), row_number()) %>%
  # tidy up factors and remove unneeded cols
  mutate(across(where(is.factor), fct_drop)) %>%
  select(processid, name, nucleotides.trim, slen, 
         family_name, genus_name, species_name)
#
glimpse(polychaeta.CO1.sample)
# summary(polychaeta.CO1.sample)
rm(polychaeata.notaligned)

# Make 'seq' a DNA stringset
polychaeta.CO1.sample$seq <-
  Biostrings::DNAStringSet(polychaeta.CO1.sample$nucleotides.trim)
# Attach labels with taxonomy + id
names(polychaeta.CO1.sample$seq) <- polychaeta.CO1.sample$name


a <- ggplot(polychaeta.CO1.sample, aes(x = slen)) +
  geom_density() + 
  xlab('Sequence length') +
  geom_rangeframe() + 
  theme_tufte()
b <- polychaeta.CO1.sample %>%
  mutate(`% GC` = str_count(seq, '[GC]')/slen*100) %>%
  ggplot(aes(x=`% GC`)) + 
  geom_density() +
  geom_rangeframe() + 
  theme_tufte()
a/b
rm(a,b)

# keep only ids, names and seqs for simplicity
polychaeta.seqs <- polychaeta.CO1.sample %>%
  dplyr::select(!contains('nucleotides')) 

rm(polychaeta.CO1, polychaeta.CO1.sample, poorly_aligned)



## Alignment -----

# This is the seq which others will be oriented to:
# polychaeta.seqs[which(polychaeta.seqs$slen == max(polychaeta.seqs$slen)),]

# Reorient any sequences that might align better as rev-complement...
polychaeta.seqs$seq <- OrientNucleotides(polychaeta.seqs$seq)

# Perform alignment with decipher defaults
polychaeta.align <- DECIPHER::AlignSeqs(polychaeta.seqs$seq, verbose = TRUE)


# translated alignment performs poorly
# polychaeta.align2 <- DECIPHER::AlignTranslation(polychaeta.seqs$seq, verbose = TRUE)
# BrowseSeqs(polychaeta.align2)

# transform decipher alignment into DNAbin object
polychaeta.bin <- as.DNAbin(polychaeta.align)
rm(polychaeta.align)


## Clustering ----
# names(polychaeta.bin) <- word(names(polychaeta.bin), 1)

# NJ clustering by K2P distance metric
cluster.COI.seq <- function(model, threshold, algo){
  DECIPHER::IdClusters(
    dist.dna(x = polychaeta.bin, 
             model = model,
             as.matrix = FALSE, 
             pairwise.deletion = TRUE),
    method = algo,
    cutoff = threshold,
    showPlot = FALSE,
    type = "both")
} 

JC.clusters <- cluster.COI.seq('raw', 0.2, 'UPGMA')
JC.dend <- JC.clusters[[2]] %>%
  set('labels_cex', c(0.5))
  

ggplot(JC.dend, horiz=T) + theme(plot.margin = margin(1,2,1,1, 'cm'))


clusters.k80 <- cluster.COI.seq('k80', 0.2, 'UPGMA')
# plot(clusters.k80[[2]])
clusters.tn93 <-  cluster.COI.seq('tn93', 0.2, 'UPGMA')

## Visually compare phylogenetic inferences among models (JC, K2P, TN93)

# JC vs k80

clusters.k80[[2]] <- clusters.k80[[2]] %>%
  dendextend::match_order_by_labels(clusters.JC[[2]])
COI.dendrograms <- dendlist(clusters.JC[[2]], clusters.k80[[2]])
dendextend::tanglegram(COI.dendrograms, common_subtrees_color_branches = FALSE) 

# JC vs tn93
clusters.tn93[[2]] <- clusters.tn93[[2]] %>%
  dendextend::match_order_by_labels(clusters.JC[[2]])
COI.dendrograms2 <- dendlist(clusters.JC[[2]], clusters.tn93[[2]])
dendextend::tanglegram(COI.dendrograms2, common_subtrees_color_branches = FALSE)

# k80 vs tn93
clusters.tn93[[2]] <- clusters.tn93[[2]] %>%
  dendextend::match_order_by_labels(clusters.k80[[2]])
COI.dendrograms3 <- dendlist(clusters.k80[[2]], clusters.tn93[[2]])
dendextend::tanglegram(COI.dendrograms3, common_subtrees_color_branches = FALSE)




## Maximum-likelihood clustering
clusters.ML <- DECIPHER::IdClusters(
  polychaeta.align,
  myDistMatrix = dist.dna(
    x = polychaeta.bin, model = 'k80',
    as.matrix = FALSE, pairwise.deletion = TRUE),
  method = 'ML',
  processors = 4,
  cutoff = 0.2,
  showPlot = TRUE,
  type = "both",
  verbose = TRUE)
write_rds(clusters.ML, './data/clusters.ML.rds')


# Molecular evolution analysis ----

# create a df with all the distance measures as columns
dist.df <- tibble(model = c('raw','N', 'k80', 'tn93', 'gg95', 'TS', 'TV')) %>%
  mutate(data = lapply(model, function(x) distVector(x, polychaeta.bin))) %>%
  pivot_wider(names_from = 'model', values_from = 'data') %>%
  unnest(cols = raw:TV)

glimpse(dist.df)


# plot p-distance ~ k80 distance
raw_k80.plot <- dist.df %>%
  ggplot(aes(x=k80, y=raw)) + 
  geom_point(alpha = 0.5, size = 0.5, pch = 1) +
  # add diagonal line
  geom_path(
    data = data.frame(
      x = seq(0, 1, 0.05),
      y = seq(0, 1, 0.05),
      group = factor('group')
    ),
    aes(x = x, y = y, group = group),
    alpha = 0.7,
    lty = 2
  ) +
  theme_classic() +
  xlim(c(0, 0.72)) +
  ylim(c(0, 0.45)) +
  labs(x = 'K2P distance', y = 'p-distance',
       title = 'Saturation plot: distances of Polychaete CO1 sequences',
       subtitle = 'p-distance underestimates greater distances',
       caption = paste0('data from BOLD. (', nrow(polychaeta.bin),
                        '), 2020/10/31'))
raw_k80.plot

# plot transitions and transversions ~ k80 distance
ts_tv.plot <- dist.df %>%
  select(k80, TS, TV) %>%
  # manipulate table into long-form
  pivot_longer(cols=c(TS, TV),
               names_to = 'Substitution', values_to = 'Proportion') %>%
  # create scatterplot w smooth fits
  ggplot(aes(x = k80, y = Proportion, color = Substitution, group = Substitution)) +
  geom_point(size=0.04, alpha = 0.3) +
  geom_smooth(method = 'loess', formula = 'y~x', se=FALSE) +
  theme_classic() +
  scale_color_carto_d('Type') +
  guides(
    colour = guide_legend(override.aes = list(size=5, alpha=1, pch=15))
    ) +
  labs(
    x = 'K2P distance',
    y = 'Substiutions',
    title = '',
    subtitle = 'Transitions and transverions by K2P distance in Polychaeta COI sequences',
    caption = paste0('Data from BOLD; accessed on 2020/10/31')
    )
# takes a long time to plot
ts_tv.plot
beep(2)

## do same plot with other distances (JC, )

rm(raw_k80.plot, ts_tv.plot, dist.df)




# just trying to figure out what I can do with this data? in comparing trees?
x <- data.frame(
  clusters = as.integer(clusters.COI[[1]]$cluster),
  names = rownames(clusters.COI[[1]])
  ) %>%
  mutate(family = word(names, 1),
         genus = word(names, 2)) %>%
  arrange(clusters)

glimpse(x)
nrow(x)
length(unique(x$cluster))
rownames(clusters.COI[[1]]) <- paste(x$family, x$genus)
rownames(clusters.COI[[2]]) <- paste(x$family, x$genus)
glimpse(x)



# NJ clustering by K2P distance metric
clusters.COI <- DECIPHER::IdClusters(
  dist.dna(x = polychaeta.bin, 
           model = 'k80',
           as.matrix = FALSE, 
           pairwise.deletion = TRUE),
  method = 'NJ',
  processors = 4,
  cutoff = 0.3,
  showPlot = FALSE,
  type = "both",
  verbose = TRUE)
plot(clusters.COI[[2]])
# clusters.COI[[2]]

## Maximum-likelihood clustering
clusters.ML <- DECIPHER::IdClusters(
  polychaeta.align,
  myDistMatrix = dist.dna(
    x = polychaeta.bin, model = 'k80',
    as.matrix = FALSE, pairwise.deletion = TRUE),
  method = 'ML',
  processors = 4,
  cutoff = 0.2,
  showPlot = TRUE,
  type = "both",
  verbose = TRUE)
write_rds(clusters.ML, './data/clusters.ML.rds')

# onstructing initial neighbor-joining tree:
# JC69:     -ln(L)=147278, AICc=290910, BIC=300410
# JC69+G4:  -ln(L)=154643, AICc=305648, BIC=315147
# K80:      -ln(L)=107559, AICc=211480, BIC=220979
# K80+G4:   -ln(L)=111017, AICc=218404, BIC=227901
# F81:      -ln(L)=163321, AICc=323019, BIC=332514
# F81+G4:   -ln(L)=167780, AICc=331945, BIC=341439
# HKY85:    -ln(L)=163201, AICc=322787, BIC=332281
# HKY85+G4: -ln(L)=167783, AICc=331960, BIC=341453
# T92:      -ln(L)=163730, AICc=323829, BIC=333326
# T92+G4:   -ln(L)=169946, AICc=336269, BIC=345765
# TN93:     -ln(L)=163201, AICc=322796, BIC=332289
# TN93+G4:  -ln(L)=167735, AICc=331872, BIC=341363
# 
# The selected model was:  K80
# 
# Maximizing Likelihood of Tree:
#   -ln(Likelihood) = 96246 (10.52% improvement), 315 NNIs 
# 
# Transition rates = 0.600
# Transversion rates = 1
# Time difference of 687.45 secs

## Plot dendrogram
clusters.ML[[2]]
plot(clusters.ML[[2]], main = 'ML tree')

length(unique(clusters.ML[[1]]$cluster))
glimpse(clusters.ML[[1]])
class(clusters.ML[[2]])
# glimpse(clusters.ML)




