
###############################################################################
# INTRO
# LeFSE analysis; top differentially abundant taxa in fecal samples of FMT recipients day 0 of C.difficile challenge (Fig 3D)
# 12-26-21
#Madeline R Barron
###############################################################################

# Goal: Generate a plot of the top differentially abundant bacterial taxa between FMT recipients that did or do not become colonized
#by C. difficile (fecal samples day of C.difficile challenge)

###############################################################################
# (01) LIBRARIES
###############################################################################
library(ggtext)
library(RColorBrewer)
library(tidyverse)


###############################################################################
# (02)READ DATA
###############################################################################

#   (A)Read in metadata file
meta <-read.csv("data/antip40_exp1_exp2_exp3_fmt_metadata.csv")

#   (B) Read in lefse summary file output by mothur; depicts OTUs that differ significantly between samples susceptiblevs resistant to C.difficile
lefse <- read.csv("data/fmt_d0_lefse.0.03.lefse_summary.csv") %>%
  rename_all(tolower) #Read in lefse summary file (mothur output)

#   (C)Read in shared file indicating counts of OTUs in samples on D0 C.diff challenge (7 days post-FMT) used in lefse

lefse_shared <- read.csv("data/fmt_d0_lefse.shared.csv")%>%
  select(-label, -numOtus) %>%
  pivot_longer(cols=-group, names_to="otu", values_to="count") %>%
  rename_all(tolower)

#   (D) Read in taxonomy file

tax <- read.csv("data/final.taxonomy.csv") %>% #Reading in taxonomy info
  rename_all(tolower) %>% # Formatting colnames
  mutate(taxonomy = str_replace_all(taxonomy, "(\\(\\d+\\)|\\;$)", "")) %>% # Removing confidence scores and trailing ';'
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";") %>% # Breaking tax into levels
  gather(key = level, value = classification, -otu, -size)

# (E) Read in file denoting which partition (community type) samples fall into

design <- read.csv("data/fmt_d0_lefse_design.csv")

design_meta <- left_join(meta, design) %>%
              filter(day == 0,
                     mouse_genotype.x != "il10",
                     treatment != "No_treatment",
                     fmt_donor_recipient != "donor")


###############################################################################
# (03) MERGE LEFSE AND META DATA
###############################################################################

#     (A) Merge lefse summary file with shared, filtering so only OTUs where LDA >= 3 are shown
lefse_shared_2 <- inner_join(lefse_shared, lefse) %>%
  filter(lda >= 3)

#     (B) Merge meta and lefse summary/shared data
lefse_meta <- inner_join(lefse_shared_2,meta)

###############################################################################
# (03) CALCULATE RELATIVE ABUNDANCE
###############################################################################

#     (A)  Relative abundance
rel_abund <- lefse_meta%>%
  group_by(group) %>% # Group by sample col
  mutate(total_count = sum(count)/6, # Calculating total reads per sample and appending to df
         rel_abund = count/total_count) %>%  # Calculating rel abund for each Otu and appending to df
  ungroup() #The result is a table containing relative abundance of each OTU for each sample

smallest_non_zero <- rel_abund %>%
    filter(rel_abund > 0) %>%
    slice_min(rel_abund) %>%
    pull(rel_abund) %>% .[1]
rel_abund <- rel_abund %>%
    mutate(rel_abund_c = rel_abund + smallest_non_zero / 10)
###############################################################################
# (04) MERGE TAXONOMY AND SAMPLE INFORMATION
###############################################################################

#     (A) Join taxonomy with relative abundance dataframe
rel_abund_tax <- inner_join(rel_abund, tax)

# (B) Join sample community type with relative abundance/taxonomy info

rel_abund_design <- left_join(rel_abund_tax, design) #%>%
#select(group, rel_abund, class, partition, level, classification,
#otu, treatment_no_cd, pvalue, lda) #grouping by these parameters

# (B)
agg_otu_rel_abund <- rel_abund_design %>%
  filter(level == "genus") %>%
  group_by(group, rel_abund_c, ever_pos_cdiff,
           otu, treatment, pvalue, classification,lda) %>% # Aggregates classifications based on treatment
  summarize(agg_rel_abund_c = sum(rel_abund_c)) %>% # Calculating aggregated rel abund
  ungroup()


#    (C) Finding tax classifications for a given tax level
top_otu_tax <- agg_otu_rel_abund %>%
  group_by(classification, group, ever_pos_cdiff) %>%
  summarize(mean_rel_abund_c = mean(agg_rel_abund_c)) %>%
  mutate(tax_label = case_when(
      str_detect(classification, '_') ~ str_replace(classification, "(^\\w+)_(.*)", "_\\1_ \\2"),
                               TRUE ~ paste0('_', classification, '_'))
         )

###############################################################################
# (05) PLOTTING
###############################################################################

#   (A) Plotting log10 relative abundance in each sample identified as significant (LDA >= 2)
# in lefse analysis. We are grouping by the C.diff susceptibility of each sample.

colors <- c("no" = "#404040",
            "yes" = "#999999")

abund <- top_otu_tax %>%
  ggplot(aes(x = tax_label, y = mean_rel_abund_c, color = ever_pos_cdiff)) +
  geom_boxplot(lwd = 1) +
  geom_hline(yintercept = smallest_non_zero, linetype = 'dashed') +
  labs(y = expression("Relative Abundance ("*Log[10]*")"),
       color = "Ever Positive _C.diff_") +
  scale_color_manual(values = colors)+
  scale_y_log10() +
  #stat_summary(aes(group=ever_pos_cdiff),fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",
               #width=0.5,colour="black",position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_markdown(size = 11, colour = "black", angle = 60, hjust=1),
        axis.title.x = element_blank(),
        legend.title = element_markdown(),
        text = element_text(colour = "black", size = 11))

abund
ggsave('figures/raw/Fig4_rel-abun.pdf', plot = abund)
#Note: Figure labels and legend modified in Adobe Illustrator.
