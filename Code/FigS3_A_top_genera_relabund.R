###############################################################################
# Relative Abundance of Top 10 Genera DMM Analysis (Figure)
#01/17/22
#Madeline R. Barron
###############################################################################

# Goals: Plot relative abundance of top 10 bacterial genera in samples belonging to enterotype 1,
#2, and 3, as determined by DMM.  

############################################
# (01) LIBRARIES
###############################################################################
library(tidyverse)
library(RColorBrewer)
install.packages("rstatix")
library(rstatix)

############################################
# (02) READ DATA
###############################################################################

#Reading in taxonomy file and tidying up
tax <- read_csv(file = "data/final.taxonomy.csv") %>% 
  rename_all(tolower) %>% # Formatting colnames
  mutate(taxonomy = str_replace_all(taxonomy, "(\\(\\d+\\)|\\;$)", "")) %>% # Removing confidence scores and trailing ';'
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";") %>% # Breaking tax into levels
  gather(key = level, value = classification, -otu, -size)

#Reading shared file of samples and grouping by sample. Calculating relative abundances 
shared_2 <- read_csv("data/sample.final.shared.csv", col_types=cols(group=col_character())) %>% 
  select(-label, -numOtus) %>% #Get rid of label columns and numOTUs columns 
  pivot_longer(cols=-group, names_to="otu", values_to="count")

design <- read.csv("data/genus_cdiff_only_dmm.design.csv")

#Reading metadata
meta <- read_csv("data/dmm_metadata.csv") 

############################################
# (03) JOINING DATA
###############################################################################

#Joining taxonomy & OTU data frames
tax_otu_join <- inner_join(tax, shared_2)

#Joining taxonomy & OTU combined dataframe with metadata (uses joining by group (i.e. sample) by default)
agg_data <- inner_join(tax_otu_join, meta)

agg_data <- inner_join(agg_data, design)

############################################
# (04) CALCULATING RELATIVE ABUNDANCE 
###############################################################################

#Relative abundance 
rel_abund <- agg_data %>% 
  group_by(group) %>% # Groups by sample col
  mutate(total_count = sum(count)/6, # Calculating total reads per sample and appending to df
         rel_abund = count/total_count) %>%  # Calculating rel abund for each Otu and appending to df
  ungroup()

smallest_non_zero <- rel_abund %>%
  filter(rel_abund > 0) %>%
  slice_min(rel_abund) %>%
  pull(rel_abund) %>% .[1]

rel_abund <- rel_abund %>%
  mutate(rel_abund_c = rel_abund + smallest_non_zero / 10)


agg_otu_rel_abund <- rel_abund %>% 
  filter(level == "genus") %>% 
  group_by(group,                                     
           partition, 
           pos_cdiff_d1, 
           exp_id_number,
           treatment,
           classification) %>% # Aggregates classifications based on treatment
  summarize(agg_rel_abund_c = sum(rel_abund_c)) %>% # Calculating aggregated rel abund
  ungroup()

#Finding top tax classifications for a given tax level
top_otu_tax <- agg_otu_rel_abund %>% 
  group_by(classification) %>%
  summarize(mean_rel_abund_c = mean(agg_rel_abund_c)) %>% # Calculating mean aggregated relative abundance across all samples
  top_n(n = 10, wt = mean_rel_abund_c) %>% 
  filter(mean_rel_abund_c > 0) %>%  # Only choosing Otus with non-zero mean relative abundances
  pull(classification) #Creating vector of Otu tax classifications

# Reassigning any tax classification as "Other" if not in the top_otu_tax list
# If the df consists of samples (no controls)...
top_otu_data <- agg_otu_rel_abund %>% 
  mutate(classification = case_when(!(classification %in% top_otu_tax) ~ "Other",
                                    TRUE ~ classification)) %>% 
  group_by(group,                                     
           partition, 
           pos_cdiff_d1, 
           exp_id_number,
           treatment,
           classification) %>% # Need to collapse all of the newly labelled "Other" tax into single entries
  summarize(agg_rel_abund_c = sum(agg_rel_abund_c)) %>% 
  mutate(mean_rel_abund_c = log10(agg_rel_abund_c + 1)) %>%# Adds up all entries for "Other" classification
 filter(classification != "Other") %>% 
   ungroup()

###########################################
# (04) PLOTTING 
###############################################################################

abund<- top_otu_data %>% 
  ggplot(aes(x=mean_rel_abund_c,
                              y= classification,color = partition)) +
  geom_boxplot(lwd =0.75) +
  geom_vline(xintercept = smallest_non_zero, linetype = 'dashed', color = "gray") +
  labs(x = expression("Relative Abundance ("*Log[10]*")"),
       y = "OTU",
       color = "Community type") +
  scale_color_brewer(palette = "Set1")+
  scale_x_log10() + 
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 11, colour = "black", angle = 90, hjust=1),
        axis.title.x = element_blank(),
        legend.title = element_text(),
        text = element_text(colour = "black", size = 11))

abund + coord_flip()


###########################################
# (04) Statistics
###############################################################################

#Determining differences between enterotypes (partitions) for reach genus.

#Akkermansia

akker <- filter(top_otu_data, classification == "Akkermansia")


kruskal.test(mean_rel_abund_c ~ partition, data = akker)
dunn_test(mean_rel_abund_c~ partition, p.adjust.method = "bonferroni", data = akker)


#Alistipes


alis <- filter(top_otu_data, classification == "Alistipes")


kruskal.test(mean_rel_abund_c ~ partition, data = alis)
dunn_test(mean_rel_abund_c ~ partition, p.adjust.method = "bonferroni", data = alis)

#Bacteroides

bacteroides <- filter(top_otu_data, classification == "Bacteroides")
kruskal.test(mean_rel_abund_c ~ partition, data = bacteroides)
dunn_test(mean_rel_abund_c ~ partition, p.adjust.method = "bonferroni", data = bacteroides)

#Bifidobacterium

bifido <- filter(top_otu_data, classification == "Bifidobacterium")
kruskal.test(mean_rel_abund_c ~ partition, data = bifido)
dunn_test(mean_rel_abund_c ~ partition, p.adjust.method = "bonferroni", data = bifido)

#Erysipelotrichaceae

erysip <- filter(top_otu_data, classification == "Erysipelotrichaceae_unclassified")
kruskal.test(mean_rel_abund_c ~ partition, data = erysip)
dunn_test(mean_rel_abund_c ~ partition, p.adjust.method = "bonferroni", data = erysip)

#Eneterobaceriaceae

entero <- filter(top_otu_data, classification == "Enterobacteriaceae_unclassified")
kruskal.test(mean_rel_abund_c ~ partition, data = entero)
dunn_test(mean_rel_abund_c ~ partition, p.adjust.method = "bonferroni", data = entero)


#Lachnospiraceae

lachno <- filter(top_otu_data, classification == "Lachnospiraceae_unclassified")
kruskal.test(mean_rel_abund_c ~ partition, data = lachno)
dunn_test(mean_rel_abund_c ~ partition, p.adjust.method = "bonferroni", data = lachno)

#Lactobacillus 

lacto <- filter(top_otu_data, classification == "Lactobacillus")
kruskal.test(mean_rel_abund_c ~ partition, data = lacto)

#Porphyromonodaceae

porphy <- filter(top_otu_data, classification == "Porphyromonadaceae_unclassified")
kruskal.test(mean_rel_abund_c ~ partition, data = porphy)
dunn_test(mean_rel_abund_c ~ partition, p.adjust.method = "bonferroni", data = porphy)

#Ruminococcacae 

rumino <- filter(top_otu_data, classification == "Ruminococcaceae_unclassified")
kruskal.test(mean_rel_abund_c ~ partition, data = rumino)
dunn_test(mean_rel_abund_c ~ partition, p.adjust.method = "bonferroni", data = rumino)

