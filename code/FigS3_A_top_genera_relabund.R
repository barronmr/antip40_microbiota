###############################################################################
# Relative Abundance of Top 10 Genera DMM Analysis (Figure S3)
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
         rel_abund = count/total_count*100) %>%  # Calculating rel abund for each Otu and appending to df
  ungroup()


agg_otu_rel_abund <- rel_abund %>% 
  filter(level == "genus") %>% 
  group_by(group,                                     
           partition, 
           pos_cdiff_d1, 
           exp_id_number,
           treatment,
           classification) %>% # Aggregates classifications based on treatment
  summarize(agg_rel_abund = sum(rel_abund)) %>% # Calculating aggregated rel abund
  ungroup()

#Finding top tax classifications for a given tax level
top_otu_tax <- agg_otu_rel_abund %>% 
  group_by(classification) %>%
  summarize(mean_rel_abund = mean(agg_rel_abund)) %>% # Calculating mean aggregated relative abundance across all samples
  top_n(n = 10, wt = mean_rel_abund) %>% 
  filter(mean_rel_abund > 0) %>%  # Only choosing Otus with non-zero mean relative abundances
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
  summarize(agg_rel_abund = sum(agg_rel_abund)) %>% 
  mutate(log_mean_rel_abund = log10(agg_rel_abund + 1)) %>%# Adds up all entries for "Other" classification
 filter(classification != "Other") %>% 
   ungroup()

###########################################
# (04) PLOTTING 
###############################################################################

abund<- top_otu_data %>% 
  ggplot(aes(x=log_mean_rel_abund,group = partition,
                              y= classification,color = partition)) +
  geom_jitter(position=position_jitterdodge(dodge.width = 0.7,jitter.width = 0.08), size=2) +
  labs(x = "Log10 (Relative Abundance +1)",
       y = "OTU",
       color = "Community type") +
  scale_color_brewer(palette = "Set1")+
  scale_x_continuous(limits = c(0, 3)) + 
  stat_summary(aes(group=partition),fun=mean,fun.min=mean,fun.max=mean,geom="crossbar",
               width=0.5,colour="black",position=position_dodge(width=0.7)) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 11, colour = "black", angle = 90),
        text = element_text(colour = "black", size = 11))

abund + coord_flip()


###########################################
# (04) Statistics
###############################################################################

#Determining differences between enterotypes (partitions) for reach genus.

#Akkermansia

akker <- filter(top_otu_data, classification == "Akkermansia")


kruskal.test(log_mean_rel_abund ~ partition, data = akker)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = akker)

.#y.                group1      group2    n1    n2 statistic       p   p.adj p.adj.signif
#* <chr>              <chr>       <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
 # 1 log_mean_rel_abund Partition_1 Parti…    20    18     -3.58 3.46e-4 0.00104 **          
 # 2 log_mean_rel_abund Partition_1 Parti…    20     4     -3.18 1.45e-3 0.00435 **          
 # 3 log_mean_rel_abund Partition_2 Parti…    18     4     -1.05 2.92e-1 0.877   ns    

#Alistipes

alis <- filter(top_otu_data, classification == "Alistipes")

kruskal.test(log_mean_rel_abund ~ partition, data = alis)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = alis)

#  .y.                group1      group2    n1    n2 statistic       p   p.adj p.adj.signif
#* <chr>              <chr>       <chr>  <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
#  1 log_mean_rel_abund Partition_1 Parti…    20    18      4.53 5.92e-6 1.77e-5 ****        
#  2 log_mean_rel_abund Partition_1 Parti…    20     4     -1.85 6.39e-2 1.92e-1 ns          
#3 log_mean_rel_abund Partition_2 Parti…    18     4     -4.50 6.85e-6 2.05e-5 **** 

#Bacteroides

bacteroides <- filter(top_otu_data, classification == "Bacteroides")
kruskal.test(log_mean_rel_abund ~ partition, data = bacteroides)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = bacteroides)

#Not significant 

#Bifidobacterium

bifido <- filter(top_otu_data, classification == "Bifidobacterium")
kruskal.test(log_mean_rel_abund ~ partition, data = bifido)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = bifido)

#.y.                group1      group2         n1    n2 statistic       p  p.adj p.adj.signif
#* <chr>              <chr>       <chr>       <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
#  1 log_mean_rel_abund Partition_1 Partition_2    20    18   -2.60   0.00945 0.0283 *           
 # 2 log_mean_rel_abund Partition_1 Partition_3    20     4    0.0957 0.924   1      ns          
#3 log_mean_rel_abund Partition_2 Partition_3    18     4    1.62   0.105   0.316  ns  

#Erysipelotrichaceae

erysip <- filter(top_otu_data, classification == "Erysipelotrichaceae_unclassified")
kruskal.test(log_mean_rel_abund ~ partition, data = erysip)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = erysip)

# A tibble: 3 × 9
#.y.                group1      group2         n1    n2 statistic        p    p.adj p.adj.signif
#* <chr>              <chr>       <chr>       <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
  #1 log_mean_rel_abund Partition_1 Partition_2    20    18    -3.67  0.000239 0.000716 ***         
  #2 log_mean_rel_abund Partition_1 Partition_3    20     4    -3.01  0.00265  0.00794  **          
  #3 log_mean_rel_abund Partition_2 Partition_3    18     4    -0.819 0.413    1        ns  

#Eneterobaceriaceae

entero <- filter(top_otu_data, classification == "Enterobacteriaceae_unclassified")
kruskal.test(log_mean_rel_abund ~ partition, data = entero)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = entero)

#* <chr>              <chr>       <chr>       <int> <int>     <dbl>      <dbl>      <dbl> <chr>       
 # 1 log_mean_rel_abund Partition_1 Partition_2    20    18     -4.71 0.00000247 0.00000741 ****        
 # 2 log_mean_rel_abund Partition_1 Partition_3    20     4      1.56 0.120      0.360      ns          
#3 log_mean_rel_abund Partition_2 Partition_3    18     4      4.31 0.0000163  0.0000490  ****        

#Lachnospiraceae

lachno <- filter(top_otu_data, classification == "Lachnospiraceae_unclassified")
kruskal.test(log_mean_rel_abund ~ partition, data = lachno)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = lachno)

#.y.                group1      group2         n1    n2 statistic         p     p.adj p.adj.signif
#* <chr>              <chr>       <chr>       <int> <int>     <dbl>     <dbl>     <dbl> <chr>       
 # 1 log_mean_rel_abund Partition_1 Partition_2    20    18      3.96 0.0000758 0.000228  ***         
 # 2 log_mean_rel_abund Partition_1 Partition_3    20     4     -1.89 0.0588    0.176     ns          
#3 log_mean_rel_abund Partition_2 Partition_3    18     4     -4.20 0.0000269 0.0000806 ****   

#Lactobacillus 

lacto <- filter(top_otu_data, classification == "Lactobacillus")
kruskal.test(log_mean_rel_abund ~ partition, data = lacto)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = lacto)

#.y.                group1      group2         n1    n2 statistic       p   p.adj p.adj.signif
#* <chr>              <chr>       <chr>       <int> <int>     <dbl>   <dbl>   <dbl> <chr>       
 # 1 log_mean_rel_abund Partition_1 Partition_2    20    18     -2.63 0.00846 0.0254  *           
 # 2 log_mean_rel_abund Partition_1 Partition_3    20     4      1.48 0.139   0.416   ns          
#3 log_mean_rel_abund Partition_2 Partition_3    18     4      3.01 0.00257 0.00771 ** 

#Porphyromonodaceae

porphy <- filter(top_otu_data, classification == "Porphyromonadaceae_unclassified")
kruskal.test(log_mean_rel_abund ~ partition, data = porphy)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = porphy)

# A tibble: 3 × 9
#.y.                group1      group2         n1    n2 statistic        p    p.adj p.adj.signif
#* <chr>              <chr>       <chr>       <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
#  1 log_mean_rel_abund Partition_1 Partition_2    20    18      1.85 0.0638   0.191    ns          
#2 log_mean_rel_abund Partition_1 Partition_3    20     4     -2.60 0.00920  0.0276   *           
 # 3 log_mean_rel_abund Partition_2 Partition_3    18     4     -3.67 0.000242 0.000727 *** 

#Ruminococcacae 

rumino <- filter(top_otu_data, classification == "Ruminococcaceae_unclassified")
kruskal.test(log_mean_rel_abund ~ partition, data = rumino)
dunn_test(log_mean_rel_abund ~ partition, p.adjust.method = "bonferroni", data = rumino)

# A tibble: 3 × 9
.#y.                group1      group2         n1    n2 statistic        p    p.adj p.adj.signif*# <chr>              <chr>       <chr>       <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
 # 1 log_mean_rel_abund Partition_1 Partition_2    20    18      3.73 0.000188 0.000565 ***         
 # 2 log_mean_rel_abund Partition_1 Partition_3    20     4     -1.29 0.195    0.586    ns          
#3 log_mean_rel_abund Partition_2 Partition_3    18     4     -3.48 0.000505 0.00152  **  
