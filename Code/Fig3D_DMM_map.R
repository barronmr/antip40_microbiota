
###############################################################################
# DMM Heatmap (Figure 3)
#01/21/21
#Roberto J. Cieza and Madeline R. Barron
###############################################################################

# Goals: Generate a heatmap illustrating relative abundances of top taxa in fecal
#samples on D0 C. difficile challenge in mice without IBD and those with active
#vs treated IBD. Annotate with treatment group, community type (as determined)
#via DMM modelling, and whether mice go on to become positive for C.difficile one day post-challenge.
#(Figure 3D)

###############################################################################
# (01) Packages

###############################################################################
############################################
# (01) LIBRARIES
###############################################################################
suppressPackageStartupMessages({
  # CRAN packages
  #######################
  library(tidyverse)            # Packages to use in everyday data analyses
  library(RColorBrewer)         # Color palettes for R
  # Bioconductor packages
  #######################
  library(ComplexHeatmap)       # Efficient to visualize associations between different sources of data sets and reveal potential patterns
})     # Efficient to visualize associations between different sources of data sets and reveal potential patterns

###############################################################################
# (02) INPUT DATA & PARAMETERS
###############################################################################

# INPUT:
# Contains METADATA: group, mouse_id, day, sex, mouse_genotype
inp_metadata  <- "data/dmm_metadata.csv"
# Contains EXPERIMENTAL information
inp_design    <- "data/genus_cdiff_only_dmm.design.csv"     # Partition information
inp_OTU_res   <- "data/genus_cdiff_only_dmm.shared.csv"         # OTU results information
# Contains TAXONOMY results
inp_taxonomy  <- "data/final.taxonomy.csv"                    # Taxonomic classification

###############################################################################
# (03) READ DATA
###############################################################################

# (03.1)
# -----------------------------------------------------------------------------
# Goal: 
# - Read CSV files containing experimental DATA

# (A) Read metadata
df_metadata <- inp_metadata %>% read_csv()

# (B) Read DESIGN data
df_design  <- inp_design  %>% read_csv()

# (C) Reading shared file of samples (those included in DMM analysis) and grouping by sample. 
df_OTU_res <- inp_OTU_res %>%
  read_csv(col_types = cols(group = col_character())) %>% 
  select(-label, -numOtus) %>%      # Get rid of label columns and numOTUs columns 
  pivot_longer(cols      = -group, 
               names_to  = "genus", 
               values_to ="count")

# (D) Read taxonomy file 
df_taxonomy <- inp_taxonomy %>% 
  read_csv() %>% 
  rename_all(tolower) %>%                                                   # Formatting COLUMN names
  mutate(taxonomy = str_replace_all(taxonomy, "(\\(\\d+\\)|\\;$)", "")) %>% # Removing confidence scores and trailing ';'
  separate(taxonomy, into = c("kingdom",                                    # Breaking tax into levels
                              "phylum",  
                              "class", 
                              "order", 
                              "family", 
                              "genus"), sep = ";") %>%                    
  gather(key   = level, 
         value = classification, 
         -otu, 
         -size)  

# (03.2)
# -----------------------------------------------------------------------------
# Goal: 
# - Merge DATA FRAMES

# (A)  Joining shared and DMM design data frames 

df_taxonomy_OTU_design <- df_OTU_res %>% 
  inner_join(df_design)

# (C) Joining taxonomy & OTU combined data frames with metadata (uses joining by group (i.e. sample) by default)
df_aggregated <- df_taxonomy_OTU_design %>% 
  inner_join(df_metadata)

df_design_meta <- inner_join(df_design, df_metadata)



###############################################################################
# (04) RELATIVE ABUNDANCE CALCULATION AND TOP OTU IDENTIFICATION
###############################################################################

# (04.1)
# -----------------------------------------------------------------------------
# Goal: 
# - Calculate relative abundances

# (A) Relative abundance 
res_relAbundance <- df_aggregated %>% 
  group_by(group) %>%                            # Groups by sample COLUMN
  mutate(total_count = sum(count)/6,             # Calculating total reads per sample and appending to df
         rel_abund = count/total_count*100) %>%  # Calculating rel abundance for each OTU and appending to df
  ungroup()

# (B) Add a COLUMN with LOG transformed relative abundance results
res_relAbundance_log <- res_relAbundance %>% 
  mutate(log_rel_abund = log10(rel_abund + 1))

# (C) Creating df of aggregated relative abundances
res_aggregated <- res_relAbundance_log %>% 
  # Aggregates classifications based on treatment
  group_by(group,                                     
           partition, 
           pos_cdiff_d1, 
           exp_id_number,
           treatment,
           genus) %>% 
  summarize(agg_log_rel_abund = sum(log_rel_abund)) %>% # Calculating aggregated rel abund, finding log 
  ungroup()

# (04.2)
# -----------------------------------------------------------------------------
# Goal: 
# - Find top taxonomic classifications

# (A) Finding top tax classifications for a given tax level
list_top_OTU_Taxo <- res_aggregated %>% 
  group_by(genus) %>%
  summarize(mean_log_rel_abund = mean(agg_log_rel_abund)) %>% # Calculating mean aggregated relative abundance across all samples
  top_n(n = 13, wt = mean_log_rel_abund) %>% 
  filter(mean_log_rel_abund > 0) %>%                     # Only choosing Otus with non-zero mean relative abundances
  pull(genus)

# (B) Generate a DATA FRAME with results for the top OTUs by Taxonomu
res_top_OTU_Taxo <- res_aggregated %>% 
  mutate(genus = case_when(!(genus %in% list_top_OTU_Taxo) ~ "Other",
                           TRUE ~ genus)) %>% 
  group_by(group,
           genus, # Need to collapse all of the newly labeled "Other" tax into single entries
           partition,
           treatment, 
           pos_cdiff_d1)%>% 
  summarize(agg_log_rel_abund = sum(agg_log_rel_abund)) %>% # Adds up all entries for "Other" classification
  ungroup()

res_top_OTU_Taxo %>% distinct(genus)
# (C) Filter out "other" classifications
res_top_OTU_Taxo <- res_top_OTU_Taxo %>% 
  filter(genus != "Other")
###############################################################################
# (05) VISUALIZATION BY HEATMAP
###############################################################################

# (05.1)
# -----------------------------------------------------------------------------
# Goal: 
# - DATA formatting for plotting

# (A) First, format data 
res_top_OTU_final_1 <- res_top_OTU_Taxo %>%
  pivot_wider(names_from  = genus, 
              values_from = agg_log_rel_abund)


# (B) Save final CSV matrix used for plotting
res_top_OTU_final_1 %>% write_csv("results/matrix_top_genus_final.csv")


res_top_OTU_final_ordered <- res_top_OTU_final_1[order(res_top_OTU_final_1$treatment),]


# (C) Place COLUMN group as rownames
res_top_OTU_rownames <- res_top_OTU_final_ordered %>% 
  column_to_rownames("group")


# (D) Creating matrix of numeric data (rel_abund of taxa)-- required for heatmap to plotclustering_distance_rows = "euclidean"
mat_top_OTU <- as.matrix(res_top_OTU_final_ordered[,5:70]) 
rownames(mat_top_OTU) = rownames(res_top_OTU_rownames) #matching rownames (groups)

# (05.2)
# -----------------------------------------------------------------------------
# Goal: 
# - Generate annotation DATA FRAME for heatmap

# (A) Creating annotation data frames 
df_annotation <- res_top_OTU_rownames %>% # creating a df with treatment matched to rows
  select(treatment,
         partition,
         pos_cdiff_d1,
         
  )
mat_top_OTU
# (05.3)
# -----------------------------------------------------------------------------
# Goal: 
# - Plot heatmap

# (A) Color_palette for all heatmaps
num_cols_mat <- max(mat_top_OTU)                                             # Maximum value in the matrix for Relative abundance
my_color_palette <- colorRampPalette(c("#E5E4E2", "#1F45FC"))(num_cols_mat)  # Palette from grays to blues

# (B) Annotation colors
annot_cols <- rowAnnotation(df  = df_annotation,
                            col = list(treatment          = c("Hh+anti-p40" = "#E7B800",
                                                              "Hh+control"   = "#00AFBB", 
                                                              "no_treatment" = "#FC4E07"),
                                       
                                       partition         = c("Partition_1" = "maroon3",   
                                                             "Partition_2"  = "#56B4E9",
                                                             "Partition_3" = "green"),
                                       
                                       pos_cdiff_d1       = c("no"      = "#332288",   
                                                              "yes"      = "#D55E00"),
                                       na_col = "black"),
                            simple_anno_size = unit(8, "mm"),
                            gp = gpar(col = "black", lty = 0.9))



annot_cols

# (C) Generating heatmap, using Euclidean distances to cluster and total df to add annotations to each sample
plot_heatmap_OTU <- 
  mat_top_OTU %>% 
  Heatmap(col                      = my_color_palette,
          right_annotation         = annot_cols,
          name                     = "log10(Relative abundance)", 
          na_col                   = "grey10",
          width                    = ncol(mat_top_OTU)*unit(4, "mm"), 
          height                   = nrow(mat_top_OTU)*unit(4, "mm"),
          border_gp                = gpar(col = "black",  lty = 1),
          rect_gp                  = gpar(col = "grey90", lwd = 2),
          row_order                = (order()),
          #row_split                = 2,
          show_row_names           = TRUE, 
          cluster_columns          = FALSE)

#clustering_distance_rows   = "euclidean") 

plot_heatmap_OTU


# (D) Save plot as PDF
pdf_dir <- file.path("results/heatmap_top_OTUs.pdf")
pdf_dir %>% pdf(width = 10, height = 14)
plot_heatmap_OTU
dev.off() 

#Note: For final figure, heatmap portion was deleted in Adobe illustrator to leave only the 
#annotation frames.Colors and legends were modified as well. 

###############################################################################
# (06) STATISTICAL ANALYSIS
###############################################################################

# (06.1)
# -----------------------------------------------------------------------------
# (A)  Is enterotype membership associated with colonization? 

fisher.test(res_top_OTU_rownames$partition, 
           res_top_OTU_rownames$pos_cdiff_d1)

# (B) What about treatment? 
fisher.test(res_top_OTU_rownames$partition, 
           res_top_OTU_rownames$treatment)


###############################################################################
# (06) CLEAN GLOBAL ENVIRONMENT
###############################################################################

# Comment out when testing script

# Remove DATA and VALUES from Global Environment
rm(list = apropos("inp_"))
rm(list = apropos("res_"))
rm(list = apropos("df_"))
rm(list = apropos("annot_"))
rm(list = apropos("plot_"))

