###############################################################################
# INTRO 
# Code to generate PcoA plots of theta YC distances; Figure 3
# 12/21/21
#Madeline R Barron 
###############################################################################

# Goal: Generate PcoA plots of theta YC distances using 16S sequences from fecal samples collec†ed at baseline, 
#14 days post -Hh colonization, and D0 C. difficile challenge (3 weeks post-mAb).

###############################################################################
# (01) LIBRARIES
###############################################################################

library(tidyverse)


###############################################################################
# (02)READING AND FILTERING DATA 
###############################################################################

#Filter out timepoints of interest for use in mothur to generate PcoA axes/loadings

#   (A) shared file from whole mothur run

    shared <- read.csv('sample.final.shared.csv')


#   (B) meta_data, removing treatments and timepoints we don't want
  meta_filter<- read.csv('antip40_exp1_exp2_exp3_fmt_metadata.csv') %>% 
  filter(exp_id_number != "fmt_1",
         exp_id_number != "fmt_2",
         treatment != "HhWT +  sterile water  + anti-p40",
         treatment != "HhWT + anti-p40 ",
         treatment != "HhWT + anti-p40 mAb",
         treatment != "HhWT + anti-p40  mAb",
         treatment != "HhWT + control mAb",
         treatment != "HhWT + isotype control mAb",
         treatment != "HhWT + sterile water  + control mAb",
         treatment != "Sterile broth  + VPI spores + anti-p40",
         treatment != "Sterile broth +  VPI spores + control mAb",
         sample_type.x == "feces")


#     (C) Combining shared and meta-filtered file
  
      shared_meta_filter <- inner_join(shared,meta_filter)


#     (D) Subsetting timepoints of interest 
      
      dn35 <- shared_meta_filter %>% 
      filter(day == '-35')
      write.csv(dn35, file = "dn35_thetayc.csv")
      
  
      dn21 <- shared_meta_filter %>% 
      filter(day=="-21")
      write.csv(dn21, file = "dn21_thetayc.csv")

      d0 <- shared_meta_filter %>% 
      filter(day=='0')
      write.csv(d0, file = "d0_thetayc.csv")
      
dist <- read.csv("sample_dn35_dn21_d0.csv")

dist_dn35 <- left_join(dn35, dist)
write.csv(dist_dn35,"dn35.dist.csv")

#For mothur, copied and pasted all group_ids (sample names) into text editor, replaced spaces with dashes.

#After mothur has generated necessary files:

#    (A) Read in file denoting PcA axes designations (generated in mothur) 
axes <- read.csv("thetayc_feces_pcoa_axes.csv")

#   (B) Read in metadata
meta <- read.csv("antip40_exp1_exp2_exp3_fmt_metadata.csv")


###############################################################################
# (03) MERGE DATA
###############################################################################

#   (A) Join metadata with axis file 

meta_axes <- inner_join(meta, axes)

###############################################################################
# (03) SUBSET DATA (Create dataframe for each timepoint)
###############################################################################

#     (A) Baseline 
n35 <-  filter(meta_axes, day == "-35") 
#write_csv(n35, "feces_thetayc_dn35.csv")

#     (B) 2 weeks post-Hh inoculation 
n21 <-  filter(meta_axes, day == "-21") 
#write_csv(n21, "feces_thetayc_dn21.csv")

#     (C) D0 C. difficile challenge (3 wks of mAb treatment--> only showing groups that were then challenged with C.difficile to match samples included in DMM analysis)
d0 <-  filter(meta_axes, day == "0")

#write_csv(d0, "feces_thetayc_d0.csv")

###############################################################################
# (04) PLOTTING 
###############################################################################

#     (A) PcoA plot of baseline samples 

n35_plot <- ggplot(n35, aes(axis1, axis2, colour=treatment_no_cd)) +
  geom_point(size = 6) +
  scale_color_manual(values= c("#00AFBB", "#E7B800", "#FC4E07")) +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6)+
  labs(x = "PCo1 (21.18%)", y = "PCo2 (6.48%)") +
  theme_classic()

n35_plot

#Save as PDF

n21_pdf <- ggsave("dn35_beta.pdf", plot = n35_plot, width = 6, height = 4, units = "in")


#    (B) 2 weeks post-Hh

n21_plot <- ggplot(n21, aes(axis1, axis2, colour=treatment_no_cd)) +
  geom_point(size = 6) +
  scale_color_manual(values= c("#00AFBB", "#E7B800", "#FC4E07")) +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6)+
  labs(x = "PCo1 (21.18%)", y = "PCo2 (6.48%)") +
  theme_classic()

n21_plot

#Save as PDF
n21_pdf <- ggsave("dn21_beta.pdf", plot = n21_plot, width = 6, height = 4, units = "in")

?ggsave

#    (B) D0 C. difficile spore challenge 
d0_plot <- ggplot(d0, aes(axis1, axis2, colour=treatment_no_cd,label=mouse_id)) +
  geom_point(size = 6) +
  geom_text(aes(label=mouse_id), colour = "black")+
  scale_color_manual(values= c("#00AFBB", "#E7B800", "#FC4E07")) +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6)+
  labs(x = "PCo1 (21.18%)", y = "PCo2 (6.48%)") +
  theme(legend.position = "none") +
  theme_classic()

d0_plot

#Save as PDF
n21_pdf <- ggsave("d0_beta.pdf", plot = d0_plot, width = 6, height = 4, units = "in")

#Note: Cross bars and legend were added and modified in Adobe illustrator, respectively.

dev.off()
