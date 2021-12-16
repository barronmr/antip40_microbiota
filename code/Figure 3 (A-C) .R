###############################################################################
# INTRO 
# Code to generate PcoA plots of theta YC distances using fecal 16S rRNA sequences; Figure 3
#(Note: mothur analyses were completed using the following pipeline: https://github.com/barronmr/mothurPipeline) 
# 06/28/21
#Madeline R Barron 
###############################################################################

# Goals: Generate PcoA plots of theta YC distances using 16S sequences from fecal samples collec†ed at baseline, 
#14 days post -Hh colonization, and D0 C. difficile challenge (3 weeks post-mAb).

###############################################################################
# (01) LIBRARIES
###############################################################################

library(tidyverse)

###############################################################################
# (02)READ DATA 
###############################################################################

#   (A) #Read in file denoting PcA axes designations (generated in mothur) 
    axes <- read.csv("sample.feces.0.03.pick.thetayc.0.03.lt.pcoa.axes.csv")
    
#   (B) Read in metadata
    meta <- read.csv("antip40_exp1_exp2_metadata.csv")


###############################################################################
# (03) MERGE DATA
 ###############################################################################

#   (A) Join metadata with axis file 
      
    meta_axes <- inner_join(meta, axes)

###############################################################################
 # (03) SUBSET DATA (Create dataframe for each timepoint)
###############################################################################
    
#     (A) Baseline 
      n35 <-  filter(meta_axes, day == "-35",) 
      #write_csv(n35, "feces_thetayc_dn35.csv")
    
#     (B) 2 weeks post-Hh inoculation 
      n21 <-  filter(meta_axes, day == "-21") 
      #write_csv(n21, "feces_thetayc_dn21.csv")
      
#     (C) D0 C. difficile challenge (3 wks of mAb treatment--> only showing groups that were then challenged with C.difficile to match samples included in DMM analysis)
      d0 <-  filter(meta_axes, day == "0",
                    treatment != "HhWT + isotype control mAb",
                    treatment != "HhWT + anti-p40  mAb",
                    treatment != "HhWT + anti-p40  mAb",
                    treatment != "HhWT +  sterile water  + anti-p40",
                    treatment != "Sterile broth  + VPI spores + anti-p40",
                    treatment != "Sterile broth +  VPI spores + control mAb",
                    treatment != "HhWT + sterile water  + control mAb",
                    treatment != "HhWT + anti-p40 mAb")
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
          labs(x = "PCo1 (7.97%)", y = "PCo2 (4.42%)") +
          theme_classic()

n35_plot

#Save pdf of plot 
ggsave("dn35_beta.pdf", width = 4, height = 4, units = "in")

#    (B) 2 weeks post-Hh

      n21_plot <- ggplot(n21, aes(axis1, axis2, colour=treatment_no_cd)) +
        geom_point(size = 6) +
        scale_color_manual(values= c("#00AFBB", "#E7B800", "#FC4E07")) +
        xlim(-0.6, 0.6) +
        ylim(-0.6, 0.6)+
        labs(x = "PCo1 (7.97%)", y = "PCo2 (4.42%)") +
        theme_classic()

n21_plot

ggsave("dn21_beta.pdf", width = 4, height = 4, units = "in")


#    (B) D0 C. difficile spore challenge 
      d0_plot <- ggplot(d0, aes(axis1, axis2, colour=treatment_no_cd,label=mouse_id)) +
        geom_point(size = 6) +
        geom_text(aes(label=mouse_id), colour = "black")+
        scale_color_manual(values= c("#00AFBB", "#E7B800", "#FC4E07")) +
        xlim(-0.6, 0.6) +
        ylim(-0.6, 0.6)+
        labs(x = "PCo1 (7.97%)", y = "PCo2 (4.42%)") +
        theme(legend.position = "none") +
        theme_classic()

d0_plot

ggsave("d0_beta.pdf", width = 4, height = 4, units = "in")

#Note: Cross bars and legend were added and modified in Adobe illustrator, respectively. 
