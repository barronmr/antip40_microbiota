###############################################################################
# INTRO 
# Code to generate PcoA plots of theta YC distances, cecal contents; Figure S4
# 06/23/22
#Madeline R Barron 
###############################################################################

# Goal: Generate PcoA plots of theta YC distances using 16S sequences from cecal contents harvested
# on the equivalent of D0 C. difficile challenge (3 weeks post-mAb).

###############################################################################
# (01) LIBRARIES
##############################################################################
library(tidyverse)

###############################################################################
# (02)READING AND FILTERING DATA 
###############################################################################
axes <- read.csv("cecal.thetayc.pcoa.axes.csv")
meta <- read.csv("antip40_exp1_exp2_exp3_fmt_metadata.csv")

#Join metadata with axis file 
meta_axes <- left_join(meta, axes)

#Filter out cecal samples of mice harvested after 3 weeks mAb treatment (equivalent to D0 C.difficile challenge)
cecal_g1_axes <- filter(meta_axes, sample_type.x == "cecal.contents",
                treatment != "HhWT + VPI spores + control mAb",
                treatment != "Sterile broth  + VPI spores + anti-p40",
                treatment != "HhWT + VPI spores + anti-p40",
                treatment != "HhWT +  sterile water  + anti-p40",
                treatment != "HhWT + sterile water  + control mAb",
                treatment != "Sterile broth +  VPI spores + control mAb",
                treatment != "HhWT + VPI spores + anti-p40 ",
                treatment != "HhWT + VPI spores +  anti-p40 mAb",
                exp_id_number == "p40_2",
                day != "9")

#PLOTTING
cecal_plot <- ggplot(cecal_g1_axes, aes(axis1, axis2, colour=treatment_no_cd)) +
  geom_point(size = 6) +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6) +
  scale_color_manual(values = c("#009E73","#E69F00","#56B4E9"))+
  print(labs(y= "PCoA2 (17.1%)", x = "PCoA1 (40.7%)")) +
  print(labs(colour = "Treatment")) +
  theme(axis.text = element_text(colour = "black"))+
  theme_classic()
  
cecal_plot
