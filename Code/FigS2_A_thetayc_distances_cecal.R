#PCOA PLOTS OF THETAYC VALUES (CALCULATED WITH MOTHUR)
#Madeline Barron
#06/23/21

#Purpose: Generate PcoA plots of thetayc values.

#Load pre-reqs
library(tidyverse)


#Read in file denoting axes designations for each sample and metadata file
axes <- read.csv("cecal.g1.sample.final.0.03.pick.thetayc.0.03.lt.pcoa.axes.csv")
meta <- read.csv("antip40_exp1_exp2_metadata.csv")

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
                day != "9")

#PLOTTING
cecal_plot <- ggplot(cecal_g1_axes, aes(axis1, axis2, colour=treatment_no_cd)) +
  geom_point(size = 6) +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6) +
  scale_color_manual(values = c("#E69F00","#56B4E9","#009E73"))+
  print(labs(y= "PCoA2 (17.1%)", x = "PCoA1 (40.7%)")) +
  print(labs(colour = "Treatment")) +
  theme(axis.text = element_text(colour = "black"))+
  theme_classic()
  
cecal_plot

#colon contents (same procedure as above)

col_axes <- read.csv("colon.g1.sample.final.0.03.pick.thetayc.0.03.lt.pcoa.axes.csv")
meta <- read.csv("antip40_exp1_exp2_metadata.csv")

#Join metadata with axis file 
meta_col_axes <- left_join(meta,col_axes)

#Filter out cecal samples of mice harvested after 3 weeks mAb treatment (equivalent to D0 C.difficile challenge)
colon_g1_axes <- filter(meta_col_axes, sample_type == "colon.contents",
                        treatment != "HhWT + VPI spores + control mAb",
                        treatment != "Sterile broth  + VPI spores + anti-p40",
                        treatment != "HhWT + VPI spores + anti-p40",
                        treatment != "HhWT +  sterile water  + anti-p40",
                        treatment != "HhWT + sterile water  + control mAb",
                        treatment != "No_treatment_controls (sterile broth + water + mAb vehicle)",
                        treatment != "Sterile broth +  VPI spores + control mAb",
                        treatment != "HhWT + VPI spores + anti-p40 ",
                        treatment != "HhWT + VPI spores +  anti-p40 mAb")

#PLOTTING
colon_plot <- ggplot(colon_g1_axes, aes(axis1, axis2, colour=treatment, shape = treatment)) +
  geom_point(size = 6) +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5) +
  theme_classic() 
colon_plot
        
