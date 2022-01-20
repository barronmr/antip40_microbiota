################################################################################
# INTRO 
# Code to generate PcoA plots of theta YC distances in FMT recipients  (Fig 4C)
# 12/21/21
#Madeline R Barron 
###############################################################################

# Goal: Generate PcoA plots of theta YC distances using 16S sequences from fecal samples collec†ed from donors on 
#d from recipients 7 days post-FM, colored based on whether become positive or negative for 
#C.difficile

##############################################################################
# (01) LIBRARIES
###############################################################################

library(tidyverse)

##############################################################################
# (02)READING AND JOINING DATA
###############################################################################

meta <- read.csv("data/antip40_exp1_exp2_exp3_fmt_metadata.csv")
d0_fmt <- read.csv("data/recipients_d0.thetayc.0.03.pcoa.axes.csv")
d0_fmt_meta <- inner_join(d0_fmt, meta)


##############################################################################
# (03) PLOTTING
###############################################################################

colors_2 <- c("yes" = "black", 
              "no" = "grey")

d0_plot <- ggplot(d0_fmt_meta, aes(axis1, axis2, color  = ever_pos_cdiff)) +
  geom_point(size = 6, shape = 17) +
  scale_color_manual(values=colors_2)+
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6)+
  labs(x = "PCo1 (32.41%)", y = "PCo2 (18.53%)") +
  theme_classic()

d0_plot

ggsave("results/fmt_d0.pdf", height = 5, width = 8, units = "in")
