library(ggtext)
library(here)
library(tidyverse)

dat <- read_csv(here('data', 'cfu_meta_filtered.csv')) %>%
    filter(treatment != "No_treatment_controls (sterile broth + water + mAb vehicle)")

# COLORS
# No treatment: #99004C
# Hh + control mAb + Cd: #00AFBB
# Hh + anti-p40 mAb + Cd: #E7B800

treatment_colors <- c("#00AFBB", "#E7B800")
names(treatment_colors) <- c("Hh + control mAb + Cd", "Hh + anti-p40 mAb + Cd")

mouse_order <- dat %>%
    arrange(by = 'treatment') %>%
    pull(mouse) %>%
    unique()

yrange <- c(0,2,4,6,8)

dat %>%
    mutate(mouse = factor(mouse, levels = mouse_order)) %>%
    mutate(treatment = case_when(
        treatment == "HhWT + VPI spores + anti-p40" ~ "Hh + anti-p40 mAb + Cd",
        treatment == "HhWT + VPI spores + control mAb" ~ "Hh + control mAb + Cd",
        treatment == "No_treatment_controls (sterile broth + water + mAb vehicle)" ~ "No treatment",
        TRUE ~ NA_character_
    )) %>%
    ggplot(aes(x = day, y = cfu, color = treatment)) +
    geom_line() +
    geom_point() +
    facet_wrap("mouse") +
    scale_color_manual(values = treatment_colors) +
    scale_x_continuous(breaks = c(1,2,3,5,7,9)) +
    scale_y_log10(breaks = 10^yrange,
                  labels = yrange) +
    labs(x = "Days post *C. difficile* challenge",
         y = "*C. difficile* CFU/g<br>fecal or cecal content (log<sub>10</sub>)"
         ) +
    theme_minimal() +
    theme(legend.position = 'top',
          legend.title = element_blank(),
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          plot.margin = margin(t = 0, r = 2, b = 0, l = 2, unit = "pt"))

ggsave(here('figures','FigureS6.pdf'),
       device = 'pdf', width = 6.875, height = 9.0625)
