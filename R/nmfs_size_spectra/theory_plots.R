# Hypothetical Syze Spectra Plots:
# About:
"
Building Plots to represent hypothetical scenarios in
size spectrum theory.

"

library(tidyverse)
library(gmRi)
library(scales)
library(geomtextpath)


# Set ggplot theme for figures
theme_set(theme_classic() + 
            theme(
              # Titles
              plot.title = element_text(hjust = 0, face = "bold", size = 14),
              plot.subtitle = element_text(size = 12),
              plot.caption = element_text(size = 7.2, margin = margin(t = 20), color = "gray40"),
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 7.5),
              # Axes
              axis.line.y = element_line(color = "black"),
              axis.ticks.y = element_line(), 
              axis.line.x = element_line(color = "black"),
              axis.ticks.x = element_line(), 
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              rect = element_rect(fill = "transparent", color = "black"),
              # Facets
              strip.text = element_text(color = "white", 
                                        face = "bold",
                                        size = 11),
              strip.background = element_rect(
                color = "#00736D", 
                fill = "#00736D", 
                size = 1, 
                linetype="solid"))
)





####  Figure 1. Base Theory  ####


# Slope Data
dat_1 <- data.frame(
  x = c(1:5),
  y = c(7,6,5,4,3),
  group = c("1"))

# Annotation Arrows
arrows <- 
  tibble(
    x1 = c(4.5, 2.25),
    x2 = c(3.8, 3.5),
    y1 = c(5.55, 1.5), 
    y2 = c(4.7, 2)
  )


# Text for annotations
top_text <- "The slope of the size spectrum\nresults from the joint change\nin abundance and size of organisms"
bot_text <- "Abundance declines with size"


# Figure
dat_1 %>% 
  ggplot(aes(x, y)) +
  geom_labelline(aes(linetype = group), label = "Size Spectrum Relationship", 
                 size = 4.5, show.legend = FALSE, label.padding = unit(0.1, "lines"),
                 linewidth = 1, hjust = 0.5) +
  annotate(x = 4.25, y = 6.15, geom = "text", label = top_text, size = 4.5, color = "gray20") +
  geom_textline(data = data.frame(x = c(1.75, 3.25), y = c(5.25, 3.75)),
                aes(x, y), linetype = 1, color = "gray20",
                label = bot_text, vjust = 1.5,
                arrow = arrow(length = unit(0.08, "inch")), size = 4.5) +
  geom_curve(
    data = arrows[1,], aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
    color = "gray20", curvature = -0.3) +
  
  # geom_curve(
  #   data = arrows[2,], aes(x = x1, y = y1, xend = x2, yend = y2),
  #   arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
  #   color = "gray20", curvature = 0.3) +
  scale_x_continuous(labels = scales::label_math(), expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(labels = scales::label_math()) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x =  element_blank(),
        axis.text.y =  element_blank(),
        axis.ticks.y =  element_blank()) +
  labs(x = "log( Organism Size )", y = "log( Abundance )")


####  Figure 2. Changing Slope:  ####

# 1. Same Intercept, Different Slope

# Slope Data
s1_dat <- data.frame(
  x = c(1:5),
  y = c(7,6,5,4,3),
  group = c("1"))

s2_dat <- data.frame(
  x = c(1:5),
  y = c(7,5.5,4,2.5,1),
  group = c("2"))

# Annotation Arrows
arrows <- 
  tibble(
    x1 = c(4.5, 2.5),
    x2 = c(4.5, 3.5),
    y1 = c(5.5, 1.5), 
    y2 = c(4, 2)
  )


# Text for annotations
top_text <- "Original Population:\nLarger individuals are abundant"
bot_text <- "Fewer larger individuals:\nSlope steepens"


# Figure
bind_rows(list(s1_dat, s2_dat)) %>% 
  ggplot(aes(x, y)) +
  geom_line(aes(linetype = group), size = 1, show.legend = FALSE) +
  annotate(x = 4.25, y = 6.15, geom = "text", label = top_text, color = "gray20", size = 4.5) +
  annotate(x = 1.75, y = 1.65, geom = "text", label = bot_text, color = "gray20", size = 4.5) +
  geom_curve(
    data = arrows[2,], aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
    color = "gray20", curvature = 0.3) +
  geom_curve(
    data = arrows[1,], aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
    color = "gray20", curvature = -0.3) +
  scale_x_continuous(labels = scales::label_math(), expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(labels = scales::label_math()) +
  labs(x = "Individual Body Mass", y = "Abundance",
       title = "Changes in Size Spectrum Slope",
       subtitle = "Relative numbers of small & large individuals shifts")




####  Figure 3: Changing Intercept:  ####



# Slope Data
int1_dat <- data.frame(
  x = c(1:5),
  y = c(7,6,5,4,3),
  group = c("1"))

int2_dat <- data.frame(
  x = c(1:5),
  y = c(5,4,3,2,1),
  group = c("2"))

# Annotation Arrows
arrows <- 
  tibble(
    x1 = c(3.65, 2.25),
    x2 = c(2.25, 3.5),
    y1 = c(6, 1.5), 
    y2 = c(6.3, 2)
  )


# Text for annotations
top_text <- "Original Population:\nHigher productivity"
bot_text <- "Lower productivity:\nLess overall biomass "


# Figure
bind_rows(list(int1_dat, int2_dat)) %>% 
  ggplot(aes(x, y)) +
  geom_line(aes(linetype = group), size = 1, show.legend = FALSE) +
  annotate(x = 4.25, y = 6, geom = "text", label = top_text, size = 4.5) +
  annotate(x = 1.5, y = 1.65, geom = "text", label = bot_text, size = 4.5) +
  geom_curve(
    data = arrows[2,], aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
    color = "gray20", curvature = 0.3) +
  geom_curve(
    data = arrows[1,], aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
    color = "gray20", curvature = 0.3) +
  scale_x_continuous(labels = scales::label_math(), expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(labels = scales::label_math()) +
  labs(x = "Individual Body Mass", y = "Abundance",
       title = "Changes in Size Spectrum Intercept",
       subtitle = "System-level changes, energy flow unchanged")





