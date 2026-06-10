library(dplyr)
library(scatterpie)
library(ggplot2)

name_map <- c(A549 = "A549", H1 = "H1-ESC", Helas3 = "Hela", LoVo = "LoVo",
              GM12878 = "GM12878", H9 = "H9-ESC", HepG2 = "HepG2", MCF7 = "MCF7",
              GP5d = "GP5d", HCT116 = "HCT116", K562 = "K562")

ord <- c("A549", "HepG2", "HCT116", "MCF7", "K562", "GM12878",
         "LoVo", "GP5d", "H1-ESC", "H9-ESC", "Hela")

overlap<-read.csv("RAD21_cellline_overlap.csv",header = F)
colnames(overlap)<-c("File1","File2","Unique1","Unique2","Overlap")
colors <- c(Unique1 = "#BEBADA", Unique2 = "#6DCBC7", Overlap = "#FFFF83")


df <- overlap %>%
  mutate(
    cl1 = name_map[sub("_RAD21.bed", "", File1, fixed = TRUE)],
    cl2 = name_map[sub("_RAD21.bed", "", File2, fixed = TRUE)],
    p1  = match(cl1, ord),
    p2  = match(cl2, ord),
    x   = pmin(p1, p2),
    y   = 12 - pmax(p1, p2),
    Horizontal = ifelse(p1 < p2, Unique1, Unique2),
    Vertical   = ifelse(p1 < p2, Unique2, Unique1),
    Common     = Overlap
  )

cprop <- df$Common / (df$Vertical + df$Horizontal + df$Common)
common_lab <- sprintf("Common (%.1f%% \u00b1 %.1f%%)", mean(cprop) * 100, sd(cprop) * 100)

p <- ggplot() +
  geom_scatterpie(aes(x = x, y = y, r = 0.45), data = df,
                  cols = c("Vertical", "Horizontal", "Common"), color = "black") +
  scale_fill_manual(
    values = c(Vertical = "#B7AED3", Horizontal = "#3FBFB0", Common = "#F4F26B"),
    labels = c(Vertical = "Vertical", Horizontal = "Horizontal", Common = common_lab)
  ) +
  scale_x_continuous(breaks = 1:10, labels = ord[1:10]) +
  scale_y_continuous(breaks = 1:10, labels = rev(ord[2:11])) +
  coord_equal() +
  labs(title = "Cross overlapping analysis of RAD21 peaks",
       x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid  = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5))

ggsave("rad21_overlap_pie.pdf", p, width = 8, height = 8)









