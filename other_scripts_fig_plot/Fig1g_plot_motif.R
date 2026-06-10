# Motif enrichment figure (bubble plot + average bar plot) from motif_data.txt.
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(ggtext)
library(patchwork)

motif_data <- read.table("motif_data.txt", header = TRUE, sep = "\t",
                         quote = "", comment.char = "")

motif_levels <- c(
  "CTCF(C2H2_ZF)", "OCT4-SOX2-TCF-NANOG(POU,Homeobox,HMG)", "Brn1(POU,Homeobox)",
  "PU.1:IRF8(ETS:IRF)", "ERE(NR),IR3", "Oct6(POU,Homeobox)", "IRF8(IRF)",
  "Oct4(POU,Homeobox)", "Zic2(C2H2_ZF)", "RUNX(Runt)", "Fosl2(bZIP)",
  "RUNX1(Runt)", "Fra2(bZIP)", "JunB(bZIP)", "Fra1(bZIP)", "Etv2(ETS)",
  "Foxa2(Forkhead)", "Fos(bZIP)", "Fox:Ebox(Forkhead,bHLH)", "BATF(bZIP)",
  "Atf3(bZIP)", "ETS1(ETS)", "FOXM1(Forkhead)", "Sox21(HMG)", "Fli1(ETS)",
  "Sox3(HMG)", "FOXA1(Forkhead)"
)

cellline_levels <- c("HepG2", "MCF7", "A549", "lovo", "ramos", "ESC-H1",
                     "H9", "GP5D", "GM12878", "HCT116", "HLS554P")

motif_data <- motif_data %>%
  mutate(name1 = factor(name1, levels = motif_levels),
         group = factor(group, levels = cellline_levels)) %>%
  filter(!is.na(name1))

family_colors <- c(
  Forkhead = "#E8413A", POU = "#3B7FC4", Runt = "#F08A24", bZIP = "#2CA089",
  ETS = "#EB4FB0", HMG = "#21B0C2", IRF = "#F8766D", C2H2_ZF = "#E5B80B",
  NR = "#333333", Other = "#000000"
)

get_family <- function(x) {
  case_when(
    str_detect(x, "Forkhead") ~ "Forkhead",
    str_detect(x, "POU")      ~ "POU",
    str_detect(x, "Runt")     ~ "Runt",
    str_detect(x, "bZIP")     ~ "bZIP",
    str_detect(x, "ETS")      ~ "ETS",
    str_detect(x, "HMG")      ~ "HMG",
    str_detect(x, "IRF")      ~ "IRF",
    str_detect(x, "Zf|C2H2")  ~ "C2H2_ZF",
    str_detect(x, "NR")       ~ "NR",
    TRUE                      ~ "Other"
  )
}

make_label <- function(x) {
  base  <- str_extract(x, "^[^(]+")
  paren <- str_remove(x, "^[^(]+")
  paste0(base, "<span style='color:", family_colors[get_family(x)], ";'>",
         paren, "</span>")
}

p_bubble <- ggplot(motif_data, aes(x = group, y = name1)) +
  geom_point(aes(size = Targets_Sequences, color = logP)) +
  scale_color_viridis(option = "D") +
  scale_y_discrete(labels = make_label) +
  scale_size_continuous(range = c(1, 7)) +
  labs(title = "Motif enrichment (Top5)", x = NULL, y = NULL,
       size = "% Targets", color = "-logP") +
  theme_bw() +
  theme(panel.grid  = element_blank(),
        axis.text.y = element_markdown(size = 9),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        plot.title  = element_text(hjust = 0.5),
        plot.margin = margin(5, 2, 5, 5))

bar_data <- motif_data %>%
  group_by(name1) %>%
  summarise(avg = mean(Targets_Sequences), .groups = "drop")

p_bar <- ggplot(bar_data, aes(x = avg, y = name1)) +
  geom_col(fill = "#2E8B8B", color = "black", width = 0.6) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Average", x = "% Targets", y = NULL) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title   = element_text(hjust = 0.5),
        plot.margin  = margin(5, 5, 5, 2))

final_plot <- p_bubble + p_bar +
  plot_layout(widths = c(4, 1), guides = "collect")
final_plot

ggsave("motif_enrichment.pdf", final_plot, width = 10, height = 8)
ggsave("motif_enrichment.png", final_plot, width = 10, height = 8, dpi = 300)
