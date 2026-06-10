library(dplyr)
library(ggplot2)
library(ggpubr)

all_distances <- read.table("FX1_dist_to_CTCF.txt", header = TRUE, sep = "\t")
names(all_distances)[names(all_distances) != "Group"] <- "dist_to_CTCF"

remove_outliers <- function(x) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x, na.rm = TRUE)
  x[x < qnt[1] - H | x > qnt[2] + H] <- NA
  x
}

dat <- all_distances %>%
  mutate(value = log2(dist_to_CTCF + 1)) %>%
  filter(value != 0) %>%
  group_by(Group) %>%
  mutate(value = remove_outliers(value)) %>%
  ungroup() %>%
  filter(!is.na(value)) %>%
  mutate(Group = factor(Group, levels = c("Common", "R219S_specific", "WT_specific")))

n_lab <- dat %>% count(Group) %>% mutate(label = paste0("n=", n))

cols <- c(Common = "#A0A0A0", R219S_specific = "#E87D96", WT_specific = "#7C97C9")

my_comparisons <- list(c("Common", "R219S_specific"),
                       c("Common", "WT_specific"),
                       c("R219S_specific", "WT_specific"))

p <- ggplot(dat, aes(Group, value, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.3) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = cols) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  geom_text(data = n_lab, aes(x = Group, y = 12, label = label),
            inherit.aes = FALSE, size = 4) +
  labs(x = NULL, y = "Log2(min.distance to CTCF sites + 1)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("FX1_dist_to_CTCF_boxplot.pdf", p, width = 4.5, height = 5.5)
