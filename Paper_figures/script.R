library(readxl)
library(dplyr)
library(ggplot2)
library(gridExtra)

# CLUSTERS

clusters = read.csv("clust-ent-data/clusters_combined.csv", header = T, check.names = F)
clusters_melt = reshape2::melt(clusters)
colnames(clusters_melt) = c("file_name", "data_perc", "clusters")

# #ENTITIES

entities = read.csv("clust-ent-data/entities_combined.csv", header = T, check.names = F)
entities_melt = reshape2::melt(entities)
colnames(entities_melt) = c("file_name", "data_perc", "entities")

clusters_melt$entities = entities_melt$entities 

clusts_ents = reshape::melt(clusters_melt)

# BOXPLOT

clust_ent_boxplot = ggplot(data = clusts_ents, aes(x = data_perc, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  stat_summary(fun = mean, aes(group = variable), geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  scale_fill_manual(values=c("#2D86F0", "#F02D5A"), name = "", labels = c("Clusters", "Entities")) +
  xlab("Data Percentage") +
  ylab("Number of clusters or entities") +
  ggtitle("A") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  geom_hline(yintercept= 5, linetype="dashed") +
  geom_hline(yintercept= 11, linetype="dashed") ; clust_ent_boxplot
  #guides(fill=guide_legend(title=""))

# OVERSPLITTING RATIOS

oversplitting = read.csv("match_data/oversplitting_ratio.csv", header = T, check.names = F)
oversplitting_melt = reshape2::melt(oversplitting)
colnames(oversplitting_melt) = c("file_name", "data_perc", "oversplitting_ratio")

oversplitting_excl_sing = read.csv("match_data/oversplitting_excl_singles.csv", header = T, check.names = F)
oversplitting_excl_sing_melt = reshape2::melt(oversplitting_excl_sing)
colnames(oversplitting_excl_sing_melt) = c("file_name", "data_perc", "oversplitting_ratio_excl")

oversplitting_melt$oversplitting_excl = oversplitting_excl_sing_melt$oversplitting_ratio_excl

oversplitting_combo = reshape::melt( oversplitting_melt )

# BOXPLOT

oversplitting_boxplot = ggplot(data = oversplitting_combo, aes(x = data_perc, y = value)) + 
  geom_boxplot(aes(fill = variable)) +
  stat_summary(fun = mean, aes(group = variable), geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
  scale_fill_manual(values=c("#2D86F0", "#F02D5A"), name = "", labels = c("+ singletons", "- singletons")) +
  xlab("Data Percentage") +
  ylab("Oversplitting ratio") +
  ggtitle("B") +
  theme_classic() +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(title="")) +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size=16, face = "bold")) ; oversplitting_boxplot

# PERCENTAGE SINGLETONS

singletons = read.csv("match_data/percentage_single_sample_GMYC_species.csv", header = T, check.names = F)
singletons_melt = reshape2::melt(singletons)
colnames(singletons_melt) = c("file_name", "data_perc", "perc_singletons")

singletons_boxplot = ggplot(data = singletons_melt, aes(x = data_perc, y = perc_singletons)) +
  geom_boxplot(fill = "#B5D5FA") +
  stat_summary(fun = mean, geom="point", shape=17, size=4, color="black", position=position_dodge(0.77)) +
  xlab("Data Percentage") +
  ylab("Percentage singletons") +
  ggtitle("C") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(plot.title = element_text(size=16, face = "bold")) ; singletons_boxplot

# OVERSPLITTING PER SPECIES

oversplits_per_group = read.csv("oversplits/mean_oversplits_per_group_full_100.csv", header = T, check.names = F)

oversplitting_bar = ggplot(data = oversplits_per_group, aes(x = predef_unique, y = Freq)) +
  geom_bar(stat = "summary", fill = "lightgrey", col = "black") +
  xlab("Morphospecies") +
  ylab("Mean oversplitting ratio per morphospecies group, \n on the full dataset (100%)") +
  theme_classic() +
  ggtitle("D") +
  theme(axis.title.y = element_text(size = 12, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(size = 12, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(plot.title = element_text(size=16, face = "bold")) +
  geom_hline(yintercept= 1, linetype="dashed") 

combo_plots = gridExtra::grid.arrange(clust_ent_boxplot, oversplitting_boxplot, singletons_boxplot, oversplitting_bar)
ggsave(plot = combo_plots, width = 25, height = 25, dpi = 350, filename = "combo_plots.pdf", units = "cm")
