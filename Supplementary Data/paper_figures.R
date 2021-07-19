library(readxl)
library(dplyr)
library(ggplot2)
library(gridExtra)

########################################################################################################
# READ IN ALL DATA FIRST
########################################################################################################


########################################################################################################
                                                    #PATHD8
########################################################################################################

##########################################################
# OVERSPLITTING RATIO PER SPECIES GROUP ON FULL DATASET
##########################################################

oversplitting_full_data = read.csv("PATHD8/mean_oversplits_per_group_100_percent_data.csv")
oversplitting_full_data$predefined_group = as.factor(oversplitting_full_data$predefined_group)

oversplitting_pathd8 = ggplot2::ggplot(data = oversplitting_full_data, aes(x = predefined_group, y = mean)) + 
  geom_bar(stat = "identity", fill = "lightblue", colour = "black") +
  theme_classic() +
  ggtitle("A)", subtitle = "PATHD8") +
  ylab("Mean oversplitting ratio") +
  xlab("Predefined morphological species group") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(breaks = seq(0, 18, by = 2), limits = c(0, 18)) 
  #theme(legend.position = "none")
  

################################
# PERCENTAGE MATCHES EXCL SINGLES
################################

per_match_pathd8 = read.csv("PATHD8/percentage_match_excl_singles.csv", header = T, check.names = F)
per_match_pathd8 = reshape2::melt(per_match_pathd8)
per_match_pathd8$variable = as.factor(per_match_pathd8$variable)
str(per_match_pathd8)

###################################
# PERCENTAGE MATCHES INCL SINGLES
###################################

per_match_inc_singles_pathd8 = read.csv("PATHD8/percentage_match.csv", check.names = F)
per_match_inc_singles_pathd8 = reshape2::melt(per_match_inc_singles_pathd8)
per_match_inc_singles_pathd8$variable = as.factor(per_match_inc_singles_pathd8$variable)

################################
# ENTITIES
################################

entities_pathd8 = read.csv("PATHD8/entities.csv", header = T, check.names = FALSE)
entities_pathd8 = reshape2::melt(entities_pathd8)
entities_pathd8$variable = as.factor(entities_pathd8$variable)

################################
# CLUSTERS
################################

clusters_pathd8 = read.csv("PATHD8/clusters.csv", header = T, check.names = FALSE)
clusters_pathd8 = reshape2::melt(clusters_pathd8)
clusters_pathd8$variable = as.factor(clusters_pathd8$variable)

###############################
# PERCENTAGE SINGLETON GMYC SPECIES
##############################

singletons_pathd8 = read.csv("PATHD8/percentage_single_sample_GMYC_species.csv", header = T, check.names = FALSE)
singletons_pathd8 = reshape2::melt(singletons_pathd8)
singletons_pathd8$variable = as.factor(singletons_pathd8$variable)

########################################################################################################
                                                # CHRONOS (APE)
########################################################################################################
##########################################################
# OVERSPLITTING RATIO PER SPECIES GROUP ON FULL DATASET
##########################################################

#lambda 1

oversplitting_full_data_chronos_lambda1 = read.csv("chronos/mean_oversplits_per_group_100_lambda1.csv")
oversplitting_full_data_chronos_lambda1$predefined_group = as.factor(oversplitting_full_data_chronos_lambda1$predefined_group)

oversplitting_chronos_lambda1 = ggplot2::ggplot(data = oversplitting_full_data_chronos_lambda1, aes(x = predefined_group, y = mean)) + 
  geom_bar(stat = "identity", fill = "lightblue", colour = "black") +
  theme_classic() +
  ggtitle("C)", subtitle = "chronos, λ = 1") +
  ylab("Mean oversplitting ratio") +
  xlab("Predefined morphological species group") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(breaks = seq(0, 18, by = 2), limits = c(0, 18)) 
#theme(legend.position = "none")

#lambda 0

oversplitting_full_data_chronos_lambda0 = read.csv("chronos/mean_oversplits_per_group_100_lambda0.csv")
oversplitting_full_data_chronos_lambda0$predefined_group = as.factor(oversplitting_full_data_chronos_lambda0$predefined_group)

oversplitting_chronos_lambda0 = ggplot2::ggplot(data = oversplitting_full_data_chronos_lambda0, aes(x = predefined_group, y = mean)) + 
  geom_bar(stat = "identity", fill = "lightblue", colour = "black") +
  theme_classic() +
  ggtitle("B)", subtitle = "chronos, λ = 0") +
  ylab("Mean oversplitting ratio") +
  xlab("Predefined morphological species group") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(breaks = seq(0, 18, by = 2), limits = c(0, 18)) 

#lambda 10

oversplitting_full_data_chronos_lambda10 = read.csv("chronos/mean_oversplits_per_group_100_lambda10.csv")
oversplitting_full_data_chronos_lambda10$predefined_group = as.factor(oversplitting_full_data_chronos_lambda10$predefined_group)

oversplitting_chronos_lambda10 = ggplot2::ggplot(data = oversplitting_full_data_chronos_lambda10, aes(x = predefined_group, y = mean)) + 
  geom_bar(stat = "identity", fill = "lightblue", colour = "black") +
  theme_classic() +
  ggtitle("D)", subtitle = "chronos, λ = 10") +
  ylab("Mean oversplitting ratio") +
  xlab("Predefined morphological species group") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(breaks = seq(0, 18, by = 2), limits = c(0, 18)) 

oversplitting_bars = gridExtra::grid.arrange(oversplitting_pathd8, 
                                             oversplitting_chronos_lambda0, 
                                             oversplitting_chronos_lambda1, 
                                             oversplitting_chronos_lambda10,
                                             ncol = 4)

ggsave(plot = oversplitting_bars, width = 35, height = 15, dpi = 350, filename = "oversplitting_bars2.svg", units = "cm")

###################################
# PERCENTAGE MATCHES EXCL SINGLES
###################################

sht = excel_sheets("chronos/percentage_matches_excl_singles.xlsx")
chronos_per_match = lapply(setNames(sht, sht), function(s) read_excel("chronos/percentage_matches_excl_singles.xlsx", sheet=s))
chronos_per_match = bind_rows(chronos_per_match, .id="Sheet")
chronos_per_match = reshape2::melt(chronos_per_match)


###################################
# PERCENTAGE MATCHES INCL SINGLES
###################################

sht_inc_singles = excel_sheets("chronos/percentage_matches.xlsx")
chronos_per_match_inc = lapply(setNames(sht_inc_singles, sht_inc_singles), function(s) read_excel("chronos/percentage_matches.xlsx", sheet=s))
chronos_per_match_inc = bind_rows(chronos_per_match_inc, .id="Sheet")
chronos_per_match_inc = reshape2::melt(chronos_per_match_inc)

################################
# ENTITIES
################################

sht_entities = excel_sheets("chronos/entities.xlsx")
chronos_entities = lapply(setNames(sht_entities, sht_entities), function(s) read_excel("chronos/entities.xlsx", sheet=s))
chronos_entities = bind_rows(chronos_entities, .id="Sheet")
chronos_entities = reshape2::melt(chronos_entities)

################################
# CLUSTERS
################################

sht_clusters = excel_sheets("chronos/clusters.xlsx")
chronos_clusters = lapply(setNames(sht_clusters, sht_clusters), function(s) read_excel("chronos/clusters.xlsx", sheet=s))
chronos_clusters = bind_rows(chronos_clusters, .id="Sheet")
chronos_clusters = reshape2::melt(chronos_clusters)

################################
# PERCENTAGE GMYC SINGLETONS
################################

sht_singletons = excel_sheets("chronos/percentage_gmyc_singleton_species.xlsx")
chronos_singletons = lapply(setNames(sht_singletons, sht_singletons), function(s) read_excel("chronos/percentage_gmyc_singleton_species.xlsx", sheet=s))
chronos_singletons = bind_rows(chronos_singletons, .id="Sheet")
chronos_singletons = reshape2::melt(chronos_singletons)


########################################################################################################
                                                 # PLOTS
########################################################################################################


# COMBINATION LINE PLOT FOR PERCENTAGE MATCHES, EXLUDING SINGLETONS, ACROSS DATA PERCENTAGES AND 
# ULTRAMETRIC METHODS


per_match_plot = ggplot(per_match_pathd8, aes(x = as.numeric(variable)*10, y = value)) +
  geom_smooth(method = "lm", alpha = 0.1, aes(colour = "PATHD8")) +
  # add in lambda values:
  geom_smooth(  data=subset(chronos_per_match,Sheet=="lambda_0"), aes(colour = "λ = 0"), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_per_match,Sheet=="lambda_1"), aes(colour = "λ = 1"), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_per_match,Sheet=="lambda_10"), aes(colour = "λ = 10"), alpha = 0.1, method = "lm" ) +
  scale_colour_manual(name="Ultrametric method", values=c("blue", "red", "green", "orange")) +
  xlab("Data percentage") +
  ylab("Percentage match excluding singletons") + 
  ggtitle("The percentage match (exluding singletons) for different ultrametric methods") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  theme(legend.position = "bottom") ;per_match_plot

########################################################################################################

# COMBINATION LINE PLOT FOR PERCENTAGE MATCHES, INLUDING SINGLETONS, ACROSS DATA PERCENTAGES AND 
# ULTRAMETRIC METHODS


per_match_plot_inc = ggplot(per_match_inc_singles_pathd8, aes(x = as.numeric(variable)*10, y = value)) +
  geom_smooth(method = "lm", alpha = 0.1, aes(colour = "PATHD8")) +
  # add in lambda values:
  geom_smooth(  data=subset(chronos_per_match_inc,Sheet=="lambda_0"), aes(colour = "λ = 0"), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_per_match_inc,Sheet=="lambda_1"), aes(colour = "λ = 1"), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_per_match_inc,Sheet=="lambda_10"), aes(colour = "λ = 10"), alpha = 0.1, method = "lm" ) +
  scale_colour_manual(name="Ultrametric method", values=c("blue", "red", "green", "orange")) +
  xlab("Data percentage") +
  ylab("Percentage match including singletons") + 
  ggtitle("The percentage match (inluding singeltons) for different ultrametric methods") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  theme(legend.position = "bottom") ;per_match_plot_inc

########################################################################################################

# COMBINATION LINE PLOT FOR THE NUMBER OF CLUSTERS ACROSS DATA PERCENTAGES AND 
# ULTRAMETRIC METHODS

clusters_pathd8_plot = ggplot(clusters_pathd8, aes(x = variable, y = value)) +
  geom_smooth(method = "lm", alpha = 0.1, aes(colour = "PATHD8", group=1)) +
  # add in lambda values:
  geom_smooth(  data=subset(chronos_clusters,Sheet=="lambda_0"), aes(colour = "λ = 0", group=1), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_clusters,Sheet=="lambda_1"), aes(colour = "λ = 1", group=1), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_clusters,Sheet=="lambda_10"), aes(colour = "λ = 10", group=1), alpha = 0.1, method = "lm" ) +
  scale_colour_manual(name="Ultrametric method", values=c("blue", "red", "green", "orange")) +
  geom_hline(yintercept= 5, linetype="dashed", color = "red") +
  geom_hline(yintercept= 11, linetype="dashed", color = "blue") +
  xlab("Data percentage") +
  ylab("Clusters") + 
  ggtitle("A") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 30, by = 2)) +
  coord_cartesian(ylim = c(0, 25)) +
  theme(legend.position = "bottom") ; clusters_pathd8_plot

########################################################################################################

# COMBINATION LINE PLOT FOR THE NUMBER OF ENTITIES, EXLUDING SINGLETONS, ACROSS DATA PERCENTAGES AND 
# ULTRAMETRIC METHODS

entities_pathd8_plot = ggplot(entities_pathd8, aes(x = variable, y = value)) +
  geom_smooth(method = "lm", alpha = 0.1, aes(colour = "PATHD8", group=1)) +
  # add in lambda values:
  geom_smooth(  data=subset(chronos_entities,Sheet=="lambda_0"), aes(colour = "λ = 0", group=1), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_entities,Sheet=="lambda_1"), aes(colour = "λ = 1", group=1), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_entities,Sheet=="lambda_10"), aes(colour = "λ = 10", group=1), alpha = 0.1, method = "lm" ) +
  scale_colour_manual(name="Ultrametric method", values=c("blue", "red", "green", "orange")) +
  geom_hline(yintercept= 5, linetype="dashed", color = "red") +
  geom_hline(yintercept= 11, linetype="dashed", color = "blue") +
  xlab("Data percentage") +
  ylab("Entities") + 
  ggtitle("B") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 60, by = 5)) +
  coord_cartesian(ylim = c(0, 60)) +
  theme(legend.position = "bottom") ; entities_pathd8_plot

p1 = gridExtra::grid.arrange(clusters_pathd8_plot, entities_pathd8_plot, ncol = 2)
ggsave(plot = p1, width = 25, height = 15, dpi = 350, filename = "clusts_ents.png", units = "cm")


#############################################################################
# PLOT FOR PERCENTAGE MATCH DATA WITH AND WITHOUT SINGLETONS ALL ON ONE GRAPH
############################################################################
# add group = 1 to get the geom_smooth lines to show!!

per_match_plot_all = ggplot(per_match_pathd8, aes(x = variable, y = value)) +
  #geom_point() +
  geom_smooth(method = "lm", alpha = 0.1, aes(colour = "PATHD8", group=1)) +
  geom_smooth(data = per_match_inc_singles_pathd8, method = "lm", alpha = 0.1, aes(colour = "PATHD8", group=1), linetype = "dashed") +
  # add in lambda values:
  geom_smooth(  data=subset(chronos_per_match,Sheet=="lambda_0"), aes(colour = "Chronos λ = 0", group=1), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_per_match,Sheet=="lambda_1"), aes(colour = "Chronos λ = 1", group=1), alpha = 0.1, method = "lm" ) +
  geom_smooth(  data=subset(chronos_per_match,Sheet=="lambda_10"), aes(colour = "Chronos λ = 10", group=1), alpha = 0.1, method = "lm" ) +
  
  geom_smooth(  data=subset(chronos_per_match_inc,Sheet=="lambda_0"), aes(colour = "Chronos λ = 0", group=1), alpha = 0.1, method = "lm", linetype = "dashed" ) +
  geom_smooth(  data=subset(chronos_per_match_inc,Sheet=="lambda_1"), aes(colour = "Chronos λ = 1", group=1), alpha = 0.1, method = "lm", linetype = "dashed" ) +
  geom_smooth(  data=subset(chronos_per_match_inc,Sheet=="lambda_10"), aes(colour = "Chronos λ = 10", group=1), alpha = 0.1, method = "lm", linetype = "dashed" ) +
  
  scale_colour_manual(name="Ultrametric method", values=c("blue", "red", "green", "orange")) +
  xlab("Data percentage") +
  ylab("Percentage match") +
  labs(subtitle = "                ____ Excluding singletons      _ _ _ Including singletons") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  theme(legend.position = "top") ;per_match_plot_all

ggsave(plot = per_match_plot_all, width = 20, height = 15, dpi = 350, filename = "per_match_plot_all.png", units = "cm")


########################################################################
                      # SUBSET JUST 100% OF THE DATA
########################################################################

#########################################
# PERCENTAGE MATCHES EXCLUDING SINGLETONS
#########################################

hundred_perc_data = c() 

hundred_perc_data$pathd8 = subset(per_match_pathd8, variable == "100")$value
min(hundred_perc_data$pathd8)
max(hundred_perc_data$pathd8)

h1 = hist(hundred_perc_data$pathd8, main = "", xlab = "", ylab = "", ylim = c(0,100), xlim = c(50,100))


hundred_perc_data$chronos_lambda0 = subset(chronos_per_match, Sheet=="lambda_0" & variable == "100")$value
min(hundred_perc_data$chronos_lambda0)
max(hundred_perc_data$chronos_lambda0)

h2 = hist(hundred_perc_data$chronos_lambda0, main = "", xlab = "", ylab = "", ylim = c(0,100), xlim = c(50,100))

hundred_perc_data$chronos_lambda1 = subset(chronos_per_match, Sheet=="lambda_1" & variable == "100")$value
min(hundred_perc_data$chronos_lambda1)
max(hundred_perc_data$chronos_lambda1)

h3 = hist(hundred_perc_data$chronos_lambda1, main = "", xlab = "", ylab = "", ylim = c(0,100), xlim = c(50,100))

hundred_perc_data$chronos_lambda10 = subset(chronos_per_match, Sheet=="lambda_10" & variable == "100")$value
min(hundred_perc_data$chronos_lambda10 )
max(hundred_perc_data$chronos_lambda10 )

h4 = hist(hundred_perc_data$chronos_lambda10, main = "", xlab = "", ylab = "", ylim = c(0,100), xlim = c(50,100))


hundred_perc_data_melt = reshape::melt(hundred_perc_data)

#########################################
# PERCENTAGE MATCHES INCLUDING SINGLETONS
#########################################

hundred_perc_data_inc_singles = c() 

hundred_perc_data_inc_singles$pathd8 = subset(per_match_inc_singles_pathd8, variable == "100")$value

min(hundred_perc_data_inc_singles$pathd8)
max(hundred_perc_data_inc_singles$pathd8)

h1 = hist(hundred_perc_data_inc_singles$pathd8, main = "", xlab = "", ylab = "", ylim = c(0,100), xlim = c(50,100))

hundred_perc_data_inc_singles$chronos_lambda0 = subset(chronos_per_match_inc, Sheet=="lambda_0" & variable == "100")$value
min(hundred_perc_data_inc_singles$chronos_lambda0)
max(hundred_perc_data_inc_singles$chronos_lambda0)

h2 = hist(hundred_perc_data_inc_singles$chronos_lambda0, main = "", xlab = "", ylab = "",  ylim = c(0,100), xlim = c(50,100))

hundred_perc_data_inc_singles$chronos_lambda1 = subset(chronos_per_match_inc, Sheet=="lambda_1" & variable == "100")$value
min(hundred_perc_data_inc_singles$chronos_lambda1)
max(hundred_perc_data_inc_singles$chronos_lambda1)

h3 = hist(hundred_perc_data_inc_singles$chronos_lambda1, main = "", xlab = "", ylab = "", ylim = c(0,100), xlim = c(50,100))

hundred_perc_data_inc_singles$chronos_lambda10 = subset(chronos_per_match_inc, Sheet=="lambda_10" & variable == "100")$value
min(hundred_perc_data_inc_singles$chronos_lambda10)
max(hundred_perc_data_inc_singles$chronos_lambda10)

h3 = hist(hundred_perc_data_inc_singles$chronos_lambda10, main = "", xlab = "", ylab = "", ylim = c(0,100), xlim = c(50,100))

hundred_perc_data_inc_singles_melt = reshape::melt(hundred_perc_data_inc_singles)

#########################################
# MERGE THE PERCENTAGE MATCH DATA
#########################################

merged_per_data = c()
merged_per_data$inc_singles = hundred_perc_data_inc_singles_melt$value
merged_per_data$ex_singles = hundred_perc_data_melt$value
merged_per_data$method = hundred_perc_data_melt$L1
merged_per_data = as.data.frame(merged_per_data)
merged_per_data_melt = reshape::melt(merged_per_data)

#########################################
# CLUSTERS AND ENTITIES
#########################################

#########################################
# ENTITIES
#########################################

hundred_perc_data_entities = c()

hundred_perc_data_entities$pathd8 = subset(entities_pathd8, variable == "100")$value
max(hundred_perc_data_entities$pathd8)
min(hundred_perc_data_entities$pathd8)

h1 = hist(hundred_perc_data_entities$pathd8, main = "", xlab = "", ylab = "", ylim = c(0,80), xlim = c(0,100))

hundred_perc_data_entities$chronos_lambda0 = subset(chronos_entities, Sheet=="lambda_0" & variable == "100")$value
min(hundred_perc_data_entities$chronos_lambda0)
max(hundred_perc_data_entities$chronos_lambda0)

h2 = hist(hundred_perc_data_entities$chronos_lambda0, main = "", xlab = "", ylab = "", ylim = c(0,80), xlim = c(0,100))

hundred_perc_data_entities$chronos_lambda1 = subset(chronos_entities, Sheet=="lambda_1" & variable == "100")$value
min(hundred_perc_data_entities$chronos_lambda1)
max(hundred_perc_data_entities$chronos_lambda1)

h3 = hist(hundred_perc_data_entities$chronos_lambda1, main = "", xlab = "", ylab = "", ylim = c(0,80), xlim = c(0,100))

hundred_perc_data_entities$chronos_lambda10 = subset(chronos_entities, Sheet=="lambda_10" & variable == "100")$value
min(hundred_perc_data_entities$chronos_lambda10)
max(hundred_perc_data_entities$chronos_lambda10)

h4 = hist(hundred_perc_data_entities$chronos_lambda10, main = "", xlab = "", ylab = "", ylim = c(0,80), xlim = c(0,100))

par(mfrow=c(2,2))

hundred_perc_data_entities_melt = reshape::melt(hundred_perc_data_entities)

#########################################
# CLUSTERS 
#########################################

hundred_perc_data_clusters = c() # excludes singletons

hundred_perc_data_clusters$pathd8 = subset(clusters_pathd8, variable == "100")$value
min(hundred_perc_data_clusters$pathd8)
max(hundred_perc_data_clusters$pathd8)

h1 = hist(hundred_perc_data_clusters$pathd8 , ylim = c(0,80), xlim = c(0,50), main = "PATHD8", xlab = "Distribution of delimited GMYC species", ylab = "Frequency")

hundred_perc_data_clusters$chronos_lambda0 = subset(chronos_clusters, Sheet=="lambda_0" & variable == "100")$value
min(hundred_perc_data_clusters$chronos_lambda0)
max(hundred_perc_data_clusters$chronos_lambda0)

h2 = hist(hundred_perc_data_clusters$chronos_lambda0, ylim = c(0,80), xlim = c(0,50), main = "chronos, λ = 0", xlab = "", ylab = "") 

hundred_perc_data_clusters$chronos_lambda1 = subset(chronos_clusters, Sheet=="lambda_1" & variable == "100")$value
min(hundred_perc_data_clusters$chronos_lambda1)
max(hundred_perc_data_clusters$chronos_lambda1)

h3 = hist(hundred_perc_data_clusters$chronos_lambda1, ylim = c(0,80), xlim = c(0,50), main = "chronos, λ = 1", xlab = "", ylab = "") 

hundred_perc_data_clusters$chronos_lambda10 = subset(chronos_clusters, Sheet=="lambda_10" & variable == "100")$value
min(hundred_perc_data_clusters$chronos_lambda10)
max(hundred_perc_data_clusters$chronos_lambda10)

h4 = hist(hundred_perc_data_clusters$chronos_lambda10, ylim = c(0,80), xlim = c(0,50), main = "chronos, λ = 10", xlab = "", ylab = "") 

hundred_perc_data_clusters_melt = reshape::melt(hundred_perc_data_clusters )
hundred_perc_data_clusters_melt =  as.data.frame(hundred_perc_data_clusters_melt)

#########################################
# MERGE CLUSTERS AND ENTITIES
#########################################

merged_clust_ent_data = c()
merged_clust_ent_data$clusters = hundred_perc_data_clusters_melt$value
merged_clust_ent_data$entities = hundred_perc_data_entities_melt$value
merged_clust_ent_data$method = hundred_perc_data_melt$L1
merged_clust_ent_data = as.data.frame(merged_clust_ent_data)
merged_clust_ent_data = reshape::melt(merged_clust_ent_data)



#########################################
# BOXPLOTS
#########################################
new_x_labs = c("chronos λ = 0", "chronos λ = 1", "chronos λ = 10", "PATHD8")
#########################################
# CLUSTERS AND ENTITIES
#########################################

clust_ent_plot = ggplot(merged_clust_ent_data, aes(x = method, y = value, fill = variable)) +
  scale_fill_manual(values = c("lightpink", "lightblue"), name = "", labels = c("Clusters", "Entities")) +
  geom_boxplot() +
  xlab("Ultrametric method") +
  ylab("Number of clusters or entities") + 
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits=c(0,100)) +
  ggtitle("A") +
  labs(subtitle = "100% of the data") +
  theme_classic() +
  #theme(plot.title = element_text(size=16, face = "bold"), axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept= 5, linetype="dashed", color = "red") +
  geom_hline(yintercept= 11, linetype="dashed", color = "blue") +
  scale_x_discrete(labels= new_x_labs) ; clust_ent_plot


#########################################
# PERCENTAGE MATCHES
#########################################

per_match_compare1 = ggplot(merged_per_data_melt, aes(x = method, y = value, fill = variable)) +
  scale_fill_manual(values = c("lightpink", "lightblue"), name = "", labels = c("Including singletons", "Excluding singletons")) +
  geom_boxplot() +
  xlab("Ultrametric method") +
  ylab("Percentage Match") + 
  ggtitle("B") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits=c(0,100)) +
  labs(subtitle = "100% of the data") +
  theme_classic() +
  #theme(plot.title = element_text(size=16, face = "bold"), axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels= new_x_labs) ; per_match_compare1 



combo_boxplots = gridExtra::grid.arrange(clust_ent_plot, per_match_compare1, ncol = 2)
ggsave(plot = combo_boxplots, width = 25, height = 15, dpi = 350, filename = "combo_boxplots.svg", units = "cm")

#########################################
# INDIVIDUAL BOXPLOTS
#########################################

# Percentage match excluding singletons for PATHD8-----------------------------------------

per_match_pathd8_boxplot = ggplot(per_match_pathd8, aes(x = variable, y = value)) +
  geom_boxplot() +
  xlab("Data percentage") +
  ylab("Percentage Match") + 
  labs(subtitle = "Percentage Match Excluding Singletons, PATHD8") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  #ggtitle("A") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") 

# Percentage match including singletons for PATHD8-----------------------------------------

per_match_inc_pathd8_boxplot = ggplot(per_match_inc_singles_pathd8, aes(x = variable, y = value)) +
  geom_boxplot() +
  xlab("Data percentage") +
  ylab("Percentage Match") + 
  labs(subtitle = "Percentage Match Including Singletons, PATHD8") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  #ggtitle("A") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;per_match_inc_pathd8_boxplot

gridExtra::grid.arrange(per_match_pathd8_boxplot, per_match_inc_pathd8_boxplot)

# Number of entities for PATHD8-----------------------------------------

entities_pathd8_boxplot = ggplot(entities_pathd8, aes(x = variable, y = value)) +
  geom_boxplot() +
  xlab("Data percentage") +
  ylab("Number of entities") + 
  labs(subtitle = "Number of entities, PATHD8") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 60, by = 5)) +
  #ggtitle("A") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;entities_pathd8_boxplot

# Number of clusters for PATHD8-----------------------------------------

clusters_pathd8_boxplot = ggplot(clusters_pathd8, aes(x = variable, y = value)) +
  geom_boxplot() +
  xlab("Data percentage") +
  ylab("Number of clusters") + 
  labs(subtitle = "Number of clusters, PATHD8") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  #ggtitle("A") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;clusters_pathd8_boxplot

# Percentage GMYC singletons for PATHD8-----------------------------------------

singletons_pathd8_boxplot = ggplot(singletons_pathd8, aes(x = variable, y = value)) +
  geom_boxplot() +
  xlab("Data percentage") +
  ylab("Percentage GMYC singletons") + 
  labs(subtitle = "A) PATHD8") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  expand_limits(x = 0, y = c(0,100)) +
  #ggtitle("A") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;singletons_pathd8_boxplot



# Percentage match excluding singletons for chronos -----------------------------------------

per_matches_chronos_plot_0 = ggplot(chronos_per_match, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_per_match, Sheet=="lambda_0")) +
  xlab("Data percentage") + 
  ylab("Percentage matches") + 
  #ggtitle("C") + 
  labs(subtitle = "Percentage Match Excluding Singletons, chronos(), lambda = 0") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(legend.position = "none") ; per_matches_chronos_plot_0


per_matches_chronos_plot_1 = ggplot(chronos_per_match, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_per_match, Sheet=="lambda_1")) +
  xlab("Data percentage") + 
  ylab("Percentage matches") + 
  #ggtitle("E") + 
  labs(subtitle = "Percentage Match Excluding Singletons, chronos(), lambda = 1") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(legend.position = "none") ; per_matches_chronos_plot_1


per_matches_chronos_plot_10 = ggplot(chronos_per_match, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_per_match, Sheet=="lambda_10")) +
  xlab("Data percentage") + 
  ylab("Percentage matches") + 
  #ggtitle("G") + 
  labs(subtitle = "Percentage Match Excluding Singletons, chronos(), lambda = 10") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(legend.position = "none") ;per_matches_chronos_plot_10


gridExtra::grid.arrange(per_matches_chronos_plot_0, per_matches_chronos_plot_1, per_matches_chronos_plot_10)


# Percentage match including singletons for chronos -----------------------------------------

per_matches_chronos_plot_0_inc = ggplot(chronos_per_match_inc, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_per_match_inc, Sheet=="lambda_0")) +
  xlab("Data percentage") + 
  ylab("Percentage matches") + 
  #ggtitle("C") + 
  labs(subtitle = "Percentage Match including Singletons, chronos(), lambda = 0") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(legend.position = "none") ;per_matches_chronos_plot_0_inc


per_matches_chronos_plot_1_inc = ggplot(chronos_per_match_inc, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_per_match_inc, Sheet=="lambda_1")) +
  xlab("Data percentage") + 
  ylab("Percentage matches") + 
  #ggtitle("E") + 
  labs(subtitle = "Percentage Match including Singletons, chronos(), lambda = 1") +
  theme_classic() +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;per_matches_chronos_plot_1_inc


per_matches_chronos_plot_10_inc = ggplot(chronos_per_match_inc, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_per_match_inc, Sheet=="lambda_10")) +
  xlab("Data percentage") + 
  ylab("Percentage matches") + 
  #ggtitle("G") + 
  labs(subtitle = "Percentage Match including Singletons, chronos(), lambda = 10") +
  theme_classic() +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;per_matches_chronos_plot_10_inc

# Percentage GMYC singletons for chronos -----------------------------------------

singletons_chronos_plot_0 = ggplot(chronos_singletons, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_singletons, Sheet=="lambda_0")) +
  xlab("Data percentage") + 
  ylab("Percentage GMYC singletons") + 
  #ggtitle("C") + 
  labs(subtitle = "B) chronos(), λ = 0") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  expand_limits(x = 0, y = c(0,100)) +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(legend.position = "none") ;singletons_chronos_plot_0


singletons_chronos_plot_1 = ggplot(chronos_singletons, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_singletons, Sheet=="lambda_1")) +
  xlab("Data percentage") + 
  ylab("Percentage GMYC singletons") + 
  #ggtitle("C") + 
  labs(subtitle = "C) chronos(), λ = 1") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  expand_limits(x = 0, y = c(0,100)) +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(legend.position = "none") ;singletons_chronos_plot_1


singletons_chronos_plot_10 = ggplot(chronos_singletons, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_singletons, Sheet=="lambda_10")) +
  xlab("Data percentage") + 
  ylab("Percentage GMYC singletons") + 
  #ggtitle("C") + 
  labs(subtitle = "D) chronos(), λ = 10") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  expand_limits(x = 0, y = c(0,100)) +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  theme(legend.position = "none") ;singletons_chronos_plot_10

singleton_boxplots = gridExtra::grid.arrange(singletons_pathd8_boxplot, singletons_chronos_plot_0, singletons_chronos_plot_1, singletons_chronos_plot_10)
ggsave(plot = singleton_boxplots, width = 20, height = 18, dpi = 350, filename = "singleton_boxplots.png", units = "cm")

# CREATE A COMBINED GRAPH FOR ALL

gridExtra::grid.arrange(per_matches_chronos_plot_0, per_matches_chronos_plot_1, per_matches_chronos_plot_10,
                        per_matches_chronos_plot_0_inc, per_matches_chronos_plot_1_inc, per_matches_chronos_plot_10_inc, 
                        per_match_pathd8_boxplot, per_match_inc_pathd8_boxplot, ncol = 2)

# Chronos number of entities  -----------------------------------------

entities_chronos_0 = ggplot(chronos_entities, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_entities,Sheet=="lambda_0")) +
  xlab("Data percentage") + 
  ylab("Entities") + 
 # ggtitle("D") + 
  labs(subtitle = "Number of entities, chronos(), lambda = 0") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 100, by = 5)) +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;entities_chronos_0 


entities_chronos_1 = ggplot(chronos_entities, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_entities,Sheet=="lambda_1")) +
  xlab("Data percentage") + 
  ylab("Entities") + 
  #ggtitle("F") + 
  labs(subtitle = "Number of entities, chronos(), lambda = 1") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 85, by = 5)) +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;entities_chronos_1

entities_chronos_10 = ggplot(chronos_entities, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_entities,Sheet=="lambda_10")) +
  xlab("Data percentage") + 
  ylab("Entities") + 
 #ggtitle("H") + 
  labs(subtitle = "Number of entities, chronos(), lambda = 10") +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 100, by = 5)) +
  theme(legend.position = "none") ;entities_chronos_10


# Chronos number of clusters  -----------------------------------------

clusters_chronos_0 = ggplot(chronos_clusters, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_clusters,Sheet=="lambda_0")) +
  xlab("Data percentage") + 
  ylab("Clusters") + 
  # ggtitle("D") + 
  labs(subtitle = "Number of clusters, chronos(), lambda = 0") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;clusters_chronos_0 

clusters_chronos_1 = ggplot(chronos_clusters, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_clusters,Sheet=="lambda_1")) +
  xlab("Data percentage") + 
  ylab("Clusters") + 
  # ggtitle("D") + 
  labs(subtitle = "Number of clusters, chronos(), lambda = 1") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 40, by = 5)) +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;clusters_chronos_1

clusters_chronos_10 = ggplot(chronos_clusters, aes(x = variable, y = value)) + 
  geom_boxplot(data=subset(chronos_clusters,Sheet=="lambda_10")) +
  xlab("Data percentage") + 
  ylab("Clusters") + 
  # ggtitle("D") + 
  labs(subtitle = "Number of clusters, chronos(), lambda = 10") +
  stat_summary(fun = mean, geom="point", shape=18, size=3, color="red", fill="red") +
  scale_y_continuous(breaks = seq(0, 40, by = 5)) +
  theme_classic() +
  theme(plot.title = element_text(size=10, face = "bold")) +
  theme(legend.position = "none") ;clusters_chronos_10

gridExtra::grid.arrange(clusters_chronos_0, clusters_chronos_1, clusters_chronos_10)


###########################################
# GET STATISTICAL SUMMARY VALUES
###########################################

###########################################
                # PATHD8
###########################################

#PATHD8 PERCENTAGE MATCHES EXC. SINGLETONS
stats.per_match_PATHD8 = Rmisc::summarySE(per_match_pathd8, measurevar = "value", groupvars = c("variable"))

#PATHD8 PERCENTAGE MATCHES INC. SINGLETONS
stats.per_match_inc_PATHD8 = Rmisc::summarySE(per_match_inc_singles_pathd8, measurevar = "value", groupvars = c("variable"))

#PATHD8 ENTITIES 
stats.entities.PATHD8 = Rmisc::summarySE(entities_pathd8, measurevar = "value", groupvars = c("variable"))

#PATHD8 CLUSTERS
stats.clusters.PATHD8 = Rmisc::summarySE(clusters_pathd8, measurevar = "value", groupvars = c("variable"))

#PATHD8 SINGLETONS
stats.singeltons.PATHD8 = Rmisc::summarySE(singletons_pathd8, measurevar = "value", groupvars = c("variable"))

ggplot(data = stats.singeltons.PATHD8, aes(x = variable, y = value, group = 1)) + 
  geom_line(colour = "grey") + 
  geom_point() +  
  theme_classic() 

###########################################
              # CHRONOS
###########################################

# CHRONOS PERCENTAGE MATCHES EXC. SINGLETONS
stats.chronos_per_match.0 = Rmisc::summarySE(data=subset(chronos_per_match,Sheet=="lambda_0"), measurevar = "value", groupvars = c("variable"))
stats.chronos_per_match.1 = Rmisc::summarySE(data=subset(chronos_per_match,Sheet=="lambda_1"), measurevar = "value", groupvars = c("variable"))
stats.chronos_per_match.10 = Rmisc::summarySE(data=subset(chronos_per_match,Sheet=="lambda_10"), measurevar = "value", groupvars = c("variable"))

stats.chronos_per_match.0[10,]
stats.chronos_per_match.1[10,]
stats.chronos_per_match.10[10,]

# CHRONOS PERCENTAGE MATCHES INC. SINGLETONS
stats.chronos_per_match_inc.0 = Rmisc::summarySE(data=subset(chronos_per_match_inc,Sheet=="lambda_0"), measurevar = "value", groupvars = c("variable"))
stats.chronos_per_match_inc.1 = Rmisc::summarySE(data=subset(chronos_per_match_inc,Sheet=="lambda_1"), measurevar = "value", groupvars = c("variable"))
stats.chronos_per_match_inc.10 = Rmisc::summarySE(data=subset(chronos_per_match_inc,Sheet=="lambda_10"), measurevar = "value", groupvars = c("variable"))

stats.chronos_per_match_inc.0[10,]
stats.chronos_per_match_inc.1[10,]
stats.chronos_per_match_inc.10[10,]

#write.csv(stats.chronos_per_match.1, file = "stats.chronos_per_match.1.csv")

# CHRONOS ENTITIES
stats.chronos_entities.0 = Rmisc::summarySE(data=subset(chronos_entities,Sheet=="lambda_0"), measurevar = "value", groupvars = c("variable"))
stats.chronos_entities.1 = Rmisc::summarySE(data=subset(chronos_entities,Sheet=="lambda_1"), measurevar = "value", groupvars = c("variable"))
stats.chronos_entities.10 = Rmisc::summarySE(data=subset(chronos_entities,Sheet=="lambda_10"), measurevar = "value", groupvars = c("variable"))

stats.chronos_entities.0[10,]
stats.chronos_entities.1[10,]
stats.chronos_entities.10[10,]

# CHRONOS CLUSTERS
stats.chronos_clusters.0 = Rmisc::summarySE(data=subset(chronos_clusters,Sheet=="lambda_0"), measurevar = "value", groupvars = c("variable"))
stats.chronos_clusters.1 = Rmisc::summarySE(data=subset(chronos_clusters,Sheet=="lambda_1"), measurevar = "value", groupvars = c("variable"))
stats.chronos_clusters.10 = Rmisc::summarySE(data=subset(chronos_clusters,Sheet=="lambda_10"), measurevar = "value", groupvars = c("variable"))

stats.chronos_clusters.0[10,]
stats.chronos_clusters.1[10,]
stats.chronos_clusters.10[10,]

# CHRONOS GMYC SINGLETONS
stats.chronos_singletons.0 = Rmisc::summarySE(data=subset(chronos_singletons,Sheet=="lambda_0"), measurevar = "value", groupvars = c("variable"))
stats.chronos_singletons.1 = Rmisc::summarySE(data=subset(chronos_singletons,Sheet=="lambda_1"), measurevar = "value", groupvars = c("variable"))
stats.chronos_singletons.10 = Rmisc::summarySE(data=subset(chronos_singletons,Sheet=="lambda_10"), measurevar = "value", groupvars = c("variable"))
