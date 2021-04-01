#load packages
library(readr)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(ggsignif)
library(ggiraphExtra)
library(ggpubr)
library(psych)
library(forcats)

#load data
PG1_16S_D40 <- read_csv("PG1_16s_D40.csv", col_types = cols(X1 = col_skip()))
PG1_16S_D40_short = filter(PG1_16S_D40, family != "Unknown") %>% dplyr::select(-`#OTU ID`) %>%spread(family, relative_abundance)
PG1_16S_D40_short = PG1_16S_D40_short[,-(which(colSums(PG1_16S_D40_short[,3:62]) <0.001) + 2)]


PG1_16S_D60 <- read_csv("PG1_16s_D60.csv", col_types = cols(X1 = col_skip()))
PG1_16S_D60_short = filter(PG1_16S_D60, family != "Unknown") %>% dplyr::select(-`#OTU ID`) %>%spread(family, relative_abundance)
PG1_16S_D60_short = PG1_16S_D60_short[,-(which(colSums(PG1_16S_D60_short[,3:62]) <0.001) + 2)]

PG1_16S_D40_D60 = rbind(PG1_16S_D40, PG1_16S_D60)
PG1_16S_D40_D60_short = rbind(PG1_16S_D40_short, PG1_16S_D60_short)

family_short <- read_csv("family_short.csv", col_types = cols(X1 = col_skip()))

species_short <- read_csv("species_short.csv", col_types = cols(X1 = col_skip()))

final_cyt_LOD_short <- read_csv("final_cyt_LOD_short.csv")

final_cyt <- gather(final_cyt_LOD_short, Assay, concentration, 4:23) %>%
  mutate(Treatment2 = factor(ifelse(Treatment == "EVC001" & Day == 6, " EVC001 Baseline", ifelse(Treatment == "Control", " Control", Treatment))))

final_cyt_16s_short = merge(final_cyt_LOD_short, PG1_16S_D40_D60_short %>% rename(Subject = Sample), by = c("Subject", "Day"), all = TRUE)

#######################################################################################
#plot data

##heatmap
cyts.to.use = grepl("IFN-B|12p70|13|17A|21|23|27|31|33|4|MIP", colnames(final_cyt_16s_short))
cor.mat_D60_adj = corr.test(filter(final_cyt_16s_short, Day == "60")[,cyts.to.use], filter(final_cyt_16s_short, Day == "60")[,c(24:48)], method = "spearman", adjust = "fdr")

r.mat = data.frame(cor.mat_D60_adj$r) %>%
  rownames_to_column(var = "cytokine") %>%
  gather(family, spearman_correlation,2:ncol(.)) %>%
  mutate(Day = 60) 

p.mat = data.frame(cor.mat_D60_adj$p) %>%
  rownames_to_column(var = "cytokine") %>%
  gather(family, p_value,2:ncol(.)) %>%
  mutate(Day = 60, symbol = ifelse(p_value <= 0.0001, "****", 
                                   ifelse(p_value > 0.0001 & p_value <= 0.001, "***",
                                          ifelse(p_value > 0.001 & p_value <= 0.01, "**",
                                                 ifelse(p_value > 0.01 & p_value <= 0.05, "*",NA)))))

heatmap_correlation_stats <- merge(r.mat, p.mat, by = c("cytokine", "family", "Day"))
#write.csv(heatmap_correlation_stats, "heatmap_correlation_stats.csv", row.names = FALSE)

#heatmap plot - vertical
ggplot(data = r.mat, aes(x = cytokine, y = reorder(family, desc(family)))) +
  ##heatmap tiles
  geom_tile(aes(fill = spearman_correlation)) +
  #add p-value labels
  geom_text(data = p.mat, aes(label = symbol)) +
  scale_fill_gradient2(low="#24A79C", mid = "white", high="deeppink", name = "Spearman\nCorrelation", limits = c(-1,1), na.value="grey") +
  theme_classic() +
  #change labels
  ylab("Bacterial Family") +
  xlab("Cytokine") +
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5))

#heatmap plot - horizontal
ggplot(data = r.mat, aes(x = cytokine, y = reorder(family, desc(family)))) +
  ##heatmap tiles
  geom_tile(aes(fill = spearman_correlation)) +
  #add p-value labels
  geom_text(data = p.mat, aes(label = symbol)) +
  scale_fill_gradient2(low="#24A79C", mid = "white", high="deeppink", name = "Spearman\nCorrelation", limits = c(-1,1), na.value="grey") +
  #make horizontal
  coord_flip() +
  theme_classic() +
  #change labels
  ylab("Bacterial Family") +
  xlab("Cytokine") +
  theme(text = element_text(size = 18), axis.text.x = element_text(angle = 90, vjust = 0.5))


##radar plot
###calculate median data
PG1_cyt_radar = final_cyt %>%
  group_by(Day, Treatment, Assay) %>%
  summarize(med = median(concentration, na.rm = TRUE)) %>%
  spread(Assay, med)

PG1_cyt_radar$DayLabel   =  paste("Day ", PG1_cyt_radar$Day)
PG1_cyt_radar$DayLabel = as.factor(PG1_cyt_radar$DayLabel)
levels(PG1_cyt_radar$DayLabel) <- c("Day 40","Day 6\n(Baseline)","Day 60")
PG1_cyt_radar$DayLabel = factor(PG1_cyt_radar$DayLabel, levels(PG1_cyt_radar$DayLabel)[c(2,1,3)])
levels(PG1_cyt_radar$DayLabel)
PG1_cyt_radar$Treatment = as.factor(PG1_cyt_radar$Treatment)
levels(PG1_cyt_radar$Treatment) <- c("Control","EVC001")
PG1_cyt_radar = PG1_cyt_radar[,c(1,2,23,17,18,7,10,19,6,4,22,20,9,21,8,15,13,11,16,14,3,5)]
PG1_cyt_radar$Treatment2 = ifelse(PG1_cyt_radar$Treatment == "EVC001" & PG1_cyt_radar$Day == "6", "EVC001 Baseline", paste(PG1_cyt_radar$Treatment))
PG1_cyt_radar$Treatment2 = as.factor(PG1_cyt_radar$Treatment2)
PG1_cyt_radar$Treatment2 = factor(PG1_cyt_radar$Treatment2, levels(PG1_cyt_radar$Treatment2)[c(1,3,2)])

radar <- filter(PG1_cyt_radar, Day != "40")[,-c(1,5,7,8,10,11:13,22)]

ggRadar(group.point.size = 4, data = radar, aes(color = Treatment2, facet = DayLabel), rescale = TRUE, size = 1) +
  theme_minimal() +
  scale_fill_manual(values = c("#686D70", "grey10", "#24A79C")) +
  scale_color_manual(values =  c("#686D70", "grey10", "#24A79C")) +
  scale_x_discrete(labels = c('IL-4', 'IL-13', 'IL-12p7', expression('MIP-3'*alpha), 'IL-17A', 'IL-31', "IL-23", "IL-21", "IL-33", "IL-27", expression('IFN-'*beta))) +
  labs(color = "Treatment", fill = "Treatment") +
  geom_text(aes(y = 1.1,label = ""), show.legend = FALSE) +
  theme(text = element_text(size=12), legend.position = "bottom") 

cyt_med_summary = 
  merge(
      filter(final_cyt, Day != 40 & grepl("IFN-B|12p70|13|17A|21|23|27|31|33|4|MIP", Assay)) %>%
      group_by(Day, Treatment, Assay) %>%
      dplyr::summarize(median = round(median(concentration, na.rm = TRUE), 2), sd = round(sd(concentration, na.rm = TRUE), 2)) %>%
      mutate(title = paste("Day", Day, Treatment), value = paste(median, " (", sd, ")", sep = "")) %>%
      ungroup() %>%
      dplyr::select(Assay, title, value) %>%
      spread(title, value) %>%
      rename(Cytokine = Assay)
    ,
      compare_means(concentration ~ Treatment, data = filter(final_cyt, Day != 40 & grepl("IFN-B|12p70|21|23|27|31|33|MIP|13|17A|4", Assay)), group.by = c("Day", "Assay"), p.adjust.method = "fdr") %>%
      mutate(title = paste("Day", Day, "P-value")) %>%
      dplyr::select(Assay, title, p.adj) %>%
      spread(title, p.adj) %>%
      rename(Cytokine = Assay),
    by = "Cytokine", all.x = TRUE)

#write.csv(cyt_med_summary, "cyt_med_summary.csv", row.names = FALSE)


