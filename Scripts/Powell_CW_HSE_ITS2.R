###Author: Maya Powell
###Last updated: Jan 8th 2025
#This script describes ITS2 processing of Symbiodiniaceae samples
#These samples are from a heat stress experiment completed in July 2023 in Curaçao

#load libraries
#install.packages("BiocManager") #use BiocManager to install packages if not available through CRAN
#BiocManager::install("phyloseq") #replace with any packages to install as needed
library("phyloseq"); packageVersion("phyloseq") #‘1.54.0’
library('ggplot2'); packageVersion("ggplot2") #‘4.0.1’
library('Rmisc'); packageVersion("Rmisc") #1.5.1’
library("cowplot"); packageVersion("cowplot") #‘1.2.0’
library("ggpubr"); packageVersion("ggpubr") #‘0.6.2’
library("vegan"); packageVersion("vegan") #‘2.7.2’
library("tidyverse"); packageVersion("tidyverse") #‘2.0.0’
library("stats"); packageVersion("stats") #‘4.5.2’
library("microbiome"); packageVersion("microbiome") #1.32.0
#remotes::install_github("KarstensLab/microshades")
#library("microshades"); packageVersion("microshades") #only needed for microshades section
#library("remotes")
#remotes::install_github("mikemc/speedyseq")
library("speedyseq"); packageVersion("speedyseq") #‘0.5.3.9021’
#remotes::install_github("david-barnett/microViz")
library("microViz"); packageVersion("microViz") #‘0.12.7’
library("nnet"); packageVersion("nnet") #‘7.3.20’
library("car"); packageVersion("car") #‘3.1.3’
library("emmeans"); packageVersion("emmeans") #‘2.0.0’
library("here"); packageVersion("here") #‘1.0.2’
#here() starts at /Users/mayapowell/Documents/Castillo_Lab/CW_HSE/CW_HSE_Repo

###INITIAL DATA PROCESSING####
#Commented out because only needed once initially - read in all processed data below
#Check for low read samples
counts <- read.csv(here("Data","ITS2",'symportal_profile_counts.csv'),header=TRUE,row.names=1,check.names=FALSE)
plot(rowSums(counts))
counts$sum <- rowSums(counts)
print(counts$sum == "0")
#4 are zeros
counts.no0 <- subset(counts, sum != 0)
print(counts.no0$sum == "0")
counts.no0 <- counts.no0 %>% select(-sum)
plot(rowSums(counts.no0))
#removed A01, G08, G12, H07
#H07 is negative control, A01, G08, G12 are samples, will be removed when making phyloseq

#sample dataframe
samdf<-read.csv(here("Data","ITS2","its2_metadata.csv"))
head(samdf)
#samdf <- samdf %>% select(-X)
rownames(samdf) <- samdf$id

#Make all phyloseq objects
#ran this once and then read them in later
#import taxa info
taxa <- read.csv(here("Data","ITS2","symportal_taxa.csv"),header=TRUE)
rownames(taxa) <- as.factor(taxa$ITS2.type.profile)
mtaxa <- as.matrix(taxa)
# import counts (absolute abundance from its2 type profiles)
mcounts <- as.matrix(counts.no0)
# Construct phyloseq object
ps <- phyloseq(otu_table(mcounts, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(mtaxa))
ps
#22 taxa and 87 samples
#checked and all 4 samples with no counts are correctly removed from this ps object

counts.all <- read.csv(here("Data","ITS2","symportal_all_counts.csv"),header = T, row.names=1)
plot(rowSums(counts.all))
counts.all$sum <- rowSums(counts.all)
print(counts.all$sum == "0")
counts.all.no0 <- subset(counts.all, sum != 0)
print(counts.all.no0$sum == "0")
counts.all.no0 <- counts.all.no0 %>% select(-sum)
plot(rowSums(counts.all.no0))

# import counts (absolute abundance from its2 post med)
mcounts.all <- as.matrix(counts.all.no0)

taxa.all <- read.csv(here("Data","ITS2","symportal_taxa_all.csv"),header=TRUE)
rownames(taxa.all) <- as.factor(taxa.all$ITS2.Type)
mtaxa.all <- as.matrix(taxa.all)
# Construct phyloseq object
ps.all <- phyloseq(otu_table(mcounts.all, taxa_are_rows=FALSE),
                   sample_data(samdf),
                   tax_table(mtaxa.all))
ps.all
#125 taxa and 91 samples
#THIS KEEPS NEGATIVE CONTROL BUT ROWSUMS SHOWS VERY LITTLE IN THOSE SAME 4 SAMPLES - NEED TO MANUALLY REMOVE LATER XOXO

##additional cleanup to add nice labels and remove duplicate samples

# #add variable to make specific site labels
# ps.all <- ps.all %>% ps_mutate(site_nice =
#                       case_when(site == "SMB" ~ "Santa Martha Bay",
#                                 site == "SMR" ~ "Santa Martha Reef",
#                                 site == "SWB" ~ "Spaanse Water Bay",
#                                 site == "DB" ~ "Spaanse Water Reef"))

##save data
saveRDS(ps, here("Data","ITS2","ps.its2.RDS"))

saveRDS(ps.all,here("Data","ITS2","ps.all.its2.RDS"))

#####READ IN DATA####
ps <- readRDS(here("Data","ITS2","ps.its2.RDS")) #includes ALL DATA
ps.all <- readRDS(here("Data","ITS2","ps.all.its2.RDS")) #includes ALL DATA

#create relative dataframes using the no-duplicate (nd) dataframes
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps.all.rel <- transform_sample_counts(ps.all, function(x) x / sum(x))

#subset everything by species

#relative by species
ps.ss.rel <- subset_samples(ps.rel, species=="Siderastrea siderea")
ps.pp.rel <- subset_samples(ps.rel, species=="Branching Porites sp.")

ps.ss.all.rel <- subset_samples(ps.all.rel, species=="Siderastrea siderea")
ps.pp.all.rel <- subset_samples(ps.all.rel, species=="Branching Porites sp.")

######BAR PLOTS#####
#taxa <- as.data.frame(ps.rel@tax_table)
#unique(taxa$Majority_ITS2_sequence)


#c("A4","A4ch","B19","B19aq","B19ao","B19ap","B5","B2","C46","C3","C42ef/C42eg","C46/C3","C46/C1","C47a","C46/C42.2","C42eg/C42ef","C1","C45l","C45a","C3af","C42.2","D1")
# maj_its2_colors = c("A4"= "#ffaabb","A4ch" = "#ffaaff",
#                     "B19"= "#aaaa00","B19aq" = "#aaad77","B19ao" = "#baaa88","B19ap" = "#caaa99",
#                     "B5"= "#44bb99","B2"="#225522",
#                     "C46"="#222255","C3"="#77aadd","C42ef/C42eg"="#225544",
#                     "C46/C3"= "#222588","C46/C1"= "#225599",
#                     "C47a"="#99ddff","C46/C42.2"= "#228899","C42eg/C42ef"="#224456",
#                     "C1"="#eedd88","C45l"="#ee8860","C45a"="#ee8880",
#                     "C3af"="#77aabb","C42.2"="#225555","D1"="#994455")
# 
# ps.all.its2 <- aggregate_rare(ps.all.rel,
#                               level = "ITS2_type",
#                               detection = 0.1,   # 10% RA threshold
#                               prevalence = 0,      # keep all genera, just lump low-RA
#                               aggregate_name = "Other")
# #tax <- as.data.frame(ps.all.its2@tax_table)
# #unique(tax$ITS2_type)
# #c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other")
# its2_cols = c("A_A4" = "#ffaabb",     
#               "B_B19"="#aaaa00", "B_B5" = "#44bb99", 
#               "C_C1"="#eedd88","C_C1c"="lightyellow2","C_C3"="#77aadd",  
#               "C_C42.2"="#225555","C_C42ef"="darkgreen","C_C42eg"="steelblue4",
#               "C_C45a"="#ee8860","C_C45l"="coral3","C_C46"="#222255",
#               "C_C47a"="#99ddff","C_C47c"="lightblue3",  
#               "D_D1"="#994455","D_D4"="darkred","Other"="darkgray")
# 
# ps.ss.all.its2 <- subset_samples(ps.all.its2, host_species=="siderea")
# ps.sr.all.its2 <- subset_samples(ps.all.its2, host_species=="radians")
# ps.pp.all.its2 <- subset_samples(ps.all.its2, host_species=="porites")

##QUICK GLANCES AT DATA##
#all samples ITS2 type profile raw abundance

all.raw <- plot_bar(ps, x = "sample_id", fill="Majority.ITS2.sequence") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Raw Abundance")+
  facet_wrap(treatment~group, scales = "free") +
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))
all.raw
ggsave(all.raw,file=here("Output","ITS2","raw.abundance.group.treatment.species.pdf"),h=30,w=30)

all.relative <- plot_bar(ps.rel, x = "sample_id", fill="Majority.ITS2.sequence") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(treatment~group, scales = "free") +
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))
all.relative
ggsave(all.relative,file=here("Output","ITS2","relative.abundance.group.treatment.species.pdf"),h=30,w=30)


#ITS2 profile by species
gg.bar <- plot_bar(ps.rel,"sample_id",fill="ITS2.type.profile")+
  geom_bar(stat="identity")+
  theme_classic()+
  facet_grid(~species, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  xlab("Species")+
  ylab("Relative Abundance")+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())
gg.bar

#siderastrea siderea - treatment vs control
ss.treat <- plot_bar(ps.ss.rel, x = "sample_id", fill="Majority.ITS2.sequence") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(~treatment, scales = "free") +
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss.treat
ggsave(ss.treat,file=here("Output","ITS2","ss.treatment.maj.its2.barplot.pdf"),h=10,w=20)

#siderastrea siderea - across RTE group
ss.group <- plot_bar(ps.ss.rel, x = "sample_id", fill="Majority.ITS2.sequence") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(~group, scales = "free") +
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss.group
ggsave(ss.group,file=here("Output","ITS2","ss.rte.group.maj.its2.barplot.pdf"),h=15,w=15)

#siderastrea siderea - across RTE group & treatment
ss.group.treat <- plot_bar(ps.ss.rel, x = "sample_id", fill="Majority.ITS2.sequence") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(treatment~group, scales = "free", ncol = 4, nrow = 2) +
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss.group.treat
ggsave(ss.group.treat,file=here("Output","ITS2","ss.rte.group.treatment.maj.its2.barplot.pdf"),h=15,w=20)

#siderastrea siderea - treatment vs control
pp.treat <- plot_bar(ps.pp.rel, x = "sample_id", fill="Majority.ITS2.sequence") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(~treatment, scales = "free") +
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))
pp.treat
ggsave(pp.treat,file=here("Output","ITS2","pp.treatment.maj.its2.barplot.pdf"),h=10,w=20)

#siderastrea siderea - across RTE group
pp.group <- plot_bar(ps.pp.rel, x = "sample_id", fill="Majority.ITS2.sequence") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(~group, scales = "free") +
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))
pp.group
ggsave(pp.group,file=here("Output","ITS2","pp.rte.group.maj.its2.barplot.pdf"),h=15,w=15)

#siderastrea siderea - across RTE group & treatment
pp.group.treat <- plot_bar(ps.pp.rel, x = "sample_id", fill="Majority.ITS2.sequence") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(treatment~group, scales = "free", ncol = 4, nrow = 2) +
  #scale_fill_manual(name = "ITS2 Type", values = its2_cols, labels = c("A4","B19","B5","C1","C1c","C3","C42.2","C42ef","C42eg","C45a","C45l","C46","C47a","C47c","D1","D4","Other"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))
pp.group.treat
ggsave(pp.group.treat,file=here("Output","ITS2","pp.rte.group.treatment.maj.its2.barplot.pdf"),h=15,w=20)

#put plots together - just time for this section
sym_all_time <- ggarrange(ss.site.time,
                          ggarrange(sr.site.time, pp.site.time, ncol = 2, 
                                    labels = c("B", "C"), font.label = list(size = 50)), 
                          #heights = c(1,2), widths = c(2,1), 
                          nrow = 2,
                          #legend = "right", common.legend = TRUE, 
                          labels = c("A"),font.label = list(size = 50))
sym_all_time
ggsave(sym_all_time, file=here("ITS2","its2.type.plot.time.all.pdf"),width=25,height=30)


#### Dominant ITS2 Majority Types and DIVs ####

#commented out initial section needed to create majority types dataset
# 
# #generate sequence table
# #seq.all <- data.frame(ps.all.rel@otu_table)
# seq <- data.frame(ps.rel@otu_table)
# samdf <- data.frame(ps.rel@sam_data)
# 
# # change column names to majority its2 sequence
# sym_taxa = read.csv(here("ITS2","symportal_taxa.csv"), header = TRUE)
# 
# sym_taxa$Majority_ITS2_sequence #you need to do it this way because you need to add the 1s etc to not have duplicates
# colnames(seq) =  c("A4.1","A4.2","A4.3","A4ch","A4.4",
#                    "A4.5","B19.1","B19.2","B19aq","B19.3",
#                    "B19ao","B19ap","B5","B19.4","B19.5",
#                    "B2","C46.1","C3.1","C42ef.C42eg","C3.2",
#                    "C46.C3","C46.C1","C47a.1","C46.C42.2","C42eg.C42ef",
#                    "C47a.2","C1.1","C1.2","C3.3","C1.3",
#                    "C1.4","C47a.3","C45l","C1.5","C3.4",
#                    "C45a","C1.6","C3.5","C1.7","C1.8",
#                    "C3af","C1.9","C42.2.1","C46.2","C42.2.2",
#                    "D1.1","D1.2","D1.3","D1.4","D1.5","D1.6","D1.7")
# 
# # make new data frame and sum columns with the same majority its2 sequence
# seqtab.rel.sums <- seq %>%
#   mutate(A4_sum = rowSums(select(., starts_with("A4")))) %>%
#   #mutate(A4ch_sum = A4ch) %>%
#   mutate(B19_sum = rowSums(select(., starts_with("B19")))) %>%
#   #mutate(B19aq_sum = B19aq) %>%
#   #mutate(B19ao_sum = B19ao) %>%
#   #mutate(B19ap_sum = B19ap) %>%
#   mutate(B5_sum = B5) %>%
#   mutate(B2_sum = B2) %>%
#   mutate(C46_sum = rowSums(select(., starts_with("C46.")))) %>%
#   mutate(C3_sum = rowSums(select(., starts_with("C3")))) %>%
#   #mutate(C42ef.C42eg_sum = C42ef.C42eg) %>%
#   #mutate(C46.C3_sum = C46.C3) %>%
#   #mutate(C46.C1_sum = C46.C1) %>%
#   mutate(C47a_sum = rowSums(select(., starts_with("C47a.")))) %>%
#   #mutate(C46.C42.2_sum =C46.C42.2) %>%
#   #mutate(C42eg.C42ef_sum = C42eg.C42ef) %>%
#   mutate(C1_sum = rowSums(select(., starts_with("C1.")))) %>%
#   mutate(C45_sum = rowSums(select(., starts_with("C45")))) %>%
#   #mutate(C45l_sum = C45l) %>%
#   #mutate(C45a_sum = C45a) %>%
#   #mutate(C3af_sum = C3af) %>%
#   mutate(C42_sum = rowSums(select(., starts_with("C42")))) %>%
#   #mutate(C42.2_sum = rowSums(select(., starts_with("C42.2.")))) %>%
#   mutate(D1_sum = rowSums(select(., starts_with("D1.")))) %>%
#   rownames_to_column(var = "id") %>%
#   select(id, contains("_sum"))
# 
# its2.sums.rel <- left_join(seqtab.rel.sums, samdf, by = "id")
# write.csv(its2.sums.rel, file = here("ITS2","its2.sums.rel.csv"), row.names = FALSE)
# 
# #use seqtab rel sum to make graphs
# #"A4","B19","B2","C46","C3","C47a","C1","C45","C42","D1"
# taxa.maj <- data.frame(colnames(seqtab.rel.sums))
# taxa.maj <- data.frame(taxa.maj[-c(1), ])
# colnames(taxa.maj)[1] <- 'maj_its2'
# taxa.maj$genus <- c("A","B","B","B","C","C","C","C","C","C","D")
# rownames(taxa.maj) <- taxa.maj$maj_its2
# 
# taxa.maj$maj_its2 = as.factor(taxa.maj$maj_its2)
# taxa.maj$genus = as.factor(taxa.maj$genus)
# 
# #seqtab rel sums
# rownames(seqtab.rel.sums) <- seqtab.rel.sums$id
# seqtab.rel.sums.maj <- seqtab.rel.sums %>% select(-id)
# #col names of seqtab.rel.sums.maj match rownames of taxa.maj
# #rownames of seqtab.rel.sums.maj match rownames of samdf
# 
# #make ps object of majority type sums
# #MADE WITH ND DATASET - no duplicates
# ps.its2.maj <- phyloseq(sample_data(samdf),
#                        otu_table(seqtab.rel.sums.maj,taxa_are_rows=FALSE),
#                        tax_table(as.matrix(taxa.maj)))
# ps.its2.maj
# # phyloseq-class experiment-level object
# # otu_table()   OTU Table:          [ 11 taxa and 143 samples ]:
# # sample_data() Sample Data:        [ 143 samples by 33 sample variables ]:
# # tax_table()   Taxonomy Table:     [ 11 taxa by 2 taxonomic ranks ]:
# # taxa are columns
# 
# ps.its2.maj <- ps.its2.maj %>% ps_mutate(tp_suffix = case_when(
#   m_y == "March_2020"    ~ "a",
#   m_y == "November_2020" ~ "b",
#   m_y == "November_2021" ~ "c")) %>%
#   ps_mutate(number_time = paste0(number,tp_suffix))
# saveRDS(ps.its2.maj, here("ITS2", "ps.its2.maj.nd.RDS"))

##READ IN MAJORITY DATA
#this is saved as a relative abundance phyloseq
ps.its2.maj <- readRDS(here("ITS2", "ps.its2.maj.nd.RDS"))

#subset each species
ps.ss.maj <- subset_samples(ps.its2.maj, host_species=="siderea")
ps.sr.maj <- subset_samples(ps.its2.maj, host_species=="radians")
ps.pp.maj <- subset_samples(ps.its2.maj, host_species=="porites")

#quick glance at data
plot_bar(ps.its2.maj, x="id",fill="maj_its2")+
  theme_classic()

#assign colors
its2_colors = c("A4_sum" = "#ffaabb", "C47a_sum" = "#99ddff","C42_sum" = "#225555",
                "B19_sum" =  "#aaaa00","B5_sum" = "#44bb99", "B2_sum" = "#225522",
                "C46_sum" = "#222255","C3_sum" = "#77aadd",
                  "C1_sum" = "#eedd88",
                     "C45_sum" = "#ee8860", "D1_sum" = "#994455")
#"A4","B19","B2","B5","C1","C3","C42","C45","C46","C47a","D1"

###SIDERASTREA SIDEREA###
ss.site <- plot_bar(ps.ss.maj, x="number_time", fill="maj_its2") +
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone.1, scales = "free")+
  scale_fill_manual(name = "Majority ITS2", values = its2_colors, labels = c("A4","B19","B2","B5","C1","C3","C42","C45","C46","C47a","D1")) +
  theme(axis.text.x=element_text(angle = 90, size = 15),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss.site

#full plot - to put in supplementary, with all data across sites and time
#numbers are individual samples here
ss.site.time <- plot_bar(ps.ss.maj, x = "number_time", fill="maj_its2") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone.1, scales = "free") +
  scale_fill_manual(name = "Majority ITS2", values = its2_colors, labels = c("A4","B19","B5","B2","C46","C3","C47a","C1","C45","C42","D1")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  ggtitle(substitute(paste(italic("Siderastrea siderea"))))
ss.site.time

###Siderastrea radians
sr.site <- plot_bar(ps.sr.maj, x="number_time", fill="maj_its2") +
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone.1, scales = "free")+
  scale_fill_manual(name = "Majority ITS2", values = its2_colors, labels = c("A4","B19","B5","B2","C46","C3","C47a","C1","C45","C42","D1")) +
  theme(axis.text.x=element_text(angle = 90, size = 15),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
sr.site

#full plot - to put in supplementary, with all data acrosr sites and time
#numbers are individual samples here
sr.site.time <- plot_bar(ps.sr.maj, x = "number_time", fill="maj_its2") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone.1, scales = "free", nrow = 3, ncol = 2) +
  scale_fill_manual(name = "Majority ITS2", values = its2_colors, labels = c("A4","B19","B5","B2","C46","C3","C47a","C1","C45","C42","D1")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  ggtitle(substitute(paste(italic("Siderastrea radians"))))
sr.site.time

#####Branching Porites sp.
pp.site <- plot_bar(ps.pp.maj, x="number_time", fill="maj_its2") +
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~site_zone.1, scales = "free")+
  scale_fill_manual(name = "Majority ITS2", values = its2_colors, labels = c("A4","B19","B5","B2","C46","C3","C47a","C1","C45","C42","D1")) +
  theme(axis.text.x=element_text(angle = 90, size = 15),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))
pp.site

#full plot - to put in supplementary, with all data acropp sites and time
#numbers are individual samples here
pp.site.time <- plot_bar(ps.pp.maj, x = "number_time", fill="maj_its2") + 
  theme_classic(base_size = 30)+
  xlab("Coral Genotype")+
  ylab("Relative Abundance")+
  facet_wrap(m_y~site_zone.1, scales = "free", nrow = 2, ncol = 2) +
  scale_fill_manual(name = "Majority ITS2", values = its2_colors, labels = c("A4","B19","B5","B2","C46","C3","C47a","C1","C45","C42","D1")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "bottom",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30))+
  xlab("Sample ID") +
  labs(title = expression(paste("Branching ", italic("Porites "), "sp.")))
pp.site.time

###put plots together
sym_pp_sr_plot <- ggarrange(sr.site,pp.site, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 40))
sym_all_plot <- ggarrange(ss.site,sym_pp_sr_plot, nrow=1, ncol = 2, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 50))
sym_all_plot
ggsave(sym_all_plot, file=here("ITS2","its2.plot.all.pdf"),width=25,height=12)

#put plots together
sym_all_time <- ggarrange(ss.site.time,
                          ggarrange(sr.site.time, pp.site.time, ncol = 2, 
                                    labels = c("B", "C"), font.label = list(size = 50)), 
                          #heights = c(1,2), widths = c(2,1), 
                          nrow = 2,
                          #legend = "right", common.legend = TRUE, 
                          labels = c("A"),font.label = list(size = 50))
sym_all_time
ggsave(sym_all_time, file=here("ITS2","its2.plot.time.all.pdf"),width=25,height=30)


#####STATS###

#commenting out initial data processing - only needed once, use its2.dom.rel dataframe below
# its2.sums.rel <- read.csv(here("ITS2", "its2.sums.rel.csv"))
# 
# #convert to factors and numeric as needed
# its2.sums.rel = its2.sums.rel %>%
#   mutate_at(c(2:12), as.numeric)
# 
# its2.sums.rel = its2.sums.rel %>%
#   mutate_at(c(13:35), as.factor)
# 
# # add in dominant and minor distinctions to use in other plots
# its2.dom.rel <- its2.sums.rel %>%
#   mutate(dominant_type = case_when(A4_sum >= 0.7 ~ "A4",
#                                    B19_sum >= 0.7 ~ "B19",
#                                    B5_sum >= 0.7 ~ "B5",
#                                    B2_sum >= 0.7 ~ "B2",
#                                    C46_sum >= 0.7 ~ "C46",
#                                    C3_sum >= 0.7 ~ "C3",
#                                    C42_sum >= 0.7 ~ "C42",
#                                    C47a_sum >= 0.7 ~ "C47a",
#                                    C1_sum >= 0.7 ~ "C1",
#                                    C45_sum >= 0.7 ~ "C45",
#                                    D1_sum >= 0.7 ~ "D1")) %>%
#   mutate(minor_type = case_when(A4_sum < 0.7 & A4_sum > 0.0 ~ "A4",
#                                 B19_sum < 0.7 & B19_sum > 0.0 ~ "B19",
#                                 B5_sum < 0.7 & B5_sum > 0.0 ~ "B5",
#                                 B2_sum < 0.7 & B2_sum > 0.0 ~ "B2",
#                                 C46_sum < 0.7 & C46_sum > 0.0 ~ "C46",
#                                 C3_sum < 0.7 & C3_sum > 0.0 ~ "C3",
#                                 C42_sum < 0.7 & C42_sum > 0.0 ~ "C42",
#                                 C47a_sum < 0.7 & C47a_sum > 0.0 ~ "C47a",
#                                 C1_sum < 0.7 & C1_sum > 0.0 ~ "C1",
#                                 C45_sum < 0.7 & C45_sum > 0.0 ~ "C45",
#                                 D1_sum < 0.7 & D1_sum > 0.0 ~ "D1"))
# 
# #filter out Porites astreoides bc not using in this analysis
# its2.dom.rel <- filter(its2.dom.rel, host_species!="astreoides")
# 
# write.csv(its2.dom.rel, file = here("ITS2","ITS2.dominanttype.CW_2020.csv"), row.names = FALSE)


####STATS#####
its2.dom.rel <- read.csv(here("ITS2","ITS2.dominanttype.CW_2020.csv"))

ss.its2.rel.dom <- filter(its2.dom.rel, host_species=="siderea")
sr.its2.rel.dom <- filter(its2.dom.rel, host_species=="radians")
pp.its2.rel.dom <- filter(its2.dom.rel, host_species=="porites")

# summarize proportion of individuals with more than 70% of the sym types
its2.dom.spp <- its2.dom.rel %>%
  group_by(host_species) %>%
  summarize(n = n(),
            n_maj_A4 = sum(A4_sum >= 0.7),
            #n_maj_B19 = sum(B19_sum >= 0.7),
            #n_maj_B5 = sum(B5_sum >= 0.7),
            #n_maj_B2 = sum(B2_sum >= 0.7),
            n_maj_C46 = sum(C46_sum >= 0.7),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C42 = sum(C42_sum >= 0.7),
            n_maj_C47a = sum(C47a_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            #n_maj_C45 = sum(C45_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_A4 = n_maj_A4/n,
            #p_maj_B19 = n_maj_B19/n,
            #p_maj_B5 = n_maj_B5/n,
            #p_maj_B2 = n_maj_B2/n,
            p_maj_C46 = n_maj_C46/n,
            p_maj_C3 = n_maj_C3/n,
            p_maj_C42 = n_maj_C42/n,
            p_maj_C47a = n_maj_C47a/n,
            p_maj_C1 = n_maj_C1/n,
            #p_maj_C45 = n_maj_C45/n,
            p_maj_D1 = n_maj_D1/n)

#based on this - removing ITS2 type summary calculations from the species that do not have any majority of certain sym types
#pp keep = A4, C42, C47a
#sr keep = C46, C1, D1
#ss keep = C3, C1, D1
#all remove = B19, B5, B2, C45

####SSID######

#timepoint
ss.its2.rel.dom %>%
  group_by(m_y) %>%
  summarize(n = n(),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C3 = n_maj_C3/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   m_y               n n_maj_C3 n_maj_C1 n_maj_D1 p_maj_C3 p_maj_C1 p_maj_D1
#   <fct>         <int>    <int>    <int>    <int>    <dbl>    <dbl>    <dbl>
# 1 March 2020       24        6        3       12    0.25    0.125     0.5  
# 2 November 2020    24        6        1       13    0.25    0.0417    0.542
# 3 November 2021    30        7        3       18    0.233   0.1       0.6  

#reef vs bay
ss.its2.rel.dom %>%
  group_by(reef_bay) %>%
  summarize(n = n(),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C3 = n_maj_C3/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   reef_bay     n n_maj_C3 n_maj_C1 n_maj_D1 p_maj_C3 p_maj_C1 p_maj_D1
#   <fct>    <int>    <int>    <int>    <int>    <dbl>    <dbl>    <dbl>
# 1 bay         40        0        3       33      0      0.075    0.825
# 2 reef        38       19        4       10      0.5    0.105    0.263

#site
ss.its2.rel.dom %>%
  group_by(site_zone) %>%
  summarize(n = n(),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C3 = n_maj_C3/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   site_zone     n n_maj_C3 n_maj_C1 n_maj_D1 p_maj_C3 p_maj_C1 p_maj_D1
#   <fct>     <int>    <int>    <int>    <int>    <dbl>    <dbl>    <dbl>
# 1 SM_bay       22        0        2       17    0       0.0909    0.773
# 2 SM_reef      21        8        4        6    0.381   0.190     0.286
# 3 SW_bay       18        0        1       16    0       0.0556    0.889
# 4 SW_reef      17       11        0        4    0.647   0         0.235

########SRAD######

#timepoint
sr.its2.rel.dom %>%
  group_by(m_y) %>%
  summarize(n = n(),
            n_maj_C46 = sum(C46_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C46 = n_maj_C46/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   m_y               n n_maj_C46 n_maj_C1 n_maj_D1 p_maj_C46 p_maj_C1 p_maj_D1
#   <chr>         <int>     <int>    <int>    <int>     <dbl>    <dbl>    <dbl>
# 1 March_2020       12        12        0        0     1        0       0     
# 2 November_2020    14        10        2        1     0.714    0.143   0.0714
# 3 November_2021    10        10        0        0     1        0       0     

#site
sr.its2.rel.dom %>%
  group_by(site_zone) %>%
  summarize(n = n(),
            n_maj_C46 = sum(C46_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C46 = n_maj_C46/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_D1 = n_maj_D1/n)
#   site_zone     n n_maj_C46 n_maj_C1 n_maj_D1 p_maj_C46 p_maj_C1 p_maj_D1
#   <chr>     <int>     <int>    <int>    <int>     <dbl>    <dbl>    <dbl>
# 1 SM_bay       17        15        2        0     0.882    0.118   0     
# 2 SW_bay       19        17        0        1     0.895    0       0.0526

####PPOR######

#timepoint
pp.its2.rel.dom %>%
  group_by(m_y) %>%
  summarize(n = n(),
            n_maj_A4 = sum(A4_sum >= 0.7),
            n_maj_C42 = sum(C42_sum >= 0.7),
            n_maj_C47a = sum(C47a_sum >= 0.7),
            p_maj_A4 = n_maj_A4/n,
            p_maj_C42 = n_maj_C42/n,
            p_maj_C47a = n_maj_C47a/n)
#   m_y               n n_maj_A4 n_maj_C42 n_maj_C47a p_maj_A4 p_maj_C42 p_maj_C47a
#   <chr>         <int>    <int>     <int>      <int>    <dbl>     <dbl>      <dbl>
# 1 November_2020    11        6         2          1    0.545     0.182     0.0909
# 2 November_2021    12        6         2          4    0.5       0.167     0.333 

#site
pp.its2.rel.dom %>%
  group_by(site_zone) %>%
  summarize(n = n(),
            n_maj_A4 = sum(A4_sum >= 0.7),
            n_maj_C42 = sum(C42_sum >= 0.7),
            n_maj_C47a = sum(C47a_sum >= 0.7),
            p_maj_A4 = n_maj_A4/n,
            p_maj_C42 = n_maj_C42/n,
            p_maj_C47a = n_maj_C47a/n)
#   site_zone     n n_maj_A4 n_maj_C42 n_maj_C47a p_maj_A4 p_maj_C42 p_maj_C47a
#   <chr>     <int>    <int>     <int>      <int>    <dbl>     <dbl>      <dbl>
# 1 SM_bay       12       12         0          0        1     0          0    
# 2 SM_reef      11        0         4          5        0     0.364      0.455


######MULTINOMIAL MODELS#####
#####STATS#####

###Trying multinomial models using nnet

#coral species
spp_model <-multinom(dominant_type ~ coral_species, data = its2.dom.rel)
Anova(spp_model, type = 3)
summary(spp_model)
# Response: dominant_type
#               LR Chisq Df Pr(>Chisq)    
# coral_species   226.95 12  < 2.2e-16 ***
coeffs<-tidy(spp_model)
results<-glance(spp_model) 

spp_em <- emmeans(spp_model, ~ coral_species | dominant_type, type = "response")
pairs(spp_em)

#siderastrea siderea
ss_null <- multinom(dominant_type ~ 1, data = ss.its2.rel.dom)
ss_habitat <- multinom(dominant_type ~ reef_bay, data = ss.its2.rel.dom)
ss_hab_loc <- multinom(dominant_type ~ reef_bay*area, data = ss.its2.rel.dom)
ss_site <- multinom(dominant_type ~ site_nice, data = ss.its2.rel.dom)
ss_time <- multinom(dominant_type ~ m_y, data = ss.its2.rel.dom)
ss_add <- multinom(dominant_type ~ site_nice + m_y, data = ss.its2.rel.dom)
ss_int <- multinom(dominant_type ~ reef_bay * area * m_y, data = ss.its2.rel.dom)

anova(ss_null, ss_habitat, ss_hab_loc, ss_time, ss_add, ss_int)
#lowest p = int, but second lowest = habitat
AIC(ss_null, ss_habitat, ss_hab_loc, ss_time, ss_add, ss_int)
#lowest AIC = habitat

Anova(ss_habitat, type = 3)

ss_t_em <- emmeans(ss_habitat, ~ reef_bay | dominant_type)
pairs(ss_t_em)


#siderastrea radians
sr_null <- multinom(dominant_type ~ 1, data = sr.its2.rel.dom)
sr_site <- multinom(dominant_type ~ site_nice, data = sr.its2.rel.dom)
sr_time <- multinom(dominant_type ~ m_y, data = sr.its2.rel.dom)
sr_add <- multinom(dominant_type ~ site_nice + m_y, data = sr.its2.rel.dom)
sr_int <- multinom(dominant_type ~ site_nice * m_y, data = sr.its2.rel.dom)

anova(sr_null, sr_site, sr_time, sr_add, sr_int)
#no significant p-value for any
AIC(sr_null, sr_site, sr_time, sr_add, sr_int)

Anova(sr_time, type = 3)

#branching porites sp.
pp_null <- multinom(dominant_type ~ 1, data = pp.its2.rel.dom)
pp_site <- multinom(dominant_type ~ site_nice, data = pp.its2.rel.dom)
pp_time <- multinom(dominant_type ~ m_y, data = pp.its2.rel.dom)
pp_add <- multinom(dominant_type ~ site_nice + m_y, data = pp.its2.rel.dom)
pp_int <- multinom(dominant_type ~ site_nice * m_y, data = pp.its2.rel.dom)

anova(pp_null, pp_site, pp_time, pp_add, pp_int)
#pp_site = best, smallest p = 5.918183e-07
AIC(pp_null, pp_site, pp_time, pp_add, pp_int)
#lowest AIC = pp_site  df 4 AIC 20.36732

Anova(pp_time, type =3)
# Response: dominant_type
#           LR Chisq Df Pr(>Chisq)    
# site_nice    28.68  2  5.918e-07 ***  

pp_t_em <- emmeans(pp_site, ~ site_nice | dominant_type)
pairs(pp_t_em)

# dominant_type = A4:
#  contrast                             estimate      SE df t.ratio p.value
#  Santa Martha Bay - Santa Martha Reef    1.000 0.00294  4 339.947  <.0001
# 
# dominant_type = C42:
#  contrast                             estimate      SE df t.ratio p.value
#  Santa Martha Bay - Santa Martha Reef   -0.444 0.16600  4  -2.683  0.0551
# 
# dominant_type = C47a:
#  contrast                             estimate      SE df t.ratio p.value
#  Santa Martha Bay - Santa Martha Reef   -0.556 0.16600  4  -3.354  0.0285

####OTHER THINGS I DID AND TRIED THAT MIGHT BE HELPFUL TO YOU!!

# #####K-W STATS ON PROPORTIONS#####
# # Kruskal-Wallis (non-parametric alternative to one-way anova) to see if there is a difference
# # in mean proportion D based on timepoint, reef vs bay, and site
# 
# ##CORAL SPECIES####
# kruskal.test(C47a_sum ~ host_species, data = its2.rel.dom)
# dunnTest(C47a_sum ~ host_species, data = its2.rel.dom, method="bonferroni")
# # data:  C3_sum by host_species
# # Kruskal-Wallis chi-squared = 30.343, df = 2, p-value = 2.577e-07
# # Comparison         Z      P.unadj        P.adj
# # 1 porites - radians  0.000000 1.000000e+00 1.000000e+00
# # 2 porites - siderea -4.367233 1.258306e-05 3.774918e-05
# # 3 radians - siderea -4.543210 5.540404e-06 1.662121e-05
# 
# #data:  D1_sum by host_species
# #Kruskal-Wallis chi-squared = 87.675, df = 2, p-value < 2.2e-16
# #         Comparison           Z      P.unadj        P.adj
# # 1 porites - radians  0.01162176 9.907274e-01 1.000000e+00
# # 2 porites - siderea -7.41652703 1.202314e-13 3.606941e-13
# # 3 radians - siderea -7.72931509 1.081267e-14 3.243802e-14
# 
# # data:  C1_sum by host_species
# # Kruskal-Wallis chi-squared = 8.8846, df = 2, p-value = 0.01177
# # Comparison          Z     P.unadj      P.adj
# # 1 porites - radians -0.7596603 0.447457678 1.00000000
# # 2 porites - siderea -2.7481119 0.005993955 0.01798186
# # 3 radians - siderea -1.9476548 0.051456277 0.15436883
# 
# # data:  C46_sum by host_species
# # Kruskal-Wallis chi-squared = 110, df = 2, p-value < 2.2e-16
# # Comparison         Z      P.unadj        P.adj
# # 1 porites - radians -8.092023 5.868194e-16 1.760458e-15
# # 2 porites - siderea  0.358183 7.202064e-01 1.000000e+00
# # 3 radians - siderea 10.078778 6.857379e-24 2.057214e-23
# 
# # data:  A4_sum by host_species
# # Kruskal-Wallis chi-squared = 53.756, df = 2, p-value = 2.124e-12
# # Comparison          Z      P.unadj        P.adj
# # 1 porites - radians  6.4502792 1.116443e-10 3.349329e-10
# # 2 porites - siderea  6.6946068 2.162523e-11 6.487570e-11
# # 3 radians - siderea -0.7725698 4.397771e-01 1.000000e+00
# 
# # data:  C42_sum by host_species
# # Kruskal-Wallis chi-squared = 33.21, df = 2, p-value = 6.145e-08
# # Comparison         Z      P.unadj        P.adj
# # 1 porites - radians 4.3334275 1.468056e-05 4.404169e-05
# # 2 porites - siderea 5.6528673 1.577932e-08 4.733795e-08
# # 3 radians - siderea 0.6828201 4.947205e-01 1.000000e+00
# 
# # data:  C47a_sum by host_species
# # Kruskal-Wallis chi-squared = 32.669, df = 2, p-value = 8.053e-08
# # Comparison          Z      P.unadj        P.adj
# # 1 porites - radians  4.8390231 1.304789e-06 3.914367e-06
# # 2 porites - siderea  5.3608265 8.284206e-08 2.485262e-07
# # 3 radians - siderea -0.2274366 8.200843e-01 1.000000e+00
# 
# ####SSID####
# library(FSA)
# #sym types are: D1, C1, and C3
# kruskal.test(D1_sum ~ site_zone, data = ss.its2.rel.dom)
# # data:  D1_sum by m_y
# # Kruskal-Wallis chi-squared = 0.44807, df = 2, p-value = 0.7993
# # data:  D1_sum by reef_bay
# # Kruskal-Wallis chi-squared = 21.166, df = 1, p-value = 4.212e-06
# # data:  D1_sum by site_zone
# # Kruskal-Wallis chi-squared = 23.323, df = 3, p-value = 3.458e-05
# dunnTest(D1_sum ~ site_zone, data = ss.its2.rel.dom, method="bonferroni")
# #          Comparison          Z      P.unadj        P.adj
# # 1  SM_bay - SM_reef  2.8629102 4.197695e-03 0.0251861676
# # 2   SM_bay - SW_bay -1.4455974 1.482902e-01 0.8897409920
# # 3  SM_reef - SW_bay -4.1495122 3.331846e-05 0.0001999107
# # 4  SM_bay - SW_reef  2.4418297 1.461304e-02 0.0876782289
# # 5 SM_reef - SW_reef -0.2602169 7.946964e-01 1.0000000000
# # 6  SW_bay - SW_reef  3.6900069 2.242480e-04 0.0013454879
# 
# kruskal.test(C1_sum ~ site_zone, data = ss.its2.rel.dom)
# # data:  C1_sum by m_y
# # Kruskal-Wallis chi-squared = 0.70146, df = 2, p-value = 0.7042
# # data:  C1_sum by reef_bay
# # Kruskal-Wallis chi-squared = 0.92886, df = 1, p-value = 0.3352
# # data:  C1_sum by site_zone
# # Kruskal-Wallis chi-squared = 4.7558, df = 3, p-value = 0.1906
# kruskal.test(C3_sum ~ site_zone, data = ss.its2.rel.dom)
# # data:  C3_sum by m_y
# # Kruskal-Wallis chi-squared = 0.10039, df = 2, p-value = 0.951
# # data:  C3_sum by reef_bay
# # Kruskal-Wallis chi-squared = 35.042, df = 1, p-value = 3.226e-09
# # data:  C3_sum by site_zone
# # Kruskal-Wallis chi-squared = 37.907, df = 3, p-value = 2.958e-08
# dunnTest(C3_sum ~ site_zone, data = ss.its2.rel.dom, method="bonferroni")
# #          Comparison          Z      P.unadj        P.adj
# # 1  SM_bay - SM_reef -3.4273855 6.094232e-04 3.656539e-03
# # 2   SM_bay - SW_bay  0.3825644 7.020427e-01 1.000000e+00
# # 3  SM_reef - SW_bay  3.6338211 2.792546e-04 1.675527e-03
# # 4  SM_bay - SW_reef -4.9037320 9.403274e-07 5.641965e-06
# # 5 SM_reef - SW_reef -1.6486878 9.921162e-02 5.952697e-01
# # 6  SW_bay - SW_reef -5.0417118 4.613857e-07 2.768314e-06
# 
# ####SRAD####
# #sym types are: C46, C1, D1
# kruskal.test(D1_sum ~ m_y, data = sr.its2.rel.dom)
# # data:  D1_sum by m_y
# # Kruskal-Wallis chi-squared = 1.4667, df = 2, p-value = 0.4803
# # data:  D1_sum by site_zone
# # Kruskal-Wallis chi-squared = 0.94737, df = 1, p-value = 0.3304
# kruskal.test(C1_sum ~ site_zone, data = sr.its2.rel.dom)
# # data:  C1_sum by m_y
# # Kruskal-Wallis chi-squared = 3.0148, df = 2, p-value = 0.2215
# # data:  C1_sum by site_zone
# # Kruskal-Wallis chi-squared = 2.1698, df = 1, p-value = 0.1407
# kruskal.test(C46_sum ~ m_y, data = sr.its2.rel.dom)
# # data:  C46_sum by site_zone
# # Kruskal-Wallis chi-squared = 1.052, df = 1, p-value = 0.305
# # data:  C46_sum by m_y
# # Kruskal-Wallis chi-squared = 6.9282, df = 2, p-value = 0.0313
# dunnTest(C46_sum ~ m_y, data = sr.its2.rel.dom, method="bonferroni")
# #                      Comparison          Z    P.unadj      P.adj
# # 1    March 2020 - November 2020  2.4266191 0.01524025 0.04572074
# # 2    March 2020 - November 2021  0.3378417 0.73548248 1.00000000
# # 3 November 2020 - November 2021 -1.9477617 0.05144349 0.15433046
# 
# ####PPOR######
# #sym types are A4, C42, C47a
# kruskal.test(A4_sum ~ site_zone, data = pp.its2.rel.dom)
# # data:  A4_sum by m_y
# # Kruskal-Wallis chi-squared = 0.56796, df = 1, p-value = 0.4511
# # data:  A4_sum by site_zone
# # Kruskal-Wallis chi-squared = 20.711, df = 1, p-value = 5.342e-06
# kruskal.test(C42_sum ~ m_y, data = pp.its2.rel.dom)
# # data:  C42_sum by m_y
# # Kruskal-Wallis chi-squared = 0.28407, df = 1, p-value = 0.594
# # data:  C42_sum by site_zone
# # Kruskal-Wallis chi-squared = 0.015913, df = 1, p-value = 0.8996
# kruskal.test(C47a_sum ~ site_zone, data = pp.its2.rel.dom)
# # data:  C47a_sum by m_y
# # Kruskal-Wallis chi-squared = 1.5364, df = 1, p-value = 0.2152
# # data:  C47a_sum by site_zone
# # Kruskal-Wallis chi-squared = 16.5, df = 1, p-value = 4.865e-05


####MICROSHADES GRAPHS####
# 
# #testing out microshades palette
# # Use microshades function prep_mdf to agglomerate, normalize, and melt the phyloseq object
# ps.rel <- ps.rel %>% tax_mutate(Majority_ITS2_sequence = gsub('/','', as.character(Majority_ITS2_sequence)))
# 
# mdf_prep <- prep_mdf(ps.rel, subgroup_level = "Majority_ITS2_sequence")
# # Create a color object for the specified data
# color_obj_prep_test <- create_color_dfs(mdf_prep, selected_groups = c("D","C","B","A"), group_level = "Clade", subgroup_level = "Majority_ITS2_sequence", cvd = TRUE)
# color_obj_prep_test$mdf$Majority_ITS2_sequence <- as.character(color_obj_prep_test$mdf$Majority_ITS2_sequence)
# # Extract
# mdf_test <- color_obj_prep_test$mdf
# cdf_test <- color_obj_prep_test$cdf
# 
# #subset different groups out
# mdf_ss <- mdf_test %>% filter(host_species == "siderea")
# mdf_sr <- mdf_test %>% filter(host_species == "radians")
# mdf_pp <- mdf_test %>% filter(host_species == "porites")
# 
# #now make graphs
# 
# #ssid
# ss_shade <- plot_microshades(mdf_ss, cdf_test, group_label = "Clade - Majority ITS2") + 
#   theme_classic(base_size = 30)+
#   ylab("Relative Abundance")+
#   facet_wrap(~site_nice, scales = "free")+
#   theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
#   ggtitle(substitute(paste(italic("Siderastrea siderea"))))
# ss_shade
# ggsave(ss_shade,file=here("ITS2","sym.shade.ss.site.pdf"),h=10,w=15)    
# 
# #full plot - to put in supplementary, with all data across sites and time
# #numbers are individual samples here
# ss_shade_time <- plot_microshades(mdf_ss, cdf_test, x = "number", group_label = "Clade - Majority ITS2") + 
#   theme_classic(base_size = 30)+
#   xlab("Site & Timepoint")+
#   ylab("Relative Abundance")+
#   facet_wrap(m_y~site_nice, scales = "free") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle(substitute(paste(italic("Siderastrea siderea"))))
# ss_shade_time
# ggsave(ss_shade_time,file=here("ITS2","sym.shade.ss.site.time.pdf"),h=15,w=25)    
# 
# #srad
# sr_shade <- plot_microshades(mdf_sr, cdf_test, group_label = "Clade - Majority ITS2") + 
#   theme_classic(base_size = 30)+
#   ylab("Relative Abundance")+
#   facet_wrap(~site_nice, scales = "free")+
#   theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
#   ggtitle(substitute(paste(italic("Siderastrea radians"))))
# sr_shade
# ggsave(sr_shade,file=here("ITS2","sym.shade.sr.site.pdf"),h=5,w=15)    
# 
# #full plot - to put in supplementary, with all data acrosr sites and time
# #numbers are individual samples here
# sr_shade_time <- plot_microshades(mdf_sr, cdf_test, x = "number", group_label = "Clade - Majority ITS2") + 
#   theme_classic(base_size = 30)+
#   xlab("Site & Timepoint")+
#   ylab("Relative Abundance")+
#   facet_wrap(m_y~site_nice, scales = "free", ncol = 2) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle(substitute(paste(italic("Siderastrea radians"))))
# sr_shade_time
# ggsave(sr_shade_time,file=here("ITS2","sym.shade.sr.site.time.pdf"),h=15,w=15)
# 
# #ppor
# pp_shade <- plot_microshades(mdf_pp, cdf_test, group_label = "Clade - Majority ITS2") + 
#   theme_classic(base_size = 30)+
#   ylab("Relative Abundance")+
#   facet_wrap(~site_nice, scales = "free")+
#   theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
#   ggtitle(substitute(paste(italic("Porites spp."))))
# pp_shade
# ggsave(pp_shade,file=here("ITS2","sym.shade.pp.site.pdf"),h=5,w=15)
# 
# #full plot - to put in supplementary, with all data acropp sites and time
# #numbers are individual samples here
# pp_shade_time <- plot_microshades(mdf_pp, cdf_test, x = "number", group_label = "Clade - Majority ITS2") + 
#   theme_classic(base_size = 30)+
#   xlab("Site & Timepoint")+
#   ylab("Relative Abundance")+
#   facet_wrap(m_y~site_nice, scales = "free") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle(substitute(paste(italic("Porites spp."))))
# pp_shade_time
# ggsave(pp_shade_time,file=here("ITS2","sym.shade.pp.site.time.pdf"),h=10,w=15)  
# 
# #put plots together
# sym_pp_sr_plot <- ggarrange(sr_shade,pp_shade, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 40))
# sym_all_plot <- ggarrange(ss_shade,sym_pp_sr_plot, nrow=1, ncol = 2, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 40))
# ggsave(sym_all_plot, file=here("ITS2","sym_all_plot_microshades.pdf"),width=25,height=12)
# 
# #put plots together
# sym_all_time <- ggarrange(ss_shade_time,
#                           ggarrange(sr_shade_time, pp_shade_time, ncol = 2, 
#                                     labels = c("B", "C"), font.label = list(size = 50),
#                                     legend = FALSE), 
#                           #heights = c(1,2), widths = c(2,1), 
#                           nrow = 2,
#                           legend = "right", common.legend = TRUE, 
#                           labels = c("A"),font.label = list(size = 50))
# ggsave(sym_all_time, file=here("ITS2","sym_all_time_plot_supp_microshades.pdf"),width=25,height=30)
# 
# #SUPPLEMENTARY PLOTS OF ALL POST MED SEQUENCES
# mdf_prep <- prep_mdf(ps.nd, subgroup_level = "ITS2_type")
# # Create a color object for the specified data
# color_obj_prep_test <- create_color_dfs(mdf_prep, selected_groups = c("D","C","B","A"), group_level = "Clade", subgroup_level = "ITS2_type", cvd = TRUE)
# color_obj_prep_test$mdf$ITS2_type <- as.character(color_obj_prep_test$mdf$ITS2_type)
# # Extract
# mdf_test <- color_obj_prep_test$mdf
# cdf_test <- color_obj_prep_test$cdf
# 
# #subset different groups out
# mdf_ss <- mdf_test %>% filter(host_species == "siderea")
# mdf_sr <- mdf_test %>% filter(host_species == "radians")
# mdf_pp <- mdf_test %>% filter(host_species == "porites")
# 
# #now make graphs
# 
# #ssid
# ss_shade <- plot_microshades(mdf_ss, cdf_test, group_label = "Clade - ITS2") + 
#   theme_classic(base_size = 30)+
#   ylab("Relative Abundance")+
#   facet_wrap(~site_zone, scales = "free")+
#   theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
#   ggtitle(substitute(paste(italic("Siderastrea siderea"))))
# ss_shade
# ggsave(ss_shade,file=here("ITS2","sym.shade.ss.all.site.pdf"),h=10,w=15)    
# 
# #full plot - to put in supplementary, with all data across sites and time
# #numbers are individual samples here
# ss_shade_time <- plot_microshades(mdf_ss, cdf_test, x = "number", group_label = "Clade - ITS2 type") + 
#   theme_classic(base_size = 30)+
#   xlab("Coral Genotype")+
#   ylab("Relative Abundance")+
#   facet_wrap(m_y~site_zone, scales = "free") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle(substitute(paste(italic("Siderastrea siderea"))))
# ss_shade_time
# ggsave(ss_shade_time,file=here("ITS2","sym.shade.ss.all.site.time.pdf"),h=15,w=25)    
# 
# #srad
# sr_shade <- plot_microshades(mdf_sr, cdf_test, group_label = "Clade - ITS2 type") + 
#   theme_classic(base_size = 30)+
#   ylab("Relative Abundance")+
#   facet_wrap(~site_zone, scales = "free")+
#   theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
#   ggtitle(substitute(paste(italic("Siderastrea radians"))))
# sr_shade
# ggsave(sr_shade,file=here("ITS2","sym.shade.sr.all.site.pdf"),h=5,w=15)    
# 
# #full plot - to put in supplementary, with all data acrosr sites and time
# #numbers are individual samples here
# sr_shade_time <- plot_microshades(mdf_sr, cdf_test, x = "number", group_label = "Clade - ITS2 type") + 
#   theme_classic(base_size = 30)+
#   xlab("Coral Genotype")+
#   ylab("Relative Abundance")+
#   facet_wrap(m_y~site_zone, scales = "free", ncol = 2) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle(substitute(paste(italic("Siderastrea radians"))))
# sr_shade_time
# ggsave(sr_shade_time,file=here("ITS2","sym.shade.sr.all.site.time.pdf"),h=15,w=15)    
# 
# #ppor
# pp_shade <- plot_microshades(mdf_pp, cdf_test, group_label = "Clade - ITS2 type") + 
#   theme_classic(base_size = 30)+
#   ylab("Relative Abundance")+
#   facet_wrap(~site_zone, scales = "free")+
#   theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
#   ggtitle(substitute(paste(italic("Porites sp."))))
# pp_shade
# ggsave(pp_shade,file=here("ITS2","sym.shade.pp.all.site.pdf"),h=5,w=15)    
# 
# #full plot - to put in supplementary, with all data acropp sites and time
# #numbers are individual samples here
# pp_shade_time <- plot_microshades(mdf_pp, cdf_test, x = "number", group_label = "Clade - ITS2 type") + 
#   theme_classic(base_size = 30)+
#   xlab("Coral Genotype")+
#   ylab("Relative Abundance")+
#   facet_wrap(m_y~site_zone, scales = "free") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle(substitute(paste(italic("Porites spp."))))
# pp_shade_time
# ggsave(pp_shade_time,file=here("ITS2","sym.shade.pp.all.site.time.pdf"),h=10,w=15)    
# 
# #put plots together
# sym_pp_sr_plot <- ggarrange(sr_shade,pp_shade, nrow=2, ncol = 1, legend = "none", labels = c("B","C"),font.label = list(size = 40))
# sym_all_plot <- ggarrange(ss_shade,sym_pp_sr_plot, nrow=1, ncol = 2, legend = "right", common.legend = TRUE, labels = c("A"),font.label = list(size = 40))
# ggsave(sym_all_plot, file=here("ITS2","sym_all_plot_ALL_postmed_microshades.pdf"),width=25,height=12)
# 
# #put plots together
# sym_all_time <- ggarrange(ss_shade_time,
#                           ggarrange(sr_shade_time, pp_shade_time, ncol = 2, 
#                                     labels = c("B", "C"), font.label = list(size = 50),
#                                     legend = FALSE), 
#                           #heights = c(1,2), widths = c(2,1), 
#                           nrow = 2,
#                           legend = "right", common.legend = TRUE, 
#                           labels = c("A"),font.label = list(size = 50))
# ggsave(sym_all_time, file=here("ITS2","sym_all_time_plot_supp_ALL_postmed_microshades.pdf"),width=25,height=30)


####BETA DIVERSITY PLOTS#####
# #type profiles by site
# 
# #siderastrea siderea
# ps.ss.ord <- ordinate(ps.ss,"NMDS",distance="bray")
# gg.ss.ord <- plot_ordination(ps.ss, ps.ss.ord, color ="reef_bay", shape="area")+
#   geom_point(alpha=0.5)+
#   stat_ellipse(aes(linetype=area))+
#   theme_cowplot()
# gg.ss.ord
# ggsave(gg.ss.ord,file=here("ITS2","gg.ss.ord.pdf"),h=10,w=10) 
# 
# #siderastrea radians
# ps.sr.ord <- ordinate(ps.sr,"PCoA",distance="bray")
# gg.sr.ord <- plot_ordination(ps.sr, ps.sr.ord, color ="reef_bay", shape="area")+
#   geom_point(alpha=0.5)+
#   stat_ellipse(aes(linetype=area))+
#   theme_cowplot()
# gg.sr.ord
# ggsave(gg.sr.ord,file=here("ITS2","gg.sr.ord.pdf"),h=10,w=10) 
# 
# #porites porites
# ps.pp.ord <- ordinate(ps.pp,"PCoA",distance="bray")
# gg.pp.ord <- plot_ordination(ps.pp, ps.pp.ord, color ="reef_bay", shape="area")+
#   geom_point(alpha=0.5)+
#   stat_ellipse(aes(linetype=area))+
#   theme_cowplot()
# gg.pp.ord
# ggsave(gg.pp.ord,file=here("ITS2","gg.pp.ord.pdf"),h=10,w=10) 
# 
# 
# #porites astreoides
# ps.pa.ord <- ordinate(ps.pa,"PCoA",distance="bray")
# gg.pa.ord <- plot_ordination(ps.pa, ps.pa.ord, color ="reef_bay", shape="area")+
#   geom_point(alpha=0.5)+
#   stat_ellipse(aes(linetype=area))+
#   theme_cowplot()
# gg.pa.ord
# ggsave(gg.pa.ord,file=here("ITS2","gg.pa.ord.pdf"),h=10,w=10) 

# #####BETA DIVERSITY STATS####
# ###Stats###
# library(vegan)
# #remotes::install_github("Jtrachsel/funfuns")
# library(funfuns)
# library(dplyr)
# #BiocManager::install("edgeR")
# library(edgeR)
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# library(pairwiseAdonis)
# 
# #separate by species and only look at site - not timepoint
# 
# #siderastrea siderea
# seq.ss.ps <- data.frame(ps.ss@otu_table)
# samdf.ss.ps <- data.frame(ps.ss@sam_data)
# dist.ss.ps <- vegdist(seq.ss.ps, method = "bray")
# bet.ss.ps <- betadisper(dist.ss.ps,samdf.ss.ps$site)
# anova(bet.ss.ps) #significant
# set.seed(123)
# permutest(bet.ss.ps,pairwise=TRUE,permutations=999)
# plot(bet.ss.ps)
# 
# adonis2(seq.ss.ps ~ site_zone, data=samdf.ss.ps, permutations=999) #p<001***
# pairwise.adonis2(seq.ss.ps ~ site_zone, data = samdf.ss.ps, permutations=999)
# #               pairs   F.Model         R2 p.value p.adjusted
# #1  SW_bay vs SW_reef 18.560955 0.34018751   0.001      0.001
# #2   SW_bay vs SM_bay  4.861360 0.11342057   0.003      0.003
# #3  SW_bay vs SM_reef 11.402678 0.22183042   0.001      0.001
# #4  SW_reef vs SM_bay  6.496329 0.16041772   0.002      0.002
# #5 SW_reef vs SM_reef  5.516496 0.13287479   0.001      0.001
# #6  SM_bay vs SM_reef  3.704101 0.08881862   0.002      0.002
# 
# #siderastrea radians
# seq.sr.ps <- data.frame(ps.sr@otu_table)
# samdf.sr.ps <- data.frame(ps.sr@sam_data)
# dist.sr.ps <- vegdist(seq.sr.ps)
# bet.sr.ps <- betadisper(dist.sr.ps,samdf.sr.ps$site)
# anova(bet.sr.ps) #ns
# set.seed(123)
# permutest(bet.sr.ps,pairwise=TRUE,permutations=999)
# #       SMB   SWB
# #SMB         0.189
# #SWB 0.18582
# plot(bet.sr.ps)
# adonis2(seq.sr.ps ~ site, data=samdf.sr.ps, permutations=999) #ns
# pairwise.adonis2(seq.sr.ps, factors=samdf.sr.ps$site, permutations=999) #
# #pairs  F.Model         R2 p.value p.adjusted
# #1 SMB vs SWB 1.696757 0.04623724   0.128      0.128
# 
# 
# #porites porites
# seq.pp.ps <- data.frame(ps.pp@otu_table)
# samdf.pp.ps <- data.frame(ps.pp@sam_data)
# dist.pp.ps <- vegdist(seq.pp.ps)
# bet.pp.ps <- betadisper(dist.pp.ps,samdf.pp.ps$site)
# anova(bet.pp.ps) #very sig!!
# set.seed(123)
# permutest(bet.pp.ps,pairwise=TRUE,permutations=999)
# #         SMB   SMR
# #SMB            0.001
# #SMR 8.5674e-05   
# plot(bet.pp.ps)#+theme(legend.margin=0.01)
# adonis2(seq.pp.ps ~ site, data=samdf.pp.ps, permutations=999) #p<001***
# pairwise.adonis(seq.pp.ps, factors=samdf.pp.ps$site, permutations=999)
# #       pairs  F.Model        R2 p.value p.adjusted
# #1 SMR vs SMB 11.18058 0.2650646   0.001      0.001
# 
# ###Now look at timepoint as a factor & split them up by site
# 
# #siderastrea siderea
# ps.ss.SWB <- subset_samples(ps.ss, site_zone=="SW_bay")
# ps.ss.SWR <- subset_samples(ps.ss, site_zone=="SW_reef")
# ps.ss.SMB <- subset_samples(ps.ss, site_zone=="SM_bay")
# ps.ss.SMR <- subset_samples(ps.ss, site_zone=="SM_reef")
# 
# seq.ss.site.ps <- data.frame(ps.ss.SWB@otu_table)
# samdf.ss.site.ps <- data.frame(ps.ss.SWB@sam_data)
# dist.ss.site.ps <- vegdist(seq.ss.site.ps)
# bet.ss.site.ps <- betadisper(dist.ss.site.ps,samdf.ss.site.ps$m_y)
# anova(bet.ss.site.ps) #SWB (.), SWR (ns), SMB (ns), SMR (ns)
# set.seed(123)
# permutest(bet.ss.site.ps,pairwise=TRUE,permutations=999)
# 
# plot(bet.ss.site.ps)
# adonis2(seq.ss.site.ps ~ m_y, data=samdf.ss.site.ps, permutations=999) #SWB (.) SWB (ns) SMB (ns) SMR (ns)
# set.seed(123)
# pairwise.adonis(seq.ss.site.ps, factors=samdf.ss.site.ps$m_y, permutations=999)
# 
# #siderastrea radians
# ps.sr.SWB <- subset_samples(ps.sr, site_zone=="SW_bay")
# ps.sr.SMB <- subset_samples(ps.sr, site_zone=="SM_bay")
# 
# seq.sr.site.ps <- data.frame(ps.sr.SMB@otu_table)
# samdf.sr.site.ps <- data.frame(ps.sr.SMB@sam_data)
# dist.sr.site.ps <- vegdist(seq.sr.site.ps)
# bet.sr.site.ps <- betadisper(dist.sr.site.ps,samdf.sr.site.ps$m_y)
# anova(bet.sr.site.ps) #SWB (*), SMB (ns)
# set.seed(123)
# permutest(bet.sr.site.ps,pairwise=TRUE,permutations=999)
# 
# plot(bet.sr.ps)
# adonis2(seq.sr.site.ps ~ m_y, data=samdf.sr.site.ps, permutations=999) #SWB (ns) SMB ()
# pairwise.adonis(seq.sr.site.ps, factors=samdf.sr.site.ps$m_y, permutations=999) #
# #pairs  F.Model         R2 p.value p.adjusted
# #1 SMB vs SWB 1.696757 0.04623724   0.128      0.128
# 
# #porites porites
# ps.pp.SMR <- subset_samples(ps.pp, site_zone=="SM_reef")
# ps.pp.SMB <- subset_samples(ps.pp, site_zone=="SM_bay")
# 
# seq.pp.site.ps <- data.frame(ps.pp.SMR@otu_table)
# samdf.pp.site.ps <- data.frame(ps.pp.SMR@sam_data)
# dist.pp.site.ps <- vegdist(seq.pp.site.ps)
# bet.pp.site.ps <- betadisper(dist.pp.site.ps,samdf.pp.site.ps$m_y)
# anova(bet.pp.site.ps) #SMR ns
# set.seed(123)
# permutest(bet.pp.site.ps,pairwise=TRUE,permutations=999)
# plot(bet.pp.site.ps)
# adonis2(seq.pp.site.ps ~ m_y, data=samdf.pp.site.ps, permutations=999) #SMR ns
# pairwise.adonis(seq.pp.site.ps, factors=samdf.pp.site.ps$m_y, permutations=999)
# 
