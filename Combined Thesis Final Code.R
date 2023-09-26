#*************** loading packages ******************
setwd("G:/My Drive/UHasselt/Second Year/Semester Two/Master Thesis/Work/THESIS")
library(microbiome)
library(gridExtra)
library(tidyverse)
library(OTUtable)
library(magrittr)
library(phyloseq)
library(MMUPHin)
library(ANCOMBC)
library(miaTime)
library(scater)
library(readxl)
library(miaViz)
library(readr)
library(vegan)
library(caret)
library(mia)
library(DT)

#*************** Data loading *************************************
#******************************************************************
df <- read.csv("Data/Deos_file.csv", header = TRUE)
case_control <- read.csv("Data/preex_pcos.tsv", header = TRUE, sep = '')
SweMaMi_24FEB23 <- read_excel("Data/SweMaMi_24FEB23.xlsx")
fecal0 <- read.csv("Data/Fecal_T1_taxa.csv")
dim(fecal0) # 1446 1320


# Obtaining the taxon data and meta data
df_meta <- df[, c(2:11, 322:379)]
df_taxon0 <- df[, c(2, 12:321)]
rownames(df_taxon0) <- df_meta$Studienummer
rownames(df_meta) <- df_meta$Studienummer
dim(df_taxon0); # 3690  311 
dim(df_meta); # 3690 68
dim(case_control) # 329 3
dim(SweMaMi_24FEB23) # 5449  590

#*************** Analysis for vaginal microbiome *****************
#*****************************************************************
# Subset the microbiome data including only our case controls
case_control_2 <- case_control
case_control_2$extra <- c(1:329)
case_control_3 <- gather(case_control_2,"extra","key", 1:3)

case_control_taxa <- merge(case_control_3, df_taxon0,
                           by.x = "key", by.y = "Studienummer",
                           all.x = TRUE, all..y = FALSE)

case_control_taxa <- na.omit(case_control_taxa)
case_control_taxa$extra <- sapply(case_control_taxa$extra,
                                  function(x) gsub("1", "", x))
case_control_taxa$extra <- sapply(case_control_taxa$extra,
                                  function(x) gsub("2", "", x))
dim(case_control_taxa); # 915 312
df_taxon1 <- t(case_control_taxa)
dim(df_taxon1) # 312 915

# Filtering the samples
# To make a table containing only OTUs with at least 0.1% abundance 
# in at least one sample that were observed
# (at any abundance) in at least 10% of samples:
df1 <- df_taxon1[-c(1:2), ]
df2 <- matrix(as.numeric(df1), ncol = 915)
rownames(df2) <- rownames(df1)
df_taxon_filt <- filter_taxa(df2, abundance = 0.1, persistence = 10)
df_taxas <- as.data.frame(rbind(df_taxon1[c(1:2),], df_taxon_filt))
colnames(df_taxas) <- case_control_taxa$key
dim(df_taxas) # 176 915

# Meta Data for the case/control samples
df_meta_2 <- df_meta %>% filter(Studienummer %in% case_control_taxa$key)
df_meta_3 <- as.data.frame(cbind(case_control= df_taxon1[2,], df_meta_2))
rownames(df_meta_3) <- df_meta_2$Studienummer
dim(df_meta_3)
swema <- SweMaMi_24FEB23
swema1 <- swema %>% dplyr::filter(Studienummer %in% rownames(df_meta_3))
df_meta_3 <- as.data.frame(cbind(df_meta_3, 
                                 Gest_week = as.numeric(swema1$Gest_week_C), 
                                 Age = as.numeric(swema1$Q1_Age)))
dim(df_meta_3) # 915  71

# Zero cells (0.5813897)- 58.1%
sum(df_taxas == 0)/(nrow(df_taxas)*ncol(df_taxas))

# Checking missing data in PCOS variable
# filtering out missing values
sum(is.na(df_meta_3$pcos_score)) # 8
df_metadata <- df_meta_3 %>% filter(!is.na(pcos_score), !is.na(Gest_week),
                                    # pcos_score %in% c('Yes', 'No'),
                                    Age_groups_2 %in% c(1,2,3,4,5),
                                    parity %in% c('Multiparous', 'Nulliparous'),
                                    Pregnancy_week %in% c('Under 8','8 to 11',
                                                          '12 to 15','16 to 20')) 
df_metadata$Age_groups_2 <- as.factor(df_metadata$Age_groups_2)
df_metadata$Subfertility_ <- as.factor(df_metadata$Subfertility_)
df_metadata <- df_metadata %>% mutate(PCOS = pcos_score)
# The taxon data based on samples in metadata
df_taxon <- df_taxas[-c(1:2), rownames(df_metadata)]

# taxon data and metadata dimension
dim(df_taxon) # 174 778
dim(df_metadata) # 778  72

# Total sum of 778 samples
df_taxon2 <- matrix(as.numeric(unlist(df_taxon)), ncol = 778)
total_sum778 <- apply(df_taxon2, MARGIN = 2, sum)

# phyloseq data type
base::colnames(df_taxon2) <- base::rownames(df_metadata)
base::rownames(df_taxon2) <- base::rownames(df_taxon)
OTU2  <- otu_table(df_taxon2, taxa_are_rows = TRUE)
meta2 <- sample_data(df_metadata)
df_phyloseq <- merge_phyloseq(OTU2, meta2)

# Create the tse object
assays = S4Vectors::SimpleList(counts = df_taxon2)
meta = S4Vectors::DataFrame(df_metadata)
tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, 
                                                         colData = meta)
# Add clr-transformation on samples
tse <- transformCounts(tse, method = "relabundance")
tse <- transformCounts(tse, assay_name = "relabundance",
                       method = "clr", pseudocount = 0.000000001)
df_reabund <- assay(tse,"relabundance")
df_clr <- assay(tse,"clr")


#*************** Analysis for Faecal microbiome *****************
#*****************************************************************
# Subset the microbiome data including only our case controls
fecal_col_name <- colnames(fecal0)
fecal_col_names <- data.frame(Fecal_names = fecal_col_name[-c(1:4)])
fecal_col_names_modified <- lapply(fecal_col_names,
                                   function(x) str_split(x,"__"))

names <- c()
for (i in 1:1316) {
  if(!is.na(fecal_col_names_modified$Fecal_names[[i]][8])){
    names[i] <- fecal_col_names_modified$Fecal_names[[i]][8]} else
    {names[i] <- 'unknown'}
}

col_names_fecal <- c(fecal_col_name[1:4], names)
colnames(fecal0) <- col_names_fecal
#view(fecal0)

fecal <- fecal0 %>% select(-unknown)

#***************** Data For The Analysis******************************************
#*****************************************************************
# saveRDS(df_phyloseq5, 'Data/df_phyloseq5.rds')
# Vaginal Microbiome
df_metadata <- readRDS('Data/df_metadata.rds')
df_taxon2 <- readRDS('Data/df_taxon2.rds')
df_reabund <- readRDS('Data/df_reabund.rds')
df_phyloseq <- readRDS('Data/df_phyloseq.rds')
df_phyloseq5 <- readRDS('Data/df_phyloseq5.rds')
tse <- readRDS('Data/tse.rds')
tse4 <- readRDS('Data/tse4.rds')
# Faecal Microbiome
df_metadata <- readRDS('Data/df_metadata33.rds')
df_taxon2 <- readRDS('Data/df_taxon33.rds')
df_reabund <- readRDS('Data/df_reabund33.rds')
df_phyloseq <- readRDS('Data/df_phyloseq33.rds')
df_phyloseq5 <- readRDS('Data/df_phyloseq35.rds')
tse <- readRDS('Data/tse33.rds')
tse4 <- readRDS('Data/tse34.rds')

## Top Features
top_taxa_faecal <- getTopTaxa(tse, top = 5)
top_taxa_vaginal <- getTopTaxa(tse, top = 6)
#****************** Descriptive statistics for study variables *************
#***************************************************************************
library(table1)

df2 <- df_metadata
df2 <- df2 %>% mutate(Age_groups_2 = factor(Age_groups_2, levels = c(1,2,3,4,5), 
                                            labels = c('<= 24', '25-29','30-34', 
                                                       '35-39', '≥40')),
                      Subfertility_ = factor(Subfertility_, levels = c(0, 1), 
                                             labels = c("No", "Yes")))

# Set up template for tables
rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.default(x))
  what <- switch(name,
                 Age  = "Mean (SD)",
                 Gest_week = "Mean (SD)")
  parse.abbrev.render.code(c("", what))(x)
}

#----- Table 1: Participant’s characteristics
strata <- c(split(df2, df2$PCOS), list(Total = df2))

labels <- list(
  variables = list(Age = "Age (year)",
                   parity = "Parity",
                   BMI_category = "BMI category",
                   Subfertility_ = "Fertility treatment",
                   Age_groups_2 = "Age groups",
                   Gest_week = "Gestational length",
                   Conception_method = "Conception method"),
  groups = list("PCOS", "", ""))

table1 <- table1(strata, labels, groupspan = c(2, 2, 1), render = rndr, 
                 caption = "Descriptive statistics for study variables")
print(table1)
xtable::xtable(as.data.frame(table1), digits = 4)

# Statistical Analysis of characteristics on the two groups (PCOS)
test_age <- wilcox.test(df_metadata$Age~df_metadata$PCOS) 
test_age # p-value = 0.2981
tab_parity <- table(df_metadata$parity,df_metadata$PCOS)
chisq.test(tab_parity) # p-value = 0.3573
tab_BMI <- table(df_metadata$BMI_category,df_metadata$PCOS)
chisq.test(tab_BMI) # p-value = 0.2916
tab_fert <- table(df_metadata$Subfertility_,df_metadata$PCOS)
chisq.test(tab_BMI) # p-value = 0.2916
tab_age <- table(df_metadata$Age_groups_2, df_metadata$PCOS)
chisq.test(tab_age) # p-value = 0.9991
test_Gest <- wilcox.test(df_metadata$Gest_week~df_metadata$PCOS)
test_Gest # p-value = 0.7246
tab_Conc <- table(df_metadata$Conception_method, df_metadata$PCOS)
chisq.test(tab_Conc) # p-value = 6.052e-08

#************ Alpha Diversity ****************************************
#*********************************************************************
# Calculate alpha diversity 
df_metadata$simpson <- vegan::diversity(df_taxon2, MARGIN = 2,index = "simpson")
df_metadata %<>% mutate(shannon2 = microbiome::alpha(df_taxon2, index = "shannon"),
                        chao12 = microbiome::alpha(df_reabund, index = "chao1"),
                        pielou2 = microbiome::alpha(df_taxon2, index = "pielou"),
                        simpson2 = microbiome::alpha(df_taxon2, index = "simpson"))
dim(df_taxon2) # 174 778
dim(df_metadata) # 778  77

# Histogram of the two groups
hist(df_metadata$shannon2$diversity_shannon[df_metadata$PCOS == 'Yes'])
hist(df_metadata$shannon2$diversity_shannon[df_metadata$PCOS == 'No'])

# Homogeneity Test within the groups
# We conduct the Bartlett’s test of homogeneity of variances
bartlett.test(df_metadata$shannon2$diversity_shannon, df_metadata$PCOS) 
# p-value = 0.2229

# We now use Fligner-Killeen test to check the homoscedasticity. 
fligner.test(df_metadata$shannon2$diversity_shannon, df_metadata$PCOS) 
# p-value = 0.5046
mod <-lm(shannon2$diversity_shannon ~ PCOS+BMI_category+Gest_week, df_metadata) 
summary(mod)
# Do-Wilcoxon of Alpha diversity in groups
# Do-t.test of Alpha Diversity
t.test(simpson2$evenness_simpson ~ PCOS, df_metadata)
wilcox.test(simpson2$evenness_simpson ~ PCOS, df_metadata)
t.test(shannon2$diversity_shannon ~ PCOS, df_metadata) 
wilcox.test(shannon2$diversity_shannon ~ PCOS, df_metadata)
t.test(pielou2$evenness_pielou ~ PCOS, df_metadata) 
wilcox.test(pielou2$evenness_pielou ~ PCOS, df_metadata)
t.test(chao12$chao1 ~ PCOS, df_metadata) 
wilcox.test(chao12$chao1 ~ PCOS, df_metadata)
wilcox.test(shannon2$diversity_shannon ~ parity, df_metadata)
wilcox.test(shannon2$diversity_shannon ~ Subfertility_, df_metadata)
summary(aov(shannon2$diversity_shannon ~ Age_groups_2, df_metadata))
aov(shannon2$diversity_shannon ~ BMI_category, df_metadata)
kruskal.test(shannon2$diversity_shannon ~ Age_groups_2, df_metadata)
kruskal.test(shannon2$diversity_shannon ~ BMI_category, df_metadata)

#********************* Top Bacterial Species Analysis **********************
#***************************************************************************
top_taxa <- getTopTaxa(tse, top = 6)
relative_abund_table <- decostand(t(df_taxon2), method = "total")
p <- 0.000000001
df_metadata %<>% mutate(L.crispatus=log10(df_taxon2["Lactobacillus.crispatus",] + p),
                        L.iners=log10(df_taxon2["Lactobacillus.iners",] + p),
                        Gardnerella.vaginalis=log10(df_taxon2["Gardnerella.vaginalis",] + p),
                        L.jensenii=log10(df_taxon2["Lactobacillus.jensenii",] + p),
                        L.gasseri=log10(df_taxon2["Lactobacillus.gasseri",] + p),
                        Peptoniphilus.lacrimalis=log10(df_taxon2["Peptoniphilus.lacrimalis", ] + p))

grid.arrange(ggplot(df_metadata, aes(x = PCOS, y=L.crispatus, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -4.5, label = "P = 0.826", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'A) L.crispatus') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=L.iners, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -8, label = "P = 0.470", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'B) L.iners') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=Gardnerella.vaginalis, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -4.5, label = "P = 0.439", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'C) G.vaginalis') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=L.jensenii, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -7.5, label = "P = 0.570", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'D) L.jensenii') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=L.gasseri, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -7.5, label = "P = 0.579", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'E) L.gasseri') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=Peptoniphilus.lacrimalis, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -7.5, label = "P = 0.562", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'F) P.lacrimalis') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) , ncol = 3)

w1 <- wilcox.test(L.crispatus ~ PCOS, df_metadata) # p = 0.8262
w2 <- wilcox.test(L.iners ~ PCOS, df_metadata) # p = 0.47
w3 <- wilcox.test(Gardnerella.vaginalis ~ PCOS, df_metadata) # p = 0.4389
w4 <- wilcox.test(L.jensenii ~ PCOS, df_metadata) # p = 0.5704
w5 <- wilcox.test(L.gasseri ~ PCOS, df_metadata) # p = 0.5787
w6 <- wilcox.test(Peptoniphilus.lacrimalis ~ PCOS, df_metadata) # p = 0.5616

#***********************************************************************
#*        PCOS Vs Chao1 - Pielou - Shannon Visualization
#***********************************************************************
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x){
                m <- mean(x)
                ymin <- m-sd(x)
                ymax <- m+sd(x)
                return(c(y = m, ymin = ymin, ymax = ymax))}

A1 <- df_metadata %>% ggplot(aes(x = PCOS, y = chao12$chao1, col = PCOS)) +
      geom_violin(linewidth = 1.2) + theme_bw() + theme_classic() +
      labs(y = 'Chao1', x = 'PCOS', title = 'Chao1') +
      stat_summary(fun.data = data_summary) +
      geom_text(aes(label = 'P = 0.152', x = 1.5, y = 146),
                col = 'black', size = 7) +
      theme(legend.position="bottom", legend.justification="center",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.text = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 32, face = "bold")) 
A2 <- df_metadata %>% ggplot(aes(x = PCOS, y = shannon2$diversity_shannon, col = PCOS)) +
      geom_violin(linewidth = 1.2) + theme_bw()  + theme_classic() +
      labs(y = 'Shannon', x = 'PCOS', title = 'Shannon') +
      stat_summary(fun.data = data_summary) +
      geom_text(aes(label = 'P = 0.188', x = 1.5, y = 2.8),
                col = 'black', size = 7) +
      theme(legend.position="bottom", legend.justification="center",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.text = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold")) 
A3 <- df_metadata %>% ggplot(aes(x = PCOS, y = pielou2$evenness_pielou, col = PCOS)) +
      geom_violin(linewidth = 1.2) + theme_bw()  + theme_classic() +
      labs(y = 'Pielou', x = 'PCOS', title = 'Pielou') +
      stat_summary(fun.data = data_summary) +
      geom_text(aes(label = 'P = 0.831', x = 1.5, y = 0.8),
                col = 'black', size = 7) +
      theme(legend.position="bottom", legend.justification="center",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.text = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold")) 

Alpha_vaginal <- gridExtra::grid.arrange(A1, A3, A2, ncol = 3)

#************ Beta Diversity- ****************************************
#**********************************************************************
# Using Euclidean distance 
set.seed(2023)
#* Others # centered log-ratio (CLR) --- Considered
df_phyloseq_clr <- microbiome::transform(df_phyloseq, "clr")
taxa_bray <- phyloseq::distance(df_phyloseq_clr, method = "euclidean")
df_meta2 <- data.frame(sample_data(df_phyloseq))

# Homogeneity of dispersion test
beta <- betadisper(taxa_bray, df_meta2$pcos_score)
xtable::xtable(anova(beta), digits = 4)
tabb <- permutest(beta)
anova(beta)

# Adonis test for pcos_score (PERMANOVA)
set.seed(2023)
perm_pcos <- adonis2(taxa_bray ~ PCOS, data = df_meta2) 
perm_pcos # p = 0.459
set.seed(2023)
perm_gestW <- adonis2(taxa_bray ~ Gest_week, data = df_meta2) 
perm_gestW # p = 0.721
set.seed(2023)
perm_age <- adonis2(taxa_bray ~ Age_groups_2, data = df_meta2) 
perm_age # p = 0.166 (significant)
set.seed(2023)
perm_fert <- adonis2(taxa_bray ~ Subfertility_, data = df_meta2) 
perm_fert # p =  0.001 (significant)
set.seed(2023)
perm_BMI_c <- adonis2(taxa_bray ~ factor(BMI_category), data = df_meta2) 
perm_BMI_c # p = 0.001 (significant)
set.seed(2023)
perm_parity <- adonis2(taxa_bray ~ parity, data = df_meta2) 
perm_parity # p = 0.001 (significant)

#******************************************************************
#*          Data Exploration -- Beta Diversity
#******************************************************************
library(ComplexHeatmap)
library(microbiomeutilities)
library(viridis)
library(RColorBrewer)

ps_ord <- ordinate(df_phyloseq, method = "PCoA", distance = "bray")
PCoA11 <- plot_ordination(df_phyloseq, ps_ord, type = "samples", color = "PCOS")

ps_ordCCA <- ordinate(df_phyloseq, method = "CCA", distance = "bray")
CCA <- plot_ordination(df_phyloseq, ps_ordCCA, type = "samples", 
                       color = "PCOS")
ps_ordDCA <- ordinate(df_phyloseq, method = "DCA", distance = "bray")
DCA <- plot_ordination(df_phyloseq, ps_ordDCA, type = "samples", color = "PCOS")

CCA <-  CCA + geom_vline(xintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#ef3b2c", "#0570b0")) +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        theme_bw() + theme_classic() +
        labs(title = "Visualization of Vaginal Microbiome")+
        theme(legend.position="bottom", legend.justification="center",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
              legend.text = element_text(size = 14, face = "bold"),
              axis.text = element_text(size = 12, face = "bold"),
              axis.title = element_text(size = 16, face = "bold"),
              legend.title = element_text(size = 16, face = "bold"),
              text = element_text(size = 16, face = "bold"))

DCA <- DCA + geom_vline(xintercept = 0, linetype = "dashed") +
      scale_color_manual(values = c("#ef3b2c", "#0570b0")) +
      geom_hline(yintercept = 0, linetype = "dashed") + 
      theme_bw() + theme_classic() +
      theme(legend.position="bottom", legend.justification="center",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.text = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold"))  


#PCA via phyloseq
df_phyloseq_ra <- microbiome::transform(df_phyloseq,  'compositional')
ord_clr <- phyloseq::ordinate(df_phyloseq_clr, "RDA")
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

P11 <- phyloseq::plot_ordination(df_phyloseq, ord_clr,
                                 type = "samples", color="PCOS") + 
      theme_bw() + theme_classic() +
      geom_point(size = 2) +
      coord_fixed(clr2 / clr1) +
      stat_ellipse(aes(group = PCOS), linetype = 2) +
      geom_hline(yintercept = 0, linetype = 'dashed') +
      geom_vline(xintercept = 0, linetype = 'dashed') +
      theme(legend.position="bottom", legend.justification="center",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.text = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold"))

# Using MDS Ordination plots using tse data
tse4 <- runNMDS(tse4, FUN = vegan::vegdist, name = "NMDS")
tse4 <- runNMDS(tse4, FUN = vegan::vegdist, name = "PCA")
tse4 <- runMDS(tse4, FUN = vegan::vegdist, method = "bray",
               name = "MDS", exprs_values = "counts")

library(patchwork)
plots11 <- lapply(c("MDS", "NMDS", "PCA"), plotReducedDim,
                  object = tse4, colour_by = "PCOS")

PCoA11 + labs(title = 'A) PCoA', x = 'PCoA1 [43.5%]', y = 'PCoA2 [23.8%]') + 
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_color_manual(values = c("#ef3b2c", "#0570b0")) +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        theme_bw() + theme_classic() +
        theme(legend.position="bottom", legend.justification="center",
              plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
              legend.text = element_text(size = 14, face = "bold"),
              axis.text = element_text(size = 12, face = "bold"),
              axis.title = element_text(size = 16, face = "bold"),
              legend.title = element_text(size = 16, face = "bold"),
              text = element_text(size = 16, face = "bold")) +
plots11[[2]] + labs(title = 'B) NMDS', color = 'PCOS') +
              geom_vline(xintercept = 0, linetype = "dashed") +
              geom_hline(yintercept = 0, linetype = "dashed") +
              theme(legend.position="bottom", legend.justification="center",
                    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                    legend.text = element_text(size = 14, face = "bold"),
                    axis.text = element_text(size = 12, face = "bold"),
                    axis.title = element_text(size = 16, face = "bold"),
                    legend.title = element_text(size = 16, face = "bold"),
                    text = element_text(size = 16, face = "bold"))  + 
              scale_color_manual(values = c("#ef3b2c", "#0570b0")) +
DCA + labs(title = 'C) DCA', color = 'PCOS') 

#*******************************************************************
#*           RCM ------------ RCM for faecal microbiome
#*******************************************************************
library(RCM)
PCOS.rcm <- RCM(df_phyloseq, round = T, prevCutOff = 0)
PCOS2.rcm <- RCM(df_phyloseq, confounders = "PCOS", round = TRUE)
PCOS3.rcm <- RCM(df_phyloseq, covariates = "BMI_category", round = TRUE)
plot(PCOS3.rcm, samColour = "PCOS")

plot(PCOS.rcm, samColour = "PCOS") + plot(PCOS2.rcm, samColour = "PCOS")
gridExtra::grid.arrange(plot(PCOS.rcm, samColour = "PCOS"),
                        plot(PCOS2.rcm, samColour = "PCOS"),
                        ncol = 1)

#*****************************************************************
#*              Heatmap using tse data
#*****************************************************************
# Computing relative abundance
tse <- relAbundanceCounts(tse)
# Getting top taxa on a Phylum level
tse <- transformFeatures(tse, assay_name = "clr",
                         method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse, top = 20)
tse1 <- tse
tse1 <- tse1[top_taxa, ]
mat <- assay(tse1, "clr_z")

# Creates the heatmap
library(pheatmap)
library(ape)
library(ggtree)

# Hierarchical clustering
taxa_hclust <- hclust(dist(mat), method = "ward.D2")
# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)
# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin = margin(0,0,0,0)) 
# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)
taxa_tree

# Creates clusters
taxa_clusters <- cutree(tree = taxa_hclust, k = 7)
# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)
# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE]
# Prints taxa and their clusters
rowData(tse1)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters),
                                                    rownames(tse1))),]
# Hierarchical clustering
sample_hclust <- hclust(dist(t(mat)), method = "ward.D2")
# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)
# Plot sample tree
sample_tree <- ggtree(sample_tree) + 
  layout_dendrogram() +
  theme(plot.margin = margin(0,0,0,0)) 
# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))


# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 7))
# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)
# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE]
# Order data based on
tse1 <- tse1[, rownames(sample_data)]
# Add sample type data
sample_data$PCOS <- colData(tse1)$PCOS 
sample_data$BMI <- colData(tse1)$BMI_category
sample_data$Parity <- colData(tse1)$parity
cluster_data <- data.frame(clusters = sample_clusters)
sample_data %<>% select(-clusters)

# Determines the scaling of colorss
# Scale colors
breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))),
              length.out = ifelse(max(abs(mat))>5,
                                  2*ceiling(max(abs(mat))), 10))
colors <- colorRampPalette(c("darkblue", "blue", "white", "red",
                             "darkred"))(length(breaks)-1)
# Vaginal
col = list(PCOS = c("No" = "green","Yes" = "purple"),
           BMI = c("Normal weight" = "yellow", "Obese" = "darkblue",
                  "Overweight weight" = "darkred", "Underweight" = "purple"),
           Parity = c("Multiparous" = "blue", "Nulliparous" = "red") )
# faecal
col = list(PCOS = c("No" = "green","Yes" = "purple"),
           BMI = c("Normal weight" = "yellow", "Obese" = "darkblue",
                   "Overweight weight" = "darkred", "Underweight" = "purple"),
           `Age groups` = c("≤24" = "blue", "25-29" = "red", "30-34" = "purple",
                            "35-39" = "green", "40+" = "yellow") )
# CSTs
col = list(CST = c("I" = "blue", "II" = "red", "III" = "purple", 
                   "IV-B" = "green", "IV-C" = "yellow", "V" = "black") )

pheatmap(mat, annotation_row = taxa_clusters,
         annotation_col = sample_data,
         annotation_colors = col,
         breaks = breaks, color = colors,
         cutree_rows = 7, fontsize = 13.5,
         fontsize_row = 15, show_colnames = F,
         annotation_names_row = F)

#************Differential Abundance********************************
#******************************************************************
#*                   ANCOM-BC Differential Abundance 
PCOS.abc <- ancombc2(data = df_phyloseq, assay_name = "counts", 
                      fix_formula  = "PCOS + Subfertility_ + BMI_category",
                      p_adj_method = "BH",  prv_cut = 0,
                      lib_cut = 0, group = "PCOS", 
                      struc_zero = TRUE, neg_lb = TRUE, 
                      em_control = list(tol = 1e-5, max_iter = 20),
                      alpha = 0.05, global = TRUE)
results_New1 <- PCOS.abc$res
PCOS.abc2 <- ancombc2(data = df_phyloseq5, assay_name = "counts", 
                      fix_formula  = "PCOS", tax_level = "Species",
                      p_adj_method = "BH",  prv_cut = 0,
                      lib_cut = 0, group = "PCOS", 
                      struc_zero = TRUE, neg_lb = TRUE, 
                      em_control = list(tol = 1e-5, max_iter = 20),
                      alpha = 0.05, global = TRUE)
results_ABC2 <- PCOS.abc2$res
results_ABC3 <- data.frame(Species = rownames(results_ABC2),
                           lfc = results_ABC2$lfc_PCOSYes, 
                           se = results_ABC2$se_PCOSYes,
                           w = results_ABC2$W_PCOSYes,
                           Adj.Pval.PCOS = results_ABC2$q_PCOSYes)
xtable::xtable(results_ABC3, digits = 5)

#****************************************************************************
#*                   Wilcoxon Differential Abundance (Try other mtds)
OTU.TSS <- otu_table(df_phyloseq)
group <- meta(df_phyloseq)$PCOS

results.WMW <- data.frame(OTU = rownames(OTU.TSS), stat = NA, p_value = NA)

for(i in 1:nrow(OTU.TSS)) {
      tmp <- wilcox.test(
        x = as.numeric(OTU.TSS[i, group == "Yes"]),
        y = as.numeric(OTU.TSS[i, group == "No"]))
      
      results.WMW$stat[i] <- tmp$statistic
      results.WMW$p_value[i] <- tmp$p.value
}
diff_abund_Wil <- results.WMW %>%
                  arrange(p_value) %>%
                  mutate(p_Adjust = p.adjust(p_value, "BH")) %>%
                  filter(p_Adjust < 0.05)
nrow(diff_abund_Wil) # 0

#***************************************************************************
#*             Negative Binomial Differential Abundance (Try other mtds)
results.NB <- data.frame(OTU = rownames(OTU.TSS), lfc = NA, p_value = NA)

db <- data.frame(N=NA,
                 group = group,
                 LogLibSize = log(colSums(OTU.TSS)))

for(i in 1:nrow(OTU.TSS)) {
                db$N <- as.numeric(OTU.TSS[i,])
                m <- MASS::glm.nb(N ~ group + offset(db$LogLibSize), data = db)
                
                results.NB$lfc[i] <- coef(m)[2]
                results.NB$p_value[i] <- summary(m)$coef[2,4]
}
diff_abund_NB <- results.NB %>%
                arrange(p_value) %>%
                mutate(p_Adjust = p.adjust(p_value, "BH")) %>%
                filter(p_Adjust < 0.05)
nrow(diff_abund_NB) # 0


#**************************************************************************
#*               Community State Types (CSTs)
#**************************************************************************
test_df1 <- t(df_taxon2)
row_sum <- apply(test_df1, MARGIN = 1, sum)
test_df2 <- data.frame(sampleID = rownames(test_df1),
                       read_count = row_sum,
                       test_df1)
dim(test_df2) # 778 176

summarize_columns <- function(df){
  # Replace periods with underscores in column names
  names(df) <- gsub("\\.", "_", names(df))
  
  # Get a list of unique column prefixes
  prefixes <- unique(sub("^(.*?)_.*", "\\1", names(df)))
  main_taxa <- c("sampleID", "read", "Lactobacillus","Gardnerella",
                 "Prevotella", "Atopobium", "Sneathia")
  
  # Loop over prefixes and create a new column with the sum
  for (prefix in prefixes) {
    if (prefix %in% main_taxa) {
      # Leave the column as it is
    } else {
      cols_to_sum <- grep(paste0("^", prefix, "_"), names(df), value = TRUE)
      if(length(cols_to_sum) > 1){
        df[, paste0("g_", prefix)] <- rowSums(df[, cols_to_sum])
      } else{
        df[, paste0("g_", prefix)] <- df[, cols_to_sum]
      }
      df <- df[, !names(df) %in% cols_to_sum]
    }}
  
  # Return the summarized dataframe
  return(df)
}

cst_data <- summarize_columns(test_df2)

dim(cst_data) # 778 122
table(cst_data$CST)
table(cst_data$subCST)
top_taxa <- getTopTaxa(tse, top = 25)
df_taxon_long <- df_taxon2[top_taxa,]
df_taxon_long1 <- data.frame(t(df_taxon_long), CST = cst_data$CST)

data_long <- df_taxon_long1 %>%
            pivot_longer(cols = -CST, names_to = "Variable", 
                         values_to = "Value")
data_long1 <- data_long %>% group_by(Variable) %>% 
              summarise(Total = sum(Value)) 

# Plot the density of all variables against the categorical variable
ggplot(data_long, aes(x = CST, y = Value, fill = Variable)) +
  geom_bar(stat = "identity", width = 0.7) + theme_classic() +
  scale_fill_viridis(discrete=TRUE, name = "Species") +
  labs(x = "CSTs", y = "Total Abundance") + theme_bw()  + 
  theme(legend.position="bottom", legend.justification="center",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        text = element_text(size = 16, face = "bold")) + coord_flip()

#*********************************************************************
#*                      Heatmap of CSTs and species features
#*********************************************************************                      
taxon <- cst_data[, c(3:106)]
cst_taxon1 <- cst_data[, c(107:122)]
taxon1 <- t(taxon)
# taxon1 <- scale(taxon1)
colnames(taxon1) <- cst_data$sampleID
rownames(cst_taxon1) <- cst_data$sampleID
assays = S4Vectors::SimpleList(counts = df_taxon2)
meta = S4Vectors::DataFrame(cst_taxon1)
tse1 = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, 
                                                          colData = meta)
tse1 <- transformCounts(tse1, method = "relabundance")
tse1 <- transformCounts(tse1, assay_name = "relabundance",
                        method = "clr", pseudocount = 0.000000001)
tse1 <- transformCounts(tse1, assay_name = "relabundance",
                        method = "rclr", pseudocount = 0.000000001)
tse1 <- transformFeatures(tse1, assay_name = "clr",
                          method = "z", name = "clr_z")


top_taxa <- getTopTaxa(tse1, top = 25)
tse1 <- tse1[top_taxa, ]
taxon1 <- assay(tse1, "clr_z")
# Hierarchical clustering
taxa_hclust <- hclust(dist(taxon1), method = "ward.D2")
taxa_tree <- as.phylo(taxa_hclust)
taxa_tree <- ggtree(taxa_tree) + 
            theme(plot.margin = margin(0,0,0,0)) 
taxa_ordered <- get_taxa_name(taxa_tree)
taxa_clusters <- cutree(tree = taxa_hclust, k = 6)
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE]
rowData(tse1)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters),
                                                    rownames(tse1))),]
# Hierarchical clustering
sample_hclust <- hclust(dist(t(taxon1)), method = "ward.D2")
sample_tree <- as.phylo(sample_hclust)
sample_tree <- ggtree(sample_tree) + 
              layout_dendrogram() +
              theme(plot.margin = margin(0,0,0,0)) 
samples_ordered <- rev(get_taxa_name(sample_tree))


# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 6))
sample_data <- data.frame(clusters = sample_clusters)
sample_data <- sample_data[samples_ordered, , drop = FALSE]
tse1 <- tse1[, rownames(sample_data)]
sample_data$CST <- sort(colData(tse1)$CST )
sample_data %<>% select(-clusters)

# Annotation/coloring
breaks <- seq(-ceiling(max(abs(taxon1))), ceiling(max(abs(taxon1))),
              length.out = ifelse(max(abs(taxon1))>5,2*ceiling(max(abs(taxon1))), 10))
colors <- colorRampPalette(c("darkblue", "blue", "white", "red",
                             "darkred"))(length(breaks)-1)
pheatmap(taxon1, annotation_row = taxa_clusters,
         annotation_col = sample_data, breaks = breaks, 
         color = colors, cutree_rows = 4, display_numbers = F,
         fontsize = 13.5, fontsize_row = 15,
         annotation_legend = T)

#**********************************************************************
#*                Box-plot of Shannon Diversity Vs subCSTs
#*                
# Shannon Diversity Vs subCSTs
rownames(cst_taxon1) <- cst_data$sampleID
metaaa <- as.data.frame(cbind(df_metadata,cst_taxon1))

metaaa %>% ggplot(aes(y = shannon2$diversity_shannon,x = subCST, fill = subCST)) +
          geom_boxplot(show.legend = F) + theme_bw() + theme_classic() +
          labs( y = "Shannon diversity") +
          theme(legend.position="bottom", legend.justification="center",
                plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                legend.text = element_text(size = 14, face = "bold"),
                axis.text = element_text(size = 12, face = "bold"),
                axis.title = element_text(size = 16, face = "bold"),
                legend.title = element_text(size = 16, face = "bold"),
                text = element_text(size = 32, face = "bold"))  +
          coord_flip()

#****************** Descriptive statistics for CSTs **********************
#*************************************************************************
library(table1)

df2 <- metaaa
# Set up template for tables
rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.default(x))
  what <- switch(name,)
  parse.abbrev.render.code(c("", what))(x)
}

#----- Table 1: Participant’s characteristics
strata <- c(split(df2, df2$PCOS), list(Total = df2))

labels <- list(
  variables = list(CST = "CSTs"),
  groups = list("PCOS", "", ""))

table2 <- table1(strata, labels, groupspan = c(2, 2, 1), render = rndr, 
                 caption = "Descriptive statistics for CSTs")
print(table2)
xtable::xtable(as.data.frame(table2), digits = 4)


#******************************************************************
#*               Genus level Analysis for faecal microbiome
#******************************************************************
genus_level <- function(df) {
  prefixes <- unique(sub("^(.*?)_.*", "\\1", colnames(df)))
  
  # Loop over prefixes and create a new column with the sum
  for (prefix in prefixes) {
    cols_to_sum <- grep(paste0("^", prefix, "_"), colnames(df), value = TRUE)
    if(length(cols_to_sum) > 1){
      df[, paste0("g_", prefix)] <- rowSums(df[, cols_to_sum])
    } 
    else{
      df[, paste0("g_", prefix)] <- df[, cols_to_sum]
    }
    df <- df[, !colnames(df) %in% cols_to_sum]
  }
  
  # Return the summarized dataframe
  return(df)
}

df_taxon_g <- t(df_taxon)
df_taxon_gg <- read.csv("Data/df_taxon_g.csv")
rownames(df_taxon_gg) <- df_taxon_gg$X
df_taxon_gg %<>% select(-X)
df_genus <- genus_level(df_taxon_gg)

dim(df_taxon_g) # 738 537 
dim(df_genus) # 738 277
df_genus <- readRDS("Data/df_genus.rds")

# phyloseq data type
df_genus2 <- t(df_genus)
OTU2  <- otu_table(df_genus2, taxa_are_rows = TRUE)
meta2 <- sample_data(df_metadata)
df_phyloseq_genus <- merge_phyloseq(OTU2, meta2)

#*ANCOM-BC2 Differential Abundance 
PCOS.abc_genus <- ancombc2(df_phyloseq_genus, fix_formula  = "PCOS", 
                           p_adj_method = "BH",  prv_cut = 0,
                           lib_cut = 0, group = "PCOS", 
                           struc_zero = TRUE, neg_lb = TRUE, 
                           em_control = list(tol = 1e-5, max_iter = 20),
                           alpha = 0.05, global = TRUE)
results_ABC_genus <- PCOS.abc_genus$res
results_ABC3 <- data.frame(Genera = rownames(results_ABC_genus),
                           lfc = results_ABC_genus$lfc_PCOSYes, 
                           se = results_ABC_genus$se_PCOSYes,
                           w = results_ABC_genus$W_PCOSYes,
                           Adj.P.value = results_ABC_genus$q_PCOSYes)

xtable::xtable(head(results_ABC3, 6), digits = 5)












































