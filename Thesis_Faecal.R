setwd("G:/My Drive/UHasselt/Second Year/Semester Two/Master Thesis/Work/THESIS")
library(microbiome)
library(gridExtra)
library(tidyverse)
library(OTUtable)
library(magrittr)
library(GUniFrac)
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
library(RCM)
library(mia)
library(DT)

#*************** Data loading ******************
df <- read.csv("Data/Deos_file.csv", header = TRUE)
case_control <- read.csv("Data/preex_pcos.tsv", header = TRUE, sep = '')
SweMaMi_24FEB23 <- read_excel("Data/SweMaMi_24FEB23.xlsx")
fecal0 <- read.csv("Data/Fecal_T1_taxa.csv")
dim(fecal0) # 1446 1320

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

# Obtaining the taxon data and meta data
df_meta <- df[, c(2:11, 322:379)]
rownames(df_meta) <- df_meta$Studienummer
dim(df_meta); # 3690 68
dim(case_control) # 329 3
dim(SweMaMi_24FEB23) # 5449  590
dim(fecal) # 1446 1248


#*************** Analysis ******************
# Subset the microbiome data including only our case controls
case_control_2 <- case_control
case_control_2$extra <- c(1:329)
case_control_3 <- gather(case_control_2,"extra","key", 1:3)

case_control_fecal <- merge(case_control_3, fecal,
                            by.x = "key", by.y = "Studienummer",
                            all.x = TRUE, all..y = FALSE)
dim(case_control_fecal) # 1017 1249

case_control_fecal <- na.omit(case_control_fecal)
dim(case_control_fecal) # 967 1249
case_control_fecal$extra <- sapply(case_control_fecal$extra,
                                   function(x) gsub("1", "", x))
case_control_fecal$extra <- sapply(case_control_fecal$extra,
                                   function(x) gsub("2", "", x))
case_control_fecal <- as.data.frame(case_control_fecal)
dim(case_control_fecal); # 967 1249
df_taxon11 <- t(case_control_fecal)
dim(df_taxon11) # 1249  967

# Filtering the samples
# To make a table containing only OTUs with at least 0.1% abundance 
# in at least one sample that were observed
# (at any abundance) in at least 10% of samples:
df11 <- df_taxon11[-c(1:5), ]
df22 <- matrix(as.numeric(df11), ncol = 967)
rownames(df22) <- rownames(df11)
df_taxon_filt11 <- filter_taxa(df22, abundance = 0.1, 
                               persistence = 10)
df_taxas11 <- as.data.frame(rbind(df_taxon11[c(1:5),], df_taxon_filt11))
colnames(df_taxas11) <- case_control_fecal$key
dim(df_taxas11) # 542  967

# Meta Data for the case/control samples
df_meta_22 <- df_meta %>% dplyr::filter(Studienummer %in% colnames(df_taxas11))
df_taxas11 <- df_taxas11[, df_meta_22$Studienummer]
df_meta_33 <- as.data.frame(cbind(case_control= as.matrix(df_taxas11)[2,], df_meta_22))
rownames(df_meta_33) <- df_meta_22$Studienummer
dim(df_meta_33)
swema <- SweMaMi_24FEB23
swema11 <- swema %>% dplyr::filter(Studienummer %in% rownames(df_meta_33))
df_meta_33 <- as.data.frame(cbind(df_meta_33, 
                                  Gest_week = as.numeric(swema11$Gest_week_C), 
                                  Age = as.numeric(swema11$Q1_Age)))
df_meta_33 <- df_meta_33 %>% mutate(PCOS = pcos_score)
dim(df_meta_33) # 871  72

# Zero cells (0.5571998)- 55.7%
sum(df_taxas11 == 0)/(nrow(df_taxas11)*ncol(df_taxas11))

# Checking missing data in PCOS variable
# filtering out missing values
sum(is.na(df_meta_33$pcos_score)) # 8
df_metadata33 <- df_meta_33 %>% filter(!is.na(pcos_score), !is.na(Gest_week),
                                       # pcos_score %in% c('Yes', 'No'),
                                       Age_groups_2 %in% c(1,2,3,4,5),
                                       parity %in% c('Multiparous', 'Nulliparous'),
                                       Pregnancy_week %in% c('Under 8','8 to 11',
                                                             '12 to 15','16 to 20'))
df_metadata33$Age_groups_2 <- as.factor(df_metadata33$Age_groups_2)
df_metadata33$Subfertility_ <- as.factor(df_metadata33$Subfertility_)

# The taxon data based on samples in metadata
df_taxon11 <- df_taxas11[-c(1:5), rownames(df_metadata33)]

# df_taxon <- df_taxon11
# df_metadata <- df_metadata33
dim(df_taxon11) # 537 738
dim(df_metadata33) # 738  72
df_taxon <- df_taxon11
df_metadata <- df_metadata33

# Vaginal Microbiome
# Maybe = 10    No = 516  Yes = 262

# Faecal Microbiome
# Maybe = 12    No = 483  Yes = 255

# European Helicobacter and Microbiota Study Group; 5-7 September, Payment.

# Total sum of 738 samples
df_taxon <- matrix(as.numeric(unlist(df_taxon)), ncol = 738)
total_sum738 <- apply(df_taxon, MARGIN = 2, sum)

# phyloseq data type
colnames(df_taxon) <- rownames(df_metadata)
rownames(df_taxon) <- rownames(df_taxon11)
OTU2  <- otu_table(df_taxon, taxa_are_rows = TRUE)
meta2 <- sample_data(df_metadata)
df_phyloseq <- merge_phyloseq(OTU2, meta2)

# Create the tse object
assays = S4Vectors::SimpleList(counts = df_taxon)
meta = S4Vectors::DataFrame(df_metadata)
tse = TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, 
                                                         colData = meta)
# Add clr-transformation on samples
tse <- transformCounts(tse, method = "relabundance")
tse <- transformCounts(tse, assay_name = "relabundance",
                       method = "clr", pseudocount = 0.000000001)
tse <- transformCounts(tse, assay_name = "relabundance",
                       method = "rclr", pseudocount = 0.000000001)
df_reabund <- assay(tse,"relabundance")

# phyloseq with species features
OTU2  <- otu_table(df_taxon, taxa_are_rows = TRUE)
meta2 <- sample_data(df_metadata33)
taxa_data <- data.frame(Kingdom = rownames(df_taxon),
                        Phylum = rownames(df_taxon),
                        Class = rownames(df_taxon),
                        Order = rownames(df_taxon),
                        Family = rownames(df_taxon),
                        Genus = rownames(df_taxon),
                        Species = rownames(df_taxon),
                        row.names = rownames(df_taxon))
taxa_dat <- as.matrix(taxa_data)
taxa_tab2 <- tax_table(taxa_dat)
df_phyloseq35 <- merge_phyloseq(OTU2, meta2, taxa_tab2)


#*****************Data************************
# saveRDS(tse4, 'Data/tse34.rds')
df_metadata <- readRDS('Data/df_metadata33.rds')
df_taxon <- readRDS('Data/df_taxon33.rds')
df_reabund <- readRDS('Data/df_reabund33.rds')
df_phyloseq2 <- readRDS('Data/df_phyloseq33.rds')
df_phyloseq5 <- readRDS('Data/df_phyloseq35.rds')
tse <- readRDS('Data/tse33.rds')
#tse4 <- readRDS('Data/tse34.rds')

#******************Descriptive statistics for study variables*************
#*************************************************************************
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
                   Gest_week = "Gestation weeks"),
  groups = list("PCOS", "", ""))

table1 <- table1(strata, labels, groupspan = c(2, 2, 1), render = rndr, 
                 caption = "Descriptive statistics for study variables")
print(table1)
xtable::xtable(as.data.frame(table1), digits = 4)

#************Alpha Diversity***************************************
#*******************************************************************
# Calculate alpha diversity (Shannon-index)
# df_taxon <- matrix(as.numeric(unlist(df_taxon)), ncol = 778)
df_metadata$simpson <- vegan::diversity(df_taxon, MARGIN = 2,index = "simpson")
df_metadata %<>% mutate(shannon2 = microbiome::alpha(df_taxon, index = "shannon"),
                        chao12 = microbiome::alpha(df_taxon, index = "chao1"),
                        pielou2 = microbiome::alpha(df_taxon, index = "pielou"),
                        simpson2 = microbiome::alpha(df_taxon, index = "simpson"))
dim(df_taxon) # 537 738
dim(df_metadata) # 738  77

# Do-ANOVA
anova_result_pcos1 <- aov(simpson2$evenness_simpson ~ PCOS, df_metadata) 
summary(anova_result_pcos1) # P = 0.0296 (significant)
anova_result_pcos2 <- aov(shannon2$diversity_shannon ~ PCOS, df_metadata) 
summary(anova_result_pcos2) # P = 0.136
anova_result_pcos3 <- aov(pielou2$evenness_pielou ~ PCOS, df_metadata) 
summary(anova_result_pcos3) # P = 0.813
anova_result_pcos4 <- aov(chao12$chao1 ~ PCOS, df_metadata) 
summary(anova_result_pcos4) # P =  0.128
anova_result_Age <- aov(shannon2$diversity_shannon ~ Age_groups_2, df_metadata)
summary(anova_result_Age) # P = 0.818 
anova_result_BMI <- aov(shannon2$diversity_shannon ~ BMI_category, df_metadata)
summary(anova_result_BMI) # P = 0.000776 (significant)
anova_result_parity <- aov(shannon2$diversity_shannon ~ parity, df_metadata)
summary(anova_result_parity) # P = 0.746 

Alph_F <- data.frame(Parameters = c('Simpson', 'Shannon','Pielou', 'Chao1'),
                     F_value1 = c(4.752, 2.229, 0.056, 2.321),
                     P_value1 = c(0.0296, 0.136, 0.813, 0.128),
                     F_value2 = c(4.752, 2.229, 0.056, 2.321),
                     P_value2 = c(0.0296, 0.136, 0.813, 0.128))
# Shannon: median = 4.013933
# Simpson: median = 0.1248887
# Pielou: median = 0.7317871
# Chao1: median = 242

#********************* Top Bacterial Species Analysis **********************
#***************************************************************************
top_taxa <- getTopTaxa(tse, top = 6)
relative_abund_table <- decostand(t(df_taxon2), method = "total")
p <- 0.000000001
df_metadata %<>% mutate(L.crispatus=log10(df_taxon2["Faecalibacterium_prausnitzii",] + p),
                        L.iners=log10(df_taxon2["Eubacterium_rectale",] + p),
                        Gardnerella.vaginalis=log10(df_taxon2["Clostridia_bacterium",] + p),
                        L.jensenii=log10(df_taxon2["Ruminococcus_bromii",] + p),
                        L.gasseri=log10(df_taxon2["Prevotella_copri_clade_A",] + p),
                        Peptoniphilus.lacrimalis=log10(df_taxon2["Blautia_wexlerae", ] + p))

grid.arrange(ggplot(df_metadata, aes(x = PCOS, y=L.crispatus, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -8, label = "P = 0.3820", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'A) F.prausnitzii') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=L.iners, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -8, label = "P = 0.6191", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'B) E.rectale') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=Gardnerella.vaginalis, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -8, label = "P = 0.3393", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'C) C.bacterium') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=L.jensenii, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -7.5, label = "P = 0.7095", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'D) R.bromii') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=L.gasseri, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = 0, label = "P = 0.1090", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'E) P.copri_clade_A') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) ,
             ggplot(df_metadata, aes(x = PCOS, y=Peptoniphilus.lacrimalis, col = PCOS)) +
               geom_boxplot(show.legend = F,linewidth = 1.2) + theme_bw() + theme_classic() +
               geom_label(x = 1.5, y = -7.5, label = "P = 0.3496", size = 5,
                          fill = "white", color = "black", label.padding = unit(0.5, "lines"),
                          label.r = unit(0.15, "lines")) +
               labs(y = 'log10(Abundance)', x = 'PCOS', title = 'F) B.wexlerae') +
               theme(legend.position="bottom", legend.justification="center",
                     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                     legend.text = element_text(size = 14, face = "bold"),
                     axis.text = element_text(size = 12, face = "bold"),
                     axis.title = element_text(size = 16, face = "bold"),
                     legend.title = element_text(size = 16, face = "bold"),
                     text = element_text(size = 16, face = "bold")) , ncol = 3)

w1 <- wilcox.test(L.crispatus ~ PCOS, df_metadata) # p-value = 0.382
w2 <- wilcox.test(L.iners ~ PCOS, df_metadata) # p-value = 0.6191
w3 <- wilcox.test(Gardnerella.vaginalis ~ PCOS, df_metadata) # p-value = 0.3393
w4 <- wilcox.test(L.jensenii ~ PCOS, df_metadata) # p-value = 0.7095
w5 <- wilcox.test(L.gasseri ~ PCOS, df_metadata) # p-value = 0.109
w6 <- wilcox.test(Peptoniphilus.lacrimalis ~ PCOS, df_metadata) # p-value = 0.3496
#***********************************************************************
#*        PCOS Vs Chao1 - Pielou - Shannon Visualization
#*        
# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x){
                m <- mean(x)
                ymin <- m-sd(x)
                ymax <- m+sd(x)
                return(c(y = m, ymin = ymin, ymax = ymax))}

A11 <- df_metadata %>% ggplot(aes(x = PCOS, y = chao12$chao1, col = PCOS)) +
      geom_violin(linewidth = 1.2) + theme_bw()  + theme_classic() +
      labs(y = 'Chao1', x = 'PCOS', title = 'Chao1') +
      stat_summary(fun.data = data_summary) +
      geom_text(aes(label = 'P = 0.417', x = 1.5, y = 355),
                col = 'black', size = 7) +
      theme(legend.position="bottom", legend.justification="center",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.text = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 32, face = "bold")) 
      

A22 <- df_metadata %>% ggplot(aes(x = PCOS, y = shannon2$diversity_shannon, col = PCOS)) +
      geom_violin(linewidth = 1.2) + theme_bw() + theme_classic() +
      labs(y = 'Shannon', x = 'PCOS', title = 'Shannon') +
      stat_summary(fun.data = data_summary) +
      geom_text(aes(label = 'P = 0.762', x = 1.5, y = 4.56),
                col = 'black', size = 7) +
      theme(legend.position="bottom", legend.justification="center",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.text = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 32, face = "bold")) 

A33 <- df_metadata %>% ggplot(aes(x = PCOS, y = pielou2$evenness_pielou, col = PCOS)) +
      geom_violin(linewidth = 1.2) + theme_bw()  + theme_classic() +
      labs(y = 'Pielou', x = 'PCOS', title = 'Pielou') +
      stat_summary(fun.data = data_summary) +
      geom_text(aes(label = 'P = 0.565', x = 1.5, y = 0.92),
                col = 'black', size = 7) +
      theme(legend.position="bottom", legend.justification="center",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            legend.text = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 32, face = "bold")) 

Alpha_fecal <- gridExtra::grid.arrange(A11, A33, A22, ncol = 3)

#************Beta Diversity****************************************
#******************************************************************
# Using Euclidean distance 
set.seed(2023)
#* Others # centered log-ratio (CLR) --- Considered
df_phyloseq_clr2 <- microbiome::transform(df_phyloseq2, "clr")
taxa_euclidean <- phyloseq::distance(df_phyloseq_clr2, method = "euclidean")
df_meta2 <- data.frame(sample_data(df_phyloseq2))

# Homogeneity of dispersion test
beta <- betadisper(taxa_euclidean, df_meta2$PCOS)
xtable::xtable(anova(beta), digits = 4)
tabb <- permutest(beta) # p = 0.082
anova(beta)

# Adonis test for pcos_score
set.seed(2023)
perm_pcos <- adonis2(taxa_euclidean ~ PCOS, data = df_meta2) 
perm_pcos # p = 0.459
set.seed(2023)
perm_gestW <- adonis2(taxa_euclidean ~ Gest_week, data = df_meta2) 
perm_gestW # p = 0.08
set.seed(2023)
perm_age <- adonis2(taxa_bray ~ Age_groups_2, data = df_meta2) 
perm_age # p = 0.002 (significant)
set.seed(2023)
perm_fert <- adonis2(taxa_euclidean ~ Subfertility_, data = df_meta2) 
perm_fert # p =  0.03 (significant)
set.seed(2023)
perm_BMI_c <- adonis2(taxa_euclidean ~ factor(BMI_category), data = df_meta2) 
perm_BMI_c # p = 0.001 (significant)
set.seed(2023)
perm_parity <- adonis2(taxa_euclidean ~ parity, data = df_meta2) 
perm_parity # p = 0.003 (significant)

Beta_F <- data.frame(Parameters = c('PCOS', 'Gestation weeks','Age groups', 'Fertility Treatment', 'BMI categories', 'Parity'),
                     F_value1 = c(0.9924, 1.2056, 1.3484, 1.3468, 2.3646, 1.6327),
                     R2 = c(0.00135, 0.00164, 0.0073, 0.00183, 0.00957, 0.00221),
                     P_value1 = c(0.459, 0.08, 0.002, 0.03, 0.001, 0.003))

set.seed(2023)
perm_pcos <- adonis2(t(assay(tse,"rclr")) ~ PCOS,
                     by = "margin", #each term (here only
                     # ↪ 'Group') analyzed individually
                     data = colData(tse),
                     method = "euclidean",
                     permutations = 999)
perm_gestW <- adonis2(t(assay(tse,"rclr")) ~ Gest_week,
                      by = "margin", 
                      data = colData(tse),
                      method = "euclidean",
                      permutations = 999)
perm_age <- adonis2(t(assay(tse,"rclr")) ~ factor(Age_groups_2),
                    by = "margin",
                    data = colData(tse),
                    method = "euclidean",
                    permutations = 999)
perm_fert <- adonis2(t(assay(tse,"rclr")) ~ factor(Subfertility_),
                     by = "margin", 
                     data = colData(tse),
                     method = "euclidean",
                     permutations = 999)
perm_BMI_c <- adonis2(t(assay(tse,"rclr")) ~ factor(BMI_category),
                      by = "margin", 
                      data = colData(tse),
                      method = "euclidean",
                      permutations = 999) 
perm_parity <- adonis2(t(assay(tse,"clr")) ~ parity,
                       by = "margin",                 
                       data = colData(tse),
                       method = "euclidean",
                       permutations = 999) 


full_data_fecal <- as.data.frame(rbind(perm_pcos, perm_gestW, perm_age, 
                                       perm_fert, perm_BMI_c, perm_parity))
xtable::xtable(full_data_fecal, digits = 4)

#******************************************************************
#******************************************************************
#*          Data Exploration -- Beta Diversity
#******************************************************************
library(ComplexHeatmap)
library(microbiomeutilities)
library(viridis)
library(RColorBrewer)

# Using MDS Ordination plots
df_phyloseq_ra <- microbiome::transform(df_phyloseq,  'compositional')
OTU2  <- otu_table(df_taxon2[,-676], taxa_are_rows = TRUE)
meta2 <- sample_data(df_metadata)
df_phyloseq1 <- merge_phyloseq(OTU2, meta2)
ps_ord2 <- ordinate(df_phyloseq, method = "MDS", distance = "bray")
PCoA <- plot_ordination(df_phyloseq, ps_ord2, type = "samples", color = "PCOS")

ps_ordCCA <- ordinate(df_phyloseq_ra, method = "CCA", distance = "bray")
CCA <- plot_ordination(df_phyloseq, ps_ordCCA, type = "samples", color = "PCOS")

#PCA via phyloseq
ord_clr2 <- phyloseq::ordinate(df_phyloseq_clr2, "RDA")
clr11 <- ord_clr2$CA$eig[1] / sum(ord_clr2$CA$eig)
clr22 <- ord_clr2$CA$eig[2] / sum(ord_clr2$CA$eig)

P1 <- phyloseq::plot_ordination(df_phyloseq2, ord_clr2,
                          type = "samples", color="PCOS") + 
          theme_bw() +
          geom_point(size = 2) +
          coord_fixed(clr22 / clr11) +
          stat_ellipse(aes(group = PCOS), linetype = 2) +
          geom_hline(yintercept = 0, linetype = 'dashed') +
          geom_vline(xintercept = 0, linetype = 'dashed') +
          # geom_text(aes(label = paste(' P = 0.459'), 
          #               x = 3, y = 4.0),col = 'black') +
          theme(legend.position="bottom", legend.justification="center",
                plot.title = element_text(hjust = 0.5, face = "bold")) +
          theme_classic()

# Using MDS Ordination plots using tse data
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "NMDS")
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "PCA")
tse <- runMDS(tse, FUN = vegan::vegdist, method = "bray",
               name = "MDS", exprs_values = "counts")

library(patchwork)
plots <- lapply(c("MDS", "NMDS", "PCA"), plotReducedDim,
                object = tse, colour_by = "PCOS")
PCoA + labs(title = 'A) MDS', x = 'MDS1 [10.2%]', y = 'MDS2 [7.2%]') + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_bw() + theme_classic() +
  theme(legend.position="bottom", legend.justification="center",
        plot.title = element_text(hjust = 0.5, face = "bold"))  + 
  plots[[2]] + labs(title = 'B) NMDS') +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position="bottom", legend.justification="center",
        plot.title = element_text(hjust = 0.5, face = "bold"))  + 
  P1 + plot_layout(guides = "keep") + labs(title = 'C) PCA') +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position="bottom", legend.justification="center",
        plot.title = element_text(hjust = 0.5, face = "bold"))

#*****************************************************
#*          Clustered Heatmap --- Using phyloseq data
#*****************************************************                  
plot_taxa_heatmap(df_phyloseq2,
                  subset.top = 20,
                  VariableA = "PCOS",
                  heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                  transformation = "clr",
                  main = 'Clustered Heatmap',
                  cutree_rows = 6) +
  theme(legend.position = "bottom")

#*******************************************************************
#*           MDS1 & MDS2
# MDS1 and MDS2 for vaginal microbiome
df_taxon2 <- df_taxon
beta_dist <- vegdist(t(df_taxon2), method = "bray")
mds <- metaMDS(beta_dist)
mds_data <- as.data.frame(mds$points)

mds_data$Studienummer <- df_metadata$Studienummer
# df_meta <- df_meta %>% select(-X)
mds_data <- dplyr::left_join(mds_data, df_metadata)

# plot our two new dimensions and color them by Study-id
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = pcos_score)) +
  geom_point()

#*******************************************************************
#*           RCM ------------ RCM
library(RCM)
# phy <- prune_taxa(df_phyloseq)
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
mat <- assay(tse1, "clr_z") #'Sd8cdb3f' '676'
dim(mat)
mat <- mat[,-676]

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
taxa_clusters <- cutree(tree = taxa_hclust, k = 5)
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
sample_clusters <- factor(cutree(tree = sample_hclust, k = 5))
# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)
# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE]
# Order data based on
tse1 <- tse1[, rownames(sample_data)]
# Add sample type data
sample_data$PCOS <- colData(tse1)$PCOS 
sample_data$BMI <- colData(tse1)$BMI_category
sample_data$`Age groups` <- colData(tse1)$Age_groups_2
sample_data$`Age groups` <- as.character(sample_data$`Age groups`)
sample_data$`Age groups`[sample_data$`Age groups`==1]<-'≤24'
sample_data$`Age groups`[sample_data$`Age groups`==2]<-'25-29'
sample_data$`Age groups`[sample_data$`Age groups`==3]<-'30-34'
sample_data$`Age groups`[sample_data$`Age groups`==4]<-'35-39'
sample_data$`Age groups`[sample_data$`Age groups`==5]<-'40+'
sample_data$`Age groups` <- as.factor(sample_data$`Age groups`)
# Determines the scaling of colorss
# Scale colors 
breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))),
              length.out = ifelse(max(abs(mat))>5,2*ceiling(max(abs(mat))), 10))
colors <- colorRampPalette(c("darkblue", "blue", "white", "red",
                             "darkred"))(length(breaks)-1)
pheatmap(mat, annotation_row = taxa_clusters,
         annotation_col = sample_data,
         breaks = breaks, color = colors,
         cutree_rows = 7, fontsize = 13.5,
         fontsize_row = 15)

#************Differential Abundance********************************
#******************************************************************
#*                   ANCOM-BC Differential Abundance 
PCOS.abc <- ancombc2(data = tse, assay_name = "counts", 
                     fix_formula  = "PCOS", 
                     p_adj_method = "BH",  prv_cut = 0,
                     lib_cut = 0, group = "PCOS", 
                     struc_zero = TRUE, neg_lb = TRUE, 
                     em_control = list(tol = 1e-5, max_iter = 20),
                     alpha = 0.05, global = TRUE)
results_ABC <- PCOS.abc$res
results_ABC2 <- data.frame(Species = rownames(results_ABC),
                           lfc = results_ABC$lfc_PCOSYes, 
                           diff.PCOS = results_ABC$diff_PCOSYes,
                           Adj.Pval.PCOS = results_ABC$p_PCOSYes,
                           row.names = rownames(results_ABC))
diff_abund_AB <- results_ABC2 %>% filter(diff.PCOS == 'TRUE') 
nrow(diff_abund_AB) # 0

# Confounders
# Subfertility_, BMI_category, parity*********************************
PCOS_others.abc  <- ancombc2(data = tse, assay_name = "counts", 
                             fix_formula  = "PCOS + Subfertility_ + BMI_category + parity", 
                             p_adj_method = "BH",  prv_cut = 0,
                             lib_cut = 0, group = "PCOS", 
                             struc_zero = TRUE, neg_lb = TRUE, 
                             em_control = list(tol = 1e-5, max_iter = 20),
                             alpha = 0.05, global = TRUE)

results_other_ABC <- PCOS_others.abc$res
diff_abund_AB2 <- results_other_ABC %>% dplyr::filter(diff_PCOSYes == 'TRUE') 
nrow(diff_abund_AB2) # 0

#BMI_category, parity*********************************
PCOS_others.abc2  <- ancombc2(data = tse, assay_name = "counts", 
                              fix_formula  = "PCOS + BMI_category + parity", 
                              p_adj_method = "BH",  prv_cut = 0,
                              lib_cut = 0, group = "PCOS", 
                              struc_zero = TRUE, neg_lb = TRUE, 
                              em_control = list(tol = 1e-5, max_iter = 20),
                              alpha = 0.05, global = TRUE)

results_other_ABC2 <- PCOS_others.abc2$res
diff_abund_AB3 <- results_other_ABC2 %>% dplyr::filter(diff_PCOSYes == 'TRUE') 
nrow(diff_abund_AB3) # 0

# parity*********************************
PCOS_others.abc3  <- ancombc2(data = tse, assay_name = "counts", 
                              fix_formula  = "PCOS + parity", 
                              p_adj_method = "BH",  prv_cut = 0,
                              lib_cut = 0, group = "PCOS", 
                              struc_zero = TRUE, neg_lb = TRUE, 
                              em_control = list(tol = 1e-5, max_iter = 20),
                              alpha = 0.05, global = TRUE)

results_other_ABC3 <- PCOS_others.abc3$res
diff_abund_AB3 <- results_other_ABC3 %>% dplyr::filter(diff_PCOSYes == 'TRUE') 
nrow(diff_abund_AB3) # 0

# age category*********************************
PCOS_others.abc4 <- ancombc2(data = tse, assay_name = "counts", 
                              fix_formula  = "PCOS + Age_groups_2", 
                              p_adj_method = "BH",  prv_cut = 0,
                              lib_cut = 0, group = "PCOS", 
                              struc_zero = TRUE, neg_lb = TRUE, 
                              em_control = list(tol = 1e-5, max_iter = 20),
                              alpha = 0.05, global = TRUE)

results_other_ABC4 <- PCOS_others.abc4$res
diff_abund_AB4 <- results_other_ABC4 %>% dplyr::filter(diff_PCOSYes == 'TRUE') 
nrow(diff_abund_AB4) # 0

#******************************************************************
#*                   Wilcoxon Differential Abundance 
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

#******************************************************************
#*             Negative Binomial Differential Abundance 
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

#******************************************************************
#*
#*               Community Detection
#*
## Community Detection in the data Neighboring 
library(bluster)
library(patchwork)
library(scater)
# Performing and storing UMAP
tse <- transformCounts(tse, method = "rclr")
tse <- runUMAP(tse, name="UMAP", exprs_values="rclr")

k <- c(2,3,6,10)
ClustAndPlot <- function(x) {
  # Creating the graph and running the short random walks algorithm
  graph_clusters <- clusterRows(t(assays(tse)$rclr), NNGraphParam(k=x))
  
  # Results of the clustering as a color for each sample
  plotUMAP(tse, colour_by = I(graph_clusters)) +
    labs(title = paste0("k = ", x))
}

# Applying the function for different k values
plots <- lapply(k, ClustAndPlot)
# Displaying plots in a grid
(plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])

#******************************************************************
#*
#*               Genus level Analysis
# view(df_taxon_g)
dim(df_taxon) # 537 738

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
view(df_genus)

# phyloseq data type
df_genus2 <- t(df_genus)
OTU2  <- otu_table(df_genus2, taxa_are_rows = TRUE)
meta2 <- sample_data(df_metadata)
df_phyloseq_genus <- merge_phyloseq(OTU2, meta2)

#*ANCOM-BC Differential Abundance 
PCOS.abc_genus <- ancombc2(df_phyloseq_genus, fix_formula  = "PCOS", 
                     p_adj_method = "BH",  prv_cut = 0,
                     lib_cut = 0, group = "PCOS", 
                     struc_zero = TRUE, neg_lb = TRUE, 
                     em_control = list(tol = 1e-5, max_iter = 20),
                     alpha = 0.05, global = TRUE)
results_ABC_genus <- PCOS.abc_genus$res
# saveRDS(results_ABC, "Data/results_ABC_genus.rds")
results_ABC_genus <- readRDS("Data/results_ABC_genus.rds")
results_ABC3 <- data.frame(Genera = rownames(results_ABC_genus),
                           lfc = results_ABC_genus$lfc_PCOSYes, 
                           se = results_ABC_genus$se_PCOSYes,
                           w = results_ABC_genus$W_PCOSYes,
                           Adj.P.value = results_ABC_genus$q_PCOSYes)

results_ANC1 <- results_ABC3 %>%  arrange(Adj.P.value) %>% head(7)
xtable::xtable(results_ANC1, digits = 5)

diff_abund_AB <- results_ABC2 %>% filter(diff.PCOS == 'TRUE') 
nrow(diff_abund_AB) # 0










