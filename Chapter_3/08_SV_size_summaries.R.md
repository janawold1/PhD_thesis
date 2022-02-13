---
title: 'Chapter 4: Short-read SV discovery in Kākāpō'
author: "Jana R. Wold"
date: "31 January 2022"
output:
  html_document:
    self_contained: false
---

Here is the outline of charting differences among structural variants (SVs) identified using 3 discovery tools: 1) Delly; 2) Manta; and 3) Smoove.  

Below we assess variation in the number of SVs discovered by type, size and chromosome for each of these three tools. We also assess the impacts of filtering parameters. 

```{r setup, include=FALSE}
library(ggplot2)
library(ggridges)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(tidyverse)

#install.packages("devtools")
#devtools::install_github("G-Thomson/Manu")
library(Manu)

pal <- get_pal("Kakapo")
get_pal("Kakapo")
print_pal(pal)

knitr::opts_chunk$set(dev = c("svg", "png"),
                      dpi = 300,
                      echo = FALSE,
                      cache = TRUE)

setwd("G:/My Drive/Data/Kakapo/Short-read SV Chapter Data/bwa/SV_summaries")

# Data from each discovery/genotyping tool
delly <-read.table("inputs/delly_summary.tsv", sep = "\t", header = TRUE)
manta <-read.table("inputs/manta_summary.tsv", sep = "\t", header = TRUE)
smoove <-read.table("inputs/smoove_summary.tsv", sep = "\t", header = TRUE)
svs <- rbind(delly, manta, smoove)

# Sorting into filtered and unfiltered data sets
unfiltered <- svs %>%
                filter(across(data, ~grepl('unfiltered' , .)))
filtered <- svs %>%
                filter(across(data, ~grepl('SVfiltered' , .)))

# Data for generation comparisons
smoove_gen <- read.table("inputs/smoove_generations.tsv", sep = "\t", header = TRUE)
delly_gen <- read.table("inputs/delly_generations.tsv", sep = "\t", header = TRUE)
manta_joint_gen <-read.table("inputs/manta_joint_generations.tsv", sep = "\t", header = TRUE)
manta_joint_gen$size <- manta_joint_gen$end - manta_joint_gen$start

manta_batch_gen <- read.table("inputs/manta_batch_generations.tsv", sep = "\t", header = TRUE)
manta_batch_gen$size <- manta_batch_gen$end - manta_batch_gen$start

smoove_indiv_size <- smoove_gen %>%
                        group_by(indiv) %>%
                        summarise(min = min(size),
                                  max = max(size),
                                  median = median(size),
                                  mean = mean(size),
                                  sd = sd(size))
smoove_indiv_size$data <- "Smoove"

delly_indiv_size <- delly_gen %>%
                        group_by(indiv) %>%
                        summarise(min = min(size),
                                  max = max(size),
                                  median = median(size),
                                  mean = mean(size),
                                  sd = sd(size))
delly_indiv_size$data <- "Delly"
indiv_size <- rbind(smoove_indiv_size, delly_indiv_size)

delly_size_relation <- delly_gen %>% 
                          select(indiv,size) %>% 
                          group_by(indiv) %>% 
                          summarise(indiv = indiv, 
                                    count = n(),
                                    median = median(size),
                                    mean = mean(size)) %>% 
                          unique()
smoove_size_relation <- smoove_gen %>% 
                          select(indiv,size) %>% 
                          group_by(indiv) %>% 
                          summarise(indiv = indiv, 
                                    count = n(),
                                    median = median(size),
                                    mean = mean(size)) %>% 
                          unique()

# Data for testing differences in SV size and type among SV discovery tools
filtered_del <- filtered [ filtered$SV_type == "Deletion" ,]
filtered_dup <- filtered [ filtered$SV_type == "Duplication" ,]
filtered_ins <- filtered [ filtered$SV_type == "Insertion" ,]
filtered_inv <- filtered [ filtered$SV_type == "Inversion" ,]
counts <- filtered %>%
            group_by(data) %>%
            count(vars = SV_type)
counts <- spread(counts, vars, n)
```

Here I explored the number of structural variants per data set before and after filtering. 

To begin we have structural variants discovered using Delly, Manta and Smoove. Two strategies were used for structural variant discovery with Manta, 1) a batched approach; and 2) joint variant calling. 

```{r type counts, echo=FALSE, fig.align='center'}
ggplot(unfiltered, aes(x=SV_type, fill=data)) +
  geom_histogram(stat = "count", position = position_dodge()) +
  labs(x = "Structural Variant Type", y = "Count", title = "Count of all Structural Variants") +
  scale_fill_manual(values = c("delly_unfiltered" = "#7D9D33", "manta_batch_unfiltered" = "#CD8862", "manta_joint_unfiltered" = "#BCA888", "smoove_unfiltered" = "#DCC949")) +
  theme_light()

ggplot(filtered, aes(x=SV_type, fill=data)) +
  geom_histogram(stat = "count", position = position_dodge()) +
  labs(x = "Structural Variant Type", y = "Count", title = "Count of Filtered Structural Variants") +
  scale_fill_manual(values = c("delly_SVfiltered" = "#7D9D33", "manta_batch_SVfiltered" = "#CD8862", "manta_joint_SVfiltered" = "#BCA888", "smoove_SVfiltered" = "#DCC949")) +
  theme_light()
```

Here, we compared the size distribution of SVs called by all three tools. Due to the experimental set up violating many assumptions of independece, comparisons of variance between the tools were not made. Instead differences are in the context of means and given the highly skewed size distribution across SV types in all data sets, the natural log of SV size was used to graph the size distribution and estimate a mean (geometric mean).

```{r size_distribution, echo=FALSE, fig.align = 'center'}
filtered$logsize <- log(filtered$size)
head(filtered)

# Not transformed
ggplot(filtered, aes(x=size, fill=data)) +
  geom_density(alpha = 0.7) +
  xlim(0,10000) +
  labs(x = "Structural Variant Size", y = "Proportion", title = "Size Distribution of Filtered Structural Variants") +
  theme_light() +
  facet_wrap(~SV_type, scales = "free_y") +
  scale_fill_manual(values = c("delly_SVfiltered" = "#7D9D33", "manta_batch_SVfiltered" = "#CD8862", "manta_joint_SVfiltered" = "#BCA888", "smoove_SVfiltered" = "#DCC949"))

ggplot(filtered, aes(x=size, fill=data)) +
  geom_histogram(alpha = 0.7, position = "dodge") +
  xlim(0,10000) +
  labs(x = "Structural Variant Size", y = "Proportion", title = "Size Distribution of Filtered Structural Variants") +
  theme_light() +
  facet_wrap(~SV_type, scales = "free_y") +
  scale_fill_manual(values = c("delly_SVfiltered" = "#7D9D33", "manta_batch_SVfiltered" = "#CD8862", "manta_joint_SVfiltered" = "#BCA888", "smoove_SVfiltered" = "#DCC949"))

# Size transformed

ggplot(filtered, aes(x=logsize, fill=data)) +
  geom_density(alpha = 0.7) +
  labs(x = "Structural Variant Size", y = "Proportion", title = "Size Distribution of Filtered Structural Variants") +
  theme_light() +
  facet_wrap(~SV_type, scales = "free_y") +
  scale_fill_manual(values = c("delly_SVfiltered" = "#7D9D33", "manta_batch_SVfiltered" = "#CD8862", "manta_joint_SVfiltered" = "#BCA888", "smoove_SVfiltered" = "#DCC949"))

ggplot(filtered, aes(x=logsize, fill=data)) +
  geom_histogram(alpha = 0.7, position = "dodge") +
  labs(x = "Structural Variant Size", y = "Proportion", title = "Size Distribution of Filtered Structural Variants") +
  theme_light() +
  facet_wrap(~SV_type, scales = "free_y") +
  scale_fill_manual(values = c("delly_SVfiltered" = "#7D9D33", "manta_batch_SVfiltered" = "#CD8862", "manta_joint_SVfiltered" = "#BCA888", "smoove_SVfiltered" = "#DCC949"))

```

And finally, we assessed whether some chromosomes were more 'rich' in SVs than others and the relative proportion of each chromosome impacted by SVs.

```{r Chromosome Comparisons, echo=FALSE, fig.align='center'}
filtered$proportion <- filtered$size / filtered$chrom_size

# The absolute count of unfiltered SVs per chromosome
ggplot(unfiltered, aes(x = chrom, fill = SV_type)) + 
  geom_histogram(stat = "count", position = "stack") +
  labs(x = "Chromosome", y = "Count", title = "Total Number of Structural Variants per Chromosome") +
  theme_light() + 
  facet_wrap(~data, scales = "free_y") +
  scale_fill_manual(values = c("Deletion" = "#7D9D33", "Duplication" = "#CD8862", "Insertion" = "#BCA888", "Inversion" = "#DCC949"))

# The absolute count of filtered SVs per chromosome
ggplot(filtered, aes(x = chrom, fill = SV_type)) + 
  geom_histogram(stat = "count", position = "stack") +
  labs(x = "Chromosome", y = "Count", title = "Number of Filtered Structural Variants per Chromosome") +
  theme_light() + 
  facet_wrap(~data, scales = "free_y") +
  scale_fill_manual(values = c("Deletion" = "#7D9D33", "Duplication" = "#CD8862", "Insertion" = "#BCA888", "Inversion" = "#DCC949"))


# The proportion of the chromosome impacted by each filtered data set
filtered %>%
  ggplot(aes(x = chrom, y = proportion, fill = SV_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Chromosome", y = "Proportion of impacted base pairs", title = "Relative proportion of chromosome impacted by structural variants") +
  theme_light() +
  facet_wrap(~data, scales = "free_y") +
  scale_fill_manual(values = c("Deletion" = "#7D9D33", "Duplication" = "#CD8862", "Insertion" = "#BCA888", "Inversion" = "#DCC949"))

```

We then compared a rough count of SVs carried by individual kākāpō by type and putative 'origin'.

```{r Genotype Comparisons, echo=FALSE, fig.align='center'}
ggplot(delly_gen, aes(x = indiv, fill=SV_type)) +
  geom_histogram(stat = 'count') +
  labs(x = "Individual", y = "Structural Variant Count", title = "Number of Delly SVs carried by individual kākāpō by population") +
  theme_light() +
  scale_fill_manual(values = c("Deletion" = "#7D9D33", "Duplication" = "#CD8862", "Insertion" = "#BCA888", "Inversion" = "#DCC949"))

ggplot(smoove_gen, aes(x = indiv, fill=SV_type)) +
  geom_histogram(stat = 'count') +
  labs(x = "Individual", y = "Structural Variant Count", title = "Number of Smoove SVs carried by individual kākāpō by population") +
  theme_light() +
  scale_fill_manual(values = c("Deletion" = "#7D9D33", "Duplication" = "#CD8862", "Insertion" = "#BCA888", "Inversion" = "#DCC949"))

ggplot(manta_batch_gen, aes(x = indiv, fill=SV_type)) +
  geom_histogram(stat = 'count') +
  ylim(0,350) +
  labs(x = "Individual", y = "Structural Variant Count", title = "Number of Manta - Batch SVs carried by individual kākāpō by population") +
  theme_light() +
  scale_fill_manual(values = c("Deletion" = "#7D9D33", "Duplication" = "#CD8862", "Insertion" = "#BCA888", "Inversion" = "#DCC949"))

ggplot(manta_joint_gen, aes(x = indiv, fill=SV_type)) +
  geom_histogram(stat = 'count') +
  ylim(0,350) +
  labs(x = "Individual", y = "Structural Variant Count", title = "Number of Manta - Joint SVs carried by individual kākāpō by population") +
  theme_light() +
  scale_fill_manual(values = c("Deletion" = "#7D9D33", "Duplication" = "#CD8862", "Insertion" = "#BCA888", "Inversion" = "#DCC949"))

```
```{r Individual_Size_Comparisons, echo=FALSE, fig.align='center'}
ggplot(delly_size_relation, aes(x = indiv, y = mean)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "Mean structural variant size", title = "Mean Delly SV size per individual kākāpō") +
  theme_light()

ggplot(smoove_size_relation, aes(x = indiv, y = mean)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "Mean structural variant size", title = "Mean Smoove SV size per individual kākāpō") +
  theme_light()

ggplot(delly_size_relation, aes(x = count, y = mean, label = indiv)) +
  geom_point() +
  labs(x = "Number of SVs Carried", y = "Mean structural variant size", title = "Relationship between number of Delly SVs carried and mean size of SVs") +
  geom_text(aes(label=indiv), hjust = 0, vjust = 0) +
  theme_light()

ggplot(smoove_size_relation, aes(x = count, y = mean, label = indiv)) +
  geom_point() +
  labs(x = "Number of SVs Carried", y = "Mean structural variant size", title = "Relationship between number of Smoove SVs carried and mean size of SVs") +
  geom_text(aes(label=indiv), hjust = 0, vjust = 0) +
  theme_light()

ggplot(delly_size_relation, aes(x = count, y = median, label = indiv)) +
  geom_point() +
  labs(x = "Number of SVs Carried", y = "Median structural variant size", title = "Relationship between number of Delly SVs carried and median size of SVs") +
  geom_text(aes(label=indiv), hjust = 0, vjust = 0) +
  theme_light()

ggplot(smoove_size_relation, aes(x = count, y = median, label = indiv)) +
  geom_point() +
  labs(x = "Number of SVs Carried", y = "Median structural variant size", title = "Relationship between number of Smoove SVs carried and median size of SVs") +
  geom_text(aes(label=indiv), hjust = 0, vjust = 0) +
  theme_light()

delly_gen %>%
  summarise(mu = mean(size),
            popvar = var(size))

ggplot(smoove_size_relation, aes(x = mean)) +
  geom_density()
```


```{r Generational_Comparisons, echo=FALSE, fig.align='center'}
delly_gen %>%
  count(indiv, gen) %>%
  ggplot(aes(x = gen, y = n, fill = gen)) +
    geom_violin() +
    labs(x = "Ancestry", y = "Total Number of Delly SV's", title = "Counts of Delly SV's by generation ") +
    theme(axis.text = element_text(size = 14)) +
    scale_fill_manual(values = c("RH0" = "#AD6B17","RH1" = "#CCB62F", "RH2" = "#D0C471", "SI0" = "#3E4331", "SI1" = "#66743B", "SI2" = "#BAC4C2")) +
    theme_light()

smoove_gen %>%
  count(indiv,gen) %>%
  ggplot(aes(x = gen, y = n, fill = gen)) +
    geom_violin() +
    labs(x = "Ancestry", y = "Total Number of Smoove SV's", title = "Counts of Smoove SV's by generation ") +
    theme(axis.text = element_text(size = 14)) +
    scale_fill_manual(values = c("RH0" = "#AD6B17","RH1" = "#CCB62F", "RH2" = "#D0C471", "SI0" = "#3E4331", "SI1" = "#66743B", "SI2" = "#BAC4C2")) +
    theme_light()

manta_batch_gen %>%
  count(indiv,gen) %>%
  ggplot(aes(x = gen, y = n, fill = gen)) +
    geom_violin() +
    ylim(0,350) +
    labs(x = "Ancestry", y = "Total Number of SV's", title = "Counts of Manta - Batch SV's by generation ") +
    theme(axis.text = element_text(size = 14)) +
    scale_fill_manual(values = c("RH0" = "#AD6B17","RH1" = "#CCB62F", "RH2" = "#D0C471", "SI0" = "#3E4331", "SI1" = "#66743B", "SI2" = "#BAC4C2")) +
    theme_light()

manta_joint_gen %>%
  count(indiv,gen) %>%
  ggplot(aes(x = gen, y = n, fill = gen)) +
    geom_violin() +
    ylim(0,350) +
    labs(x = "Ancestry", y = "Total Number of SV's", title = "Counts of Manta - Joint SV's by generation ") +
    theme(axis.text = element_text(size = 14)) +
    scale_fill_manual(values = c("RH0" = "#AD6B17","RH1" = "#CCB62F", "RH2" = "#D0C471", "SI0" = "#3E4331", "SI1" = "#66743B", "SI2" = "#BAC4C2")) +
    theme_light()
```
