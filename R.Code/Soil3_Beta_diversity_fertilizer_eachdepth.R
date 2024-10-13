library(phyloseq)
library(tidyverse)
library(ggplot2)
library(vegan)
library(reshape2)
#import ASV table and clean and format
ASV=read.table("data/ASV_table.txt",header = TRUE,sep = "\t")
ASV.t=as.data.frame(t(ASV))
colnames(ASV.t) <- ASV.t[1,]
ASV.t <- ASV.t[-1,]
ASV.t[] <- lapply(ASV.t,as.numeric) #this is because later to convert to phyloseq object with error. it was noted becasue of non-numeric value, but in fact all are numeric

#import tax table and clean and format
tax=read.table("data/Taxonomy_table.txt",header = TRUE,sep = "\t")
row.names(tax)=tax$X
tax.clean=tax[,-which(names(tax) %in% c("X"))]
tax.clean$Species=gsub("Undetermined Undetermined","Undetermined",tax.clean$Species)

#import metadata
meta=read.table("data/metadata.txt",header = TRUE,sep = "\t")
row.names(meta)=meta$Sample_Name
meta=meta[,-which(names(meta) %in% c("Sample_Name"))]

ASV.UF=otu_table(as.matrix(ASV.t),taxa_are_rows = T)
tax.UF=tax_table(as.matrix(tax.clean))
meta.UF=sample_data(as.data.frame(meta))

soil3 <- phyloseq(ASV.UF,tax.UF,meta.UF)
soil3 <- phyloseq::prune_taxa(taxa_sums(soil3)>0,soil3)

# Extract sample data
sample_data_df <- data.frame(sample_data(soil3))

# Split sample data by Depth_range
sample_data_list <- split(sample_data_df, sample_data_df$Depth_range)

# Function to subset phyloseq object and transform to relative abundances
transform_to_relative <- function(physeq, sample_data_subset) {
  # Subset phyloseq object
  physeq_subset <- prune_samples(rownames(sample_data_subset), physeq)
  
  # Transform to relative abundances
  physeq_subset_rel <- transform_sample_counts(physeq_subset, function(counts) counts / sum(counts))
  
  return(physeq_subset_rel)
}

# Apply the function to each subset
physeq_rel_abundance_list <- lapply(sample_data_list, function(subset) transform_to_relative(soil3, subset))

#here 0-30______________________________________________________________________________________________________
physeq.30 <- physeq_rel_abundance_list$`0-30`
bray.30 <- phyloseq::distance(otu_table(physeq.30),method = "bray")
bray.30.matrix <- as.matrix(bray.30)
head(bray.30.matrix)[,1:6]
dim(bray.30.matrix)

bray.30.melt <- melt(bray.30.matrix) 
meta.30 <- sample_data(physeq.30)
meta.30 <- data.frame(meta.30)
meta.30$Sample <- rownames(meta.30)
dim(meta.30)
bray.30.df <- merge(bray.30.melt, meta.30[,c("Sample","Fertilizer","Group")], 
                    by.x="Var1",by.y="Sample",all.x=T)
head(bray.30.df)
bray.30.df <- bray.30.df[bray.30.df$value !=0, ]

ggplot(bray.30.df, aes(x=Fertilizer, y=value, color=Fertilizer))+
  geom_boxplot(alpha=0.6)+
  stat_boxplot(aes(ymin= ..lower.., ymax= ..upper..),outlier.shape = NA)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..))+
  geom_jitter() +
  labs(x = "", y = "Bray-Curtis Distance") +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        #legend.position = c(1,0.15),
        #legend.justification = c(1.5,0)
        )
  
# res.30 <- aov(value ~ Fertilizer, data = bray.30.df)
# summary(res.30)

kruskal.test(value ~ Fertilizer, data = bray.30.df)

#here 30-40______________________________________________________________________________________________________
physeq.40 <- physeq_rel_abundance_list$`30-40`
bray.40 <- phyloseq::distance(otu_table(physeq.40),method = "bray")
bray.40.matrix <- as.matrix(bray.40)
head(bray.40.matrix)[,1:6]
dim(bray.40.matrix)

bray.40.melt <- melt(bray.40.matrix) 
meta.40 <- sample_data(physeq.40)
meta.40 <- data.frame(meta.40)
meta.40$Sample <- rownames(meta.40)
dim(meta.40)
bray.40.df <- merge(bray.40.melt, meta.40[,c("Sample","Fertilizer", "Group")], 
                    by.x="Var1",by.y="Sample",all.x=T)
head(bray.40.df)
bray.40.df <- bray.40.df[bray.40.df$value !=0, ]

ggplot(bray.40.df, aes(x=Fertilizer, y=value, color=Fertilizer))+
  geom_boxplot(alpha=0.6)+
  stat_boxplot(aes(ymin= ..lower.., ymax= ..upper..),outlier.shape = NA)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..))+
  geom_jitter() +
  labs(x = "", y = "Bray-Curtis Distance") +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        #legend.position = c(1,0.15),
        #legend.justification = c(1.5,0)
  )

res.40 <- aov(value ~ Fertilizer, data = bray.40.df)
summary(res.40)

kruskal.test(value ~ Fertilizer, data = bray.40.df)

pairwise.wilcox.test(bray.40.df$value, bray.40.df$Fertilizer, p.adjust.method= "bonferroni")

#here 40-50______________________________________________________________________________________________________
physeq.50 <- physeq_rel_abundance_list$`40-50`
bray.50 <- phyloseq::distance(otu_table(physeq.50),method = "bray")
bray.50.matrix <- as.matrix(bray.50)
head(bray.50.matrix)[,1:6]
dim(bray.50.matrix)

bray.50.melt <- melt(bray.50.matrix) 
meta.50 <- sample_data(physeq.50)
meta.50 <- data.frame(meta.50)
meta.50$Sample <- rownames(meta.50)
dim(meta.50)
bray.50.df <- merge(bray.50.melt, meta.50[,c("Sample","Fertilizer", "Group")], 
                    by.x="Var1",by.y="Sample",all.x=T)
head(bray.50.df)
bray.50.df <- bray.50.df[bray.50.df$value !=0, ]

ggplot(bray.50.df, aes(x=Fertilizer, y=value, color=Fertilizer))+
  geom_boxplot(alpha=0.6)+
  stat_boxplot(aes(ymin= ..lower.., ymax= ..upper..),outlier.shape = NA)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..))+
  geom_jitter() +
  labs(x = "", y = "Bray-Curtis Distance") +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        #legend.position = c(1,0.15),
        #legend.justification = c(1.5,0)
  )

res.50 <- aov(value ~ Fertilizer, data = bray.50.df)
summary(res.50)

kruskal.test(value ~ Fertilizer, data = bray.50.df)

pairwise.wilcox.test(bray.50.df$value, bray.50.df$Fertilizer, p.adjust.method= "bonferroni")


#here 50-70______________________________________________________________________________________________________
physeq.70 <- physeq_rel_abundance_list$`50-70`
bray.70 <- phyloseq::distance(otu_table(physeq.70),method = "bray")
bray.70.matrix <- as.matrix(bray.70)
head(bray.70.matrix)[,1:6]
dim(bray.70.matrix)

bray.70.melt <- melt(bray.70.matrix) 
meta.70 <- sample_data(physeq.70)
meta.70 <- data.frame(meta.70)
meta.70$Sample <- rownames(meta.70)
dim(meta.70)
bray.70.df <- merge(bray.70.melt, meta.70[,c("Sample","Fertilizer", "Group")], 
                    by.x="Var1",by.y="Sample",all.x=T)
head(bray.70.df)
bray.70.df <- bray.70.df[bray.70.df$value !=0, ]

ggplot(bray.70.df, aes(x=Fertilizer, y=value, color=Fertilizer))+
  geom_boxplot(alpha=0.6)+
  stat_boxplot(aes(ymin= ..lower.., ymax= ..upper..),outlier.shape = NA)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..))+
  geom_jitter() +
  labs(x = "", y = "Bray-Curtis Distance") +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        #legend.position = c(1,0.15),
        #legend.justification = c(1.5,0)
  )

# res.70 <- aov(value ~ Fertilizer, data = bray.70.df)
# summary(res.70)

kruskal.test(value ~ Fertilizer, data = bray.70.df)

pairwise.wilcox.test(bray.70.df$value, bray.70.df$Fertilizer, p.adjust.method= "bonferroni")

#here 70-100______________________________________________________________________________________________________
physeq.100 <- physeq_rel_abundance_list$`70-100`
bray.100 <- phyloseq::distance(otu_table(physeq.100),method = "bray")
bray.100.matrix <- as.matrix(bray.100)
head(bray.100.matrix)[,1:6]
dim(bray.100.matrix)

bray.100.melt <- melt(bray.100.matrix) 
meta.100 <- sample_data(physeq.100)
meta.100 <- data.frame(meta.100)


meta.100$Sample <- rownames(meta.100)
dim(meta.100)
bray.100.df <- merge(bray.100.melt, meta.100[,c("Sample","Fertilizer", "Group")], 
                    by.x="Var1",by.y="Sample",all.x=T)
head(bray.100.df)
bray.100.df <- bray.100.df[bray.100.df$value !=0, ]

ggplot(bray.100.df, aes(x=Fertilizer, y=value, color=Fertilizer))+
  geom_boxplot(alpha=0.6)+
  stat_boxplot(aes(ymin= ..lower.., ymax= ..upper..),outlier.shape = NA)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..))+
  geom_jitter() +
  labs(x = "", y = "Bray-Curtis Distance") +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        #legend.position = c(1,0.15),
        #legend.justification = c(1.5,0)
  )
res.100 <- aov(value ~ Fertilizer, data = bray.100.df)
summary(res.100)

kruskal.test(value ~ Fertilizer, data = bray.100.df)

pairwise.wilcox.test(bray.100.df$value, bray.100.df$Fertilizer, p.adjust.method= "bonferroni")

#overall_______________________________________________________________________________________________________________
data <-rbind(bray.30.df, bray.40.df,bray.50.df,bray.70.df, bray.100.df)
data <- data %>%
  separate(Var1, into=c("Number","Depth"), sep="-",remove=F)
data$Depth <- recode(data$Depth, "30"="0-30", "40"="30-40","50"="40-50","70"="50-70", "100"="70-100")
data$Depth <- factor(data$Depth, levels = c("0-30","30-40","40-50","50-70","70-100"))
data$Group <- factor(data$Group, levels = c("Topsoil","Subsoil"))
data$Fertilizer <- factor(data$Fertilizer, levels = c("Control","NK","NP","PK","NPK"))
p_fer_eachdepth <- ggplot(data, aes(x=Fertilizer, y=value, color=Fertilizer))+
  facet_wrap(~Depth) +
  geom_boxplot(linetype="dashed", width=0.1)+
  stat_boxplot(aes(ymin= ..lower.., ymax= ..upper..),outlier.shape = NA,width=0.5)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.5)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.5)+
  #geom_jitter(alpha=0.3) +
  labs(x = "", y = "Bray-Curtis Distance") +
  theme_bw(base_size = 8) +
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(1,0.2),
        legend.key.size = unit(10,"pt"),
        legend.justification = c(1.5,0)
  )
p_fer_eachdepth
ggsave("figures/p_fer_eachdepth.tiff",width = 5,height = 4)

p_fer_ts <- ggplot(data, aes(x=Fertilizer, y=value, color=Fertilizer))+
  facet_wrap(~Group) +
  geom_boxplot(linetype="dashed")+
  stat_boxplot(aes(ymin= ..lower.., ymax= ..upper..),outlier.shape = NA)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..))+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..))+
  #geom_jitter(alpha=0.6) +
  labs(x = "", y = "Bray-Curtis Distance") +
  theme_bw(base_size = 14) +
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        #legend.position = c(1,0.15),
        #legend.justification = c(1.5,0)
  )
p_fer_ts 
ggsave("figures/p_fer_ts.tiff",width = 6,height = 5)

data.sub <- data %>%
  filter(Group=="Subsoil")
res.sub <- aov(value ~ Fertilizer, data = data.sub)
summary(res.sub)

kruskal.test(value ~ Fertilizer, data = data.sub)

pairwise.wilcox.test(data.sub$value, data.sub$Fertilizer,p.adjust.method= "bonferroni")

data.top <- data %>%
  filter(Group=="Topsoil")
res.top <- aov(value ~ Fertilizer, data = data.top)
summary(res.top)

kruskal.test(value ~ Fertilizer, data = data.top)


