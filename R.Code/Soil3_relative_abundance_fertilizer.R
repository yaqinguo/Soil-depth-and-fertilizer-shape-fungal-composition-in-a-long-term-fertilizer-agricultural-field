library(phyloseq)
library(tidyverse)
library(ggplot2)
library(vegan)
library(scales)
library(ggpubr)
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
sample_data_list <- split(sample_data_df, sample_data_df$Fertilizer)

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

# Combine the transformed subsets back into a single phyloseq object
soil3_rel_abundance <- do.call(merge_phyloseq, physeq_rel_abundance_list)

# Aggregate at the Phylum level
soil3_phylum <- tax_glom(soil3_rel_abundance, taxrank = "Phylum")

# Melt the data for ggplot2
soil3_melt <- psmelt(soil3_phylum)

# Calculate the relative abundance per Phylum within each Depth_range
phyla_summary <- soil3_melt %>%
  group_by(Fertilizer, Phylum) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  mutate(RelativeAbundance = TotalAbundance / sum(TotalAbundance) * 100)

# Identify Phyla with relative abundance < 5% and replace them with "Other"
phyla_to_collapse <- phyla_summary %>%
  filter(RelativeAbundance < 1) %>%
  pull(Phylum) %>%
  unique()

# Update the Phylum column in the melted data
soil3_melt <- soil3_melt %>%
  mutate(Phylum = ifelse(Phylum %in% phyla_to_collapse, "Other", Phylum))
soil3_melt$Fertilizer=factor(soil3_melt$Fertilizer,levels = c("Control","NK","NP","PK","NPK"))
# Create the plot with y-axis as percentage
p.bar.Phylum.Fer <- ggplot(soil3_melt, aes(x = Fertilizer, y = Abundance,fill=Phylum)) +
  geom_bar(stat = "identity", position = "fill") + # Use 'position = "fill"' to stack bars to 100%
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) + # Ensure y-axis shows percentage
  labs(x = "Depth Range", y = "Relative Abundance (%)") +
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Set1")+
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.margin = margin(0,0,0,0),
        legend.key.size = unit(0.5,"cm"),
        legend.spacing = unit(-2,"cm"),
        #legend.box.margin = margin(10,10,100,10),
        #plot.margin = margin(0,0,0,0),
        #axis.text.x.bottom = element_text(angle = -90,size=6),
        #legend.position = "bottom"
  )
p.bar.Phylum.Fer
ggsave("figures/p.bar.Phylum.Fer.tiff",width = 5,height = 4)

soil3_melt_mean <- soil3_melt %>%
  group_by(Fertilizer,Phylum)%>%
  summarise(mean=mean(Abundance))

p.line.phylum.fer <- ggplot(soil3_melt_mean,aes(x = Fertilizer, y = mean,group=Phylum))+
  theme_bw(base_size = 14)+
  geom_line(size=1,aes(color=Phylum))+
  geom_point(aes(color=Phylum),size=3)+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Depth Range", y = "Relative Abundance (%)") +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.margin = margin(0,0,0,0),
        legend.key.size = unit(0.5,"cm"),
        legend.spacing = unit(-2,"cm"),
        #legend.box.margin = margin(10,10,100,10),
        #plot.margin = margin(0,0,0,0),
        #axis.text.x.bottom = element_text(angle = -90,size=6),
        #legend.position = "bottom"
        )
p.line.phylum.fer
ggsave("figures/p.line.phylum.fer.tiff",width = 5,height = 4)

#here i did anova analysis for each Phylum at different fertilizer
soil3_melt1 <- soil3_melt %>%
  filter(Phylum=="Ascomycota")

kruskal.test(Abundance ~ Fertilizer, data=soil3_melt1)

soil3_melt2 <- soil3_melt %>%
  filter(Phylum=="Basidiomycota")

kruskal.test(Abundance ~ Fertilizer, data=soil3_melt2)

soil3_melt3 <- soil3_melt %>%
  filter(Phylum=="Mortierellomycota")

kruskal.test(Abundance ~ Fertilizer, data=soil3_melt3)

soil3_melt4 <- soil3_melt %>%
  filter(Phylum=="Glomeromycota")

kruskal.test(Abundance ~ Fertilizer, data = soil3_melt4)

pairwise.wilcox.test(soil3_melt4$Abundance, soil3_melt4$Fertilizer,p.adjust.methods="bonferroni")

soil3_melt4.mean <- soil3_melt4 %>%
  group_by(Fertilizer) %>%
  summarise(mean=mean(Abundance))

  
p.fer.Glo<- ggboxplot(soil3_melt4,x="Fertilizer", y="Abundance",color="Fertilizer",
                      width = 0.5)+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..,color=Fertilizer),width=0.2,show.legend = F)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..,color=Fertilizer),width=0.2)+
  stat_compare_means(label.x = 1,label.y = 0.5,size=2)+
  stat_compare_means(ref.group="Control", label="p.signif", label.y = 0.2, size=2)+
  theme()+
  theme_classic(base_size = 8)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  labs(y="Relative Abundance")+
  scale_color_brewer(palette = "Set2",guide=F)

p.fer.Glo
ggsave("figures/p.fer.Glo.tiff",width = 3.2,height = 2.5)

??stat_compare_means
soil3_melt5 <- soil3_melt %>%
  filter(Phylum=="Other")
aov5 <- aov(Abundance ~ Fertilizer,data=soil3_melt5)
summary(aov5)

#here Guild_______________________________________________________________________________________________________
# Aggregate at the Guild level
#remove unassigned at the beginning 

tax.guild <- rownames(tax.clean[tax.clean$Guild != "Unassigned",]) 
soil3.guild <- prune_taxa(tax.guild, soil3)

# Extract sample data
sample_data_df <- data.frame(sample_data(soil3.guild))

# Split sample data by Depth_range
sample_data_list <- split(sample_data_df, sample_data_df$Fertilizer)

# Function to subset phyloseq object and transform to relative abundances
transform_to_relative <- function(physeq, sample_data_subset) {
  # Subset phyloseq object
  physeq_subset <- prune_samples(rownames(sample_data_subset), physeq)
  
  # Transform to relative abundances
  physeq_subset_rel <- transform_sample_counts(physeq_subset, function(counts) counts / sum(counts))
  
  return(physeq_subset_rel)
}

# Apply the function to each subset
physeq_rel_abundance_list <- lapply(sample_data_list, function(subset) transform_to_relative(soil3.guild, subset))

# Combine the transformed subsets back into a single phyloseq object
soil3_rel_abundance <- do.call(merge_phyloseq, physeq_rel_abundance_list)


soil3_Guild <- tax_glom(soil3_rel_abundance, taxrank = "Guild")

# Melt the data for ggplot2
soil3_melt <- psmelt(soil3_Guild)

soil3_melt <- soil3_melt %>%
  select(Abundance,Fertilizer,Guild) %>%
  # filter(Guild != "Unassigned") %>%
  #mutate(Abundance=Abundance*100) %>%
  mutate(Fertilizer=factor(Fertilizer, levels = c("Control","NK","NP","PK","NPK"))) %>%
  mutate(Guild=factor(Guild, levels = c("Saprotroph","Pathotroph","Arbuscular Mycorrhizal")))

soil3_melt.mean <- soil3_melt %>%
  group_by(Fertilizer,Guild) %>%
  summarise(mean=mean(Abundance))

# Create the plot with y-axis as percentage
p.bar.Guild.Fer <- ggplot(soil3_melt, aes(x = Fertilizer, y = Abundance,fill=Guild)) +
  #facet_wrap(~Guild)+
  #geom_area(alpha=0.6,size=1,color="black")+
  geom_bar(stat = "identity", width=0.2, position = "fill") + # Use 'position = "fill"' to stack bars to 100%
  labs(x = "Fertilizer", y = "Relative Abundance (%)") +
  theme_bw(base_size = 8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Accent")+
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        #legend.position = "none",
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8),
        legend.margin = margin(0,0,0,0),
        legend.key.size = unit(0.5,"cm"),
        legend.spacing = unit(-2,"cm"),
        strip.text = element_text(size=6),
        #legend.box.margin = margin(10,10,100,10),
        #plot.margin = margin(0,0,0,0),
        axis.text.x.bottom = element_text(angle = 0,size=6),
        #legend.position = "bottom"
  )
p.bar.Guild.Fer
ggsave("figures/p.bar.Guild.Fer.tiff",width = 3.2,height = 2)


soil3_melt1 <- soil3_melt %>%
  filter(Guild=="Arbuscular Mycorrhizal")

kruskal.test(Abundance ~ Fertilizer, data = soil3_melt1)

pairwise.wilcox.test(soil3_melt1$Abundance, soil3_melt1$Fertilizer)

soil3_melt2 <- soil3_melt %>%
  filter(Guild=="Pathotroph")

kruskal.test(Abundance ~ Fertilizer, data = soil3_melt2)
pairwise.wilcox.test(soil3_melt2$Abundance, soil3_melt2$Fertilizer)

soil3_melt3 <- soil3_melt %>%
  filter(Guild=="Saprotroph")

kruskal.test(Abundance ~ Fertilizer, data = soil3_melt3)
pairwise.wilcox.test(soil3_melt3$Abundance, soil3_melt3$Fertilizer)






