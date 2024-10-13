library(phyloseq)
library(tidyverse)
library(ggplot2)
library(vegan)
library(scales)
library(ggpubr)
library(RColorBrewer)
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

# Combine the transformed subsets back into a single phyloseq object
soil3_rel_abundance <- do.call(merge_phyloseq, physeq_rel_abundance_list)

# Aggregate at the Phylum level
soil3_phylum <- tax_glom(soil3_rel_abundance, taxrank = "Phylum")

# Melt the data for ggplot2
soil3_melt <- psmelt(soil3_phylum)

# Calculate the relative abundance per Phylum within each Depth_range
phyla_summary <- soil3_melt %>%
  group_by(Depth_range, Phylum) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  mutate(RelativeAbundance = TotalAbundance / sum(TotalAbundance) * 100)

# Identify Phyla with relative abundance < 1% and replace them with "Other"
phyla_to_collapse <- phyla_summary %>%
  filter(RelativeAbundance < 1) %>%
  pull(Phylum) %>%
  unique()

# Update the Phylum column in the melted data
soil3_melt <- soil3_melt %>%
  mutate(Phylum = ifelse(Phylum %in% phyla_to_collapse, "Other", Phylum))

# soil3_melt_mean <- soil3_melt %>%
#   group_by(Depth_range,Phylum)%>%
#   summarise(mean=mean(Abundance))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

soil3_melt2 <- data_summary(soil3_melt, varname = "Abundance",
                            groupnames = c("Depth_range","Phylum"))

soil3_melt2$Depth_range=as.factor(soil3_melt2$Depth_range)

p.line.phylum.depth <- ggplot(soil3_melt2,aes(x = Depth_range, y =Abundance,group=Phylum,color=Phylum))+
  theme_classic(base_size = 8)+
  #coord_flip()+
  geom_line(size=0.5,aes(color=Phylum))+
  geom_errorbar(aes(ymin=Abundance-sd, ymax=Abundance+sd), width=0.2)+
  geom_point(aes(color=Phylum),size=1)+
  scale_color_brewer(palette = "Accent")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Depth (cm)", y = "Relative Abundance (%)") +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8),
        legend.margin = margin(0,0,0,0),
        legend.key.size = unit(0.3,"cm"),
        legend.spacing = unit(-1.5,"cm"),
        #aspect.ratio = 1,
        #legend.box.margin = margin(10,10,100,10),
        #plot.margin = margin(0,0,0,0),
        #axis.text.x.bottom = element_text(angle = -90,size=6),
        #legend.position = c(0.2,0.8)
        )
p.line.phylum.depth
ggsave("figures/p.line.phylum.depth.tiff",width = 3.2,height = 2.5)

#here Guild_______________________________________________________________________________________________________
# Aggregate at the Guild level
#remove unassigned at the beginning 

tax.guild <- rownames(tax.clean[tax.clean$Guild != "Unassigned",]) 
soil3.guild <- prune_taxa(tax.guild, soil3)

# Extract sample data
sample_data_df <- data.frame(sample_data(soil3.guild))

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
physeq_rel_abundance_list <- lapply(sample_data_list, function(subset) transform_to_relative(soil3.guild, subset))

# Combine the transformed subsets back into a single phyloseq object
soil3_rel_abundance <- do.call(merge_phyloseq, physeq_rel_abundance_list)


soil3_Guild <- tax_glom(soil3_rel_abundance, taxrank = "Guild")

# Melt the data for ggplot2
soil3_melt <- psmelt(soil3_Guild)

soil3_melt <- soil3_melt %>%
  select(Abundance,Depth_range,Guild) %>%
  mutate(Depth_range=factor(Depth_range, levels = c("0-30","30-40","40-50","50-70","70-100"))) %>%
  mutate(Guild=factor(Guild, levels = c("Saprotroph","Pathotroph","Arbuscular Mycorrhizal")))

soil3_melt_mean <- soil3_melt %>%
  group_by(Depth_range,Guild) %>%
  summarise(mean=mean(Abundance))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

soil3_melt_guild <- data_summary(soil3_melt, varname = "Abundance", groupnames =c("Depth_range","Guild"))

soil3_melt_guild$Depth_range=as.factor(soil3_melt_guild$Depth_range)


# Create the plot with y-axis as percentage
p.bar.Guild.Depth <- ggplot(soil3_melt_guild, aes(x = Depth_range, y = Abundance, fill=Guild)) +
  #facet_wrap(~Guild)+
  #geom_errorbar(aes(ymin=Abundance-sd, ymax=Abundance+sd), width=0.2)+
  geom_bar(stat = "identity", width = 0.2,show.legend = T, position = "fill") + # Use 'position = "fill"' to stack bars to 100%
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + # Ensure y-axis shows percentage
  labs(x = "Depth (cm)", y = "Relative Abundance (%)") +
  theme_bw(base_size = 8) +
  scale_fill_brewer(palette = "Accent")+
  #scale_fill_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
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
p.bar.Guild.Depth

ggsave("figures/p.bar.Guild.Depth.tiff",width = 3.2,height = 2)


#here i did anova analysis for each Phylum at different Depth 
soil3_melt1 <- soil3_melt %>%
  filter(Guild=="Arbuscular Mycorrhizal")

kruskal.test(Abundance ~ Depth_range,data=soil3_melt1)



soil3_melt2 <- soil3_melt %>%
  filter(Guild=="Pathotroph")

kruskal.test(Abundance ~ Depth_range,data=soil3_melt2)

soil3_melt3 <- soil3_melt %>%
  filter(Guild=="Saprotroph")

kruskal.test(Abundance ~ Depth_range,data=soil3_melt3)
