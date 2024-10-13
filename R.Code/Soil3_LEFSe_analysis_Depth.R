library(microeco)
library(magrittr)
library(patchwork)
library(tidyverse)

#import ASV table and clean and format
ASV=read.table("data/ASV_table.txt",header = TRUE,sep = "\t")
ASV.t=as.data.frame(t(ASV))

colnames(ASV.t) <- ASV.t[1,]
ASV.t <- ASV.t[-1,]
ASV.t[] <- lapply(ASV.t,as.numeric)

meta=read.table("data/metadata.txt",header = TRUE,sep = "\t")
row.names(meta)=meta$Sample_Name

#import tax table and clean and format
tax=read.table("data/Taxonomy_table.txt",header = TRUE,sep = "\t")
row.names(tax)=tax$X
tax.clean=tax[,-which(names(tax) %in% c("X"))]
tax.clean$Species=gsub("Undetermined Undetermined","Undetermined",tax.clean$Species)
tax.clean <- tax.clean %>%
  mutate(Kingdom=paste0("k_",Kingdom)) %>%
  mutate(Phylum=paste0("p_",Phylum)) %>%
  mutate(Class=paste0("c_",Class)) %>%
  mutate(Order=paste0("o_",Order)) %>%
  mutate(Family=paste0("f_",Family)) %>%
  mutate(Genus=paste0("g_",Genus)) %>%
  mutate(Species=paste0("s_", Species)) %>%
  select(-Kingdom,-Resolution,-Guild)

#each depth_______________________________________________________________________________________

meta <- meta %>%
  select(-Depth, -Group, -Fertilizer,-Sample_Name) %>%
  mutate(sample=rownames(meta)) %>%
  rename(group=Depth_range)

dataset <- microtable$new(sample_table = meta, otu_table = ASV.t, tax_table = tax.clean)
dataset 

lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = "group",
                        lefse_subgroup = NULL)

#results: 6 Phylum;19 class;35 order;Family 73; Genus 101

head(lefse$res_diff)

data.taxa.sig <- lefse$res_diff

xlsx::write.xlsx(data.taxa.sig, file = "data.taxa.sig.eachdepth.xlsx", col.names = T, row.names = T)

p.LDA.eachdepth <- lefse$plot_diff_bar(use_number = 1:30,width=0.2,
                    group_order = c("0-30","40-50","50-70","70-100"))+
  theme_classic(base_size = 8)+
  scale_fill_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  scale_color_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"), guide=F)+
  labs(fill="Depth (cm)")+
  theme(
    plot.background = element_blank(),
    plot.margin = margin(5.5, 0, 50, 5.5),
    legend.position = c(0.2,-0.15),
        # legend.justification = "left",
        legend.key.height = unit(0.2, "cm"),
        # legend.text = element_text(size = 6)
        )+
  guides(fill=guide_legend(nrow = 1))

p.LDA.eachdepth
ggsave("figures/p.LDA.eachdepth.tiff",width = 3,height = 4)