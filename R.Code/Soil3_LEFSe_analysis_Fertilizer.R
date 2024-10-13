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

# row.names(ASV)=ASV$X
# ASV=ASV[,-which(names(ASV) %in% c("X"))]

meta=read.table("data/metadata.txt",header = TRUE,sep = "\t")
row.names(meta)=meta$Sample_Name
# meta=meta[,-which(names(meta) %in% c("Sample_Name"))]

#meta$group=factor(meta$group, levels = c("Topsoil", "Sbusoil"))

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
  select(-Resolution, -Guild,-Kingdom)

###fertilizer______________________________________________________________________________

meta <- meta %>%
  select(-Depth,-Sample_Name, -Group, -Depth_range) %>%
  mutate(sample=rownames(meta)) %>%
  rename(group=Fertilizer)

dataset <- microtable$new(sample_table = meta, otu_table = ASV.t, tax_table = tax.clean)
dataset

lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = "group", lefse_subgroup = NULL, 
                        #taxa_level = "Species"
                        )

head(lefse$res_diff)

# data.taxa.sig.fer <- lefse$res_diff
# 
# xlsx::write.xlsx(data.taxa.sig.fer, file = "data.taxa.sig.fer.xlsx", col.names = T, row.names = T)

p.LDA.fer <- lefse$plot_diff_bar(use_number = 1:40,width=0.15, group_order = c("Control","NK","NP","PK","NPK")) +
  theme_classic(base_size = 8) +
  scale_fill_brewer(palette = "Set2", guide=F)+
  scale_color_brewer(palette = "Set2", guide=F)+
  labs(fill="Fertilizer")+
  theme(
    plot.background = element_blank(),
    plot.margin = margin(5.5, 0, 50, 5.5),
    legend.position = c(0.1,-0.15),
    legend.key.height = unit(0.2, "cm"),
  )+
  guides(fill=guide_legend(nrow = 1))

p.LDA.fer

ggsave("figures/p.LDA.fer.tiff",width = 3.5,height = 4)


