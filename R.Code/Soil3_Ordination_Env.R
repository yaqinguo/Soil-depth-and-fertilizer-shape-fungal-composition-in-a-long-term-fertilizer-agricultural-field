library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
#import ASV table and clean and format
ASV=read.table("data/ASV_table.txt",header = TRUE,sep = "\t")
row.names(ASV)=ASV$X
ASV=ASV[,-which(names(ASV) %in% c("X"))]

rowSums(ASV)
ASV.Norm <- decostand(ASV,"hellinger")
rowSums(ASV.Norm)

meta=read.table("data/metadata.txt",header = TRUE,sep = "\t")
row.names(meta)=meta$Sample_Name
meta=meta[,-which(names(meta) %in% c("Sample_Name"))]
meta$Fertilizer=factor(meta$Fertilizer,levels=c("Control","NK","NP","PK","NPK"))

Env=read.table("data/Environmental_data.txt",header = TRUE,sep = "\t")
row.names(Env)=Env$Sample_Name
Env=Env[,-which(names(Env) %in% c("Sample_Name"))]

# Env %>%
#   group_by(Depth) %>%
#   summarise(texture_clay.mean=mean(texture_clay),
#             TC.mean=mean(TC),
#             CN.mean=mean(CN),
#             BD.mean=mean(bulk_density))


Env=Env %>% mutate(CN=TC/TN) %>%
  select(bulk_density, Moisture, TN, CN, TC, pH, EC,texture_clay, texture_silt, texture_sand)

#normalize numeric variables 
env_log <- Env %>%
  mutate(across(where(is.numeric),log1p))


dist <- vegdist(ASV.Norm,method = "bray")
set.seed(1)
#here i remove TN
db.rda.1 <- dbrda(dist ~ bulk_density + Moisture + CN + TC 
                  + pH + EC + texture_clay + texture_silt + texture_sand, data=env_log)

set.seed(629)
anova(db.rda.1, by = "term", permutations = 999)
anova(db.rda.1, by = "axis", permutations = 999)

db.rda.res <- summary(db.rda.1)
db.rda.res

df1 <- data.frame(db.rda.res$sites[,1:2])
df1 <- cbind(df1, meta)
df1$Depth <- factor(df1$Depth)
df1$Group <- factor(df1$Group)
df2 <- data.frame(db.rda.res$biplot[,1:2])
df2$Factor <- row.names(df2)
df2 <- df2 %>% filter(Factor !="texture_silt" & Factor !="texture_sand" & Factor !="Moisture")
df2 <- df2 %>% mutate(Factor = recode(Factor,
                                      'texture_clay' ="Clay",
                                      'TC'="TOC",
                                      'bulk_density'="BD",
                                      'CN' = "C/N"))

 db.rda.plot <- ggplot() +
   theme_bw(base_size = 5)+
   geom_vline(xintercept = 0, color='grey',size=0.2, linetype=2)+
   geom_hline(yintercept = 0, color='grey', size=0.2,linetype=2)+
   labs(color="Depth (cm)", x = "dbRDA1 (7.23%)", y = "dbRDA2 (3.07%)")+
   geom_point(data = df1, aes(x=dbRDA1, y=dbRDA2,color=Depth,fill=Depth), size=1)+
   scale_color_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"),labels=c("0-30","30-40","40-50","50-70","70-100"))+
   scale_fill_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"),labels=c("0-30","30-40","40-50","50-70","70-100"))+
   theme(panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.box = "vertical",
          legend.spacing.y = unit(0,"cm"),
          legend.key.height = unit(0.3,"cm"),
          legend.key.width = unit(0,"cm"),
          plot.background = element_blank())+
   geom_segment(data = df2, color="black", size=0.2, 
                aes(x=0, y=0, xend=dbRDA1, yend=dbRDA2),
                arrow = arrow(length = unit(0.01, "npc")))+
   geom_text(data = df2,aes(x=dbRDA1, y=dbRDA2, label=Factor,
                            hjust=0.7, vjust=0.6*(1-sign(dbRDA2))), size=1.2)+
   guides(fill=F)
   
db.rda.plot

ggsave("figures/db.rda.plot.tiff",width = 2.5,height = 2)

















