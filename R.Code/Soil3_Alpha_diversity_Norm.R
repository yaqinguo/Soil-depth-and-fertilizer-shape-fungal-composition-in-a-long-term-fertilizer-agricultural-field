library(phyloseq)
library(tidyverse)
library(ggplot2)
library(vegan)
library(car) 
#import ASV table and clean and format
ASV=read.table("data/ASV_table.txt",header = TRUE,sep = "\t")
row.names(ASV)=ASV$X
ASV=ASV[,-which(names(ASV) %in% c("X"))]

rowSums(ASV)
ASV.Norm <- decostand(ASV,method = "total")
rowSums(ASV.Norm)

meta=read.table("data/metadata.txt",header = TRUE,sep = "\t")
row.names(meta)=meta$Sample_Name
meta=meta[,-which(names(meta) %in% c("Sample_Name"))]
meta$Fertilizer=factor(meta$Fertilizer,levels=c("NK","NP","PK","NPK","Control"))
meta$richness=specnumber(ASV.Norm)
meta$Shannon=diversity(ASV.Norm,index = "shannon")
meta$Invsimpson=diversity(ASV.Norm,index = "invsimpson")

table(meta$Depth_range,meta$Fertilizer) #to check balance
meta$Depth_range=as.factor(meta$Depth_range)

# statistic test for richness
richness.aov <- aov(richness ~ Fertilizer+Depth_range,data=meta)
summary(richness.aov)

#Depth_range P=6.2e-14 ***  Fertilizer P=0.77  
#only depth affect the richness but not fertilizer

#to check homegen
#this is visiluation to check it 
plot(richness.aov,1)

#this is to use statistical method to check 
leveneTest(richness ~ Fertilizer*Depth_range,data=meta ) #P=0.686
#check residuals normal or not 
#this is visiluation to check it 
plot(richness.aov,2)
#this is to use statistical method to check 
richness_residuals <- residuals(object=richness.aov)
shapiro.test(x=richness_residuals)

P_observd<-ggplot(meta, aes(x=Depth_range,y=richness,color=Depth_range))+
  theme_classic(base_size = 10)+
  stat_boxplot(geom = "errorbar",width=0.1,show.legend = F)+
  geom_boxplot(aes(color=Depth_range),position = position_dodge(0.9),width=0.2,show.legend = F,outlier.shape = 2)+
  geom_point(cex=0.3,shape=21,color="black", fill="white",aes(fill=Depth_range),show.legend = F)+
  scale_color_manual(values = c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  scale_fill_manual(values = c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  scale_y_continuous(limits = c(0,480))+
  labs(x="Depth (cm)",y="Observed ASV")+
  theme(legend.position = c(0.88,0.88),
        legend.background = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.background = element_blank())
P_observd
ggsave("figures/P_observd_Norm.tiff",width = 3,height = 2.5)


P_observd_fer<-ggplot(meta, aes(x=Depth_range,y=richness,fill=Fertilizer))+
  theme_bw(base_size = 10)+
  stat_boxplot(position = position_dodge(0.9),geom = "errorbar",width=0.3,show.legend = F)+
  geom_boxplot(position = position_dodge(0.9),width=0.3)+
  scale_fill_brewer(palette="Set2",labels=c("Control","NP","NK","PK","NPK"))+
  labs(x="Depth (cm)",y="Observed ASV")+
  theme(legend.position = c(0.88,0.88),
        legend.background = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        #legend.title = element_blank(),
        panel.grid = element_blank(),plot.background = element_blank())
P_observd_fer
ggsave("figures/P_observd_fer.tiff",width = 5,height = 4)


# statistic test for shannon
shannon.aov <- aov(Shannon ~ Fertilizer+Depth_range,data=meta)
summary(shannon.aov )

#Depth_range 1.23e-07 ***  Fertilizer P=0.796 
#only depth affect the shannon diversity but not fertilizer

#to check homegen
#this is visiluation to check it 
plot(shannon.aov,1)
#this is to use statistical method to check 
leveneTest(Shannon ~ Fertilizer*Depth_range,data=meta ) #P=0.2088
#check residuals normal or not 
#this is visiluation to check it 
plot(shannon.aov,2)
#this is to use statistical method to check 
shannon_residuals <- residuals(object=shannon.aov)
shapiro.test(x=shannon_residuals)

library(multcomp)
summary(glht(richness.aov, linfct = mcp(Depth_range = "Tukey")))


#here plot shannon 
P_shannon<-ggplot(meta, aes(x=Depth_range,y=Shannon,color=Depth_range))+
  theme_classic(base_size = 10)+
  stat_boxplot(geom = "errorbar",width=0.1,show.legend = F)+
  geom_boxplot(aes(color=Depth_range),position = position_dodge(0.9),width=0.2,show.legend = F,outlier.shape = 2)+
  geom_point(cex=0.3,shape=21,aes(fill=Depth_range),color="black",show.legend = F)+
  scale_color_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  scale_fill_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  scale_y_continuous(limits=c(0,6))+
  labs(x="Depth (cm)",y="Shannon Diversity")+
  theme(legend.position = c(0.88,0.88),
        legend.background = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.background = element_blank())
P_shannon
ggsave("figures/P_shannon_Norm.tiff",width = 3,height = 2.5)

P_shannon_fer<-ggplot(meta, aes(x=Depth_range,y=Shannon,fill=Fertilizer))+
  theme_bw(base_size = 10)+
  stat_boxplot(position = position_dodge(0.9),geom = "errorbar",width=0.3,show.legend = F)+
  geom_boxplot(position = position_dodge(0.9),width=0.3)+
  scale_fill_brewer(palette="Set2",labels=c("Control","NP","NK","PK","NPK"))+
  labs(x="Depth (cm)",y="Shannon Diversity")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.background = element_blank())
P_shannon_fer
ggsave("figures/P_shannon_fer.tiff",width = 5,height = 4)

# statistic test for shannon at each depth
data1=meta %>% filter(Depth_range=="0-30")
shannon.aov.fer <- aov(Shannon ~ Fertilizer,data=data1)
summary(shannon.aov.fer)
data2=meta %>% filter(Depth_range=="30-40")
shannon.aov.fer <- aov(Shannon ~ Fertilizer,data=data2)
summary(shannon.aov.fer)
data3=meta %>% filter(Depth_range=="40-50")
shannon.aov.fer <- aov(Shannon ~ Fertilizer,data=data3)
summary(shannon.aov.fer)
data4=meta %>% filter(Depth_range=="50-70")
shannon.aov.fer <- aov(Shannon ~ Fertilizer,data=data4)
summary(shannon.aov.fer)
data5=meta %>% filter(Depth_range=="70-100")
shannon.aov.fer <- aov(Shannon ~ Fertilizer,data=data5)
summary(shannon.aov.fer)


# statistic test for InvSimpson
Invsimpson.aov <- aov(Invsimpson ~ Fertilizer+Depth_range,data=meta)
summary(Invsimpson.aov)

#Depth_range 1.26e-05 ***  Fertilizer P=0.832 
#only depth affect the Invsimpson diversity but not fertilizer

#to check homegen
#this is visiluation to check it 
plot(Invsimpson.aov,1)
#this is to use statistical method to check 
leveneTest(Invsimpson ~ Fertilizer*Depth_range,data=meta ) #P=0.2088
#check residuals normal or not 
#this is visiluation to check it 
plot(Invsimpson.aov,2)
#this is to use statistical method to check 
Invsimpson_residuals <- residuals(object=Invsimpson.aov)
shapiro.test(x=Invsimpson_residuals)


P_InvSimpson<-ggplot(meta, aes(x=Depth_range,y=Invsimpson,color=Depth_range))+
  theme_classic(base_size = 10)+
  stat_boxplot(geom = "errorbar",width=0.1,show.legend = F)+
  geom_boxplot(aes(color=Depth_range),position = position_dodge(0.9),width=0.2,show.legend = F, outlier.shape = 2)+
  geom_point(cex=0.3,shape=21,aes(fill=Depth_range),color="black",show.legend = F)+
  scale_color_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  scale_fill_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  scale_y_continuous(limits = c(0,40))+
  #scale_y_continuous(limits=c(0,5))+
  labs(x="Depth (cm)",y="InvSimpson Diversity")+
  theme(legend.position = c(0.88,0.88),
        legend.background = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.background = element_blank())
P_InvSimpson
ggsave("figures/P_InvSimpson_Norm.tiff",width = 3,height = 2.5)


P_InvSimpson_fer<-ggplot(meta, aes(x=Depth_range,y=Invsimpson,fill=Fertilizer))+
  theme_bw(base_size = 10)+
  stat_boxplot(position = position_dodge(0.9),geom = "errorbar",width=0.3,show.legend = F)+
  geom_boxplot(position = position_dodge(0.9),width=0.3)+
  scale_fill_brewer(palette="Set2",labels=c("Control","NP","NK","PK","NPK"))+
  labs(x="Depth (cm)",y="InvSimpson Diversity")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key.size = unit(0.3,"cm"),
        legend.title = element_blank(),
        panel.grid = element_blank(),plot.background = element_blank())
P_InvSimpson_fer
ggsave("figures/P_InvSimpson_fer.tiff",width = 5,height = 4)

# statistic test for InvSimpson at each depth
data1=meta %>% filter(Depth_range=="0-30")
Inv.aov.fer <- aov(Invsimpson ~ Fertilizer,data=data1)
summary(Inv.aov.fer)
data2=meta %>% filter(Depth_range=="30-40")
Inv.aov.fer <- aov(Invsimpson ~ Fertilizer,data=data2)
summary(Inv.aov.fer)
data3=meta %>% filter(Depth_range=="40-50")
Inv.aov.fer <- aov(Invsimpson ~ Fertilizer,data=data3)
summary(Inv.aov.fer)
data4=meta %>% filter(Depth_range=="50-70")
Inv.aov.fer <- aov(Invsimpson ~ Fertilizer,data=data4)
summary(Inv.aov.fer)
data5=meta %>% filter(Depth_range=="70-100")
Inv.aov.fer <- aov(Invsimpson ~ Fertilizer,data=data5)
summary(Inv.aov.fer)
