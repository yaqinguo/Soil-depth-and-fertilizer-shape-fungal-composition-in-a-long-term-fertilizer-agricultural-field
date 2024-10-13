library(phyloseq)
library(tidyverse)
library(ggplot2)
library(vegan)
library(ape)
#import ASV table and clean and format
ASV=read.table("data/ASV_table.txt",header = TRUE,sep = "\t")
row.names(ASV)=ASV$X
ASV=ASV[,-which(names(ASV) %in% c("X"))]

ASV.Norm <- decostand(ASV,method = "total")
rowSums(ASV.Norm)

meta=read.table("data/metadata.txt",header = TRUE,sep = "\t")

#bray distance
dist_matrix_bray <- vegdist(ASV.Norm,method = "bray")

cap_model_bray <- capscale(dist_matrix_bray ~ Fertilizer*Depth_range, data = meta)
summary(cap_model_bray)
cap_scores_bray <- as.data.frame(scores(cap_model_bray,display="sites"))
cap_scores_bray <- cbind(cap_scores_bray,meta)

cap_plot_bray <- ggplot(cap_scores_bray,aes(x=CAP1,y=CAP2,shape=Fertilizer))+
  geom_point(size=1.5,aes(fill=Depth_range,color=Depth_range))+
  theme_bw(base_size = 5)+
  geom_hline(yintercept =0,linetype="dashed",color="grey")+
  geom_vline(xintercept =0,linetype="dashed",color="grey")+
  scale_shape_manual(values = c(21,22,23,24,25),labels=c("Control","NK","NP","PK","NPK"))+
  scale_color_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"),labels=c("0-30","30-40","40-50","50-70","70-100"))+
  scale_fill_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"),labels=c("0-30","30-40","40-50","50-70","70-100"))+
  labs(x="CAP1 (12.71%)", y="CAP2 (7.66%)", color="Depth (cm)")+
  theme(panel.grid = element_blank(),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical",
        legend.spacing.y = unit(0,"cm"),
        legend.key.height = unit(0.3,"cm"),
        legend.key.width = unit(0,"cm"),
        plot.background = element_blank()
  )+
  guides(fill=F)
cap_plot_bray

ggsave("figures/cap_plot_bray.tiff",width = 2.8,height = 2)


adonis(dist_matrix_bray ~ Fertilizer*Depth_range,data = meta,permutations = 10000)

#only depth__________________________________________________________________________________________________
cap_depth_bray <- capscale(dist_matrix_bray ~ Depth_range + Condition(Fertilizer), data = meta)
summary(cap_depth_bray)

depth_scores_bray <- as.data.frame(scores(cap_depth_bray,display="sites"))
depth_scores_bray <- cbind(depth_scores_bray,meta)

depth_bray_mean <- depth_scores_bray %>%
  group_by(Depth_range) %>%
  summarise(mean_CAP1=mean(CAP1),
            mean_CAP2=mean(CAP2))

depth_scores_bray <- depth_scores_bray %>%
  left_join(depth_bray_mean, by="Depth_range")

cap_depth_bray <- ggplot(depth_scores_bray,aes(x=CAP1,y=CAP2,color=Depth_range))+
  geom_point(size=3)+
  theme_bw()+
  geom_segment(aes(xend=mean_CAP1, yend=mean_CAP2),
               arrow = arrow(length = unit(0.2,"cm")), alpha=0.5, show.legend = F)+
  geom_hline(yintercept =0,linetype="dashed",color="grey")+
  geom_vline(xintercept =0,linetype="dashed",color="grey")+
  scale_color_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"))+
  labs(x="CAP1 (13.11%)", y="CAP2 (3.52%)", color="Depth (cm)")+
  theme(panel.grid = element_blank(),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical",
        #legend.spacing.y = unit(0,"cm"),
        #legend.key.height = unit(-1,"cm"),
        #legend.key.width = unit(-1.2,"cm"),
        plot.background = element_blank()
  )
cap_depth_bray 
ggsave("figures/cap_depth_bray.tiff",width = 5,height = 4)

#1 fit the model to remove the effect of fertilizer
cap_fertilizer_out <- capscale(dist_matrix_bray ~ Fertilizer, data=meta)
#2 extract th residuals from the model
resuiduals_fertilizer <- residuals(cap_fertilizer_out)
#3 perform adonis on the residuals with respect to depth alone
adonis(resuiduals_fertilizer ~ Depth_range, data=meta, permutations=1000)


#Jaccard distance
dist_matrix_jac <- vegdist(ASV.Norm,method = "jaccard")

cap_model_jac <- capscale(dist_matrix_jac ~ Fertilizer*Depth_range, data = meta)
summary(cap_model_jac)
cap_scores_jac <- as.data.frame(scores(cap_model_jac,display="sites"))
cap_scores_jac <- cbind(cap_scores_jac,meta)

cap_plot_jac <- ggplot(cap_scores_jac,aes(x=CAP1,y=CAP2,shape=Fertilizer,color=Depth_range))+
  geom_point(size=3,aes(fill=Depth_range))+
  theme_bw()+
  geom_hline(yintercept =0,linetype="dashed",color="grey")+
  geom_vline(xintercept =0,linetype="dashed",color="grey")+
  scale_shape_manual(values = c(21,22,23,24,25),labels=c("Control","NK","NP","PK","NPK"))+
  scale_color_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"),labels=c("0-30","30-40","40-50","50-70","70-100"))+
  scale_fill_manual(values=c("#543005", "#8c510a","#bf812d","#dfc27d","#f6e8c3"),labels=c("0-30","30-40","40-50","50-70","70-100"))+
  labs(x="CAP1 (8.22%)", y="CAP2 (5.33%)", color="Depth (cm)")+
  theme(panel.grid = element_blank(),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical",
        #legend.spacing.y = unit(0,"cm"),
        #legend.key.height = unit(-1,"cm"),
        #legend.key.width = unit(-1.2,"cm"),
        plot.background = element_blank()
  )+guides(fill=F)
cap_plot_jac

ggsave("figures/cap_plot_jac.tiff",width = 5,height = 4)


adonis(dist_matrix_jac ~ Fertilizer*Depth_range,data = meta,permutations = 10000)

cap_depth_jac <- capscale(dist_matrix_jac ~ Depth_range + Condition(Fertilizer), data = meta)
summary(cap_depth_jac)

depth_scores_jac <- as.data.frame(scores(cap_depth_jac,display="sites"))
depth_scores_jac <- cbind(depth_scores_jac,meta)


cap_depth_jac <- ggplot(depth_scores_jac,aes(x=CAP1,y=CAP2,color=Depth_range))+
  geom_point(size=3)+
  theme_bw()+
  geom_hline(yintercept =0,linetype="dashed",color="grey")+
  geom_vline(xintercept =0,linetype="dashed",color="grey")+
  scale_color_brewer(palette="Set1")+
  labs(x="CAP1 (8.34%)", y="CAP (2.45%)", color="Depth (cm)")+
  theme(panel.grid = element_blank(),
        #legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box = "vertical",
        #legend.spacing.y = unit(0,"cm"),
        #legend.key.height = unit(-1,"cm"),
        #legend.key.width = unit(-1.2,"cm"),
        plot.background = element_blank()
  )
cap_depth_jac 
ggsave("figures/cap_depth_jac.tiff",width = 5,height = 4)

#1 fit the model to remove the effect of fertilizer
cap_fertilizer_out <- capscale(dist_matrix_jac ~ Fertilizer, data=meta)
#2 extract th residuals from the model
resuiduals_fertilizer <- residuals(cap_fertilizer_out)
#3 perform adonis on the residuals with respect to depth alone
adonis(resuiduals_fertilizer ~ Depth_range, data=meta, permutations=1000)
