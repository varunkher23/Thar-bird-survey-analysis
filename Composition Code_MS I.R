##### Loading packages ####
library(tidyverse)
library(vegan)
library(reshape2)
library(ggfortify)
library(iNEXT)
library(VennDiagram)

#### Importing data ####
library(readxl)
a <- read_excel("A:/Thar/Final data/Data/Data.xlsx")
a$Type=factor(a$Type, levels = c("Protected Grasslands", "Rangelands","Non-irrigated Croplands","Irrigated Croplands"))
a=a%>%
  filter(Type!="Irrigated Croplands")
a$corrected_number=(a$Number/a$`Detection probability`)*(1000/a$`Truncation distance`)*(1/(a$Transect_length*3/1000))

comm_data=dcast(a,Type+Transect_Name+Season~Species, value.var = "corrected_number",sum)

#Massaging Data 
comm_data$Hab=comm_data$Type
comm_data = comm_data %>%
  mutate(Hab = as.character(Hab))%>%
  mutate(Hab = ifelse(Hab=="Protected Grasslands"|Hab== "Rangelands", "Grasslands",Hab))%>%
  mutate(Hab = ifelse(Hab=="Irrigated Croplands"|Hab== "Non-irrigated Croplands","Agriculture", Hab))
comm_data$Type=factor(comm_data$Type, levels = c("Protected Grasslands",
                                                 "Rangelands",
                                                 "Non-irrigated Croplands"))
comm_data$Season=factor(comm_data$Season, levels=c("Winter","Summer"))
group <- factor(comm_data$Type, levels = c("Protected Grasslands",
                                           "Rangelands",
                                           "Non-irrigated Croplands"))
comm_data_winter=comm_data%>%
  filter(Season=="Winter")
comm_data_summer=comm_data%>%
  filter(Season=="Summer")


####NMDS####

nmds_winter <- metaMDS(comm_data_winter[,c(-1:-3,-ncol(comm_data_winter))], k=2, trymax = 1000)
nmds_winter
nmds_summer <- metaMDS(comm_data_summer[,c(-1:-3,-ncol(comm_data_summer))], k=2, trymax = 1000)
nmds_summer

##### perMANOVA ####

library(RVAideMemoire)
library(ecodist)
set.seed(5)

adonis_winter=adonis2(comm_data_winter[,c(-1,-2,-3,-ncol(comm_data_winter))]~comm_data_winter$Type,
               na.action="na.fail", distance = "bray")
write.table(tibble(adonis_winter),"clipboard")
adonis_summer=adonis2(comm_data_summer[,c(-1,-2,-3,-ncol(comm_data_summer))]~comm_data_summer$Type,
               na.action="na.fail", distance = "bray")
write.table(tibble(adonis_summer),"clipboard")

bray_winter=bcdist(comm_data_winter[c(-1,-2,-3,-ncol(comm_data_winter))])
pairwise.adonis_winter=data.frame(pairwise.perm.manova(bray_winter, comm_data_winter$Type, 
                                                p.method = "bonferroni")[["p.value"]],
                           nperm=9999)
write.table(tibble(pairwise.adonis_winter),"clipboard")
TukeyHSD(betadisper(bray_winter, group=comm_data_winter$Type))
betadisper(bray_winter, group=comm_data_winter$Type)

bray_summer=bcdist(comm_data_summer[c(-1,-2,-3,-ncol(comm_data_summer))])
pairwise.adonis_summer=data.frame(pairwise.perm.manova(bray_summer, comm_data_summer$Type, 
                                                       p.method = "bonferroni")[["p.value"]],
                                  nperm=9999)
pairwise.adonis_summer
TukeyHSD(betadisper(bray, group=comm_data_summer$Type))
betadisper(bray, group=comm_data_summer$Type)

#### ANOSIM ####

anosim1_winter=anosim(comm_data_winter[comm_data_winter$Type!="Rangelands",c(-1,-2,-3,-ncol(comm_data_winter))],comm_data_winter$Type[comm_data_winter$Type!="Rangelands"], distance = "bray")
anosim2_winter=anosim(comm_data_winter[comm_data_winter$Type!="Non-irrigated Croplands",c(-1,-2,-3,-ncol(comm_data_winter))],comm_data_winter$Type[comm_data_winter$Type!="Non-irrigated Croplands"], distance = "bray")
anosim3_winter=anosim(comm_data_winter[comm_data_winter$Type!="Protected Grasslands",c(-1,-2,-3,-ncol(comm_data_winter))],comm_data_winter$Type[comm_data_winter$Type!="Protected Grasslands"], distance = "bray")

anosim_winter=data.frame(Type=c("Rangelands","Non-irrigated Croplands"),
                         Protected_Grasslands=c(anosim2_winter[["statistic"]],anosim1_winter[["statistic"]]),
                         Rangelands=c("NA",anosim3_winter[["statistic"]]))

anosim1_summer=anosim(comm_data_summer[comm_data_summer$Type!="Rangelands",c(-1,-2,-3,-ncol(comm_data_summer))],comm_data_summer$Type[comm_data_summer$Type!="Rangelands"], distance = "bray")
anosim2_summer=anosim(comm_data_summer[comm_data_summer$Type!="Non-irrigated Croplands",c(-1,-2,-3,-ncol(comm_data_summer))],comm_data_summer$Type[comm_data_summer$Type!="Non-irrigated Croplands"], distance = "bray")
anosim3_summer=anosim(comm_data_summer[comm_data_summer$Type!="Protected Grasslands",c(-1,-2,-3,-ncol(comm_data_summer))],comm_data_summer$Type[comm_data_summer$Type!="Protected Grasslands"], distance = "bray")

anosim_summer=data.frame(Type=c("Rangelands","Non-irrigated Croplands"),
                         Protected_Grasslands=c(anosim2_summer[["statistic"]],anosim1_summer[["statistic"]]),
                         Rangelands=c("NA",anosim3_summer[["statistic"]]))

write.table(tibble(anosim_winter),"clipboard")

rm(anosim1_summer,anosim2_summer,anosim3_summer, anosim2_winter,anosim1_winter, anosim3_winter)

##### Visualisation #####

#NMDS 

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

scores_winter <- data.frame(scores(nmds_winter))
scores_winter$sites <- rownames(nmds_winter)
scores_winter <- data.frame(scores_winter, group=comm_data_winter$Type)
NMDS.mean_winter=aggregate(scores_winter[,1:2],list(group=comm_data_winter$Type),mean)
plot(nmds_winter,type = "t",display = "sites")
ordiellipse(nmds_winter, comm_data_winter$Type, display = "sites", kind = "se", conf = 0.95)
ord_winter<-ordiellipse(nmds_winter, comm_data_winter$Type, display = "sites", kind = "se", conf = 0.95)

df_ell3_winter <- data.frame()
for(g in levels(scores_winter$group)){
  df_ell3_winter <- rbind(df_ell3_winter, cbind(as.data.frame(with(scores_winter[scores_winter$group==g,],
                                                                   veganCovEllipse(ord_winter[[g]]$cov,
                                                                                   ord_winter[[g]]$center,
                                                                                   ord_winter[[g]]$scale)))
                                                ,group=g))
}

scores_summer <- data.frame(scores(nmds_summer))
scores_summer$sites <- rownames(nmds_summer)
scores_summer <- data.frame(scores_summer, group=comm_data_summer$Type)
NMDS.mean_summer=aggregate(scores_summer[,1:2],list(group=comm_data_summer$Type),mean)
plot(nmds_summer,type="t",display = "sites")
ordiellipse(nmds_summer, comm_data_summer$Type, display = "sites", kind = "se", conf = 0.95)
ord_summer<-ordiellipse(nmds_summer, comm_data_summer$Type, display = "sites", kind = "se", conf = 0.95)

df_ell3_summer <- data.frame()
for(g in levels(scores_summer$group)){
  df_ell3_summer <- rbind(df_ell3_summer, cbind(as.data.frame(with(scores_summer[scores_summer$group==g,],
                                                                   veganCovEllipse(ord_summer[[g]]$cov,
                                                                                   ord_summer[[g]]$center,
                                                                                   ord_summer[[g]]$scale)))
                                                ,group=g))
}

scores_winter$Habitat=as.factor(comm_data_winter$Hab)
scores_summer$Habitat=as.factor(comm_data_summer$Hab)
scores_summer$Season=as.factor("Summer")
scores_winter$Season=as.factor("Winter")

scores=rbind(scores_summer,scores_winter)
scores$Season=factor(scores$Season, levels=c("Winter","Summer"))

df_ell3_summer$Season=as.factor("Summer")
df_ell3_winter$Season=as.factor("Winter")
df_ell3=rbind(df_ell3_summer,df_ell3_winter)
df_ell3$Season=factor(df_ell3$Season,levels=c("Winter","Summer"))

write.csv(scores,"Community composition/nmds.csv")
write.csv(df_ell3, "Community composition/ellipse.csv")

ggplot(scores,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(col=as.factor(group),size=0.001, pch=as.factor(Habitat)))+ # add the point markers  
  geom_path(data=df_ell3, aes(x=NMDS1, y=NMDS2,colour=group), size=1, 
            linetype=5)+
  xlim(c(-1.3,1.3))+
  ylim(c(-1.3,1.3))+
  theme_bw()+
  labs(title="NMDS ordination", col="Land-use Type")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  facet_grid(.~Season)

ggplot(scores_summer,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(col=as.factor(group),size=0.001, pch=as.factor(Habitat)))+ # add the point markers  
  geom_path(data=df_ell3_summer, aes(x=NMDS1, y=NMDS2,colour=group), size=1, 
            linetype=5)+
  xlim(c(-1.3,1.3))+
  ylim(c(-1.3,1.3))+
  theme_bw()+
  labs(title="Summer", col="Land-use Type")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#### Venn diagram #####

library(VennDiagram)

venn_winter=list(`Protected Grasslands`=unique(a$Species[a$Type=="Protected Grasslands"&a$Number>0&a$Season=="Winter"]),
          `Rangelands`=unique(a$Species[a$Type=="Rangelands"&a$Number>0&a$Season=="Winter"]),
                              `Non-irrigated Croplands`=unique(a$Species[a$Type=="Non-irrigated Croplands"&a$Number>0&a$Season=="Winter"]))
venn.diagram(venn_winter,filename = "venn_winter.png",
             output = FALSE ,
             imagetype="png" ,
             height = 480 , 
             width = 480 , 
             resolution = 300,
             compression = "lzw",
             lwd = 1,
             col=c("#440154ff", '#21908dff', '#fde725ff'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.3,
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
             rotation = 1
)
