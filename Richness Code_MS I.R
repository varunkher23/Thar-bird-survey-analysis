##### Loading packages ####
library(tidyverse)
library(vegan)
library(reshape2)
library(ggfortify)
library(iNEXT)

#### Importing data ####
library(readxl)
a <- read_excel("A:/Thar/Final data/Data/Data.xlsx")
a$Type=factor(a$Type, levels = c("Protected Grasslands", "Rangelands","Non-irrigated Croplands","Irrigated Croplands"))
a=a%>%
  filter(Type!="Irrigated Croplands")
a$corrected_number=(a$Number/a$`Detection probability`)*(1000/a$`Truncation distance`)*(1/(a$Transect_length*3/1000))

length=data.frame(Transect_Name=a$Transect_Name,Type=a$Type,Transect_length=a$Transect_length)%>%
  unique.data.frame()

species=data_frame(Species=a$Species,Diet=a$Diet,Habitat=a$Habitat,
                   `Detection probability`=a$`Detection probability`,
                   `Truncation distance`=a$`Truncation distance`)%>%
  unique.data.frame()

a_blank=data.frame(Species=rep(unique(a$Species),length(unique(a$Transect_Name))*2))%>% ### 2 is because 2 seasons
  mutate(Transect_Name=(rep(unique(a$Transect_Name),each=length(unique(a$Species))*2)))%>% 
  arrange(Species,Transect_Name)%>%
  mutate(Season=rep(unique(a$Season),length(unique(a$Transect_Name))*length(unique(a$Species))))%>%
  mutate(Number=0)%>%
  merge(length,all=T)%>%
  merge(species,all=T)  
a=merge(a,a_blank,all=T)

comm_data_specpool=dcast(a,Type+Transect_Name+Season~Species, value.var = "Number",sum)

#### Species pool - Naive ####
Richness_naive=a%>%
  group_by(Type)%>%
  filter(Number>0)%>%
  summarise(Richness_naive_summer=n_distinct(Species[Season=="Summer"]),
            Richness_naive_winter=n_distinct(Species[Season=="Winter"]),
            Richness_naive_overall=n_distinct(Species))

#### Species Richness - Local ####

localspecr_overall=a%>%
  group_by(Type, Transect_Name, Season)%>%
  summarise(specr=n_distinct(Species[Number>0]), Diet="Overall",Habitat="Overall")

localspecr_Diet=a%>%
  group_by(Type, Transect_Name, Season, Diet)%>%
  summarise(specr=n_distinct(Species[Number>0]),Habitat="Overall")

localspecr_Habitat=a%>%
  group_by(Type, Transect_Name, Season, Habitat)%>%
  summarise(specr=n_distinct(Species[Number>0]),Diet="Overall")

localspecr=rbind(localspecr_overall,localspecr_Diet,localspecr_Habitat)

write.csv(localspecr, file="Results/Localspecr/localspecr.csv")

#### Mean Species Richness - Normal distribution ####

meanspecr_overall <- localspecr%>%
  filter(Diet=="Overall")%>%
  filter(Habitat=="Overall")%>%
  group_by(Type,Season)%>%
  summarise(mean=mean(specr), se=sd(specr)/sqrt(length(unique(Transect_Name))),
            Habitat="Overall",Diet="Overall",Classification="Overall")

meanspecr_Diet=localspecr%>%
  filter(Habitat=="Overall")%>%
  group_by(Type,Season,Diet)%>%
  summarise(mean=mean(specr), se=sd(specr)/sqrt(length(unique(Transect_Name))),
            Habitat="Overall",Classification="Diet")

meanspecr_Habitat=localspecr%>%
  filter(Diet=="Overall")%>%
  group_by(Type,Season,Habitat)%>%
  summarise(mean=mean(specr), se=sd(specr)/sqrt(length(unique(Transect_Name))),
            Diet="Overall",Classification="Habitat")

meanspecr=rbind(meanspecr_overall,meanspecr_Diet,meanspecr_Habitat)

### Statistical test ###

library(MASS)
library(lme4)
library(car)
library(MuMIn)

localspecr.glm.winter=glm(specr~Type, family = "poisson", data = localspecr,
                          subset = (Season=="Winter"& Diet=="Overall" & Habitat=="Overall") )
localspecr.null.winter=glm(specr~1, family = "poisson", data = localspecr,
                          subset = (Season=="Winter"& Diet=="Overall" & Habitat=="Overall") )
model.sel(localspecr.glm.winter,localspecr.null.winter)


localspecr.glm.summer=glm(specr~Type, family = "poisson", data = localspecr,
                          subset = (Season=="Summer"& Diet=="Overall" & Habitat=="Overall"))
localspecr.null.summer=glm(specr~1, family = "poisson", data = localspecr,
                          subset = (Season=="Summer"& Diet=="Overall" & Habitat=="Overall") )
model.sel(localspecr.glm.summer,localspecr.null.summer)

newdat=data.frame(Type=unique(factor(localspecr$Type)))
localspecr.glm.winter.predict=cbind(newdat,as.data.frame(predict(localspecr.glm.winter, newdata = newdat, type = "response", se.fit = T)),Season="Winter")
localspecr.glm.summer.predict=cbind(newdat,as.data.frame(predict(localspecr.glm.summer, newdata = newdat, type = "response", se.fit = T)),Season="Summer")

localspecr.glm.predict=rbind(localspecr.glm.winter.predict,localspecr.glm.summer.predict)

localspecr.glm.predict$C=exp(1.96*sqrt(log(1+((localspecr.glm.predict$se.fit/localspecr.glm.predict$fit)^2))))
localspecr.glm.predict$LL=localspecr.glm.predict$fit/localspecr.glm.predict$C
localspecr.glm.predict$UL=localspecr.glm.predict$fit*localspecr.glm.predict$C
localspecr.glm.predict

write.csv(localspecr.glm.predict, file  = "Results/Localspecr/meanspecr.csv")

#### Diet-wise GLM ####

#### Species accumulation ####

library(BiodiversityR)

C.estspecr_summer <- specaccum(comm_data_specpool[comm_data_specpool$Type=="Protected Grasslands"&comm_data_specpool$Season=="Summer",c(-1,-2,-3, -ncol(comm_data_specpool))], method = "rarefaction", xvar=c("sites","individuals"))
C.estspecr_summer=data.frame(cbind
                       (richness=as.numeric(as.character(
                         C.estspecr_summer[["richness"]]))),
                       sites=as.numeric(as.character(C.estspecr_summer[["sites"]])), 
                       sd=as.numeric(as.character(C.estspecr_summer[["sd"]])),
                       indivuduals=as.numeric(
                         as.character(C.estspecr_summer[["individuals"]])),
                       Type="Protected Grasslands")

GR.estspecr_summer <- specaccum(comm_data_specpool[comm_data_specpool$Type=="Rangelands"&comm_data_specpool$Season=="Summer",c(-1,-2,-3, -ncol(comm_data_specpool))], method = "rarefaction", xvar=c("sites","individuals"))
GR.estspecr_summer=data.frame(cbind
                        (richness=as.numeric(as.character(
                          GR.estspecr_summer[["richness"]]))),
                        sites=as.numeric(as.character(GR.estspecr_summer[["sites"]])), 
                        sd=as.numeric(as.character(GR.estspecr_summer[["sd"]])),
                        indivuduals=as.numeric(
                          as.character(GR.estspecr_summer[["individuals"]])),
                        Type="Rangelands")
NIA.estspecr_summer <- specaccum(comm_data_specpool[comm_data_specpool$Type=="Non-irrigated Croplands"&comm_data_specpool$Season=="Summer",c(-1,-2,-3, -ncol(comm_data_specpool))], method = "rarefaction", xvar=c("sites","individuals"))
NIA.estspecr_summer=data.frame(cbind
                      (richness=as.numeric(as.character(
                        NIA.estspecr_summer[["richness"]]))),
                      sites=as.numeric(as.character(NIA.estspecr_summer[["sites"]])), 
                      sd=as.numeric(as.character(NIA.estspecr_summer[["sd"]])),
                      indivuduals=as.numeric(
                        as.character(NIA.estspecr_summer[["individuals"]])),
                      Type="Non-irrigated Croplands")
Specr_accum_summer=data.frame(rbind(C.estspecr_summer, GR.estspecr_summer, NIA.estspecr_summer),Season="Summer")
Specr_accum_summer$sites=round(Specr_accum_summer$sites)
Specr_accum_summer=filter(Specr_accum_summer, sd!=0)

C.estspecr_winter <- specaccum(comm_data_specpool[comm_data_specpool$Type=="Protected Grasslands"&comm_data_specpool$Season=="Winter",c(-1,-2,-3, -ncol(comm_data_specpool))], method = "rarefaction", xvar=c("sites","individuals"))
C.estspecr_winter=data.frame(cbind
                             (richness=as.numeric(as.character(
                               C.estspecr_winter[["richness"]]))),
                             sites=as.numeric(as.character(C.estspecr_winter[["sites"]])), 
                             sd=as.numeric(as.character(C.estspecr_winter[["sd"]])),
                             indivuduals=as.numeric(
                               as.character(C.estspecr_winter[["individuals"]])),
                             Type="Protected Grasslands")

GR.estspecr_winter <- specaccum(comm_data_specpool[comm_data_specpool$Type=="Rangelands"&comm_data_specpool$Season=="Winter",c(-1,-2,-3, -ncol(comm_data_specpool))], method = "rarefaction", xvar=c("sites","individuals"))
GR.estspecr_winter=data.frame(cbind
                              (richness=as.numeric(as.character(
                                GR.estspecr_winter[["richness"]]))),
                              sites=as.numeric(as.character(GR.estspecr_winter[["sites"]])), 
                              sd=as.numeric(as.character(GR.estspecr_winter[["sd"]])),
                              indivuduals=as.numeric(
                                as.character(GR.estspecr_winter[["individuals"]])),
                              Type="Rangelands")
NIA.estspecr_winter <- specaccum(comm_data_specpool[comm_data_specpool$Type=="Non-irrigated Croplands"&comm_data_specpool$Season=="Winter",c(-1,-2,-3, -ncol(comm_data_specpool))], method = "rarefaction", xvar=c("sites","individuals"))
NIA.estspecr_winter=data.frame(cbind
                               (richness=as.numeric(as.character(
                                 NIA.estspecr_winter[["richness"]]))),
                               sites=as.numeric(as.character(NIA.estspecr_winter[["sites"]])), 
                               sd=as.numeric(as.character(NIA.estspecr_winter[["sd"]])),
                               indivuduals=as.numeric(
                                 as.character(NIA.estspecr_winter[["individuals"]])),
                               Type="Non-irrigated Croplands")
Specr_accum_winter=data.frame(rbind(C.estspecr_winter, GR.estspecr_winter, NIA.estspecr_winter),Season="Winter")
Specr_accum_winter$sites=round(Specr_accum_winter$sites)
Specr_accum_winter=filter(Specr_accum_winter, sd!=0)

Specr_accum=rbind(Specr_accum_summer,Specr_accum_winter)
Specr_accum$Season=factor(Specr_accum$Season,levels = c("Winter","Summer"))

write.csv(Specr_accum,file="Species pool/Specr_accum.csv")

#### Estimated species richness ####

c.estS <- specpool(c.comm_data[,c(-1,-2,-3)])
c.estS$Type <- list("Protected Grasslands")
GR.estS <- specpool(GR.comm_data[,c(-1,-2,-3)])
GR.estS$Type <- list("Rangelands")
NIA.estS <- specpool(NIA.comm_data[,c(-1,-2,-3)])
NIA.estS$Type <- list("Non-irrigated Croplands")
IA.estS <- specpool(IA.comm_data[,c(-1,-2,-3)])
IA.estS$Type <- list("Irrigated Croplands")
estS <- rbind(c.estS, GR.estS, NIA.estS, IA.estS, RA.estS)

##### Visualisation #####

#Overall site level richness
meanspecr$Season=factor(meanspecr$Season, levels=c("Winter","Summer"))
bar_plot_overall=ggplot(subset(meanspecr,Classification=="Overall"), aes(x=Type, y=mean, fill=Type))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"), 
        axis.line.x.top = element_line(colour = "black"))+
  ggtitle(label = "Local species richness in each land cover type")+
  geom_errorbar(aes(ymin=mean-se, 
                    ymax=mean+se), 
                width=.2,position=position_dodge(.9))+
  ylim(0.0,15)+
  xlab(label = "Type of Land-use")+
  ylab(label = "Mean Species Richness")+
  facet_grid(~Season)
plot(bar_plot)

#Diet-wise site level richness

bar_plot=ggplot(subset(meanspecr,Classification=="Diet"&Diet!="FruiNect"&Diet!="Overall"&Diet!="VertFishScav"), aes(x=Type, y=mean, fill=Type))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"), 
        axis.line.x.top = element_line(colour = "black"))+
  ggtitle(label = "Local species richness in each land cover type")+
  geom_errorbar(aes(ymin=mean-se, 
                    ymax=mean+se), 
                width=.2,position=position_dodge(.9))+
  xlab(label = "Type of Land-use")+
  ylab(label = "Mean Species Richness")+
  facet_grid(Diet~Season,scales = "free_y", space = "free_y")
plot(bar_plot)

bar_plot=ggplot(subset(meanspecr,Classification=="Habitat"&Habitat!="Overall"&Habitat!="Forest"&Habitat!="Wetland"), aes(x=Type, y=mean, fill=Type))+
  geom_bar(stat = "identity")+
  theme_void()+
  theme(axis.line = element_line(colour = "black"), 
        axis.line.x.top = element_line(colour = "black"))+
  ggtitle(label = "Local species richness in each land cover type")+
  geom_errorbar(aes(ymin=mean-se, 
                    ymax=mean+se), 
                width=.2,position=position_dodge(.9))+
  xlab(label = "Type of Land-use")+
  ylab(label = "Mean Species Richness")+
  facet_grid(Habitat~Season,scales = "free_y", space = "free_y")


localspecr$Season=factor(localspecr$Season, levels=c("Winter","Summer"))
jitter_plot=ggplot(subset(localspecr,Habitat=="Overall"&Diet=="Overall"), aes(x=Type, y=specr,col=Type))+
    theme_minimal()+
  theme(axis.line = element_line(colour = "black"), 
        axis.line.x.top = element_line(colour = "black"),
        panel.grid.minor = element_blank())+
  geom_dotplot(aes(y=specr,x=Type,fill=Type,alpha=0.2),
               width = 0.10, alpha=0.3,
               binaxis = "y",stackdir = "center")+
  ggtitle(label = "Local species richness in each land cover type")+
  geom_errorbar(inherit.aes=F,data=localspecr.glm.predict,
                aes(x=Type,ymin=LL,ymax=UL,
                    width=.1))+
  geom_point(inherit.aes = F,data=localspecr.glm.predict,
             aes(x=Type,y=fit),col="Black")+
  xlab(label = "Type of Land-use")+
  ylab(label = "Mean Species Richness")+
  ylim(0,20)+
  facet_grid(.~Season)

plot(jitter_plot)

ggsave(path = "Plots/",
       filename = "Specr_local_overall.tiff", 
       plot = jitter_plot, device = "tiff", dpi = 300, 
       width = 14, height=7 )

#Species Pool

ggplot(Specr_accum, aes(x=sites,y=richness, color=Type))+
  geom_line(size=2)+
  geom_errorbar(aes(ymin=richness-sd,ymax=richness+sd),width=.2 )+
  ylim(0,40)+
  xlim(0,26)+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"), 
        axis.line.x.top = element_line(colour = "black"))+
  ggtitle(label = "Rarified species pool in each land cover type")+
  xlab(label = "Sites")+
  ylab(label = "Estimated Species Richness")+
  facet_grid(~Season)
