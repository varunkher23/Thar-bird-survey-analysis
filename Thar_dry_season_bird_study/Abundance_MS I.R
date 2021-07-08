##### Loading packages ####
library(tidyverse)
library(vegan)
library(reshape2)
library(ggfortify)

#### Importing data ####
library(readxl)
a <- read_excel("D:/WII - M.Sc/Thar/Final data/Data/Data.xlsx")

a=a%>%
  filter(Type!="Irrigated Croplands")


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

a=merge(a,a_blank, all=T)

##### Data Massaging #####

abundance=a%>%
  filter(Perp_dist<=`Truncation distance`)
abundance=merge(abundance,a_blank,all=T)
abundance$Corrected_Count_per_sqkm=(abundance$Number/abundance$`Detection probability`)*(1000/abundance$`Truncation distance`)*(1/(abundance$Transect_length*3/1000))
abundance$Type=factor(abundance$Type, levels = c("Protected Grasslands", "Rangelands","Non-irrigated Croplands"))
abundance$Season=factor(abundance$Season, levels=c("Winter","Summer"))

##### Abundance - Local ####

localabundance_overall=abundance%>%
  group_by(Type, Transect_Name, Season)%>%
  summarise(abundance=sum(Corrected_Count_per_sqkm),Diet="Overall",Habitat="Overall")

localabundance_Diet=abundance%>%
  group_by(Type, Transect_Name, Season, Diet)%>%
  summarise(abundance=sum(Corrected_Count_per_sqkm),Habitat="Overall")

localabundance_Habitat=abundance%>%
  group_by(Type, Transect_Name, Season,Habitat)%>%
  summarise(abundance=sum(Corrected_Count_per_sqkm),Diet="Overall")

localabundance=rbind(localabundance_overall,localabundance_Diet,localabundance_Habitat)

meanabundance_overall <- localabundance%>%
  filter(Diet=="Overall")%>%
  filter(Habitat=="Overall")%>%
  group_by(Type,Season)%>%
  summarise(mean=mean(abundance), se=sd(abundance)/sqrt(length(unique(Transect_Name))),
            Habitat="Overall",Diet="Overall",Classification="Overall")

meanabundance_Habitat <- localabundance%>%
  filter(Diet=="Overall")%>%
  group_by(Type,Season,Habitat)%>%
  summarise(mean=mean(abundance), se=sd(abundance)/sqrt(length(unique(Transect_Name))),
            Diet="Overall",Classification="Diet")

meanabundance_Diet <- localabundance%>%
  filter(Habitat=="Overall")%>%
  group_by(Type,Season,Diet)%>%
  summarise(mean=mean(abundance), se=sd(abundance)/sqrt(length(unique(Transect_Name))),
            Habitat="Overall",Classification="Habitat")

meanabundance=rbind(meanabundance_Diet,meanabundance_Habitat,meanabundance_overall)

meanabundance$C=exp(1.96*sqrt(log(1+((meanabundance$se/meanabundance$mean)^2))))
meanabundance$LL=meanabundance$mean/meanabundance$C
meanabundance$UL=meanabundance$mean*meanabundance$C

write.csv(localabundance,"Localabundance/localabundance.csv")

###Statistical test###

abundance.lm.winter=lm(log(abundance)~Type,data = localabundance,subset = (Season=="Winter"& Diet=="Overall"& Habitat=="Overall") )
abundance.null.winter=lm(abundance~1,data = localabundance,subset = (Season=="Winter"& Diet=="Overall"& Habitat=="Overall") )
model.sel(abundance.lm.winter,abundance.null.winter)
summary(abundance.lm.winter)
Confint(abundance.lm.winter)
par(mfrow=(c(2,2)))
plot(abundance.lm.winter)

abundance.lm.summer=lm(log(abundance)~Type, data = localabundance,subset = (Season=="Summer"&Diet=="Overall"&Habitat=="Overall"))
summary(abundance.lm.summer)
Confint(abundance.lm.summer)
plot(abundance.lm.summer)
abundance.null.summer=lm(abundance~1,data = localabundance,subset = (Season=="Summer"& Diet=="Overall"& Habitat=="Overall") )
model.sel(abundance.lm.summer,abundance.null.summer)

newdat=data.frame(Type=unique(factor(localabundance$Type)))

localabundance.glm.winter.predict=cbind(newdat,as.data.frame(predict(abundance.lm.winter, newdata = newdat, type = "response", se.fit = T)),Season="Winter")
localabundance.glm.summer.predict=cbind(newdat,as.data.frame(predict(abundance.lm.summer, newdata = newdat, type = "response", se.fit = T)),Season="Summer")

localabundance.glm.predict=rbind(localabundance.glm.winter.predict,localabundance.glm.summer.predict)

localabundance.glm.predict$C=exp(1.96*sqrt(log(1+((localabundance.glm.predict$se.fit/localabundance.glm.predict$fit)^2))))
localabundance.glm.predict$LL=localabundance.glm.predict$fit/localabundance.glm.predict$C
localabundance.glm.predict$UL=localabundance.glm.predict$fit*localabundance.glm.predict$C

localabundance.glm.predict$fit=exp(localabundance.glm.predict$fit)
localabundance.glm.predict$LL=exp(localabundance.glm.predict$LL)
localabundance.glm.predict$UL=exp(localabundance.glm.predict$UL)
localabundance.glm.predict

write.table(localabundance.glm.predict,"Localabundance/meanabundance.csv")

Confint(abundance.lm.summer)

#### GlM- Diet and Habitat

###Diet
abundance.insectivore.lm.winter=lm(log(abundance+1)~Type,data = localabundance_Diet,subset = (Season=="Winter"& Diet=="Invertebrate"& Habitat=="Overall") )
abundance.granivore.lm.winter=lm(log(abundance+1)~Type,data = localabundance_Diet,subset = (Season=="Winter"& Diet=="PlantSeed"& Habitat=="Overall") )
abundance.omnivore.lm.winter=lm(log(abundance+1)~Type,data = localabundance_Diet,subset = (Season=="Winter"& Diet=="Omnivore"& Habitat=="Overall") )

summary(abundance.granivore.lm.winter)
exp(coef(abundance.insectivore.lm.winter))
par(mfrow=(c(2,2)))
plot(abundance.granivore.lm.winter)
shapiro.test(residuals(abundance.granivore.lm.winter))

predict.abundance.insectivore.winter=cbind(newdat,Diet="Insectivore",as.data.frame(predict(abundance.insectivore.lm.winter, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.granivore.winter=cbind(newdat,Diet="Granivore",as.data.frame(predict(abundance.granivore.lm.winter, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.omnivore.winter=cbind(newdat,Diet="Omnivore",as.data.frame(predict(abundance.omnivore.lm.winter, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.diet.winter=data.frame(rbind(predict.abundance.insectivore.winter,predict.abundance.granivore.winter,predict.abundance.omnivore.winter),Season="Winter",Classification="Diet")

abundance.insectivore.lm.summer=lm(log(abundance+1)~Type,data = localabundance_Diet,subset = (Season=="Summer"& Diet=="Invertebrate"& Habitat=="Overall"))
abundance.granivore.lm.summer=lm(log(abundance+1)~Type,data = localabundance_Diet,subset = (Season=="Summer"& Diet=="PlantSeed"& Habitat=="Overall") )
abundance.omnivore.lm.summer=lm(log(abundance+1)~Type,data = localabundance_Diet,subset = (Season=="Summer"& Diet=="Omnivore"& Habitat=="Overall") )

summary(abundance.granivore.lm.summer)

predict.abundance.insectivore.summer=cbind(newdat,Diet="Insectivore",as.data.frame(predict(abundance.insectivore.lm.summer, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.granivore.summer=cbind(newdat,Diet="Granivore",as.data.frame(predict(abundance.granivore.lm.summer, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.omnivore.summer=cbind(newdat,Diet="Omnivore",as.data.frame(predict(abundance.omnivore.lm.summer, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.diet.summer=data.frame(rbind(predict.abundance.insectivore.summer,predict.abundance.granivore.summer,predict.abundance.omnivore.summer),Season="Summer",Classification="Diet")

predict.abundance.diet=rbind(predict.abundance.diet.winter,predict.abundance.diet.summer)
predict.abundance.diet$Habitat="Overall"

predict.abundance.diet

#Habitat

abundance.grassland.lm.winter=lm(log(abundance+1)~Type,data = localabundance_Habitat,subset = (Season=="Winter"& Diet=="Overall"& Habitat!="Wetland" & Habitat!="Overall" & Habitat!="Generalist" & Habitat!="Forest") )
abundance.generalist.lm.winter=lm(log(abundance+1)~Type,data = localabundance_Habitat,subset = (Season=="Winter"& Diet=="Overall"& Habitat!="Wetland" & Habitat!="Overall" & Habitat!="Grassland" & Habitat!="Forest" & Habitat!="Desert") )
abundance.other.lm.winter=lm(log(abundance+1)~Type,data = localabundance_Habitat,subset = (Season=="Winter"& Diet=="Overall"& Habitat!="Grassland" & Habitat!="Overall" & Habitat!="Generalist" & Habitat!="Desert") )

predict.abundance.grassland.winter=cbind(newdat,Habitat="Grassland",as.data.frame(predict(abundance.grassland.lm.winter, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.generalist.winter=cbind(newdat,Habitat="Generalist",as.data.frame(predict(abundance.generalist.lm.winter, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.other.winter=cbind(newdat,Habitat="Others",as.data.frame(predict(abundance.other.lm.winter, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.habitat.winter=data.frame(rbind(predict.abundance.grassland.winter,predict.abundance.generalist.winter,predict.abundance.other.winter),Season="Winter",Classification="Habitat")

summary(abundance.generalistt.lm.winter)

abundance.grassland.lm.summer=lm(log(abundance+1)~Type,data = localabundance_Habitat,subset = (Season=="Summer"& Diet=="Overall"& Habitat!="Wetland" & Habitat!="Overall" & Habitat!="Generalist" & Habitat!="Forest") )
abundance.generalist.lm.summer=lm(log(abundance+1)~Type,data = localabundance_Habitat,subset = (Season=="Summer"& Diet=="Overall"& Habitat!="Wetland" & Habitat!="Overall" & Habitat!="Grassland" & Habitat!="Forest" & Habitat!="Desert") )
abundance.other.lm.summer=lm(log(abundance+1)~Type,data = localabundance_Habitat,subset = (Season=="Summer"& Diet=="Overall"& Habitat!="Grassland" & Habitat!="Overall" & Habitat!="Generalist" & Habitat!="Desert") )

summary(abundance.grassland.lm.Summer)

predict.abundance.grassland.summer=cbind(newdat,Habitat="Grassland",as.data.frame(predict(abundance.grassland.lm.summer, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.generalist.summer=cbind(newdat,Habitat="Generalist",as.data.frame(predict(abundance.generalist.lm.summer, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.other.summer=cbind(newdat,Habitat="Others",as.data.frame(predict(abundance.other.lm.summer, newdata = newdat, type = "response", se.fit = T)))
predict.abundance.habitat.summer=data.frame(rbind(predict.abundance.grassland.summer,predict.abundance.generalist.summer,predict.abundance.other.summer),Season="Summer",Classification="Habitat")

predict.abundance.habitat=rbind(predict.abundance.habitat.winter,predict.abundance.habitat.summer)
predict.abundance.habitat$Diet="Overall"

#Combining diet and habitat guild data

predict.abundance.guild=rbind(predict.abundance.diet,predict.abundance.habitat)

predict.abundance.guild$C=exp(1.96*sqrt(log(1+((predict.abundance.guild$se.fit/predict.abundance.guild$fit)^2))))
predict.abundance.guild$LL=predict.abundance.guild$fit/predict.abundance.guild$C
predict.abundance.guild$UL=predict.abundance.guild$fit*predict.abundance.guild$C

predict.abundance.guild$fit=exp(predict.abundance.guild$fit)
predict.abundance.guild$LL=exp(predict.abundance.guild$LL)
predict.abundance.guild$UL=exp(predict.abundance.guild$UL)
predict.abundance.guild

write.csv(predict.abundance.guild,"Results/Guild structure/guildstructure.csv")

##### Visialisation #####

meanabundance$Type=as.factor(meanabundance$Type)
meanabundance$Diet=as.factor(meanabundance$Diet)

bar_plot_diet=ggplot(subset(meanabundance,Diet!="Overall"&Diet!="FruiNect"&Diet!="VertFishScav"), aes(x=Diet, y=mean, fill=Type))+
  geom_bar(stat = "identity",position=position_dodge())+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"), 
        axis.line.x.top = element_line(colour = "black"))+
  ggtitle(label = "Diet wise local abundance")+
  geom_errorbar(aes(ymin=LL, 
                    ymax=UL), 
                width=.2,position=position_dodge(.9))+
  ylim(0,350)+
  xlab(label = "Type of Land-use")+
  ylab(label = "Mean abundance")+
  facet_grid(~Season)
plot(bar_plot_diet)

ggsave(path = "Plots/",
       filename = "Abundance_local_diet.tiff", 
       plot = bar_plot_diet, device = "tiff", dpi = 300, 
       width = 14, height=7 )

bar_plot_habitat=ggplot(subset(meanabundance,Habitat!="Overall"&Habitat!="Forest"&Habitat!="Wetland"), aes(x=Habitat, y=mean, fill=Type))+
  geom_bar(stat = "identity",position=position_dodge())+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"), 
        axis.line.x.top = element_line(colour = "black"))+
  ggtitle(label = "Habitat-wise abundance")+
  geom_errorbar(aes(ymin=LL, 
                    ymax=UL), 
                width=.2,position=position_dodge(.9))+
  ylim(0,400)+
  xlab(label = "Type of Land-use")+
  ylab(label = "Mean Species Richness")+
  facet_grid(~Season)
plot(bar_plot_habitat)

ggsave(path = "Plots/",
       filename = "Abundance_local_habitat.tiff", 
       plot = bar_plot_habitat, device = "tiff", dpi = 300, 
       width = 14, height=7 )

localabundance$Season=factor(localabundance$Season, levels=c("Winter","Summer"))
jitter_plot_abundance=ggplot(subset(localabundance,Habitat=="Overall"&Diet=="Overall"), aes(x=Type, y=abundance,col=Type))+
  theme_minimal()+
  theme(axis.line = element_line(colour = "black"), 
        axis.line.x.top = element_line(colour = "black"),
        panel.grid.minor = element_blank())+
  geom_dotplot(aes(y=abundance,x=Type,fill=Type,alpha=0.2),
               width = 0.10, alpha=0.3,
               binaxis = "y",stackdir = "center")+
  ggtitle(label = "Local species richness in each land cover type")+
  geom_errorbar(inherit.aes=F,data=localabundance.glm.predict,
                aes(x=Type,ymin=LL,ymax=UL,
                    width=.1))+
  geom_point(inherit.aes = F,data=localabundance.glm.predict,
             aes(x=Type,y=fit),col="Black")+
  xlab(label = "Type of Land-use")+
  ylab(label = "Mean Species Richness")+
  ylim(0,1500)+
  facet_grid(.~Season)
plot(jitter_plot_abundance)

ggsave(path = "Plots/",
       filename = "Abundance_local_overall.tiff", 
       plot = jitter_plot_abundance, device = "tiff", dpi = 300, 
       width = 14, height=7 )
