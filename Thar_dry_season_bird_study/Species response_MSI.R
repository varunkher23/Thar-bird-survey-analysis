##### Loading packages ####
library(tidyverse)
library(reshape2)

#### Importing data ####
library(readxl)
a <- read_excel("A:/Thar/Final data/Data/Data.xlsx")

a=a%>%
  filter(Type!="Irrigated Croplands")

##### Adding blank records

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

a_species=merge(a,a_blank, all=T)
a_species$Encounter_rate=(a_species$Number)/(a_species$Transect_length*3/1000) ###Encounter rate per square km

##### Data Massaging #####

localabundance_species=a_species%>%
  group_by(Type,Season,Transect_Name,Species)%>%
  summarise(Encounter_rate=sum(Encounter_rate))

meanabundance_species=localabundance_species%>%
  group_by(Type,Season,Species)%>%
  summarise(mean=mean(Encounter_rate), 
            se=sd(Encounter_rate)/sqrt(length(Transect_Name)))

meanabundance_species$C=exp(1.645*sqrt(log(1+((meanabundance_species$se/meanabundance_species$mean)^2)))) ### 95% CI 
meanabundance_species$LL= meanabundance_species$mean/meanabundance_species$C 
meanabundance_species$UL= meanabundance_species$mean*meanabundance_species$C 

meanabundance_species=meanabundance_species%>%
  arrange(Species,Season,Type)

ER_mean=dcast(meanabundance_species,formula = Species+Season~Type, value.var = "mean")
colnames(ER_mean)=c("Species","Season","NIC_Mean","PG_Mean","RL_Mean")
ER_mean$NICdiff=(ER_mean$PG_Mean-ER_mean$NIC_Mean)/ER_mean$PG_Mean
ER_mean$RLdiff=(ER_mean$PG_Mean-ER_mean$RL_Mean)/ER_mean$PG_Mean

ER_LL=dcast(meanabundance_species,formula = Species+Season~Type, value.var = "LL")
ER_UL=dcast(meanabundance_species,formula = Species+Season~Type, value.var = "UL")

ER_CI=cbind(ER_LL,ER_UL)
ER_CI=ER_CI[,c(-6,-7)]
colnames(ER_CI)=c("Species","Season","NIC_LL","PG_LL","RL_LL","NIC_UL","PG_UL","RL_UL")

ER_CI=ER_CI%>%
  mutate(NICoverlap="NA")%>%
  mutate(RLoverlap="NA")%>%
  mutate(NICoverlap=ifelse(NIC_LL>PG_UL | NIC_UL<PG_LL,"No","Yes"))%>%
  mutate(RLoverlap=ifelse(RL_LL>PG_UL | RL_UL<PG_LL,"No","Yes"))

ER=merge(ER_CI,ER_mean)

ER=ER%>%
  mutate(NICSignificance="NA")%>%
  mutate(RLSignificance="NA")%>%
  mutate(NICSignificance=ifelse(NICoverlap=="No"& NICdiff>0.20,"Significant decrease","Not significant"))%>%
  mutate(NICSignificance=ifelse(NICoverlap=="No"& NICdiff<(-0.20),"Significant increase",NICSignificance))%>%
  mutate(NICSignificance=ifelse(NICdiff>0.999999&NICdiff<1.00001,"Absent",NICSignificance))%>%
  mutate(NICSignificance=ifelse(NICdiff=="-Inf","Presence",NICSignificance))%>%
  mutate(RLSignificance=ifelse(RLoverlap=="No"& RLdiff>0.20,"Significant decrease","Not significant"))%>%
  mutate(RLSignificance=ifelse(RLoverlap=="No"& RLdiff<(-0.20),"Significant increase",RLSignificance))%>%
  mutate(RLSignificance=ifelse(RLdiff>0.999999&RLdiff<1.00001,"Absent",RLSignificance))%>%
  mutate(RLSignificance=ifelse(RLdiff=="-Inf","Presence",RLSignificance))

ER=merge(ER,species)

ndetections=a_species%>%
  filter(is.na(ID)==FALSE)%>%
  group_by(Species,Season)%>%
  summarise(ndetections=n_distinct(ID))

ER=merge(ER,ndetections,all.y = T)

ER%>%
  filter(Season=="Winter")%>%
  filter(NICSignificance!="NA"| RLSignificance!="NA")%>%
  write.csv("Speciesresponse/Species_response_Winter.csv")

ER%>%
  filter(Season=="Summer")%>%
  filter(NICSignificance!="NA"| RLSignificance!="NA")%>%
  write.csv("Speciesresponse/Species_response_Summer.csv")
