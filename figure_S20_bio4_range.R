#The R script written by Gaurav Baruah for the paper "Temperature tolerance and adaptive rewiring jointly predict plant–pollinator network architecture globally"

#figure S20- network architecture and temperature variability


#load libraries
library(bipartite)
library(ggplot2)
library(tidyverse)
library(colorBlindness)
library(ggfortify)
library(ggpmisc)
library(corrplot)
library(factoextra)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
d<-colorBlindness::displayAllColors(viridis::viridis(6))
#climatic data from the longitudes and latitudes of all the networks collated, from 
climatic_data<-read.csv("points_with_climate_data.csv")
head(climatic_data)
str(climatic_data)

climatic_data <- as.data.frame(climatic_data)
colnames(climatic_data)[colnames(climatic_data) == "wc2.1_2.5m_bio_1"] <- "MAT" #only interested in MAT, which is mean annual temperature!
colnames(climatic_data)[colnames(climatic_data) == "wc2.1_2.5m_bio_10"] <- "BIO10" 
#mean temperature of warmest quarter
colnames(climatic_data)[colnames(climatic_data) == "wc2.1_2.5m_bio_4"] <- "BIO4" 

#Correlation of MAT and mean temperature of the warmest quater

cor(climatic_data$MAT,climatic_data$BIO4,method = "pearson")
#-0.75

climatic_data$BIO4<-climatic_data$BIO4/100 #dividing by 100 to give SD in similar scales
#high network size --> high H2, more species less to 
conn<-climatic_data %>%
  ggplot(aes(y =Connectance , x = BIO4)) +
  geom_point(size = 3, alpha = 0.75,color="#01665e")  +
  xlab("SD of monthly temperatures")+
  ylab("Network Connectance") +
  theme_classic() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey2")+
  ggpmisc::stat_poly_eq(
    aes(label = paste (..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE,rsquared.conf.level = 0.95,
    coef.digits = 3, rr.digits = 2,
    f.digits = 3,
    p.digits = 3,label.x = "right", label.y = "top"
  ) #..rr.label..,
conn

nodf<-climatic_data  %>% 
  ggplot(aes(x=BIO4,y=NODF))+
  geom_point(size=3,alpha=0.75,color="#8c510a")+
  xlab("SD of monthly temperatures")+
  ylab("Network nestedness(NODF)")+
  theme_classic()+
  geom_smooth(method = "lm",formula = 'y~x',se=T, color ="grey2")+
  ggpmisc::stat_poly_eq(
    aes(label = paste( ..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE, label.x = "right", label.y = "top"
  )#..rr.label.., ..eq.label..,


h2<-climatic_data  %>%  
  ggplot(aes(x=BIO4,y=H2))+
  geom_point(size=3,alpha=0.75, color ="#b2182b")+
  xlab("SD of monthly temperatures")+
  ylab("Network specialisation (H2)")+
  theme_classic()+
  geom_smooth(method = "lm",formula = 'y~x',se=T, color ="grey2")+
  ggpmisc::stat_poly_eq(
    aes(label = paste(..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE, label.x = "right", label.y = "top"
  )#..eq.label..,

mod<-climatic_data  %>%  
  ggplot(aes(x=BIO4,y=modularity))+
  geom_point(size=3,alpha=0.75, color ="#D55E00")+
  xlab("SD of monthly temperatures")+
  ylab("Modularity")+
  theme_classic()+
  geom_smooth(method = "lm",formula = 'y~x',se=T, color ="grey2")+
  ggpmisc::stat_poly_eq(
    aes(label = paste(..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE, label.x = "right", label.y = "top"
  )#..eq.label..,
mod
ggpubr::ggarrange(conn,nodf, mod, h2,  
                  labels = c("A","B","C","D"),nrow=1,ncol=4)

