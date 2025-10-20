#R script to reproduce Figure 5 in main-text, figure S2 of appendix, figure S14 of the appendix etc. The R script written by Gaurav Baruah for the paper "Adpative
#rewiring and temperature tolerance shapes the architecture of plant-pollinator networks globally"

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

#high network size --> high H2, more species less to 
conn<-climatic_data %>%
  ggplot(aes(y =H2 , x = MAT)) +
  geom_point(size = 3, alpha = 0.75,color="#01665e")  +
  xlab("MAT(°C)") +
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
  ggplot(aes(x=MAT,y=NODF))+
  geom_point(size=3,alpha=0.75,color="#8c510a")+
  xlim(c(-1,30))+
  xlab("MAT(°C)")+
  ylab("Network nestedness(NODF)")+
  theme_classic()+
  geom_smooth(method = "lm",formula = 'y~x',se=T, color ="grey2")+
  ggpmisc::stat_poly_eq(
    aes(label = paste( ..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE, label.x = "right", label.y = "top"
  )#..rr.label.., ..eq.label..,


h2<-climatic_data  %>%  
  ggplot(aes(x=MAT,y=H2))+
  geom_point(size=3,alpha=0.75, color ="#b2182b")+
  xlab("MAT(°C)")+
  ylab("Network specialisation (H2)")+
  theme_classic()+
  geom_smooth(method = "lm",formula = 'y~x',se=T, color ="grey2")+
  ggpmisc::stat_poly_eq(
    aes(label = paste(..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE, label.x = "right", label.y = "top"
  )#..eq.label..,

mod<-climatic_data  %>%  
  ggplot(aes(x=MAT,y=modularity))+
  geom_point(size=3,alpha=0.75, color ="#D55E00")+
  xlab("MAT (°C)")+
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


##################### Supplementary figure  S2- PCA #################




#some variables in the data is not named properly such as `nestedness`
princpdata<-climatic_data %>% filter(MAT > -10) %>%  select(Connectance,modularity,NODF,Species) 

#principal component analysis
res.pca <- prcomp(princpdata, scale = TRUE)
print(res.pca)
eig.val<-get_eigenvalue(res.pca)
eig.val
pc2<-fviz_eig(res.pca, col.var="blue")

# Color by cos2 values: quality on the factor map
PCA_PLOT<-fviz_pca_var(res.pca, col.var = "cos2",
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                       repel = TRUE # Avoid text overlapping
)

#this is for figure S1 inn the appendix
ggpubr::ggarrange(PCA_PLOT,pc2)


# scores of the pca axis that are not uncorrelated now also known as loadings
PC1_scores<-res.pca$x[,1] 
PC2_scores<-res.pca$x[,2]

cor(PC2_scores,princpdata$Species)
cor(PC1_scores,princpdata$Species)
cor(PC1_scores,princpdata$modularity)
cor(PC1_scores,princpdata$NODF)
cor(PC1_scores,princpdata$Connectance)

#cor(PC2_scores,princpdata$NODF)

climatic_data<-climatic_data %>%
  dplyr::filter(MAT > -10) %>% 
  dplyr::mutate(PC1 = PC1_scores,
                PC2 = PC2_scores)

#diagnoistic plots --- looks good!
climatic_data|>
  lm(PC1_scores ~ MAT, data = _) |>
  autoplot()
climatic_data|>
  lm(PC2_scores ~ MAT, data = _) |>
  autoplot()



#PC1 ---- positively correlated with network connectance + nestedness, negatively correlated with modularity
pc1<-climatic_data %>%
  ggplot(aes(y =PC1_scores , x = MAT)) +
  geom_point(size = 3, alpha = 0.75,color="#01665e")  +
  xlab("MAT(°C)") +
  ylab("PC1 scores") +
  theme_classic() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey2")+
  ggpmisc::stat_poly_eq(
    aes(label = paste (..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE,rsquared.conf.level = 0.95,
    coef.digits = 3, rr.digits = 2,
    f.digits = 3,
    p.digits = 3,label.x = "right", label.y = "top"
  ) #..rr.label..,
pc1

ggpubr::ggarrange(PCA_PLOT,pc1,labels = c("A","B"))



### relationship between PC1 and H2########### and reproduce FIGURE 5
pc1_h2<-climatic_data %>%
  ggplot(aes(y =H2 , x = MAT)) +
  geom_point(size = 3, alpha = 0.75,color ="#b2182b")  +
  xlab("MAT(°C)") +
  ylab("Network specialisation (H2')") +
  theme_classic() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "grey2")+
  ggpmisc::stat_poly_eq(
    aes(label = paste (..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE,rsquared.conf.level = 0.95,
    coef.digits = 3, rr.digits = 2,
    f.digits = 3,
    p.digits = 3,label.x = "right", label.y = "top"
  ) #..rr.label..,

pc1_h2

ggpubr::ggarrange(pc1,pc1_h2,labels = c("A","B"))



# Read coordinate data
coords <- read.csv("coordinates.csv")  # Must have Latitude and Longitude columns for all plant-poll networks compiled
head(coords)
# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot
r1<-ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
  geom_point(data = coords, aes(x = Longitude, y = Latitude),
             color = "darkred", size = 2.5, alpha = 0.7) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(title = "Plant–Pollinator Network Sites",
       x = "Longitude", y = "Latitude")


#Figure 5 of main-text
ggpubr::ggarrange(
  r1,
  ggpubr::ggarrange(pc1,pc1_h2, labels = c("B", "C"), ncol = 2),
  labels = "A",
  ncol = 1,
  nrow = 2,
  heights = c(1, 1)  )



## this is the figure without the PCA, i.e., FIGURE S14 in the appendix
ggpubr::ggarrange(
  r1,
  ggpubr::ggarrange(conn,nodf,h2,mod, labels = c("B", "C","D","E"), ncol = 4),
  labels = "A",
  ncol = 1,
  nrow = 2,
  heights = c(1, 1)  
)



#### Figure S3- relationship between PC1 and H2'

pc1_h23 <- climatic_data %>%
  dplyr::filter(MAT > -10) %>%
  ggplot(aes(y = H2, x = PC1_scores)) +
  geom_point(size = 3, alpha = 0.75, color = "#01665e") +
  xlab("PC1") +
  ylab("Network specialisation, H2'") +
  theme_classic() +
  geom_smooth(method = "glm", formula = y ~ x, se = TRUE, color = "grey2") +
  ggpmisc::stat_poly_eq(
    aes(label = paste(..p.value.label.., sep = "*\"; \"*")),
    formula = y ~ x, parse = TRUE, rsquared.conf.level = 0.95,
    coef.digits = 3, rr.digits = 2,
    f.digits = 3, p.digits = 3,
    label.x = "right", label.y = "top"
  ) +
  annotate("text", x = max(PC1_scores, na.rm = TRUE) - 0.5,
           y = min(climatic_data$H2, na.rm = TRUE) + 0.02,
           label = "(low network size, \n high connectance)", hjust = 1, size = 4) +
  annotate("text", x = min(PC1_scores, na.rm = TRUE) + 0.5,
           y = min(climatic_data$H2, na.rm = TRUE) + 0.02,
           label = "(high network size, \n low connectance)", hjust = 0, size = 4)

pc1_h23
