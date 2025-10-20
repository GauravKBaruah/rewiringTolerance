#R script to reproduce Figure 3 and Figure 4 in main-text.The R script written by Gaurav Baruah for the paper "Adpative
#rewiring and temperature tolerance shapes the architecture of plant-pollinator networks globally" - authors: Gaurav Baruah, Meike Wittmann

#loading libraries need for this R script
source('01_functions_rewiring.R')
require(deSolve) 
require(cowplot) 
library(dplyr)
library(readr)
library(beepr)
library(ggplot2)
library(viridis)
library(igraph)
library(bipartite)
library(sna)
library(network)
library(GGally)

#load simulated data from for 20 species networks
load("network_metrics_20species.RData")

#stretch the data to see whats there
str(sp_data)  
unique(sp_data$width_competition)
unique(sp_data$a)

#relabel factors and columns
a_labels <- c( "0"="a=0", "0.1"="a=0.1",  
              "0.125"= "a=0.125")

#network metrics- connectance gower variance
conn_gv<-sp_data %>% 
  ggplot(aes(x=Temperature_shift,y=Connectance_gowervariance_unweighted, color=factor(width_competition)))+
  geom_point(size=2, alpha=0.3,shape=1)+
  stat_smooth(method = "lm",
              formula = "y~x",
              se = TRUE) + theme_classic()+
  labs(color="Competition \n width")+
  scale_color_manual(values = c("#01665e", "#01305e", "#CC79A7")) + # Okabe-Ito colors
  xlab(expression(paste("Temperature (",degree ~C, ")")))+
  ylab("Connectance ")+
  facet_wrap(.~a, scales = "free",
             labeller = labeller(
               a = as_labeller(a_labels)))

conn_gv

#nestedness by gowervariance metric
nodf_gv<-sp_data %>%  
  ggplot(aes(x=Temperature_shift,y=nestedness_gowervariance_unweighted, color=factor(width_competition)))+
  geom_point(size=2, alpha=0.3,shape=1)+
  # smooth curve via quasibinomial GLM
  stat_smooth(method = "lm",formula = "y~x",
              se = TRUE)+ theme_classic()+
  labs(color="Competition \n width")+
  scale_color_manual(values = c("#01665e", "#01305e", "#CC79A7")) + # Okabe-Ito colors
 # scale_color_manual(values = c("#8c510a", "#D55E00", "#CC79A7")) + # Okabe-Ito colors
  xlab(expression(paste("Temperature (",degree ~C, ")")))+
  ylab("Nestedness \n (NODF) ")+
  facet_wrap(.~a, scales = "free",
             labeller = labeller(
               a = as_labeller(a_labels)))
nodf_gv


#network metrics- H2 by gowervariance metric
h2_gv<-sp_data %>% 
  ggplot(aes(x=Temperature_shift,y=H2_gowervariance_int, color=factor(width_competition)))+
  geom_point(size=2, alpha=0.3,shape=1)+
  # smooth curve via quasibinomial GLM
  stat_smooth(method = "lm",formula = "y~x",
              se = TRUE)+ theme_classic()+
  labs(color="Competition \n width")+
  scale_color_manual(values = c("#01665e", "#01305e", "#CC79A7")) + # Okabe-Ito colors
  
  #scale_color_manual(values = c("#b2182b","#a1502a", "#CC79A7")) + # Okabe-Ito colors
  xlab(expression(paste("Temperature (",degree ~C, ")")))+
  ylab("Specialisation \n (H2') ")+
  facet_wrap(.~a, scales = "free",
             labeller = labeller(
               a = as_labeller(a_labels)))
h2_gv

#modularity estimation from the gowervariance metric
mod_gv<-sp_data %>%   
  ggplot(aes(x=Temperature_shift,y=modularity_gowervariance, color=factor(width_competition)))+
  geom_point(size=2, alpha=0.3,shape=1)+
  # smooth curve via quasibinomial GLM
  stat_smooth(method = "lm",formula = "y~x",
              se = TRUE) + theme_classic()+
  labs(color="Competition \n width")+
  scale_color_manual(values = c("#01665e", "#01305e", "#CC79A7")) + # Okabe-Ito colors
  
 # scale_color_manual(values = c("#0072B2", "#D55E00", "#CC79A7")) + # Okabe-Ito colors
  xlab(expression(paste("Temperature (",degree ~C, ")")))+
  ylab("Modularity ")+
  facet_wrap(.~a, scales = "free",
             labeller = labeller(
               a = as_labeller(a_labels)))
mod_gv


#plotting all together
ggpubr::ggarrange(conn_gv,nodf_gv, mod_gv,
                  h2_gv, labels = c("A","B","C","D"),
                  common.legend = TRUE,legend = "bottom",nrow=2,ncol=2)



############################# REWIRING DATA ############################
load("rewiring_data_20species.RData") #rewiring data for 20 species simulated plant-pollinator network 
df <- net_data %>%
  mutate(time = time * 50) #in ODE solver, time units were 50 time points apart, so this multiplication is a scaling factor

#arranging the dataframe from the simulations
df2 <- df %>%
  arrange(a, network_size, starting_temperature, competition, time) %>%
  group_by(a, network_size, starting_temperature, competition) %>%
  mutate(replicate = cumsum(time == 50)) %>%
  ungroup()

# creating a dataframe of mean rewiwirng
df_mean <- df %>%
  group_by(a, network_size, starting_temperature, competition, time) %>%
  summarise(
    mean_rewiring = mean(Rewiring, na.rm = TRUE), #median rewiring across reps
    sd  = sd(Rewiring,   na.rm = TRUE),
    n   = dplyr::n(),
    se  = sd / sqrt(40), #40 is the replicaties
    lower = mean_rewiring -  1.96*se,
    upper = mean_rewiring +  1.96*se,
    .groups = "drop"
  ) %>%
  # if n==1, CI is undefined; fall back to the mean (flat ribbon)
  mutate(
    lower = ifelse(is.finite(lower), lower, mean_rewiring),
    upper = ifelse(is.finite(upper), upper, mean_rewiring)
  )

#raw rewiring over time, all reps plotteed
rewiring_dyn<-ggplot(df2, 
                     aes(x     = time,
                         y     = Rewiring,
                         group = replicate,
                         colour = factor(starting_temperature))) +
  geom_point(alpha = 0.5) +
  facet_grid(a ~ competition,
             labeller = labeller(a = label_both,
                                 competition = label_both)) +
  scale_colour_viridis_d(option = "magma",name = expression(paste("Temperature (",degree~C,")"))) +
  theme_cowplot() +
  labs(x = "Time", y = "Rewiring")

rewiring_dyn


# mean rewiring across replicates, so mean and and upper and lower intervals .e.e sd error plotted over time
rewiring_trend <- ggplot(
  df_mean,
  aes(x = time,
      y = mean_rewiring,
      colour = factor(starting_temperature),
      fill   = factor(starting_temperature))
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20, colour = NA) +
  geom_line(size = 1) +
  facet_grid(a ~ competition,
             labeller = labeller(a = label_both,
                                 competition = label_both)) +
  scale_colour_viridis_d(option = "magma",
                         name = expression(paste("Temperature (", degree*C, ")"))) +
  scale_fill_viridis_d(option = "magma", guide = "none") +
  theme_classic() +
  labs(x = "Time", y = "Mean rewiring  \n across replicates")

rewiring_trend #rewiring over time with mean and sd over replicates

#reconstruct dataframe
df2 <- net_data %>%
  arrange(a, starting_temperature, competition, time) %>%
  group_by(a, starting_temperature, competition) %>%
  mutate(replicate = cumsum(time == 1)) %>%
  ungroup()

# compute median rewiring over time per replicate per temperature 
df_med <- df2 %>% 
  group_by(a, competition, starting_temperature, replicate) %>% #mutate(final_rewiring)
  summarise(median_rewiring = median(Rewiring, na.rm=TRUE),
            .groups = "drop")

#summarise(median_rewiring = max(Rewiring, na.rm=TRUE),
#         .groups = "drop")

#summarise across replicates at each temp/ treatment combo
df_trend <- df_med %>%
  group_by(a, competition, starting_temperature) %>%
  summarise(
    mean_rw = mean(median_rewiring),
    sd_rw   = sd(median_rewiring),
    n_rep   = n(),
    se_rw   = sd_rw / sqrt(n_rep),
    lower   = mean_rw - se_rw,
    upper   = mean_rw + se_rw,
    .groups = "drop"
  )


rewiring_median<-ggplot(df_trend,
                        aes(x = starting_temperature,
                            y = mean_rw,
                            colour = factor(a),
                            fill   = factor(a),
                            group  = factor(a))) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.3, colour = NA) +
  geom_point()+
  geom_line(size = 1.5) +
  facet_wrap(~ competition,
             labeller = labeller(competition = label_both)) +
  #  scale_colour_viridis_d(option= "plasma",name = "a") +
  # scale_fill_viridis_d(option="plasma",name = "a") +
  scale_color_brewer(palette = "Set1", name = "a")+
  scale_fill_brewer(palette = "Set1", name = "a")+
  theme_classic() +
  labs(x = expression(paste("Temperature (", degree~C, ")")),
       y = "Median rewiring \n per temperature",
       colour = "a",
       group = "a")
rewiring_median

ggpubr::ggarrange(rewiring_trend,rewiring_median,nrow=2,
                  labels = c("A","B"))

