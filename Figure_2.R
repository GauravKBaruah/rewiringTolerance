rm(list = ls())
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
library(scales)
library(network)
library(GGally)
set.seed(123)


g <- matrix(1, nrow = 11, ncol = 10) #ork web names
Aspecies <- dim(g)[2] # no of animal species
Plantspecies <- dim(g)[1] # no of plant species
degree.animals <- degree.plants <- numeric()


#parameters for modelling: intialisation of the model
## initial species densities
Temperature <- 28
N <-  rep(1, each = (Aspecies + Plantspecies))
na <-N[1:Aspecies] 
np <- N[(Aspecies + 1):(Aspecies + Plantspecies)]#
muA <-  runif((Aspecies), (Temperature - 4.5), (Temperature +4.5))# mu[1:Aspecies]#initial mean phenotypic optimum trait values
muP <-runif((Plantspecies), (Temperature - 4.5), (Temperature +4.5))
tmax <- 5000
time_range <- c(0, tmax)
nestedness <- nestedness_NODF(g)
C <- Connectance(g)
dganimals <- degree.animals #degree animals - all species have same degree
dgplants <- degree.plants #degree plants - all species have same degree
alpha_ii <- 1 #density dependence term a_ii

# trait variances 
sig <- runif((Aspecies + Plantspecies), 0.1, 0.25)
envA <- runif(Aspecies, 0.009, 0.01) #environmental variances for animals
envP <- runif(Plantspecies, 0.009, 0.01) #env variances for plants
sa <- runif((Aspecies), 0.1, 0.25) # gen variances for animals
sp <- runif(( Plantspecies), 0.1, 0.25) # gen variances for plants
Pa <- envA + sa #phenotypic variances animals
Pp <- envP + sp # phenotypic variances plants
trait_evolve <- 1
  
#tolerance parameters of temperature curves
bw <- 4.5 #width of the tolerance curve
gi <- 1.5 # peak of tolerance temperature curve
w <- 1 # mutualistic gaussian kernel width
a <- 0.1 # tolerance temperature curve parameter that shapes the temperature curve
w_c <- 0.05 # inter specific competition width
opTemp<- Temperature # environmental temperature
ic <-c(na, np, muA, muP,  opTemp ) #inital state values 
mut.strength <- 1 #avg. mutualistic strength
params <- #list of parameters
  list(
    time = time,
    matrix = g,
    pref_matrix=g,
    Temp = Temperature,
    alpha_ii = 1,
    w_c = w_c,
    ki = 0.1,
    gi = gi,
    bw = bw,
    a = a,
    tmax = tmax,
    sig = sig,
    w = w,
    trait_evolve = trait_evolve,
    mut.strength = mut.strength,
    C = C,
    gvarA=sa,
    gvarP=sp,
    nestedness = nestedness,
    dganimals = degree.animals,
    degree = degree,
    dgplants = degree.plants,
    envA = envA,
    envP = envP
  )

#   initial snapshot of the trait distribution of species in a community
t1<- plot_snapshot(Na=na, Np=np, m = c(muA,muP), 
                   sigma = c(Pa,Pp), Temp=Temperature, limits = c(22,33),res = 1001)+
  annotate("text",x=23, y =5, label="time == 0",parse=T)  

t1

#calculating the overlap between plants and pollinators and then creating the incidence matrix
tem_mets<-quantitative_network(muA = muA,muP = muP, 
                     Na = na, Np = np, sa = (Pa),
                     sp = (Pp),w = 1,thr = 0.2)

#Gower variance incidence matrix
incidence_mat<-tem_mets$web_matrix_gowervariance
tem_mets$connectance_true_est
# plotting the network initially based on overlap and similarity
#degree of plants and anichmals
for (i in 1:Plantspecies) {
  degree.plants[i] <- sum(incidence_mat[i, ])
} # degree of plants
for (j in 1:Aspecies) {
  degree.animals[j] <- sum(incidence_mat[, j]) # degree of animals
}

degree <- c(degree.animals, degree.plants)

# network graph creation based on the gower variance incidence matrix
net1 <- network(
  incidence_mat,
  bipartite    = TRUE,
  directed     = FALSE,
  matrix.type  = "bipartite",
  ignore.eval  = FALSE,        
  names.eval   = "weight"   
)

#names of vertices
names<-network.vertex.names(net1)
net1 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net1 %v% "color" = ifelse(net1 %v% "groups" == "plants","#b2182b", "#E69F00" )

#plotting the network
deg<-c(degree.plants,degree.animals)
(w1<-ggnet2(net1, mode="circle", node.size = 3,
            color ="color", edge.alpha = 1, 
            legend.position = "",label = "A"))

w1


# simulating eco-evo dynamics of the network
tmax <- tmax #total time points for simulation
#dynamics from ode solve
OUT<-ode(func=eqs_type2_new, y=ic, parms=params,
         times=seq(0, tmax, by=tmax/100),rtol = 1e-9, atol = 1e-12) %>% 
  organize_results_rewire(pars = params) 

#getting the state variables out for maximum timepoint
sol_1<-  OUT %>% filter(time == tmax)
Na1<- (sol_1 %>% filter(type=="N"))$v #density of animals
Np1<-(sol_1 %>% filter(type =="P"))$v #density of plants
ma1<-(sol_1 %>% filter(type == "muA"))$v #trait values of animals
mp1<-(sol_1 %>% filter(type == "muP"))$v #trait values of plants

#phenotypic trait variances
sa1<-params$gvarA
sp1<-params$gvarP
Pa1<-sa1 + params$envA
Pp1<-sp1 +params$envP

#plotting trait distributions at final timepoint
t2<- plot_snapshot(Na=Na1, Np=Np1, m = c(ma1,mp1), 
                   sigma = c(Pa1,Pp1), Temp=28, limits = c(23,31),res = 1001)+
    annotate("text",x=29.5, y =1.85, label="time == 10000",parse=T)  

t2

#getting the network metrics out from the overlap matrices
q_metrics<-quantitative_network(muA = ma1,muP = mp1,Na = Na1,Np = Np1,sa = Pa1,sp = Pp1, w=1 ,thr = 0.2)

# adjancency matrix based on gower variances
incidence_matrix<-q_metrics$web_matrix_gowervariance

 
#degree of plants and anichmals
for (i in 1:Plantspecies) {
  degree.plants[i] <- sum(incidence_matrix[i, ])
} # degree of plants
for (j in 1:Aspecies) {
  degree.animals[j] <- sum(incidence_matrix[, j]) # degree of animals
}



# network graph creation based on the gower variance incidence matrix
net2 <- network(
  incidence_matrix,
  bipartite    = TRUE,
  directed     = FALSE,
  matrix.type  = "bipartite",
  ignore.eval  = FALSE,        # KEEP the edge values
  names.eval   = "weight"      # store them in an edge attribute called "weight"
)


names<-network.vertex.names(net2)
net2 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(incidence_matrix)[1])) ), "plants", "animals")
net2 %v% "color" = ifelse(net2 %v% "groups" == "plants","#b2182b", "#E69F00" )

deg<-c(degree.plants,degree.animals)
(w2<-ggnet2(net2, mode="circle", node.size = 3,
            color ="color", edge.alpha = 5, 
            legend.position = "",label = "A"))

w2
net_h_temp<-ggpubr::ggarrange(w1,t1,t2,w2, nrow=1, 
                  labels = c("B","C","D","E"),
                  widths = c(1, 1.75, 1.75, 1))




################################# For 2 degree C temperature ##############################
Temperature <- 2 # env. temperature
na <-N[1:Aspecies]  # densities of animals
np <- N[(Aspecies + 1):(Aspecies + Plantspecies)] #densities of plants

muA <-  runif((Aspecies), (Temperature - 4.5), (Temperature +4.5))# mu[1:Aspecies]#initial mean phenotypic optimum trait values
muP <-runif((Plantspecies), (Temperature - 4.5), (Temperature +4.5))

#tolerance temperature parameters
bw <- 4.5 #width of the tolerance curve
gi <- 1.5 # peak of tolerance temperature curve
w <- 1 # mutualistic gaussian kernel width
a <- 0.1 # tolerance temperature curve parameter that shapes the temperature curve
w_c <- 0.05 # inter specific competition width
opTemp<- Temperature # environmental temperature
ic <-c(na, np, muA, muP,  opTemp ) #inital state values 
mut.strength <- 1 #avg. mutualistic strength
#list of parameters
params <-
  list(
    time = time,
    matrix = g,
    pref_matrix=g,
    Temp = Temperature,
    alpha_ii = 1,
    w_c = w_c,
    ki = 0.1,
    gi = gi,
    bw = bw,
    a = a,
    tmax = tmax,
    sig = sig,
    w = w,
    trait_evolve = trait_evolve,
    mut.strength = mut.strength,
    C = C,
    gvarA=sa,
    gvarP=sp,
    nestedness = nestedness,
    dganimals = degree.animals,
    degree = degree,
    dgplants = degree.plants,
    envA = envA,
    envP = envP
  )


#plotting trait distributions at initial timepoint
t3<- plot_snapshot(Na=na, Np=np, m = c(muA,muP), 
                   sigma = c(Pa,Pp), Temp=2, limits = c(-5,9),res = 1001)+
  annotate("text",x=2, y =5, label="time == 0",parse=T)  

t3


#calculating the overlap between plants and pollinators and then creating the incidence matrix
tem_mets3<-quantitative_network(muA = muA,muP = muP, 
                               Na = na, Np = np, sa = (Pa),
                               sp = (Pp),w = 1,thr = 0.2)

#Gower incidence matrix
incidence_mat3<-tem_mets3$web_matrix_gowervariance

# plotting the network initially based on overlap and similarity
#degree of plants and anichmals
for (i in 1:Plantspecies) {
  degree.plants[i] <- sum(incidence_mat3[i, ])
} # degree of plants
for (j in 1:Aspecies) {
  degree.animals[j] <- sum(incidence_mat3[, j]) # degree of animals
}

degree <- c(degree.animals, degree.plants)


# network graph creation based on the gower variance incidence matrix
net3 <- network(
  incidence_mat3,
  bipartite    = TRUE,
  directed     = FALSE,
  matrix.type  = "bipartite",
  ignore.eval  = FALSE,        # KEEP the edge values
  names.eval   = "weight"      # store them in an edge attribute called "weight"
)


names<-network.vertex.names(net3)
net3 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(incidence_mat3)[1])) ), "plants", "animals")
net3 %v% "color" = ifelse(net3 %v% "groups" == "plants","#b2182b", "#E69F00" )

deg<-c(degree.plants,degree.animals)
(w3<-ggnet2(net3, mode="circle", node.size = 3,
            color ="color", edge.alpha = 1, 
            legend.position = "",label = "A"))

w3

# simulating eco-evo dynamics of the network
tmax <- tmax #total time points for simulation
#dynamics from ode solve
OUT2<-ode(func=eqs_type2_new, y=ic, parms=params,
         times=seq(0, tmax, by=tmax/1000), rtol = 1e-9, atol = 1e-12) %>% 
  organize_results_rewire(pars = params) 


#getting the state variables out for maximum timepoint
sol_2<-  OUT2 %>% filter(time == tmax)
Na2<- (sol_2 %>% filter(type=="N"))$v
Np2<-(sol_2 %>% filter(type =="P"))$v
ma2<-(sol_2 %>% filter(type == "muA"))$v
mp2<-(sol_2 %>% filter(type == "muP"))$v

#phenotypic trait variances
sa2<-params$gvarA+params$envA
sp2<-params$gvarP+params$envP
Pa2<-sa2
Pp2<-sp2

#plotting trait distributions at final timepoint
t4<- plot_snapshot(Na=Na2, Np=Np2, m = c(ma2,mp2), 
                   sigma = c(sa2,sp2), Temp=Temperature, limits = c(-5,9),res = 1001)+
  annotate("text",x=5, y =.8, label="time == 10000",parse=T)  


t4
#network metrics from inter guild overlaps of trait values
q_metrics_2<-quantitative_network(muA = ma2,muP = mp2,Na = Na2,Np = Np2,sa = Pa2,sp = Pp2,w=1,thr = 0.2)

#incidence matrix from overlap of species traits
incidence_matrix_3<-q_metrics_2$web_matrix_gowervariance



#degree of plants and anichmals
for (i in 1:Plantspecies) {
  degree.plants[i] <- sum(incidence_matrix_3[i, ])
} # degree of plants
for (j in 1:Aspecies) {
  degree.animals[j] <- sum(incidence_matrix_3[, j]) # degree of animals
}

degree <- c(degree.animals, degree.plants)

#network from adjancency matrix derived from gower-variances 
net4 <- network(  incidence_matrix_3,
  bipartite    = TRUE,
  directed     = FALSE,
  matrix.type  = "bipartite",
  ignore.eval  = FALSE,       
  names.eval   = "weight")      

#plotting
names<-network.vertex.names(net4)
net4 %v% "groups" = ifelse( names[1:sum(dim(g))] %in% c( as.character(seq(1:dim(g)[1])) ), "plants", "animals")
net4 %v% "color" = ifelse(net4 %v% "groups" == "plants","#b2182b", "#E69F00" )

deg<-c(degree.plants,degree.animals)*0.25

(w4<-ggnet2(net4, mode="circle", node.size = 3,
            color ="color", edge.alpha = 5, 
            legend.position = "",label = "A"))

w4


net_low_temp<-ggpubr::ggarrange(w3,t3,t4,w4, nrow=1, 
                                labels = c("F","G","H","I"),
                                widths = c(1, 1.75, 1.75, 1))

ggpubr::ggarrange(net_h_temp,net_low_temp, nrow=2)




################################ function pgr ################################

# 1) define your parameters
muA   <- c(0, 7, 14, 21,28)
varA  <- rep(0.5, length(muA))
Temp  <- seq(-2, 32, 0.5)

ki <- 0.1
gi <- 1.5
a  <- 0.1
bw <- 4.5

# 2) build the data.frame of all Temp × muA combos
df <- expand.grid(
  Temp = Temp,
  muA  = muA
)
df$varA <- 0.25  # same length as muA, constant here

# 3) compute b for each row
df$b <- with(df,
             gi / (bw - a * muA) *
               (bw - a * muA) / sqrt((bw - a * muA)^2 + varA) *
               exp( - (Temp - muA)^2 / (2 * ( (bw - a * muA)^2 + varA )) ) -
               ki
)

# 4) plot
r1<-ggplot(df, aes(x = Temp, y = b, group = muA)) +
  geom_ribbon(aes(ymin = 0, ymax = b, fill = muA), alpha = 0.2, colour = NA) +
  geom_line(aes(color = muA), size = 1) +
  ylim(c(-0.15, 0.75)) +
  scale_color_gradient(name = expression(mu[A]),low = "skyblue", high = "salmon") +
  scale_fill_gradient(name = expression(mu[A]),low = "skyblue", high = "salmon") +
  labs(
    x = "Temperature",
    y = "Growth rate",
    title = ""
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()

r1

blankPlot <- ggplot() + theme_void()


ss<-ggpubr::ggarrange(blankPlot,r1,blankPlot, ncol=3)

ggpubr::ggarrange(ss,net_h_temp,net_low_temp,nrow = 3,labels = "A",
                  
                  heights = c(0.5, 1, 1) )
  
