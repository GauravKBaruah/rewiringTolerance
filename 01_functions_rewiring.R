
#R script functions to be used to reproduce figures and simulate dynbamics. The R script written by Gaurav Baruah for the paper "Adpative
#rewiring and temperature tolerance shapes the architecture of plant-pollinator networks globally"
require(statmod)

#cutoff to remove negative densities from odesolveer  
cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)

  adj.mat<-function(data){
    #dat <- paste('network.csv',sep='')
    d <- read.csv(file=data,header=FALSE )
    dat<-as.matrix(d)
    dat[dat > 0] = 1
    dat<-apply(dat,2,as.numeric)
    return(dat)}

# To solve the double integral of the type-2 curve of plant-animal interaction is not possible analytically. So we employ a numerical approach with 
# a powerful solver such as a Gaussian Quadrature approach. And the function below does that.
  
# type 2 functional curve estimation from Gaussian quadrature
#input: m - trait values
# sigma - trait variances
# w - mutualistic gaussian interaction width
# mut.strenght- avg. strenght o plant-pollinator interaction
#a_index- index value of which animal species  is in the loop of species interaction.
#points - no. of points used to approximate the normal distributions: 6 is perfect! below that is a bit shaky!
# mat, pref_matrix - both incidence matrix of which plant species interacts with which pollinator species, which in our case is all 1's
# degree.animal - all species have same degree because potentiall any species can interact with anybody based on the trait values they possess and how similar they are to others.
#Pa -  phenotypic variance of plants
type_2_animals<-function(m,sigma,w,h,np,na,mut.strength,a_index,
                           points,mat,pref_matrix,degree.animal,Pa){
    temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma[a_index], sigma =sqrt(sigma$sa[a_index]))$nodes #z of animals
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$ma[a_index],sigma =sqrt(sigma$sa[a_index]))$weights #p_i(z)
    
    z2<-matrix(0,nrow=length(np),ncol=points)
    w2<-matrix(0,nrow=length(np),ncol=points)
    numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(np))
    N_strength<-m_strength<-g_strength<-numeric()
    
    for(j in 1:points){  
      for(k in 1:length(np)){
        z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$mp[k], sigma =sqrt(sigma$sp[k]))$nodes #z''
        
        #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
        w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                                mu=m$mp[k],sigma =sqrt(sigma$sp[k]))$weights #pj(z'')
        
        numer_a[j,k]<- np[k]*mat[k,a_index]*pref_matrix[k,a_index]*(mut.strength/degree.animal)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        denom_a[j,k]<- np[k]*mat[k,a_index]*pref_matrix[k,a_index]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        
        numer_m[j,k]<- np[k]*mat[k,a_index]*pref_matrix[k,a_index]*(mut.strength/degree.animal)*sum((z1[j]-m$ma[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        denom_m[j,k]<- np[k]*mat[k,a_index]*pref_matrix[k,a_index]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        
        numer_g[j,k]<- np[k]*mat[k,a_index]*pref_matrix[k,a_index]*1/Pa[a_index]^2*(mut.strength/degree.animal)*sum(((z1[j]-m$ma[a_index])^2-Pa[a_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        denom_g[j,k]<- np[k]*mat[k,a_index]*pref_matrix[k,a_index]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        
      }
      
      N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
      m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
      g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
    }
    
    
    G = sum(N_strength)
    B = sum(m_strength) 
    V = sum(g_strength)
    
    
    
    return(list(G= G, B = B,V=V))
    
  }

# type 2 functional curve estimation from Gaussian quadrature
#input: m - trait values
# sigma - trait variances
# w - mutualistic gaussian interaction width
# mut.strenght- avg. strenght o plant-pollinator interaction
#p_index- index value of which plant species  is in the loop of species interaction.
#points - no. of points used to approximate the normal distributions: 6 is perfect! below that is a bit shaky!
# mat, pref_matrix - both incidence matrix of which plant species interacts with which pollinator species, which in our case is all 1's
# degree.plant - all species have same degree because potentially any species can interact with anybody based on the trait values they possess and how similar they are to others.
#Pa -  phenotypic variance of plants
  type_2_plants<-function(m,sigma,w,h,np,na,mut.strength,p_index,
                          points,mat,pref_matrix,degree.plant,Pa){
    temp2<-dat2<-x2<-x3<-gvar<-j1<-array(dim=c(points))
    
    z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp[p_index], sigma =sqrt(sigma$sp[p_index]))$nodes #z of animals
    w1<-gauss.quad.prob(points, dist = "normal", 
                        mu=m$mp[p_index],sigma =sqrt(sigma$sp[p_index]))$weights #p_i(z)
    
    z2<-matrix(0,nrow=length(na),ncol=points)
    w2<-matrix(0,nrow=length(na),ncol=points)
    numer_a<-denom_a<-numer_m<-denom_m<-numer_g<-denom_g<-matrix(0,nrow=points,ncol=length(na))
    N_strength<-m_strength<-g_strength<-numeric()
  
    for(j in 1:points){  
      for(k in 1:length(na)){
        z2[k,]<-gauss.quad.prob(points, dist = "normal", mu=m$ma[k], sigma =sqrt(sigma$sa[k]))$nodes #z''
        
        #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
        w2[k,]<-gauss.quad.prob(points, dist = "normal", 
                                mu=m$ma[k],sigma =sqrt(sigma$sa[k]))$weights #pj(z'')
        
        numer_a[j,k]<- na[k]*mat[p_index,k]*pref_matrix[p_index,k]*(mut.strength/degree.plant)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        denom_a[j,k]<- na[k]*mat[p_index,k]*pref_matrix[p_index,k]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        
        numer_m[j,k]<- na[k]*mat[p_index,k]*pref_matrix[p_index,k]*(mut.strength/degree.plant)*sum((z1[j]-m$mp[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        denom_m[j,k]<- na[k]*mat[p_index,k]*pref_matrix[p_index,k]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        
        numer_g[j,k]<- na[k]*mat[p_index,k]*pref_matrix[p_index,k]*1/Pa[p_index]^2*(mut.strength/degree.plant)*sum(((z1[j]-m$mp[p_index])^2-Pa[p_index])*exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        denom_g[j,k]<- na[k]*mat[p_index,k]*pref_matrix[p_index,k]*(mut.strength)*sum(exp(-(z1[j]- z2[k,])^2/w^2)*w2[k,])
        
      }
      
      N_strength[j] <-(sum(numer_a[j,])/(1+h*sum(denom_a[j,])))*w1[j]
      m_strength[j] <-(sum(numer_m[j,])/(1+h*sum(denom_m[j,])))*w1[j]
      g_strength[j] <-(sum(numer_g[j,])/(1+h*sum(denom_g[j,])))*w1[j]
    }
    

    
    G = sum(N_strength)
    B = sum(m_strength) 
    V = sum(g_strength)
    
    
    
    return(list(G= G, B = B,V=V))
    
  }
  
  # function for dynamical eco-evolutionary dynamics of mutualistic networks; the main dynamical function to simulate eco-evo dynamics - eq 5 and 6 of the main text.
  #time : time points
  #state : state variables of plants, animals, densities and trait values
  #pars: parameter values
  eqs_type2_new<- function(time, state, pars) {
    A <- dim(pars$matrix)[2]  ## number of animal species
    P <-dim(pars$matrix)[1]
    Np<-state[(A+1):(A+P)]
    Na<-state[1:A]
    s <- pars$sig ## species' trait standard deviations
    ma<-state[(A+P+1):(A+P+A)]
    mp<-state[(A+P+A+1):(A+P+A+P)]
    opTemp<-state[(A+P+A+P+1):(A+P+A+P+1)]
    gvarA<-pars$gvarA
    gvarP<-pars$gvarP
    pgr<-numeric()
    s <- pars$sig ## species' trait standard deviations
    alpha.a<-pars$Amatrix ## alpha matrix
    alpha.p<-pars$Pmatrix
    
    aij<-bij<-vij<-matrix(0, nrow=A,ncol=P) 
    aji<-bji<-vji<-matrix(0, nrow=P,ncol=A) 
    muA<-ma
    muP<-mp
    aj<-bj<-ai<-bi<-vj<-vi<-numeric()

    varA <- ( (gvarA + pars$envA))
    varP <- ( (gvarP + pars$envP))
    
    dmA <- outer(muA, muA, FUN="-") ## difference matrix of trait means
    dmP <- outer(muP, muP, FUN="-") ## difference matrix of trait means
    
    svA <- outer(varA, varA, FUN="+") ## sum matrix of trait variances
    svP <- outer(varP, varP, FUN="+") ## sum matrix of trait variances
    
    #alpha matrix in eq 6
    alphaA <- exp(-dmA^2/(2*svA+pars$w_c^2))*pars$w_c/sqrt(2*svA+pars$w_c^2) ## alpha matrix
    alphaP <-  exp(-dmP^2/(2*svP+pars$w_c^2))*pars$w_c/sqrt(2*svP+pars$w_c^2) ## alpha matrix
    diag(alphaA)<- pars$alpha_ii #*diag(alphaA)  #intraspecific competition which was fixed
    diag(alphaP)<- pars$alpha_ii #*diag(alphaP) #intraspecific competition which was fixed
    
    
    # beta matrix in eq 6
    betaA <- alphaA*2*varA*(-dmA)/(2*svA+pars$w_c^2)^1.5 ## beta matrix
    betaP <- alphaP*2*varP*(-dmP)/(2*svP+pars$w_c^2)^1.5 ## beta matrix
    diag(betaA)<- 0 #diag(betaA)
    diag(betaP)<- 0 #diag(betaP)
 
    # a- value that captures the shape of temperature-tolerance curve.
    a<-pars$a
    #environmental temperature
    Temp<-pars$Temp  
    
    #growth rate terms i.e., b for density changes over time -  bi(T) in eq 5
    ba<- pars$gi/(pars$bw-a*muA)*(pars$bw-a*muA)/(sqrt((pars$bw-a*muA)^2+varA))*exp(-(Temp- muA)^2/(2*(pars$bw-a*muA)^2+varA)) - pars$ki
    bp<- pars$gi/(pars$bw-a*muP )*(pars$bw-a*muP)/(sqrt((pars$bw-a*muP)^2+ varP))*exp(-(Temp- muP)^2/(2*(pars$bw-a*muP)^2+varP)) - pars$ki
    
    #rate of change in growth due to tolerance to temperature and its impact on mean trait - \bar{bi} in equation in eq 6
    bar_ba<- pars$gi/(pars$bw-a*muA)*exp(-(Temp- muA)^2/(2*(pars$bw-a*muA)^2+varA))*(varA)*(pars$bw-a*muA)*(Temp -muA)/(2*(pars$bw-a*muA)^2 + varA)^1.5
    bar_bp<-pars$gi/(pars$bw-a*muP)*exp(-(Temp- muP)^2/(2*(pars$bw-a*muP)^2+varP))*(varP)*(pars$bw-a*muP)*(Temp - muP)/((pars$bw-a*muP)^2 + varP)^1.5 
    
    #looping over Animal species to caluclate the type-2 functional curve estimation of mutualistic species interaction by Gaussian quadrature
    for(r in 1:A){
        m.temp<-list(ma=muA,mp=muP)
        sigma1<-list(sa=varA,sp=varP)
        temp1<-type_2_animals(m=m.temp,sigma=sigma1,w=pars$w,h=pars$H,np=Np,na=Na,mut.strength=pars$mut.strength,a_index=r,
                       points=6,mat=pars$matrix,pref_matrix=pars$pref_matrix,degree.animal=pars$dganimals[r],Pa=varA)
      ai[r]<-temp1$G
      bi[r]<-temp1$B
      vi[r]<-temp1$V
      
    }
    #looping over plant species to caluclate the type-2 functional curve estimation of mutualistic species interaction by Gaussian quadrature
    for(k in 1:P){
        m2.temp<-list(ma=muA,mp=muP)
        sigma2<-list(sa=varA,sp=varP)
        temp2<-type_2_plants(m=m2.temp,sigma=sigma2,w=pars$w,h=pars$H,np=Np,na=Na,
                               mut.strength=pars$mut.strength,p_index=k,
                               points=6,mat=pars$matrix, pref_matrix=pars$pref_matrix,
                               degree.plant =pars$dgplants[k],
                               Pa=varP)
  
      aj[k]<-temp2$G
      bj[k]<-temp2$B
      vj[k]<-temp2$V
    }
    
    
    dndt_a<- Na*(ba-alphaA%*%Na+ai)*cutoff(Na/(1e-6)) #  na*(ba-alpha.a%*%na+ai)*cutoff(na/(1e-8)) #population dynamics
    dndt_p<- Np*(bp-alphaP%*%Np+aj)*cutoff(Np/(1e-6))  #population dynamics
    dudt_A<- pars$trait_evolve*gvarA/varA*(bar_ba- betaA%*%Na+ bi) #mean trait dynamics
    dudt_P<- pars$trait_evolve*gvarP/varP*(bar_bp -betaP%*%Np+ bj) #mean trait dynamics
    dtemp<-0.0
    ## return equations by first flattening them back into a single vector
    return(list(c(dndt_a, dndt_p,dudt_A,dudt_P,dtemp)))
  }
  
  
  
 #function to estimate rewiring over time: again used in the simulations for the Theobiota Cluster
  #input:
  #dat: data frame from the ODE solver
  #pars: list of parameters
  #output: uses functions from the bipartite package such as"betalinkr_multi() to output rewiring (OS), and turnover(ST) over time 
  rewiring_estimate<-function(dat,pars){
    nestedness_ot<-connectance_ot<-connectance_weighted<-numeric()
    timepoints<-unique(dat$time)

    matrix_time<-NULL
    
    #looping over time points
    for(t in timepoints){
    
    #for each time point get out the species densities
    density_n <-(dat %>%filter(type=="N",time == t) %>% select(species,v))$v 
    density_p <-(dat %>%filter(type=="P",time == t) %>% select(species,v))$v 
    #for each time point get out the trait values
    trait_n <-(dat %>%filter(type=="muA",time == t) %>% select(species,v))$v 
    trait_p <-(dat %>%filter(type=="muP",time == t) %>% select(species,v))$v 
    #trait variances
    va_r<-pars$gvarA + pars$envA
    vp_r<-pars$gvarP + pars$envP
    
    # differences in mean trait values of species
    d_ma <- outer(trait_n, trait_p, FUN="-") 
    #addition of mean trait variances of species
    sva<- outer(va_r,vp_r, FUN="+")
    
    #modified gowervariance matrix from the main-text
    Gowers_variance_matrix<-  abs(d_ma)/sqrt(sva) 
    Gowers_variance_matrix<-  1- pmin(Gowers_variance_matrix,1)
    
    Plant_n<-density_p
    Poll_n<-density_n
    #threshold density cut-off, any species below 1e-4 are considered not present
    Poll_n[Poll_n<1e-4]<-0
    Plant_n[Plant_n<1e-4]<-0
    # all species are either present or absent, 1 or 0 based on the density cut-off
    Plant_n[Plant_n>0]<-1
    Poll_n[Poll_n>0]<-1
    
    #estimating presence-based gower variance matrix with threshodl of 0.2
    N_mat<-outer(Poll_n,Plant_n,FUN="*")
    Gowers_variance_matrix<-Gowers_variance_matrix*N_mat
    #here 0.2 is the threshold cut-off of whether a species is interacting with others or not. If its not, interaction is 0
    Gowers_variance_matrix[which(Gowers_variance_matrix>=0.2)] <-1 
    Gowers_variance_matrix[which(Gowers_variance_matrix<0.2)] <-0
    
    matrix_time[t]<-list(Gowers_variance_matrix)
    #calculate netstedness and connectance over time based on the gowervariance matrix
    nestedness_ot[t] <- nestedness_NODF(Gowers_variance_matrix)
    connectance_ot[t] <- Connectance(Gowers_variance_matrix)
    }
    matrix_time<-matrix_time[!sapply(matrix_time,is.null)]
  
    #prepare dataframe for rewiring estimate calcuation over time
    Ls<-c(letters,LETTERS, seq(500,800,1))
    bp.pois_OS<-bp.pois_ST<-numeric()
    for(i in 1:(length(matrix_time)-1)){
      dimnames(matrix_time[[i]])<-list(Ls[1:dim(matrix_time[[1]])[1]], Ls[1:dim(matrix_time[[1]])[2]])
      dimnames(matrix_time[[i+1]])<-list(Ls[1:dim(matrix_time[[1]])[1]], Ls[1:dim(matrix_time[[1]])[2]])
      }
    
    #Rewiring OS based on Poisot's paper
    bp.pois_OS<-betalinkr_multi(webarray = webs2array(matrix_time),partitioning="poisot")
    bp.pois_ST<-betalinkr_multi(webarray = webs2array(matrix_time),partitioning="poisot")
    
    OS<-bp.pois_OS$OS[1:(length(matrix_time)-1)]
    ST<-bp.pois_OS$ST[1:(length(matrix_time)-1)]
    
    
    return(list(OS=OS,ST=ST,connectance_ot=connectance_ot,nestedness_ot=nestedness_ot))
  }
  
  
  # function to simulate the dynamics of differential equations of the networks, as well as quantify all network metrics (used for the Cluster simulations not for the main-texts figures)
  #params: parameter list
  #ic: state variables
  #tmax: max time points
  #This function was used for simulating dynamics at the Theobiota cluster in Bielefeld University.
  cluster_run_func<-function(params, ic, tmax ){
     
      OUT<-ode(func=eqs_type2_new, y=ic, parms=params,
             times=seq(0, tmax, by=tmax/1000),rtol = 1e-9, atol = 1e-12) %>% 
      organize_results_rewire(pars = params) 
    
    (Na_r<-(OUT  %>% filter(time == tmax,type=="N"))$v)
    (Np_r<-(OUT %>% filter(time == tmax,type=="P"))$v)
    mua_r<- (OUT  %>% filter(time > 0,type=="muA") %>% group_by(species) %>%   summarise(
      mean_v = mean(v, na.rm = TRUE),
      n_timepoints = n(),
      .groups = "drop"
    ))$mean_v
    mup_r<-(OUT %>% filter(time> 0,type=="muP") %>% group_by(species) %>% summarise(
      mean_v = mean(v, na.rm = TRUE),
      n_timepoints = n(),
      .groups = "drop"
    ))$mean_v
    va_r<-params$gvarA
    vp_r<-params$gvarP  
    
    Na_r[which(Na_r<0)]<-0
    Np_r[which(Np_r<0)]<-0
  
     A<-length(Na_r)
     P<-length(Np_r)
     na<-ic[1:length(Na_r)]
     np<-ic[(length(Na_r)+1) : (length(Na_r)+length(Np_r) )]
    p_i <- c(Na_r,Np_r)/sum(c(Na_r,Np_r))
    inv_simpson<-1/sum(p_i^2)
  
    N_initial<-c(na,np)   
    final_density<-c(Na_r,Np_r) #final density at tmax
    richness<-length(which(final_density>1e-4)) #final richness not needed in iur model
    g<-params$matrix #adjacency matrix all 1's
    Aspecies<- dim(g)[2] 
    Plantspecies<- dim(g)[1] 
    d_ma <- outer(mup_r, mua_r, FUN="-") ## difference matrix of trait means
    svA <- outer( (vp_r+params$envP), (va_r+params$envA), FUN="+") ## sum matrix of trait variances
    multi_va<-outer( (vp_r+params$envP), (va_r+params$envA), FUN = "*") #multiplication of variances

   network_metrics<-quantitative_network(muA = mua_r,muP = mup_r, 
                                         Na = Na_r, Np = Np_r, sa = (va_r+params$envA),
                                         sp = (vp_r+params$envP),w = params$w,thr = 0.2)
   
   # H2 index 
   H2_bhatd_int<-network_metrics$H2_bhat
   H2_gowervariance_int<-network_metrics$H2_gowervariance
   H2_true_est<-network_metrics$H2_true_est
  
   # network modularity
   modularity_bhat <- (computeModules(network_metrics$incidence_matrix_bhat))@likelihood
   modularity_gowervariance <-  (computeModules(network_metrics$web_matrix_gowervariance))@likelihood
   modularity_true_est<- (computeModules(network_metrics$incidence_matrix_true_estimate))@likelihood
 
   #connectance 
   connectance_bhat_unweighted<-network_metrics$connectance_bhat
   connectance_gowervariance_unweighted<- network_metrics$connectance_gowervariance
   connectance_true_est_unweighted<- network_metrics$connectance_true_est
  
   #nestedness NODF                                         
   nestedness_bhat_unweighted<-nestedness_NODF(network_metrics$incidence_matrix_bhat)
   nestedness_gowervariance_unweighted<-nestedness_NODF(network_metrics$web_matrix_gowervariance)
   nestedness_true_est_unweighted<- nestedness_NODF(network_metrics$incidence_matrix_true_estimate)
   
     dd<-OUT
    dat_temp<-rewiring_estimate(dd,pars=params)
    rewiring_overtime <-dat_temp$OS
    turnover_time <- dat_temp$ST
    connectance_time<-dat_temp$connectance_ot
    nestedness_time<-dat_temp$nestedness_ot

    connectance_time<-connectance_time[!is.na(connectance_time)]
    nestedness_time<-nestedness_time[!is.na(nestedness_time)]
    
    connectance_bhat_unweighted<-connectance_bhat_unweighted
    connectance_true_est_unweighted <- connectance_true_est_unweighted
    nestedness_bhat_unweighted<- nestedness_bhat_unweighted
    
    
    community_biomass<- sum(final_density)
    network_size <- (Aspecies+Plantspecies)
  
  
    species_index<-params$species_index
    
    timeseries_data<-data.frame(rewiring_overtime,
                                   turnover_time,
  				rep(as.character(params$variation), each=length(rewiring_overtime)),
                                rep(as.character(params$rewiring),each=length(rewiring_overtime)),
                                   seq(1,length(turnover_time),1),
                                   rep(as.character(params$web.name),each=length(rewiring_overtime)),
                                  rep(params$a, each=length(rewiring_overtime)),
                                rep(params$H, each=length(rewiring_overtime)),
                                rep(network_size, each=length(rewiring_overtime)),
                                   rep(params$w_c, each=length(rewiring_overtime)),
                                   rep(params$h2,each=length(rewiring_overtime)),
                                   rep(params$Temp,each=length(rewiring_overtime)),
                                   connectance_time[1:length(rewiring_overtime)],
                                   nestedness_time[1:length(rewiring_overtime)])

    colnames(timeseries_data)<-c("Rewiring","Turnover","variation","percent_rewire", "time","webname","a","Handling_time", "network_size", "competition",
                                 "h2","starting_temperature","connectance","Nestedness")

    species_level_data <-data.frame(
      rep(seq(1,(nrow(g)+ncol(g)),1)), #number of species
      c(mua_r, mup_r),
      c(va_r,vp_r),
      c(params$envA,params$envP),
      as.numeric(c(Na_r,Np_r )),
      params$mut.strength[1],
      params$ki,
      params$Temp,
      params$H,
      params$w,
      params$w_c,
      rep(as.character(params$web.name),each=(Aspecies+Plantspecies)),
      as.character(params$variation),
     as.numeric(rep(nestedness_bhat_unweighted, each=((Aspecies+Plantspecies)) )),
     as.numeric(rep(nestedness_gowervariance_unweighted, each=((Aspecies+Plantspecies)) )),
     as.numeric(rep(nestedness_true_est_unweighted, each=((Aspecies+Plantspecies)) )),
  as.numeric(rep(connectance_bhat_unweighted, each=((Aspecies+Plantspecies)) )), 
  as.numeric(rep(connectance_gowervariance_unweighted, each=((Aspecies+Plantspecies)) )), 
  as.numeric(rep(connectance_true_est_unweighted, each=((Aspecies+Plantspecies)) )), 
 
  as.numeric(rep(modularity_bhat,each=((Aspecies+Plantspecies)))),
  as.numeric(rep(modularity_gowervariance,each=((Aspecies+Plantspecies)))),
  as.numeric(rep(modularity_true_est,each=((Aspecies+Plantspecies)))), 

  as.numeric(rep(H2_bhatd_int, each=((Aspecies+Plantspecies)) )),
  as.numeric(rep(H2_gowervariance_int,each=((Aspecies+Plantspecies)) )),
as.numeric(rep(H2_true_est,each=((Aspecies+Plantspecies)))),
as.numeric(rep( (Aspecies+Plantspecies),each=((Aspecies+Plantspecies)))),
      as.character(params$h2),
      as.numeric(rep( community_biomass,each=((Aspecies+Plantspecies)))),
  as.numeric(rep(richness,each=((Aspecies+Plantspecies)))),
  as.numeric(rep(inv_simpson,each=((Aspecies+Plantspecies)))),
  
    as.numeric(rep( params$a,each=((Aspecies+Plantspecies)))))
    # 
     colnames(species_level_data)<-c("Species",
                                     "trait",
                                     "variance",
                                     "environmental_variance",
                                     "density", 
                                     "mutualism_strength", 
                                     "mortality",
                                     "Temperature_shift",
                                     "Handling_time",
                                     "width_mutualism",
                                    "width_competition",
                                     "webname",
                                     "variation",
                                    "nestedness_bhat_unweighted",
                                    "nestedness_gowervariance_unweighted",
                                    "nestedness_true_estimate",
                                    "Connectance_bhat_unweighted",
                                    "Connectance_gowervariance_unweighted",
                                    "Connectance_true_estimate",
				                            "modularity_bhat",
				                           "modularity_gowervariance",
				                           "modularity_true_est",
                                    "H2_bhat_int",
                                    "H2_gowervariance_int",
                                    "H2_true_estimate",
                                     "Network_size",
                                     "h2",
                                      "biomass",
                                    "richness",
                                    "inverse_simpson",
                                    "a")
     
    # final_output<- list(species_level_data=species_level_data,timeseries_data=timeseries_data)
     
     final_output<- list(timeseries_data=timeseries_data,species_level_data=species_level_data)
    return(final_output)
    
  }
  
  #function to estimate network metrics from trait overlaps 
  #muA, muP: trait values
  #Na,Np: densities at final timepoint
  #sa,sp: trait variances
  #w : width of mutualistic interaction
  #thr: threshold interaction
  quantitative_network<-function(muA,muP,Na,Np,sa,sp,w,thr){
    
    Plant_n<-Np
    Poll_n<-Na
    Poll_n[Poll_n<0]<-0
    Plant_n[Plant_n<0]<-0
    Poll_n[Poll_n<1e-5]<-0
    Plant_n[Plant_n<1e-5]<-0
    Poll_n[Poll_n>=1e-5]<-1
    Plant_n[Plant_n>=1e-5]<-1
    
    
    N_mat<-outer(Plant_n,Poll_n,"*")
    
    d_ma <- outer(muP, muA, FUN="-") ## difference matrix of trait means
    svA <- outer( sp, sa, FUN="+") ## sum matrix of trait variances
    multi_va<-outer(sp, sa, FUN = "*") #multiplication of variances
    
    #gower matrix
    
    gower_matrix<-1- abs(d_ma)/(max(c(muP,muA))-min(c(muP,muA)))
    adjacency_gower<-gower_matrix*N_mat
    adjacency_gower[which(adjacency_gower <=thr)]<-0
    adjacency_gower[which(adjacency_gower >thr)]<-1
    web_matrix_gower<-adjacency_gower
    
    #gowers trait variance weighted
    gowers_matrix_variance <- abs(d_ma)/sqrt(svA)
    gowers_matrix_variance <- 1- pmin(gowers_matrix_variance,1)
    
    incidence_matrix2 <- gowers_matrix_variance*N_mat
    incidence_matrix2[which(incidence_matrix2<= thr)] <- 0
    incidence_matrix2[which(incidence_matrix2> thr)] <- 1
    web_matrix_gowervariance<-incidence_matrix2
    

    standardised_mp1<- (muP-min(muP))/(max(muP)-min(muP))
    standardised_ma1<- (muA-min(muA))/(max(muA)-min(muA))
    standardised_d_ma <- outer(standardised_mp1, standardised_ma1, FUN="-") ## difference matrix of trait means
    standardised_va<-outer(sp, sa, FUN = "+")
    t_com_mat<-abs(standardised_d_ma) < 0.5*((standardised_va))
    t_com_mat[which(t_com_mat == TRUE)]<-1
    t_com_mat[which(t_com_mat == FALSE)] <-0
    web_matrix_complimentary<-t_com_mat*N_mat
    
    #BC coefficient
    bhat_distance <- 1/sqrt((svA)/(2*sqrt(multi_va)))*exp(-d_ma^2/(4*svA))
    bhat_distance <- (bhat_distance - min(bhat_distance))/(max(bhat_distance)-min(bhat_distance))
    #bhat_distance <-pmax(bhat_distance,0.05)
    incidence_matrix_bhat<-bhat_distance*N_mat
    incidence_matrix_bhat[which(incidence_matrix_bhat>thr)]<-1
    incidence_matrix_bhat[which(incidence_matrix_bhat<=thr)]<-0   
    
    #Actual estimate 
    true_estimate <- w/sqrt(svA+w)*exp(-d_ma^2/(svA+w))
    true_estimate<- (true_estimate - min(true_estimate))/(max(true_estimate)-min(true_estimate))
    incidence_matrix_true_estimate<-true_estimate*N_mat
    incidence_matrix_true_estimate[which(incidence_matrix_true_estimate>thr)]<-1
    incidence_matrix_true_estimate[which(incidence_matrix_true_estimate<=thr)]<-0   
    
    
    H2_complimentary <- (bipartite::H2fun((t_com_mat*N_mat),H2_integer = FALSE))[1]
    H2_bhat <- (bipartite::H2fun((incidence_matrix_bhat),H2_integer = FALSE))[1]
    H2_gowervariance <-(bipartite::H2fun((incidence_matrix2),H2_integer = FALSE))[1]
    H2_gower<-(bipartite::H2fun((web_matrix_gower),H2_integer = FALSE))[1]
    H2_true_est<-(bipartite::H2fun((incidence_matrix_true_estimate),H2_integer = FALSE))[1]
    
    connectance_complimentary <-Connectance(t_com_mat*N_mat)
    connectance_bhat<-Connectance(incidence_matrix_bhat)
    connectance_gowervariance<-Connectance(incidence_matrix2)
    connectance_gower<-Connectance(web_matrix_gower)
    connectance_true_est<-Connectance(incidence_matrix_true_estimate)
    
    return(list(connectance_complimentary=connectance_complimentary,connectance_gowervariance=connectance_gowervariance,connectance_gower=connectance_gower,
                connectance_bhat=connectance_bhat,H2_bhat=H2_bhat,H2_complimentary=H2_complimentary,H2_gowervariance=H2_gowervariance,H2_gower=H2_gower,
                incidence_matrix_bhat=incidence_matrix_bhat,web_matrix_gowervariance=web_matrix_gowervariance,
                web_matrix_complimentary=web_matrix_complimentary,web_matrix_gower=web_matrix_gower,
                incidence_matrix_true_estimate=incidence_matrix_true_estimate,connectance_true_est=connectance_true_est,
                H2_true_est=H2_true_est))
  }
  
  
  #organise the results from the output of ODE solver. It takes the output from the ode solver and a list of parameters
  #to organise them in a nice tibble
  #inputs are : sol -  solution from the ODE solver funciton
  # pars - list of parameters
  organize_results_rewire <- function(sol, pars) {
    S <- length(pars$sig) ## number of species
    A<-dim(pars$matrix)[2] # no. of animals
    P<-dim(pars$matrix)[1] # no. of plants
    temp<- sol %>% as.data.frame %>% as_tibble ## convert to tibble
    ## name the first column "time"
    temp<- temp 
    names(temp)[1] <- "time"
    names(temp)[2:(A+1)] <- paste0("N_", 1:(A)) ## name abundance columns (n_k)
    
    names(temp)[(A+2):(A+1+P)] <- paste0("P_", 1:P) ## name trait mean columns
    names(temp)[(A+P+2):(A+P+A+1)] <-  paste0("muA_", 1:(A))
    names(temp)[(A+P+A+2):(A+P+A+P+1)] <- paste0("muP_", 1:(P))
    names(temp)[(A+P+A+P+2):(A+P+A+P+2)] <- paste0("Op_", 1:1)
    
    temp <- temp %>%
      tidyr::gather("variable", "v", 2:ncol(temp)) %>% ## normalize the data
      tidyr::separate(variable, c("type", "species"), sep="_") %>%
      #spread(type, v) %>% ## separate columns for animal densities n and plant densities m
      dplyr::select(time, type, species,v) %>% ## rearrange columns
      mutate(species=as.integer(species), w=pars$gamma,
             Nestedness=pars$nestedness, Connectance=pars$C,
             theta=pars$theta,Web.name=pars$web.name) ## add params
    return(as_tibble(temp))
  }
  
  
  
  ## Plot time series of densities, time series of trait values, and
  ## snapshot of the trait distributions at time = moment
  ## - dat: data generated by organize_results()
  ## - moment: time at which trait distribution should be plotted
  ## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
  ## - res: number of evenly spaced sampling points along the trait axis
  ##               for the trait distribution plot
  ## Output:
  ## - a ggplot2 plot with three panels in one column: abundance time series,
  ##   trait value time seties, and snapshot of trait distribution
  plot_all <- function(dat, moment=0, limits=c(-1, 1), res=1001) {
    plot_grid(plot_density(dat), ncol=1, align="hv") %>%
      return
  }
  
  
  
  ## Plot species densities through time
  ## Input:
  ## - dat: data generated by organize_results()
  ## Output:
  ## - a ggplot2 plot
  ## used to produce figure 2.
  plot_density<- function(dat) {
    dat %>%
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=1.25) +
      scale_y_continuous(name="population density") +
      theme_classic()+
      scale_colour_viridis_d(option = "plasma")+
      theme(legend.position="none") + facet_wrap(.~type,scales = "free") %>%
      return
  }
  
  ## Plot species densities through time
  ## Input:
  ## - dat: data generated by organize_results()
  ## Output:
  ## - a ggplot2 plot
  ## used to produce figure 2.
  plot_density_timeseries<- function(dat) {
    dat %>% filter(type == c("N","P")) %>% 
      ggplot +
      geom_line(aes(x=time, y=v, colour = factor(species)),size=1.25) +
      scale_y_continuous(name="population density") +
      theme_cowplot()+
      theme(legend.position="none") + 
      facet_wrap(.~type,scales = "free") %>%
      return
  }
  

  
  #function to plot trait density distribution over various time snapshots, used to produce figure 2
  # Na,Np- densities of animals and plants
  #m - vector of trait values for that particular time point
  #sigma - vector of trait variances of all species
  # Temp -  environmental temperature
  #res - resolution
  #returns a figure 
  plot_snapshot <- function(Na, Np, m, sigma, Temp, limits=c(-1, 1), res=1001) {
    S_a <-length(Na) ## number of species
    S_p <- length(Np)
    ma<- m[1:(S_a)]
    mp<- m[(S_a+1):(S_a+S_p)]
    sigma_a <-sigma[1:(S_a)]
    sigma_p <- sigma[(S_a+1):(S_a+S_p)]
    traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
    #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
    traits_a <- expand.grid(species=1:S_a, trait=traitaxis) %>% as_tibble ## trait table
    traits_p <- expand.grid(species=(S_a+1):(S_a+S_p), trait=traitaxis) %>% as_tibble ## trait table
    
    traits_a["density"] <- 0 ## add column for population densities
    traits_p["density"] <- 0
    
    for (i in 1:S_a) {
      #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
      traits_a$density[(traits_a$species==i)] <- Na[i]*
        dnorm(traits_a$trait[(traits_a$species==i)], ma[i], (sigma_a[i])) ## times density
    }
    traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
    
    for (i in 1:S_p) {
      #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
      traits_p$density[(traits_p$species==((S_a)+i))] <- Np[i]*dnorm(traits_p$trait[(traits_p$species==((S_a)+i))], 
                                                                       mp[i], (sigma_p[i])) ## times density
    }
    traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
    
    
    
    landscape <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
      mutate(r= (1.5/(4.5-0.12*(trait))*exp(-(Temp-trait)^2/(2*(4.5-0.12*(trait))^2))-0.1)) %>%      
      mutate(r=ifelse(r<=0, NA, r)) %>%
      mutate(r=r*max(c(traits_a$density,traits_p$density),na.rm = T))
    
 
    traits<-data.frame(rbind(traits_a,traits_p), 
                       species_group=c(rep("Animals", nrow(traits_a)),
                                       rep("Plants", nrow(traits_p))))
    n_animals <- length(unique(traits$species[traits$species_group == "Animals"]))
    n_plants  <- length(unique(traits$species[traits$species_group == "Plants"]))
    
    # generate gradient palettes
    orange_pal <- gradient_n_pal(c("#FFE5B4", "#FF8C00"))(seq(0,1,length.out=n_animals))
    purple_pal <- gradient_n_pal(c("#F4A6A6", "#b2182b"))(seq(0, 1, length.out = n_plants))  # light red → firebrick

    # assign colors by species
    species_levels <- sort(unique(traits$species))
    species_groups <- traits$species_group[match(species_levels, traits$species)]
    species_colors <- c(
      setNames(orange_pal, species_levels[species_groups == "Animals"]),
      setNames(purple_pal, species_levels[species_groups == "Plants"])
    )
    
    
    ggplot(traits) + ## generate plot
      geom_line(aes(x=trait, y=density, colour=factor(species)), alpha=1, na.rm=TRUE) +
      geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                  alpha=0.5, colour=NA)+
      facet_wrap(.~species_group, nrow = 2)+
      theme_classic()+
      theme(legend.title = element_text(size = 20), 
            legend.position = "right", panel.background = element_blank(), 
            axis.text = element_text(colour = "black", size = 20), 
            axis.title = element_text(size = 20), 
            legend.text = element_text(size =20 ), legend.key = element_blank(),
            strip.text.x = element_text(size= 20),
            panel.spacing = unit(1.5, "lines"))+
      geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
                colour="black", alpha=1, na.rm=TRUE) +
      scale_x_continuous(name=expression("Temperature (" ~ degree*C ~ ")"), limits=limits) +
      scale_y_continuous(name="Trait density", limits=c(0, NA)) + 
      scale_fill_manual(values = species_colors) +
      scale_color_manual(values = species_colors)+
      theme(legend.position="none") %>%
      return 
  }
  
  
 
  #computes the raw NODF taken from Song et al 2017 J. Animal Ecology
  #input: web = mutualistic network
  #output: raw NODF of the given network
  nestedness_NODF <- function(web){
    web[web > 0] = 1
    SA <- nrow(web) #no of animal species
    SP <- ncol(web) #no of plant species
    N <- t(web) %*% web
    num <- N
    num[lower.tri(num,diag=TRUE)]=1
    den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
    dele <- den - t(den)
    dele[lower.tri(dele,diag=TRUE)] <- 1
    num[dele == 0] <- 0
    den <- pmin(den,t(den))
    den[lower.tri(den,diag=TRUE)] = 1
    nes <- num/den
    nes[lower.tri(nes,diag=TRUE)] = 0
    nes[is.na(nes)] <- 0
    n1 <- sum(nes)
    
    N <- web %*% t(web)
    num <- N
    num[lower.tri(num,diag=TRUE)]=1
    den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
    dele <- den - t(den)
    dele[lower.tri(dele,diag=TRUE)] <- 1
    num[dele ==0 ] <- 0
    den <- pmin(den,t(den))
    den[lower.tri(den,diag=TRUE)]=1
    nes <- num/den
    nes[lower.tri(nes,diag=TRUE)] = 0
    nes[is.na(nes)] <- 0
    n2 <- sum(nes)
    out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
    return(out)
  }
  
  
  
  
  # function for competition coefficients within a guild.
  # competitive interactions  were scaled by the total number of species within a guild as Dakos & Bascompte 2014 PNAS.
  # matrix: network of interactions which are 0 or 1. 
  # strength: average competition str
  # measures connectance of a web network
  Connectance<-function(web)
  {
    return(sum(web)/(ncol(web)*nrow(web)))}
  
  
  

