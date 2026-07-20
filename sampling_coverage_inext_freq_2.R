# R script with metadata of sampling coverage of 118 networks used to investigate
# how sampling effort varied across latitudes and whether it impacted
#network architecture such as connectance, nestedeness, H2, and modularity


#uses the metadata extracted from original sources of 118 networks which had sampling frequencies
temp <- read.csv("subset_network_data_sampling_coverage.csv",sep=",")

head(temp)
plot(temp$wc2.1_2.5m_bio_1, temp$coverage)

#regressions with sampling coverage of 118 networks
summary(lm(connectance  ~ wc2.1_2.5m_bio_1 + coverage , data = temp))
summary(lm(NODF      ~ wc2.1_2.5m_bio_1 + coverage, data = temp))
summary(lm(modularity  ~ wc2.1_2.5m_bio_1 + coverage , data = temp))
summary(lm(H2          ~ wc2.1_2.5m_bio_1 + coverage, data = temp))
