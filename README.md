# Github repository for the paper titled " Adaptive rewiring and temperature tolerance shapes the architecture of plant-pollinator networks globally"
### Authors: Gaurav Baruah, Meike J. Wittmann
### Affiliations: Theoretical Biology, Faculty of Biology, University of Bielefeld, Germany.

##### Operating System: Ubuntu 22.04.5 LTS or Windows 11 ,MacOS 10.13
##### R version: R 4.5.2 (2025-10-31)

##### Required R packages:
 1. deSolve (integrating differential equations),statmod (statistical analysis).
 2. cowplot (visualisation graphics),dplyr (efficient data structuring),readr (data analysis), beepr (data arranging), ggplot2 (plotting graphics), viridis(colors for plotting),igraph (network graphics), bipartite (network analysis),sna, network, GGally (network plotting)

##### Required non-standard hardware: NONE.
On a standard desktop system, installation is immediate if R is present; otherwise, it depends on the time required to install R and its dependencies.

The repository below includes all project files, as well as the `README.md` and the GNU General Public License v3.0 in `copyingLicense.txt`.

## Empirical and model-simulated data
1. Information of the plant-pollinator network data were compiled from www.web-of-life.es database, which is an open-acces data based on species interactions. And the rest of the data was compiled from Marjakangas et al 2025, Eco. Letts. paper.
2. `coordinates.csv` contains the latitude and longitude coordinates for all the plant-pollinator networks compiled. This was used to create a map of the locations of where these networks were found, and as well to extract climate data. Used to generate Fig. 5 in the main-text and used in the R script named `Figure_5.R`.
3. `points_with_climate_data.csv`- is a csv file containing all the climatic variables extracted from the locations of the plant-pollinator networks, only mean annual temperature from this .csv file was used. Alongwith it , the .csv data file contains information of network nestedness, connectance, latitude,longitude, specialisation etc, which was also extracted from the information of plant-pollinator networks compiled. Used in the R script `Figure_5.R`.
4. `rewiring_data_20species.RData` and `network_metrics_20species.RData` - is an example data from the model simulated data of artificial plant-pollinator network of 20 species. These two data were used to produce figure 3 and figure 4 in main-text, and used in the R script `Figure_3-4.R`.
5. `subset_network_data_sampling_coverage.csv` -  is a metadata and a csv file containing climatic variables of 118 plant-pollinator networks from the original database which had interaction sampling frequencies, and compiled to have sampling coverage estiamted for each of these 118 networks. This file is used to get the table S2 in the supplementary information, and used in `sampling_coverage_inext_freq_2.R`. 

## R scripts
1. `01_functions_rewiring.R`- is an R script containing all the functions that are either used to simulate the eco-evolutionary dynamics of plant-pollinator networks or used to plot figures in the main-text, or used to estimate network metrics such as connectance, nestedness, etc. Specifically used in producing figure 2 in main-text and in R script `Figure_2.R`.
2. `Figure_2.R`-  the R script to produce figure 2 of the main-text.
3. `Figure_3-4.R` and `Figure_5.R` - are R scripts to produce figure 3, figure 4, and figure 5 of the main-text, and some other in the supplementary appendix.
4. `sampling_coverage_inext_freq_2.R` - R script for evaluating regressions of how sampling coverage and sampling effort impacts network architecture, for table S2.
5. `figure_S20_bio4_range.R` -  R script for estimating temperature variability and its relationship with network architecture for figure S20.
6. `figure_S19_BIO10_analysis.R` - R script for estimating quaterly active temperature and its relationship with network architecture for figure S19.

The R scripts and functions were used in Ubuntu (linux-gnu), and R version was 4.5.1.
