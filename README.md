# Github repository for the paper titled " Adaptive rewiring and temperature tolerance shapes the architecture of plant-pollinator networks globally"
## Authors: Gaurav Baruah, Meike J. Wittmann
## Affiliations: Theoretical Biology, Faculty of Biology, University of Bielefeld, Germany.


## Empirical and model-simulated data
1. Information of the plant-pollinator network data were compiled from www.web-of-life.es database, which is an open-acces data based on species interactions. And the rest of the data was compiled from Marjakangas et al 2025, Eco. Letts. paper.
2. `coordinates.csv` contains the latitude and longitude coordinates for all the plant-pollinator networks compiled. This was used to create a map of the locations of where these networks were found, and as well to extract climate data. Used to generate Fig. 5 in the main-text and used in the R script named `Figure_5.R`.
3. `points_with_climate_data.csv`- is a csv file containing all the climatic variables extracted from the locations of the plant-pollinator networks, only mean annual temperature from this .csv file was used. Alongwith it , the .csv data file contains information of network nestedness, connectance, latitude,longitude, specialisation etc, which was also extracted from the information of plant-pollinator networks compiled. Used in the R script `Figure_5.R`.
4. `rewiring_data_20species.RData` and `network_metrics_20species.RData` - is an example data from the model simulated data of artificial plant-pollinator network of 20 species. These two data were used to produce figure 3 and figure 4 in main-text, and used in the R script `Figure_3-4.R`.

## R scripts
1. `01_functions_rewiring.R`- is an R script containing all the functions that are either used to simulate the eco-evolutionary dynamics of plant-pollinator networks or used to plot figures in the main-text, or used to estimate network metrics such as connectance, nestedness, etc. Specifically used in producing figure 2 in main-text and in R script `Figure_2.R`.
2. `Figure_2.R`-  the R script to produce figure 2 of the main-text.
3. `Figure_3-4.R` and `Figure_5.R` - are R scripts to produce figure 3, figure 4, and figure 5 of the main-text, and some other in the supplementary appendix.

The R scripts and functions were used in Ubuntu (linux-gnu), and R version was 4.5.1.
