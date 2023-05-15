# Simulation Study: Efficiently Selecting Representative Audit Samples
This repository contains all the necessary files to replicate the simulation study for the manuscript “Efficiently Selecting Representative Audit Samples”. The goal of the simulation study is to demonstrate application and investigate performance of the sample size approach, a method for efficient audit sample selection elaborately introduced in the manuscript.

Ethical approval for the simulation study was granted by the FETC (FETC-22-1861). 

The data used for this study are generated within the simulation study scripts. Because all scripts are seeded, one can generate the exact same data as was used to execute the study. By executing simulation_part1.R, and then simulation_part2.R, the results presented in the manuscript can be exactly replicated. Here follows an overview of all files and their purpose in the research archive:

|Files/Folders|	Description|
|-------------|------------|
|simulation_part1.R|	R-script to execute the first part of the simulation study. In the first part, simulation conditions are created and the relation between bias and deviance is determined for every condition.|
|simulation_part2.R|	R-script to execute the second part of the simulation study. In the second part, we generate data sets for every condition, and apply the sample size approach to these data sets (using the relation between bias and deviance obtained from the first part of the simulation study)  |
|process_results.R| R-script to extract the relevant features from the results of the entire simulation study, and generate the tables and plots presented in the manuscript. |
|/functions|	Folder containing the functions used throughout the simulation study. |
|/fictonal_example|	Folder containing the R-scripts for the fictional example presented in the manuscript. We provide an example data set (example_data.csv), a script to determine the deviance threshold for this example data set, and a script to subsequently apply the sample size to the example data set. It is important to emphasize that the provided example data is made up, purely for illustrative purposes. |

Software requirements are [R](http://www.r-project.org}) and the following R-packages:
- [‘nloptr’](https://cran.r-project.org/web/packages/nloptr/index.html) (version 2.0.3) was used to solve the constrained optimization model.
- ['hitandrun'](https://cran.r-project.org/web/packages/hitandrun/index.html) (version 0.5-6) was used to generate random probabilities.
- ['extraDistr'](https://cran.r-project.org/web/packages/extraDistr/index.html) (version 1.9.1) was used to sample from a multivariate hypergeometric distribution. 
- ['ggplot2'](https://cran.r-project.org/web/packages/ggplot2/index.html) (version 3.4.0) was used for plotting.
- ['dplyr'](https://cran.r-project.org/web/packages/dplyr/index.html) (version 1.0.10) was used for data manipulation. 
- ['ggpubr'](https://cran.r-project.org/web/packages/ggpubr/index.html) (version 0.5.0) was used for getting the plots in a grid. 

The entire research archive is pubicly available on Github, under responsibility of Lotte Mensink. For any help with the files in this archive, please contact lottemensink1998@gmail.com. 
