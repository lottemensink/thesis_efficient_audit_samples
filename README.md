# Simulation Study: Efficiently Selecting Representative Audit Samples
This repository contains all the necessary files to replicate the simulation study for the manuscript “Efficiently Selecting Representative Audit Samples”.

|Files/Folders|	Description|
|-------------|------------|
|simulation_part1.R|	R-script to execute the first part of the simulation study. In the first part, simulation conditions are created and the relation between bias and deviance is determined for every condition.|
|simulation_part2.R|	R-script to generate data sets for every condition, and apply the sample size approach to these data sets (using the relation between bias and deviance obtained from the first step)  |
|process_results.R| R-script to extract the relevant features from the results, and generate the tables and plots presented in the manuscript. |
|/functions|	Folder containing the functions used throughout the R scripts. |
|/fictonal_example|	Folder containing the R-scripts for the fictional example presented in the manuscript. We provide an example data set, a script to determine the deviance threshold for this example script, and a script to apply the sample size to the example data set. |

Software requirements are [R](http://www.r-project.org}) and the following R-packages:
- [‘nloptr’](https://cran.r-project.org/web/packages/nloptr/index.html) (version 2.0.3) was used to solve the constrained optimization model.
- ['hitandrun'](https://cran.r-project.org/web/packages/hitandrun/index.html) (version 0.5-6) was used to generate random probabilities.
- ['ggplot2'](https://cran.r-project.org/web/packages/ggplot2/index.html) (version 3.4.0) was used for plotting.
- ['dplyr'](https://cran.r-project.org/web/packages/dplyr/index.html) (version 1.0.10) was used for data manipulation. 

For any help with the files in this archive, please contact Lotte Mensink (lottemensink1998@gmail.com)
