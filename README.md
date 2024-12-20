# ABM-Spaces
Agent-based modelling of infectious disease using the ABM Spaces model

## What is it?
ABM Spaces is an agent-based model used to model infectious disease interventions. The term ‘spaces’ refers to the fact that all agents in the model are associated with settings  such as specific households or workplaces. Since transmission risk and social contact is mediated through these ‘spaces’, the model is well-suited to the modelling of airborne infections such as tuberculosis. It is not an appropriate model for infections that are more directly transmitted from person to person, such as HIV.
ABM Spaces is coded in C++. Some R packages and Makefiles are included to assist with analysis of model outputs.
There are three sub-directories, each containing a different variant of ABM Spaces:
•	ABM-Spaces-TB/ contains a version of the model used to explore the impact of five different case-finding strategies on tuberculosis incidence and mortality.
•	ABM-Spaces-Tests/ contains a version of the model used to explore the impact of variable test sensitivity and test frequency on TB incidence and mortality.
•	ABM-Spaces-COVID/ contains a version of the model used to explore the impact of test turnaround times and contact tracing on COVID-19 incidence and mortality.

## How to use it
In each of the three directories there is a C++ file with name starting "ABM" and extention .cc. This file contains the entire model. At the top of each of these files there are instructions on how to compile and how to run the model with various witches. For example, when compiling to a.out you can run "a.out num_runs=10" to run the model 10 times. At the core of the model is a loop of timesteps. You can modify the behaviour of the model by adding and removing functions from this loop. 
The model writes output to the terminal in CSV format. To capture the out put in a CSV file you can use the ">" operator. You could for example run "a.out num_runs=10 > a.csv" to capture the output of 10 runs of the model in a file called "a.csv".

## Makefiles and R scripts
Each of the three directories contains a Makefile with several options. You can open the Makefile in a text editor to see the various options. Some options such as "Make csv" simply runs the model and then renders the output as an HTML table for quick testing. Other options run the model many times with different model configurations and then generates HTML output with key statistics and graphs. (Note that after running some of these Makefile options there will be new CSV, HTML and in some cases PNG files generated in your working directory.) You will need to have Make installed for this to work.
Most of the Makefile options make use of R and Rmarkdown to do the analysis and render output. You will thus need to have R installed for it to work plus the R packages, Rmarkdown, dplyr, htmlTable, and ggplot2. Rmarkdown also requires that you have Pandoc installed. We are confident that the various Rmarkdown files produce accurate analysis, but tightening up the coding and better annotation is on our to do list.

## Who is behind ABM Spaces?
ABM Spaces was developed by Dr Marcus Low – with support from Dr Nathan Geffen. It was developed as part of Low's PhD in the Department of Computer Science at the University of Cape Town. The PhD was supervised by Professor Michelle Kuttel and co-supervised by Geffen. The thesis, titled 'A complex, high-performance agent-based model used to explore tuberculosis and COVID-19 case-finding interventions', is online here: https://open.uct.ac.za/server/api/core/bitstreams/dc1b12de-c9e2-4b10-ae1a-d9b5d096b60c/content
