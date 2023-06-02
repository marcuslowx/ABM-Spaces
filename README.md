# ABM-Spaces
Agent-based modelling of infectious disease using the ABM Spaces model

What is it?
ABM Spaces is an agent-based model used to model infectious disease interventions. The term ‘spaces’ refers to the fact that all agents in the model are associated with settings  such as specific households or workplaces. Since transmission risk and social contact is mediated through these ‘spaces’, the model is well-suited to the modelling of airborne infections such as tuberculosis. It is not an appropriate model for infections that are more directly transmitted from person to person, such as HIV.
ABM Spaces is coded in C++. Some R packages and Makefiles are included to assist with analysis of model outputs.
There are three sub-directories, each containing a different variant of ABM Spaces:
•	ABM-Spaces-TB/ contains a version of the model used to explore the impact of five different case-finding strategies on tuberculosis incidence and mortality.
•	ABM-Spaces-Tests/ contains a version of the model used to explore the impact of variable test sensitivity and test frequency on TB incidence and mortality.
•	ABM-Spaces-COVID/ contains a version of the model used to explore the impact of test turnaround times and contact tracing on COVID-19 incidence and mortality. An earlier version of this work has been published as a preprint at https://www.medrxiv.org/content/10.1101/2020.10.06.20207761v1.
 More details on how to use each version of the model will be shared in readme files inside each directory.
ABM Spaces was developed by Marcus Low – with support from Dr Nathan Geffen – as part of studies at the Department of Computer Science at the University of Cape Town.

