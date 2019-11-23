# Bioecon_Fisheries.
Code used in "Integrating Economic dynamics into Ecological Networks: The case of fishery sustainability"
Example Code Tutorial:
Description of the example code used in the set up and simulation of each treatment type.
Valentin Cocco (valentin.cocco@ens.fr) and Paul Glaum (prglaum@umich.edu)

1)	webgeneration.m: Generation of the network structures to be used in simulations.

•	Inputs:

-	n: number of networks wanted
-	spe: number of trophic species
-	con: connectance
-	err: authorized error with connectance

•	Function called: NicheModel.m

•	Outputs: the n foodwebs created are saved in two different ways:

-	in the WebStock folder: each foodweb is described by a .csv file (WebXXX.csv)

2.) fishingfreewebsdriver.m: Using the generated food-web networks from the webgeneration.m script, complete the fishing free initial period. In this example code, a random web is chosen to simulate.

•	Inputs: 
o	Webs1100.mat: Matlab .mat file containing the 1100 webs generated and used in this study. 
o	Spe: number of species. Keep this the same unless the user generates larger webs first

•	Functions called:
o	Webproperties.m: takes in adjacency matrix of the food web and produces an array of 18 different web properties and the trophic level of each species. 
o	Setup.m: takes in the adjacency matrix of the food web, the trophic levels of each speicis, and the logical array delineating fish species. Returns biological parameters of the ATN model; economic parameters of fishing (included the choice of the harvested species); parameters of the numerical integration
o	Differential.m: Takes in time span, initial biomasses, and parameters. Runs the simulation for tspan returning derivatives of biomass (and effort in Open Access) according to ATN model. 

•	Output:
o	X: matrix of time series of each species’ biomass across the simulation. 

3.) fixedEffortwebdriverExample.m: Using the “conserved” webs from the fishing free simulations, a random conserved web is chosen and in that web a random fish is chosen for harvesting. It is fished at one of the Effort levels used in the Fixed Effort treatments. 

•	Inputs: 
o	SimCons.mat: Matlab .mat file containing the “conserved” webs (see Methods) used in this study. These webs are used in the fishing treatments. 

•	Functions called:
o	Webproperties.m: takes in adjacency matrix of the food web and produces an array of 18 different web properties and the trophic level of each species. 
o	Setup.m: takes in the adjacency matrix of the food web, the trophic levels of each speicis, and the logical array delineating fish species. Returns biological parameters of the ATN model; economic parameters of fishing (included the choice of the harvested species); parameters of the numerical integration
o	Differential.m: Takes in time span, initial biomasses, and parameters. Runs the simulation for tspan returning derivatives of biomass (and effort in Open Access) according to ATN model. 

•	Output:
o	X: matrix of time series of each species’ biomass and the fishing effort level across the simulation. Fishing Effort does not change in the Fixed Effort simulations. 

4.) openAccesswebdriverExample.m: Using the “conserved” webs from the fishing free simulations, a random conserved web is chosen and in that web a random fish is chosen for harvesting. Under the Open Access fishing treatment, the fish is harvested under the linear pricing model with an initial Effort level of 1. 

•	Inputs: 
o	SimConsLin.mat: Matlab .mat file containing the “conserved” webs (see Methods) used in this study. These webs are used in the fishing treatments. These are the same webs in SimCons.mat but SimConsLin.mat has additional categories for any data stored regarding the variable Effort level.

•	Functions called:
o	Webproperties.m: takes in adjacency matrix of the food web and produces an array of 18 different web properties and the trophic level of each species. 
o	Setup.m: takes in the adjacency matrix of the food web, the trophic levels of each speicis, and the logical array delineating fish species. Returns biological parameters of the ATN model; economic parameters of fishing (included the choice of the harvested species); parameters of the numerical integration
o	Differential.m: Takes in time span, initial biomasses, and parameters. Runs the simulation for tspan returning derivatives of biomass (and effort in Open Access) according to ATN model. 

•	Output:
o	X: matrix of time series of each species’ biomass and the fishing effort level across the simulation. 

