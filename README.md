University Lake Reservoir Gradient Project
==========

This repository contains open-source code, data, & text files for the Reservoir Gradient Project.
This project is the product of the 2013 SSRP work done with N. Nelson at Indiana University.

## Project Goals

* **Aim 1)** Determine if there is an underlying stream to lake gradient in small reservoirs.

* **Aim 2)** Determine if microbial communities change along a stream to lake gradient.

* **Aim 3)** Determine if microbial activity changes along a stream to lake gradient.

### Repo Contents


* **data:**
	* *res.grad.metab.txt*: Metabolism data (bacterial production and respiration) along the gradient.
	* *ResGrad_EnvDat.csv*: Environmental data - Contains data collected from the hydrolab and calculated environmental data.
	This file includes the following:
		1. UTM lat/long; WGS84
		2. distance to dam (m)
		3. water temperature (°C)
		4. Specific conductivity (SpC; ms/cm)
		5. Dissolved Oxygen (DO; mg/L)
		6. pH
		7. Total Phosphorus (TP; µg P/L)
		8. Respiration Rate (µM O2 Hr-1)
		9. Chlorphyll a ()
	* *ul_resgrad.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.tax.summary*: Taxonomy summary from mothur
	* *ul_resgrad.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy*: Taxonomy from mothur
	* *ul_resgrad.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared*: Site-by-OTU matrix from mothur
	* *ul_resgrad.tree*: Phylogeny 
	* *UL.design.txt*: Design corresponding to site by OTU matrix


* **bin:**
	* *mothur_tools.R*: An R script written by M. Muscarella and modified by N. Wisnoski containing functions used in the analysis of mothur raw files.
	* *Chao_functions.R*: R functions from Chiu and Chao to estimate richness in the presence of sequencing errors

* **figures:**
  	* Includes a number of intermediate and final figures included in the manuscript. 

* **intermediate-data:**
	* Includes intermediate output used in the main analysis file.

* **maps:**
  	* Includes shapefiles of university lake

* **mothur:**
	* batch file and list of raw sequence files for sequence processing

* Miscellaneous cache and output files from knitr processing of the main `ReservoirGradient.Rmd` analysis file.
* `ReservoirGradient.Rmd`: main analysis file. 

## Contributors

[Nathan Wisnoski](https://nwisnoski.github.io/).

[Mario Muscarella](http://mmuscarella.github.io/).

[Megan Larsen](http://meganllarsen.wordpress.com).

[Ariane Peralta](https://www.peraltalab.com/).

[Dr. Jay Lennon](http://www.indiana.edu/~microbes/people.php): Principle Investigator, Professor, Department of Biology, Indiana University, Bloomington. Head of the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).

