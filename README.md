University Lake Reservoir Gradient Project
==========

This repository contains open-source code, data, & text files for the Reservoir Gradient Project.
This project is the product of the 2013 SSRP work done with Nelson.

## Project Goals

* **Aim 1)** Determine if there is an underlying stream to lake gradient in small reservoirs.

* **Aim 2)** Determine if microbial communities change along a stream to lake gradient.

* **Aim 3)** Determine if microbial functions change along a stream to lake gradient.

### Repo Contents

* **analyses:**
	* *ULReservoirGradient_MothurSummary.R*: An R script to test mothur output and return basic information including sequence coverage and species richness.
	Used to determine if mothur pipeline returned logical output.


* **data:**
	* *res.grad.metab.txt*:
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
	* *UL.bac.final.0.03.tax.summary*
	* *UL.bac.final.0.03.taxonomy*:
	* *UL.bac.final.0.05.tax.summarY*:
	* *UL.bac.final.0.05.taxonomy*:
	* *UL.bac.final.shared*: Mothur output shared file.
	Used to generate the site-by-species matrix.
	* *UL_euks.final.shared*: Mothur output shared file for eukaryotes.


* **bin:**
	* *MothurTools.R*: An R script written by Mario Muscarella (Indiana University) containing functions used in the analysis of mothur raw files.


* **figures:**
  * *Figure1.png*: Microbial metabolism figure.


* **maps:**


* **mothur:**
  * *UL.raw.seqs.sh*: Shell script to download raw illumina data from ()  for Reservoir Gradient Project (work in progress).  
  * *UL.files*: Mothur input file that organizes raw sequence files into read 1 and read 2 by sample design.  
  * *UL.Bacteria_A.Batch*: Mothur batch script for the initial sequence processing.
	This script goes through chimera checking and the output is used by all subsequent batch scripts.
  * *UL.Bacteria_B.Batch*: Mothur batch script for analysis of bacterial sequences.
  * *UL.Cyano.Batch*: Mothur batch script for analysis of cyanobacteria only
  * *UL.Euks.Batch*: Mothur batch script for anlysis of eukaryotes based on plastid sequences
  * *UL.Bacteria_A.sh*: Torque script to execute UL.Bacteria_A.Batch on Mason (IU).   
	Ends up submitting all subsequent torque scripts.
  * *UL.Bacteria_B.sh*: Torque script to execute UL.Bacteria_B.sh on Mason (IU).
  * *UL.Cyano.sh*: Torque script to execute UL.Cyano.Batch on Mason (IU).
  * *UL.Euks.sh*: Torque script to execute UL.Euks.Batch on Mason (IU).


## Contributors

[Nathan Wisnoski](https://nwisnoski.github.io/).

[Mario Muscarella](http://mmuscarella.github.io/).

[Megan Larsen](http://meganllarsen.wordpress.com).

[Ariane Peralta](https://www.peraltalab.com/).

[Dr. Jay Lennon](http://www.indiana.edu/~microbes/people.php): Principle Investigator, Professor, Department of Biology, Indiana University, Bloomington. Head of the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).

