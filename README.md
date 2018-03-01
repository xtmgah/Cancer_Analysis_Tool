# Cancer Analysis Tool
M. Bayati, H.R. Rabiee, et al., H. Alinejad-Rokny, “CANCERSIGN: a user-friendly and robust tool for identification and classification of mutational signatures and patterns in cancer genomes”, preparing for submission.

Table of contents
=================
  * [Prerequisites](#prerequisites)
    * [Platform](#platform)
    * [R](#r)
  * [Usage](#usage)
    * [General overview](#general-overview)
    * [Preprocessing tab](#preprocessing-tab)
    * [3-mer signatures tab](#3-mer-signatures-tab)
    * [5-mer signatures tab](#5-mer-signatures-tab)
    * [Clustering tab](#clustering-tab)
    * [Simulation tab](#simulation-tab)
  	* [Terminate the tool](#terminate-the-tool)


Prerequisites
=============

Platform
-------------------
This tool is currently applicable in **Linux** and **MacOSX**.

R
-
You need to have **R** installed on your machine. In addition, the following **R packages** are required to be installed in R default library path:
* BSgenome.Hsapiens.UCSC.hg19
* doParallel
* ggplot2
* ggfortify
* plotly
* grid
* gridExtra
* reshape2
* shiny
* shinyjs
* rhandsontable

Usage
=====

General overview
----------------
The overall procedure of using this tool is as follows: 
1. Download the whole project folder (**_Cancer_Analysis_Tool_** folder) and open it.
2. Put the mutational catalog file (your data file) in the **_data_** folder.
3. Double click on the **_Run_** file and use the tool!

See the screenshot below:

![1](https://user-images.githubusercontent.com/16561858/30983872-8e4c0496-a498-11e7-959c-af6baab07e7e.png) 

_The first 3 folders in the project folder (see the above image) are interactively used between the user and the software. In this manual, we refer to them  as "**data** folder", "**output** folder" and "**results** folder" respectively._
 
Once the user-interface is opened, you will see the **Preprocessing tab** where you can specify the input data file and get the tool to perform preprocessing necessary for analyzing the data. After this step, you will be shown the main tabs of the tool: **3-mer signatures tab**, **5-mer signatures tab**, **Clustering tab** and **Simulation tab**. We will describe all of the mentioned tabs in detail in the following sections.

Preprocessing tab
-----------------
This is the first thing you will see when the user-interface is opened. It is assumed that you have put your input file in the **_data_** folder and now you need to introduce it to the tool by writing its name (with extension) in the field at the top of the page. Then you should specify the input file format.

The file format can be **ICGC**, **TCGA** or **Custom**. Make sure that your input file is well formatted in each case:
* For **ICGC** files, make sure that the table has the following columns:
	* _icgc_sample_id_
    * _chromosome_
    * _chromosome_start_
    * _chromosome_end_
    * _reference_genome_allele_
    * _mutated_to_allele_
 
* For **TCGA** files, make sure that the table has the following columns:
	* _?_
    * _?_
    * _?_
    * _?_
    * _?_
    * _?_

* For **Custom** files, make sure that the table has the following columns **_with exact names_**:
	* _sample_id_   
    * _chromosome_  
    * _position_    
    * _reference_   
    * _mutated_to_ 
    
After specifying the input file, you can click on the green button to start preprocessing. When this step is performed, all previous results stored in the **_result_** folder and **_output_** folder will be erased. Note that for a specific data file, the preprocessing is only needed to be performed for one time and its results will be stored in **_result_** folder. You will normally want to run the tool for a specific input data for many other times to perform different analyses, so you wish to skip this step in each subsequent run. This can be done by clicking on the orange button which causes the tool to bypass the preprocessing step, and gives you the options for preventing the tool from erasing the contents of the **_output_** folder. See the screenshot below:

![2 with arrows](https://user-images.githubusercontent.com/16561858/30988310-ed844d16-a4a6-11e7-9cbe-aafbb1dddb09.png)


3-mer signatures tab
-----------------
This tab is dedicated to extracting 3-mer mutational signatures from mutational catalogue of input data. First, you are asked 

Do not close the terminal window.
To restrart the tool, just double click on the run file again.


*
*
*
*
*
*
This manual will be completed soon...
=====================================
