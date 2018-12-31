


# CANCERSIGN manual
M. Bayati, H.R. Rabiee, et al., and H. Alinejad-Rokny, “**CANCERSIGN: a user-friendly and robust tool for identification and classification of mutational signatures and patterns in cancer genomes**”, preparing for submission.

Table of contents
=================
  * [Prerequisites](#prerequisites)
    * [Operating System](#operating-system)
    * [R](#r)
  * [Usage](#usage)
    * [Overview](#overview)
    * [Preprocessing tab](#preprocessing-tab)
    * [Terminate or Restart](#terminate-or-restart)


Prerequisites
=============

Operating System
-------------------
The recommended operating systems for using this tool are **Linux** and **MacOSX**.

R
-
You need to have **R** installed on your machine. In addition, the following **R packages** are required to be installed in R default library path:
* BSgenome.Hsapiens.UCSC.hg19
* doParallel
* data.table
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

Overview
----------------
The overall procedure of using this tool is as follows: 
1. Download the whole project folder (**_CANCERSIGN_** folder) and open it.
2. Put the mutational catalog file (your data file) in the **_data_** folder.
3. Double click on the **_Run_** file.

![1](https://user-images.githubusercontent.com/16561858/30983872-8e4c0496-a498-11e7-959c-af6baab07e7e.png) 

_The "**data** folder", "**output** folder" and "**result** folder" are interactively manipulated by the user and the tool._

Once the user-interface is opened, you will see the **Preprocessing tab** where you can specify the input data file and get the tool to perform preprocessing necessary for analyzing the data. After this step, you will be shown the main tabs of the tool: **3-mer signatures tab**, **5-mer signatures tab** and **Clustering tab**. We will describe them in more detail below.

Preprocessing tab
-----------------
This is the first thing you will see when the user-interface is opened. It is assumed that you have put your input file in the **_data_** folder and now you need to introduce it to the tool by writing its name (with extension) in the field at the top of the page. 

The valid input file for CANCERSIGN must contain the following fields:

1. **sample_id**
2. **chromosome**
3. **position**
4. **reference**
5. **mutated_to**
    
After specifying the input data file, you can click on the green button to start preprocessing. When this step is performed, the results are stored in the **_result_** folder and its old files are deleted. The contents of the **_output_** folder are also cleared. However, the user may want to repeat running the tool to perform various analyses for the same data. In such cases, the preprocessing results can remain in the **_result_** folder and the user can bypass the preprocessing step by clicking on the orange button. The user is then asked whether he wants to preserve the contents of the **_output_** folder or let them be erased for subsequent runs.

![2 with arrows](https://user-images.githubusercontent.com/16561858/30988310-ed844d16-a4a6-11e7-9cbe-aafbb1dddb09.png)


Terminate or Restart
-----------------
To terminate the tool, close the corresponding terminal window.  
To restart the tool, just double click on the **_Run_** file again.
<!--stackedit_data:
eyJoaXN0b3J5IjpbOTk1MzU2Mzg4LDY3OTYzNDk4Nl19
-->
