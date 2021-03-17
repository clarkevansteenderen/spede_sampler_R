<img src="https://github.com/CJMvS/spede-sampler/blob/main/www/spede_sampler_R_logo.png" height="200" align="right">

# SPEDE-SAMPLER (GMYC)

*Created by:*

*Clarke van Steenderen*

*Department of Zoology and Entomology*

[*The Centre for Biological Control*](https://www.ru.ac.za/centreforbiologicalcontrol/)

*Rhodes University, Grahamstown, Eastern Cape, South Africa*

*2021*

*e-mail:* vsteenderen@gmail.com

---

This R Shiny App is the final part of the SPEDE-SAMPLER program, and allows the user to run GMYC analyes on multiple phylogenies that have been created by randomly resampling the sequences in an aligned Fasta file.

To run this app through R, type the following into the console:

`install.packages("shiny") # install the shiny package` 

`library(shiny) # load up the shiny library` 

`shiny::runGitHub("spede-sampler", "CJMvS", ref="main") # run the app`

# **USER GUIDE**

**OVERVIEW**

This R Shiny App is the final step of the analysis pipeline following from the SPEDE-SAMPLER Python program (see the diagram below). The application requires the user to input a folder directory containing all the tree files created by FastTree or RAxML. If the user wishes to run the analysis on only one tree, this tree file needs to be saved into a folder first, which can then be selected.
Using the "ape" package, each tree is opened and converted to become fully dichotomous (**multi2di()** function) and ultrametric (**chronos()** function).
The GMYC species delimitation algorithm is then run on each tree using the R "splits" package. The number of clusters and entities for each tree is recorded in a dataframe, and can be downloaded and/or plotted in the application under the "Plot Results" tab.
Each GMYC clustering tree can be viewed and downloaded.
If the user has predefined grouping data for their samples, this can be uploaded as an Excel csv file. These predefined groups are then compared to the groups estimated by the GMYC analysis, and a percentage match is calculated.
 
---

<img src="spede_sampler_overview.png" alt="drawing" width="850"/>

**USAGE**

Select a folder, or manually input the file path containing the files created by either FastTree or RAxML. 
Select the approapriate radio button to indicate which ML program was used to create your tree files.
The file path will display on the screen as confirmation of your choice.

To upload predefined grouping information, browse for the relevant .csv file. 

The csv file needs one column for sample names, and another for their corresponding predefined groups. For example:


| sample_id | group |
|-----------|-------|
| MN1234    | sp1   |
| MN1235    | sp1   |
| MN1236    | sp3   |

<br />
Where the **group** column could be, for example, morphospecies.
Select which column is the grouping, and which is the sample name column from the dropdown menus.

Click the "Run" button to start the GMYC analysis. 

Once the analysis is complete, the "View Data" tab allows for the viewing and downloading of data tables for the estimated number of clusters and entities. When downloading, add the .csv extension manually to the file name.

The "Plot Results" tab allows for the viewing, editing (line and point colours), and downloading of plots for:

- The number of clusters vs the number of entities as a scatter plot
- A boxplot for the numbers of clusters and entities
- The number of clusters and entities estimated by each ML iteration file

When downloading these plots, add the .svg extension to the file name.

The "Plot trees" tab allows for the plotting of any individual GMYC tree with support values for detected species groups from those originally uploaded. 

The "Percentage matches" tab allows the user to select any input tree file, and view 1) the GMYC results table with the appended predefined groups (as uploaded by the user), the percentage matches, percentage of single-sample GMYC species, and the oversplitting ratio. 


