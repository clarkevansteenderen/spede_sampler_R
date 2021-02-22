## **SPEDE-SAMPLER: GMYC ANALYSIS** 

### **USER GUIDE**

---

*Created by:*

*Clarke van Steenderen*

*Department of Zoology and Entomology*

[*The Centre for Biological Control*] (https://www.ru.ac.za/centreforbiologicalcontrol/)

*Rhodes University, Grahamstown, Eastern Cape, South Africa*

*2021*

*e-mail:* vsteenderen@gmail.com

---

**OVERVIEW**

This R Shiny App is the final step of the analysis pipeline following from the SPEDE-SAMPLER Python program. The application requires the user to input a folder directory containing all the tree files created by FastTree or RAxML.
Using the "ape" package, each tree is opened and converted to a fully dichotomous (multi2di() function) and ultrametric (chronos() function) phylogeny.
A GMYC species delimitation algorithm is then run on each tree using the R "splits" package. The number of clusters and entities for each tree is recorded in a dataframe, and can be downloaded and/or plotted in the application under the "Plot Results" tab.

---

**USAGE**

Please either select a folder, or manually input the file path containing the files created by either FastTree or RAxML. 
Select the approapriate radio button to indicate which ML program was used to create your tree files.
The file path will display on the screen as confirmation of your choice.

Click the "Run" button to start the GMYC analysis. 

Once the analysis is complete, the "View Data" tab allows you to view and download data tables for the estimated number of clusters and entities. When downloading, add the .csv extension manually to the file name.


The "Plot Results" tab allows you to view, edit (line and point colours), and download plots for:

- The number of clusters vs the number of entities as a scatter plot
- A boxplot for the numbers of clusters and entities
- The number of clusters and entities estimated by each ML iteration file

When downloading these plots, add the .svg extension to the file name.



