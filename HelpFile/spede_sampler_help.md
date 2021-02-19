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
Using the "ape" package, each tree is opened and converted to a fully dichotomous (multi2di() function) and ultrametric (chronopl() function) phylogeny.
A GMYC species delimitation algorithm is then run on each tree using the R "splits" package. The number of clusters and entities for each tree is recorded in a dataframe, and can be downloaded and/or plotted in the application under the "Plot Results" tab.

---

**USAGE**

Please select a folder containing the files created by either FastTree or RAxML, and select the approapriate radio button to indicate which program was used.
The file path will display on the screen as confirmation of your choice.

Click the "Run" button to start the GMYC analysis. Once it has completed, you can download the number of clusters and entities for each phylogeny as a .csv file (please add the .csv extension after the name before saving) and/or plot the results in the application.
This plot can be downloaded as a .svg file (again, add the .svg after the name before saving). 

---

