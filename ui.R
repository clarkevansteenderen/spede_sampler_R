
mypackages <- c("ape", "shiny", "shinyhelper", "magrittr", "shinyFiles", "shinythemes", "splits")
checkpkg <- mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(checkpkg)) install.packages(checkpkg, dependencies = TRUE)

library(shiny)
library(ape)
library(splits)
library(shinyhelper)
library(magrittr) # allows you to use %>%
library(shinythemes)
library(shinyFiles)


ui <- fluidPage(theme = shinytheme("cyborg"),
                
                titlePanel(strong(h1("SPEDE-SAMPLER-GMYC"))),
                titlePanel(h3("Run GMYC analyses on resampled ML trees", )),
                #img(src='logo.png',  height = '150px', width = '200px'),
                titlePanel(h4("Created by Clarke van Steenderen")),
                br(),
                
                tabsetPanel(tabPanel("Home",
                    br(), br(),
                    #sidebarLayout(
                        sidebarPanel(strong("Click to view the help file:", style = "color:lightgreen")
                            %>% helper(type = "markdown", content = "spede_sampler_help", colour = "lightgreen", icon = "question-circle"),
                            br(),
                            shinyDirButton('directory', 'Folder select', 'Please select a folder', style="color: black; background-color: white; border-color: black"),
                            br(), br(),
                            textOutput('folder_path'),
                            #verbatimTextOutput('wd'),
                            br(), br(), br(),
                            radioButtons("data_type", label="Maximum Likelihood program used to produce tree files:", choices = c("FastTree", "RAxML")),
                            actionButton("run_gmyc", label = strong("Run"), style="color: black; background-color: white; border-color: black", icon("check")),
                            br(),br(),
                            htmlOutput("complete"), tags$head(tags$style("#complete{color: lightgreen; font-size: 25px;}")),
                            br(),
                            htmlOutput("warn_files"), tags$head(tags$style("#warn_files{color: red; font-size: 20px;}")),
                            br(),br(),
                            downloadButton("download_clust_ent_data", label = "Download Data", style="color: black; background-color: white; border-color: black")
                                                           ),
                    
                            mainPanel() ),
         
                            tabPanel("Plot Results", 
                                     br(), br(),
                                     actionButton("plot_clusts", label = "Plot", style="color: black; background-color: white; border-color: black"), 
                                     br(), br(),
                                     downloadButton("download_clust_plot", label = "Download Plot", style="color: black; background-color: white; border-color: black"),
                                     plotOutput("clust_ent_plot", height = "650px")
                                     )
                            
                )
)

#)