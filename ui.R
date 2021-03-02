
mypackages <- c("ape", "shiny", "shinyhelper", "magrittr", "shinyFiles", "shinythemes", "shinyalert", "splits", "phytools", "reshape2", "devtools")
checkpkg <- mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(checkpkg)) install.packages(checkpkg, dependencies = TRUE)

library(shiny)
library(ape)
library(splits)
library(shinyhelper)
library(magrittr) # allows you to use %>%
library(shinythemes)
library(shinyFiles)
library(reshape2)
library(shinyalert)
library(phytools)

ui <- fluidPage(theme = shinytheme("flatly"),
                
                useShinyalert(), # set up the shinyalert package
                
                titlePanel(strong(h1("SPEDE-SAMPLER-GMYC"))),
                titlePanel(h3("Run GMYC analyses on resampled ML trees", )),
                #img(src='logo.png',  height = '150px', width = '200px'),
                titlePanel(h4("Created by Clarke van Steenderen")),
                br(),
                
                tabsetPanel(tabPanel(strong("Home: multiple ML trees", style = "color:darkblue"),
                    br(), br(),
                    #sidebarLayout(
                        sidebarPanel(strong("Click to view the help file:", style = "color:green; font-size: 18px")
                            %>% helper(type = "markdown", content = "spede_sampler_help", colour = "green", icon = "question-circle"),
                            br(),
                            h5(strong('Select the folder containing your tree files:')),
                            shinyDirButton('directory', 'Folder select', 'Please select a folder containing your tree files', style="color: black; background-color: white; border-color: black"),
                            br(), br(),
                            h5("You have selected the folder path: "),
                            br(),
                            textOutput('folder_path'),
                            br(), br(),
                            h3('OR'),
                            br(),
                            textInput(inputId = 'raw_file_path', label = 'Manually insert a file path: '),
                            br(),
                            radioButtons("data_type", label="Maximum Likelihood program used to produce tree files:", choices = c("FastTree", "RAxML")),
                            actionButton("run_gmyc", label = strong("RUN"), style="color: black; background-color: lightgreen; border-color: green", icon("check-square")),
                            br(),br()
                                                           ),
                    
                            mainPanel() ),
                    
                            tabPanel(strong("View Data", style = "color:darkblue"),
                                     br(), br(),
                                     actionButton('all_data', label=strong('Show all data'), style="color: black; background-color: lightgreen; border-color: green", icon("edit")),
                                     downloadButton("download_clust_ent_data", label = strong("Download all data"), style="color: black; background-color: lightgreen; border-color: green"),
                                     br(), br(),
                                     actionButton('summary_data', label = strong('Show summary table'), style="color: black; background-color: lightblue; border-color: darkblue", icon("edit")),
                                     downloadButton("download_stat_summary", label = strong("Download summary table"), style="color: black; background-color: lightblue; border-color: darkblue"),
                                     br(), br(),
                                     tableOutput("data_table"),
                                     br()
                               
                            ),
         
                            tabPanel(strong("Plot Results", style = "color:darkblue"), 
                                     br(), br(),
                                     actionButton("plot_clusts", label = strong("Plot clusters vs entities"), style="color: black; background-color: lightgreen; border-color: green", icon("drafting-compass")), 
                                     downloadButton("download_clust_plot", label = strong("Download"), style="color: black; background-color: lightgreen; border-color: green"),
                                     selectInput("clust_vs_ent_plot_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen")),
                                     br(), 
                                     actionButton("plot_boxplot", label = strong("Plot boxplot"), style="color: black; background-color: lightblue; border-color: darkblue", icon("drafting-compass")), 
                                     downloadButton("download_boxplot", label = strong("Download"), style="color: black; background-color: lightblue; border-color: darkblue"),
                                     br(), br(), br(),
                                     actionButton("plot_clusts_vs_iterations", label = strong("Plot clusters vs iteration file"), style="color: black; background-color: lightyellow; border-color: black", icon("drafting-compass")), 
                                     downloadButton("download_clusts_vs_iterations", label = strong("Download"), style="color: black; background-color: lightyellow; border-color: black"),
                                     selectInput("plot_clusts_vs_iterations_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen")),
                                     selectInput("plot_clusts_vs_iterations_line_colour", "Line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen")),
                                     br(), br(),
                                     actionButton("plot_ents_vs_iterations", label = strong("Plot entities vs iteration file"), style="color: black; background-color: lightpink; border-color: black", icon("drafting-compass")), 
                                     downloadButton("download_ents_vs_iterations", label = strong("Download"), style="color: black; background-color: lightpink; border-color: black"),
                                     selectInput("plot_ents_vs_iterations_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen")),
                                     selectInput("plot_ents_vs_iterations_line_colour", "Line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen")),
                                     br(), br(),
                                     plotOutput("clust_ent_plot", height = "600px"),
                                     br(), br(),
                                     
                                     )
                            
                )
)

#)