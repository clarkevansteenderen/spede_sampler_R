
mypackages <- c("ape", "shiny", "shinyhelper", "magrittr", "shinyFiles", "shinythemes", "shinyalert", "splits", "phytools", "reshape2", "devtools", "ggplot2")
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
library(ggplot2)

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
                            br(),
                            checkboxInput("set_seed", label = "Set a seed?", value = FALSE),
                            
                            checkboxInput("group_info", label="Check this box if you wish to upload predefined grouping information. If yes, upload a .csv file, and select the columns containing your grouping information and sample names from the dropdown menu.", value = FALSE),
                            br(),
                            fileInput("predefined_groups", label="Upload a .csv file containing predefined groups for your samples:", accept = ".csv"),
                            selectInput("col.group", "Select Group Column:", choices=NULL),
                            selectInput("sample_names", "Select Sample Name Column:", choices=NULL),
                            br(),
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
                                     
                                     ),
                    
                          tabPanel(strong("Plot Trees", style = "color:darkblue"),
                                    br(), br(),
                                    selectInput("select_tree", label = "Select a tree to plot", choices = NULL),
                                    sliderInput("tip_label_size", label = "Tip label size: ", value = 0.8, min = 0.1, max = 2),
                                    sliderInput("support_value_size", label = "Support value size: ", value = 1, min = 0.1, max = 5),
                                    sliderInput("line_width", label = "Branch line width: ", value = 1, min = 0.1, max =5),
                                    selectInput("support_value_col", label = "Support value colour: ", choices = c("grey", "lightblue", "salmon", "lightgreen", "lightyellow"), selected = "lightgreen"),
                                    selectInput("support_value_frame", label = "Support value frame:", choices = c("none", "circle", "rect"), selected = "circle"),
                                    selectInput("branch_col", label = "Branch colour: ", choices = c("black", "blue", "lightblue", "red", "green", "orange"), selected = "blue"),
                                    br(),
                                    actionButton("plot_gmyc_tree", label = strong("Plot GMYC tree result"), style="color: black; background-color: lightpink; border-color: black", icon("drafting-compass")),
                                    br(), br(),
                                    plotOutput("gmyc_tree", height = "600px")
                          ),
                      
                          tabPanel(strong("Percentage Matches", style = "color:darkblue"),
                                    br(), br(),
                                    actionButton("view_match_data", label = strong("View Matches"), style="color: black; background-color: lightgreen; border-color: green", icon("edit")),
                                    actionButton("view_summary_match_data", label = strong("View Matches Summary"), style="color: black; background-color: lightblue; border-color: darkblue", icon("edit")),
                                    br(), br(),
                                    tableOutput("matches")
                          
                        )
                            
                )
)

#)