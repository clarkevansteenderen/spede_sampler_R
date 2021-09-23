
mypackages <- c("ape", "data.table", "dplyr", "shiny", "shinyWidgets", "ips", "Rmisc", "shinyhelper", "gtools", "magrittr", "shinyFiles", "shinythemes", "shinyalert", "splits", "phytools", "reshape2", "devtools", "ggplot2", "tools")
checkpkg <- mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(checkpkg)) install.packages(checkpkg, dependencies = TRUE)

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("Biostring")

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
library(gtools)
library(Rmisc)
library(dplyr)
library(data.table)
library(shinyWidgets)
library(ips)
library(tools)
library(Biostrings)

#devtools::install_github("GuillemSalazar/FastaUtils")
#library(FastaUtils)

ggthemes = list("Classic" = theme_classic(),
                "Dark" = theme_dark(),
                "Minimal" = theme_minimal(),
                "Grey" = theme_grey(),
                "Light" = theme_light(),
                "Black/White" = theme_bw(),
                "Void" = theme_void())

ui <- fluidPage(
  
  # this code suppresses error warnings on the shiny console (e.g. "Warning: Error in [<-.data.frame: new columns would leave holes after existing columns")
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  
  
  tags$head(tags$style(".shiny-notification {font-size: 25px; background-color: #ffffff;")),
  bootstrapPage('',
                
                tags$style(type = 'text/css',
                           HTML('.navbar {background-color: #79aaaa; font-size: 18px;}
                           .navbar-default .navbar-brand {color: #ffffff; font-size: 25px;}
                           .navbar-default .navbar-nav > .active > a, 
                           .navbar-default .navbar-nav > li > a {color:black;}
                           .navbar-default .navbar-nav > .active > a:focus, 
                           .navbar-default .navbar-nav > .active > a:hover {color: grey; background-color: #ffffff;}
                           .nav-tabs>li>a {color: black; background-color:#79aaaa; border-color: white}
                           .well {background-color: #E7F0F1;}'
                           )
                )),
  
                navbarPage(title = "SPEDE-SAMPLER-GMYC", id = "tabset",
                          tabPanel("Home",
                          wellPanel(align = "justify", style = "background: white",
                          br(), br(),                
                          div(img(src='spede_sampler_gmyc_logo.png',  height = 'auto', width = '70%'), style="display: block; margin-left: auto; margin-right: auto; text-align: center;"),
                          br(), br(),
                          fluidRow(actionButton("app", strong("BEGIN"), style = 'font-size:150%; color: black; background-color: #A6C3C6; border-color: black'), align = "center"),
                          br(), br(),
                          div(img(src='cover_insects.png',  height = 'auto', width = '80%'), style="display: block; margin-left: auto; margin-right: auto; text-align: center;"),
                          
                          # br(), br(),
                          # HTML("SPEDE-SAMPLER offers the following functionalities: 
                          #                                         (1) Runs GMYC analyses on multiple Maximum Likelihood tree files and records the number of clusters and entities, 
                          #                                           (2) Calculates the match percentage between user-defined grouping information and species estimated by the GMYC method,
                          #                                           (3) Calculates to what extent the GMYC method oversplits species relative to predefined groupings,
                          #                                           (4) Extracts desired output data columns from separate runs to amalgamate and plot on one graph,
                          #                                           and (5) downloads figures and Excel data produced by the application.
                          #                                           View plots and/or data at the bottom of each tab window.
                          #                                       "),      
                         
                                    
                                    )
                                    ),
                          
                           tabPanel("App", value = "app",
                                    
                                    tabsetPanel(
                                      
                                      tabPanel(strong("Random resampling"),
                                               br(),
                                               img(src="cornops.png", align = "left", height="20%", width="20%"),
                                               br(), br(), br(), br(),
                                               h3(strong("Upload a Fasta file to randomly sample")),
                                               strong("Randomly resample a desired percentage of sequences in an uploaded multiple sequence alignment file, without replacement"),
                                               br(), br(),
                                        wellPanel(
                                          fileInput("fasta_file", label="Upload a .fas file:", accept = c(".fas", ".fasta")),
                                          br(),
                                          textOutput("num_seqs"),
                                          br(),
                                          textInput("resampled_fasta_folder_name", label = 'Output fasta folder name: ', value = "RESAMPLED_FASTA_FILES"),
                                          br(),
                                          textInput("resampled_fasta_file_name", label = 'Output fasta file name: ', value = "resampled"),
                                          br(),
                                          numericInput("fasta_subsample_percent", "Percentage of sequences to resample:", value = 10,  min = 1, step = 1, width = "250px"),
                                          br(),
                                          numericInput("fasta_resample_iterations", "Number of iterations: ", value = 2, min = 1, step = 1, width = "250px"),
                                          br(),
                                          checkboxInput("set_seed_resampling", label = strong("Set a seed?"), value = FALSE),
                                         fluidRow( actionButton("resample_fastas", strong("Resample"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black"), align = "center")
                                          
                                        )       
                                        
                                      ),
                                      
                                      tabPanel(strong("BEAST XML Files"),
                                               br(),
                                               img(src="beast.png", align = "left", height="10%", width="10%"),
                                               br(), br(),
                                               h3(strong("Create XML files")),
                                               strong("Use functionality from BEAUti to create an input .xml file for a BEAST analysis"),
                                               br(), br(), br(),
                                               wellPanel(
                                                 textInput(inputId = 'resampled_fasta_file_path', label = 'Manually insert a file path to the folder containing your resampled Fasta files: '),
                                                 br(),
                                                 fluidRow(
                                                   column(width = 2,
                                                 selectInput("beast_site_model", "Site model", choices = c("GTR", "HKY", "JC69", "TN93"), width = "250px"),
                                                   ),
                                                 column(width = 4,
                                                 selectInput("beast_clock_model", "Clock model", choices = c("Strict", "Relaxed lognormal"), width = "250px"),
                                                 ),
                                                 column(width = 2,
                                                        numericInput("beast_clock_rate", "Clock rate", value = 1, min = 0, width = "250px"),
                                                 ),
                                                 ),
                                                 fluidRow(
                                                   
                                                 column(width = 4,
                                                 selectInput("beast_tree_prior", "Tree prior", choices = c("Birth-death", "Coalescent Bayesian skyline", "Coalescent constant-population", "Coalescent exponential-population", "Yule"), selected = "Yule", width = "250px"),
                                                 ),
                                                 
                                            conditionalPanel(
                                              condition = "input.beast_tree_prior == 'Birth-death'",
                                                 
                                                 column(width = 4,
                                                 selectInput("distr_b", "Birth rate distribution", choices = c("Beta", "Exponential", "Gamma", "Inverse gamma", "Laplace", "Log-normal", "Normal", "1/X", "Poisson", "Uniform"), selected = "Uniform", width = "250px"),
                                                 ),
                                                
                                                 column(width = 4,
                                                 selectInput("distr_d", "Death rate distribution", choices = c("Beta", "Exponential", "Gamma", "Inverse gamma", "Laplace", "Log-normal", "Normal", "1/X", "Poisson", "Uniform"), selected = "Uniform", width = "250px"),
                                                 ),
                                              
                                            ), # end conditionpanel
                                            
                                            conditionalPanel(
                                              condition = "input.beast_tree_prior == 'Coalescent Bayesian skyline'",
                                              
                                              column(width = 3,
                                                     numericInput("cbs_group_sizes_dim", "Group size dimension", value = 5, min = 0),
                                              ),
                                              
                                            ), # end conditionpanel
                                            
                                            conditionalPanel(
                                              condition = "input.beast_tree_prior == 'Coalescent constant-population'",
                                              fluidRow(
                                              column(width = 3,
                                                     selectInput("distr_ccp", "Rate distribution", choices = c("Beta", "Exponential", "Gamma", "Inverse gamma", "Laplace", "Log-normal", "Normal", "1/X", "Poisson", "Uniform"), selected = "1/X", width = "250px"),
                                              ),
                                              column(width = 3,
                                                     numericInput("pop_size_distribution", "Population distribution", value = 0.3, min = 0),
                                              ),
                                              ),
                                              
                                            ), # end conditionpanel
                                            
                                            conditionalPanel(
                                              condition = "input.beast_tree_prior == 'Coalescent exponential-population'",
                                              
                                              column(width = 4,
                                                     selectInput("distr_cep_pop", "Population distribution", choices = c("Beta", "Exponential", "Gamma", "Inverse gamma", "Laplace", "Log-normal", "Normal", "1/X", "Poisson", "Uniform"), selected = "1/X", width = "250px"),
                                              ),
                                              
                                              column(width = 4,
                                                     selectInput("distr_cep_gr", "Growth rate distribution", choices = c("Beta", "Exponential", "Gamma", "Inverse gamma", "Laplace", "Log-normal", "Normal", "1/X", "Poisson", "Uniform"), selected = "Laplace", width = "250px"),
                                              ),
                                              
                                            ), # end conditionpanel
                                            
                                            conditionalPanel(
                                              condition = "input.beast_tree_prior == 'Yule'",
                                              
                                              column(width = 3,
                                                selectInput("distr_yule", "Birth rate distribution", choices = c("Beta", "Exponential", "Gamma", "Inverse gamma", "Laplace", "Log-normal", "Normal", "1/X", "Poisson", "Uniform"), selected = "Uniform", width = "250px"),   
                                              ),
                                              
                                            ), # end conditionpanel
                                            
                                                 column(width = 2,
                                                 numericInput("beast_mcmc", "MCMC:", value = 1000, min = 1, step = 1, width = "250px"),
                                                 ),
                                                 
                                                 column(width = 2,
                                                 numericInput("beast_store_every", "Store every:", value = 1000, min =1000, width = "250px"),
                                                 ),
                                                 
                                                 ), # end of fluid row
                                                 
                                                 textInput("xml_folder_name", label = 'Output XML folder name: ', value = "XML_FILES"),
                                                 fluidRow(actionButton("create_xml_files", strong("Generate"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black"), align = "center")
                                               ),
                                               
                                               ),
                                              
                                      tabPanel(strong("Run BEAST"),
                                               br(),
                                               img(src="beast.png", align = "left", height="10%", width="10%"),
                                               br(), br(), br(),
                                               h3(strong("Run BEAST2")),
                                               strong("Use the generated XML files to run BEAST2, via the R package 'beastier'. Consider using the CIPRES server as an alternative for large sequence alignments."),
                                               br(), br(),
                                               wellPanel(
                                                 textInput(inputId = 'xml_file_path', label = 'Manually insert a file path to the folder containing your XML files: '),
                                                 checkboxInput("run_beagle", strong("Run BEAGLE?"), value = FALSE),
                                                 #textInput("beast_output_folder_name", label = 'BEAST output folder name: ', value = "BEAST_OUTPUT"),
                                                 fluidRow(actionButton("run_beast", strong("Run BEAST2"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black"), align = "center")
                                               ),
                                      ),
                                      
                                      tabPanel(strong("Run LogCombiner"),
                                               br(),
                                               img(src="pillbug.png", align = "left", height="15%", width="15%"),
                                               br(), br(), br(),
                                               h3(strong("Run LogCombiner")),
                                               strong("Use LogCombiner to optionally reduce the size of the .trees file generated by BEAST by resampling states at a lower frequency. If trees were sampled at every 1000 states, these can be thinned to 2000 or 3000, for example. This will ensure that TreeAnnotator does not run out of memory when processing trees."),
                                               br(), br(),
                                               strong("* Make sure that the LogCombiner.exe file is in the same directory as the input files"),
                                               br(), br(),
                                               wellPanel(
                                                 textInput(inputId = 'logcombiner_file_path', label = 'Manually insert a file path to the folder containing your .trees files: '),
                                                 numericInput("resample_freq", "Resampling frequency:", value = 2000, width = "250px"),
                                                 textInput("logcombiner_folder_results", "Name of the output folder:", value = "LogCombiner_output"),
                                                 fluidRow(actionButton("run_logcombiner", strong("Run"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black"), align = "center")
                                               )
                                      ),
                                      
                                      tabPanel(strong("Run TreeAnnotator"),
                                               br(),
                                               img(src="beast.png", align = "left", height="10%", width="10%"),
                                               br(), br(),
                                               h3(strong("Run TreeAnnotator")),
                                               strong("* Ensure that a copy of the TreeAnnotator.exe application is in the same directory as the data files"),
                                               br(), br(),
                                               wellPanel(
                                                 textInput(inputId = 'beast_trees_file_path', label = 'Manually insert a file path to the folder containing your BEAST .trees files: '),
                                                 fluidRow(
                                                   column(width = 3,
                                                 numericInput("treeannotator_burnin", "Burnin (% trees):", value = 0, min = 0, width = "250px"),
                                                 checkboxInput("treeannotator_low_mem", "Low Memory?", value = TRUE),
                                                   ),
                                                 column(width = 2,
                                                 selectInput("treeannotator_heights", "Heights:", choices = c("median", "keep", "mean", "ca")),
                                                 ),  
                                                 ), # end of fluidrow
                                                 fluidRow(actionButton("run_treeannotator", strong("Run"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black"), align = "center")
                                               )
                                               
                                      ),
                                      
                                      tabPanel(strong("Tracer"),
                                        br(),
                                        img(src="beast.png", align = "left", height="10%", width="10%"),
                                        br(), br(),
                                        h3(strong("Tracer")),
                                        strong("Check that ESS scores > 200, and for MCMC convergence"),
                                        br(), br(),
                                        wellPanel(
                                          textInput("beast_log_files_path", "Manually insert a file path to the folder containing your BEAST log files: "),
                                          fluidRow(actionButton("load_log_files", strong("Load"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black"), align = "center"),
                                          numericInput("tracer_burnin", "Burnin:", value = 0.1, min = 0, width = "250px"),
                                          numericInput("tracer_sample_interval", "Sample interval:", value = 1000, min = 0, width = "250px"),
                                          selectInput("select_log_file", label = "Select a log file:", choices = NULL),
                                          fluidRow(actionButton("table_log_files", strong("Results"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black"), align = "center")
                                        ),
                                        column(4,
                                        tableOutput("tracer_ess"),
                                        ),
                                        column(8,
                                        plotOutput("tracer_plot"),
                                        )
                                        
                                      ),
                                      
                                      tabPanel(strong("GMYC Analysis"),
                                               br(),
                                               img(src="beast_octopus.png", align = "center", height="15%", width="15%"),
                                               br(), 
                                               h3(strong("Upload multiple tree files")),
                                               strong("Run GMYC analyses on all the tree files in a designated folder"),
                                               br(), br(),
                                               wellPanel(
                                               #              h5(strong('Select the folder containing your tree files:')),
                                               #              shinyDirButton('directory', 'Folder select', 'Please select a folder containing your tree files', style="color: black; background-color: white; border-color: black"),
                                               #              br(), br(),
                                               #              
                                               #              conditionalPanel(
                                               #                condition = "input.directory",
                                               #                h5(strong("You have selected the folder path: ")),
                                               #                br(),
                                               #                textOutput('folder_path'),
                                               # ),
                                               #              
                                               #              br(), br(),
                                               #              h3('OR'),
                                               #              br(),
                                                            textInput(inputId = 'raw_file_path', label = 'Manually insert a file path to the folder containing your tree files: '),
                                                            br(),
                                                            selectInput("gmyc_threshold", "GMYC threshold:", choices = c("single", "multiple"), width = "250px"),
                                                            br(),
                                                            selectInput("tree_tool", "Select phylogenetic method:", choices = c("BEAST", "Non-ultrametric: PATHD8"), width = "250px"),
                                                  
                                               conditionalPanel(
                                                 condition = "input.tree_tool == 'Non-ultrametric: PATHD8'",
                                                 numericInput("seqlength", "Sequence length of alignment (bp):", value = 1,  min = 0, step = 1, width = "250px"),
                                                 radioButtons("data_type", label="Maximum Likelihood program used to produce tree files:", choices = c("FastTree", "RAxML"))
                                               ),
                                               
                                               
                                                            br(),
                                               ),
                                               
                                               h3(strong("[Optional] Upload a csv file with prior grouping information")),
                                               wellPanel(
                                                            fileInput("predefined_groups", label="Upload a .csv file:", accept = ".csv"),
                                                            actionButton("group_data_uploaded", strong("Confirm file"), style = 'font-size:120%; color: black; background-color: #A6C3C6; border-color: black', icon("thumbs-o-up")),
                                                            br(), br(),
                                                            selectInput("col.group", "Select Group Column:", choices="", width = "350px"),
                                                            selectInput("sample_names", "Select Sample Name Column:", choices="", width = "350px"),
                                               ),
                                               wellPanel(
                                                 tags$style("input[type=checkbox] {transform: scale(2);}"),
                                                 
                                                            tags$div(title="Check if you want the results to be reproducible if you run this again with the same data",
                                                                     checkboxInput("set_seed", label = strong("Set a seed?"), value = FALSE), style = "font-size: 16px",
                                                            ),
                                               ),  
                                                            br(),
                                                            fluidRow(actionButton("run_gmyc", strong("RUN GMYC"), style = 'font-size:150%; color: black; background-color: #A6C3C6; border-color: black'), align = "center"),
                                                            #actionButton("run_gmyc", label = strong("RUN"), style="color: black; background-color: lightgreen; border-color: black", icon("check-square")),
                                                            br(),br()
                                               ),
                                               
                                               #mainPanel(), 
                                      
                                      tabPanel(strong("View Data"),
                                               br(),
                                               img(src="iris_flea_beetle.png", align = "left", height="15%", width="15%"),
                                               br(), br(), br(),
                                               textOutput("sp.numbers"),
                                               br(),
                                               h3(strong("View cluster and entity data")),
                                               strong("Assess the number of clusters and entities estimated by the GMYC analyses"),
                                               br(), br(),
                                               wellPanel(align = "justify",
                                               h4(strong("Print all data to the screen")),
                                               br(),
                                               actionButton('all_data', label=strong('Print'), style="color: black; background-color: #A6C3C6; border-color: black", icon("edit")),
                                               downloadButton("download_clust_ent_data", label = strong("Download all data"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ),
                                               br(),
                                               wellPanel(align = "justify",
                                               h4(strong("Print summary statistics to the screen")),
                                               br(),
                                               actionButton('summary_data', label = strong('Print'), style="color: black; background-color: #A6C3C6; border-color: black", icon("edit")),
                                               downloadButton("download_stat_summary", label = strong("Download summary table"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ),
                                               br(),
                                               tableOutput("data_table"),
                                               br()
                                               
                                      ),
                                      
                                      tabPanel(strong("Plot Results"), 
                                               br(), br(),
                                               img(src="parasitoid.png", align = "left", height="20%", width="20%"),
                                               br(), br(), br(),
                                               h3(strong("Plot cluster-entity results")),
                                               br(),
                                               wellPanel(align = "justify",
                                               selectInput("ggtheme_plots", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"], width = "150px"),
                                               ),
                                               h3(strong("Plots")),
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("CLUSTERS VS ENTITIES")),
                                                         br(),
                                               fluidRow(
                                                 column(width = 3,
                                               selectInput("clust_vs_ent_plot_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen"), width = "150px"),
                                                 ),
                                               column(width = 3,
                                                sliderInput("clust_vs_ent_plot_point_size", label = "Point size: ", value = 1, min = 1, max = 10, step = 0.5),
                                               ),
                                               column(width = 3,
                                                      selectInput("clust_vs_ent_plot_point_shape", label = "Point shape: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18)),     
                                               ),
                                               ),
                                               actionButton("plot_clusts", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")),
                                               hr(),
                                               h4(strong("DOWNLOAD")),
                                               br(), 
                                               textInput("file_name_clustvsent", "File name: ", "clusters_entities"),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                 selectInput("plot_format_clustvsent", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                        ),
                                                 column(width = 3,
                                                        textInput("w_plot_clustvsent", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_clustvsent", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_clustvsent", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_clustvsent == 'png'",
                                                          textInput("res_plot_clustvsent", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_clust_plot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(),
                                               ), br(),
                                               
                                              
                                               wellPanel(align = "justify",
                                                         h4(strong("BOXPLOT")),
                                                         br(),
                                               actionButton("plot_boxplot", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")), 
                                               hr(),
                                               h4(strong("DOWNLOAD")),
                                               br(), 
                                               textInput("file_name_clust_ent_box", "File name: ", "cluster_entity_boxplot"),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_clust_ent_box", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_clust_ent_box", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_clust_ent_box", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_clust_ent_box", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_clust_ent_box == 'png'",
                                                          textInput("res_plot_clust_ent_box", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_boxplot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ), br(),
                                               
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("CLUSTERS VS ITERATION FILE")),
                                                         br(),
                                               fluidRow( 
                                                column(width = 3,
                                               selectInput("plot_clusts_vs_iterations_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen"), width = "150px"),
                                                ),
                                               column(width = 3,
                                                      sliderInput("plot_clusts_vs_iterations_point_size", label = "Point size: ", value = 1, min = 1, max = 10, step = 0.5),
                                               ),
                                               column(width = 3,
                                                      selectInput("plot_clusts_vs_iterations_point_shape", label = "Point shape: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18)),     
                                               ),
                                               column(width = 3,
                                               selectInput("plot_clusts_vs_iterations_line_colour", "Line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen", "white"), width = "150px"),
                                               ),
                                               ),
                                               actionButton("plot_clusts_vs_iterations", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")), 
                                               hr(),
                                               h4(strong("DOWNLOAD")),
                                               br(), 
                                               
                                               textInput("file_name_clustvsiter", "File name: ", "clusters_vs_iterations"),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_clustvsiter", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_clustvsiter", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_clustvsiter", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_clustvsiter", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_clustvsiter == 'png'",
                                                          textInput("res_plot_clustvsiter", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_clusts_vs_iterations", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(),
                                               ), br(),
                                               
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("ENTITIES VS ITERATION FILE")),
                                                         br(),
                                               fluidRow(  
                                                 column(width = 3,
                                               selectInput("plot_ents_vs_iterations_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen"), width = "150px"),
                                                 ),
                                               column(width = 3,
                                                      sliderInput("plot_ents_vs_iterations_point_size", label = "Point size: ", value = 1, min = 1, max = 10, step = 0.5),
                                               ),
                                               column(width = 3,
                                                      selectInput("plot_ents_vs_iterations_point_shape", label = "Point shape: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18)),     
                                               ),
                                               column(width = 3,
                                               selectInput("plot_ents_vs_iterations_line_colour", "Line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen"), width = "150px"),
                                               ),
                                               ),
                                               actionButton("plot_ents_vs_iterations", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")), 
                                               hr(),
                                               h4(strong("DOWNLOAD")),
                                               br(), 
                                               textInput("file_name_entvsiter", "File name: ", "entities_vs_iterations"),
                                              
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_entvsiter", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_entvsiter", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_entvsiter", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_entvsiter", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_entvsiter == 'png'",
                                                          textInput("res_plot_entvsiter", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_ents_vs_iterations", label = strong("Download"), style="color: black; background-color: lightpink; border-color: black"),
                                               br(),
                                               ), 
                                               h3(strong("Plot output:")),
                                               br(), br(),
                                               plotOutput("clust_ent_plot", height = "600px"),
                                               br(), br(),
                                               
                                      ),
                                      
                                      tabPanel(strong("Plot Trees"),
                                               br(), br(),
                                               img(src="phenrica.png", align = "left", height="25%", width="25%"),
                                               br(), br(), 
                                               h3(strong("Plot phylogenies")),
                                               br(),
                                               
                                               selectInput("select_tree", label = "Select a tree to plot", choices = NULL),
                                               sliderInput("tip_label_size", label = "Tip label size: ", value = 0.8, min = 0.1, max = 2),
                                               sliderInput("support_value_size", label = "Support value size: ", value = 1, min = 0.1, max = 5),
                                               sliderInput("line_width", label = "Branch line width: ", value = 1, min = 0.1, max =5),
                                               br(),
                                               wellPanel(align = "justify",
                                                         h4(strong("TREE PLOT TWEAKS")),
                                                         br(),
                                               fluidRow(
                                                 column(width = 3,
                                               selectInput("support_value_type", label = "Support value type: ", choices = c("GMYC estimate"), selected = "GMYC estimate", width = "150px"),
                                                 ),
                                               column(width = 3,
                                               selectInput("support_value_col", label = "Support value colour: ", choices = c("grey", "lightblue", "salmon", "lightgreen", "lightyellow", "white"), selected = "lightgreen", width = "150px"),
                                               ),
                                               column(width = 3,
                                               selectInput("support_value_frame", label = "Support value frame:", choices = c("none", "circle", "rect"), selected = "rect", width = "150px"),
                                               ),
                                               column(width = 3,
                                               selectInput("branch_col", label = "Branch colour: ", choices = c("black", "blue", "lightblue", "red", "green", "orange"), selected = "blue", width = "150px"),
                                               ),
                                               ),
                                               br(),
                                               actionButton("plot_gmyc_tree", label = strong("Plot GMYC tree result"), style="color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")),
                                               ),
                                               br(),
                                              
                                               wellPanel(align = "justify",
                                                         h4(strong("DOWNLOAD")),
                                                         br(),
                                                         textInput("file_name_tree", "File name: ", "phylogeny"),
                                                         
                                                         fluidRow(
                                                           column(width = 3,
                                                                  selectInput("plot_format_tree", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  textInput("w_plot_tree", "Width: ", 20, width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  textInput("h_plot_tree", "Height: ", 15, width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  selectInput("unit_plot_tree", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  conditionalPanel(
                                                                    condition = "input.plot_format_tree == 'png'",
                                                                    textInput("res_plot_tree", "Res (dpi): ", 300), width = "150px")
                                                           ),
                                                         ),
                                               
                                               
                                               downloadButton("download_gmyc_tree", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ),
                                               h3(strong("Tree output:")),
                                               br(), br(),
                                               plotOutput("gmyc_tree", height = "1250px"),
                                               br(), br()
                                      ),
                                      
                                      tabPanel(strong("Percentage Matches"),
                                               br(), br(),
                                               img(src="neochetina.png", align = "left", height="15%", width="15%"),
                                               br(), br(), br(),
                                               h3(strong("View percentage matches")),
                                               strong("Compare your predefined grouping information to the assessments made by the GMYC"),
                                               br(), br(),
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("View GMYC species list")),
                                                         br(),
                                                         selectInput("select_tree_speclist", label = "Select a tree file", choices = NULL),
                                                         br(),
                                               actionButton("view_gmyc_spec", label = strong("View"), style="color: black; background-color: #A6C3C6; border-color: black", icon("edit")), 
                                               downloadButton("download_gmyc_spec", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ), 
                                               
                                               br(),
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("View matches")),
                                                         br(),
                                               actionButton("view_match_data", label = strong("View"), style="color: black; background-color: #A6C3C6; border-color: black", icon("edit")),
                                               downloadButton("download_match_data", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ), br(),
                                               
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("View Matches Summary")),
                                                         br(),
                                               actionButton("view_summary_match_data", label = strong("View"), style="color: black; background-color: #A6C3C6; border-color: black", icon("edit")),
                                               downloadButton("download_match_data_summary", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"), 
                                               ),br(),
                                               tableOutput("matches"),
                                               br(), br()
                                               
                                      ),
                                      
                                      tabPanel(strong("Plot Percentage matches"),
                                               br(), br(),
                                               img(src="lysathia.png", align = "left", height="15%", width="15%"),
                                               br(), br(), br(),
                                               h3(strong("Plot percentage matches")),
                                               br(), 
                                               wellPanel(align = "justify",
                                                         h4(strong("PLOT TWEAKS")),
                                                         br(),
                                               selectInput("plot_matches_type", label = "Select match type: ", choices = c("Including single-sample GMYC species"=1, "Excluding single-sample GMYC species"=2), width = "350px"),
                                               fluidRow(
                                               column(width = 3,
                                               selectInput("plot_matches_point_colours", "Point colour: ", choices = c("black", "blue", "red", "darkgreen"), width = "150px"),
                                               ),
                                               column(width = 3,
                                               selectInput("plot_matches_line_colour", "Line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen", "white"), width = "150px"),
                                               ),
                                               column(width = 3,
                                               selectInput("ggtheme_matches", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"], width = "150px"),
                                               ),
                                               ),
                                               br(),
                                               actionButton("plot_matches", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")),
                                               ),
                                               wellPanel(align = "justify",
                                                         h4(strong("DOWNLOAD")),
                                                         br(),
                                                         textInput("file_name_percmatches", "File name: ", "percentage_matches"),
                                                         
                                                         fluidRow(
                                                           column(width = 3,
                                                                  selectInput("plot_format_percmatches", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  textInput("w_plot_percmatches", "Width: ", 20, width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  textInput("h_plot_percmatches", "Height: ", 15, width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  selectInput("unit_plot_percmatches", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  conditionalPanel(
                                                                    condition = "input.plot_format_percmatches == 'png'",
                                                                    textInput("res_plot_percmatches", "Res (dpi): ", 300), width = "150px")
                                                           ),
                                                         ),
                                               
                                               downloadButton("download_match_plot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ),
                                               h3(strong("Plot:")),
                                               br(), br(),
                                               plotOutput("match_plot", height = "600px"),
                                               br(), br()
                                               
                                      ),
                                      
                                      tabPanel(strong("GMYC Oversplitting"),
                                               br(), br(),
                                               img(src="bactrocera.png", align = "left", height="15%", width="15%"),
                                               br(), br(), br(),
                                               h3(strong("GMYC oversplitting")),
                                               strong("Assess which taxa are split into multiple groups, indicating potential diversity"),
                                               br(),br(),
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("View output:")),
                                                         br(),
                                               actionButton("GMCY_oversplit_full_table", label = strong("View full"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               actionButton("GMYC_oversplit_table_view", label = strong("View summary"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(), br(),
                                               downloadButton("GMYC_oversplit_full_table_download", label = strong("Download full"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               downloadButton("GMYC_oversplit_table_download", label = strong("Download summary"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ),
                                               br(), 
                                               
                                               wellPanel(align = "justify",
                                               selectInput("GMYC_oversplit_ggtheme", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"], width = "150px"),
                                               ),
                                               br(), 
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("BOXPLOT")),
                                                         br(),
                                               actionButton("GMYC_oversplit_boxplot", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               hr(),
                                               h4(strong("DOWNLOAD")),
                                               br(), 
                                               textInput("file_name_oversplit", "File name: ", "gmyc_oversplitting_boxplot"),
                                              
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_oversplit", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_oversplit", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_oversplit", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_oversplit", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_oversplit == 'png'",
                                                          textInput("res_plot_oversplit", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("GMYC_oversplit_boxplot_download", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ),
                                               br(),
                                               
                                               wellPanel(align = "justify",
                                                h4(strong("BARPLOT")),
                                                br(),
                                                fluidRow(
                                                  column(width = 3,
                                                selectInput("GMYC_barchart_fill", "Fill: ", choices = c("black", "lightgrey", "white", "lightblue", "lightgreen"), selected = "white", width = "150px"),
                                                  ),
                                                column(width = 3,
                                                selectInput("GMYC_barchart_outline", "Outline: ", choices = c("black", "steelblue", "white", "darkgreen"), selected = "black", width = "150px"),
                                                ),
                                                ),         
                                               actionButton("GMYC_oversplit_barplot", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: darkgreen"),
                                               hr(),
                                               h4(strong("DOWNLOAD")),
                                               br(), 
                                               textInput("file_name_oversplitbar", "File name: ", "gmyc_oversplit_barplot"),
                                              
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_oversplitbar", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_oversplitbar", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_oversplitbar", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_oversplitbar", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_oversplitbar == 'png'",
                                                          textInput("res_plot_oversplitbar", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("GMYC_oversplit_barplot_download", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: darkgreen"),
                                               ),
                                               br(), 
                                               h3(strong("Table/Plot output:")),
                                               br(), br(),
                                               tableOutput("GMYC_oversplit_table"),
                                               br(), br(),
                                               plotOutput("GMYC_oversplit_plot"),
                                               br(), br()
                                               
                                      ),
                                      
                                      tabPanel(strong("Amalgamate"),
                                               br(), br(),
                                               img(src="tetramesa.png", align = "left", height="15%", width="15%"),
                                               br(), br(), br(),
                                               h3(strong("Amalgamate data")),
                                               br(), 
                                               wellPanel(align = "justify",
                                                         br(),
                                               fileInput("amalgamate_multiple", label = "Upload multiple .csv files with multiple columns of output data:", accept = ".csv", multiple = TRUE),
                                               ),
                                               
                                               selectInput("amalg_col", "Select column data:", choices=NULL),
                                               actionButton("view_amalg", strong("View"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               downloadButton("download_amalg", strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(), br(),
                                               tableOutput("amalgamate_table"),
                                               br(), br()
                                               
                                      ),
                                      
                                      tabPanel(strong("Plot for multiple-column data"),
                                               br(), br(),
                                               h3(strong("Plot multiple-column data")),
                                               br(), 
                                               wellPanel(align = "justify",
                                               selectInput("ggtheme_multiple", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"], width = "150px"),
                                               ),
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("BOXPLOT")),
                                                         br(),
                                                         actionButton("multiple_input_boxplot", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")),
                                                         br(), br(),
                                                         textInput("file_name_multi_box", "File name: ", "multiple_data_boxplot"),
                                                        
                                                         fluidRow(
                                                           column(width = 3,
                                                                  selectInput("plot_format_multi_box", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  textInput("w_plot_multi_box", "Width: ", 20, width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  textInput("h_plot_multi_box", "Height: ", 15, width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  selectInput("unit_plot_multi_box", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  conditionalPanel(
                                                                    condition = "input.plot_format_multi_box == 'png'",
                                                                    textInput("res_plot_multi_box", "Res (dpi): ", 300), width = "150px")
                                                           ),
                                                         ),
                                                         
                                                         downloadButton("download_multiple_input_boxplot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")),
                                               ),
                                               
                                               fluidRow(
                                                 h4(strong("LINE PLOT")),
                                                 column(6, wellPanel(align = "justify",
                                                        h4(strong("LINE 1")),
                                                        br(),
                                                        fileInput("multiple_input", label = "Upload a .csv file with multiple columns of output data:", accept = ".csv"),
                                                        
                                                        selectInput("multiple_input_line_type", label = 'Select a line type: ', choices = c("blank" = 0, "solid" = 1, "dashed" = 2, "dotted" = 3), selected = 1),
                                                        selectInput('multiple_input_line_col', label = "Select a line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen")),
                                                        sliderInput("multiple_input_line_width", label = "Line width: ", value = 1, min = 0.5, max = 5, step = 0.5),
                                                        sliderInput("multiple_input_point_size", label = "Point size: ", value = 1, min = 1, max = 10, step = 0.5),
                                                        selectInput("multiple_input_point_shape", label = "Point shape: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18)),
                                                        selectInput("multiple_input_point_colour", label = "Point colour: ", choices = c("blue", "lightblue", "red", "green", "forestgreen", "black", "grey"), selected = "black"),
                                                 ),
                                                 ), # wellpanel end
                                                 
                                                 
                                                 column(6, wellPanel(align = "justify",
                                                        h4(strong("LINE 2")),
                                                        br(),
                                                        fileInput("multiple_input2", label = "Upload a .csv file with multiple columns of output data:", accept = ".csv"),
                                                        
                                                        #selectInput("error_bar_type2", label = "Select what error bars should represent: ", choices = c("ci", "sd", "se")),
                                                        selectInput("multiple_input_line_type2", label = 'Select a line type: ', choices = c("blank" = 0, "solid" = 1, "dashed" = 2, "dotted" = 3), selected = 1),
                                                        selectInput('multiple_input_line_col2', label = "Select a line colour: ", choices = c("black", "grey", "lightblue", "salmon", "lightgreen")),
                                                        sliderInput("multiple_input_line_width2", label = "Line width: ", value = 1, min = 0.5, max = 5, step = 0.5),
                                                        sliderInput("multiple_input_point_size2", label = "Point size: ", value = 1, min = 1, max = 10, step = 0.5),
                                                        selectInput("multiple_input_point_shape2", label = "Point shape: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18)),
                                                        selectInput("multiple_input_point_colour2", label = "Point colour: ", choices = c("blue", "lightblue", "red", "green", "forestgreen", "black", "grey"), selected = "black"),
                                                 )
                                                 ),
                                                 ), 
                                              
                                               wellPanel(align = "justify",
                                               fluidRow(
                                               column(width = 3,
                                               selectInput("error_bar_type", label = "Error bar type: ", choices = c("ci", "sd", "se"), width = "150px"),
                                               ),
                                               column(width = 3,
                                               selectInput("multiple_input_error_bar_colour", label = "Error bar colour: ", choices = c("grey", "lightblue", "black", "white"), width = "150px"),
                                               ),
                                               ),
                                               textInput("title_multiple_input", label = "Title: ", value = "Title", width = "600px"),
                                               textInput("x_lab_multiple_input", label = "X-axis label: ", value = "Resampled data (%)", width = "600px"),
                                               textInput("y_lab_multiple_input", label = "Y-axis label: ", value = "Measure variable", width = "600px"),
                                               #textInput("y_interval_multiple_input", label = "Y-axis tick-mark interval: ", value = "10", width = "200px"),
                                               #numericInput("y_interval_multiple_input", label = "Y-axis tick-mark interval: ", value = 1, min = 1),
                                               checkboxInput("include_line2", label = strong("Include second line on line plot?"), value = FALSE),
                                               
                                               conditionalPanel(
                                                 condition = "input.include_line2",
                                                 textInput("line1_lab", "Line 1 label: ", value = "Line 1"),
                                                 textInput("line2_lab", "Line 2 label: ", value = "Line 2"),
                                               ),
                                               
                                               ),
                                               
                                               fluidRow(actionButton("plot_multiple_input", label = strong("Plot"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black", icon("drafting-compass")), align = "center"),
                                               br(), br(),
                                               
                                               wellPanel(align = "justify",
                                               h4(strong("DOWNLOAD")),
                                               br(),
                                               textInput("file_name_multi_line", "File name: ", "multiple_data_lineplot"),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_multi_line", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_multi_line", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_multi_line", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_multi_line", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_multi_line == 'png'",
                                                          textInput("res_plot_multi_line", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_multiple_input_plot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ),
                                               
                                               br(), br(),
                                               
                                               plotOutput("multiple_input_plot", height = "600px"),
                                               br(), br()
                                               
                                      )
                                    )
                                    
                                    
                                    ),
                                    
                           tabPanel("About",
                                    wellPanel(align = "justify",
                                              HTML("<h1 align = 'center'>SPEDE-SAMPLER <i>1.0.0</i> </h1>"),
                                              p(Sys.Date(), align = "center"),
                                              HTML("<p align = 'center'><img src = 'GitHub.png' width = '20px' height = 'auto'> <a target='_blank' rel='noopener noreferrer' href='https://github.com/CJMvS/spede-sampler'> GitHub Link </a></p>"),
                                              img(src="spede_authors.png", align = "center", height="75%", width="75%", style="display: block; margin-left: auto; margin-right: auto;"),
                                              br(), br(),
                                              strong("Image credits, with thanks to David Taylor at the Center for Biological Control, Department of Zoology and Entomology, Rhodes University:"),
                                              br(), br(),
                                              HTML("Random resampling tab: <i>Cornops aquaticum</i>"),
                                              br(),
                                              HTML("View Data tab: Iris flea beetle"),
                                              br(),
                                              HTML("Plot Results tab: Parasitoid wasp on <i>Eragrostis curvula</i>, Eastern Cape"),
                                              br(),
                                              HTML("Plot Trees tab: <i>Phenrica guerni</i>"),
                                              br(),
                                              HTML("Percentage Matches tab: <i>Neochetina</i>"),
                                              br(),
                                              HTML("Plot percentage matches tab: <i>Lysathia</i>"),
                                              br(),
                                              HTML("GMYC Oversplitting tab: <i>Bactrocera</i>"),
                                              br(),
                                              HTML("Amalgamate tab: <i>Tetramesa</i>")
                                              #HTML("<p><b>Cite the application:</b> van Steenderen, C.J.M. & Sutton, G.F.S. (2021) SPEDE-SAMPLER:  assess GMYC sampling effects on SPecies DElimitation </i>, 11: 1526-1534 <a target='_blank' rel='noopener noreferrer' href='https://doi.org/10.1002/ece3.6928'>https://doi.org/10.1002/ece3.6928</a></p>")
                                              
                                    )
                                           
                                    )
                ),
                useShinyalert() # set up the shinyalert package
          
                # titlePanel(h3("Run GMYC analyses on resampled ML trees for SPEcies DElimitation", )),
                # 
                # titlePanel(h4("Created by Clarke van Steenderen")),
                # br(),
                
                
                     
                
)

#)
