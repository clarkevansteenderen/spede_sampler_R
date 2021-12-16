
mypackages <- c("ape", "remotes", "here", "gridExtra", "data.table", "dplyr", "shiny", "shinyWidgets", "ips", "Rmisc", "shinyhelper", "gtools", 
                "magrittr", "shinyFiles", "shinythemes", "shinyalert", "phytools", "reshape2", "devtools", "ggplot2", "tools", "devtools", "beastier", 
                "beautier", "tracerer", "rJava")

checkpkg <- mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(checkpkg)) install.packages(checkpkg, dependencies = TRUE)

install.packages("splits", repos="http://R-Forge.R-project.org")

library(devtools)
library(shiny)
library(remotes)
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
library(here)
library(gridExtra)
library(beastier)
library(tracerer)
library(beautier)
library(rJava)

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
  
                navbarPage(title = "SPEDE-sampler", id = "tabset",
                          tabPanel("Home",
                          wellPanel(align = "justify", style = "background: white",
                          br(), br(),                
                          div(img(src='spede_sampler_gmyc_logo.png',  height = 'auto', width = '70%'), style="display: block; margin-left: auto; margin-right: auto; text-align: center;"),
                          br(), br(),
                          fluidRow(actionButton("app", strong("BEGIN"), style = 'font-size:150%; color: black; background-color: #A6C3C6; border-color: black'), align = "center"),
                          br(), br(),
                          div(img(src='cover_insects.png',  height = 'auto', width = '80%'), style="display: block; margin-left: auto; margin-right: auto; text-align: center;"),
                          br(),
                          strong("Software created by: Clarke J.M. van Steenderen"),
                          br(),
                          strong("Contributing author to the associated publication: Guy F. Sutton"),
                          br(),
                          strong("Acknowledgements: The Centre for Biological Control (CBC) and Rhodes University for funding")
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
                                          textOutput("num_seqs"),
                                          br(),
                                          textInput("resampled_fasta_folder_name", label = 'Output fasta folder name: ', value = "RESAMPLED_FASTA_FILES"),
                                          br(),
                                          textInput("resampled_fasta_file_name", label = 'Output fasta file name: ', value = "resampled"),
                                          br(),
                                          fluidRow(
                                            column(width = 5,
                                          numericInput("fasta_subsample_percent", "Percentage to resample:", value = 10,  min = 1, step = 1, width = "250px"),
                                          ),
                                          column(width = 4,
                                          numericInput("fasta_resample_iterations", "Number of iterations: ", value = 2, min = 1, step = 1, width = "250px"),
                                          ),
                                          ),
                                          checkboxInput("set_seed_resampling", label = strong("Set a seed?"), value = TRUE),
                                          selectInput("resampling_approach", strong("Select resampling approach:"), choices = c("Random resampling", "Keep at least one representative sequence per predefined group"), width = "450px"),
                                          conditionalPanel(
                                            condition = "input.resampling_approach == 'Keep at least one representative sequence per predefined group'",
                                            br(),
                                            fileInput("resampling_groupings", label="Upload a .csv file with sequence names and morphospecies groups:", accept = ".csv"),
                                             actionButton("resampling_groupings_uploaded", strong("Confirm file"), style = 'font-size:120%; color: black; background-color: #A6C3C6; border-color: black', icon("thumbs-o-up")),
                                             br(), br(),
                                             selectInput("resampling_group_col", "Select Group Column:", choices="", width = "350px"),
                                             selectInput("resampling_sample_name_col", "Select Sample Name Column:", choices="", width = "350px")
                                          ),
                                          br(),
                                         fluidRow( actionButton("resample_fastas", strong("Resample"), style="font-size:150%; color: black; background-color: #A6C3C6; border-color: black"), align = "center")
                                          
                                        )       
                                        
                                      ),
                                      
                                      # tabPanel(strong("Genetic Divergence"),
                                      #   br(),
                                      #   img(src="fruit_piercing_moth.png", align = "left", height="20%", width="20%"),
                                      #   br(), br(),
                                      #   h3(strong("Genetic Divergence Estimates")),
                                      #   strong("Generate average within-group genetic divergence estimates across morphospecies"),
                                      #   br(), br(), br(),
                                      #   wellPanel(
                                      #     fileInput("genetic_divergence_groupings", label="Upload a .csv file with sequence names and morphospecies groups:", accept = ".csv"),
                                      #     actionButton("genetic_divergence_groupings_uploaded", strong("Confirm file"), style = 'font-size:120%; color: black; background-color: #A6C3C6; border-color: black', icon("thumbs-o-up")),
                                      #     br(), br(),
                                      #     selectInput("genetic_divergence_group_col", "Select Group Column:", choices="", width = "350px"),
                                      #     selectInput("genetic_divergence_sample_name_col", "Select Sample Name Column:", choices="", width = "350px"),
                                      #     br(),
                                      #     textInput(inputId = 'fasta_file_path_genetic_divergence', label = 'Manually insert a file path to the folder containing your Fasta files:'),
                                      #     br(),
                                      #     actionButton("calculate_genetic_divergences", strong("Calculate"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                      #     actionButton("view_genetic_divergences", strong("View"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                      #     
                                      #     br(), br(),
                                      #     downloadButton("download_genetic_dists", strong('Download Data'), style="color: black; background-color: #A6C3C6; border-color: black")
                                      #   ),
                                      #   tableOutput("genetic_divergence_table"),
                                      #   
                                      #   # wellPanel(
                                      #   #   strong("Genetic divergence per morphospecies"),
                                      #   #   br(), br(),
                                      #   #   downloadButton("download_divergences_per_group", strong('Download Data'), style="color: black; background-color: #A6C3C6; border-color: black")
                                      #   # )
                                      # ),
                                      
                                      tabPanel(strong("BEAST XML Files"),
                                               br(),
                                               img(src="beast.png", align = "left", height="10%", width="10%"),
                                               br(), br(),
                                               h3(strong("Create XML files")),
                                               strong("Use functionality from BEAUti to create an input .xml file for a BEAST analysis"),
                                               br(), br(), br(),
                                               wellPanel(
                                                 textInput(inputId = 'resampled_fasta_file_path', label = 'Manually insert a file path to the folder containing your resampled Fasta files:'),
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
                                               strong("Obtain maximum clade credibility trees, and specify a burnin value"),
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
                                                            selectInput("tree_tool", "Select phylogenetic method:", choices = c("BEAST"), width = "250px"),
                                                  
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
                                                                     checkboxInput("set_seed", label = strong("Set a seed?"), value = TRUE), style = "font-size: 16px",
                                                                     # checkboxInput("download_all_results", strong("Download all results upon completion?"), value = FALSE), style = "font-size: 16px",
                                                                     # conditionalPanel(condition = "input.download_all_results == true",
                                                                     #                  textInput("all_results_folder_name", "Folder name to which all results are written to:"),
                                                                     # ),
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
                                               br(), br(),
                                               #strong("These GMYC species are recorded as merges/undersplits:"),
                                               htmlOutput("merges"),
                                               tags$head(tags$style("#merges{color: darkcyan;
                                                       font-size: 16px;
                                                       }"
                                                          )
                                                      ),
                                               
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
                                      
                                      tabPanel(strong("Plot Percentage Matches"),
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
                                      
                                      tabPanel(strong("GMYC Splitting"),
                                               br(), br(),
                                               img(src="bactrocera.png", align = "left", height="15%", width="15%"),
                                               br(), br(), br(),
                                               h3(strong("GMYC splitting")),
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
                                               textInput("file_name_oversplit", "File name: ", "gmyc_splitting_boxplot"),
                                              
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
                                      
                                      tabPanel(strong("GMYC Exact Matches"),
                                               br(), br(),
                                               img(src="cornops_nymph.png", align = "left", height="20%", width="20%"),
                                               br(), br(), br(),
                                               h3(strong("GMYC Exact Matches")),
                                               strong("Assess which taxa are exact GMYC matches"),
                                               br(),br(),
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("View output:")),
                                                         br(), 
                                                         htmlOutput("percent_exact_matches"),
                                                         tags$head(tags$style("#percent_exact_matches{color: darkcyan;
                                                             font-size: 16px;
                                                             }"
                                                         )
                                                         ),
                                                         br(), br(),
                                                         actionButton("GMCY_exact_match_full_table", label = strong("View full"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                                         actionButton("GMYC_exact_match_table_view", label = strong("View summary"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                                         br(), br(),
                                                         downloadButton("GMYC_exact_match_full_table_download", label = strong("Download full"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                                         downloadButton("GMYC_exact_match_table_download", label = strong("Download summary"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               ),
                                               br(), 
                                               
                                               wellPanel(align = "justify",
                                                         selectInput("GMYC_exact_match_ggtheme", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"], width = "150px"),
                                               ),
                                               br(), 
                                               
                                               wellPanel(align = "justify",
                                                         h4(strong("BARPLOT")),
                                                         br(),
                                                         fluidRow(
                                                           column(width = 3,
                                                                  selectInput("exact_match_barchart_fill", "Fill: ", choices = c("black", "lightgrey", "white", "lightblue", "lightgreen"), selected = "white", width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  selectInput("exact_match_barchart_outline", "Outline: ", choices = c("black", "steelblue", "white", "darkgreen"), selected = "black", width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  numericInput("x_axis_angle_exact_matches", "x-axis label angle:", value = 0, min = 0, step = 1, width = "150px"),
                                                                  )
                                                         ),  
                                                         
                                                         actionButton("GMYC_exact_match_barplot", label = strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: darkgreen"),
                                                         hr(),
                                                         h4(strong("DOWNLOAD")),
                                                         br(), 
                                                         textInput("file_name_exact_match_bar", "File name: ", "gmyc_exact_match_barplot"),
                                                         
                                                         fluidRow(
                                                           column(width = 3,
                                                                  selectInput("plot_format_exact_match", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  textInput("w_plot_exact_match", "Width: ", 20, width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  textInput("h_plot_exact_match", "Height: ", 15, width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  selectInput("unit_plot_exact_match", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                           ),
                                                           column(width = 3,
                                                                  conditionalPanel(
                                                                    condition = "input.plot_format_exact_match == 'png'",
                                                                    textInput("res_plot_exact_match", "Res (dpi): ", 300), width = "150px")
                                                           ),
                                                         ),
                                                         
                                                         downloadButton("GMYC_exact_match_barplot_download", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: darkgreen"),
                                               ),
                                               br(), 
                                               h3(strong("Table/Plot output:")),
                                               br(), br(),
                                               tableOutput("GMYC_exact_match_table"),
                                               br(), br(),
                                               plotOutput("GMYC_exact_match_plot"),
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
                                      
                                      tabPanel(strong("Summary Plots"),
                                               br(), br(),
                                               img(src="hydrellia.png", align = "left", height="15%", width="15%"),
                                               br(), br(), br(),
                                               h3(strong("Summary Plots")),
                                               br(), br(),
                                               selectInput("ggtheme_summary_plots", "Select ggplot Theme:", choices = names(ggthemes), selected = ggthemes["Classic"], width = "250px"),
                                               
                                               h3(strong("Clusters and Entities")),
                                               wellPanel(
                                                 fluidRow(
                                                   column(width = 4,
                                                   fileInput("summary_clusters", "Cluster data"),
                                                       ),
                                                   column(width = 4,
                                                   fileInput("summary_entities", "Entities data"),
                                                   ),
                                                   column(width = 4,
                                                   numericInput("no_morphospecies", "Number of morphospecies:", value = 0, min = 0, step = 1, width = "250px"),
                                                   ),
                                                   
                                               ), # end of fluidrown
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                  selectInput("clust_ent_summary_line_col", "Line colour:", choices = c("black", "blue", "brown", "darkblue", "darkgreen", "forestgreen", "grey", "lightblue", "orange", "pink", "purple", "red", "turquoise", "white"), selected = "lightblue"),      
                                                 ),
                                                 column(width = 3,
                                                  selectInput("clust_ent_summary_ci_col", "CI band colour:", choices = c("black", "grey", "lightblue", "salmon", "lightgreen", "white", "yellow"), selected = "grey"),      
                                                 ),
                                                 column(width = 3,
                                                  numericInput("clust_ent_summary_ci_alpha", "Alpha value:", value = 0.2, min = 0, step = 0.1),      
                                                 ),
                                                 column(width = 3,
                                                  numericInput("clust_ent_summary_y_interval", "Y-axis interval:", value = 10, min = 1, step = 1),      
                                                 )
                                               ),
                                               
                                               actionButton("plot_summary_clusts", strong("Plot clusters"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               actionButton("plot_summary_ents", strong("Plot entities"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(), br(),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_summary_clusts_ents", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_summary_clusts_ents", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_summary_clusts_ents", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_summary_clusts_ents", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_summary_clusts_ents == 'png'",
                                                          textInput("res_plot_summary_clusts_ents", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_clust_summary_plot", label = strong("Download cluster plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               downloadButton("download_ent_summary_plot", label = strong("Download entities plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br()
                                               
                                               ), # end of wellpanel
                                               
                                               
                                               h3(strong("Splitting Ratios")),
                                               wellPanel(
                                                 fluidRow(
                                                   
                                               column(width = 6,
                                               fileInput("summary_splitting_ratio_inc", "Splitting ratio data including singletons"),
                                               ),
                                               column(width = 6,
                                                      fileInput("summary_splitting_ratio_exc", "Splitting ratio data excluding singletons"),
                                               ),
                                               ), # end of fluidrow
                                               
                                               fluidRow(
                                                 column(width = 6,
                                                        selectInput("splitting_inc_col", "Splitting ratio including singletons colour:", choices = c("grey", "hotpink", "lightblue", "lightyellow", "salmon", "lightgreen", "red", "royalblue", "white"), selected = "lightblue"),      
                                                 ),
                                                 column(width = 6,
                                                        selectInput("splitting_exc_col", "Splitting ratio excluding singletons colour:", choices = c("grey", "hotpink", "lightblue", "lightyellow", "salmon", "lightgreen", "red", "royalblue", "white"), selected = "lightyellow"),      
                                                 ),
                                                 column(width = 4,
                                                        selectInput("splitting_plot_point_shape", label = "Point shape of means: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18), selected = 17),  
                                                        ),
                                                 column(width = 4,
                                                        selectInput("splitting_plot_point_col", label = "Colour of mean points: ", choices = c("black", "white", "red", "grey", "blue", "darkblue")),  
                                                 ),
                                                 column(width = 3,
                                                        sliderInput("splitting_plot_point_size", label = "Point size: ", value = 3, min = 0, max = 10, step = 0.5),
                                                 ),
                                                 
                                               ),
                                               
                                               fluidRow(
                                                 column(width = 4,
                                               numericInput("splitting_plot_legend_x_pos", "Legend position x coordinate:", value = 0.9, min = 0, max = 1, step = 0.1, width = "250px"), 
                                                 ),
                                               column(width = 4,
                                                      numericInput("splitting_plot_legend_y_pos", "Legend position y coordinate:", value = 0.1, min = 0, max = 1, step = 0.1, width = "250px"), 
                                                 ),
                                               ),
                                               
                                               actionButton("plot_summary_splits", strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(), br(),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_summary_splits", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_summary_splits", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_summary_splits", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_summary_splits", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_summary_splits == 'png'",
                                                          textInput("res_plot_summary_splits", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_summary_splits_plot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(),
                                               
                                               ), # end of wellpanel
                                               
                                               h3(strong("Percentage Matches")),
                                               wellPanel(
                                                 fluidRow(
                                                   column(width = 6,
                                               fileInput("summary_percentage_match_inc", "Percentage match data including singletons"),
                                                   ),
                                               column(width = 6,
                                               fileInput("summary_percentage_match_exc", "Percentage match data excluding singletons"),
                                               ),
                                               
                                               ),# end of fluidrow
                                               
                                               fluidRow(
                                                 column(width = 6,
                                                        selectInput("percentage_match_inc_col", "Percentage match including singletons colour:", choices = c("grey", "hotpink", "lightblue", "lightyellow", "salmon", "lightgreen", "red", "royalblue", "white"), selected = "lightblue"),      
                                                 ),
                                                 column(width = 6,
                                                        selectInput("percentage_match_exc_col", "Percentage match excluding singletons colour:", choices = c("grey", "hotpink", "lightblue", "lightyellow", "salmon", "lightgreen", "red", "royalblue", "white"), selected = "lightyellow"),      
                                                 ),
                                                 column(width = 4,
                                                        selectInput("percentage_match_plot_point_shape", label = "Point shape of means: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18), selected = 17),  
                                                 ),
                                                 column(width = 4,
                                                        selectInput("percentage_match_plot_point_col", label = "Colour of mean points: ", choices = c("black", "white", "red", "grey", "blue", "darkblue")),  
                                                 ),
                                                 column(width = 3,
                                                        sliderInput("percentage_match_plot_point_size", label = "Point size: ", value = 3, min = 0, max = 10, step = 0.5),
                                                 ),
                                                 
                                               ),
                                               
                                               fluidRow(
                                                 column(width = 4,
                                                        numericInput("percentage_match_plot_legend_x_pos", "Legend position x coordinate:", value = 0.9, min = 0, max = 1, step = 0.1, width = "250px"), 
                                                 ),
                                                 column(width = 4,
                                                        numericInput("percentage_match_plot_legend_y_pos", "Legend position y coordinate:", value = 0.1, min = 0, max = 1, step = 0.1, width = "250px"), 
                                                 ),
                                               ),
                                               
                                               actionButton("plot_summary_percentage_matches", strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               
                                               br(), br(),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_summary_percentage_matches", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_summary_percentage_matches", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_summary_percentage_matches", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_summary_percentage_matches", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_summary_percentage_matches == 'png'",
                                                          textInput("res_plot_summary_percentage_matches", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_summary_percentage_matches_plot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(),
                                               
                                               ), # end of wellpanel
                                               
                                               h3(strong("Percentage Singletons")),
                                               wellPanel(
                                               fileInput("summary_percentage_singletons", "Percentage singletons data"),
                                               
                                               fluidRow(
                                                 
                                                 column(width = 3,
                                               selectInput("percentage_singletons_col", "Barplot colour:", choices = c("grey", "lightblue", "lightgreen", "lightyellow", "salmon", "red", "black", "royalblue", "forestgreen", "white"), selected = "lightblue",  width = "250px"), 
                                                 ),
                                               column(width = 3,
                                               selectInput("percentage_singletons_plot_point_shape", label = "Point shape of means: ", choices = c("Round filled" = 16, "Round open" = 1, "+" = 3, "X" = 4, "Square" = 15, "Triangle" = 17, "Diamond" = 18), selected = 17, width = "250px"),
                                               ),
                                               column(width = 3,
                                               selectInput("percentage_singletons_plot_point_col", label = "Colour of mean points: ", choices = c("black", "white", "red", "grey", "blue", "darkblue"), width = "250px"),  
                                               ),
                                               column(width = 3,
                                                      sliderInput("percentage_singletons_plot_point_size", label = "Point size: ", value = 3, min = 0, max = 10, step = 0.5),
                                               ),
                                               
                                               ), # end of fluidrow
                                               
                                               actionButton("plot_summary_singletons", strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               
                                               br(), br(),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_summary_percentage_singletons", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_summary_percentage_singletons", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_summary_percentage_singletons", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_summary_percentage_singletons", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_summary_percentage_singletons == 'png'",
                                                          textInput("res_plot_summary_percentage_singletons", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_summary_percentage_singletons_plot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(),
                                               
                                               ),
                                               
                                               h3(strong("Splitting Ratios per Morphospecies Group")),
                                               wellPanel(
                                               fileInput("summary_splitting_morphospecies", "Mean splitting data per morphospecies"),
                                               selectInput("splitting_morphospecies_col", "Barplot colour:", choices = c("grey", "lightblue", "royalblue", "black", "red", "lightyellow", "salmon", "lightgreen", "white"), selected = "lightblue",  width = "250px"), 
                                               numericInput("x_axis_angle", "x-axis label angle:", value = 0, min = 0, step = 1, width = "150px"),
                                               actionButton("plot_summary_oversplits_per_morphospecies", strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               
                                               br(), br(),
                                               
                                               fluidRow(
                                                 column(width = 3,
                                                        selectInput("plot_format_summary_splitting_morphospecies", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("w_plot_summary_splitting_morphospecies", "Width: ", 20, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        textInput("h_plot_summary_splitting_morphospecies", "Height: ", 15, width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        selectInput("unit_plot_summary_splitting_morphospecies", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                 ),
                                                 column(width = 3,
                                                        conditionalPanel(
                                                          condition = "input.plot_format_summary_splitting_morphospecies == 'png'",
                                                          textInput("res_plot_summary_splitting_morphospecies", "Res (dpi): ", 300), width = "150px")
                                                 ),
                                               ),
                                               
                                               downloadButton("download_summary_splitting_morphospecies_plot", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               br(),
                                               
                                               ),
                                               
                                               # h3(strong("Genetic Divergences")),
                                               # wellPanel(
                                               #   
                                               #    fileInput("genetic_divergence_data", "Genetic Divergence Data"),
                                               #  
                                               #   fluidRow(
                                               #     column(width = 3,
                                               #            selectInput("genetic_divergence_line_col", "Line colour:", choices = c("black", "blue", "brown", "darkblue", "darkgreen", "forestgreen", "grey", "lightblue", "orange", "pink", "purple", "red", "turquoise", "white"), selected = "lightblue"),      
                                               #     ),
                                               #     column(width = 3,
                                               #            selectInput("genetic_divergence_ci_col", "CI band colour:", choices = c("black", "grey", "lightblue", "salmon", "lightgreen", "white", "yellow"), selected = "grey"),      
                                               #     ),
                                               #     column(width = 3,
                                               #            numericInput("genetic_divergence_ci_alpha", "Alpha value:", value = 0.2, min = 0, step = 0.1),      
                                               #     ),
                                               #     column(width = 3,
                                               #            numericInput("genetic_divergence_y_interval", "Y-axis interval:", value = 0.005, min = 1, step = 1),      
                                               #     )
                                               #   ),
                                               #   
                                               #   textInput("genetic_divergence_y_label", "y-axis label:", value = "Mean genetic divergence"),
                                               #   br(),
                                               #   actionButton("plot_genetic_divergence", strong("Plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               #   
                                               #   br(), br(),
                                               #   
                                               #   fluidRow(
                                               #     column(width = 3,
                                               #            selectInput("plot_format_genetic_divergence", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                               #     ),
                                               #     column(width = 3,
                                               #            textInput("w_plot_genetic_divergence", "Width: ", 20, width = "150px"),
                                               #     ),
                                               #     column(width = 3,
                                               #            textInput("h_plot_genetic_divergence", "Height: ", 15, width = "150px"),
                                               #     ),
                                               #     column(width = 3,
                                               #            selectInput("unit_plot_genetic_divergence", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                               #     ),
                                               #     column(width = 3,
                                               #            conditionalPanel(
                                               #              condition = "input.plot_format_genetic_divergence == 'png'",
                                               #              textInput("res_plot_genetic_divergence", "Res (dpi): ", 300), width = "150px")
                                               #     ),
                                               #   ),
                                               #   
                                               #   downloadButton("download_genetic_divergence_plot", label = strong("Download plot"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                               #   br()
                                               #   
                                               # ), # end of wellpanel
                                               
                                               h3(strong("Multi-Plots")),
                                               wellPanel(
                                                 numericInput("gridextra_ncol", "Number of columns for plots:", value = 2, min = 1, step = 1, width = "250px"),
                                                 br(), 
                                                 actionButton("plot_multiple_summary_plots", strong("Plot All"),  style="color: black; background-color: #A6C3C6; border-color: black"),
                                                 
                                                 br(), br(),
                                                 
                                                 fluidRow(
                                                   column(width = 3,
                                                          selectInput("plot_format_multiple_summary_plots", "Image format:", choices = c("pdf", "png", "svg"), width = "150px"),
                                                   ),
                                                   column(width = 3,
                                                          textInput("w_plot_multiple_summary_plots", "Width: ", 20, width = "150px"),
                                                   ),
                                                   column(width = 3,
                                                          textInput("h_plot_multiple_summary_plots", "Height: ", 15, width = "150px"),
                                                   ),
                                                   column(width = 3,
                                                          selectInput("unit_multiple_summary_plots", "Unit: ", choices=c("cm", "in"), width = "150px"),
                                                   ),
                                                   column(width = 3,
                                                          conditionalPanel(
                                                            condition = "input.plot_format_multiple_summary_plots == 'png'",
                                                            textInput("res_plot_multiple_summary_plots", "Res (dpi): ", 300), width = "150px")
                                                   ),
                                                 ),
                                                 
                                                 downloadButton("download_multiple_summary_plots", label = strong("Download"), style="color: black; background-color: #A6C3C6; border-color: black"),
                                                 br(),
                                                 
                                                 ), # end of wellpanel
                                               
                                               br(), br(),
                                               
                                               plotOutput("summary_plot", height = "650px")
                                               
                                      )
                                      
                                    )
                                    
                                    
                                    ),
                                    
                           tabPanel("About",
                                    wellPanel(align = "justify",
                                              HTML("<h1 align = 'center'>SPEDE-SAMPLER <i>1.0.0</i> </h1>"),
                                              p(Sys.Date(), align = "center"),
                                              HTML("<p align = 'center'><img src = 'GitHub.png' width = '20px' height = 'auto'> <a target='_blank' rel='noopener noreferrer' href='https://github.com/CJMvS/spede-sampler'> GitHub Link </a></p>"),
                                              #img(src="clarke.png", align = "center", height="30%", width="30%", style="display: block; margin-left: auto; margin-right: auto;"),
                                              br(), br(),
                                              strong("Clarke van Steenderen is a PhD student, and Guy Sutton is a post-doctoral researcher with the the Centre for Biological Control in the Department of Zoology and Entomology at Rhodes University, South Africa"),
                                              br(), br(),
                                              strong("SPEDE-sampler is currently under review. Please cite as: van Steenderen, C.J.M., and Sutton, G.F. 2021. SPEDE-sampler: an R Shiny application to assess how methodological choices and taxon-sampling can affect GMYC output and interpretation."),
                                              br(), br(),
                                              strong("Image credits with thanks to David Taylor at the Centre for Biological Control (CBC):"),
                                              br(), br(),
                                              HTML("Random resampling tab: <i>Cornops aquaticum</i> adult"),
                                              br(),
                                              HTML("View Data tab: Iris flea beetle"),
                                              br(),
                                              HTML("Plot Results tab: Parasitoid wasp on <i>Eragrostis curvula</i>"),
                                              br(),
                                              HTML("Plot Trees tab: <i>Phenrica guerni</i>"),
                                              br(),
                                              HTML("Percentage Matches tab: <i>Neochetina</i>"),
                                              br(),
                                              HTML("Plot percentage matches tab: <i>Lysathia</i>"),
                                              br(),
                                              HTML("GMYC splitting tab: <i>Bactrocera</i>"),
                                              br(),
                                              HTML("GMYC Exact Matches tab: <i>Cornops aquaticum</i> nymph"),
                                              br(),
                                              HTML("Amalgamate tab: <i>Tetramesa</i>"),
                                              br(),
                                              HTML("Summary Plots tab: <i>Hydrellia</i>")
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
