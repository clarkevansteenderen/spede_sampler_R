
library(shiny)


server = function(input, output, session) {
    
    # tell the shinyhelper package what the file name of the help file is
    observe_helpers(help_dir = "HelpFile")
   
    
    volumes =  getVolumes()() # this presents all the folders on the user's PC
    
    path = reactive({ 
        shinyDirChoose(input, 'directory', roots= volumes, session=session) 
        return(parseDirPath(volumes, input$directory))
        })
    
    output$folder_path = renderText( path() )  # quick check to see if the directory is stored as 'path'
    
    
    observeEvent(input$run_gmyc, {
        
        manual_file_path = input$raw_file_path
        
        ################################################################################################
        # check whether the user has neglected to select a folder or input a file path
        ################################################################################################
        
        if (length(path()) == 0 && manual_file_path == "") shinyalert::shinyalert("No folder or path input", "You have not selected a folder or pasted in a folder path.", type = "error")
        else if (length(path()) > 0 && manual_file_path != "") shinyalert::shinyalert("Issue: Two folder inputs", "You have selected a folder and inputted a file path. Please remove the manual file path and ensure that the desired folder has been selected.", type = "error")
        
        ################################################################################################
        # If they have inputted one or the other, proceed:
        ################################################################################################
        
        else{# else1
            
        output$warn_files = renderText("")
        
        ################################################################################################
        # check whether the user has inserted only a manual file path. 
        # If so, set that path as the working directory
        ################################################################################################
            
        if (length(path()) == 0 && manual_file_path != "") {
                
            path = manual_file_path
            setwd(path)
        }
            
        ################################################################################################
        # Otherwise use the path stored from the folder input button
        ################################################################################################
            
        else setwd(path())
        # output$wd = renderPrint({getwd()}) # this checks if the wd has been set accordingly
        
        ml_input_type = input$data_type
        
        if (ml_input_type == "FastTree") files = list.files(pattern = "\\.tre$")
        else if (ml_input_type == "RAxML") files = list.files(pattern = "bestTree") 
        
        ################################################################################################
        # check whether there are files in the folder that contain bestTree or .tre in their names
        ################################################################################################
        
        #if (length(files) == 0) output$warn_files = renderText('<b>There are no appropriate tree files in the folder you have selected.')
        if (length(files) == 0) shinyalert::shinyalert("No .tre or bestTree files", "There are no appropriate tree files in the folder you have selected.", type = "error")

        ################################################################################################
        # If there are, then proceed:
        ################################################################################################
        
        else{ # else2
            
        output$warn_files = renderText("") # clear the warning if there was one before
        
        clust_ent = data.frame(matrix(nrow = length(files), ncol = 3))
        colnames(clust_ent) = c("filename", "clusters", "entities")
            
        ################################################################################################
        # Loop through the files in the chosen directory to create an ultrametric tree for each,
        # then run the GMYC analysis,
        # and store the entities and clusters
        ################################################################################################
        
        withProgress(message = 'Running GMYC', value = 0, {
          
        
        for(i in seq(along=files)) {
            
            treex = ape::read.tree(files[i])
            treex.ultra = phytools::force.ultrametric(treex) # converts it to an ultrametric tree. This is an alternative to creating the ultrametric tree in BEAST first, and then running this GMYC analysis
            #tree.ultra = phytools::force.ultrametric(treex)
            treex.ultra2 = multi2di(treex.ultra, random = T) # makes the tree fully dichotomous
            treex.gmyc = gmyc(treex.ultra2, quiet = F, method = "multiple")

            # these three lines below write the files to the same folder:
            # uncomment to write them
            
            #sink(paste(files[i], c("gmyc_out.txt"), sep = ""))
            #summary.gmyc(treex.gmyc)
            #sink()
               
            clust_ent[i,1] = files[i] # populate file name
            clust_ent[i,2] = treex.gmyc$cluster[which.max(treex.gmyc$likelihood)] # extract the number of clusters
            clust_ent[i,3] = treex.gmyc$entity[which.max(treex.gmyc$likelihood)] # extract the number of entities
            
            
            
            # Increment the progress bar, and update the detail text.
            incProgress(1/length(files), detail = paste("Processing tree file ", i, " of ", length(files)))
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            
        } # end of for loop
          
            
        }) # end of progress bar bracket
        
        #output$complete = renderText('<b>COMPLETE. Please click on the tabs at the top of the window to view the results.')
        shinyalert::shinyalert("Complete", "Please click on the tabs at the top of the window to view the results.", type = "success")
        
        ################################################################################################
        # Get summary statistics for the number of clusters and entities
        ################################################################################################
        
        # stats for clusters
        avg_clust = mean(clust_ent$clusters)
        stdev_clust = sd(clust_ent$clusters)
        min_clust = min(clust_ent$clusters)
        max_clust = max(clust_ent$clusters)
        clust_summary = round( c(avg_clust, stdev_clust, min_clust, max_clust) ,2 )
        
        # stats for entities
        avg_ent = mean(clust_ent$entities)
        stdev_ent = sd(clust_ent$entities)
        min_ent = min(clust_ent$entities)
        max_ent = max(clust_ent$entities)
        ent_summmary = round( c(avg_ent, stdev_ent, min_ent, max_ent) ,2)
        
        # create a data frame housing summary statistics, called 'stats.df'
        stats.df = data.frame(matrix(nrow=4, ncol = 3)) 
        colnames(stats.df) = c("statistic", "clusters", "entities")
        categories = c("Average", "Standard deviation", "Minimum", "Maximum")
        stats.df$statistic = categories
        stats.df$clusters = clust_summary
        stats.df$entities = ent_summmary
        
        
        observeEvent(input$all_data, {
            output$data_table =  renderTable(clust_ent, rownames = FALSE, colnames = TRUE, digits = 0)
        })
        
        observeEvent(input$summary_data, {
            output$data_table =  renderTable(stats.df, rownames = FALSE, colnames = TRUE, digits = 2)
        })
       
        
        ################################################################################################
        # Plot entities vs clusters
        ################################################################################################
        
        observeEvent(input$plot_clusts, {
            
            output$clust_ent_plot = renderPlot({
                
                plot(clust_ent$clusters, 
                     clust_ent$entities, 
                     pch = 16, 
                     ylab = "Entities",
                     xlab = "Clusters",  
                     main = "Number of Clusters vs Entities for GMYC",
                     col=input$clust_vs_ent_plot_point_colours)
                
            })
        })
        
        ################################################################################################
        # Download clusters vs entities plot
        ################################################################################################
        
        output$download_clust_plot <- downloadHandler(
            
            filename = function (){paste("clust_ent_plot", "svg", sep = '.')},
            
            content = function (file){svg(file)
                plot(clust_ent$clusters, 
                     clust_ent$entities, 
                     pch = 16, 
                     ylab = "Entities", 
                     xlab = "Clusters",  
                     main = "Number of Clusters vs Entities for GMYC",
                     col=input$clust_vs_ent_plot_point_colours)
                dev.off()
            }
        )
        
        ################################################################################################
        # Plot entities and clusters boxplot
        ################################################################################################
        
        observeEvent(input$plot_boxplot, {
            
            gg_clust_ent = reshape2::melt(clust_ent) # get into the format the ggplot can work with
            colnames(gg_clust_ent) = c("filename", "gmyc_cat", "count")
            gg_clust_ent$gmyc_cat %>% as.factor()
            
            output$clust_ent_plot = renderPlot({
                
                ggplot(gg_clust_ent, aes(x=gmyc_cat, y=count)) + geom_boxplot() +
                    theme_classic() + 
                    xlab("GMYC Category") + 
                    ylab("Count") 
            })
        })
        
        ################################################################################################
        # Download entities and clusters boxplot
        ################################################################################################
        
        output$download_boxplot <- downloadHandler(
            
            filename = function (){paste("clust_ent_boxplot", "svg", sep = '.')},
            
            content = function (file){
                ggsave(file, ggplot( gg_clust_ent, aes(x=gmyc_cat, y=count)) + geom_boxplot() +
                    theme_classic() + 
                    xlab("GMYC Category") + 
                    ylab("Count") )
            }
        )
        
        ################################################################################################
        # Plot clusters vs iteration/ML file number
        ################################################################################################
        
        observeEvent(input$plot_clusts_vs_iterations, {
            
            output$clust_ent_plot = renderPlot({
            ggplot(clust_ent, aes(x=1:nrow(clust_ent), y=clusters)) + 
                geom_line(color = input$plot_clusts_vs_iterations_line_colour) + 
                geom_point(color = input$plot_clusts_vs_iterations_point_colours)  + 
                theme_classic() +
                xlab("ML tree file/iteration number") + 
                ylab("Number of GMYC clusters") + 
                ggtitle("GMCY number of clusters vs ML tree/iteration file")
            })
        })
        
        ################################################################################################
        # Download clusters vs iteration/ML file number plot
        ################################################################################################
        
        output$download_clusts_vs_iterations <- downloadHandler(
            
            filename = function (){paste("clusters_vs_iteration", "svg", sep = '.')},
            
            content = function (file){
                ggsave(file, ggplot( clust_ent, aes(x=1:nrow(clust_ent), y=clusters)) + 
                    geom_line(color = input$plot_clusts_vs_iterations_line_colour) + 
                    geom_point(color = input$plot_clusts_vs_iterations_point_colours)  + 
                    theme_classic() +
                    xlab("ML tree file/iteration number") + 
                    ylab("Number of GMYC clusters") + 
                    ggtitle("GMCY number of clusters vs ML tree/iteration file") )
            }
        )
        
        ################################################################################################
        # Plot entities vs iteration/ML file number
        ################################################################################################
        
        observeEvent(input$plot_ents_vs_iterations, {
            
            output$clust_ent_plot = renderPlot({
                ggplot(clust_ent, aes(x=1:nrow(clust_ent), y=entities)) + 
                    geom_line(color = input$plot_ents_vs_iterations_line_colour) + 
                    geom_point(color = input$plot_ents_vs_iterations_point_colours)  + 
                    theme_classic() +
                    xlab("ML tree file/iteration number") + 
                    ylab("Number of GMYC entities") + 
                    ggtitle("GMCY number of entities vs ML tree/iteration file")
            })
        })
        
        ################################################################################################
        # Download entities vs iteration/ML file number plot
        ################################################################################################
        
        output$download_ents_vs_iterations <- downloadHandler(
            
            filename = function (){paste("entities_vs_iterations", "svg", sep = '.')},
            
            content = function (file){
                ggsave(file, ggplot( clust_ent, aes(x=1:nrow(clust_ent), y=entities)) + 
                    geom_line(color = input$plot_ents_vs_iterations_line_colour) + 
                    geom_point(color = input$plot_ents_vs_iterations_point_colours)  + 
                    theme_classic() +
                    xlab("ML tree file/iteration number") + 
                    ylab("Number of GMYC entities") + 
                    ggtitle("GMCY number of entities vs ML tree/iteration file") )
            }
        )
        
        ################################################################################################
        # Download entities and clusters data
        ################################################################################################
        
        
        output$download_clust_ent_data = downloadHandler(
            
            filename = function (){paste('clust_ent_data', 'csv', sep = '.')},
            content = function (file){write.csv(clust_ent, file, row.names = FALSE)}
        )
        
        
        ################################################################################################
        # Download stat summary for entities and clusters data
        ################################################################################################
        
        
        output$download_stat_summary = downloadHandler(
            
            filename = function (){paste('clust_ent_summary', 'csv', sep = '.')},
            content = function (file){write.csv(stats.df, file, row.names = FALSE)}
        )
        
        ################################################################################################
        
        
        
        }# end of else 1
        
        }# end of else 2
        
    })
    
    
    
   
    
}

