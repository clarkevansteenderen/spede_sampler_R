
library(shiny)


server = function(input, output, session) {
    
    # tell the shinyhelper package what the file name of the help file is
    observe_helpers(help_dir = "HelpFile")
   
    
    volumes =  getVolumes()() # this presents all the folders on the user's PC
    
    path = reactive({ 
        shinyDirChoose(input, 'directory', roots= volumes, session=session) 
        return(parseDirPath(volumes, input$directory))
        })
    
    output$folder_path = renderText( path() ) # quick check to see if the directory is stored as 'path' 
    
    
    observeEvent(input$run_gmyc, {
        
        if (length(path()) == 0) output$warn_files = renderText('<b>You have not selected a folder.')
        
        else{# else1
            
        setwd(path())
        # output$wd = renderPrint({getwd()}) # this checks if the wd has been set accordingly
        
        ml_input_type = input$data_type
        
        if (ml_input_type == "FastTree") files = list.files(pattern = "\\.tre$")
        else if (ml_input_type == "RAxML") files = list.files(pattern = "bestTree") 
        
        ################################################################################################
        # check whether there are files in the folder that contain bestTree or .tre in their names
        ################################################################################################
        
        if (length(files) == 0) output$warn_files = renderText('<b>There are no appropriate tree files in the folder you have selected.')
        
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
            treex.ultra = ape::chronopl(treex, 1) # converts it to an ultrametric tree. This is an alternative to creating the ultrametric tree in BEAST first, and then running
            # this GMYC analysis
            treex.ultra2 = multi2di(treex.ultra) # makes the tree fully dichotomous
            treex.gmyc = gmyc(treex.ultra2, quiet = T)

            # these three lines below write the files to the same folder:
            # uncomment to write them
            
            #sink(paste(files[i], c("gmyc_out.txt"), sep = ""))
            #summary.gmyc(treex.gmyc)
            #sink()
               
            clust_ent[i,1] = files[i]
            clust_ent[i,2] = treex.gmyc$cluster[which.max(treex.gmyc$likelihood)] # extract the number of clusters
            clust_ent[i,3] = treex.gmyc$entity[which.max(treex.gmyc$likelihood)] # extract the number of entities
            
            # Increment the progress bar, and update the detail text.
            incProgress(1/length(files), detail = paste("Processing tree file ", i))
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            
        } # end of for loop
          
            
        }) # end of progress bar bracket
        
        output$complete = renderText('<b>COMPLETE')
        
        
        ################################################################################################
        # Plot entities and clusters
        ################################################################################################
        
        observeEvent(input$plot_clusts, {
            
            output$clust_ent_plot = renderPlot({
                
                plot(clust_ent$clusters, clust_ent$entities, pch = 16, ylab = "Entities", xlab = "Clusters",  main = "Number of Clusters and Entities for GMYC")
                
            })
        })
        
        ################################################################################################
        # Download plot
        ################################################################################################
        
        output$download_clust_plot <- downloadHandler(
            
            filename = function (){paste("clust_ent_plot", "svg", sep = '.')},
            
            content = function (file){svg(file)
                plot(clust_ent$clusters, clust_ent$entities, pch = 16, ylab = "Entities", xlab = "Clusters",  main = "Number of Clusters and Entities for GMYC")
                dev.off()
            }
        )
        
        
        ################################################################################################
        # Download entities and clusters data
        ################################################################################################
        
        
        output$download_clust_ent_data = downloadHandler(
            
            filename = function (){paste('clust_ent_data', 'csv', sep = '.')},
            content = function (file){write.csv(clust_ent, file)}
        )
        
        ################################################################################################
        
        
        
        }# end of else 1
        
        }# end of else 2
        
    })
    
    
    
   
    
}

