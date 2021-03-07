
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
    
    ################################################################################################
    # Read in an optional csv file with predefined grouping information
    ################################################################################################
    observe({
    
        infile = input$predefined_groups
        if (is.null(infile)) {
            predefined_groups_uploaded <<- FALSE
            return(NULL)
        } 
        
        groups <<- read.csv(infile$datapath)
        predefined_groups_uploaded <<- TRUE
        
        
        updateSelectInput(session,"col.group", choices=colnames(groups))
        updateSelectInput(session,"sample_names", choices=colnames(groups))
        
        
    }) # end of observe
    ################################################################################################
    
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
        
        if (ml_input_type == "FastTree") files = gtools::mixedsort( list.files(pattern = "\\.tre$") ) # the gtools::mixedsort() function preserves the original order of the files, instead of ordering them 1, 10, 100 etc.
        else if (ml_input_type == "RAxML") files = gtools::mixedsort( list.files(pattern = "bestTree") )
        
        # populate the dropdown list to show all the files inputted, so that the user can plot one if they wish
        updateSelectInput(session,"select_tree", choices=files)
        updateSelectInput(session,"select_tree_speclist", choices=files)
        
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
        
        # initialise the dataframes to store cluster and entity info
        clust_ent = data.frame(matrix(nrow = length(files), ncol = 3))
        colnames(clust_ent) = c("filename", "clusters", "entities")
        
        # this initialises the dataframe below, even if the user doesn't input grouping info. Doesn't seem to work when these three lines are within the if statemetn for group_info
        record = data.frame(matrix(nrow = length(files), ncol = 3))
        colnames(record) = c("filename", "percentage_match", "percentage_single_sample_GMYC_species")
        record$filename = files
        
        # create a list to which each gmyc tree is stored
        tree_container = vector("list", length(files))
        gmyc_spec_container = vector("list", length = length(files))
            
        ################################################################################################
        # Loop through the files in the chosen directory to create an ultrametric tree for each,
        # then run the GMYC analysis,
        # and store the entities and clusters
        ################################################################################################
        
        withProgress(message = 'Running GMYC', value = 0, {
          
        
        for(i in seq(along=files)) {
            
            treex = ape::read.tree(files[i])
            treex.ultra = chronos(treex, lambda = 0) # converts it to an ultrametric tree. This is an alternative to creating the ultrametric tree in BEAST first, and then running this GMYC analysis
            #tree.ultra = phytools::force.ultrametric(treex)
            treex.ultra2 = multi2di(treex.ultra, random = T) # makes the tree fully dichotomous
            
            if (input$set_seed == TRUE) set.seed(1234)
            
            # Run the GMYC analysis
            treex.gmyc = gmyc(treex.ultra2, quiet = F, method = "multiple")
            
            # store the gmyc in the tree_container listC:\Users\s1000334\Documents\PhD_Tetramesa\Experiment Fasta Files\Iterations_10\FastTree_output
            tree_container[[i]] = treex.gmyc

            # these three lines below write the files to the same folder:
            # uncomment to write them
            
            #sink(paste(files[i], c("gmyc_out.txt"), sep = ""))
            #summary.gmyc(treex.gmyc)
            #sink()
               
            clust_ent[i,1] = files[i] # populate file name
            clust_ent[i,2] = treex.gmyc$cluster[which.max(treex.gmyc$likelihood)] # extract the number of clusters
            clust_ent[i,3] = treex.gmyc$entity[which.max(treex.gmyc$likelihood)] # extract the number of entities
            
            ########################################################################################################
            # Do this bit only if the user uploaded a csv file with grouping data and checked the "YES" radio button
            ########################################################################################################
            
            if (predefined_groups_uploaded == TRUE){ #input$group_info == TRUE
                
                groups_col = as.name( input$col.group )
                sample_names_col = as.name (input$sample_names )
                
                gmyc.spec = splits::spec.list(treex.gmyc)
                gmyc.spec['ids'] = NA # add an empty column that will later have the predefined IDS as uploaded by the user
                
               
                for (j in (1:nrow(gmyc.spec))){
                    
                    for (k in (1:nrow(groups))) {
                        
                        if (gmyc.spec$sample_name[j] == groups[[sample_names_col]][k])
                            gmyc.spec$ids[j] = groups[[groups_col]][k]
                    }
                }
                
                gmyc_spec_container[[i]] = gmyc.spec
                
                # split the data by the species group numbers assigned by the GMYC
                
                spec_split = split(gmyc.spec, gmyc.spec$GMYC_spec)
                gmyc.spec$GMYC_spec = as.factor(gmyc.spec$GMYC_spec)
                
                matches.df = data.frame(matrix(nrow = length(levels(gmyc.spec$GMYC_spec)), ncol = 2)) 
                colnames(matches.df) = c("GMYC_spec", "match(y/n)") 
                matches.df$GMYC_spec = 1:nrow(matches.df) 
                
                y_count = 0
                n_count = 0
                single_sample_count = 0
                
                for (m in 1:length(spec_split)){
                    
                    len_unique = length(unique(spec_split[[m]]$ids))
                    len_values = length(spec_split[[m]]$ids)
                    
                    if (len_unique == 1 && len_values > 1){
                        matches.df$`match(y/n)`[m] = "y"
                        y_count = y_count + 1
                    }
                    
                    else if (len_values == 1){
                        matches.df$`match(y/n)`[m] = "single_sample"
                        single_sample_count = single_sample_count + 1
                    }
                    
                    else {
                        matches.df$`match(y/n)`[m] = "n"
                        n_count = n_count + 1
                    }
                }
                
                total = y_count + n_count 
                
                success = round(y_count/total * 100, 2) # this is only taking the matches for GMYC species groups that had more than one sample in it
                prop_single_samples = round(single_sample_count/nrow(matches.df) * 100, 2)
                
                record[i,2] = success
                record[i,3] = prop_single_samples
                
                
            }# end of if statement
           
            ################################################################################################
            
            ################################################################################################
            
            # Increment the progress bar, and update the detail text.
            incProgress(1/length(files), detail = paste("Processing tree file ", i, " of ", length(files)))
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            
        } # end of for loop
          
            
        }) # end of progress bar bracket
        
        
        shinyalert::shinyalert("Complete", "Please click on the tabs at the top of the window to view the results.", type = "success")
        
        ################################################################################################
        # View the GMYC species table with the predefined grouping information appended as the last column
        ################################################################################################
        
        observeEvent(input$view_gmyc_spec, {
            
            gmyc_spec_to_show = input$select_tree_speclist # get the tree file name selected from the drop down menu
            file_i = which(files == gmyc_spec_to_show) # get the index of that file to match up with the trees stored in the spec_list_container list
            
            output$matches = renderTable(gmyc_spec_container[[file_i]], rownames = FALSE, colnames = TRUE, digits = 0)
            
            
            ################################################################################################
            # Download GMYC species table
            ################################################################################################
            
            output$download_gmyc_spec = downloadHandler(
                
                filename = function (){paste('gmyc_spec_data', 'csv', sep = '.')},
                content = function (file){write.csv(gmyc_spec_container[[file_i]], file, row.names = FALSE)}
            )
            
        })
        
        ################################################################################################
        # View the match data for each file
        ################################################################################################
        
        observeEvent(input$view_match_data, {
            output$matches = renderTable(record, rownames = FALSE, colnames = TRUE, digits = 0)
            
            ################################################################################################
            # Download match data for each file
            ################################################################################################
            
            output$download_match_data = downloadHandler(
                
                filename = function (){paste('match_data', 'csv', sep = '.')},
                content = function (file){write.csv(record, file, row.names = FALSE)}
            )
            
        })
        
        ################################################################################################
        # Plot percentage match data for all files as a ggplot
        ################################################################################################
        
        observeEvent(input$plot_matches, {
            
            output$match_plot = renderPlot({
                ggplot(record, aes(x=1:nrow(record), y=percentage_match)) + 
                    geom_line(color = input$plot_matches_line_colour) + 
                    geom_point(color = input$plot_matches_point_colours)  + 
                    theme_classic() +
                    xlab("ML tree file/iteration number") + 
                    ylab("Percentage match") + 
                    ggtitle("GMCY species match to predefined groups vs ML tree/iteration file")
            })
        })
        
        
        ################################################################################################
        # view the match data summary for all files (mean, stdev etc.)
        ################################################################################################
        
        observeEvent(input$view_summary_match_data, {
            
            ################################################################################################
            # Get summary statistics for the success record and number of single samples assigned to a GMYC species
            ################################################################################################
            
            # stats for percentage matches
            avg_match = mean(record$percentage_match)
            stdev_match = sd(record$percentage_match)
            min_match = min(record$percentage_match)
            max_match = max(record$percentage_match)
            match_summary = round( c(avg_match, stdev_match, min_match, max_match) ,2 )
            
            # stats for percentage of single-ssample gmyc species
            avg_single = mean(record$percentage_single_sample_GMYC_species)
            stdev_single = sd(record$percentage_single_sample_GMYC_species)
            min_single = min(record$percentage_single_sample_GMYC_species)
            max_single = max(record$percentage_single_sample_GMYC_species)
            single_summary = round( c(avg_single, stdev_single, min_single, max_single) ,2 )
            
            # create a data frame housing summary statistics, called 'match_stats.df'
            match_stats.df = data.frame(matrix(nrow=4, ncol = 3)) 
            colnames(match_stats.df) = c("statistic", "percentage_matches", "percentage_single_sample_GMYC_species")
            categories = c("Average", "Standard deviation", "Minimum", "Maximum")
            match_stats.df$statistic = categories
            match_stats.df$percentage_matches = match_summary
            match_stats.df$percentage_single_sample_GMYC_species = single_summary
            
            output$matches =  renderTable(match_stats.df, rownames = FALSE, colnames = TRUE, digits = 2)
            
            ################################################################################################
            # Download summary match data for all files
            ################################################################################################
            
            output$download_match_data_summary = downloadHandler(
                
                filename = function (){paste('match_data_summary', 'csv', sep = '.')},
                content = function (file){write.csv(match_stats.df, file, row.names = FALSE)}
            )
            
            
        })
        
        
        observeEvent(input$all_data, {
            output$data_table =  renderTable(clust_ent, rownames = FALSE, colnames = TRUE, digits = 0)
        })
        
        observeEvent(input$summary_data, {
            
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
        # Plot a user-selected GMYC tree from the dropdown list
        ################################################################################################
        
        observeEvent(input$plot_gmyc_tree, {
            
            output$gmyc_tree = renderPlot({
                
            gmyc_to_plot = input$select_tree # get the tree file name selected from the drop down menu
            file_index = which(files == gmyc_to_plot) # get the index of that file to match up with the trees stored in the tree_container list
            
            tip_label_size = input$tip_label_size
            support_value_size = input$support_value_size
            branch_width = input$line_width
            support_value_col = input$support_value_col
            support_value_frame = input$support_value_frame
            branch_col = input$branch_col
                
            
            ####################################################################################################################################
            # Function to plot the gmyc as an alternative to plot_gmyc()
            # this allows you to change cex and colours, and also avoids the "hit enter to view next plot" issue
            # https://tmfujis.wordpress.com/2016/01/26/a-better-plot-function-for-gmyc-result/
            ####################################################################################################################################
            
            plot.result.gmyc <- function(res, cex=0.5, edge.width=1, no.margin=F, show.tip.label=T, label.offset=0) {
                plot.tree <- function(tr, mrca, cex=0.5, edge.width=1, no.margin=F, show.tip.label=T, label.offset=0) {
                    traverse.subtree <- function(tr, n=1) {
                        numtip <- length(tr$tip.label)
                        sub.node <- tr$edge[which(tr$edge[,1] == n+numtip), 2]
                        
                        res.node <- c()
                        for (sn in sub.node) {
                            res.node <- c(res.node, sn)
                            res.node <- c(res.node, traverse.subtree(tr, sn-numtip))
                        }
                        return (res.node)
                    }
                    
                    numtip <- length(tr$tip.label)
                    br.col <- rep(1, length(tr$edge.length))
                    
                    for (i in mrca) {
                        for (j in traverse.subtree(tr, i-numtip)) {
                            br.col[which(tr$edge[,2]==j)] <- branch_col
                            
                        }
                    }
                    
                    plot(tr, edge.color=br.col, show.tip.label=show.tip.label, cex=cex, edge.width=edge.width, no.margin=no.margin, label.offset=label.offset)
                }
                
                plot.tree(res$tree, res$MRCA[[which.max(res$likelihood)]]+length(res$tree$tip.label), cex=cex, edge.width=edge.width, no.margin=no.margin, show.tip.label=show.tip.label, label.offset=label.offset)
            }
            
            ####################################################################################################################################
            
            support = splits::gmyc.support(tree_container[[file_index]]) 
            is.na(support[support == 0]) = TRUE
            
            plot.result.gmyc(tree_container[[file_index]], cex = tip_label_size, edge.width = branch_width)
            nodelabels(round(support, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
            
            ################################################################################################
            # Download GMYC tree
            ################################################################################################
            
            
            output$download_gmyc_tree <- downloadHandler(
                
                filename = function (){paste(files[file_index], "svg", sep = '.')},
                
                content = function (file){svg(file)
                    plot.result.gmyc(tree_container[[file_index]], cex = tip_label_size, edge.width = branch_width)
                    nodelabels(round(support, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
                    dev.off()
                }
            )
            
            
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

    #})# end of first observe
    
}

