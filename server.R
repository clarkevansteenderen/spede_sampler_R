
ggthemes = list("Classic" = theme_classic(),
                "Dark" = theme_dark(),
                "Minimal" = theme_minimal(),
                "Grey" = theme_grey(),
                "Light" = theme_light(),
                "Black/White" = theme_bw(),
                "Void" = theme_void())

server = function(input, output, session) {
  
  home_dir = getwd()
  treeannotator_dir = paste(home_dir, "/TreeAnnotator.exe", sep = "")
  logcombiner_dir = paste(home_dir, "/LogCombiner.exe", sep = "")
  
  #access to the app from the homepage link
  observeEvent(input$app, updateTabsetPanel(session = session, inputId = "tabset", selected = "app"))
  
################################################################################################
# RANDOM RESAMPLING OF FASTA FILE SEQUENCES
################################################################################################
  
  # optionally read in an excel file of sequence names and predefined groups:
  
  groups_resampling = reactive({
    
    infile = input$resampling_groupings
    
    if (is.null(infile)) {
      return(NULL)
    }
    
    else{
      resample_groups <<- read.csv(infile$datapath)
      return(resample_groups)
    }
    
  }) # end of reactive
  
  # select which column is the grouping and which is the sample names
  observeEvent(input$resampling_groupings_uploaded, {
    updateSelectInput(session,"resampling_group_col", choices=colnames(groups_resampling()))
    updateSelectInput(session,"resampling_sample_name_col", choices=colnames(groups_resampling()))
  })
  
  
  # Read in the fasta file:
  
observe({
    
    in_seqs = input$fasta_file
    if (is.null(in_seqs)) {
      return(NULL)
    }
    
  seqs = ape::read.FASTA(in_seqs$datapath) 
  
  # get the number of sequences in the alignment file
  output$num_seqs = renderText({paste("There are ", length(seqs), " sequences in your alignment.")})
  
  
observeEvent(input$resample_fastas, {
  
  if(input$resampling_approach == "Random resampling"){
  
        # get the number of sequences to resample
        num_fasta_seqs = floor( (input$fasta_subsample_percent/100) * length(seqs) )
        
        # create an empty folder to store resampled fasta files
        dir.create(input$resampled_fasta_folder_name)
        
        # ITERATE THIS N TIMES
        out_file = paste(input$resampled_fasta_folder_name, "/", input$resampled_fasta_file_name, sep = "")
        
       if(input$set_seed_resampling == TRUE) set.seed(123)
        
        # function from the FastaUtils github package, Guillem Salazar
        fasta.sample<-function(infile=NULL,nseq=NULL,file.out=NULL,replacement=FALSE){
          seqs<- Biostrings::readDNAStringSet(infile)
          selected<-seqs[sample(1:length(seqs),nseq,replace=replacement)]
          Biostrings::writeXStringSet(selected,filepath=file.out)}
        
        withProgress(message = 'Resampling...', value = 0, {
          
        for (i in 1:input$fasta_resample_iterations){
          
          fasta.sample(infile = in_seqs$datapath, nseq = num_fasta_seqs, file.out = paste(out_file, "_", i, ".fasta", sep = ""), 
                                   replacement = FALSE)
        
          # Increment the progress bar, and update the detail text.
          incProgress(1/input$fasta_resample_iterations, detail = paste("Resampling file ", i, " of ", input$fasta_resample_iterations, "[ ", round(i/input$fasta_resample_iterations*100, 0), "% ]"))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          
        }
          
        }) # end of progressbar
        
        shinyalert::shinyalert("Complete", "Successfully resampled", type = "success")
  
   } # end of if statement re input$resampling_approach
  
    else{ # if the user wishes to keep at least one sequence representative per predefined group
      
      if(input$set_seed_resampling == TRUE) set.seed(123)
      
      id_col = as.name(input$resampling_sample_name_col)
      grp_col = as.name(input$resampling_group_col)
      groups_resampling_df = groups_resampling()
      groups_resampling_df[[grp_col]] = as.factor( groups_resampling_df[[grp_col]] )
      
      num_fasta_seqs = floor( (input$fasta_subsample_percent/100) * length(seqs) )
      
      num_groups = length( levels( groups_resampling_df[[grp_col]] ) )
      
      if(num_fasta_seqs < num_groups) shinyalert::shinyalert("Warning", paste("Please select a larger percentage to resample. You need at least ", ceiling( num_groups/length(seqs) * 100 ), "%"), type = "warning")
      
      else{
              out_folder = paste(input$resampled_fasta_folder_name)
              dir.create(out_folder)
              
            withProgress(message = 'Resampling...', value = 0, { 
              
            for (k in 1:input$fasta_resample_iterations){
                  
                  # subset each predefined group, and store them all in a list called "subset_groups"
                  subset_groups = c()
                  
                  for(i in levels( groups_resampling_df[[grp_col]] )){
                    subset_groups[[i]] = subset(groups_resampling_df, groups_resampling_df[[grp_col]] == i)
                  }
                  
                  # randomly select one sequence from each predefined group, and store it in a list called "extracted_samples"
                  extracted_samples = c()
                  
                  for(j in 1:length(subset_groups)){
                    rand_ind = sample(nrow(subset_groups[[j]]), 1, replace = FALSE) 
                    extracted_samples[[j]] = subset_groups[[j]][rand_ind,]
                    subset_groups[[j]] = subset_groups[[j]][-rand_ind, ]
                  }
                  
                  # convert these lists into dataframes
                  extracted_samples = dplyr::bind_rows(lapply(extracted_samples, as.data.frame.list))
                  subset_groups = dplyr::bind_rows(lapply(subset_groups, as.data.frame.list))
                  
                  # get the number of sequences to resample, taking into account the sequences that were already extracted to serve as representative sequences for each group
                  num = num_fasta_seqs - nrow(extracted_samples)
                  
                  # randomly select from the remaining sequences
                  # get random indices (or one index, depending):
                  resampled_ind = sample(nrow(subset_groups), num, replace = FALSE) 
                  # subset the associated rows based on those indices
                  resampled = subset_groups[resampled_ind,] 
                  # bind the randomly selected sequences to the ones that were extracted in the beginning (one representative sequence for each group)
                  final = rbind(extracted_samples, resampled) 
                  # extract the actual sequences in the fasta file based on the randomly selected names
                  seqs_subsetted = subset(seqs, labels(seqs) %in% final[[id_col]]) 
                  # write the fasta file
                  write.dna(seqs_subsetted, paste(out_folder, "/", input$resampled_fasta_file_name, "_", k, ".fasta", sep = ""), format = "fasta")
                  
                  # Increment the progress bar, and update the detail text.
                  incProgress(1/input$fasta_resample_iterations, detail = paste("Resampling file ", k, " of ", input$fasta_resample_iterations, "[ ", round(k/input$fasta_resample_iterations*100, 0), "% ]"))
                  # Pause for 0.1 seconds to simulate a long computation.
                  Sys.sleep(0.1)
            
            } # end of for loop
              
            }) # end of progressbar
            
            shinyalert::shinyalert("Complete", "Successfully resampled", type = "success")
        
      } # end of else
      
    } # end of else
  
}) # end of observeEvent
  
}) # end of observe
  

  ################################################################################################
  # GENETIC DIVERGENCE CALCULATION
  ################################################################################################
  
  # read in the groupings file
  # groups_genetic_divergence = reactive({
  #   
  #   infile = input$genetic_divergence_groupings
  #   
  #   if (is.null(infile)) {
  #     return(NULL)
  #   }
  #   
  #   else{
  #     genetic_divergence_groups <<- read.csv(infile$datapath)
  #     return(genetic_divergence_groups)
  #   }
  #   
  # }) # end of reactive
  # 
  # # select which column is the grouping and which is the sample names
  # observeEvent(input$genetic_divergence_groupings_uploaded, {
  #   updateSelectInput(session,"genetic_divergence_group_col", choices=colnames(groups_genetic_divergence()))
  #   updateSelectInput(session,"genetic_divergence_sample_name_col", choices=colnames(groups_genetic_divergence()))
  # })
  # 
  # # input the file path to the fasta files
  # 
  # observeEvent(input$calculate_genetic_divergences,{
  # 
  # fasta_file_path_gen_diverg = input$fasta_file_path_genetic_divergence
  # 
  # if (!dir.exists(fasta_file_path_gen_diverg)) 
  #   shinyalert::shinyalert("Error", "File path does not exist", type = "error")
  # 
  # else{
  #   
  #   setwd(fasta_file_path_gen_diverg)
  #   gen_diverg_fasta_files = gtools::mixedsort( list.files(pattern = "\\.fas") ) 
  #   w_group_dists = c()
  #   b_group_dists = c()
  #   o_group_dists = c()
  #   
  #   spp_dists = list()
  #   
  #   seq_name_col = as.name(input$genetic_divergence_sample_name_col)
  #   morphogroup_col = as.name(input$genetic_divergence_group_col)
  #   
  #   withProgress(message = 'Calculating genetic divergences...', value = 0, {
  #     
  #     for(i in seq(along = gen_diverg_fasta_files)){
  #       
  #       target = ape::read.FASTA(gen_diverg_fasta_files[i])
  #       dists = ape::dist.dna(target)
  #       
  #       seq_names = names(target)
  #       
  #       new_dat = data.frame(matrix(nrow = length(seq_names), ncol =2))
  #       colnames(new_dat) = c("seq_name", "group")
  #       new_dat$seq_name = seq_names
  #       
  #       for(m in (1:length(seq_names))){
  #         for(n in (1:nrow(genetic_divergence_groups))){
  #           if(seq_names[m] == genetic_divergence_groups[[seq_name_col]][n])
  #             new_dat$group[m] = genetic_divergence_groups[[morphogroup_col]][n]
  #         }
  #       }
  #       
  #       tryCatch({
  #         
  #       MD = vegan::meandist(dists, new_dat$group)
  #       MD_summary = summary(MD)
  #       w_group_dists[i] = MD_summary$W
  #       b_group_dists[i] = MD_summary$B
  #       o_group_dists[i] = MD_summary$D
  #       
  #       MD_df = as.data.frame(MD)
  #       vals = c()
  #       
  #       for(k in 1:nrow(MD_df)){
  #         vals[k] = MD_df[k,k]
  #       }
  #       
  #       spp_dists[[i]] = vals
  #       
  #       }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #       
  #       # Increment the progress bar, and update the detail text.
  #       incProgress(1/length(gen_diverg_fasta_files), detail = paste("Processing Fasta file ", i, " of ", length(gen_diverg_fasta_files), "[ ", round(i/length(gen_diverg_fasta_files)*100, 0), "% ]"))
  #       # Pause for 0.1 seconds to simulate a long computation.
  #       Sys.sleep(0.1)
  #       
  #     } # end of for loop
  #     
  #     group_dists_df = data.frame(matrix(nrow = length(gen_diverg_fasta_files), ncol = 4))
  #     colnames(group_dists_df) = c("filename", "intra_dist", "inter_dist", "overall_dist")
  #     group_dists_df$filename = gen_diverg_fasta_files
  #     group_dists_df$intra_dist = w_group_dists
  #     group_dists_df$inter_dist = b_group_dists
  #     group_dists_df$overall_dist = o_group_dists
  #     
  #     # tryCatch({
  #     # spp_dists = as.data.frame(do.call(rbind, spp_dists))
  #     # colnames(spp_dists) = levels(as.factor(genetic_divergence_groups[[morphogroup_col]]))
  #     # rownames(spp_dists) = gen_diverg_fasta_files 
  #     # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #     
  #   }) # end of progress bar
  #   
  #   shinyalert::shinyalert("Complete", "Divergences Successfully Calculated", type = "success")
  #   
  #   
  # } # end of else
  # 
  # observeEvent(input$view_genetic_divergences, {
  #   output$genetic_divergence_table = renderTable(group_dists_df, rownames = TRUE, colnames = TRUE, digits = 3)
  # })
  # 
  # 
  # output$download_genetic_dists = downloadHandler(
  #   
  #   filename = function (){paste('genetic_dists', 'csv', sep = '.')},
  #   content = function (file){write.csv(group_dists_df, file, row.names = FALSE)}
  # )
  # 
  # # output$download_divergences_per_group = downloadHandler(
  # #   
  # #   filename = function (){paste('genetic_dists_per_morphogroup', 'csv', sep = '.')},
  # #   content = function (file){write.csv(spp_dists, file, row.names = FALSE)}
  # # )
  # 
  # }) # end of observeEvent
  
################################################################################################
  # GENERATE XML FILES FOR BEAST
################################################################################################

  observeEvent(input$create_xml_files, {
  
    resampled_fasta_file_path = input$resampled_fasta_file_path
    
    if (!dir.exists(resampled_fasta_file_path)) 
      shinyalert::shinyalert("Error", "File path does not exist", type = "error")
    
    else{
        setwd(resampled_fasta_file_path)
        resampled_fasta_files = gtools::mixedsort( list.files(pattern = "\\.fas") ) 
        
        # different rate distribution options
        rate_option = function(distribution_prior, n){ 
          
          if(missing(n)){
          switch(distribution_prior, 
                 "Beta" = beautier::create_beta_distr(),
                 "Exponential" = beautier::create_exp_distr(),
                 "Gamma" = beautier::create_gamma_distr(),
                 "Inverse gamma" = beautier::create_inv_gamma_distr(),
                 "Laplace" = beautier::create_laplace_distr(),
                 "Log-normal" = beautier::create_log_normal_distr(),
                 "Normal" = beautier::create_normal_distr(),
                 "1/X" = beautier::create_distr_one_div_x(),
                 "Poisson" = beautier::create_poisson_distr(),
                 "Uniform" = beautier::create_uniform_distr())
          }
          
          else {
            switch(distribution_prior, 
                   "Beta" = beautier::create_beta_distr(value = n),
                   "Exponential" = beautier::create_exp_distr(value = n),
                   "Gamma" = beautier::create_gamma_distr(value = n),
                   "Inverse gamma" = beautier::create_inv_gamma_distr(value = n),
                   "Laplace" = beautier::create_laplace_distr(value = n),
                   "Log-normal" = beautier::create_log_normal_distr(value = n),
                   "Normal" = beautier::create_normal_distr(value = n),
                   "1/X" = beautier::create_distr_one_div_x(value = n),
                   "Poisson" = beautier::create_poisson_distr(value = n),
                   "Uniform" = beautier::create_uniform_distr(value = n))
          } 
      
    }
    
    
        # Set up the inference model 
        inference_model = beautier::create_inference_model(
          
          site_model = switch(input$beast_site_model, "GTR" = beautier::create_gtr_site_model(), 
                              "HKY" = beautier::create_hky_site_model(),
                              "JC69" = beautier::create_jc69_site_model(),
                              "TN93" = beautier::create_tn93_site_model()),
          
          clock_model = switch(input$beast_clock_model, "Strict" = beautier::create_strict_clock_model(clock_rate_param = input$beast_clock_rate), 
                               "Relaxed lognormal" = beautier::create_rln_clock_model(clock_rate_param = input$beast_clock_rate)),
         
          tree_prior = switch(input$beast_tree_prior, "Birth-death" = beautier::create_bd_tree_prior(birth_rate_distr = rate_option(input$distr_b), death_rate_distr = rate_option(input$distr_d)),
                              "Coalescent Bayesian skyline" = beautier::create_cbs_tree_prior(group_sizes_dimension = input$cbs_group_sizes_dim),
                              "Coalescent constant-population" = beautier::create_ccp_tree_prior(pop_size_distr = rate_option(input$distr_ccp, n = input$pop_size_distribution)),
                              "Coalescent exponential-population" = beautier::create_cep_tree_prior(pop_size_distr = rate_option(input$distr_cep_pop), growth_rate_distr = rate_option(input$distr_cep_gr)),
                              "Yule" = beautier::create_yule_tree_prior(birth_rate_distr = rate_option(input$distr_yule))),
          
          mcmc = beautier::create_mcmc(chain_length = input$beast_mcmc, store_every = input$beast_store_every)
    
        )
        
        # create xml files
        withProgress(message = 'Creating XML files...', value = 0, {
          
        for(i in seq(along = resampled_fasta_files)){
          
          beautier::create_beast2_input_file_from_model(
            input_filename =  resampled_fasta_files[i],
            inference_model = inference_model,
            output_filename = paste(input$xml_folder_name,"/", tools::file_path_sans_ext(resampled_fasta_files[i]), ".xml", sep = "" )
          )
          
          # Increment the progress bar, and update the detail text.
          incProgress(1/length(resampled_fasta_files), detail = paste("Processing fasta file ", i, " of ", length(resampled_fasta_files), "[ ", round(i/length(resampled_fasta_files)*100, 0), "% ]"))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          
        }
          
        }) # end of progress bar
        
        shinyalert::shinyalert("Complete", "Successfully created XML files", type = "success")
  
    } # end of else regarding setwd()
  
})
  
################################################################################################
# RUN BEAST ON THE GENERATED XML FILES
################################################################################################

  observeEvent(input$run_beast, {
    
    if (!dir.exists(input$xml_file_path)) 
      shinyalert::shinyalert("Error", "File path does not exist", type = "error")
    
    else{
    
    setwd(input$xml_file_path)
    
    xml_files = gtools::mixedsort( list.files(pattern = "\\.xml") ) 
    
    withProgress(message = 'Running BEAST...', value = 0, {
    for(i in seq(along = xml_files)){
      beastier::run_beast2(
        input_filename = xml_files[i],
        verbose = TRUE,
        use_beagle = input$run_beagle
        
      )
      
      # Increment the progress bar, and update the detail text.
      incProgress(1/length(xml_files), detail = paste("Processing XML file ", i, " of ", length(xml_files), "[ ", round(i/length(xml_files)*100, 0), "% ]"))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
    }
      
    }) # end of progress bar
    
    shinyalert::shinyalert("Complete", "BEAST successfully run", type = "success")
    
    } # end of else regarding setwd()
    
  })

################################################################################################
  # RUN TREEANNOTATOR ON BEAST .TREES FILES
################################################################################################
  
  observeEvent(input$run_treeannotator, {
    
    if (!dir.exists(input$beast_trees_file_path)) 
      shinyalert::shinyalert("Error", "File path does not exist", type = "error")
    
    else{
    
          file.copy(treeannotator_dir, input$beast_trees_file_path)
          
          setwd(input$beast_trees_file_path)
          
          beast_tree_files = gtools::mixedsort( list.files(pattern = "\\.trees") ) 
          
          burnin = paste("-burnin", input$treeannotator_burnin) 
          heights = paste("-heights", input$treeannotator_heights)
          
          withProgress(message = 'Running TreeAnnotator...', value = 0, {
            
          for(i in seq(along = beast_tree_files)){
            
            in_file = beast_tree_files[i]
            out_file = paste(tools::file_path_sans_ext(beast_tree_files[i]), ".nex", sep = "")
            
            if(input$treeannotator_low_mem == TRUE)
              treeannotator_command_line = paste("TreeAnnotator.exe", burnin, heights, "-lowMem", in_file, out_file)
            else treeannotator_command_line = paste("TreeAnnotator.exe", burnin, heights, in_file, out_file)
            
            system(treeannotator_command_line)
            
            # Increment the progress bar, and update the detail text.
            incProgress(1/length(beast_tree_files), detail = paste("Processing BEAST .trees file ", i, " of ", length(beast_tree_files), "[ ", round(i/length(beast_tree_files)*100, 0), "% ]"))
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            
            # this while loop delays the main loop so that it only goes into the next iteration after the file has been written.
            # this seems to prevent the program from freezing up when treeannotator is run with very large file sizes
            
            while (!file.exists( paste(tools::file_path_sans_ext( beast_tree_files[i] ), ".nex", sep = "") ) ) {
              Sys.sleep(2)
            }
            
          }
            
          }) # end of progress bar
          
          shinyalert::shinyalert("Complete", "TreeAnnotator successfully run", type = "success")
    
    } # end of else regarding setwd()
    
  })
  
################################################################################################
  # Run TRACER
################################################################################################
 
 observeEvent(input$load_log_files, {
   
   if (!dir.exists(input$beast_log_files_path)) 
     shinyalert::shinyalert("Error", "File path does not exist", type = "error")
   
   else{
   
         setwd(input$beast_log_files_path)
         
         log_files = gtools::mixedsort( list.files(pattern = "\\.log") )
         updateSelectInput(session,"select_log_file", choices=log_files)
         
         observeEvent(input$table_log_files, {
           
           estimates <- tracerer::parse_beast_tracelog_file(input$select_log_file)
           estimates <- tracerer::remove_burn_ins(estimates, burn_in_fraction = input$tracer_burnin)
           esses <- tracerer::calc_esses(estimates, sample_interval = input$tracer_sample_interval)
           ess_table <- t(esses)
           colnames(ess_table) <- c("ESS")
           output$tracer_ess = renderTable(ess_table, rownames = TRUE, colnames = TRUE, digits = 0)
           
          output$tracer_plot = renderPlot({  ggplot2::ggplot(
             data = tracerer::remove_burn_ins(estimates, burn_in_fraction = input$tracer_burnin),
             ggplot2::aes(x = Sample)
           ) + ggplot2::geom_line(ggplot2::aes(y = posterior)) + ggplot2::ggtitle("Trace plot") + ggplot2::theme_classic()
          })
          
         })
   
   } # end of else
   
 })
  
  ################################################################################################
  # RUN LOGCOMBINER
  ################################################################################################
  
  observeEvent(input$run_logcombiner, {
    
    if (!dir.exists(input$logcombiner_file_path)) 
      shinyalert::shinyalert("Error", "File path does not exist", type = "error")
    
    else{
    
          file.copy(logcombiner_dir, input$logcombiner_file_path)
          
          setwd(input$logcombiner_file_path)
          
          # create an empty folder 
          dir.create(input$logcombiner_folder_results)
          
          logcombiner_tree_files = gtools::mixedsort( list.files(pattern = "\\.trees") ) 
          
          resampling = paste("-resample", input$resample_freq) 
          
          withProgress(message = 'Running LogCombiner...', value = 0, {
            
            for(i in seq(along = logcombiner_tree_files)){
              
              logcombiner_in_file = paste("-log", logcombiner_tree_files[i])
              logcombiner_out_file = paste("-o", paste( input$logcombiner_folder_results, "/", tools::file_path_sans_ext(logcombiner_tree_files[i]), ".trees", sep = "") )
              
              logcombiner_command_line = paste("LogCombiner.exe", logcombiner_in_file, logcombiner_out_file, resampling)
              
              system(logcombiner_command_line)
              
              # Increment the progress bar, and update the detail text.
              incProgress(1/length(logcombiner_tree_files), detail = paste("Processing LogCombiner .trees file ", i, " of ", length(logcombiner_tree_files), "[ ", round(i/length(logcombiner_tree_files)*100, 0), "% ]"))
              # Pause for 0.1 seconds to simulate a long computation.
              Sys.sleep(0.1)
              
            }
            
          }) # end of progress bar
          
          shinyalert::shinyalert("Complete", "LogCombiner successfully run", type = "success")
    
    } # end of else re. setwd()
    
  })
  
################################################################################################

################################################################################################
  
  # get the filepath for the PATHD8 executable ready, if the user chooses to select it. It is in the PATHD8 folder, in the GitHub repo
  # that is downloaded.
  pathd8_filepath = paste(getwd(), "/PATHd8/PATHd8", sep="")
  
  # tell the shinyhelper package what the file name of the help file is
  # observe_helpers(help_dir = "HelpFile")
  
  # I'm leaving this bit in, despite taking the option of a folder upload via a button out of the GUI
  
  volumes =  getVolumes()() # this presents all the folders on the user's PC
  
  path = reactive({ 
    shinyDirChoose(input, 'directory', roots= volumes, session=session) 
    return(parseDirPath(volumes, input$directory))
  })
  
  output$folder_path = renderText( path() )  # quick check to see if the directory is stored as 'path'
  
  ################################################################################################
  # Read in an optional csv file with predefined grouping information for the GMYC analysis
  ################################################################################################
  
  
  predefined_groups_uploaded = reactive({

    infile = input$predefined_groups
    
    if (is.null(infile)) {
      return(NULL)
    }
    
    else{
      groups = read.csv(infile$datapath)
      return(groups)
    }
    
  }) # end of observe of input of predefined group file
      
    ################################################################################################
    ################################################################################################
    
  observeEvent(input$group_data_uploaded, {
    updateSelectInput(session,"col.group", choices=colnames(predefined_groups_uploaded()))
    updateSelectInput(session,"sample_names", choices=colnames(predefined_groups_uploaded()))
  })
  
    ################################################################################################
    # Run the GMYC analysis
    ################################################################################################
    
    observeEvent(input$run_gmyc, {
      
      if (!dir.exists(input$raw_file_path)) 
        shinyalert::shinyalert("Error", "File path does not exist", type = "error")
      
      else{
      
            manual_file_path = input$raw_file_path
            
            ################################################################################################
            # check whether the user has neglected to select a folder or input a file path
            ################################################################################################
            
            if (length(path()) == 0 && manual_file_path == "") shinyalert::shinyalert("No folder or path input", "You have not selected a folder or pasted in a folder path.", type = "warning")
            if (length(path()) > 0 && manual_file_path != "") shinyalert::shinyalert("Issue: Two folder inputs", "You have selected a folder and inputted a file path. Please remove the manual file path and ensure that the desired folder has been selected.", type = "warning")
            if (!is.null( predefined_groups_uploaded() ) && input$col.group == "" ) shinyalert::shinyalert("Confirm columns", "Please click the Confirm button, and then select columns for grouping and sample name informaiton for your .csv file data.", type = "warning")
            
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
              
              if (input$tree_tool == "BEAST") files = gtools::mixedsort( list.files(pattern = "\\.nex") )
              
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
                record = data.frame(matrix(nrow = length(files), ncol = 6))
                colnames(record) = c("filename", "percentage_match", "percentage_match_excl_singles", "percentage_single_sample_GMYC_species", "oversplitting_ratio", "oversplitting_excl_singles")
                record$filename = files
                
                # create a list to which each gmyc tree is stored
                tree_container = vector("list", length(files))
                # create a list with gmyc.spec output for each file
                gmyc_spec_container = vector("list", length = length(files))
                # create a list of which predefined species are being oversplit by the gmyc algorithm
                oversplitting_species = vector("list")
                
                merged_gmyc_specs = vector("list")
                
                ################################################################################################
                # Loop through the files in the chosen directory to create an ultrametric tree for each,
                # then run the GMYC analysis,
                # and store the entities and clusters
                ################################################################################################
                
                withProgress(message = 'Running GMYC', value = 0, {
                  
                  
                  for(i in seq(along=files)) {
                    
                    
                    if(input$tree_tool == "BEAST") treex.ultra2 = ape::read.nexus(files[i])
                    
                    else if(input$tree_tool == "Non-ultrametric: PATHD8"){
                      treex = ape::read.tree(files[i])
                      treex_tiplabels = treex$tip.label
                      # extract the first two tip labels (could be any two labels)
                      label1 = treex_tiplabels[1]
                      label2 = treex_tiplabels[2]
                      pathd8_params = data.frame(matrix(nrow = 1, ncol = 4))
                      colnames(pathd8_params) = c("tax1", "tax2", "age_type", "age")
                      pathd8_params[1,] = c(label1, label2, "root", 1)
                      
                      tryCatch({
                      pathd8_result = ips::pathd8(phy = treex, exec = pathd8_filepath, seql = input$seqlength, calibration = pathd8_params)
                      treex.ultra = pathd8_result$mpl_tree
                      treex.ultra2 = ape::multi2di(treex.ultra, random = T) # makes the tree fully dichotomous
                      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                      
                    }
                    
                    if (input$set_seed == TRUE) set.seed(1234)
                    
                    # Run the GMYC analysis
                    # tryCatch skips through any possible errors with the gmyc function (e.g. nuclear genes that are identical)
                    
                    tryCatch({
                       treex.gmyc = splits::gmyc(treex.ultra2, quiet = F, method = input$gmyc_threshold)
                    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                    
                    # treex.gmyc = gmyc(treex.ultra2, quiet = F, method = "multiple")
                    
                    # store the gmyc in the tree_container list
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
                    # Do this bit only if the user uploaded a csv file with grouping data and clicked "CONFIRM"
                    ########################################################################################################
                    
                    if (!is.null( predefined_groups_uploaded() )){ 
                      
                    groups = predefined_groups_uploaded()
                    
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
                    
                    gmyc.spec$GMYC_spec = as.factor(gmyc.spec$GMYC_spec) # make the GMYC species a factor
                    gmyc.spec$ids = as.factor(gmyc.spec$ids) # make the user's predefined groups a factor
                    
                    matches.df = data.frame(matrix(nrow = length(levels(gmyc.spec$GMYC_spec)), ncol = 2)) 
                    colnames(matches.df) = c("GMYC_spec", "match(y/n)") 
                    matches.df$GMYC_spec = 1:nrow(matches.df) 
                    
                    # keep counts of the number of "yes", "no", and "single-sample" species:
                    
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
                    
                    # this total excludes single samples
                    total_excl_singles = y_count + n_count 
                    
                    # this total includes single samples
                    total_incl_singles = y_count + n_count + single_sample_count
                    #or just total_incl_singles = nrow(matches.df)
                    
                    # this is only taking the matches for GMYC species groups that had more than one sample in it
                    success_excl_singles = round(y_count/total_excl_singles * 100, 2) 
                    
                    # this treats all single sample GMYC species as equivalent to a "yes" match, and will therefore be a larger value than the estimate
                    # that excludes them
                    success_incl_singles = round((y_count + single_sample_count)/total_incl_singles * 100, 2)
                    
                    prop_single_samples = round(single_sample_count/nrow(matches.df) * 100, 2)
                    
                    # populate the "record" dataframe
                    
                    record[i,2] = success_incl_singles
                    record[i,3] = success_excl_singles
                    record[i,4] = prop_single_samples
                    
                    num_predefined_groups = length(levels(gmyc.spec$ids))
                    num_gmyc_groups = length(levels(gmyc.spec$GMYC_spec))
                    num_gmyc_groups_excluding_single_spp = num_gmyc_groups - single_sample_count
                    
                    # oversplitting ratio (ie. GMYC species to predefined groups)
                    record[i,5] = round( num_gmyc_groups/num_predefined_groups, 2 )
                    # oversplitting ratio excluding single species
                    record[i,6] = round( num_gmyc_groups_excluding_single_spp/num_predefined_groups, 2 )
                    
                    # find which species are being over-split:
                    
                    yes_indices = which(matches.df$`match(y/n)` == "y")
                    
                    if(length(yes_indices) > 0){
                      
                      predef_unique = c()
                      
                      for (v in 1:length(yes_indices)){
                        predef_unique[v] = unique(spec_split[[yes_indices[v]]]$ids)
                      }
                      
                      predef_unique_table = table(predef_unique)
                      predef_unique_table = predef_unique_table
                      predef_unique_table = as.data.frame(predef_unique_table)
                      
                      oversplitting_species[[i]] = predef_unique_table
                    }
                    
                    
                    ##################################
                    # added 18/10/2021:
                    # find which species are undersplit/merged:
                    
                    no_indices = which(matches.df$`match(y/n)` == "n")
                    #print(no_indices)
                    merged_gmyc_specs[[i]] = no_indices
                    
                    #if(length(no_indices) > 0){
                      
                     # output$merged_gmyc_spec = renderText(paste("Merged GMYC species: ", no_indices))
                      
                    # print(no_indices)
                    #  
                    # for(w in length(no_indices)){
                    # 
                    # print(spec_split[no_indices[w]])
                    # 
                    # }
                    
                    
                   # }
                    
                    ###################################
                    
                    
                       }# end of if statement
                    
                    
                    
                    ################################################################################################
                    
                    ################################################################################################
                    
                    # Increment the progress bar, and update the detail text.
                    incProgress(1/length(files), detail = paste("Processing tree file ", i, " of ", length(files), "[ ", round(i/length(files)*100, 0), "% ]"))
                    # Pause for 0.1 seconds to simulate a long computation.
                    Sys.sleep(0.1)
                    
                  } # end of for loop
                  
                  
                }) # end of progress bar bracket
                
                
                shinyalert::shinyalert("Complete", "Please click on the tabs at the top of the window to view the results.", type = "success")
                
                
        } # end of if re. !dir.exists()
          
          if (!is.null( predefined_groups_uploaded() )) output$sp.numbers = renderText(c("You have ", length(unique(groups[[groups_col]])), " predefined species."))
          
          gg_clust_ent = reshape2::melt(clust_ent) # get into the format the ggplot can work with
          colnames(gg_clust_ent) = c("filename", "gmyc_cat", "count")
          gg_clust_ent$gmyc_cat %>% as.factor()
          
          if (!is.null( predefined_groups_uploaded() )){ 
            
            oversplitting_species = oversplitting_species[!is.na(oversplitting_species)] # remove the NA lists
            oversplitting_species_bound = dplyr::bind_rows(oversplitting_species)
            oversplitting_species_bound = oversplitting_species_bound[oversplitting_species_bound$Freq > 1,] # remove rows with species that had a frequency of only one (ie those that were not oversplit)
            
            exact_matches = dplyr::bind_rows(oversplitting_species)
            exact_matches = exact_matches[exact_matches$Freq == 1,]
            
            mean_oversplits_per_grp = 
              oversplitting_species_bound %>% 
              dplyr::group_by(predef_unique) %>%
              dplyr::summarise(across(everything(), c(mean = mean, sd = sd, min = min, max = max)))
            
            colnames(mean_oversplits_per_grp) = c("predefined_group", "mean", "sd", "min", "max")
            
            mean_exact_matches_per_grp =
              exact_matches %>%
              dplyr::group_by(predef_unique) %>%
              dplyr::summarise(across(everything(), c(sum=sum)))

            colnames(mean_exact_matches_per_grp) = c("predefined_group", "sum")
            mean_exact_matches_per_grp$mean = round(mean_exact_matches_per_grp$sum/length(files)*100, 2)
            
            num_user_defined_groups = length(unique(groups[[groups_col]]))
            num_exact_matches = nrow(mean_exact_matches_per_grp)
            percentage_exact_matches = round(num_exact_matches/num_user_defined_groups*100, 2)
            
            output$percent_exact_matches = renderText(c("<B>", "There are ", percentage_exact_matches, "% overall exact GMYC matches."))
            
          }
          
          ################################################################################################
          # View the GMYC species table with the predefined grouping information appended as the last column
          ################################################################################################
          
          observeEvent(input$view_gmyc_spec, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
            gmyc_spec_to_show = input$select_tree_speclist # get the tree file name selected from the drop down menu
            file_i = which(files == gmyc_spec_to_show) # get the index of that file to match up with the trees stored in the spec_list_container list
            
            output$matches = renderTable(gmyc_spec_container[[file_i]], rownames = FALSE, colnames = TRUE, digits = 0)
            output$merges = renderText(c("<B>", "These GMYC species (GMYC_spec) were merged: ", merged_gmyc_specs[[file_i]]))
            }
            
            ################################################################################################
            # Download GMYC species table
            ################################################################################################
            
            output$download_gmyc_spec = downloadHandler(
              
              filename = function (){paste('gmyc_spec_data', 'csv', sep = '.')},
              content = function (file){write.csv(gmyc_spec_container[[file_i]], file, row.names = FALSE)}
            )
            
          })
          
          ################################################################################################
          # View which predefined groups were oversplit by the GMYC
          ################################################################################################
          
          observeEvent(input$GMYC_oversplit_table_view, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
            if(nrow(mean_oversplits_per_grp) > 0){
              
              output$GMYC_oversplit_table = renderTable(mean_oversplits_per_grp, rownames = FALSE, colnames = TRUE, digits = 2)
            }
            
            else shinyalert::shinyalert("No oversplits", "None of your predefined groups were oversplit by the GMYC algorithm", type = "info")
            }
          })
          
          # download the summary table
          output$GMYC_oversplit_table_download = downloadHandler(
            
            filename = function (){paste('mean_oversplits_per_group', 'csv', sep = '.')},
            content = function (file){write.csv(mean_oversplits_per_grp, file, row.names = FALSE)}
          )
          
          # View the full table data
          observeEvent(input$GMCY_oversplit_full_table,{
            output$GMYC_oversplit_table = renderTable(oversplitting_species_bound, rownames = FALSE, colnames = TRUE, digits = 2)
          })
          
          # download the full data table
          output$GMYC_oversplit_full_table_download = downloadHandler(
            
            filename = function (){paste('mean_oversplits_per_group_full', 'csv', sep = '.')},
            content = function (file){write.csv(oversplitting_species_bound, file, row.names = FALSE)}
          )
          
          ################################################################################################
          # Boxplot of which predefined groups were oversplit by the GMYC
          ################################################################################################
          
          observeEvent(input$GMYC_oversplit_boxplot, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
            if(nrow(oversplitting_species_bound) > 0){
              
              output$GMYC_oversplit_plot = renderPlot({
                ggplot(data = oversplitting_species_bound, aes(x = predef_unique, y = Freq)) + 
                  geom_boxplot() +
                  ggthemes[[input$GMYC_oversplit_ggtheme]] +
                  xlab("Predefined group") +
                  ylab("Oversplitting ratio")
              })
            }
            
            else shinyalert::shinyalert("No oversplits", "None of your predefined groups were oversplit by the GMYC algorithm", type = "info")
            } 
          })
          
          # download the boxplot
          output$GMYC_oversplit_boxplot_download = downloadHandler(
            
            filename = function (){paste(input$file_name_oversplit, input$plot_format_oversplit, sep = '.')},
            
            content = function (file){
              
              width = as.numeric(input$w_plot_oversplit) 
              height = as.numeric(input$h_plot_oversplit) 
              dpi = as.numeric(input$res_plot_oversplit)
              units = input$unit_plot_oversplit
              
              ggsave(file, width = width, height = height, dpi = dpi, units = units,
                     ggplot(data = oversplitting_species_bound, aes(x = predef_unique, y = Freq)) + 
                       geom_boxplot() +
                       ggthemes[[input$GMYC_oversplit_ggtheme]] +
                       xlab("Predefined group") +
                       ylab("Oversplitting ratio")
              )}
          )
          
          ################################################################################################
          # Barplot of which predefined groups were oversplit by the GMYC
          ################################################################################################
          observeEvent(input$GMYC_oversplit_barplot, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
            if(nrow(oversplitting_species_bound) > 0){
              
              output$GMYC_oversplit_plot = renderPlot({
                
                ggplot(data = mean_oversplits_per_grp, aes(x=predefined_group, y=mean)) + 
                  geom_bar(stat = "identity", colour = input$GMYC_barchart_outline, fill = input$GMYC_barchart_fill) + 
                  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .1) + 
                  xlab("Predefined group") +
                  ylab("Oversplitting ratio") +
                  ggthemes[[input$GMYC_oversplit_ggtheme]] 
              })
              
              
            }
            
            else shinyalert::shinyalert("No oversplits", "None of your predefined groups were oversplit by the GMYC algorithm", type = "info")
            } 
          }) 
          
          # download the barplot
          output$GMYC_oversplit_barplot_download = downloadHandler(
            
            filename = function (){paste(input$file_name_oversplitbar, input$plot_format_oversplitbar, sep = '.')},
            
            content = function (file){
              
              width = as.numeric(input$w_plot_oversplitbar) 
              height = as.numeric(input$h_plot_oversplitbar) 
              dpi = as.numeric(input$res_plot_oversplitbar)
              units = input$res_plot_oversplitbar
              
              ggsave(file, width = width, height = height, dpi = dpi, units = units,
                     ggplot(data = mean_oversplits_per_grp, aes(x=predefined_group, y=mean)) + 
                       geom_bar(stat = "identity", colour = input$GMYC_barchart_outline, fill = input$GMYC_barchart_fill) + 
                       geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .1) + 
                       xlab("Predefined group") +
                       ylab("Oversplitting ratio") +
                       ggthemes[[input$GMYC_oversplit_ggtheme]] 
              )}
          )
          
          ################################################################################################
          # View which predefined groups were exact GMYC matches
          ################################################################################################
          
          observeEvent(input$GMYC_exact_match_table_view, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
              if(nrow(mean_exact_matches_per_grp) > 0){
                
                output$GMYC_exact_match_table = renderTable(mean_exact_matches_per_grp, rownames = FALSE, colnames = TRUE, digits = 2)
              }
              
              else shinyalert::shinyalert("No exact matches", "None of your predefined groups were exact GMYC matches", type = "info")
            }
          })
          
          # download the summary table
          output$GMYC_exact_match_table_download = downloadHandler(
            
            filename = function (){paste('mean_exact_matches_per_group', 'csv', sep = '.')},
            content = function (file){write.csv(mean_exact_matches_per_grp, file, row.names = FALSE)}
          )
          
          # View the full table data
          observeEvent(input$GMCY_exact_match_full_table,{
            output$GMYC_exact_match_table = renderTable(exact_matches, rownames = FALSE, colnames = TRUE, digits = 2)
          })
          
          # download the full data table
          output$GMYC_exact_matches_full_table_download = downloadHandler(
            
            filename = function (){paste('mean_exact_matches_per_group_full', 'csv', sep = '.')},
            content = function (file){write.csv(exact_matches, file, row.names = FALSE)}
          )
          
          ################################################################################################
          # Barplot of which predefined groups were exact GMYC matches
          ################################################################################################
          observeEvent(input$GMYC_exact_match_barplot, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
              if(nrow(exact_matches) > 0){
                
                output$GMYC_exact_match_plot = renderPlot({
                  
                  ggplot(data = mean_exact_matches_per_grp, aes(x=predefined_group, y=mean)) + 
                    geom_bar(stat = "identity", colour = input$exact_match_barchart_outline, fill = input$exact_match_barchart_fill) + 
                    ylim(0,100) +
                    xlab("Predefined group") +
                    ylab("Mean exact match (%)") +
                    ggthemes[[input$GMYC_exact_match_ggtheme]] +
                    theme(axis.text.x = element_text(angle = req( input$x_axis_angle_exact_matches)) )
                })
                
                
              }
              
              else shinyalert::shinyalert("No exact matches", "None of your predefined groups were exact GMYC matches", type = "info")
            } 
          }) 
          
          # download the barplot
          output$GMYC_exact_match_barplot_download = downloadHandler(
            
            filename = function (){paste(input$file_name_exact_match_bar, input$plot_format_exact_match, sep = '.')},
            
            content = function (file){
              
              width = as.numeric(input$w_plot_exact_match) 
              height = as.numeric(input$h_plot_exact_match) 
              dpi = as.numeric(input$res_plot_exact_match)
              units = input$res_plot_exact_match
              
              ggsave(file, width = width, height = height, dpi = dpi, units = units,
                     ggplot(data = mean_exact_matches_per_grp, aes(x=predefined_group, y=mean)) + 
                       geom_bar(stat = "identity", colour = input$exact_match_barchart_outline, fill = input$exact_match_barchart_fill) + 
                       ylim(0,100) +
                       xlab("Predefined group") +
                       ylab("Mean exact match (%)") +
                       ggthemes[[input$GMYC_exact_match_ggtheme]] +
                       theme(axis.text.x = element_text(angle = input$x_axis_angle_exact_matches)) 
              )}
          )
          
          ################################################################################################
          # View the match data for each file
          ################################################################################################
          
          observeEvent(input$view_match_data, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
            output$matches = renderTable(record, rownames = FALSE, colnames = TRUE, digits = 2)
            }
           
          })
          
          ################################################################################################
          # Download match data for each file
          ################################################################################################
          
          output$download_match_data = downloadHandler(
            
            filename = function (){paste('match_data_full', 'csv', sep = '.')},
            content = function (file){write.csv(record, file, row.names = FALSE)}
          )
          
          ################################################################################################
          # Plot percentage match data for all files as a ggplot
          ################################################################################################
          
          observeEvent(input$plot_matches, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
            chosen_match_type = ""
            if( input$plot_matches_type == 1) {
              chosen_match_type = as.name("percentage_match")
              ggtit = "GMYC species match to predefined groups vs ML tree/iteration file (including single-sample species)"
            }
            else if (input$plot_matches_type == 2) {
              chosen_match_type = as.name("percentage_match_excl_singles")
              ggtit = ggtit = "GMYC species match to predefined groups vs ML tree/iteration file (excluding single-sample species)"
            }
            
            output$match_plot = renderPlot({
              ggplot(record, aes(x=1:nrow(record), y= .data[[chosen_match_type]])) + 
                geom_line(color = input$plot_matches_line_colour) + 
                geom_point(color = input$plot_matches_point_colours)  + 
                ggthemes[[input$ggtheme_matches]] +
                xlab("ML tree file/iteration number") + 
                ylab("Percentage match") + 
                ggtitle(ggtit)
            })
            }
          })
          
          ################################################################################################
          # Download Plot of percentage match data for all files as a ggplot
          ################################################################################################
          
          output$download_match_plot = downloadHandler(
            filename = function (){paste(input$file_name_percmatches, input$plot_format_percmatches, sep = '.')},
            
            content = function (file){
              
              width = as.numeric(input$w_plot_percmatches) 
              height = as.numeric(input$h_plot_percmatches) 
              dpi = as.numeric(input$res_plot_percmatches)
              units = input$res_plot_percmatches
              
              chosen_match_type = ""
              if( input$plot_matches_type == 1) {
                chosen_match_type = as.name("percentage_match")
                ggtit = "GMYC species match to predefined groups vs ML tree/iteration file (including single-sample species)"
              }
              else if (input$plot_matches_type == 2) {
                chosen_match_type = as.name("percentage_match_excl_singles")
                ggtit = ggtit = "GMYC species match to predefined groups vs ML tree/iteration file (excluding single-sample species)"
              }
              
              ggsave(file, width = width, height = height, dpi = dpi, units = units,
                     ggplot(record, aes(x=1:nrow(record), y=.data[[chosen_match_type]])) + 
                       geom_line(color = input$plot_matches_line_colour) + 
                       geom_point(color = input$plot_matches_point_colours)  + 
                       ggthemes[[input$ggtheme_matches]] +
                       xlab("ML tree file/iteration number") + 
                       ylab("Percentage match") + 
                       ggtitle(ggtit))
            }
            
          )
          
          
          ################################################################################################
          # view the match data summary for all files (mean, stdev etc.)
          ################################################################################################
          
          observeEvent(input$view_summary_match_data, {
            
            if (is.null( predefined_groups_uploaded() )) shinyalert::shinyalert("No grouping information uploaded", "You did not upload a .csv file with predefined grouping information for your sequences", type = "warning")
            else{
            ################################################################################################
            # Get summary statistics for the success record and number of single samples assigned to a GMYC species
            ################################################################################################
            
            # stats for percentage matches including single samples
            avg_match = mean(record$percentage_match)
            stdev_match = sd(record$percentage_match)
            min_match = min(record$percentage_match)
            max_match = max(record$percentage_match)
            match_summary = round( c(avg_match, stdev_match, min_match, max_match) ,2 )
            
            # stats for percentage matches excluding single samples
            avg_match_excl_singles = mean(record$percentage_match_excl_singles)
            stdev_match_excl_singles = sd(record$percentage_match_excl_singles)
            min_match_excl_singles = min(record$percentage_match_excl_singles)
            max_match_excl_singles = max(record$percentage_match_excl_singles)
            match_excl_singles_summary = round( c(avg_match_excl_singles, stdev_match_excl_singles, min_match_excl_singles, max_match_excl_singles) ,2 )
            
            # stats for percentage of single-ssample gmyc species
            avg_single = mean(record$percentage_single_sample_GMYC_species)
            stdev_single = sd(record$percentage_single_sample_GMYC_species)
            min_single = min(record$percentage_single_sample_GMYC_species)
            max_single = max(record$percentage_single_sample_GMYC_species)
            single_summary = round( c(avg_single, stdev_single, min_single, max_single) ,2 )
            
            # stats for percentage of GMYC:predefined_species
            avg_oversplit_ratio = mean(record$oversplitting_ratio)
            stdev_oversplit_ratio = sd(record$oversplitting_ratio)
            min_oversplit_ratio = min(record$oversplitting_ratio)
            max_oversplit_ratio = max(record$oversplitting_ratio)
            oversplit_ratio_summary = round( c(avg_oversplit_ratio, stdev_oversplit_ratio, min_oversplit_ratio, max_oversplit_ratio) ,2 )
            
            # stats for percentage of GMYC:predefined_species excluding all single sample representatives
            avg_oversplit_ratio_excl_single = mean(record$oversplitting_excl_singles)
            stdev_oversplit_ratio_excl_single = sd(record$oversplitting_excl_singles)
            min_oversplit_ratio_excl_single = min(record$oversplitting_excl_singles)
            max_oversplit_ratio_excl_single = max(record$oversplitting_excl_singles)
            oversplit_ratio_summary_excl_single = round( c(avg_oversplit_ratio_excl_single, stdev_oversplit_ratio_excl_single, min_oversplit_ratio_excl_single, max_oversplit_ratio_excl_single) ,2 )
            
            
            # create a data frame housing summary statistics, called 'match_stats.df'
            match_stats.df = data.frame(matrix(nrow=4, ncol = 6)) 
            colnames(match_stats.df) = c("statistic", "percentage_matches", "percentage_matches_excl_single_samples", "percentage_single_sample_GMYC_species", "oversplitting_ratio", "oversplitting_ratio_excl_single_sample_spp")
            categories = c("Average", "Standard deviation", "Minimum", "Maximum")
            match_stats.df$statistic = categories
            match_stats.df$percentage_matches = match_summary
            match_stats.df$percentage_matches_excl_single_samples = match_excl_singles_summary
            match_stats.df$percentage_single_sample_GMYC_species = single_summary
            match_stats.df$oversplitting_ratio = oversplit_ratio_summary
            match_stats.df$oversplitting_ratio_excl_single_sample_spp = oversplit_ratio_summary_excl_single
            
            output$matches =  renderTable(match_stats.df, rownames = FALSE, colnames = TRUE, digits = 2)
            }
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
          })
          
          
          ################################################################################################
          # Plot entities vs clusters
          ################################################################################################
          
          observeEvent(input$plot_clusts, {
            
            output$clust_ent_plot = renderPlot({
              
              ggplot(clust_ent, aes(x=clusters, y=entities)) + geom_point(colour = input$clust_vs_ent_plot_point_colours, shape = as.numeric(input$clust_vs_ent_plot_point_shape), size = input$clust_vs_ent_plot_point_size) +
                xlab("Clusters") + 
                ylab("Entities") +
                ggthemes[[input$ggtheme_plots]] +
                ggtitle("Number of Clusters vs Entities for GMYC")
              
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
              
              support_gmyc = splits::gmyc.support(tree_container[[file_index]]) # get the support values as estimated by the GMYC analysis
              is.na(support_gmyc[support_gmyc == 0]) = TRUE
              
              support_original_tree = as.numeric(treex.ultra2$node.label) # get the original support values straight from the ML tree read in
              
              plot.result.gmyc(tree_container[[file_index]], cex = tip_label_size, edge.width = branch_width)
              
              if (input$support_value_type == "GMYC estimate") nodelabels(round(support_gmyc, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
              
              else nodelabels(round(support_original_tree, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
              
              
              ################################################################################################
              # Download GMYC tree
              ################################################################################################
              
              
              output$download_gmyc_tree <- downloadHandler(
                
                filename = function (){paste(input$file_name_tree, input$plot_format_tree, sep = '.')},
                
                content = function (file){
                  width = as.numeric(input$w_plot_tree) * ifelse(input$unit_plot_tree == "cm", 1, 2.54)
                  height = as.numeric(input$h_plot_tree) * ifelse(input$unit_plot_tree == "cm", 1, 2.54)
                  
                  if(input$plot_format_tree == "png"){
                    png(filename = file, 
                        width = width, height = height,
                        res = as.numeric(input$res_plot), units = "cm")
                    
                    plot.result.gmyc(tree_container[[file_index]], cex = tip_label_size, edge.width = branch_width)
                    if (input$support_value_type == "GMYC estimate") nodelabels(round(support_gmyc, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
                    else nodelabels(round(support_original_tree, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
                    
                  }
                  
                  if(input$plot_format_tree == "pdf"){
                    pdf(file = file, 
                        width = width , 
                        height = height ) 
                    
                    plot.result.gmyc(tree_container[[file_index]], cex = tip_label_size, edge.width = branch_width)
                    if (input$support_value_type == "GMYC estimate") nodelabels(round(support_gmyc, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
                    else nodelabels(round(support_original_tree, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
                    
                  }
                  
                  if(input$plot_format_tree == "svg"){
                    svg(filename = file, 
                        width = width, 
                        height = height)
                    
                    plot.result.gmyc(tree_container[[file_index]], cex = tip_label_size, edge.width = branch_width)
                    if (input$support_value_type == "GMYC estimate") nodelabels(round(support_gmyc, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
                    else nodelabels(round(support_original_tree, 2), cex = support_value_size, bg = support_value_col, frame = support_value_frame)
                    
                  }
                  
                  dev.off()
                }
              )
              
            })
            
          })
          
          ################################################################################################
          # Download clusters vs entities plot
          ################################################################################################
          
          output$download_clust_plot = downloadHandler(
            
            filename = function (){paste(input$file_name_clustvsent, input$plot_format_clustvsent, sep = '.')},
            
            content = function (file){
              width.clustplot = as.numeric(input$w_plot_clustvsent) 
              height.clustplot = as.numeric(input$h_plot_clustvsent) 
              dpi.clustplot = as.numeric(input$res_plot_clustvsent)
              units.clustplot = input$unit_plot_clustvsent
              
              ggsave(file, width = width.clustplot, height = height.clustplot, dpi = dpi.clustplot, units = units.clustplot,
                     ggplot(clust_ent, aes(x=clusters, y=entities)) + geom_point(colour = input$clust_vs_ent_plot_point_colours, shape = as.numeric(input$clust_vs_ent_plot_point_shape), size = input$clust_vs_ent_plot_point_size) +
                       xlab("Clusters") + 
                       ylab("Entities") +
                       ggthemes[[input$ggtheme_plots]] +
                       ggtitle("Number of Clusters vs Entities for GMYC") )
              
            }
          )
          
          ################################################################################################
          # Plot entities and clusters boxplot
          ################################################################################################
          
          observeEvent(input$plot_boxplot, {
            
            output$clust_ent_plot = renderPlot({
              
              ggplot(gg_clust_ent, aes(x=gmyc_cat, y=count)) + geom_boxplot() +
                ggthemes[[input$ggtheme_plots]] +
                xlab("GMYC Category") + 
                ylab("Count") 
            })
          })
          
          ################################################################################################
          # Download entities and clusters boxplot
          ################################################################################################
          
          output$download_boxplot <- downloadHandler(
            
            filename = function (){paste(input$file_name_clust_ent_box, input$plot_format_clust_ent_box, sep = '.')},
            
            content = function (file){
              
              width.clustentboxplot = as.numeric(input$w_plot_clust_ent_box) 
              height.clustentboxplot = as.numeric(input$h_plot_clust_ent_box) 
              dpi.clustentboxplot = as.numeric(input$res_plot_clust_ent_box)
              units.clustentboxplot = input$unit_plot_clust_ent_box
              
              ggsave(file, width = width.clustentboxplot, height = height.clustentboxplot, dpi = dpi.clustentboxplot, units = units.clustentboxplot,
                     ggplot( gg_clust_ent, aes(x=gmyc_cat, y=count)) + geom_boxplot() +
                       ggthemes[[input$ggtheme_plots]] +
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
                geom_point(color = input$plot_clusts_vs_iterations_point_colours, size = input$plot_clusts_vs_iterations_point_size, shape = as.numeric(input$plot_clusts_vs_iterations_point_shape))  + 
                ggthemes[[input$ggtheme_plots]] +
                xlab("ML tree file/iteration number") + 
                ylab("Number of GMYC clusters") + 
                ggtitle("GMYC number of clusters vs ML tree/iteration file")
            })
          })
          
          ################################################################################################
          # Download clusters vs iteration/ML file number plot
          ################################################################################################
          
          output$download_clusts_vs_iterations <- downloadHandler(
            
            filename = function (){paste(input$file_name_clustvsiter, input$plot_format_clustvsiter, sep = '.')},
            
            content = function (file){
              
              width.clustvsiter = as.numeric(input$w_plot_clustvsiter) 
              height.clustvsiter = as.numeric(input$h_plot_clustvsiter) 
              dpi.clustvsiter = as.numeric(input$res_plot_clustvsiter)
              units.clustvsiter = input$unit_plot_clustvsiter
              
              ggsave(file, width = width.clustvsiter, height = height.clustvsiter, dpi = dpi.clustvsiter, units = units.clustvsiter,
                     ggplot( clust_ent, aes(x=1:nrow(clust_ent), y=clusters)) + 
                       geom_line(color = input$plot_clusts_vs_iterations_line_colour) + 
                       geom_point(color = input$plot_clusts_vs_iterations_point_colours, size = input$plot_clusts_vs_iterations_point_size, shape = as.numeric(input$plot_clusts_vs_iterations_point_shape))  + 
                       ggthemes[[input$ggtheme_plots]] +
                       xlab("ML tree file/iteration number") + 
                       ylab("Number of GMYC clusters") + 
                       ggtitle("GMYC number of clusters vs ML tree/iteration file") )
            }
          )
          
          ################################################################################################
          # Plot entities vs iteration/ML file number
          ################################################################################################
          
          observeEvent(input$plot_ents_vs_iterations, {
            
            output$clust_ent_plot = renderPlot({
              ggplot(clust_ent, aes(x=1:nrow(clust_ent), y=entities)) + 
                geom_line(color = input$plot_ents_vs_iterations_line_colour) + 
                geom_point(color = input$plot_ents_vs_iterations_point_colours, size = input$plot_ents_vs_iterations_point_size, shape = as.numeric(input$plot_ents_vs_iterations_point_shape))  + 
                ggthemes[[input$ggtheme_plots]] +
                xlab("ML tree file/iteration number") + 
                ylab("Number of GMYC entities") + 
                ggtitle("GMYC number of entities vs ML tree/iteration file")
            })
          })
          
          ################################################################################################
          # Download entities vs iteration/ML file number plot
          ################################################################################################
          
          output$download_ents_vs_iterations <- downloadHandler(
            
            filename = function (){paste(input$file_name_entvsiter, input$plot_format_entvsiter, sep = '.')},
            
            content = function (file){
              
              width.entvsiter = as.numeric(input$w_plot_entvsiter) 
              height.entvsiter = as.numeric(input$h_plot_entvsiter) 
              dpi.entvsiter = as.numeric(input$res_plot_entvsiter)
              units.entvsiter = input$unit_plot_entvsiter
              
              ggsave(file, width = width.entvsiter, height = height.entvsiter, dpi = dpi.entvsiter, units = units.entvsiter,
                     ggplot(clust_ent, aes(x=1:nrow(clust_ent), y=entities)) + 
                       geom_line(color = input$plot_ents_vs_iterations_line_colour) + 
                       geom_point(color = input$plot_ents_vs_iterations_point_colours, size = input$plot_ents_vs_iterations_point_size, shape = as.numeric(input$plot_ents_vs_iterations_point_shape))  + 
                       ggthemes[[input$ggtheme_plots]] +
                       xlab("ML tree file/iteration number") + 
                       ylab("Number of GMYC entities") + 
                       ggtitle("GMYC number of entities vs ML tree/iteration file") )
            }
          )
          
          
          
        }# end of else 1
        
      }# end of else 2
      
    })
    
  # }) # end of observe of input of predefined group file
  
  ################################################################################################
  # Read in csv files with summary data
  ################################################################################################
  
  multiple_summary_plots = list()
  
    # CLUSTER DATA
    
    observe({
    infile_summary_clusters = input$summary_clusters
    
    if (is.null(infile_summary_clusters)) {
      return(NULL)
    }
    
    
    observeEvent(input$plot_summary_clusts, {
      
      summary_clusters_csv = read.csv(infile_summary_clusters$datapath, check.names = F)
      summary_clusters_csv_melt = reshape2::melt(summary_clusters_csv)
      
      if(input$no_morphospecies > 0) linetype_morpho = "dashed"
      else linetype_morpho = "blank"
      colnames(summary_clusters_csv_melt) = c("file_name", "data_perc", "clusters")
      
      # CREATE PLOT
      clust_accum <<- ggplot(data = summary_clusters_csv_melt, aes(x = as.numeric( data_perc ), y = clusters)) + 
        geom_point() +
        geom_smooth(span = 3, col = input$clust_ent_summary_line_col, fill = input$clust_ent_summary_ci_col, alpha = input$clust_ent_summary_ci_alpha) +
        scale_y_continuous(breaks = seq(0, max(summary_clusters_csv_melt$clusters), by = input$clust_ent_summary_y_interval)) +
        #expand_limits(y = 0) +
        #geom_boxplot(fill = "#1EAD4B") +
        #stat_summary(fun = mean, geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
        xlab("Data Percentage") +
        ylab("Number of clusters") +
        #ggtitle(input$clust_ent_summary_plot_title) +
        ggthemes[[input$ggtheme_summary_plots]] +
        theme(legend.position = "bottom") +
        theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
        theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
        theme(axis.text.x = element_text(size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(legend.text=element_text(size=12)) +
        theme(plot.title = element_text(size=16, face = "bold")) +
        scale_x_continuous(breaks = as.numeric(summary_clusters_csv_melt$data_perc), labels = summary_clusters_csv_melt$data_perc) +
        geom_hline(yintercept= input$no_morphospecies, linetype= linetype_morpho)
      
      output$summary_plot = renderPlot(clust_accum)
      multiple_summary_plots$clust_accum <<- clust_accum 
      
    })
    
    # download the plot
    output$download_clust_summary_plot <- downloadHandler(
      
      filename = function (){paste("clust_summary_plot", input$plot_format_summary_clusts_ents, sep = '.')},
      
      content = function (file){
        
        width.clust_ent_summary = as.numeric(input$w_plot_summary_clusts_ents) 
        height.clust_ent_summary = as.numeric(input$h_plot_summary_clusts_ents) 
        dpi.clust_ent_summary = as.numeric(input$res_plot_summary_clusts_ents)
        units.clust_ent_summary = input$unit_plot_summary_clusts_ents
        
        ggsave(file, width = width.clust_ent_summary, height = height.clust_ent_summary, dpi = dpi.clust_ent_summary, units = units.clust_ent_summary,
               clust_accum)
      }
    )
    
  }) # end of observe
    
    ################################################################################################
    
  observe({
    
    # ENTITIES DATA
    
    infile_summary_entities = input$summary_entities
    
    if (is.null(infile_summary_entities)) {
      return(NULL)
    }
    
    observeEvent(input$plot_summary_ents, {
      
      summary_entities_csv = read.csv(infile_summary_entities$datapath, check.names = F)
      summary_entities_csv_melt = reshape2::melt(summary_entities_csv)
      
      if(input$no_morphospecies > 0) linetype_morpho = "dashed"
      else linetype_morpho = "blank"
      
      colnames(summary_entities_csv_melt) = c("file_name", "data_perc", "entities")
      
      ent_accum <<- ggplot(data = summary_entities_csv_melt, aes(x = as.numeric( data_perc ), y = entities)) + 
        geom_point() +
        geom_smooth(span = 3, col = input$clust_ent_summary_line_col, fill = input$clust_ent_summary_ci_col, alpha = input$clust_ent_summary_ci_alpha) +
        scale_y_continuous(breaks = seq(0, max(summary_entities_csv_melt$entities), by = input$clust_ent_summary_y_interval)) +
        #expand_limits(y = 0) +
        #geom_boxplot(fill = "#1EAD4B") +
        #stat_summary(fun = mean, geom="point", shape=17, size=3, color="black", position=position_dodge(0.77)) +
        xlab("Data Percentage") +
        ylab("Number of entities") +
        ggthemes[[input$ggtheme_summary_plots]] +
        theme(legend.position = "bottom") +
        theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
        theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
        theme(axis.text.x = element_text(size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(legend.text=element_text(size=12)) +
        theme(plot.title = element_text(size=16, face = "bold")) +
        scale_x_continuous(breaks = as.numeric(summary_entities_csv_melt$data_perc), labels = summary_entities_csv_melt$data_perc) +
        geom_hline(yintercept= input$no_morphospecies, linetype= linetype_morpho)
      
      output$summary_plot = renderPlot(ent_accum)
      multiple_summary_plots$ent_accum <<- ent_accum
      
    })
    
    # download the plot
    output$download_ent_summary_plot <- downloadHandler(
      
      filename = function (){paste("ent_summary_plot", input$plot_format_summary_clusts_ents, sep = '.')},
      
      content = function (file){
        
        width.clust_ent_summary = as.numeric(input$w_plot_summary_clusts_ents) 
        height.clust_ent_summary = as.numeric(input$h_plot_summary_clusts_ents) 
        dpi.clust_ent_summary = as.numeric(input$res_plot_summary_clusts_ents)
        units.clust_ent_summary = input$unit_plot_summary_clusts_ents
        
        ggsave(file, width = width.clust_ent_summary, height = height.clust_ent_summary, dpi = dpi.clust_ent_summary, units = units.clust_ent_summary,
               ent_accum)
      }
    )
    
  }) # end of observe
    
    ################################################################################################
    
  observe({
    
    # OVERSPLITTING RATIO DATA INC. SINGLETONS
    
    infile_summary_oversplitting_ratio_inc = input$summary_oversplitting_ratio_inc
    
    if (is.null(infile_summary_oversplitting_ratio_inc)) {
      return(NULL)
    }
    
    # PLOT
    observeEvent(input$plot_summary_oversplits, {
      
      summary_oversplitting_inc_csv = read.csv(infile_summary_oversplitting_ratio_inc$datapath, header = T, check.names = F)
      summary_oversplitting_inc_csv_melt = reshape2::melt(summary_oversplitting_inc_csv)
      
      # OVERSPLITTING RATIO DATA EXCL SINGLETONS
      infile_summary_oversplitting_ratio_exc = input$summary_oversplitting_ratio_exc
      
      if (is.null(infile_summary_oversplitting_ratio_exc)) {
        return(NULL)
      }
      
      summary_oversplitting_exc_csv = read.csv(infile_summary_oversplitting_ratio_exc$datapath, header = T, check.names = F)
      summary_oversplitting_exc_csv_melt = reshape2::melt(summary_oversplitting_exc_csv)
      
      colnames(summary_oversplitting_inc_csv_melt) = c("file_name", "data_perc", "oversplitting_ratio")
      colnames(summary_oversplitting_exc_csv_melt) = c("file_name", "data_perc", "oversplitting_ratio_excl")
      
      summary_oversplitting_inc_csv_melt$oversplitting_excl = summary_oversplitting_exc_csv_melt$oversplitting_ratio_excl
      oversplitting_combo = reshape::melt( summary_oversplitting_inc_csv_melt )
      
      oversplitting_boxplot_summary <<- ggplot(data = oversplitting_combo, aes(x = data_perc, y = value)) + 
        geom_boxplot(aes(fill = variable)) +
        expand_limits(y = 0) +
        #scale_y_continuous(breaks = seq(0, max(oversplitting_combo$value), by = input$oversplitting_plot_summary_y_interval)) +
        stat_summary(fun = mean, aes(group = variable), geom="point", shape = as.numeric(input$oversplitting_plot_point_shape), size=input$oversplitting_plot_point_size, color=input$oversplitting_plot_point_col, position=position_dodge(0.77)) +
        scale_fill_manual(values=c(input$oversplitting_inc_col, input$oversplitting_exc_col), name = "", labels = c("+ singletons", "- singletons")) +
        xlab("Data Percentage") +
        ylab("Oversplitting ratio") +
        ggthemes[[input$ggtheme_summary_plots]] +
        theme(legend.position = "bottom") +
        guides(fill=guide_legend(title="")) +
        theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
        theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
        theme(axis.text.x = element_text(size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(legend.text=element_text(size=12)) +
        geom_hline(yintercept= 1, linetype="dashed") +
        theme(plot.title = element_text(size=16, face = "bold"))
      
      output$summary_plot = renderPlot(oversplitting_boxplot_summary)
      multiple_summary_plots$oversplitting <<- oversplitting_boxplot_summary
      
    }
      
    )
    
    # download the plot
    output$download_summary_oversplits_plot <- downloadHandler(
      
      filename = function (){paste("oversplitting_summary_plot", input$plot_format_summary_oversplits, sep = '.')},
      
      content = function (file){
        
        width.oversplits_summary = as.numeric(input$w_plot_summary_oversplits) 
        height.oversplits_summary = as.numeric(input$h_plot_summary_oversplits) 
        dpi.oversplits_summary = as.numeric(input$res_plot_summary_oversplits)
        units.oversplits_summary = input$unit_plot_summary_oversplits
        
        ggsave(file, width = width.oversplits_summary, height = height.oversplits_summary, dpi = dpi.oversplits_summary, units = units.oversplits_summary,
               oversplitting_boxplot_summary)
      }
    )
    
  }) # end of observe
  
    ################################################################################################
    
    
  observe({
    
    infile_summary_percentage_match_inc = input$summary_percentage_match_inc
    
    if (is.null(infile_summary_percentage_match_inc)) {
      return(NULL)
    }
  
    # PERCENTAGE MATCH DATA INCLUDING SINGLETONS
    
    observeEvent(input$plot_summary_percentage_matches, {
      
      summary_percentage_matches_inc_csv = read.csv(infile_summary_percentage_match_inc$datapath, check.names = F)
      summary_percentage_matches_inc_csv_melt = reshape2::melt(summary_percentage_matches_inc_csv)
      
      # PERCENTAGE MATCH DATA EXCLUDING SINGLETONS
      infile_summary_percentage_match_exc = input$summary_percentage_match_exc
      
      if (is.null(infile_summary_percentage_match_exc)) {
        return(NULL)
      }
      
      summary_percentage_matches_exc = read.csv(infile_summary_percentage_match_exc$datapath, check.names = F)
      summary_percentage_matches_exc_melt = reshape2::melt(summary_percentage_matches_exc)
      
      summary_percentage_matches_inc_csv_melt$ex_sing = summary_percentage_matches_exc_melt$value
      colnames(summary_percentage_matches_inc_csv_melt) = c("file_name", "data_perc", "inc_singletons", "ex_singletons")
      percentage_matches_combo = reshape2::melt(summary_percentage_matches_inc_csv_melt)
      
      perc_match_boxplot_summary <<- ggplot(data = percentage_matches_combo, aes(x = data_perc, y = value)) + 
        geom_boxplot(aes(fill = variable)) +
        stat_summary(fun = mean, aes(group = variable), geom="point", shape=as.numeric(input$percentage_match_plot_point_shape), size=input$percentage_match_plot_point_size, color=input$percentage_match_plot_point_col, position=position_dodge(0.77)) +
        ylim(0,100) +
        scale_fill_manual(values=c(input$percentage_match_inc_col, input$percentage_match_exc_col), name = "", labels = c("+ singletons", "- singletons")) +
        xlab("Data Percentage") +
        ylab("Percentage match") +
        ggthemes[[input$ggtheme_summary_plots]] +
        theme(legend.position = "bottom") +
        guides(fill=guide_legend(title="")) +
        theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
        theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
        theme(axis.text.x = element_text(size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(legend.text=element_text(size=12)) +
        theme(plot.title = element_text(size=16, face = "bold"))
      
      output$summary_plot = renderPlot(perc_match_boxplot_summary)
      multiple_summary_plots$perc_match <<- perc_match_boxplot_summary
      
    }
    )
    
    # download the plot
    output$download_summary_percentage_matches_plot <- downloadHandler(
      
      filename = function (){paste("percentage_matches_plot", input$plot_format_summary_percentage_matches, sep = '.')},
      
      content = function (file){
        
        width.percentage_matches_summary = as.numeric(input$w_plot_summary_percentage_matches) 
        height.percentage_matches_summary = as.numeric(input$h_plot_summary_percentage_matches) 
        dpi.percentage_matches_summary = as.numeric(input$res_plot_summary_percentage_matches)
        units.percentage_matches_summary = input$unit_plot_summary_percentage_matches
        
        ggsave(file, width = width.percentage_matches_summary, height = height.percentage_matches_summary, dpi = dpi.percentage_matches_summary, units = units.percentage_matches_summary,
               perc_match_boxplot_summary)
      }
    )
  
  }) # end of observe
  
    ################################################################################################
  
  
   observe({
     
     infile_summary_percentage_singletons = input$summary_percentage_singletons
     
     if (is.null(infile_summary_percentage_singletons)) {
       return(NULL)
     }
  
    # PERCENTAGE SINGLETONS
    
    observeEvent(input$plot_summary_singletons, {
      
      summary_singletons_csv = read.csv(infile_summary_percentage_singletons$datapath, header = T, check.names = F)
      summary_singletons_csv_melt = reshape2::melt(summary_singletons_csv)
      
      colnames(summary_singletons_csv_melt) = c("file_name", "data_perc", "perc_singletons")
      
      singletons_boxplot_summary <<- ggplot(data = summary_singletons_csv_melt, aes(x = data_perc, y = perc_singletons)) +
        geom_boxplot(fill = input$percentage_singletons_col) +
        ylim(0,100) +
        stat_summary(fun = mean, geom="point", shape=as.numeric(input$percentage_singletons_plot_point_shape), size=input$percentage_singletons_plot_point_size, color=input$percentage_singletons_plot_point_col, position=position_dodge(0.77)) +
        xlab("Data Percentage") +
        ylab("Percentage singletons") +
        ggthemes[[input$ggtheme_summary_plots]] +
        theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
        theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
        theme(axis.text.x = element_text(size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(plot.title = element_text(size=16, face = "bold"))
      
      output$summary_plot = renderPlot(singletons_boxplot_summary)
      multiple_summary_plots$singletons <<- singletons_boxplot_summary
      
    })
    
    # download the plot
    output$download_summary_percentage_singletons_plot <- downloadHandler(
      
      filename = function (){paste("percentage_singletons_plot", input$plot_format_summary_percentage_singletons, sep = '.')},
      
      content = function (file){
        
        width.percentage_singletons_summary = as.numeric(input$w_plot_summary_percentage_singletons) 
        height.percentage_singletons_summary = as.numeric(input$h_plot_summary_percentage_singletons) 
        dpi.percentage_singletons_summary = as.numeric(input$res_plot_summary_percentage_singletons)
        units.percentage_singletons_summary = input$unit_plot_summary_percentage_singletons
        
        ggsave(file, width = width.percentage_singletons_summary, height = height.percentage_singletons_summary, dpi = dpi.percentage_singletons_summary, units = units.percentage_singletons_summary,
               singletons_boxplot_summary)
      }
    )
    
  }) # end of observe
  
    ################################################################################################
    
    # GENETIC DIVERGENCE DATA
    
    observe({
      infile_genetic_diversity = input$genetic_divergence_data
      
      if (is.null(infile_genetic_diversity)) {
        return(NULL)
      }
      
      
      observeEvent(input$plot_genetic_divergence, {
        
        genetic_divergences_csv = read.csv(infile_genetic_diversity$datapath, check.names = F)
        genetic_divergences_csv_melt = reshape2::melt(genetic_divergences_csv)
        
        colnames(genetic_divergences_csv_melt) = c("file_name", "data_perc", "intra_dist")
        
        # CREATE PLOT
        genetic_divergence_plot <<- ggplot(data = genetic_divergences_csv_melt, aes(x = as.numeric( data_perc ), y = intra_dist)) + 
          geom_point() +
          geom_smooth(span = 3, col = input$genetic_divergence_line_col, fill = input$genetic_divergence_ci_col, alpha = input$genetic_divergence_ci_alpha) +
          scale_y_continuous(breaks = seq(0, max(genetic_divergences_csv_melt$intra_dist), by = input$genetic_divergence_y_interval)) +
          xlab("Data Percentage") +
          ylab(input$genetic_divergence_y_label) +
          ggthemes[[input$ggtheme_summary_plots]] +
          theme(legend.position = "bottom") +
          theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
          theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
          theme(axis.text.x = element_text(size = 12)) +
          theme(axis.text.y = element_text(size = 12)) +
          theme(legend.text=element_text(size=12)) +
          theme(plot.title = element_text(size=16, face = "bold")) +
          scale_x_continuous(breaks = as.numeric(genetic_divergences_csv_melt$data_perc), labels = genetic_divergences_csv_melt$data_perc) 
          
        
        output$summary_plot = renderPlot(genetic_divergence_plot)
        multiple_summary_plots$genetic_divergence <<- genetic_divergence_plot
        
      })
      
      # download the plot
      output$download_genetic_divergence_plot <- downloadHandler(
        
        filename = function (){paste("genetic_divergence_plot", input$plot_format_genetic_divergence, sep = '.')},
        
        content = function (file){
          
          width.genetic_divergence = as.numeric(input$w_plot_genetic_divergence) 
          height.genetic_divergence = as.numeric(input$h_plot_genetic_divergence) 
          dpi.genetic_divergence = as.numeric(input$res_plot_genetic_divergence)
          units.genetic_divergence = input$unit_plot_genetic_divergence
          
          ggsave(file, width = width.genetic_divergence, height = height.genetic_divergence, dpi = dpi.genetic_divergence, units = units.genetic_divergence,
                 genetic_divergence_plot)
        }
      )
      
    }) # end of observe
    
    ################################################################################################
  
  observe({
    
    # OVERSPLITTING PER MORPHOSPECIES GROUP
    infile_summary_oversplitting_morphospecies = input$summary_oversplitting_morphospecies
    
    if (is.null(infile_summary_oversplitting_morphospecies)) {
      return(NULL)
    }
    
    observeEvent(input$plot_summary_oversplits_per_morphospecies, {
      
      summary_oversplits_per_group_csv = read.csv(infile_summary_oversplitting_morphospecies$datapath, header = T, check.names = F)
      
      oversplitting_bar_summary <<- ggplot(data = summary_oversplits_per_group_csv, aes(x = predef_unique, y = Freq)) +
        geom_bar(stat = "summary", fill = input$oversplitting_morphospecies_col, col = "black") +
        expand_limits(y = 0) +
        xlab("Predefined group") +
        ylab("Mean oversplitting ratio per predefined \n group, on the full dataset (100%)") +
        ggthemes[[input$ggtheme_summary_plots]] +
        theme(axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 20, b = 0, l = 0))) +
        theme(axis.title.x = element_text(size = 14, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
        theme(axis.text.x = element_text(size = 12)) +
        theme(axis.text.y = element_text(size = 12)) +
        theme(plot.title = element_text(size=16, face = "bold")) +
        theme(axis.text.x = element_text(angle = req( input$x_axis_angle)) ) +
        geom_hline(yintercept= 1, linetype="dashed") 
      
      output$summary_plot = renderPlot(oversplitting_bar_summary)
      multiple_summary_plots$oversplitting <<- oversplitting_bar_summary
      
    })
    
    # download the plot
    output$download_summary_oversplitting_morphospecies_plot <- downloadHandler(
      
      filename = function (){paste("oversplitting_predefined_groups_plot", input$plot_format_summary_oversplitting_morphospecies, sep = '.')},
      
      content = function (file){
        
        width.oversplitting_morphospecies_summary = as.numeric(input$w_plot_summary_oversplitting_morphospecies) 
        height.oversplitting_morphospecies_summary = as.numeric(input$h_plot_summary_oversplitting_morphospecies) 
        dpi.oversplitting_morphospecies_summary = as.numeric(input$res_plot_summary_oversplitting_morphospecies)
        units.oversplitting_morphospecies_summary = input$unit_plot_summary_oversplitting_morphospecies
        
        ggsave(file, width = width.oversplitting_morphospecies_summary, height = height.oversplitting_morphospecies_summary, dpi = dpi.oversplitting_morphospecies_summary, units = units.oversplitting_morphospecies_summary,
               oversplitting_bar_summary)
      }
    )
    
  }) # end of observe
  
  observe({
    
    req(input$gridextra_ncol)
      
  observeEvent(input$plot_multiple_summary_plots, {
    
    output$summary_plot = renderPlot(gridExtra::grid.arrange(grobs = multiple_summary_plots, ncol = input$gridextra_ncol))
    
    # download the plot
    output$download_multiple_summary_plots <- downloadHandler(
      
      filename = function (){paste("mutlple_summary_plots", input$plot_format_multiple_summary_plots, sep = '.')},
      
      content = function (file){
        
        width.multiple_summary_plots = as.numeric(input$w_plot_multiple_summary_plots) 
        height.multiple_summary_plots = as.numeric(input$h_plot_multiple_summary_plots) 
        dpi.multiple_summary_plots = as.numeric(input$res_plot_multiple_summary_plots)
        units.multiple_summary_plots = input$unit_multiple_summary_plots
        
        ggsave(file, width = width.multiple_summary_plots, height = height.multiple_summary_plots, dpi = dpi.multiple_summary_plots, units = units.multiple_summary_plots,
               gridExtra::grid.arrange(grobs = multiple_summary_plots, ncol = input$gridextra_ncol))
      }
    )
    
  })
  
  }) # end of observe
  
  #########################################################################################################################
  # Upload multiple .csv files to extract a chosen column from each, and merge all those columns into one amalgamated file
  #########################################################################################################################
  
  observe({
    
    
    infile_multiple = input$amalgamate_multiple
    if (is.null( infile_multiple)) {
      return(NULL)
    } 
    # find how many files were uploaded, to use as a guide in a loop
    nfiles = nrow(infile_multiple) 
    # create an empty list to store each file in
    csv = list()
    column_names = c()
    
    for (i in 1:nfiles)
    {
      csv[[i]] = read.csv(infile_multiple$datapath[i])
      column_names[i] = i
    }
    
    updateSelectInput(session, "amalg_col", choices=colnames(csv[[1]][-1])) # the -1 removes the filename from the columns to choose from
    
    amalg_data = reactive({
      
      # store the selected column to extract from each file
      chosen_col = as.name( input$amalg_col )
      
      # create an empty dataframe to house the extracted columns
      desired_data = data.frame(matrix(nrow = nrow(csv[[1]]), ncol = nfiles)) 
      
      
      colnames(desired_data) = column_names
      rownames(desired_data) = csv[[1]]$filename
      
      # populate the desired_data dataframe with the merged columns from each file
      for (i in 1:nfiles){
        desired_data[,i] = csv[[i]][[chosen_col]]
      }
      
      return(desired_data)
      
    })
    
    # View the output on the screen
    observeEvent(input$view_amalg, {
      
      # output the merged dataframe to the screen
      output$amalgamate_table = renderTable(amalg_data(), rownames = TRUE)
      
    })
    
    # download the merged dataframe
    output$download_amalg = downloadHandler(
      filename = function (){paste(input$amalg_col, 'csv', sep = '.')},
      content = function (file){write.csv(amalg_data(), file, row.names = TRUE)}
    )
    
  }) # end of observe 
    
    
  #########################################################################################################################
  # END OF APPLICATION
  #########################################################################################################################
  
}
