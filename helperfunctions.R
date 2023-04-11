library(tidyverse,quietly=TRUE,verbose = FALSE,warn.conflicts = FALSE)
library(bnlearn,quietly=TRUE,verbose = FALSE,warn.conflicts = FALSE)
library(pcalg,quietly=TRUE,verbose = FALSE,warn.conflicts = FALSE)
library(parallel,quietly=TRUE,verbose = FALSE,warn.conflicts = FALSE)
library(LocalFCI,quietly=TRUE,verbose = FALSE,warn.conflicts = FALSE)

# DAG Info -------------------------------------------------------------

# Build the file structure for this simulation setting and get network info
# Obtain network information and simulated data
get_network_DAG <- function(net){
  current_wd <- getwd() # Save the current working directory to return to at the end
  setwd(rds_dir) # Go to the directory with the network .rds files
  
  if (paste0(net,".rds") %in% list.files()){
    network <- readRDS(paste0(net,".rds"))
    p <- length(names(network))
    node_names <- names(network)
  } else {
    stop('Invalid Network Name')
  }
  
  true_dag <- amat(network) # Get the adjacency matrix of the true DAG
  
  network_cpdag <- cpdag(network) # Get the CPDAG of the true DAG (bnlearn)
  true_cpdag <- amat(network_cpdag) # Adjacency matrix of the true CPDAG
  
  setwd(current_wd) # Return to original working directory
  
  return(list(
    "net"=net,
    "p"=p,
    "node_names"=node_names,
    "network"=network,
    "true_dag"=true_dag,
    "cpdag"=true_cpdag
  ))
}


# Simulation Data Functions -----------------------------------------------

# Function responsible for generating the data for the simulation setting
simulation_data_creation <- function(){
  # Check to see if the data has been generated already
  go_to_dir("data")
  sims_not_created <- check_sims_created()
  # Generate Data
  if (sims_not_created){
    sims_text_output()
    gdg <- generate.data.grid(data.grid,
                              out.dir=getwd(),
                              array_num=array_num,
                              verbose=FALSE,
                              path.start=home_dir)
    
    cat("finished",file = paste0(net,"_",array_num,"/marker.txt"))
  } else {
    Sys.sleep(3)
  }
  setwd("..") # Return to the original directory
}

# Check whether or not simulations have been created so they do not have to be created again
check_sims_created <- function(){
  sim_file <- paste0(net,"_",array_num)
  if (dir.exists(sim_file)){
    files <- list.files(sim_file) # Get the list of files in this directory\
    sims_not_created <- !any(str_detect(files,"data[[:digit:]]+.txt"))
  } else {
    sims_not_created <- TRUE
  }
  return(sims_not_created)
  
}

# Provides a text output for which datasets are being created
sims_text_output <- function(){
  q1 <- format(paste("Generating Datasets for",net),width = 35,justify = "left")
  q2 <- format(paste("U.B. Variance:",data.grid$high),width = 20,justify = "centre")
  q3 <- format(paste("U.B. Coefs:",data.grid$ub),width = 20,justify = "centre")
  q4 <- format(paste("n =",data.grid$n.obs),width = 10,justify = "centre")
  cat(q1,"|",q2,"|",q3,"|",q4,"\n")
}

# Grab dataframe for network at given sample size
grab_data <- function(df_num){
  go_to_dir("data")
  go_to_dir(paste0(net,"_",array_num))
  check <- check_file(paste0("data",df_num,".txt"))
  if (!check){
    cat("The data did not pass the check. File size:",file.size(paste0("data",df_num,".txt")),"\n")
  } else {
    cat("The data did pass the check. File size:",file.size(paste0("data",df_num,".txt")),"\n")
  }
  df <- read.table(paste0("data",df_num,".txt"))
  colnames(df) <- network_info$node_names
  setwd("../..")
  return(df)
}

# Checks whether or not the file with the data for the current setting exists
check_file <- function(f_name){
  files <- list.files()
  i <- 0
  while (!(file.exists("marker.txt")) & i < 20){
    i <- i + 1
    Sys.sleep(3)
    files <- list.files()
  }
  return(f_name %in% list.files())
}

# Target functions --------------------------------------------------------

# Function to obtain targets of various sizes for the network \
check_targets_defined_get_targets <- function(net_info){
  # Determine whether or not targets have been defined
  if (!file.exists(paste0("targets_",net,".rds"))){
    targets <- get_targets(net_info$p)
    saveRDS(targets,paste0("targets_",net,".rds"))
  } else {
    targets <- readRDS(paste0("targets_",net,".rds"))
    cat("Targets for",net,"obtained from file.\n")
  }
  
  return(targets)
}

# Internal function to create the list of targets
get_targets <- function(p){
  max_targets_per_category <- 8
  if (p<max_targets_per_category){
    max_targets_per_category <- p
  }
  # Generate a list of matrices containing targets sets of various sizes
  targets <- lapply(1:max_targets,function(num){
    if (p<200){
      target_options <- combn(1:p,num)
      return(target_options[,sample(1:ncol(target_options),max_targets_per_category)])
    } else {
      targets <- matrix(nrow = num,ncol = max_targets_per_category)
      for (j in 1:max_targets_per_category){
        targets[,j] <- sample(1:p,num)
      }
      return(targets)
    }
  })
  # Create final list containing all combinations of targets
  final_target_list <- list()
  # Add individual targets first
  final_target_list <- c(final_target_list,targets[[1]]) 
  offset <- length(final_target_list)
  j <- 1
  # Add remaining sets of targets
  for (i in 2:max_targets){
    while (j < ncol(targets[[i]])){
      final_target_list[[offset+j]]<-targets[[i]][,j]
      j <- j + 1
    }
    offset <- length(final_target_list)
    j <- 1
  }
  
  return(final_target_list)
}


# Global PC -----------------------------------------------------

# Run global PC 
run_global_pc <- function(df,trial_num){
  
  # Lists to store results
  time_diff <- list()
  lmax_list <- list()
  num_tests <- list()
  
  if (file.exists(paste0("pc_",array_num,"_",trial_num,"_results.rds"))){
    cat("We are loading PC algorithm results from a saved file ... ")
    return(readRDS(paste0("pc_",array_num,"_",trial_num,"_results.rds")))
  }
  
  largest_possible_sepset <- 5
  pc_test_file <- paste0("pc_",array_num,"_tests.txt")
  sink(file = pc_test_file)
  start <- Sys.time()
  pc.fit <- as(pc(suffStat = list(C = cor(df), n = n),
                  indepTest = gaussCItest, ## indep.test: partial correlations
                  alpha=alpha, labels = network_info$node_names,
                  verbose = TRUE,m.max=largest_possible_sepset),"amat")
  end <- Sys.time()
  sink(file = NULL)
  diff <- end - start
  units(diff) <- "secs"
  time_diff$PC <- diff
  lmax_list$PC <- get_lmax(pc_test_file)
  num_tests$PC <- get_pc_test_num(pc_test_file)
  
  pc_res <- list(
    "pc"=pc.fit,
    "time_diff"=time_diff,
    "lmax"=lmax_list,
    "num_tests"=num_tests
  )
  saveRDS(pc_res,paste0("pc_",array_num,"_",trial_num,"_results.rds"))

  return(pc_res)
}




# PC Helpers --------------------------------------------------------------

# Get the maximum size of separating sets used in PC algorithm
get_lmax <- function(file){
  # Read the file with PC algorithm output
  pc_file <- read_file(file = file)
  # Only identify lines with tests, which will include a part including the
  # separating set
  sep_sets <- unlist(str_extract_all(pc_file,'S= .* :'))
  # Count the number of nodes in the separating set from each test
  counts <- sapply(sep_sets,function(s){
    res <- str_count(s,"[0-9]+")
  })
  
  if (is.list(counts)){
    counts <- unlist(counts)
    if (length(counts)==0){
      counts <- 0
    }
  }
  # Return the maximum separating set size observed in the file
  return(max(counts))
}

get_pc_test_num <- function(file){
  # Read the PC algorithm output file
  pc_file <- read_file(file = file)
  # Each p-value recorded in the file represents a test
  tests <- unlist(str_extract_all(pc_file,'pval = [0-9]+'))
  # Return the number of p-values observed
  return(length(tests))
}


# Local FCI ---------------------------------------------------------------

run_fci_target <- function(t,df,num,results_pc,algo,curr_dir){
  output_text(t,num,algo,"lfci")
  setwd(curr_dir)
  # vars <- create_target_directory(t) 
  # Run local FCI
  results <- run_local_fci(t,df,num,results_pc,algo)
  
  return(results)
}

run_pc_target <- function(t,df,num,results_pc,algo,curr_dir){
  output_text(t,num,algo,"lpc")
  setwd(curr_dir)
  # vars <- create_target_directory(t) 
  # Run local PC
  results <- run_local_pc(t,df,num,results_pc,algo)
  
  return(results)
}

# Local FCI Helpers -------------------------------------------------------

run_local_fci <- function(t,df,num,results_pc,algo){
  # Set initial values for lmax and start timer
  val_lmax <- results_pc$lmax$PC
  lmax <- max(1,val_lmax)
  results_pc$lmax[["Local FCI"]] <- lmax
  true_dag <- network_info$true_dag
  start <- Sys.time()
  # Run local FCI
  localfci_result <- localfci(data=df,
                              targets=t,
                              lmax=lmax,
                              tol=alpha,
                              mb_tol=mb_alpha,
                              method=algo,
                              verbose = FALSE)
  end <- Sys.time()
  # Record values
  diff <- end - start
  units(diff) <- "secs"
  results_pc$time_diff[["Local FCI"]] <- diff
  results_pc$num_tests[["Local FCI"]] <- localfci_result$NumTests
  check_anc_edge_markings(localfci_result$amat,
                          paste0(result_dir,"/anc_edge_errors.txt"),
                          array_num)
  results <- neighborhood_results(t,localfci_result,results_pc,num)
  return(results)
}

run_local_pc <- function(t,df,num,results_pc,algo){
  # Set initial values for max. sep. set size and start timer
  val_lmax <- results_pc$lmax$PC
  lmax <- max(1,val_lmax)
  true_dag <- network_info$true_dag
  start <- Sys.time()
  # Run local PC
  localpc_result <- localpc(data=df,
                            targets=t,
                            lmax=lmax,
                            tol=alpha,
                            mb_tol = mb_alpha,
                            method = algo,
                            verbose = FALSE)
  end <- Sys.time()
  diff <- end - start
  units(diff) <- "secs"
  localpc_result$time_diff <- diff
  localpc_result$lmax<- lmax
  results <- neighborhood_results_pc(t,localpc_result,num)
  
  return(results)
}

create_target_directory <- function(t){
  go_to_dir(paste0("Target",ifelse(length(t)>1,"s",""),"=",paste(t,collapse = ",")))
  return()
}


# Result Functions --------------------------------------------------------

get_nbhds <- function(t,local_result){
  # Gather all the neighborhood nodes (true)
  nbhd <- as.vector(sapply(t,function(targ){
    check_neighbors_retrieval(network_info$p,network_info$node_names,
                              network_info$true_dag,targ-1)+1
  }))
  ### First, we look at the results using only the true neighborhoods (narrow)
  if (is.list(nbhd)) nbhd <- unlist(nbhd)
  # Add the target to the neighborhood set
  nbhd <- sort(unique(c(t,nbhd)))
  
  ### Second, we consider the results on the union of the estimated neighborhood
  ### nodes and the true nodes
  nbhd_broad <- union(nbhd,local_result$Nodes)
  
  return(
    list(
      "narrow"=nbhd,
      "broad"=nbhd_broad
    )
  )
}

# Return the matrices we need to compare for our metrics
get_ref_mats <- function(nbhds,pc_results,local_result,get_pc=TRUE){
  # Get PC matrix around narrow neighborhood
  pc_mat_narrow<-NULL
  pc_mat_broad<-NULL
  if (get_pc){
    pc_mat_narrow <- matrix(pc_results$pc,
                            nrow = network_info$p)[nbhds$narrow,nbhds$narrow]
    pc_mat_broad <- matrix(pc_results$pc,
                           nrow = network_info$p)[nbhds$broad,nbhds$broad]
  }
  
  # Remove all ancestral edges from consideration
  local_mat <- local_result$amat
  n_dim <- nrow(local_mat)
  for (i in 1:n_dim){
    for (j in i:n_dim){
      if (local_mat[i,j]>1){
        if (local_mat[j,i]<=1){
          cat("Local matrix i,j val:",local_mat[i,j],"\n")
          cat("Local matrix j,i val:",local_mat[j,i],"\n")
          warning("Improperly set ancestral edge from local algorithm")
        }
        local_mat[i,j] <- local_mat[j,i] <- 0
      }
    }
  }
  
  # Ground truth graph from narrow neighborhood
  # subgraph of CPDAG is Ground Truth
  true_neighborhood_graph_narrow <- network_info$cpdag[nbhds$narrow,nbhds$narrow] 
  # local algorithm's narrow graph
  local_mat_narrow <- local_result$amat[nbhds$narrow,nbhds$narrow]
  local_mat_narrow_shd <- local_mat[nbhds$narrow,nbhds$narrow]
  
  ### Broad neighborhood graphs
  true_neighborhood_graph_broad <- network_info$cpdag[nbhds$broad,nbhds$broad]
  local_mat_broad <- local_result$amat[nbhds$broad,nbhds$broad]
  local_mat_broad_shd <- local_mat[nbhds$broad,nbhds$broad]
  
  ### Fix matrices if only one node
  if (length(nbhds$narrow)==1){
    if(get_pc) {pc_mat_narrow <- as.matrix(pc_mat_narrow,nrow=1,ncol=1)}
    local_mat_narrow <- as.matrix(local_mat_narrow,nrow=1,ncol=1)
    true_neighborhood_graph_narrow <- as.matrix(true_neighborhood_graph_narrow,
                                                nrow=1,ncol=1)
  }
  if (length(nbhds$broad)==1){
    if (get_pc) {pc_mat_broad <- as.matrix(pc_mat_broad,nrow=1,ncol=1)}
    local_mat_broad <- as.matrix(local_mat_broad,nrow=1,ncol=1)
    true_neighborhood_graph_broad <- as.matrix(true_neighborhood_graph_broad,
                                               nrow=1,ncol=1)
  }
  return(
    list(
      "ground_narrow"=true_neighborhood_graph_narrow,
      "ground_broad"=true_neighborhood_graph_broad,
      "pc_narrow"=pc_mat_narrow,
      "pc_broad"=pc_mat_broad,
      "local_narrow"=local_mat_narrow,
      "local_broad"=local_mat_broad,
      "local_narrow_shd"=local_mat_narrow_shd,
      "local_broad_shd"=local_mat_broad_shd
    )
  )
}

get_result_metrics <- function(targets,graphs,nbhds,method='lfci'){
  if (method=='pc'){
    g_narrow <- graphs$pc_narrow
    g_broad <- graphs$pc_broad
  } else {
    g_narrow <- graphs$local_narrow
    g_broad <- graphs$local_broad
  }
  # Get these results for local FCI
  results_narrow <- allMetrics(g_narrow,
                               graphs$ground_narrow,
                               sapply(targets,function(t) {
                                 which(nbhds$narrow==t)-1
                               }),
                               network_info$true_dag,
                               nbhds$narrow-1,
                               algo=method,which_nodes = "narrow",verbose = FALSE) %>% 
    select(contains("skel"),contains("_v_"),contains("F1"))
  results_broad <- allMetrics(g_broad,
                              graphs$ground_broad,
                              sapply(targets,function(t) {
                                which(nbhds$broad==t)-1
                              }),
                              network_info$true_dag,
                              nbhds$broad-1,
                              algo=method,which_nodes = "broad",verbose = FALSE)
  results_broad <- results_broad %>% 
    select(contains("skel"),contains("_v_"),contains("pra"),contains("ancestor"))
  
  if (method=='pc'){
    results_broad <- results_broad %>% select(-contains("ancestor"))
  }
  
  results <- cbind(results_narrow,results_broad)
  return(results)
}

# Return neighborhood information as well as Markov Blanket Recovery and FCI rules
# usage metrics
additional_metrics <- function(graphs,local_result,pc_result,
                               targets,nbhds,method='lfci'){
  if (method=='lfci') {
    ### General information about the neighborhood (metrics.cpp)
    nbhd_metrics <- getNeighborhoodMetrics(graphs$ground_narrow)
    ### Markov Blanket Recovery metrics
    # (mbEst.R) Markov Blanket Recovery TP, FN, FP
    mb_metrics <- mbRecovery(network_info$cpdag,local_result$referenceDAG,targets)
    # MB children FN, TP; MB parents FN, TP; MB spouses FN, TP; MB Total FP
    mb_metrics_add <- mbRecoveryMetrics(network_info$cpdag,
                                        local_result$referenceDAG,targets)
    # Time and Tests given by MXM function
    mb_time <- getTotalMBTime(local_result$mbList)
    mb_tests <- getTotalMBTests(local_result$mbList)
    
    ### Statistics for the amount of times each FCI rule was used
    rule_usage <- t(data.frame(local_result$RulesUsed))
    dimnames(rule_usage) <- list(NULL,paste0("rule",0:10))
    
    ### Combine results so far
    results <- cbind(nbhd_metrics,mb_metrics,mb_metrics_add,
                     mb_time,mb_tests,rule_usage)
    ### Add previously collected results about tests, timing, lmax, and nbhds
    results <- results %>% 
      mutate(pc_num_tests=pc_result$num_tests[["PC"]],
             lfci_num_tests=pc_result$num_tests[["Local FCI"]],
             pc_time=pc_result$time_diff[["PC"]],
             lfci_time=pc_result$time_diff[["Local FCI"]]) %>% 
      mutate(pc_lmax = pc_result$lmax$PC,
             lfci_lmax=pc_result$lmax[["Local FCI"]]) %>%
      mutate(totalMBEstTime=mb_time,
             totalMBEstTimeInclusive=local_result$mbEstTime,
             totalMBTests=mb_tests,
             totalSkeletonTime=local_result$totalSkeletonTime,
             targetSkeletonTimes=paste(local_result$targetSkeletonTimes,collapse = ","),
             totalcpptime=local_result$totalTime,
             nodes_narrow=paste(nbhds$narrow,collapse = ","),
             nodes_broad=paste(nbhds$broad,collapse = ","),
             nodes_local_est=paste(local_result$Nodes,collapse = ",")
      )
  } else {
    ### Simplified results data for local PC
    results <- data.frame(
      lpc_time=as.numeric(local_result$time_diff),
      lpc_lmax=local_result$lmax,
      lpc_num_tests=local_result$NumTests,
      lpc_rule0=local_result$rules_used[1],
      lpc_rule1=local_result$rules_used[2],
      lpc_rule2=local_result$rules_used[3],
      lpc_rule3=local_result$rules_used[4],
      lpc_rule4=local_result$rules_used[5],
      targetSkeletonTimes=paste(local_result$targetSkeletonTimes,collapse = ","),
      totalcpptime=local_result$totalTime
    )
    
  }
  
  return(results)
}

# Add values for number of tests and timing (everything should be in seconds),
# Maximum size of separating set, simulation setting values, targets, trial number,
# Network information, timing and testing for MB estimation, timing breakdown for steps
clean_modify_results <- function(results,num,t){
  results <- results %>% 
    mutate(sim_number=array_num,
           alpha=alpha,
           mb_alpha=mb_alpha,
           algo=algo,
           net=net,
           n=n,
           ub=ub,
           high=high,
           trial_num=num,
           num_targets=length(t),
           targets=paste(t,collapse = ","),
           p=network_info$p,
           net_edges=sum(network_info$true_dag)) %>%
    select(sim_number,net,p,net_edges,trial_num,n,alpha,mb_alpha,algo,ub,high,
           num_targets,targets,contains("nodes"),size,num_edges,
           contains("lfci"),contains("mb",ignore.case=FALSE),contains("pc"),
           contains("time"),contains("tests"),starts_with("rule"))
  return(results)
}

# TODO: Fix this 
# Idea: remove all ancestral edges and go from there
get_shd <- function(ref_mat,est_mat){
  p <- nrow(ref_mat)
  ref_g <- bnlearn::empty.graph(nodes=as.character(1:p))
  est_g <- bnlearn::empty.graph(nodes=as.character(1:p))
  colnames(ref_mat) <- rownames(ref_mat) <- as.character(1:p)
  colnames(est_mat) <- rownames(est_mat) <- as.character(1:p)
  amat(ref_g) <- ref_mat
  amat(est_g,check.cycles=FALSE) <- est_mat
  # if (!bnlearn::acyclic(est_g)){
  #    
  # }
  shd_res <- bnlearn::shd(learned = est_g,true = ref_g)
  return(shd_res)
}

# Compile all results about the simulation for local FCI
neighborhood_results <- function(t,localfci_result,pc_results,num){
  # Obtain narrow and broad neighborhoods
  nbhd <- get_nbhds(t,localfci_result)
  # Obtain narrow and broad matrices from algorithms and reference
  mats <- get_ref_mats(nbhd,pc_results,localfci_result)
  
  # Compare results
  # Returns the following:
  # Skeleton FP, FN, TP
  # V-Structure FP, FN, TP
  # PRA FP, FN, TP, Potential (Undirected edge in at least one of the graphs)
  # Ancestors Correct, Missing, Missing Orientation, Reversed Orientation, Falsely Oriented Edge, Falsely added connection
  # Overall F1 Score
  lfci_metrics <- get_result_metrics(t,mats,nbhd)
  pc_metrics <- get_result_metrics(t,mats,nbhd,'pc')
  if (length(nbhd)>0){
    pc_metrics$pc_narrow_shd <- get_shd(mats$ground_narrow,mats$pc_narrow)
    pc_metrics$pc_broad_shd <- get_shd(mats$ground_broad,mats$pc_broad)
    lfci_metrics$lfci_narrow_shd <- get_shd(mats$ground_narrow,mats$local_narrow_shd)
    lfci_metrics$lfci_broad_shd <- get_shd(mats$ground_broad,mats$local_broad_shd)
  } else {
    pc_metrics$pc_narrow_shd <- 0
    pc_metrics$pc_broad_shd <- 0
    lfci_metrics$lfci_narrow_shd <- 0
    lfci_metrics$lfci_broad_shd <- 0
  }
  
  
  add_results <- additional_metrics(mats,localfci_result,pc_results,t,nbhd)
  
  # Combine results
  results <- cbind(lfci_metrics,pc_metrics,add_results)
  results <- clean_modify_results(results,num,t)
  
  return(results)
}

neighborhood_results_pc <- function(t,localpc_result,num){
  # Obtain narrow and broad neighborhoods
  nbhd <- get_nbhds(t,localpc_result)
  # Obtain narrow and broad matrices from algorithms and reference
  mats <- get_ref_mats(nbhd,NULL,localpc_result,get_pc = FALSE)

  lpc_results <- get_result_metrics(t,mats,nbhd,method = 'lpc')
  lpc_results$lpc_narrow_shd <- get_shd(mats$ground_narrow,mats$local_narrow)
  lpc_results$lpc_broad_shd <- get_shd(mats$ground_broad,mats$local_broad)
  add_results <- additional_metrics(mats,localpc_result,pc_results,
                                    t,nbhd,method = 'lpc') %>% 
    mutate(trial_num=num,
           targets=paste(t,collapse = ","))
  results <- cbind(lpc_results,add_results)

  return(results)
}

# General Functions -------------------------------------------------------

go_to_dir <- function(file){
  # cat("Current dir:\n",list.dirs(),"\n\n")
  # cat("Objective dir:\n",list.dirs(),"\n\n")
  if (!dir.exists(file)){
    dir.create(file)
  }
  setwd(file)
}

find_neighbors <- function(t,true_dag){
  # Parents of node t
  parents <- which(true_dag[,t]==1)
  # Children of node t
  children <- which(true_dag[t,]==1)
  # Spouses share the same children
  spouses <- c()
  for (c in children){ # loop through all of t's children
    # Find all parents of current child
    potential <- which(true_dag[,c]==1)
    # Add all parents that are not t to spouses
    spouses <- c(spouses,potential[potential!=t])
  }
  return(unique(c(parents,children,spouses)))
}

find_all_neighbors <- function(t,true_dag){
  if (length(t)==1){
    return(find_neighbors(t,true_dag))
  } else {
    neighbors <- lapply(t,find_neighbors,true_dag=true_dag)
    all_neighbors <- sort(unique(unlist(neighbors)))
    return(all_neighbors)
  }
}

# Check whether these simulations have been completed
check_completion <- function(dir){
  if (!dir.exists(dir)){
    return(NULL)
  } else {
    current_wd <- getwd()
    setwd(dir)
    rds_files <- list.files(pattern = "*.rds$")
    file_l <- length(rds_files)
    if (file_l==0){
      setwd(current_wd)
      return(NULL)
    } else if (file_l==1){
      cat("Simulations Completed. Name of file:",rds_files," ")
      results <- readRDS(rds_files)
      cat("Number of rows from file:",nrow(results),"\n")
      setwd(current_wd)
      return(results)
    } else {
      stop("Something wrong with simulation completion check code.")
    }
  }
}

output_text <- function(t,num,algo,method){
  q1 <- format(paste("Tols.:",paste0(alpha,",",mb_alpha)),width = 21,justify = "centre")
  q2 <- format(paste("Net:",net),width = 25,justify = "centre")
  q3 <- format(paste("Sample N:",data.grid$n.obs),width = 15,justify = "centre")
  q4 <- format(paste("Algo:",method),width = 10,justify = "centre")
  q5 <- format(paste("MB Algo:",algo),width = 14,justify = "centre")
  q6 <- format(paste("Targets:",paste(t,collapse=",")),width = 23,justify = "left")
  q7 <- format(paste("Dataset:",num),width = 15,justify = "centre")
  cat(q1,"|",q2,"|",q3,"|",q4,"|",q5,"|",q6,"|",q7,"\n")
}

check_anc_edge_markings <- function(G,file_name,sim_num){
  for (i in 1:nrow(G)){
    for (j in 1:ncol(G)){
      if ((G[i,j] %in% c(2,3)) && (G[j,i]==1)){
        cat("Incorrect ancestral labels for",sim_num,"\n",file = file_name)
        return()
      }
    }
  }
}
