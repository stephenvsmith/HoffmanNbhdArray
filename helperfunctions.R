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
  
  network_cpdag <- cpdag(network) # Get the CPDAG of the true DAG
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

simulation_data_creation <- function(){
  # Check to see if the data has been generated already
  go_to_dir("data")
  sims_not_created <- check_sims_created(n)
  # Generate Data
  if (sims_not_created){
    sims_text_output()
    gdg <- generate.data.grid(data.grid,
                              out.dir=getwd(),
                              array_num=array_num,
                              verbose=FALSE,
                              path.start=home_dir)
    # dir.create(paste0(net,"_",array_num))
    # data_files <- list.files(paste0(net,"; n = ",n,"; c = 0"))
    # sapply(data_files,function(f){
    #   file.copy(paste0(net,"; n = ",n,"; c = 0/",f),
    #             paste0(net,"_",array_num))
    # })
    cat("finished",file = paste0(net,"_",array_num,"/marker.txt"))
    #unlink(paste0(net,"; n = ",n,"; c = 0"),recursive = TRUE)
  } else {
    Sys.sleep(60)
  }
  setwd("..") # Return to the original directory
}

# Check whether or not simulations have been created so they do not have to be created again
check_sims_created <- function(n){
  sim_file <- paste0(net,"_",array_num)
  if (dir.exists(sim_file)){
    files <- list.files(sim_file) # Get the list of files in this directory\
    sims_not_created <- !any(str_detect(files,"data[[:digit:]]+.txt"))
  } else {
    sims_not_created <- TRUE
  }
  return(sims_not_created)
  
}

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

get_targets <- function(p){
  max_targets_per_category <- 15
  if (p<25){
    max_targets_per_category <- p
  }
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
  final_target_list <- list()
  final_target_list <- c(final_target_list,targets[[1]]) # Add individual targets
  offset <- length(final_target_list)
  j <- 1
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
run_global_pc <- function(df){
  
  # Lists to store results
  time_diff <- list()
  lmax_list <- list()
  num_tests <- list()
  
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
  units(diff) <- "mins"
  time_diff$PC <- diff
  lmax_list$PC <- get_lmax(pc_test_file)
  num_tests$PC <- get_pc_test_num(pc_test_file)
  
  return(list(
    "pc"=pc.fit,
    "time_diff"=time_diff,
    "lmax"=lmax_list,
    "num_tests"=num_tests))
}




# PC Helpers --------------------------------------------------------------

# Get the maximum size of separating sets used in PC algorithm
get_lmax <- function(file){
  pc_file <- read_file(file = file)
  
  sep_sets <- unlist(str_extract_all(pc_file,'S= .* :'))
  counts <- sapply(sep_sets,function(s){
    res <- str_count(s,"[0-9]+")
  })
  
  if (is.list(counts)){
    counts <- unlist(counts)
    if (length(counts)==0){
      counts <- 0
    }
  }
  
  return(max(counts))
}

get_pc_test_num <- function(file){
  pc_file <- read_file(file = file)
  tests <- unlist(str_extract_all(pc_file,'pval = [0-9]+'))
  return(length(tests))
}


# Local FCI ---------------------------------------------------------------

run_fci_target <- function(t,df,num,results_pc,algo,curr_dir){
  output_text(t,num,algo)
  setwd(curr_dir)
  # vars <- create_target_directory(t) 
  # Run local FCI
  results <- run_local_fci(t,df,num,results_pc,algo)
  
  return(results)
}

run_pc_target <- function(t,df,num,results_pc,algo,curr_dir){
  output_text(t,num,algo)
  setwd(curr_dir)
  # vars <- create_target_directory(t) 
  # Run local PC
  results <- run_local_pc(t,df,num,results_pc,algo)
  
  return(results)
}

# Local FCI Helpers -------------------------------------------------------

run_local_fci <- function(t,df,num,results_pc,algo){
  val_lmax <- results_pc$lmax$PC
  lmax <- max(1,val_lmax)
  results_pc$lmax[["Local FCI"]] <- lmax
  true_dag <- network_info$true_dag
  start <- Sys.time()
  
  localfci_result <- localfci(data=df,
                              targets=t,
                              lmax=lmax,
                              tol=alpha,
                              mb_tol=mb_alpha,
                              method=algo,
                              verbose = FALSE)
  end <- Sys.time()
  diff <- end - start
  units(diff) <- "mins"
  results_pc$time_diff[["Local FCI"]] <- diff
  results_pc$num_tests[["Local FCI"]] <- localfci_result$NumTests
  results <- neighborhood_results(t,localfci_result,results_pc,num)
  return(results)
}

run_local_pc <- function(t,df,num,results_pc,algo){
  val_lmax <- results_pc$lmax$PC
  lmax <- max(1,val_lmax)
  true_dag <- network_info$true_dag
  start <- Sys.time()
  localpc_result <- localpc(data=df,
                            targets=t,
                            lmax=lmax,
                            tol=alpha,
                            mb_tol = mb_alpha,
                            method = algo,
                            verbose = FALSE)
  end <- Sys.time()
  diff <- end - start
  units(diff) <- "mins"
  localpc_result$time_diff <- diff
  localpc_result$lmax<- lmax
  results <- neighborhood_results_pc(t,localpc_result,num)
  
  return(results)
}

create_target_directory <- function(t){
  go_to_dir(paste0("Target",ifelse(length(t)>1,"s",""),"=",paste(t,collapse = ",")))
  return()
}

# Compile all results about the simulation
neighborhood_results <- function(t,localfci_result,pc_results,num){
  nbhd <- as.vector(sapply(t,function(targ){
    check_neighbors_retrieval(network_info$p,network_info$node_names,network_info$true_dag,targ-1)+1
  }))
  if (is.list(nbhd)) nbhd <- unlist(nbhd)
  nbhd <- sort(unique(c(t,nbhd)))
  # Zoom in on estimated and true DAGs (only the target and first-order neighbors)
  nodes_zoom <- network_info$node_names[nbhd]
  pc_mat <- matrix(pc_results$pc,nrow = network_info$p)[nbhd,nbhd]
  true_neighborhood_graph <- network_info$cpdag[nbhd,nbhd] # CPDAG is Ground Truth
  localfci_mat <- localfci_result$amat[nbhd,nbhd]
  # Compare results
  if (length(nbhd)==1){
    pc_mat <- as.matrix(pc_mat,nrow=1,ncol=1)
    localfci_mat <- as.matrix(localfci_mat,nrow=1,ncol=1)
    true_neighborhood_graph <- as.matrix(true_neighborhood_graph,nrow=1,ncol=1)
  }
  results <- allMetrics(localfci_mat,
                        true_neighborhood_graph,
                        sapply(t,function(tg) {which(nbhd==tg)-1}),
                        algo="lfci",ref="sub_cpdag",verbose = FALSE) 
  results_pc <- allMetrics(pc_mat,
                           true_neighborhood_graph,
                           sapply(t,function(tg) {which(nbhd==tg)-1}),
                           algo = "pc",ref = "sub_cpdag",verbose = FALSE)
  nbhd_metrics <- getNeighborhoodMetrics(true_neighborhood_graph)
  mb_metrics <- mbRecovery(network_info$cpdag,localfci_result$referenceDAG,t)
  mb_metrics_add <- mbRecoveryMetrics(network_info$cpdag,localfci_result$referenceDAG,t)
  mb_time <- getTotalMBTime(localfci_result$mbList)
  mb_tests <- getTotalMBTests(localfci_result$mbList)
  rule_usage <- t(data.frame(localfci_result$RulesUsed))
  dimnames(rule_usage) <- list(NULL,paste0("rule",0:10))
  results <- cbind(nbhd_metrics,results,rule_usage,
                   results_pc,mb_metrics,mb_metrics_add)
  results <- results %>% 
    mutate(pc_num_tests=pc_results$num_tests[["PC"]],
           lfci_num_tests=pc_results$num_tests[["Local FCI"]],
           pc_time=pc_results$time_diff[["PC"]],
           lfci_time=pc_results$time_diff[["Local FCI"]]) %>%
    mutate(pc_lmax = pc_results$lmax$PC,
           lfci_lmax=pc_results$lmax[["Local FCI"]]) %>%
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
    mutate(totalMBEstTime=mb_time,
           totalMBTests=mb_tests,
           totalSkeletonTime=localfci_result$totalSkeletonTime,
           targetSkeletonTimes=paste(localfci_result$targetSkeletonTimes,collapse = ","),
           totalcpptime=localfci_result$totalTime,
           nodes=paste(localfci_result$Nodes,collapse = ","),
           true_nodes=paste(nbhd,collapse = ",")
    )
  
  # saveRDS(results,file = paste0("results_df",num,".rds"))
  
  # capture.output(results %>% select(size,num_edges,contains("pc")),
  #                file = paste0("results_pc",num,".txt"))
  # capture.output(results %>% select(size,num_edges,contains("fci")),
  #                file = paste0("results_fci",num,".txt"))
  
  # write.table(localfci_result$amat,paste0("estAmat",num,".txt"))
  # write.table(localfci_result$referenceDag,"refDAG.txt")
  # saveRDS(localfci_result$S,paste0("SepSet",num,".rds"))
  # saveRDS(localfci_result$mbList,paste0("mbList",num,".rds"))
  # saveRDS(localfci_result$data_means,paste0("dataMeans",num,".rds"))
  # saveRDS(localfci_result$data_cov,paste0("dataCov",num,".rds"))
  
  # saveRDS(mb_metrics,paste0("mbMetrics",num,".rds"))
  
  
  # Testing Diagnostics
  
  # capture.output(results %>% select(size,num_edges,contains("pc")),
  #                file = paste0("results_pc",num,".txt"))
  # capture.output(results %>% select(size,num_edges,contains("fci")),
  #                file = paste0("results_fci",num,".txt"))
  
  # write.table(localfci_result$amat,paste0("estAmat",num,".txt"))
  # write.table(localfci_result$referenceDag,"refDAG.txt")
  # saveRDS(localfci_result$S,paste0("SepSet",num,".rds"))
  # saveRDS(localfci_result$mbList,paste0("mbList",num,".rds"))
  # saveRDS(localfci_result$data_means,paste0("dataMeans",num,".rds"))
  # saveRDS(localfci_result$data_cov,paste0("dataCov",num,".rds"))
  
  return(results)
}

neighborhood_results_pc <- function(t,localpc_result,num){
  nbhd <- as.vector(sapply(t,function(targ){
    check_neighbors_retrieval(network_info$p,network_info$node_names,network_info$true_dag,targ-1)+1
  }))
  if (is.list(nbhd)) nbhd <- unlist(nbhd)
  nbhd <- sort(unique(c(t,nbhd)))
  # Zoom in on estimated and true DAGs (only the target and first-order neighbors)
  nodes_zoom <- network_info$node_names[nbhd]
  true_neighborhood_graph <- network_info$cpdag[nbhd,nbhd] # CPDAG is Ground Truth
  localpc_mat <- localpc_result$amat[nbhd,nbhd]
  # Compare results
  if (length(nbhd)==1){
    localpc_mat <- as.matrix(localpc_mat,nrow=1,ncol=1)
    true_neighborhood_graph <- as.matrix(true_neighborhood_graph,nrow=1,ncol=1)
  }
  
  results <- allMetrics(localpc_mat,
                         true_neighborhood_graph,
                         sapply(t,function(tg) {which(nbhd==tg)-1}),
                         algo = "lpc",ref = "sub_cpdag",verbose = FALSE) # Need to check out the sapply here
  results <- results %>%
    mutate(lpc_num_tests=localpc_result$NumTests,
           lpc_time=localpc_result$time_diff) %>%
    mutate(lpc_lmax=localpc_result$lmax) %>% 
    mutate(trial_num=num,
           targets=paste(t,collapse = ",")) %>% 
    mutate(targetSkeletonTimes=paste(localpc_result$targetSkeletonTimes,collapse = ","),
           totalcpptime=localpc_result$totalTime,
           nodes=paste(localpc_result$Nodes,collapse = ",")
    )
  
  return(results)
}

# General Functions -------------------------------------------------------

go_to_dir <- function(file){
  cat("Current dir:\n",list.dirs(),"\n\n")
  cat("Objective dir:\n",list.dirs(),"\n\n")
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

output_text <- function(t,num,algo){
  q1 <- format(paste("Tolerances:",paste0(alpha,",",mb_alpha)),width = 28,justify = "centre")
  q2 <- format(paste("Net:",net),width = 25,justify = "centre")
  q3 <- format(paste("Sample Size:",data.grid$n.obs),width = 18,justify = "centre")
  q4 <- format(paste("MB Algorithm:",algo),width = 18,justify = "centre")
  q5 <- format(paste("Targets:",paste(t,collapse=",")),width = 23,justify = "left")
  q6 <- format(paste("Dataset:",num),width = 15,justify = "centre")
  cat(q1,"|",q2,"|",q3,"|",q4,"|",q5,"|",q6,"\n")
}
