##########################################################################
# This file provides the kernel for arrayscript{test|hoffman}.R
##########################################################################

### Simulation Setup
# Load the file with functions for generating data
source(data_gen_file)
# Define the number of trials for each simulation setting
num_trials <- 3
# Define the number of targets considered for each setting
max_targets <- 4
# Define the number of cores used for parallel processing
num_cores <- min(parallel::detectCores()-2,4)
num_cores <- 1
cat("We are using",num_cores,"core(s).\n\n")

# Set up a generic data grid for simulated data
# True values will be replace these later
data.grid <- data.frame(network = "asia",
                        data.type = "continuous",
                        n.dat = num_trials,
                        n.obs = 1000,
                        c.ratio = 0,
                        max.in.degree = Inf,
                        lb = 0.4,  # lower bound of coefficients
                        ub = 5,  # upper bound of coefficients
                        low = 0.01,  # lower bound of variances if continuous, of number of levels if discrete
                        high = 5,  # upper bound of variances if continuous, of number of levels if discrete
                        scale = FALSE,
                        stringsAsFactors = FALSE)


go_to_dir(result_dir)
### Read the table containing each simulation setting
sim_vals <- read.csv("sim_vals.csv",stringsAsFactors = FALSE)[,-1]

### TODO: MODIFY THIS SO THAT WE AREN'T WORKING ON SAME SIMULATION TWICE
# if (file.exists("completed_sims.txt")){
#   completed_sims <- read.table("completed_sims.txt")
#   completed_sims <- unlist(completed_sims)
#   names(completed_sims) <- NULL
#   completed_sims <- sort(completed_sims)
#   remaining_sims <- setdiff(seq(nrow(sim_vals)),completed_sims)
#   ind <- which.min(abs(remaining_sims-array_num))
#   # array_num <- min(array_num,length(remaining_sims))
#   array_num <- remaining_sims[ind]
#   cat("Array Number (Updated):",array_num,"\n\n")
# } else {
#   cat("Array Number:",array_num,"\n\n")
# }

### Save the settings for this simulation run
alpha <- sim_vals$alpha[array_num]
mb_alpha <- sim_vals$mb_alpha[array_num]
net <- sim_vals$net[array_num]
high <- sim_vals$high[array_num]
ub <- sim_vals$ub[array_num]
n <- sim_vals$n[array_num]
algo <- sim_vals$algo[array_num]

cat("This is row",array_num,"of our simulation grid.\n")
cat("Simulation Parameters:\n")
cat("Significance:",alpha,"\n")
cat("Markov Blanket Estimation Significance:",mb_alpha,"\n")
cat("Network:",net,"\n")
cat("high:",high,"\n")
cat("ub:",ub,"\n")
cat("n:",n,"\n")
cat("MB algo:",algo,"\n")

# Obtain network information, including the true DAG adj. mat.
network_info <- get_network_DAG(net)

# Generate/Retrieve Targets
targets <- check_targets_defined_get_targets(network_info)

# Set up for simulated data and directory for context
data.grid$network <- net
data.grid$n.obs <- n
n <- n
ub <- ub
high <- high
data.grid$ub <- ub
data.grid$high <- high

# Generate the data
simulation_data_creation()

# Grab simulated data
df_list <- lapply(1:num_trials,function(i) grab_data(i))

# Keep network directory
curr_dir <- getwd()

# Get results for each trial if they exist
results_list <- lapply(1:num_trials,function(num){
  if (file.exists(paste0("results_",array_num,"_",num,".rds"))){
    return(readRDS(paste0("results_",array_num,"_",num,".rds")))
  }
  
  # Run Global PC Algorithm
  trial_num <- num
  cat("Running Global PC for Dataset",num,"... ")  
  results_pc <- run_global_pc(df_list[[num]],trial_num)
  time_pc <- results_pc$time_diff$PC
  units(time_pc) <- "secs"
  cat("completed in",as.numeric(time_pc),units(time_pc),"\n")

  # Run Local FCI Algorithm for all sets of targets
  results_lfci_df <- mclapply(targets,
                              run_fci_target,
                              df=df_list[[num]],
                              num,
                              results_pc,
                              algo,
                              curr_dir,
                              mc.preschedule = FALSE,mc.cores = num_cores)
  
  # Run Local FCI Algorithm for all sets of targets
  results_lpc_df <- mclapply(targets,
                             run_pc_target,
                             df=df_list[[num]],
                             num,
                             results_pc,
                             algo,
                             curr_dir,
                             mc.preschedule = FALSE,mc.cores = num_cores)

  results_final_df_lfci <- data.frame(do.call(rbind,results_lfci_df))
  results_final_df_lpc <- data.frame(do.call(rbind,results_lpc_df))
  results_final_df <- results_final_df_lfci %>% 
    left_join(results_final_df_lpc,by=c("targets","trial_num"),suffix=c("_lfci","_lpc")) %>% 
    rename_with(~paste0(str_extract(.,"(lpc|lfci)$"),
                        "_",
                        str_extract(.,"^[a-z]+")),
                ends_with(c("_lpc","_lfci")))
  
  saveRDS(results_final_df,paste0("results_",array_num,"_",num,".rds"))
  setwd(curr_dir)
  return(results_final_df)
})

results_iteration <- data.frame(do.call(rbind,results_list))
saveRDS(results_iteration,paste0("results_",array_num,".rds"))
# Remove temporary files (if the program gets to this point)
for (i in 1:num_trials){
  unlink(paste0("results_",array_num,"_",i,".rds"))
}
unlink(paste0("pc_",array_num,"_results.rds"))
cat(array_num,file = "completed_sims.txt",append = TRUE)
cat("\n",file = "completed_sims.txt",append = TRUE)


