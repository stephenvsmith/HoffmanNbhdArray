### Directories (Testing)
#home_dir <- '/home/stephen'
home_dir <- '/Users/stephensmith'
result_dir <- paste0(home_dir,'/Desktop/','ResultsFullSample-',format(Sys.Date(),"%m-%y"))
rds_dir <- paste0(home_dir,'/Dropbox/Academics/Research/Code/Networks/rds')
data_gen_file <- paste0(home_dir,'/Dropbox/Academics/Research/Code/Scripts/data_gen.R')

### Setup (Testing)
source(paste0(home_dir,'/Dropbox/Academics/Research/Code/Scripts/HoffmanNbhdArray/helperfunctions.R'))
sim_vals <- read.csv(paste0(home_dir,'/Desktop/ResultsFullSample-',format(Sys.Date(),"%m-%y"),'/sim_vals.csv'))
array_num <- sample(1:nrow(sim_vals),1)
array_num <- 357
#cat("Simulation",array_num,file = "~/Desktop/lastSim.txt")

source(paste0(home_dir,'/Dropbox/Academics/Research/Code/Scripts/HoffmanNbhdArray/arraykernel.R'))
