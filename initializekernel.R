#########################################################################
# This is the kernel for initialize{test|hoffman}.R
#########################################################################

set.seed(555)

# Download the newest version of the Local FCI package
if (load_package){
  devtools::install_github("stephenvsmith/LocalFCI",upgrade="never") # ,quiet=TRUE  
}


# Set simulation settings
alpha <- c(0.01,0.05,0.1)
mb_alpha <- c(0.01,0.05,0.1)
net_names <- gsub(".rds","",list.files(rds_dir))
net <- c("alarm","barley","insurance","arth150","andes","diabetes",
         "link","mildew","hepar2","pigs","andes","munin2","munin3","munin4")
high <- c(0.25) # For variances
ub <- c(0.75) # For coefficients
n <- c(500,1000,10000)
algos <- c("MMPC","SES")

# Create a table with every combination of simulation settings
simulation_values <- expand.grid(alpha=alpha,mb_alpha=mb_alpha,net=net,
                                 high=high,ub=ub,n=n,algos=algos)

go_to_dir(result_dir)
if (file.exists("incomplete.txt")){
  cat("We are skipping making the sim_vals file")
} else {
  write.csv(simulation_values,file = "sim_vals.csv")
}
