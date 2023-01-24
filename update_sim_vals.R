### Re-create the file with simulation values

if (file.exists("completed_sims.txt")){
  completed_sims <- read.table("completed_sims.txt")
  completed_sims <- unlist(completed_sims)
  names(completed_sims) <- NULL
  completed_sims <- sort(completed_sims)
  remaining_sims <- setdiff(seq(nrow(sim_vals)),completed_sims)
  ind <- which.min(abs(remaining_sims-array_num))
  # array_num <- min(array_num,length(remaining_sims))
  array_num <- remaining_sims[ind]
  cat("Array Number (Updated):",array_num,"\n\n")
} else {
  cat("Array Number:",array_num,"\n\n")
}