# BEAST to RevGadgets Conversion
# BSP: we need skyline.popSize1 etc, skyline.groupSize1 etc, and node ages from the trees

library(ape)

BEAST_results_dir = "scripts/validation/"
validation_dir = "output/validation/"

for (s in 1:4){
  filename = paste0("simulation_fromGMRF_ntips_36_sequences_", s, "_BEAST_BSP")
  
  logfile_path = paste0(BEAST_results_dir, filename, ".log")
  
  logfile = read.table(logfile_path, header = TRUE)
  
  popSizes = logfile[, grep(pattern="popSize", colnames(logfile))]
  write.table(popSizes, paste0(validation_dir, filename, "_RevVersion_NEs.log"), quote = FALSE, row.names = FALSE, sep = "\t")
  
  treefile_path = paste0(validation_dir, filename, ".trees")
  trees = read.nexus(treefile_path)
  
  branching_points = lapply(trees, function(x) sort(branching.times(x)))
  
  groupSizes = logfile[, grep(pattern="groupSize", colnames(logfile))]
  
  intervalIndices = t(apply(groupSizes, 1, function(x) cumsum(as.numeric(x))))
  
  times = list()
  for (i in 1:length(branching_points)){
    times[[i]] = branching_points[[i]][intervalIndices[i,]]
  }
  
  times_df = t(data.frame(times))
  colnames_times = c()
  for (n in 1:dim(times_df)[2]){
    colnames_times = c(colnames_times, paste0("changepoint_", n))
  }
  
  colnames(times_df) = colnames_times
  rownames(times_df) = NULL
  write.table(times_df, paste0(validation_dir, filename, "_RevVersion_times.log"), quote = FALSE, row.names = FALSE, sep = "\t")
  
}

