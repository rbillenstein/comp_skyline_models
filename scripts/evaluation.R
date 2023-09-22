# library(devtools)
# devtools::install_github("revbayes/RevGadgets", force = TRUE, ref = "development")

library(RevGadgets)

SAMPLETIMES = "iso"
INPUTDATA = "sequences"
PIECES = "constant"
NUM_TREES = "NG"
TREESOURCE = "NG"

NUM_ITERATIONS = as.numeric(100000)
NUM_REPLICATES = 2
THINNING = 10

NUM_INTERVALS = 50

DATAMODEL = c("BSP", "GMRF")

PATH = paste0("output/", as.integer(NUM_ITERATIONS), "iterations/", INPUTDATA, "/")

burnin = 0.1
probs = c(0.025, 0.975)
summary = "median"

num_grid_points = 500

max_age = 5e5

spacing = "equal"
min_age = 0

######################

for (i in 1:length(DATAMODEL)){
  MODEL = DATAMODEL[i]
  
  if (MODEL == "GMRF") NUM_INTS = NUM_INTERVALS else NUM_INTS = "NG"
  
  filename = paste0("horses_numiterations", as.integer(NUM_ITERATIONS), "_thinning", THINNING, "_", SAMPLETIMES, "_", MODEL, "_", PIECES, "_int", NUM_INTS, "_", INPUTDATA, "_numtrees", NUM_TREES, "_from", TREESOURCE, "trees_reps", NUM_REPLICATES)
  
  population_size_log = paste0(PATH, filename, "_NEs.log")
  interval_change_points_log = paste0(PATH, filename, "_times.log")
  
  df <- processPopSizes(population_size_log, interval_change_points_log, model = PIECES,
                        burnin = burnin, probs = probs, summary = summary,
                        num_grid_points = num_grid_points, spacing = spacing,
                        max_age = max_age, min_age = min_age)
  
  write.table(df, file = paste0("evaluation/horses_iso_", MODEL,".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}
