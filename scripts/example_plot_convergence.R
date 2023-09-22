# library(devtools)
# devtools::install_github("revbayes/RevGadgets", force = TRUE, ref = "development")

library(RevGadgets)

SAMPLETIMES = "iso"
INPUTDATA = "sequences"
PIECES = "constant"
NUM_TREES = "NG"
TREESOURCE = "NG"

NUM_ITERATIONS = as.numeric(100000)
# RevBayes changes format of big numbers starting at 1e7 (1e+07),
# below that it is e.g. 1e4 -> 10000
NUM_REPLICATES = 2
THINNING = 10

NUM_INTS = "NG"

PATH = paste0("output/", as.integer(NUM_ITERATIONS), "iterations/", INPUTDATA, "/")

MODEL = "SkyfishAC"
NAME = "Skyfish"

pairs_to_check = combn(1:NUM_REPLICATES,2)
names_runs = c()
for (i in 1:dim(pairs_to_check)[2]){
  names_runs = c(names_runs, paste0("run ", pairs_to_check[1, i], " vs run ", pairs_to_check[2,i]))
}

burnin = 0.1
num_grid_points = 500

max_age = 5e5
min_age = 1e-2
model = "constant"

ymax = 0.1


filename = paste0("horses_numiterations", as.integer(NUM_ITERATIONS), "_thinning", THINNING, "_", SAMPLETIMES, "_", MODEL, "_", PIECES, "_int", NUM_INTS, "_", INPUTDATA, "_numtrees", NUM_TREES, "_from", TREESOURCE, "trees_reps", NUM_REPLICATES)

pop_size_dist = list()
for (r in 1:NUM_REPLICATES){
  population_size_log = paste0(PATH, "/", filename, "_NEs_run_", r, ".log")
  interval_change_points_log = paste0(PATH, "/", filename, "_times_run_", r, ".log")

  pop_size_dist[[r]] = processPopSizes(population_size_log, interval_change_points_log, model = model,
                                       burnin = burnin, num_grid_points = num_grid_points, spacing = "exponential",
                                       max_age = max_age, min_age = min_age, distribution = TRUE)
  
  cat(paste0("done with distribution of rep ", r, "\n"))
  
}

ks_pairs = list()
for (combs in 1:dim(pairs_to_check)[2]){
  first = pairs_to_check[1, combs]
  second = pairs_to_check[2, combs]
  ks_test_stat = list()
  for (i in 1:num_grid_points){
    ks_test_stat[[i]] = ks.test(pop_size_dist[[first]][i,], pop_size_dist[[second]][i,])$statistic
  }
  
  ks_data = as.data.frame(cbind("timepoints" = as.numeric(rownames(pop_size_dist[[1]])), "D" = unlist(ks_test_stat, use.names = FALSE)))
  ks_pairs[[combs]] = ks_data
  
  ymax = max(ymax, max(unlist(ks_test_stat)))
}

ks_p = ggplot(data = ks_pairs[[1]], mapping = aes(x = timepoints, y = D)) +
        geom_point(color = rainbow(dim(pairs_to_check)[2])[1], size = 1, shape = 19)

if (dim(pairs_to_check)[2] > 1){
  for (combs in 2:dim(pairs_to_check)[2]){
    ks_p = ks_p +
      geom_point(data = ks_pairs[[combs]], color = rainbow(dim(pairs_to_check)[2])[combs], size = 1, shape = 19)
  }
}

ks_p = ks_p +
        geom_hline(yintercept = 0.0921, col = "blue", linewidth = 1) +
        scale_x_reverse() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank()) +
        ggtitle(NAME) +
        theme(plot.title = element_text(hjust = 0.5, face = "italic"), aspect.ratio = 1) +
        theme(panel.background = element_rect(fill='white', color=NA), plot.background = element_rect(fill='transparent', color=NA))

# x axis
ks_p <- ks_p +
  scale_x_reverse(breaks = seq(0, max_age, length.out = 3), labels=seq(0, max_age, length.out = 3)/1000) +
  xlab("Age (k years)") +
  theme(axis.text.x = element_text(face="bold", size=10), axis.title=element_text(size=14))


# y axis
ks_p <- ks_p +
  theme(axis.text.y = element_text(face="bold", size=10), axis.title=element_text(size=14))


filename_figure = paste0("horses_", SAMPLETIMES, "_", NAME, "_", PIECES, "_", INPUTDATA, "_numtrees", NUM_TREES, "_from", TREESOURCE, "trees_reps", NUM_REPLICATES)

ggplot2::ggsave(paste0("figures/convergence/convergence_", filename_figure,".pdf"), ks_p,
                width = 20, height = 25, units = "cm", bg = "transparent")
