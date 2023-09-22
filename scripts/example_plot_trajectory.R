# library(devtools)
# devtools::install_github("revbayes/RevGadgets", force = TRUE, ref = "development")

library(ggplot2)
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

burnin = 0.1
probs = c(0.025, 0.975)
summary = "median"

num_grid_points = 500

max_age = 5e5
min_age = 1e-2

spacing = "exponential"

y_lim = c(1e3, 1e7)
x_lim = c(max_age, 0)


filename = paste0("horses_numiterations", as.integer(NUM_ITERATIONS), "_thinning", THINNING, "_", SAMPLETIMES, "_", MODEL, "_", PIECES, "_int", NUM_INTS, "_", INPUTDATA, "_numtrees", NUM_TREES, "_from", TREESOURCE, "trees_reps", NUM_REPLICATES)


population_size_log = paste0(PATH, filename, "_NEs.log")
interval_change_points_log = paste0(PATH, filename, "_times.log")

df <- processPopSizes(population_size_log, interval_change_points_log, model = PIECES,
                      burnin = burnin, probs = probs, summary = summary,
                      num_grid_points = num_grid_points, spacing = spacing,
                      max_age = max_age, min_age = min_age)
p <- plotPopSizes(df) + ggplot2::coord_cartesian(ylim = y_lim, xlim = x_lim) + ggplot2::ggtitle(NAME) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "italic"), aspect.ratio = 1) +
  ggplot2::theme(panel.background = element_rect(fill='white', color=NA), plot.background = element_rect(fill='transparent', color=NA))

# x axis
p <- p +
  scale_x_reverse(breaks = seq(0, max_age, length.out = 3), labels=seq(0, max_age, length.out = 3)/1000) +
  xlab("Age (k years)") +
  theme(axis.text.x = element_text(face="bold", size=10), axis.title=element_text(size=14))
  
# y axis
p <- p +
  scale_y_log10(breaks = c(1e3,1e5,1e7), labels = c(1, 100, 10000)) +
  ylab("Population size\n(k individuals)") +
  theme(axis.text.y = element_text(face="bold", size=10), axis.title=element_text(size=14))

filename_figure = paste0("popsize_horses_maxage", (max_age/1000), "k_", SAMPLETIMES, "_", NAME, "_", PIECES, "_", INPUTDATA, "_numtrees", NUM_TREES, "_from", TREESOURCE, "trees_reps", NUM_REPLICATES)

ggplot2::ggsave(paste0("figures/trajectories/", filename_figure, ".pdf"), p,
               width = 20, height = 20, units = "cm", bg = "transparent")
