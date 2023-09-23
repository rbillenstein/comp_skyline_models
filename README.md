# Comparison of Bayesian Coalescent Skyline Plot  Models for Inferring Demographic Histories

This repository contains the scripts of the manuscript "Comparison of Bayesian Coalescent Skyline Plot Models for Inferring Demographic Histories".

To run the analyses, please install `RevBayes` from the [dev-coalescent branch](https://github.com/revbayes/revbayes/tree/dev-coalescent). You can find instructions on how to do it [here](https://revbayes.github.io/compile-linux) (for the development branch).

A [tutorial](https://revbayes.github.io/tutorials/coalescent/) for running Bayesian Coalescent Skyline Plot Models is provided on the `RevBayes` website.

## Empirical Analyses

The empirical analyses for the manuscript were run on mitochondrial horse data from Vershinina *et al.* (2021).

All analyses can be performed by providing the appropriate arguments to the `mcmc_template.Rev` script in the `scripts` directory:

```bash
rb scripts/mcmc_template.Rev --args DATAFILE AGEFILE SAMPLETIMES INPUTDATA MODEL PIECES NUM_INTERVALS NUM_TREES TREESOURCE NUM_MCMC_ITERATIONS NUM_REPLICATES THINNING
```

- `DATAFILE`: sequence file from the `data` directory

- `AGEFILE`: file providing tip ages from the `data` directory or `NG` in case of isochronous data

- `SAMPLETIMES`: `iso` in the case of isochronous samples, `hetero` in the case of heterochronous samples

- `INPUTDATA`: can be `sequences`, `MAPtree` or `treesample`

- `MODEL`: the Bayesian Coalescent Skyline Plot Model chosen for the analysis (name should match the name in the scripts/models/ directory)

- `PIECES`: the per-interval demographic function - `constant` or `linear`

- `NUM_INTERVALS`: the number of intervals, should be an integer in case of coalescent event independent models (*GMRF*, *HSMRF*, all *Skyfish* models, *Skygrid*) or `NG` in all other cases

- `NUM_TREES`: the number of trees (integer) for the `treesample` analysis, should be `NG` for the two other `INPUTDATA` types

- `TREESOURCE`: the original analysis to draw the `treesample` or `MAPtree` from, can be `constant`or `skyline` for this `mcmc_template`, should be `NG` in case of a sequence based analysis

- `NUM_MCMC_ITERATIONS`: the number of MCMC iterations (integer)

- `NUM_REPLICATES`: the number of replicates (should be 2 or higher for convergence checks)

- `THINNING`: the distance in iterations between samples

Note that all scripts are run from the top directory. 

### All Models, Isochronous Data

In the scripts/models/ directory, Rev scripts for all models described in the manuscript are provided. To run the respective analyses as described in the manuscript for the isochronous dataset, the arguments are:

- `DATAFILE`: `horses_isochronous_sequences.fasta`

- `AGEFILE`: `NG`

- `SAMPLETIMES`: `iso`

- `INPUTDATA`: `sequences`

- `MODEL`: the Bayesian Coalescent Skyline Plot Model chosen for the analysis (name should match the name in the scripts/models/ directory)

- `PIECES`: `constant`

- `NUM_INTERVALS`: `50` or `NG`

- `NUM_TREES`: `NG`

- `TREESOURCE`: `NG`

- `NUM_MCMC_ITERATIONS`: `100000`

- `NUM_REPLICATES`: `2`

- `THINNING`: `10`

This results in the following command:

```bash
rb scripts/mcmc_template.Rev --args horses_isochronous_sequences.fasta NG iso sequences MODEL constant NUM_INTERVALS NG NG 100000 2 10
```

### All Models, Heterochronous Data

For the same analysis on the heterochronous dataset, only the `DATAFILE`, the `AGEFILE`, and the `SAMPLETIMES` need to be changed:

- `DATAFILE`: `horses_heterochronous_sequences.fasta`

- `AGEFILE`: `horses_heterochronous_ages.tsv`

- `SAMPLETIMES`: `hetero`

This results in the following command:

```bash
rb scripts/mcmc_template.Rev --args horses_heterochronous_sequences.fasta horses_heterochronous_ages.tsv hetero sequences MODEL constant NUM_INTERVALS NG NG 100000 2 10
```

### Comparing Joint and Sequential Inference

In the manuscript, we compare demographic inferences from sequence data (jointly inferring genealogies and population size trajectories), from a single *maximum a posteriori* (MAP) tree and from a sample of trees (sequential inference, 10 or 100 trees). We chose the trees to come from the *Constant* analysis with sequence data. All these analyses were run with a *GMRF* model with 50 intervals on the heterochronous dataset.

The following arguments are the same for all four analyses:

- `DATAFILE`: `horses_heterochronous_sequences.fasta`

- `AGEFILE`: `horses_heterochronous_ages.tsv`

- `SAMPLETIMES`: `hetero`

- `MODEL`: `GMRF`

- `PIECES`: `constant`

- `NUM_INTERVALS`: `50`

- `NUM_MCMC_ITERATIONS`: `100000`

- `NUM_REPLICATES`: `2`

- `THINNING`: `10`

The following arguments should be adjusted to the specific analysis:

- `INPUTDATA`: `sequences` or `MAPtree`  or `treesample`

- `NUM_TREES`: `NG`  (in case of `sequences` and `MAPtree`), `10`, or `100`

- `TREESOURCE`: `NG` (in case of `sequences`), `constant` (in the other cases)

### Comparing Constant and Linear Per-Interval Population Sizes

To compare constant and linear per-interval population sizes, we ran the *GMRF* analysis with 50 intervals on heterochronous sequence data for both settings. The only argument that differs is the `PIECES` argument:

- `DATAFILE`: `horses_heterochronous_sequences.fasta`

- `AGEFILE`:` horses_heterochronous_ages.tsv`

- `SAMPLETIMES`:` hetero`

- `INPUTDATA`: `sequences`

- `MODEL`:` GMRF`

- `PIECES`: `constant` or `linear`

- `NUM_INTERVALS`: `50`

- `NUM_TREES`: `NG`

- `TREESOURCE`: `NG`

- `NUM_MCMC_ITERATIONS`: `100000`

- `NUM_REPLICATES`: `2`

- `THINNING`: `10`

### Comparing Autocorrelated and Uncorrelated *Skyfish* Models

We tested three different implementations of our new *Skyfish* model on the heterochronous sequence data. An autocorrelated version (scripts/models/SkyfishAC.Rev), an uncorrelated version (scripts/models/SkyfishIID.Rev), and an uncorrelated verion with an empirically informed prior (scripts/models/SkyfishEP.Rev). For these analyses, only the `MODEL` argument differs:

- `DATAFILE`: `horses_heterochronous_sequences.fasta`

- `AGEFILE`: `horses_heterochronous_ages.tsv`

- `SAMPLETIMES`: `hetero`

- `INPUTDATA`: `sequences`

- `MODEL`: `SkyfishAC` or `SkyfishIID` or `SkyfishEP`

- `PIECES`: `constant` 

- `NUM_INTERVALS`: `NG`

- `NUM_TREES`: `NG`

- `TREESOURCE`: `NG`

- `NUM_MCMC_ITERATIONS`: `100000`

- `NUM_REPLICATES`: `2`

- `THINNING`: `10`

## Simulation Study

For our simulation study, we used the resulting median population size trajectories from the analysis with the *BSP* model and the *GMRF* model on isochronous sequence data to simulate 10 datasets of sequences. We then performed analyses of the 10 simulated datasets with both the *BSP* and the *GMRF* model.

First, for extracting the population size trajectories from the two analyses, run the `evaluation` script in the scripts directory from the top directory. In the evaluation folder, we also provide the extracted trajectories used for the manuscript.

Next, run the simulations by running the `simulate_fromBSP.Rev` and the `simulate_fromGMRF.Rev` script that are in the scripts/simulations/ directory from the top directory.

Last, run the script `mcmc_template_simstudy` that is also in the scripts/simulations/ directory from the top directory. It takes the following arguments:

- `DATAMODEL`: the model which was used for the simulation - `BSP` or `GMRF`

- `INDEX`: the index of the simulation, ranging from 1 to 10

- `MODEL`: the model that should be used for the analysis of the simulated data - we chose `BSP` or `GMRF`

- `PIECES`: the per-interval demographic function - `constant`or `linear`, we chose `constant`

- `NUM_INTS`: the number of intervals - for our analyses `NG` in the case of the *BSP* model,  `50` in the case of the *GMRF* model

- `N_TIPS `= the number of tips to be simulated, we chose `36` as in the original isochronous dataset

## Validation

For validation of our implementation, we compared results of an analysis with the *BSP* model to results generated with [BEAST v1.10.4](https://beast.community/). We adjusted our MCMC scripts to match the `BEAST` scripts as much as possible. Both software programs were run on the first four simulations generated under the resulting population size trajectory from the analysis of isochronous data with the *GMRF* model. We provide these four simulated datasets in output/simulations. The scripts for running both the `RevBayes` and the `BEAST` analyses are provided in the scipts/validation directory. The `RevBayes`scripts should be run from the top directory, but the `BEAST` scripts should be run from the scripts/validation directory. Once done, they can be converted to `RevBayes`-like output files running the `conversion_BEAST_toRB.R` script that is in the scripts/validation/ directory from the top directory. Afterwards, all (converted) ouput files will be in output/validation/ and can be compared with each other by plotting with `RevGadgets` (see next two sections).

### Plotting Results

For plotting population size trajectories, we used the `R` package `RevGadgets`, installed from the [development branch](https://github.com/revbayes/RevGadgets/tree/development) using `devtools`:

```R
devtools::install_github("revbayes/RevGadgets", ref = "development)
```

`RevGadgets` includes two functions to work with output from the Bayesian coalescent skyline plot models in `RevBayes`:

- `processPopSizes` is for post-processing the output, evaluating the population sizes and interval change-points from the `RevBayes` output on a user-defined grid. It can either return median population sizes and their credible intervals (CIs) on that grid or the whole posterior distributions at the grid points. The latter is useful for convergence assessment (see next section).

- `plotPopSizes` uses the dataframe with median population sizes and CIs generated with `processPopSizes` for plotting the population size trajectories.

For an example of how to use these functions, see scripts/example_plot_trajectory.R. You can also find a detailed description of the different options for these functions in the associated documentation in `R` (e. g. `?processPopSizes`) or in the [post processing section](https://revbayes.github.io/tutorials/coalescent/postprocessing) of the tutorial.

## Convergence Assessment

The convergence assessment of the population size trajectories in the manuscript was done using a Kolmogorov-Smirnov (ks) test. The threshold for the test statistic *D* was set to $0.0921$ (Fabreti and Höhna, 2022). We used `processPopSizes` from the `R` package `RevGadgets` (see previous section) to get posterior distributions of the population sizes from two replicates of the same analysis. These were evaluated at user-defined grid points, in our case usually 500 grid points, and then plotted. For an example script see scripts/example_plot_convergence.R. You can also find a description in the [post processing section](https://revbayes.github.io/tutorials/coalescent/postprocessing) of the tutorial.

## References

Fabreti and Höhna 2022. Convergence assessment for Bayesian phylogenetic analysis using MCMC simulation. Methods in Ecology and Evolution, 13(1): 77–90

Vershinina *et al.* 2021. Ancient horse genomes reveal the timing and extent of dispersals across the bering land bridge. Molecular Ecology, 30(23): 6144–6161

*BSP*: Drummond *et al.* 2005. Bayesian Coalescent Inference of Past Population Dynamics from Molecular Sequences. Molecular Biology and Evolution, 22(5): 1185–1192

*EBSP*: Heled and Drummond 2008. Bayesian inference of population size history from multiple loci. BMC Evolutionary Biology, 8(1): 289

*GMRF* and *HSMRF*: Faulkner *et al.* 2020. Horseshoe-based Bayesian nonparametric estimation of effective population size trajectories. Biometrics, 76(3): 677–690

*Skygrid*: Gill *et al.* 2012. Improving Bayesian Population Dynamics Inference: A Coalescent-Based Model for Multiple Loci. Molecular Biology and Evolution, 30(3): 713–724

*Skyride*: Minin *et al.* Smooth Skyride through a Rough Skyline: Bayesian Coalescent-Based Inference of Population Dynamics. Molecular Biology and Evolution, 25(7): 1459–1471
