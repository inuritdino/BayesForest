# Bayes Forest Toolbox
#### Realistic data-based stochastic modelling of tree morphologies

We offer an approach to model structure/morphology of trees by iteratively optimizing "distance" between empirical distributions describing structure of the data tree (as obtained from the laser scanning and represented as a Quantitative Structure Model, QSM) and those describing structure of a Stochastic Structure tree Model (SSM). The resulting SSM has 'optimal' parameter values to produce trees *statistically* similar to the data tree and to each other, that is the SSM sample trees are not exact copies of each other. Thus, unlimited number of sample trees can be generated.

This algorithm can be viewed as a "clone" trees generator, given that the SSM optimized to the data possesses stochastic rules of growth
and is statistically equivalent (as much as the assumptions behind the model, its parameters, and stochastic influences are plausible) 
to the data tree (which is, in fact, again a model, i.e. a reconstruction from the laser scanning point cloud data). The resulting best-fit SSM is a generator of the morphological clones (given its own stochasticity is on). Morphological clones
possess similar overall structure (coarse-grain scale), while vary in minute
details of organization (fine-grain scale).

Our approach is based upon five distinct parts:

1. Quantitative Structural Model (QSM) that is to be obtained from the Terrestrial Laser Scanning (TLS) data (see [Raumonen et al., *Remote Sensing* 2013]).
2. Stochastic Structural Model (SSM), that is an analytical tree growth model. For example, one of Functional-Structural Plant Models (FSPM) or structural models with heuristic rules for growth (SHM or procedural models).	
3. Structural data sets, that is data sets U relating different physical dimensions as well as spatial location of various parts and segments of a tree with optional sorting by the topological characteristics.
4. Distance, that is a measure of proximity between any two data sets, in other words, a value quantifying how similar the two data sets are. 
5. Optimization algorithm, that is an iterative procedure capable of finding a minimum of any given function (Newton algorithm, genetic algorithm etc.).

This approach is described in further details in the following publications:
a. "Data-based stochastic modeling of tree growth and structure formation" by
I. Potapov et al., Silva Fennica 50(1):1413.

==============================

The main function is BayesForest.m, which takes an input configuration file describing
all the details of the whole procedure. Refer to the help documentation of a 
particular function for further details (configuration file syntax is described
in *bf_process_input()* function's documentation).



