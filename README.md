# Bayes Forest Toolbox
#### Realistic data-based stochastic modelling of tree morphologies

We offer an approach to model structure/morphology of trees by iteratively optimizing "distance" between empirical distributions describing structure of the data tree (as obtained from the laser scanning) and those describing structure of a Stochastic Structural tree Model (SSM). The resulting SSM has 'optimal' parameter values to produce trees *statistically* similar to the data tree and to each other, that is the SSM sample trees are not exact copies of each other. Thus, unlimited number of sample trees can be generated.

This algorithm can be viewed as a "clone" trees generator, given that the SSM optimized to the data possesses stochastic rules of growth
and is statistically equivalent (as much as the assumptions behind the model, its parameters, and stochastic influences are plausible) 
to the data tree (which is, in fact, again a model, i.e. a reconstruction from the laser scanning point cloud data).

Our approach is based upon five distinct parts:

1.	Quantitative Structural Model (QSM) that is to be obtained from the Terrestrial Laser Scanning (TLS) data (see [Raumonen et al., 2013]).
2. Stochastic Structural Model (SSM), that is an analytical tree growth model. For example, one of Functional-Structural Plant Models (FSPM) or structural models with heuristic rules for growth (SHM or procedural models).	
3. Structural data sets, that is data sets U relating different physical dimensions as well as spatial location of various parts and segments of a tree with optional sorting by the topological characteristics.
4. Distance, that is a measure of proximity between any two data sets, in other words, a value quantifying how similar the two data sets are. 
5. Optimization algorithm, that is an iterative procedure capable of finding a minimum of any given function (Newton algorithm, genetic algorithm etc.).


