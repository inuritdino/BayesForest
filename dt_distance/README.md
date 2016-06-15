## Distance measure code

Written by Marko Järvenpää.

A quantitative measure of the distribution topology attempts to calculate the complexity
of probability distributions and to compare any given pair of distributions.

The topology is captured by calculating the empirical cumulative distribution functions
along 1D directions, which are placed randomly (uniform by angles or principal component
directions) in the space of the distribution's variates.

The pairwise comparisons are made by calculating the 1D statistic (KS, L^2 etc.) between
the two compared distributions. Then collective behaviour of these 1D statistic measures
along the directions is assessed.

See elsewhere the extensive description of the method.
