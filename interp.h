/*  denest.h: Density estimation from discrete samples.  
    Copyright (C) 2009-2010 Will M. Farr <w-farr@northwestern.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

/* interp: Density estimation from a discrete set of samples.  Written
   Dec 2009 by Will M. Farr <w-farr@northwestern.edu>.  

   Modified:
   April 2010 by Will M. Farr <w-farr@northwestern.edu>.

   Given a set of points in arbitrary dimension, this code subdivides
   the domain of the points recursively into cells such that each cell
   contains exactly one unique point (there may be many copies of the
   same point, such as the repeated samples that might come out of an
   MCMC).  It does this by cutting the domain along the dimension with
   the largest spread of points such that (approximately) half the
   points lie in each of the two resulting cells.

   The code assumes that the domain is topologically R^n; if there
   are, e.g., angular coordinates that wrap then there will be
   inaccuracies near the edges of the domain (0 and 2*pi) where points
   that are actually close are regarded by the code to be
   well-separated.  Beware!

   The code does not try to scale any of the dimensions of the domain;
   it is to your advantage to have all dimensions of approximately the
   same magnitude, since the code will choose the numerically largest
   extent of points to partition at each subdivision.

   To use this code to interpolate a PDF from a set of discrete
   samples, provide your samples to the make_density_tree function.
   Then use the sample_density function to draw output points
   distributed according to a (constant-in-cell) interpolation of your
   PDF.  This is done by first choosing (uniformly) a random sample
   point, then finding the cell that contains only this sample point,
   and finally drawing an output point uniformly in this cell.  This
   corresponds to drawing from a piecewise-constant-in-cell
   interpolation of the PDF represented by the sample points.

   The jump_probability function is useful when using sample_density
   in the context of an MCMC, where one needs to both propose random
   jumps from the current point and compute the probability of
   choosing the given jump.  Given a point (assumed to be drawn using
   the sample_density function), jump_probability returns the
   probability density that sample_density would return that point
   from among all possible choices.  (NOTE: jump_probability returns
   the actual jump probability, NOT the log of the jump probability.
   Logs are often used in MCMC simulations to prevent under- and
   over-flow of probabilities.)

   The code makes use of the qsort_r function (which is a ?BSD?
   extension to the C standard library qsort function).  It is
   thread-safe, re-entrant, and all that jazz.  The time complexity of
   constructing the tree is O(N*log(N)^2) (contrast with the minimum
   time for construction of a kD tree, O(N*log(N))---the extra log(N)
   factor comes from sorting the points to determine the split along
   each dimension rather than using a faster, O(N), median-finding
   algorithm).  The time complexity for proposing a jump and computing
   the jump probability is O(log(N)).  

*/

#ifndef __INTERP_H__
#define __INTERP_H__

#include<string.h>

/* The tree datastructure we use for computing the density. */
struct cell_struct {
  size_t ndim;
  size_t npts; /* The number of repeated points in this cell. */
  double *lower_left; /* ndim coordinates of lower left of box. */
  double *upper_right; /* ndim coordinates of upper right corner. */
  struct cell_struct *left; /* the smaller points in the coordinate number split_dim. */
  struct cell_struct *right; /* the larger points in the coordinate number split_dim. */
  /* left AND right are NULL when this cell contains only one (unique)
     point from the initial distribution. */
};

typedef struct cell_struct tree;

/* Typedef for the random number generator function used to sample
   from the density tree in sample_density.  The function should take
   a pointer to some data (e.g. RNG state), and return a uniform
   random double in [0,1]. */
typedef double (*uniform_random)(void *data);

/* Construct a density tree from the npts points stored in pts; each
   point is an array of ndim double precision numbers.  The total
   boundary of the domain is a box in R^n bounded below by lower_left,
   and above by upper_right. 

   The procedure will make its own copies of the lower_left and
   upper_right arrays, so these may be freed or altered upon return
   arbitrarily.  Similarly with the pts array. */
tree *make_interp_tree(size_t ndim, size_t npts, double **pts, double *lower_left, 
                        double *upper_right);

/* Free the memory associated with the tree t. */
void free_interp_tree(tree *t);

/* Draw a sample from the piecewise-constant interpolation of a
   distribution sampled by the npts points (each of dimension ndim) in
   the array pts using the tree data structure in tree.  rng is used
   to draw random numbers for the sampling.  The pointer output_pt
   should be to a memory area with space for ndim doubles; on output,
   it will contain the coordinates of the sample point. */
void sample_density(size_t ndim, size_t npts, double **pts, tree *tree, 
                    uniform_random rng, void *rng_data, double *output_pt);

/* Returns the probability density that a call to sample_density with
   npts number of points will return the point pt.  This probability
   is (part) of the computation of the jump_probability in the
   Metropolis-Hastings sampling algorithm for an MCMC.  */
double jump_probability(double *pt, tree *tree);

/* Returns the value of the piecewise-constant PDF in tree at the
   point pt. */
double probability_density(double *pt, tree *tree);

#endif /* __INTERP_H__ */
