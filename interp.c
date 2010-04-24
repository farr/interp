/*  denest.c: Density estimation from discrete samples. 
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


#include"interp.h"
#include<assert.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>

/* Utility functions and typedefs for manipulating arrays of objects. */

static int
all_equal_pts(size_t n, size_t ndim, double **pts) {
  size_t i;
  double *pt0 = pts[0];

  for (i = 1; i < n; i++) {
    size_t j;

    for (j = 0; j < ndim; j++) {
      if (pt0[j] != pts[i][j]) return 0;
    }
  }

  return 1;
}

void free_interp_tree(tree *t) {
  if (t == NULL) {
    return;
  } else {
    free(t->lower_left);
    free(t->upper_right);
    free_interp_tree(t->left);
    free_interp_tree(t->right);
    free(t);
    return;
  }
}

static size_t
widest_dimension(size_t npts, size_t ndim, double **pts) {
  size_t i, dim = 0;
  double delta = -1.0;

  for (i = 0; i < ndim; i++) {
    size_t j;
    double xmin = pts[0][i], xmax = pts[0][i];
    for (j = 0; j < npts; j++) {
      if (xmin > pts[j][i]) xmin = pts[j][i];
      if (xmax < pts[j][i]) xmax = pts[j][i];
    }

    if (delta < xmax - xmin) {
      delta = xmax - xmin;
      dim = i;
    }
  }
  
  return dim;
}

static void
split_along_dimension(size_t ndim, size_t dim, double x, double *low, double *high,
                      double *new_low, double *new_high) {
  size_t i;
  
  memcpy(new_low, low, ndim*sizeof(double));
  memcpy(new_high, high, ndim*sizeof(double));

  new_low[dim] = x;
  new_high[dim] = x;
}

static int
compare_on_dimension(size_t *dim, const double **pt1, const double **pt2) {
  const double *p1 = *pt1;
  const double *p2 = *pt2;
  size_t comp_dim = *dim;

  if (p1[comp_dim] == p2[comp_dim]) {
    return 0;
  } else if (p1[comp_dim] < p2[comp_dim]) {
    return -1;
  } else {
    return 1;
  }
}

static double 
find_split_point(size_t ndim, size_t npts, double **pts, size_t dim) {
  size_t imid = npts/2;
  size_t step;
  size_t i;

  qsort_r(pts, npts, sizeof(double *), &dim, 
          (int (*)(void *, const void *, const void *))compare_on_dimension);

  for (i = 0; i < npts-1; i++) {
    assert(pts[i][dim] <= pts[i+1][dim]);
  }

  for (step = 1; step <= imid || imid + step < npts; step++) {
    if (step <= imid) {
      if (pts[imid][dim] != pts[imid-step][dim]) {
        return 0.5*(pts[imid][dim] + pts[imid-step][dim]);
      }
    }

    if (imid + step < npts) {
      if (pts[imid][dim] != pts[imid+step][dim]) {
        return 0.5*(pts[imid][dim] + pts[imid+step][dim]);
      }
    }
  }

  assert(0);
}

static size_t
index_of_split(size_t ndim, size_t npts, double **pts, size_t dim, double x) {
  size_t ind = 0;

  while (pts[ind][dim] < x) ind++;
  
  return ind;
}

static tree *
internal_make_interp_tree(size_t ndim, size_t npts, double **pts, double *lower_left,
                           double *upper_right) {
  double *ll, *ur;
  tree *result;

  assert(npts > 0);

  result = malloc(sizeof(tree));
  assert(result != 0);

  ll = malloc(ndim*sizeof(double));
  assert(ll != 0);
  
  ur = malloc(ndim*sizeof(double));
  assert(ur != 0);

  memcpy(ll, lower_left, ndim*sizeof(double));
  memcpy(ur, upper_right, ndim*sizeof(double));

  if (all_equal_pts(npts, ndim, pts)) {
    result->ndim = ndim;
    result->npts = npts;
    result->lower_left = ll;
    result->upper_right = ur;
    result->left = 0;
    result->right = 0;
    return result;
  } else {
    size_t dim = widest_dimension(npts, ndim, pts);
    double x = find_split_point(ndim, npts, pts, dim);
    size_t spind = index_of_split(ndim, npts, pts, dim, x);
    double *new_ll = malloc(ndim*sizeof(double));
    double *new_ur = malloc(ndim*sizeof(double));
    
    assert(new_ll != 0);
    assert(new_ur != 0);

    split_along_dimension(ndim, dim, x, ll, ur, new_ll, new_ur);

    result->ndim = ndim;
    result->npts = npts;
    result->lower_left = ll;
    result->upper_right = ur;
    result->left = make_interp_tree(ndim, spind, pts, ll, new_ur);
    result->right = make_interp_tree(ndim, npts-spind, pts+spind, new_ll, ur);

    free(new_ll);
    free(new_ur);
    
    return result;
  }
}

tree *
make_interp_tree(size_t ndim, size_t npts, double **input_pts, double *ll, double *ur) {
  double **pts = malloc(npts*sizeof(double*));
  tree *result;
  assert(pts != 0);

  memcpy(pts, input_pts, npts*sizeof(double *));

  result = internal_make_interp_tree(ndim, npts, pts, ll, ur);

  free(pts);

  return result;
}

static int
in_bounds(size_t ndim, double *pt, double *low, double *high) {
  size_t i;
  
  for (i = 0; i < ndim; i++) {
    if (pt[i] < low[i] || pt[i] > high[i]) return 0;
  }

  return 1;
}

static tree *
do_find_cell(size_t ndim, double *pt, tree *t) {
  if (!in_bounds(ndim, pt, t->lower_left, t->upper_right)) {
    return NULL;
  } else if (t->left == NULL || t->right == NULL) {
    return t;
  } else { 
    tree *maybe_left = do_find_cell(ndim, pt, t->left);
    if (maybe_left != NULL) {
      return maybe_left;
    } else {
      return do_find_cell(ndim, pt, t->right);
    }
  }
}

static tree *
find_cell(size_t ndim, double *pt, tree *t) {
  assert(in_bounds(ndim, pt, t->lower_left, t->upper_right));
  return do_find_cell(ndim, pt, t);
}

void
sample_density(size_t ndim, size_t npts, double **pts, tree *t,
               uniform_random rng, void *rng_data, double *output_pt) {
  double *pt = pts[(size_t) (npts*rng(rng_data))];
  tree *cell = find_cell(ndim, pt, t);
  size_t i;

  for (i = 0; i < ndim; i++) {
    output_pt[i] = cell->lower_left[i] + rng(rng_data)*(cell->upper_right[i] - cell->lower_left[i]);
  }
}

static double
cell_volume(tree *t) {
  double vol = 1.0;
  size_t i;

  for (i = 0; i < t->ndim; i++) {
    vol *= (t->upper_right[i] - t->lower_left[i]);
  }

  return vol;
}

double 
jump_probability(double *pt, tree *t) {
  size_t npts = t->npts;
  size_t ndim = t->ndim;
  tree *cell = find_cell(ndim, pt, t);
  double vol = cell_volume(cell);

  return ((double)cell->npts)/(vol*npts);
}

double
probability_density(double *pt, tree *t) {
  size_t ndim = t->ndim;

  if (!in_bounds(ndim, pt, t->lower_left, t->upper_right)) {
    return 0.0;
  } else {
    return jump_probability(pt, t);
  }
}

void
bounds_of_points(size_t ndim, size_t npts, double **pts, 
                 double *lower_left, double *upper_right) {
  size_t i;

  memcpy(lower_left, pts[0], ndim*sizeof(double));
  memcpy(upper_right, pts[0], ndim*sizeof(double));

  for (i = 1; i < npts; i++) {
    size_t j;
    for (j = 0; j < ndim; j++) {
      if (lower_left[j] > pts[i][j]) lower_left[j] = pts[i][j];
      if (upper_right[j] < pts[i][j]) upper_right[j] = pts[i][j];
    }
  }
}
