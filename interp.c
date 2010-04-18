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

/* Utility functions and typedefs for manipulating arrays of objects. */

/* Compare functions should return -1 if obj1 < obj2, 0 if obj1 =
   obj2, and 1 if obj1 > obj2. */
typedef int (*compare)(void *obj1, void *obj2, void *data);

/* Predicates return zero if obj does not satisfy the predicate, and 1
   if not. */
typedef int (*predicate)(void *obj, void *data);

static void
swap(void **objs, size_t i, size_t j) {
  void *tmp = objs[i];

  objs[i] = objs[j];
  objs[j] = tmp;
}

/* Returns nsat, the number of objets that satisfy pred in the array
   of objects.  These objects are moved to the front of the array, so
   objs[0] through objs[nsat-1] satisfy the predicate; objs[nsat]
   through objs[n-1] do not satisfy the predicate. */
static size_t
partition(size_t n, void **objs, predicate pred, void *pred_data) {
  size_t lowest_false = 0, highest_true = n-1;

  while (1) {
    while (pred(objs[lowest_false], pred_data)) {
      lowest_false++;
    }
    
    while (!(pred(objs[highest_true], pred_data))) {
      highest_true--;
    }

    if (lowest_false < highest_true) {
      swap(objs, lowest_false, highest_true);
    } else {
      break;
    }
  }

  return lowest_false;  
}

typedef struct {
  compare comp;
  void *comp_data;
  void *obj;
} less_than_given_object_data;

static int
less_than_given_object(void *obj, less_than_given_object_data *data) {
  int c = data->comp(obj, data->obj, data->comp_data);

  return (c < 0);
}

static int
all_equal(size_t n, void **objs, compare comp, void *comp_data) {
  size_t i;
  void *obj0 = objs[0];

  for (i = 1; i < n; i++) {
    if (comp(obj0, objs[i], comp_data) != 0) {
      return 0;
    }
  }

  return 1;
}

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

/* Returns the nth obj in order according to the comparison function
   comp.  The array of objects is re-ordered during this operation
   (not sorted---the computation cost is O(n), not O(n*log(n))). */
static void *
find_nth(size_t n, size_t nth, void **objs, compare comp, void *comp_data) {
  assert(nth < n);
  if (n == 1 && nth == 0) {
    return objs[0];
  } else if (all_equal(n, objs, comp, comp_data)) {
    return objs[0];
  } else {
    void *partition_obj = objs[rand() % n]; /* Probably *very* not random, but good enough. */
    less_than_given_object_data ltdata;
    size_t nlt;

    ltdata.comp = comp;
    ltdata.comp_data = comp_data;
    ltdata.obj = partition_obj;

    nlt = partition(n, objs, (predicate) less_than_given_object, &ltdata);

    if (nth < nlt) {
      return find_nth(nlt, nth, objs, comp, comp_data);
    } else {
      return find_nth(n-nlt, nth-nlt, objs, comp, comp_data);
    }
  }
}

void free_density_tree(tree *t) {
  if (t == NULL) {
    return;
  } else {
    free(t->lower_left);
    free(t->upper_right);
    free_density_tree(t->left);
    free_density_tree(t->right);
    free(t);
    return;
  }
}

typedef struct {
  size_t dim;
  double x0;
} one_d_data;

static int 
lt_one_d(double *pt, one_d_data *data) {
  return pt[data->dim] < data->x0;
}

static int 
gt_one_d(double *pt, one_d_data *data) {
  return pt[data->dim] > data->x0;
}

static int
compare_along_dim(double *pt1, double *pt2, size_t *d) {
  size_t dim = *d;

  if (pt1[dim] < pt2[dim]) {
    return -1;
  } else if (pt1[dim] > pt2[dim]) {
    return 1;
  } else {
    return 0;
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

static size_t
partition_pts(size_t n, size_t ndim, size_t dim, double **pts) {
  double *median_pt = find_nth(n, n/2, (void **)pts, (compare) compare_along_dim, &dim);
  one_d_data data;
  size_t nlt;

  data.dim = dim;
  data.x0 = median_pt[dim];

  nlt = partition(n, (void **) pts, (predicate) lt_one_d, &data);
  if (nlt == 0) {
    nlt = partition(n, (void **) pts, (predicate) gt_one_d, &data);
    nlt = n-nlt;
  }

  return nlt;
}

static tree *
internal_make_density_tree(size_t ndim, size_t npts, double **pts, double *lower_left,
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
    size_t nlt = partition_pts(npts, ndim, dim, pts);
    double x = 0.5*(pts[nlt-1][dim] + pts[nlt][dim]);
    double *new_ll = malloc(ndim*sizeof(double));
    double *new_ur = malloc(ndim*sizeof(double));
    
    assert(new_ll != 0);
    assert(new_ur != 0);

    split_along_dimension(ndim, dim, x, ll, ur, new_ll, new_ur);

    result->ndim = ndim;
    result->npts = npts;
    result->lower_left = ll;
    result->upper_right = ur;
    result->left = make_density_tree(ndim, nlt, pts, ll, new_ur);
    result->right = make_density_tree(ndim, npts-nlt, pts+nlt, new_ll, ur);

    free(new_ll);
    free(new_ur);
    
    return result;
  }
}

tree *
make_density_tree(size_t ndim, size_t npts, double **input_pts, double *ll, double *ur) {
  double **pts = malloc(npts*sizeof(double*));
  tree *result;
  assert(pts != 0);

  memcpy(pts, input_pts, npts*sizeof(double *));

  result = internal_make_density_tree(ndim, npts, pts, ll, ur);

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
find_cell(size_t ndim, double *pt, tree *t) {
  if (!in_bounds(ndim, pt, t->lower_left, t->upper_right)) {
    return NULL;
  } else if (t->left == NULL || t->right == NULL) {
    return t;
  } else { 
    tree *maybe_left = find_cell(ndim, pt, t->left);
    if (maybe_left != NULL) {
      return maybe_left;
    } else {
      return find_cell(ndim, pt, t->right);
    }
  }
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

double 
jump_probability(size_t ndim, size_t npts, double *pt, tree *t) {
  tree *cell = find_cell(ndim, pt, t);
  double vol = 1.0;
  size_t i;

  for (i = 0; i < ndim; i++) {
    vol *= (cell->upper_right[i] - cell->lower_left[i]);
  }

  return ((double)cell->npts)/(vol*npts);
}
