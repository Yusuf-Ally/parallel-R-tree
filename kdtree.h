
#ifndef _KDTREE_H_
#define _KDTREE_H_

struct kdtree;

typedef double knum_t;

struct kdtree * kd_create( knum_t ** values, int dim, int size);

struct kdtree * kd_open(const char * filename);

void kd_free(struct kdtree * tree);

void kd_insert(struct kdtree * tree, knum_t * value);


int kd_cts_range_count(struct kdtree *tree, knum_t * searchrect);
int kd_range_count(struct kdtree *tree, knum_t * searchrange)

// Slow functions
int kd_par_range_count(struct kdtree *tree, knum_t * searchrange);
int kd_par_cts_range_count(struct kdtree *tree, knum_t * searchrect);

// Fast functions. Is faster than the ones above
int kd_par_range_count2(struct kdtree *tree, knum_t * searchrange);
int kd_par_cts_range_count2(struct kdtree *tree, knum_t * searchrect);

#endif
