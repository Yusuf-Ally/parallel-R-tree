/* Created by Dylan Jacobsen */

#ifndef _LSDTREE_H_
#define _LSDTREE_H_

typedef double num_t;

struct lsdtree;

struct lsdtree * lsd_create(num_t **values,   int dim, int size);

struct lsdtree * lsd_create_empty(int dim);

void lsd_insert(struct lsdtree * tree, num_t * value);

void lsd_free (struct lsdtree *tree);

// search range for 2 dimensional has to be in the format (xl,xh,yl,yh).
// for 3 dimension (xl,xh,yl,yh,zl,zh) etc for higher dimensions.
int lsd_range_count(struct lsdtree *tree, num_t * searchrange);
int lsd_cts_range_count(struct lsdtree *tree, num_t * searchrect);

int lsd_par_range_count(struct lsdtree *tree, num_t * searchrange);
int lsd_par_cts_range_count(struct lsdtree *tree, num_t * searchrect);

int lsd_par_range_count2(struct lsdtree *tree, num_t * searchrange);
int lsd_par_cts_range_count2(struct lsdtree *tree, num_t * searchrect);

#endif	/* _LDSTREE_H_ */
