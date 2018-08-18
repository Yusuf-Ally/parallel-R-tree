
#ifndef _KDTREE_H_
#define _KDTREE_H_

struct kdtree;

typedef double knum_t;

struct kdtree * kd_create( knum_t ** values, int dim, int size);

struct kdtree * kd_open(const char * filename);

void kd_free(struct kdtree * tree);

void kd_insert(struct kdtree * tree, knum_t * value);


#endif
