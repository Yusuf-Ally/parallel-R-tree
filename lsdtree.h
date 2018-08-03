#ifndef _LSDTREE_H_
#define _LSDTREE_H_

struct lsdtree;

typedef double num_t;

struct lsdtree * lsd_create(int dim, int bucketsize, const char * filename);

struct lsdtree * lsd_open(char * filename);

void lsd_flush(struct lsdtree * tree);

void lsd_free(struct lsdtree * tree);

void lsd_insert(struct lsdtree * tree, num_t * value);


#endif
