#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "kdtree.h"
#include "queue.h"
#include "queuea.h"

#define KD_NUMOFTHREADS 2

struct kdnode
{
    knum_t * point;

    struct kdnode * left, * right;

    int dim;
};

struct kdtree
{
    struct kdnode * root;

    int numofdim;

};

struct kdnodeH
{
    struct kdnode * node;
    int height;
};

void kd_sort(knum_t **values, int dim, int size);
struct kdnode * kd_create_subnodes(knum_t **values, int dim,  int size, int height);

struct kdnode * kd_create_node()
{
    struct kdnode * node;
    node = malloc(sizeof *node);
    node->point = NULL;
    node->left = NULL;
    node->right = NULL;

    return node;
}

struct kdtree * kd_create(knum_t ** values, int dim, int size)
{
    struct kdtree * tree;
    tree = malloc(sizeof *tree);

    tree->numofdim = dim;

    tree->root = kd_create_subnodes(values,dim,size,0);

    return tree;
}

struct kdnode * kd_create_subnodes(knum_t **values, int dim,  int size, int height)
{
    if (size <= 0)
        return NULL;

    int splitdim = height % dim;

	// Sort the data across the specific dimension
	kd_sort(values,splitdim,size);

	// find the datapoint to be split by
	int median = size / 2;

    knum_t * point = malloc(dim * sizeof(knum_t));
    memcpy(point, values[median], dim * sizeof(knum_t));

    struct kdnode * node = kd_create_node();

    node->point = point;

    node->left = kd_create_subnodes(values, dim, median, height+1);

    node->right = kd_create_subnodes(values+median+1, dim, (size % 2 == 0)? median - 1 : median, height+1);

    return node;
}

int kd_sort_partition(knum_t **values, int dim, int size)
{
	knum_t *pivot = values[size-1];

	int i = -1;

	for (int j = 0; j < size - 1; j++)
	{
		knum_t *val_j = values[j];

		if (val_j[dim] < pivot[dim])
		{
			i++;
			knum_t *val_i = values[i];
			values[j] = val_i;
			values[i] = val_j;
		}
	}
	knum_t * val_i = values[i+1];
	values[size-1] = val_i;
	values[i+1] = pivot;

	return i + 1;
}

void kd_sort(knum_t **values, int dim, int size)
{
	if (size > 0)
	{
		int p = kd_sort_partition(values, dim, size);

		kd_sort(values, dim, p);
		kd_sort(values + p+1, dim,  size - p - 1 );
	}
}


void kd_insert(struct kdtree * tree, knum_t * value)
{
    struct kdnode * node = tree->root;

    int depth = 0;

    struct kdnode * newnode = kd_create_node();
    newnode->point = malloc(sizeof(knum_t) * tree->numofdim);
    memcpy(newnode->point, value, sizeof(knum_t) * tree->numofdim);

    while (1)
    {
        int dim = depth % tree->numofdim;

        if (value[dim] <= node->point[dim])
        {
            if (node->left == NULL)
            {
                node->left = newnode;
                break;
            }
        }
        else
        {
            if (node->right == NULL)
            {
                node->right = newnode;
                break;
            }
        }

        depth++;
    }
}

void kd_free_node(struct kdnode * node)
{
    if (node->left != NULL)
        kd_free_node(node->left);

    if (node->right != NULL)
        kd_free_node(node->right);

    free(node->point);

    free(node);
}

void kd_free(struct kdtree * tree)
{
    kd_free_node(tree->root);
}

int kd_is_point_contained(knum_t * point, knum_t * searchrange, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        knum_t low = searchrange[i*2];
        knum_t high = searchrange[i*2+1];

        if (point[i] < low || point[i] > high)
            return 0;
    }

    return 1;
}

int kd_range_node_count(struct kdnode *node, knum_t * searchrange, int dim, int height)
{
    int ret = 0;

    if (node == NULL)
        return 0;

    if (kd_is_point_contained(node->point,searchrange, dim))
        ret += 1;

    int sdim = height % dim;
    knum_t spos = node->point[sdim];

    knum_t xlow = searchrange[sdim*2];
    knum_t xhigh = searchrange[sdim*2+1];

    if(xlow <= spos)
    {
        ret += kd_range_node_count(node->left,searchrange, dim, height + 1);
    }
    if (xhigh >= spos)
    {
        ret += kd_range_node_count(node->right, searchrange, dim, height + 1);
    }

    return ret;
}

int kd_range_count(struct kdtree *tree, knum_t * searchrange)
{
    return kd_range_node_count(tree->root, searchrange, tree->numofdim, 0);
}



int kd_is_rect_intersecting(knum_t * rect, knum_t * searchrect, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        int p = 2*i;
        knum_t amin = rect[p];
        knum_t amax = rect[p + 1];

        knum_t bmin = searchrect[p];
        knum_t bmax = searchrect[p + 1];

        if(amax < bmin || amin > bmax)
            return 0;
    }

    return 1;
}

int kd_cts_range_node_count(struct kdnode *node, knum_t * searchrect, int dim, int height)
{
    int ret = 0;

    if (node == NULL)
        return 0;

    if (kd_is_rect_intersecting(node->point, searchrect,dim))
        ret++;

    int sdim = height % (dim*2);
    knum_t spos = node->point[sdim];

    int sndim = sdim / 2;

    knum_t xlow = searchrect[sndim*2];
    knum_t xhigh = searchrect[sndim*2+1];

    if (sdim % 2 == 0)
    {
        ret += kd_cts_range_node_count(node->left,searchrect,dim, height+1);
        if (spos <= xhigh)
            ret += kd_cts_range_node_count(node->right,searchrect,dim, height+1);
    }
    else
    {
        ret += kd_cts_range_node_count(node->right,searchrect,dim, height+1);
        if (spos >= xlow)
            ret += kd_cts_range_node_count(node->left,searchrect,dim, height+1);

    }

    return ret;
}

int kd_cts_range_count(struct kdtree *tree, knum_t * searchrect)
{
    return kd_cts_range_node_count(tree->root,searchrect,tree->numofdim / 2, 0);
}

struct kdnodeH * kd_create_par_node(struct kdnode * node, int height)
{
    struct kdnodeH * ndh = malloc(sizeof *ndh);
    ndh->height = height;
    ndh->node = node;

    return ndh;
}

void kd_par_range_node_count(knum_t * searchrange, int dim, queuea * q, long * tcount, long * ncount)
{
    //int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    struct kdnode * node;

    #pragma omp critical
    node = queuea_pop(q);

    if (node == NULL)
    {
        #pragma omp atomic
        (*tcount)++;

    }

    while (*tcount < thread_count)
    {
        if (node != NULL)
        {
            int ret = 0;
            int height = node->dim;

            struct kdnode * knode = node;

            //queue * q2 = queue_create();

            while (knode != NULL)
            {
                if (kd_is_point_contained(knode->point,searchrange, dim))
                    ret += 1;

                int sdim = height % dim;

                knum_t spos = knode->point[sdim];

                knum_t xlow = searchrange[sdim*2];
                knum_t xhigh = searchrange[sdim*2+1];

                int goneleft = 0;
                struct kdnode * nnode = NULL;

                if(xlow <= spos)
                {
                    nnode = knode->left;
                    goneleft = 1;
                }

                if (xhigh >= spos)
                {
                    if (goneleft)
                    {
                        if (knode->right != NULL)
                        {
                            knode->right->dim = height+1;
                            #pragma omp critical
                            queuea_push(q,knode->right);
                        }
                    }
                    else
                        nnode = knode->right;
                }

                knode = nnode;
                height++;
            }

            #pragma omp atomic
            *ncount += ret;

            //#pragma omp critical
            //queue_merge(q,q2);
            //queue_free(q2);
            //free(node);
        }

        struct kdnode * nextnode;

        #pragma omp critical
        nextnode = queuea_pop(q);

        if (nextnode == NULL && node != NULL)
        {
            #pragma omp atomic
            (*tcount)++;
        }
        if (nextnode != NULL && node == NULL)
        {
            #pragma omp atomic
            (*tcount)--;
        }

        node = nextnode;
    }
}

int kd_par_range_count(struct kdtree *tree, knum_t * searchrange)
{
    queuea * q = queuea_create(1000000);
    tree->root->dim = 0;
    queuea_push(q, tree->root);
    long tcount = 0;
    long ncount = 0;

    #pragma omp parallel num_threads(KD_NUMOFTHREADS)
    kd_par_range_node_count(searchrange, tree->numofdim, q, &tcount, &ncount);

    queuea_free(q);

    return ncount;
}

int kd_get_max_balanced_height(struct kdtree * tree)
{
    struct kdnode * node = tree->root;
    int h = 0;

    while (node != NULL)
    {
        node = node->left;
        h++;
    }
    return h;
}


void kd_par_range_node_count2(struct kdnode *node, knum_t * searchrange, int dim, int height, int maxheight, long *ncount)
{
    //int my_rank = omp_get_thread_num();
    //int thread_count = omp_get_num_threads();

    if (node == NULL)
        return;

    if (kd_is_point_contained(node->point,searchrange, dim))
    {
        #pragma omp atomic
        (*ncount)++;
    }

    int sdim = height % dim;
    knum_t spos = node->point[sdim];

    knum_t xlow = searchrange[sdim*2];
    knum_t xhigh = searchrange[sdim*2+1];

    int gone = 0;

    if(xlow <= spos)
    {
        if (node->left != NULL)
        {
            kd_par_range_node_count2(node->left,searchrange, dim, height + 1, maxheight, ncount);
            gone = 1;
        }
    }
    if (xhigh >= spos)
    {
        if(node->right != NULL)
        {
            if (gone && height < maxheight)
            {
                #pragma omp task
                kd_par_range_node_count2(node->right, searchrange, dim, height + 1, maxheight, ncount);
            }
            else
                kd_par_range_node_count2(node->right, searchrange, dim, height + 1, maxheight, ncount);

        }
    }
}


int kd_par_range_count2(struct kdtree *tree, knum_t * searchrange)
{
    long ncount = 0;

    int maxheight = kd_get_max_balanced_height(tree) - 5;

    #pragma omp parallel num_threads(KD_NUMOFTHREADS) shared(ncount)
    {
        #pragma omp master
        kd_par_range_node_count2(tree->root,searchrange,tree->numofdim,0,maxheight,&ncount);
    }

    return ncount;
}


void kd_par_cts_range_node_count(knum_t * searchrect, int dim, queue * q, long * tcount, long * ncount)
{
    //int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    struct kdnodeH * node;

    #pragma omp critical
    node = queue_pop(q);

    if (node == NULL)
    {
        #pragma omp atomic
        (*tcount)++;

    }

    while (*tcount < thread_count)
    {
        if (node != NULL)
        {
            int ret = 0;
            int height = node->height;

            struct kdnode * knode = node->node;


            while (knode != NULL)
            {
                if (kd_is_rect_intersecting(knode->point, searchrect,dim))
                    ret++;

                int sdim = height % (dim*2);
                knum_t spos = knode->point[sdim];


                int sndim = sdim / 2;

                knum_t xlow = searchrect[sndim*2];
                knum_t xhigh = searchrect[sndim*2+1];

                struct kdnode * nnode = NULL;

                if (sdim % 2 == 0)
                {
                    nnode = knode->left;
                    if (spos <= xhigh)
                    {
                        #pragma omp critical
                        queue_push(q, kd_create_par_node(knode->right, height+1));
                    }
                }
                else
                {
                    nnode = knode->right;
                    if (spos >= xlow)
                    {
                        #pragma omp critical
                        queue_push(q, kd_create_par_node(knode->left, height+1));
                    }

                }
                height++;
                knode = nnode;
            }

            #pragma omp atomic
            *ncount += ret;

            free(node);
        }

        struct kdnodeH * nextnode;

        #pragma omp critical
        nextnode = queue_pop(q);

        if (nextnode == NULL && node != NULL)
        {
            #pragma omp atomic
            (*tcount)++;
        }
        if (nextnode != NULL && node == NULL)
        {
            #pragma omp atomic
            (*tcount)--;
        }

        node = nextnode;
    }
}

int kd_par_cts_range_count(struct kdtree *tree, knum_t * searchrect)
{
    queue * q = queue_create();
    queue_push(q, kd_create_par_node(tree->root, 0));
    long tcount = 0;
    long ncount = 0;

    #pragma omp parallel num_threads(KD_NUMOFTHREADS)
    kd_par_cts_range_node_count(searchrect, tree->numofdim / 2, q, &tcount, &ncount);

    queue_free(q);

    return ncount;
}



void kd_par_cts_range_node_count2(struct kdnode *node, knum_t * searchrect, int dim, int height, int maxheight, long * ncount)
{
    if (node == NULL)
        return;

    if (kd_is_rect_intersecting(node->point, searchrect,dim))
    {
        #pragma omp atomic
        (*ncount)++;
    }

    int sdim = height % (dim*2);
    knum_t spos = node->point[sdim];

    int sndim = sdim / 2;

    knum_t xlow = searchrect[sndim*2];
    knum_t xhigh = searchrect[sndim*2+1];

    if (sdim % 2 == 0)
    {
        kd_par_cts_range_node_count2(node->left,searchrect,dim, height+1, maxheight, ncount);

        if (spos <= xhigh)
        {
            if (height < maxheight)
            {
                #pragma omp task untied
                kd_par_cts_range_node_count2(node->right,searchrect,dim, height+1, maxheight, ncount);
            }
            else
                kd_par_cts_range_node_count2(node->right,searchrect,dim, height+1, maxheight, ncount);
        }
    }
    else
    {
        kd_par_cts_range_node_count2(node->right,searchrect,dim, height+1, maxheight, ncount);

        if (spos >= xlow)
        {
            if (height < maxheight)
            {
                #pragma omp task untied
                kd_par_cts_range_node_count2(node->left,searchrect,dim, height+1, maxheight, ncount);
            }
            else
                kd_par_cts_range_node_count2(node->left,searchrect,dim, height+1, maxheight, ncount);
        }

    }
}

int kd_par_cts_range_count2(struct kdtree *tree, knum_t * searchrect)
{
    long ncount = 0;

    int maxheight = kd_get_max_balanced_height(tree) - 7;

    #pragma omp parallel num_threads(KD_NUMOFTHREADS) shared(ncount)
    {
        #pragma omp master
        kd_par_cts_range_node_count2(tree->root,searchrect,tree->numofdim / 2, 0, maxheight, &ncount);
    }

    return ncount;
}

