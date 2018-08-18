#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "lsdtree.h"
#include "queue.h"

#define LSD_BUCKET_SIZE 100

#define LSD_NODE_BUFFER_DEPTH 10

struct lsdnode
{
	int splitdim;

	num_t splitpos;

	struct lsdnode *left, *right;

	num_t **bucket;

	int bucketsize;
};

struct lsdtree
{
	struct lsdnode *root;

	int numofdim;

	int bucketsize;

};

void lsd_create_subnodes(struct lsdnode * node, num_t **values, int dim,  int size, int bucketsize);
void lsd_sort(num_t **values, int dim, int size);
void lsd_free_node(struct lsdnode *node);


struct lsdtree * lsd_create(num_t **values,   int dim, int size)
{
	struct lsdtree *tree;

	tree = malloc(sizeof *tree);
	tree->numofdim = dim;
	tree->root = malloc(sizeof *(tree->root));

	tree->root->splitdim = 0;

	int bucketsize = LSD_BUCKET_SIZE;
	tree->bucketsize = bucketsize;

	lsd_create_subnodes(tree->root, values,  dim,size, bucketsize);

	return tree;
}

struct lsdtree * lsd_create_empty(int dim)
{
	struct lsdtree *tree;

	tree = malloc(sizeof *tree);
	tree->numofdim = dim;

	int bucketsize = LSD_BUCKET_SIZE;
	tree->bucketsize = bucketsize;

	tree->root = malloc(sizeof *(tree->root));
	tree->root->splitdim = 0;
	tree->root->bucket = malloc(bucketsize * sizeof(num_t*));
	tree->root->bucketsize = 0;
	tree->root->left = NULL;
	tree->root->right = NULL;

	return tree;
};


void lsd_create_subnodes(struct lsdnode * node, num_t **values, int dim,  int size, int bucketsize)
{
	if (size <= bucketsize)
	{
        num_t ** bucket = malloc(bucketsize * sizeof(num_t*));

        for (int i = 0; i < size; i++) bucket[i] = values[i];

        node->bucket = bucket;
		node->bucketsize = size;
		node->left = NULL;
		node->right = NULL;
		return;
	}

	// Get the dimesnion to split across
	int splitdim = node->splitdim;

	// Calculate Average across the split dimension
	num_t average = 0;
	for (int i = 0; i < size; i++)
	{
        num_t * val = values[i];
		average += val[splitdim];
	}
	average /= size;
	node->splitpos = average;

	// Sort the data across the specific dimension
	lsd_sort(values,splitdim,size);

	// find the datapoint to be split by
	int splitvpos = 0;
	for (int i = 0; i < size; i++)
	{
	    num_t *val = values[i];
		if (val[splitdim] > average)
		{
			splitvpos = i-1;
			break;
		}
	}

	node->left = malloc(sizeof *(node->left));
	node->right = malloc(sizeof *(node->right));

	int newdim = (splitdim == dim-1)? 0 : splitdim+1;

	node->left->splitdim = newdim;
	node->right->splitdim = newdim;

	lsd_create_subnodes(node->left, values, dim, splitvpos + 1, bucketsize);
	lsd_create_subnodes(node->right, values + splitvpos + 1, dim, size - splitvpos - 1, bucketsize);

}

int lsd_sort_partition(num_t **values, int dim, int size)
{
	num_t *pivot = values[size-1];

	int i = -1;

	for (int j = 0; j < size - 1; j++)
	{
		num_t *val_j = values[j];

		if (val_j[dim] < pivot[dim])
		{
			i++;
			num_t *val_i = values[i];
			values[j] = val_i;
			values[i] = val_j;
		}
	}
	num_t * val_i = values[i+1];
	values[size-1] = val_i;
	values[i+1] = pivot;

	return i + 1;
}

void lsd_sort(num_t **values, int dim, int size)
{
	if (size > 0)
	{
		int p = lsd_sort_partition(values, dim, size);

		lsd_sort(values, dim, p);
		lsd_sort(values + p+1, dim,  size - p - 1 );
	}
}

void lsd_print(struct lsdtree * tree)
{
    struct lsdnode *root = tree->root;

    while (root->left != NULL)
    {
        printf("%d %f\n", root->splitdim, root->splitpos);

        root = root->left;
    }

}

void lsd_insert(struct lsdtree * tree, num_t * value)
{
    struct lsdnode * node = tree->root;

    int maxbucketsize = tree->bucketsize;

    while (node->left != NULL)
    {
        int splitdim = node->splitdim;
        num_t val = value[splitdim];
        if (val <= node->splitpos)
            node = node->left;
        else
            node = node->right;
    }

    if (node->bucketsize == maxbucketsize)
    {
        num_t * values[maxbucketsize + 1];
        for (int i = 0; i < maxbucketsize; i++)
            values[i] = node->bucket[i];
        values[maxbucketsize] = value;

        free(node->bucket);

        lsd_create_subnodes(node,values,tree->numofdim,maxbucketsize + 1,maxbucketsize);
    }
    else
    {
        node->bucket[node->bucketsize] = value;
        (node->bucketsize)++;
    }


}

void lsd_free_node(struct lsdnode *node)
{
    if (node->left != NULL)
    {
        lsd_free_node(node->left);
        lsd_free_node(node->right);
    }
    else
    {
        free(node->bucket);
    }

    free(node);
}

void lsd_free(struct lsdtree *tree)
{
    lsd_free_node(tree->root);

    free(tree);
}

int lsd_is_point_contained(num_t * point, num_t * searchrange, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        num_t low = searchrange[i*2];
        num_t high = searchrange[i*2+1];

        if (point[i] < low || point[i] > high)
            return 0;
    }

    return 1;
}

int lsd_range_node_count(struct lsdnode *node, num_t * searchrange, int dim)
{
    int ret = 0;

    if (node->left == NULL)
    {
        for (int i = 0; i < node->bucketsize; i++)
        {
            if (lsd_is_point_contained(node->bucket[i], searchrange,dim))
                ret++;
        }
        return ret;
    }

    int sdim = node->splitdim;
    num_t spos = node->splitpos;

    num_t xlow = searchrange[sdim*2];
    num_t xhigh = searchrange[sdim*2+1];

    if(xlow < spos)
    {
        ret += lsd_range_node_count(node->left,searchrange, dim);
    }
    if (xhigh > spos)
    {
        ret += lsd_range_node_count(node->right, searchrange, dim);
    }

    return ret;
}

int lsd_range_count(struct lsdtree *tree, num_t * searchrange)
{
    return lsd_range_node_count(tree->root, searchrange, tree->numofdim);
}

int lsd_is_rect_intersecting(num_t * rect, num_t * searchrect, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        int p = 2*i;
        num_t amin = rect[p];
        num_t amax = rect[p + 1];

        num_t bmin = searchrect[p];
        num_t bmax = searchrect[p + 1];

        if(amax < bmin || amin > bmax)
            return 0;
    }

    return 1;
}

int lsd_cts_range_node_count(struct lsdnode *node, num_t * searchrect, int dim)
{
    int ret = 0;

    if (node->left == NULL)
    {
        for (int i = 0; i < node->bucketsize; i++)
        {
            if (lsd_is_rect_intersecting(node->bucket[i], searchrect,dim))
                ret++;
        }
        return ret;
    }

    int sdim = node->splitdim;
    num_t spos = node->splitpos;

    int sndim = sdim / 2;

    num_t xlow = searchrect[sndim*2];
    num_t xhigh = searchrect[sndim*2+1];

    if (sdim % 2 == 0)
    {
        ret += lsd_cts_range_node_count(node->left,searchrect,dim);
        if (spos <= xhigh)
            ret += lsd_cts_range_node_count(node->right,searchrect,dim);
    }
    else
    {
        ret += lsd_cts_range_node_count(node->right,searchrect,dim);
        if (spos >= xlow)
            ret += lsd_cts_range_node_count(node->left,searchrect,dim);

    }

    return ret;
}

int lsd_cts_range_count(struct lsdtree *tree, num_t * searchrect)
{
    return lsd_cts_range_node_count(tree->root,searchrect,tree->numofdim / 2);
}


int lsd_par_range_node_count(num_t * searchrange, int dim, queue * q, long * tcount, long * ncount)
{
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    struct lsdnode * node;

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
            while (node->left != NULL)
            {
                int sdim = node->splitdim;
                num_t spos = node->splitpos;

                num_t xlow = searchrange[sdim*2];
                num_t xhigh = searchrange[sdim*2+1];

                int goneleft = 0;
                struct lsdnode * nnode = NULL;

                if(xlow < spos)
                {
                    nnode = node->left;
                    goneleft = 1;
                }

                if (xhigh > spos)
                {
                    if (goneleft)
                    {
                        #pragma omp critical
                        queue_push(q,node->right);
                    }
                    else
                        nnode = node->right;
                }

                node = nnode;
            }

            int ret = 0;
            for (int i = 0; i < node->bucketsize; i++)
            {
                if (lsd_is_point_contained(node->bucket[i], searchrange,dim))
                    ret++;
            }

            #pragma omp atomic
            *ncount += ret;

        }

        struct lsdnode * nextnode;

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

int lsd_par_range_count(struct lsdtree *tree, num_t * searchrange)
{
    queue * q = queue_create();
    queue_push(q, tree->root);
    long tcount = 0;
    long ncount = 0;

    #pragma omp parallel num_threads(4)
    lsd_par_range_node_count(searchrange, tree->numofdim, q, &tcount, &ncount);

    queue_free(q);

    return ncount;
}
int lsd_par_cts_range_node_count(num_t * searchrect, int dim, queue * q, long * tcount, long * ncount)
{
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    struct lsdnode * node;

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
            while (node->left != NULL)
            {
                int sdim = node->splitdim;
                num_t spos = node->splitpos;

                int sndim = sdim / 2;

                num_t xlow = searchrect[sndim*2];
                num_t xhigh = searchrect[sndim*2+1];

                struct lsdnode * nnode = NULL;

                if (sdim % 2 == 0)
                {
                    nnode = node->left;
                    if (spos <= xhigh)
                    {
                        #pragma omp critical
                        queue_push(q, node->right);
                    }
                }
                else
                {
                    nnode = node->right;
                    if (spos >= xlow)
                    {
                        #pragma omp critical
                        queue_push(q, node->left);
                    }

                }

                node = nnode;
            }

            int ret = 0;
            for (int i = 0; i < node->bucketsize; i++)
            {
                if (lsd_is_rect_intersecting(node->bucket[i], searchrect,dim))
                    ret++;
            }

            #pragma omp atomic
            *ncount += ret;
        }

        struct lsdnode * nextnode;

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

int lsd_par_cts_range_count(struct lsdtree *tree, num_t * searchrect)
{
    queue * q = queue_create();
    queue_push(q, tree->root);
    long tcount = 0;
    long ncount = 0;

    #pragma omp parallel num_threads(4)
    lsd_par_cts_range_node_count(searchrect, tree->numofdim / 2, q, &tcount, &ncount);

    queue_free(q);

    return ncount;
}


