#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kdtree.h"

#define KDHEADERSIZE 128

struct kdnode
{
    knum_t * point;

    struct kdnode * left, * right;
};

struct kdtree
{
    struct kdnode * root;

    int numofdim;

    int fheadersize;

    int fnodesize;

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

struct kdtree * kd_create(int dim, knum_t ** values, int size)
{
    struct kdtree * tree;
    tree = malloc(sizeof *tree);

    tree->numofdim = dim;

    tree->fheadersize = KDHEADERSIZE;

    tree->fnodesize = 16 + dim * sizeof (knum_t);

    tree->root = kd_create_subnodes(values,dim,size,0);
}

struct kdnode * kd_create_subnodes(knum_t **values, int dim,  int size, int height)
{
    if(size <= 0)
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

int kd_nodes_count(struct kdnode * node)
{
    if (node == NULL)
        return 0;

    return 1 + kd_nodes_count(node->left) + kd_nodes_count(node->right);
}


void kd_node_make_array(struct kdnode * node, struct kdnode ** nodes, long * nodecount)
{
    if (node == NULL)
        return;

    nodes[*nodecount] = node;

    (*nodecount)++;

    kd_node_make_array(node->left,nodes,nodecount);
    kd_node_make_array(node->right,nodes,nodecount);
}

long kd_get_node_position(struct kdnode * node, struct kdnode ** nodes, long size)
{
    for (int i = 0; i < size; i++)
    {
        if (nodes[i] == node)
            return i;
    }
    return -1;
}

void kd_write_node(char * buffer, struct kdnode * node, struct kdnode ** nodes, long size, long dim)
{
    int wpos = 0;

    *(long *)(buffer + wpos) = kd_get_node_position(node->left,nodes,size);
    wpos += sizeof(long);

    *(long *)(buffer + wpos) = kd_get_node_position(node->right,nodes,size);
    wpos += sizeof(long);

    for (int i = 0; i < dim; i++)
    {
        *(knum_t *)(buffer + wpos) = node->point[i];
        wpos += sizeof(knum_t);
    }
}

void kd_save(struct kdtree * tree, const char * filename)
{
    int icount = kd_nodes_count(tree->root);

    struct kdnode ** nodes;
    nodes = malloc(sizeof (struct kdnode *) * icount);

    int nodecount = 0;
    kd_node_make_array(tree->root,nodes,&nodecount);

    int nodesize = 16 + sizeof(knum_t)*(tree->numofdim);

    int wpos = 0;
    char * header = malloc(tree->fheadersize);

    *(int *)(header + wpos) = tree->fheadersize;
    wpos += sizeof(int);

    *(int *)(header + wpos) = tree->numofdim;
    wpos += sizeof(int);

    *(int *)(header + wpos) = nodesize;
    wpos += sizeof(int);

    *(long *)(header + wpos) = icount;
    wpos += sizeof(long);

    FILE * f = fopen(filename,"wb+");

    fwrite(header,1,tree->fheadersize,f);

    int pagenodecount = 50;
    long pagesize = pagenodecount * nodesize;

    char * page = malloc(pagesize);
    int wpage = 0;

    int noofpages = icount / pagenodecount + 1;

    for (int i = 0; i < noofpages; i++)
    {
        wpage = 0;
        memset(page, 0, pagesize);

        for (int j = 0; j < pagenodecount; i++)
        {
            if (i * pagenodecount + j < icount)
            {

            }
            else
            {
                pagesize = j * nodesize;
                break;
            }
        }
    }


    free(nodes);
    free(page);
}
