#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lsdtree.h"

#define LSD_VERSION 1

#define LSD_TREE_PAGE_HEIGHT 5
#define LSD_TREE_HEADER_SIZE 1024
#define LSD_TREE_NODE_SIZE 64

struct lsdnode
{
    // The dimension to split across
	int splitdim;

	// the position to split across
	num_t splitpos;

	// The left and right sub-nodes
	struct lsdnode *left, *right;

    // If the node is at the end of a tree page define the page number of the new page
	long pagenumber;

	// the index of the bucket in the bucket page.
	long bucketindex;

    // the size of a bucket
	int bucketsize;
};

struct lsdpage
{
    struct lsdnode * root;
    long pagenumber;
};

struct lsdtree
{
    struct lsdpage * rootpage;

	int numofdim;

	int bucketsize;

    long numofpages;

    long numofbucketpages;

    int pageheight;
    // The size of the header in bytes
    int headersize;
    // the size of each page in bytes
    int pagesize;
    // the size of each node in bytes
    int nodesize;

    int version;

	char * filename;

	FILE * file;

	FILE * bfile;

};

struct lsdpage * lsd_get_treepage(struct lsdtree * tree, long pagenumber);

struct lsdnode * lsd_node_create()
{
    struct lsdnode * node;
    node = malloc(sizeof *node);
    node->splitdim = 0;
    node->splitpos = 0.0;
    node->bucketindex = -1;
    node->bucketsize = 0;
    node->left = NULL;
    node->right = NULL;
    node->pagenumber = -1;
};

struct lsdpage * lsd_page_create()
{
    struct lsdpage * page;
    page = malloc(sizeof *page);

    page->pagenumber = -1;

    //page->root = lsd_node_create();

};


int lsd_get_page_size(int pageheight)
{
    int res = 0;
    for (int i = 0; i < pageheight; i++)
        res += 1 << i;

    return res * LSD_TREE_NODE_SIZE;
}

int lsd_get_numofdims(struct lsdtree * tree)
{
    return tree->numofdim;
}

struct lsdtree * lsd_create(int dim, int bucketsize, const char * filename)
{
    struct lsdtree * tree;

    tree = malloc(sizeof *tree);

    tree->headersize = LSD_TREE_HEADER_SIZE;

    tree->numofdim = dim;

    tree->bucketsize = bucketsize;

    tree->filename = malloc(strlen(filename) + 1);
    strcpy(tree->filename, filename);

    tree->file = fopen(filename, "wb+");

    tree->numofpages = 1;

    tree->pageheight = LSD_TREE_PAGE_HEIGHT;

    tree->pagesize = lsd_get_page_size(LSD_TREE_PAGE_HEIGHT);

    tree->nodesize = LSD_TREE_NODE_SIZE;

    tree->version = LSD_VERSION;

    tree->rootpage = lsd_page_create();
    tree->rootpage->pagenumber = 0;
    tree->rootpage->root = lsd_node_create();
    tree->rootpage->root->bucketindex = 0;

    tree->numofbucketpages = 0;

    lsd_flush(tree);


    char bucketfile[strlen(tree->filename) + 10];
    strcpy(bucketfile, tree->filename);
    strcat(bucketfile, ".pgs");

    tree->bfile = fopen(bucketfile, "wb+");

    lsd_new_bucketpage(tree,NULL,0);

    return tree;
};

struct lsdtree * lsd_open(char * filename)
{
    FILE * f;

    f = fopen(filename,"rb+");

    fseek(f,0,SEEK_SET);

    struct lsdtree * tree;
    tree = malloc(sizeof *tree);

    tree->filename = malloc(strlen(filename) + 1);
    strcpy(tree->filename, filename);

    tree->file = f;

    int headersize = 0;
    int wpos = 0;

    wpos += sizeof(int);

    fread(&headersize, sizeof(int), 1, f);
    tree->headersize = headersize;

    fseek(f,0,SEEK_SET);

    char wdata[headersize];
    memset(wdata,0,headersize);

    fread(wdata, 1, headersize, f);

    int version = *(int *)(wdata+wpos);
    tree->version = version;
    wpos += sizeof(tree->version);

    if (tree->version != LSD_VERSION)
        return NULL;

    int numofdim = *(int *)(wdata+wpos);
    tree->numofdim = numofdim;
    wpos += sizeof (tree->numofdim);

    tree->bucketsize = *(int *)(wdata+wpos);
    wpos += sizeof (tree->bucketsize);

    tree->numofpages = *(long *)(wdata+wpos);
    wpos += sizeof(tree->numofpages);

    tree->numofbucketpages = *(long *)(wdata+wpos);
    wpos += sizeof(tree->numofbucketpages);

    tree->pageheight = *(int *)(wdata+wpos);
    wpos += sizeof(tree->pageheight);

    tree->pagesize = *(int *)(wdata+wpos);
    wpos += sizeof(tree->pagesize);

    tree->nodesize = *(int *)(wdata+wpos);
    wpos += sizeof(tree->nodesize);

    tree->rootpage = lsd_get_treepage(tree,0);

    char bucketfile[strlen(tree->filename) + 10];
    strcpy(bucketfile, tree->filename);
    strcat(bucketfile, ".pgs");

    tree->bfile = fopen(bucketfile, "rb+");


    return tree;
};

void lsd_flush(struct lsdtree * tree)
{
    FILE * f = tree->file;

    fseek(f, 0, SEEK_SET);

    // WRITE HEADER
    int headersize = tree->headersize;
    int wpos = 0;
    char wdata[headersize];
    memset(wdata,0,headersize);

    *(int *)(wdata+wpos) = headersize;
    wpos += sizeof(int);

    *(int *)(wdata+wpos) = tree->version;
    wpos += sizeof(int);

    *(int *)(wdata+wpos) = tree->numofdim;
    wpos += sizeof(tree->numofdim);

    *(int *)(wdata+wpos) = tree->bucketsize;
    wpos += sizeof(tree->bucketsize);

    *(long *)(wdata+wpos) = tree->numofpages;
    wpos += sizeof(tree->numofpages);

    *(long *)(wdata+wpos) = tree->numofbucketpages;
    wpos += sizeof(tree->numofbucketpages);

    *(int *)(wdata+wpos) = tree->pageheight;
    wpos += sizeof(tree->pageheight);

    *(int *)(wdata+wpos) = tree->pagesize;
    wpos += sizeof(tree->pagesize);

    *(int *)(wdata+wpos) = tree->nodesize;
    wpos += sizeof(tree->nodesize);

    fwrite(wdata, 1, headersize, f);;
    fflush(f);


    lsd_put_treepage(tree, tree->rootpage);
}



void lsd_free_node(struct lsdnode *node)
{
    if (node == NULL)
        return;

    if (node->left != NULL)
    {
        lsd_free_node(node->left);
        lsd_free_node(node->right);
    }

    free(node);
}

void lsd_free_page(struct lsdpage *page)
{
    if (page == NULL)
        return;

    lsd_free_node(page->root);

    free(page);
}

void lsd_free(struct lsdtree * tree)
{
    if (tree == NULL)
        return;

    lsd_flush(tree);

    lsd_free_page(tree->rootpage);

    fclose(tree->file);

    fclose(tree->bfile);

    free(tree->filename);

    free(tree);
}

void lsd_free_bucket(struct lsdtree *tree, num_t ** bucket)
{
    for (int i = 0; i < tree->bucketsize; i++)
    {
        free(bucket[i]);
    }
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

void lsd_create_subnodes(struct lsdtree * tree, struct lsdnode * node, num_t **values,int splitdim, int size)
{
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
		if (values[i][splitdim] > average)
		{
			splitvpos = i-1;
			break;
		}
	}

	//printf("NSplit %f\n", average);

	node->left = lsd_node_create();
	node->right = lsd_node_create();

    node->left->bucketindex = node->bucketindex;
    node->right->bucketindex = tree->numofbucketpages;

    node->left->bucketsize = splitvpos + 1;
    node->right->bucketsize = size - splitvpos - 1;

    node->bucketindex = -1;
    node->bucketsize = 0;
    node->splitdim = splitdim;

    lsd_put_bucketpage(tree,node->left->bucketindex,values,node->left->bucketsize);
    lsd_new_bucketpage(tree, &values[splitvpos+1], node->right->bucketsize);
}

void lsd_get_bucketpage(struct lsdtree * tree, int bucketindex, num_t ** bucket)
{
    FILE * f = tree->bfile;

    int nbucketsize = tree->bucketsize * tree->numofdim;
    int pagesize = ( sizeof(num_t) * nbucketsize);
    int ppos = bucketindex * pagesize;

    fseek(f,ppos, SEEK_SET);

    num_t cpage[nbucketsize];
    memset(cpage,0,pagesize);

    fread(cpage,1,pagesize,f);

    for (int i = 0; i < tree->bucketsize; i++)
    {
        num_t * val = malloc(sizeof(num_t) * tree->numofdim);
        for (int j = 0; j < tree->numofdim; j++)
            val[j] = cpage[i*(tree->numofdim) + j];

        bucket[i] = val;
    }
}

void lsd_put_bucketpage(struct lsdtree * tree, int bucketindex, num_t ** bucket, int bucketsize)
{
    FILE * f = tree->bfile;

    int nbucketsize = tree->bucketsize * tree->numofdim;
    int pagesize = ( sizeof(num_t) * nbucketsize);

    int ppos;
    if (bucketindex > -1)
        ppos = bucketindex * pagesize;
    else
    {
        ppos = tree->numofbucketpages * pagesize;
        (tree->numofbucketpages)++;
    }

    fseek(f,ppos, SEEK_SET);

    num_t cpage[nbucketsize];
    memset(cpage, 0, pagesize);


    for (int i = 0; i < bucketsize; i++)
    {
        int ipos = i * tree->numofdim;
        for (int j = 0; j < tree->numofdim; j++)
            cpage[ipos + j] = bucket[i][j];
    }

    fwrite(cpage, 1, pagesize, f);
    fflush(f);
}

void lsd_new_bucketpage(struct lsdtree * tree, num_t ** bucket, int bucketsize)
{
    lsd_put_bucketpage(tree,-1,bucket,bucketsize);
}



void lsd_insert(struct lsdtree * tree, num_t * value)
{
    struct lsdpage * page;
    struct lsdnode * parent, * node;
    int height = 1;

    int maxbucketsize = tree->bucketsize;

    page = tree->rootpage;
    node = page->root;

    while (node->left != NULL)
    {
        parent = node;

        int splitdim = node->splitdim;

        num_t val = value[splitdim];

        if (val <= node->splitpos)
            node = node->left;
        else
            node = node->right;

        height++;

        if (node->pagenumber >= 0)
        {
            long pagenumber = node->pagenumber;

            if (page->pagenumber > 0)
                lsd_free_page(page);

            page = lsd_get_treepage(tree, pagenumber);
            node = page->root;

            height = 1;
        }
    }

    num_t ** bucket = malloc(sizeof(num_t *) * (tree->bucketsize + 1));
    lsd_get_bucketpage(tree,node->bucketindex, bucket);

    int bcksize = node->bucketsize;

    if (bcksize < tree->bucketsize)
    {
        for(int i = 0; i < tree->numofdim;i++)
        {
            num_t val = value[i];
            bucket[bcksize][i] = val;
        }

        (node->bucketsize)++;

        lsd_put_bucketpage(tree,node->bucketindex, bucket, node->bucketsize);

        lsd_put_treepage(tree, page);

    }
    else
    {
        //printf("New Node\n");
        int splitdim;

        if (parent == NULL)
            splitdim = 0;
        else
            splitdim = (parent->splitdim + 1 < tree->numofdim) ? parent->splitdim + 1: 0;

        bucket[tree->bucketsize] = value;

        lsd_create_subnodes(tree,node,bucket, splitdim, tree->bucketsize + 1);

        if (height >= tree->pageheight)
        {
            struct lsdpage * rpage = lsd_page_create();
            rpage->root = node;

            struct lsdnode * rnode = lsd_node_create();
            rnode->pagenumber = tree->numofpages;

            if (node == parent->left)
                parent->left = rnode;
            else
                parent->right = rnode;

            lsd_new_treepage(tree, rpage);
        }

        lsd_put_treepage(tree, page);
    }


    if (page->pagenumber > 0)
        lsd_free_page(page);

    lsd_free_bucket(tree,bucket);
    free(bucket);
}

// fint the node in a list of nodes
int lsd_get_node_position(struct lsdnode * node, struct lsdnode ** nodes, int nodesnums)
{
    for (int i = 0; i < nodesnums; i++)
    {
        if (nodes[i] == node)
            return i;
    }
    return -1;
}

void lsd_write_node(char * buffer, struct lsdnode * node, struct lsdnode ** nodes, int nodesnums)
{
    int wpos = 0;

    num_t t = node->splitdim;

    *(int *)(buffer+wpos) = node->splitdim;
    wpos += sizeof(int);

    *(num_t *)(buffer+wpos) = node->splitpos;
    wpos += sizeof(num_t);

    *(int *)(buffer+wpos) = lsd_get_node_position(node->left, nodes, nodesnums);
    wpos += sizeof(int);

    *(int *)(buffer+wpos) = lsd_get_node_position(node->right, nodes, nodesnums);
    wpos += sizeof(int);

    *(long *)(buffer+wpos) = node->pagenumber;
    wpos += sizeof(long);

    *(long *)(buffer+wpos) = node->bucketindex;
    wpos += sizeof(long);

    *(int *)(buffer+wpos) = node->bucketsize;
    wpos += sizeof(int);
}

void lsd_read_node(char * buffer, struct lsdnode * node, struct lsdnode ** nodes)
{
    int wpos = 0;

    node->splitdim = *(int *)(buffer + wpos);
    wpos += sizeof (int);

    node->splitpos = *(num_t *)(buffer + wpos);
    wpos += sizeof (num_t);

    int lpos = *(int *)(buffer + wpos);
    node->left = (lpos == -1)? NULL : nodes[lpos];
    wpos += sizeof (int);

    int rpos = *(int *)(buffer + wpos);
    node->right = (rpos == -1)? NULL : nodes[rpos];
    wpos += sizeof (int);

    node->pagenumber = *(long *)(buffer + wpos);
    wpos += sizeof (long);

    node->bucketindex = *(long *)(buffer + wpos);
    wpos += sizeof (long);

    node->bucketsize = *(int *)(buffer + wpos);
    wpos += sizeof (int);
}


struct lsdpage * lsd_get_treepage(struct lsdtree * tree, long pagenumber)
{
    FILE * f = tree->file;

    long wpos = 0;

    int pagesize = tree->pagesize;

    wpos += tree->headersize + pagenumber * (pagesize);

    fseek(f,wpos,SEEK_SET);

    char * wpage = malloc(sizeof(char) * pagesize);

    fread(wpage,1,pagesize,f);

    int npos = 0;

    char c = wpage[0];

    int numofnodes = *(int*)wpage;
    npos += sizeof(int);

    struct lsdnode ** nodes;
    nodes = malloc(sizeof (struct lsdnode *) * numofnodes);

    for (int i = 0; i < numofnodes; i++)
        nodes[i] = lsd_node_create();

    for (int i = 0; i < numofnodes; i++)
    {
        lsd_read_node(wpage + npos,nodes[i], nodes);

        npos += tree->nodesize;
    }

    struct lsdpage * page = malloc(sizeof *page);
    page->root = nodes[0];
    page->pagenumber = pagenumber;

    free(nodes);
    free(wpage);

    return page;
};

// Counts the number of nodes up to a certain depth
int lsd_count_nodes(struct lsdnode * node, int depth)
{
    if (node->left == NULL || node->pagenumber > 0 || depth == 0)
        return 1;

    return lsd_count_nodes(node->left, depth - 1) + lsd_count_nodes(node->right, depth - 1) + 1;
}

// get a list of all nodes and subnodes up to a specific depth
void lsd_node_list(struct lsdnode * node, int depth, struct lsdnode ** result_list, int * numberofnodes)
{
    int nn = *numberofnodes;
    result_list[nn] = node;
    (*numberofnodes) += 1;

    if (node->left == NULL || node->pagenumber > 0 || depth == 0)
        return;

    lsd_node_list(node->left, depth - 1, result_list, numberofnodes);
    lsd_node_list(node->right, depth - 1, result_list, numberofnodes);
};

void lsd_put_treepage(struct lsdtree * tree, struct lsdpage * page)
{
    int inodes = lsd_count_nodes(page->root,tree->pageheight);

    char * wpage = malloc(tree->pagesize);
    memset(wpage,0,tree->pagesize);

    int pageindex = 0;

    struct lsdnode ** nodes;
    nodes = malloc((sizeof(struct lsdnode *)) * inodes);

    // Get the list of nodes to write
    int numofnodes = 0;
    lsd_node_list(page->root, tree->pageheight, nodes, &numofnodes);

    // write the number of nodes to the buffer
    *(int *)wpage = inodes;
    memcpy(wpage,&inodes, sizeof(int));
    pageindex += sizeof (int);

    for (int i = 0; i <inodes;i++)
    {
        char nodebuffer[tree->nodesize];
        memset(nodebuffer,0,tree->nodesize);

        struct lsdnode * nd = nodes[i];

        lsd_write_node(nodebuffer, nodes[i], nodes, inodes);

        memcpy((wpage + pageindex), nodebuffer, tree->nodesize);

        pageindex += tree->nodesize;
    }

    FILE * f = tree->file;

    if (page->pagenumber >= 0)
    {
        int pageposition = tree->headersize + page->pagenumber * tree->pagesize;

        fseek(f, pageposition, SEEK_SET);
        fwrite(wpage, 1, tree->pagesize, f);
        fflush(f);
    }
    else
    {
        int pageposition = tree->headersize + tree->numofpages * tree->pagesize;

        fseek(f, pageposition, SEEK_SET);
        fwrite(wpage, 1,tree->pagesize, f);
        fflush(f);


        fseek(f,pageposition,SEEK_SET);
        char rpage[tree->pagesize];
        fread(rpage,1,tree->pagesize,f);
        char cp = rpage[0];

        page->pagenumber = tree->numofpages;

        tree->numofpages++;
    }

    free(nodes);
    free(wpage);
}

void lsd_new_treepage(struct lsdtree * tree, struct lsdpage * page)
{
    page->pagenumber = -1;
    lsd_put_treepage(tree, page);
}


