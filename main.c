#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include "lsdtree.h"
#include "kdtree.h"
#include "queue.h"


int is_point_contained(num_t * point, num_t * searchrange, int dim);
int search_count(num_t ** points, num_t * searchrange,int numofpoints, int dim);

int is_rect_intersecting(num_t * rect, num_t * searchrange, int dim);
int search_cts_count(num_t ** rects, num_t * searchrange, int numofrects, int dim);

int main()
{
    kdctsrangesearch();

    //int nooft = 20;
    //#pragma omp parallel
    //mptest();
    return 0;
}

void mptest()
{
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    printf("c = %d\n", my_rank);
}

void queuetest()
{

    srand((unsigned int)time(NULL));

    queue * q = queue_create();

    int data[10];

    for (int i = 0; i < 10; i++)
    {
        data[i] = rand() % 100;
        queue_push(q,&data[i]);
        printf ("d=%d\n", data[i]);
    }

    void * item;
    while ((item = queue_pop(q)) != NULL)
    {
        printf("%d\n", *(int*)item);
    }

    queue_free(q);
}


void kdctsrangesearch()
{
    int numofvalues = 1000000;
    int numofdims = 2;
    int numofctsdims = 2 * numofdims;
    int upper = 5000;
    int maxrectsize = 2000;

    srand((unsigned int)time(NULL));

    knum_t * arrvalues = malloc(sizeof(knum_t) *numofctsdims * numofvalues);
    knum_t ** values = malloc(sizeof(knum_t*) * numofvalues);

    for (int i = 0; i < numofvalues; i++)
    {
        for (int j = 0; j < numofdims; j++)
        {
            int pos = i * numofctsdims + j * 2;
            arrvalues[pos] = rand() % upper;
            arrvalues[pos+1] = arrvalues[pos] + rand() % maxrectsize;
        }

        values[i] = arrvalues + i * numofctsdims;
    }

    struct kdtree * tree = kd_create(values,numofctsdims, numofvalues);

    printf("tree created with %d values\n", numofvalues);

    knum_t searchrect[numofdims * 2];

    searchrect[0] = 3000;
    searchrect[1] = 4000;
    searchrect[2] = 3000;
    searchrect[3] = 4000;

    clock_t start = clock();
    int count1 = kd_cts_range_count(tree, searchrect);
    clock_t diff = clock() - start;
    int msec = diff;
    printf("n = %d  t = %d\n", count1, msec);


    start = clock();
    int count2 = kd_par_cts_range_count(tree,searchrect);
    diff = clock() - start;
    msec = diff;
    printf("n = %d  t = %d\n", count2,msec);


    start = clock();
    int count3 = search_cts_count(values, searchrect, numofvalues, numofdims);
    diff = clock() - start;
    msec = diff;
    printf("n = %d  t = %d\n", count3,msec);

    free(arrvalues);
    free(values);
}

void kdrangesearch()
{
    int numofvalues = 1000000;
    int numofdims = 2;
    int upper = 10000;
    srand((unsigned int)time(NULL));

    int arrsize = numofdims * numofvalues;

    knum_t * arrvalues = malloc(arrsize * sizeof(knum_t));

    for (int i = 0; i < arrsize; i++)
        arrvalues[i] = rand() % upper;

    knum_t **values = malloc(sizeof(knum_t *) * numofvalues);

    for (int i = 0; i < numofvalues; i++)
        values[i] = arrvalues + i * numofdims;

    struct kdtree * tree = kd_create(values,numofdims, numofvalues);

    printf("tree created with %d values\n", numofvalues);

    knum_t searchrange[numofdims * 2];

    searchrange[0] = 5000;
    searchrange[1] = 8000;
    searchrange[2] = 5000;
    searchrange[3] = 8000;

    clock_t start = clock();
    int count1 = kd_range_count(tree, searchrange);
    clock_t diff = clock() - start;
    int msec = diff;
    printf("n = %d  t = %d\n", count1, msec);


    start = clock();
    int count2 = kd_par_range_count(tree,searchrange);
    diff = clock() - start;
    msec = diff;
    printf("n = %d  t = %d\n", count2,msec);


    start = clock();
    int count3 = search_count(values, searchrange, numofvalues, numofdims);
    diff = clock() - start;
    msec = diff;
    printf("n = %d  t = %d\n", count3,msec);

    kd_free(tree);
    free(arrvalues);
    free(values);
}

void lsdctsrangesearch()
{
    int numofvalues = 100000;
    int numofdims = 2;
    int numofctsdims = 2 * numofdims;
    int upper = 5000;
    int maxrectsize = 2000;

    srand((unsigned int)time(NULL));

    num_t * arrvalues = malloc(sizeof(num_t) *numofctsdims * numofvalues);
    num_t ** values = malloc(sizeof(num_t*) * numofvalues);

    for (int i = 0; i < numofvalues; i++)
    {
        for (int j = 0; j < numofdims; j++)
        {
            int pos = i * numofctsdims + j * 2;
            arrvalues[pos] = rand() % upper;
            arrvalues[pos+1] = arrvalues[pos] + rand() % maxrectsize;
        }

        values[i] = arrvalues + i * numofctsdims;
    }

    struct lsdtree * tree = lsd_create(values,numofctsdims, numofvalues);

    printf("tree created with %d values\n", numofvalues);

    num_t searchrect[numofdims * 2];

    searchrect[0] = 3000;
    searchrect[1] = 4000;
    searchrect[2] = 3000;
    searchrect[3] = 4000;

    clock_t start = clock();
    int count1 = lsd_cts_range_count(tree, searchrect);
    clock_t diff = clock() - start;
    int msec = diff;
    printf("n = %d  t = %d\n", count1, msec);


    start = clock();
    int count2 = lsd_par_cts_range_count(tree,searchrect);
    diff = clock() - start;
    msec = diff;
    printf("n = %d  t = %d\n", count2,msec);


    start = clock();
    int count3 = search_cts_count(values, searchrect, numofvalues, numofdims);
    diff = clock() - start;
    msec = diff;
    printf("n = %d  t = %d\n", count3,msec);

    free(arrvalues);
    free(values);
}

void lsdrangesearch()
{
    int numofvalues = 40000;
    int numofdims = 2;
    int upper = 10000;
    srand((unsigned int)time(NULL));

    int arrsize = numofdims * numofvalues;

    num_t * arrvalues = malloc(arrsize * sizeof(num_t));

    for (int i = 0; i < arrsize; i++)
        arrvalues[i] = rand() % upper;

    num_t **values = malloc(sizeof(num_t *) * numofvalues);

    for (int i = 0; i < numofvalues; i++)
        values[i] = arrvalues + i * numofdims;

    struct lsdtree * tree = lsd_create(values,numofdims, numofvalues);

    printf("tree created with %d values\n", numofvalues);

    num_t searchrange[numofdims * 2];

    searchrange[0] = 5000;
    searchrange[1] = 8000;
    searchrange[2] = 5000;
    searchrange[3] = 8000;

    clock_t start = clock();
    int count1 = lsd_range_count(tree, searchrange);
    clock_t diff = clock() - start;
    int msec = diff;
    printf("n = %d  t = %d\n", count1, msec);


    start = clock();
    int count2 = lsd_par_range_count(tree,searchrange);
    diff = clock() - start;
    msec = diff;
    printf("n = %d  t = %d\n", count2,msec);


    start = clock();
    int count3 = search_count(values, searchrange, numofvalues, numofdims);
    diff = clock() - start;
    msec = diff;
    printf("n = %d  t = %d\n", count3,msec);

    lsd_free(tree);
    free(arrvalues);
    free(values);
}

int is_point_contained(num_t * point, num_t * searchrange, int dim)
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

int search_count(num_t ** points, num_t * searchrange,int numofpoints, int dim)
{
    int ret = 0;

   for (int i = 0; i < numofpoints; i++)
   {
        if (is_point_contained(points[i], searchrange, dim))
            ret++;
   }

   return ret;
}

int is_rect_intersecting(num_t * rect, num_t * searchrect, int dim)
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

int search_cts_count(num_t ** rects, num_t * searchrect, int numofrects, int dim)
{
    int ret = 0;

    for (int i= 0; i < numofrects; i++)
    {
        if (is_rect_intersecting(rects[i],searchrect, dim))
            ret++;
    }

    return ret;
}

