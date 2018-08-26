#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "lsdtree.h"
#include "kdtree.h"
#include "queuea.h"
#include "queue.h"
#include "test.h"


int is_point_contained(num_t * point, num_t * searchrange, int dim);
int search_count(num_t ** points, num_t * searchrange,int numofpoints, int dim);

int is_rect_intersecting(num_t * rect, num_t * searchrange, int dim);
int search_cts_count(num_t ** rects, num_t * searchrange, int numofrects, int dim);

void queuetest();
void kdctsrangesearch();
void kdrangesearch();
void lsdctsrangesearch();
void lsdrangesearch();


int main()
{
    search_rect_size_test();

    return 0;
}

void queuetest()
{

    srand((unsigned int)time(NULL));

    queuea * q = queuea_create(1000);

    int data[10];

    for (int i = 0; i < 10; i++)
    {
        data[i] = rand() % 100;
        queuea_push(q,&data[i]);
        printf ("d=%d\n", data[i]);
    }

    void * item;
    while ((item = queuea_pop(q)) != NULL)
    {
        printf("r=%d\n", *(int*)item);
    }

    queuea_free(q);
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

    double start = omp_get_wtime();
    int count1 = kd_cts_range_count(tree, searchrect);
    double diff = omp_get_wtime() - start;
    double msec = diff;
    printf("n = %d  t = %f\n", count1, msec);


    start = omp_get_wtime();
    int count2 = kd_par_cts_range_count2(tree,searchrect);
    diff = omp_get_wtime() - start;
    msec = diff;
    printf("n = %d  t = %f\n", count2,msec);


    start = omp_get_wtime();
    int count3 = search_cts_count(values, searchrect, numofvalues, numofdims);
    diff = omp_get_wtime() - start;
    msec = diff;
    printf("n = %d  t = %f\n", count3,msec);

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

    searchrange[0] = 3000;
    searchrange[1] = 8000;
    searchrange[2] = 3000;
    searchrange[3] = 8000;

    double start = omp_get_wtime();
    int count1 = kd_range_count(tree, searchrange);
    double diff = omp_get_wtime() - start;
    double msec = diff;
    printf("n = %d  t = %f\n", count1, msec);


    start = omp_get_wtime();
    int count2 = kd_par_range_count2(tree,searchrange);
    diff = omp_get_wtime() - start;
    msec = diff;
    printf("n = %d  t = %f\n", count2,msec);


    start = omp_get_wtime();
    int count3 = search_count(values, searchrange, numofvalues, numofdims);
    diff = omp_get_wtime() - start;
    msec = diff;
    printf("n = %d  t = %f\n", count3,msec);

    kd_free(tree);
    free(arrvalues);
    free(values);
}

void lsdctsrangesearch()
{
    int numofvalues = 400000;
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

    double start = omp_get_wtime();
    int count1 = lsd_cts_range_count(tree, searchrect);
    double diff = omp_get_wtime() - start;
    double msec = diff;
    printf("n = %d  t = %f\n", count1, msec);


    start = omp_get_wtime();
    int count2 = lsd_par_cts_range_count2(tree,searchrect);
    diff = omp_get_wtime() - start;
    msec = diff;
    printf("n = %d  t = %f\n", count2,msec);


    start = omp_get_wtime();
    int count3 = search_cts_count(values, searchrect, numofvalues, numofdims);
    diff = omp_get_wtime() - start;
    msec = diff;
    printf("n = %d  t = %f\n", count3,msec);

    free(arrvalues);
    free(values);
}

void lsdrangesearch()
{
    int numofvalues = 1000000;
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

    searchrange[0] = 3000;
    searchrange[1] = 8000;
    searchrange[2] = 3000;
    searchrange[3] = 8000;

    double start = omp_get_wtime();
    int count1 = lsd_range_count(tree, searchrange);
    double diff = omp_get_wtime() - start;
    double msec = diff;
    printf("n = %d  t = %f\n", count1, msec);


    start = omp_get_wtime();
    int count2 = lsd_par_range_count2(tree,searchrange);
    diff = omp_get_wtime() - start;
    msec = diff;
    printf("n = %d  t = %f\n", count2,msec);


    start = omp_get_wtime();
    int count3 = search_count(values, searchrange, numofvalues, numofdims);
    diff = omp_get_wtime() - start;
    msec = diff;
    printf("n = %d  t = %f\n", count3,msec);

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

