#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "RTree.h"
#include "genRect.h"

int main(){

srand(time(NULL));

/// Create array of rectangle objects
struct Rect *rectList;
rectList = genClusteredRectList(NUMRECTS, 4, GRID_MAX);
//rectList = genRectList(NUMRECTS);

    //printf("\nArray : \n\n");
    //for (int n = 0;n<NUMRECTS;n++){
    //    rectList[n] = generateRect(GRID_MAX,n);
    //}
double buildStart = omp_get_wtime();
	radixSort(rectList,NUMRECTS);

    sortY(rectList);

/// Build tree
    struct Node *root;
    root = buildTree(rectList);

double buildEnd = omp_get_wtime();
double buildTime = (buildEnd - buildStart);
printf("Time to build tree of %d rects is %lfs\n",NUMRECTS,buildTime);


/// Range search

    int num_of_average = 100;

    struct Rect query;
    query = ownRect(0,GRID_MAX/100*100,0,GRID_MAX,-1);
    
    double queryTime = 0;
    int intersectingRects = 0;
    for(int i = 0; i < num_of_average; i++)
    {
        double startQuery = omp_get_wtime();
        intersectingRects = RangeSearch(root,query);
        double endQuery = omp_get_wtime();
        queryTime += (double)(endQuery - startQuery) * 1000;
    }
    queryTime /= num_of_average;


    double queryTimeq = 0;
    int intersectingRectsq = 0;
    for(int i = 0; i < num_of_average; i++)
    {
        double startQueryq = omp_get_wtime();
        intersectingRectsq = RangeSearchQ(root,query);
        double endQueryq = omp_get_wtime();
        queryTimeq += (double)(endQueryq - startQueryq) * 1000;
    }
    queryTimeq /= num_of_average;
    double speedup = queryTime / queryTimeq;

    printf("\n\nFound %d rects inside query window\nTime to search: %lfms\n\nFound %d rects inside query window\nTime to search queue: %lfms\n",intersectingRects,queryTime,intersectingRectsq,queryTimeq);
    printf("\nSpeedup - %lf\n",speedup);






}


