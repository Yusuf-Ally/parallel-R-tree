#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "RTree.h"
#include "genRect.h"

#define GRID_MAX 1000000

int main(){

srand(time(NULL));

/// Create array of rectangle objects
struct Rect *rectList;
rectList = (struct Rect*)malloc(NUMRECTS*sizeof(struct Rect));

    if (rectList == NULL){
        printf("Error allocating memory.");
        exit(0);}

    for (int n = 0;n<NUMRECTS;n++){
        rectList[n] = generateRect(100,n);
        printf("%d %d %d %d ID:%d\n",rectList[n].x1,rectList[n].x2,rectList[n].y1,rectList[n].y2,rectList[n].rectID);
        }

///Sort rects by their low-x-coordinate
    sortX(rectList);

    printf("\nArray after sorting by x value : \n\n");
    for (int n = 0;n<NUMRECTS;n++){
        printf("%d %d %d %d ID:%d\n",rectList[n].x1,rectList[n].x2,rectList[n].y1,rectList[n].y2,rectList[n].rectID);
    }

    sortY(rectList);

    printf("\nArray after sorting by y value : \n\n");
    for (int n = 0;n<NUMRECTS;n++){
        printf("%d %d %d %d ID:%d\n",rectList[n].x1,rectList[n].x2,rectList[n].y1,rectList[n].y2,rectList[n].rectID);
    }

}


