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
int numRects = 10;

/// Create array of rectangle objects
struct Rect *rectList;
rectList = (struct Rect*)malloc(numRects*sizeof(struct Rect));

    if (rectList == NULL){
        printf("Error allocating memory.");
        exit(0);}

    for (int n = 0;n<numRects;n++){
        rectList[n] = generateRect(100000,n);
        printf("%d %d %d %d ID:%d\n",rectList[n].x1,rectList[n].x2,rectList[n].y1,rectList[n].y2,rectList[n].rectID);}

///Sort rects by their low-x-coordinate
struct Rect temp;

    for (int i=0;i<numRects-1;i++){
        for (int j=i+1;j<numRects;j++){
            if (rectList[i].x1 > rectList[j].x2){
                temp = rectList[i];
                rectList[i] = rectList[j];
                rectList[j] = temp;}
                }
        }
printf("\nArray after sorting by x value : \n\n");
    for (int n = 0;n<numRects;n++){

    printf("%d %d %d %d ID:%d\n",rectList[n].x1,rectList[n].x2,rectList[n].y1,rectList[n].y2,rectList[n].rectID);}









}


