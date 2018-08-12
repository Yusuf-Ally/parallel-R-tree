#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "RTree.h"
#include "genRect.h"

#define GRID_MAX 1000000

int main(){

///Create root node, with MBR of whole grid ie (1 000 000 x 1 000 000)
struct Node *root = Root(0,GRID_MAX,0,GRID_MAX);


printf("%d",root->nodeMBR.x2);
printf("\n");
srand(time(NULL));

    struct Rect testRect = generateRect(GRID_MAX);
    printf("%d",testRect.x1);
    printf("\n");
    printf("%d",testRect.x2);
    printf("\n");
    printf("%d",testRect.y1);
    printf("\n");
    printf("%d",testRect.y2);
    printf("\n");

}

