#ifndef GENRECT_H
#define GENRECT_H

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "RTree.h"

struct Rect generateRect(int grid){

int maxLength = grid/100;    ///1% of grid length
int allowed = 0;
int rX;
int rY;
struct Rect rect;

int rWidth= 1+rand()%maxLength;
int rHeight = 1+rand()%maxLength;

while (allowed == 0){
    rX = 1+rand()%grid;
    rY = 1+rand()%grid;

    if ( (rX+rWidth < grid) && (rY+rHeight < grid) )
        {allowed = 1;}
}

rect.x2 = rX+rWidth;
rect.x1 = rX;

rect.y2 = rY+rHeight;
rect.y1 = rY;

return (rect);

}





#endif
